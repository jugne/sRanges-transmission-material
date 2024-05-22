rm(list = ls())
gc()

library(ggplot2)
library(dplyr)
library(coda)

####################################################
## Make sure the SA package is installed in BEAST ##
####################################################

##### SET PATHS ######
appLauncherPath <- "'/Applications/BEAST 2.6.7/bin/applauncher'"
logCombinerPath <- "'/Applications/BEAST 2.6.7/bin/logcombiner'"
branchRateAnalyserPath <- "~/Documents/Source/beast2.7/stratigraphic-ranges/out/artifacts/branch_rate_analyser/branchRateAnalyser.jar"
directionAnalyserPath <- 

sRangesPath <- "~/Documents/Source/beast2.7/sRanges-transmission-material/inference/bulkseq_allserial/no_removal_all_change/relaxed_clock/"
sRangesPriorPath <- "~/Documents/Source/beast2.7/sRanges-transmission-material/inference/bulkseq_allserial/no_removal_all_change/relaxed_clock/prior_runs/"
figure_path <- "~/Documents/Source/beast2.7/sRanges-transmission-material/inference/figures/"
figure_extension <- ".pdf" # e.g., "pdf", "png", "svg"


combine_logs <- function(path, pattern, extesion, burnin, seedInName=F){
  print(paste0(pattern, "_\\d.", extesion))
  if (seedInName){
    logs <- list.files(path=path, pattern = paste0(pattern, "_[0-9]\\.[0-9]*.", extesion))
  } else{
    logs <- list.files(path=path, pattern = paste0(pattern, "_[0-9]\\.", extesion))
  }
  cmd <- paste(logCombinerPath, "-b", burnin)
  for (log in logs){
    cmd <- paste(cmd, "-log", log)
  }
  cmd <- paste0(cmd, " -o ", pattern, ".combined")
  cmd <- paste0(cmd, ".", extesion)
  print(cmd)
  system(paste0("cd ",path ,";",cmd))
}

getSATable<- function(treeFile, outFile, analize=T){
  if (analize){
    cmd <- paste(appLauncherPath,'SampledAncestorTreeAnalyser -file', treeFile, ">", outFile)
    system(cmd)
  }
  # Determine the total number of lines in the file
  # total_lines <- length(readLines(outFile))
  # Skip the first two lines and read until (total_lines - 2)
  data <- read.table(outFile, sep = "\t", header = TRUE, skip = 2)#, nrows = total_lines - 3)
  
  old_names <- c("Burnside_Palaeudyptes", "Waimanu_tuatahi")
  new_names <- c("Burnside_Palaeeudyptes", "Muriwaimanu_tuatahi")
  for (i in 1:length(old_names)){
    id<-which(data$SA==old_names[i])
    if (length(id)!=0){
      data$SA[id] <- new_names[i]
    }
  }
  data$SA <- gsub("_last", "", gsub("_first", "", data$SA))
  data$SA <- gsub("_", " ", data$SA)
  n <- names(which(table(data$SA)==2))
  for (nn in n){
    ids <- which(data$SA==nn) 
    data <- data[-ids[1],]
  }
  return(data)
}


plot_clock_ranges_difference <- function(filePath, figurePath, width=10, height=10){
  
  # Open the connection to the file
  connection <- file(filePath, open = "r")
  
  # Read and ignore the first line (header)
  readLines(connection, n = 1)
  
  initial_size <- 100000
  rates_noRange <- vector("numeric", initial_size)
  weight_noRange <- vector("numeric", initial_size)
  weight_range <- vector("numeric", initial_size)
  rates_range <- vector("numeric", initial_size)
  count_nr <- 0
  count_r <- 0
  # Process each line
  while (length(line <- readLines(connection, n = 1)) > 0) {
    
    
    # Split the line by tabs
    columns <- strsplit(line, "\t")[[1]]
    
    # Assuming there are at least three columns in the file
    if (length(columns) >= 3) {
      col1 <- columns[1]
      col2 <- columns[2]
      col3 <- columns[3]
      
      if (col3=="1"){
        count_r <- count_r + 1
        # Resize the vector if needed
        if (count_r > length(weight_range)) {
          weight_range <- c(weight_range, vector("numeric", initial_size))
          rates_range <- c(rates_range, vector("numeric", initial_size))
        }
        rates_range[count_r] <- as.numeric(col1)
        weight_range[count_r] <- as.numeric(col2)
        
      } else {
        count_nr <- count_nr + 1
        if (count_nr > length(weight_noRange)) {
          weight_noRange <- c(weight_noRange, vector("numeric", initial_size))
          rates_noRange <- c(rates_noRange, vector("numeric", initial_size))
        }
        rates_noRange[count_nr] <- as.numeric(col1)
        weight_noRange[count_nr] <- as.numeric(col2)
      }
    }
  }
  close(connection)
  rates_noRange <- rates_noRange[1:count_nr]
  weight_noRange <- weight_noRange[1:count_nr]
  
  weight_range <- weight_range[1:count_r]
  rates_range <- rates_range[1:count_r]
  
  # id_range <- which(rates_range>1e-5)
  # id_noRange <- which(rates_noRange>1e-5)
  # dd <- data.frame(range=c(rep("range", length(id_range)), rep("non range", length(id_noRange))), vals=c(rates_range[id_range], rates_noRange[id_noRange]), weights=c(weight_range[id_range], weight_noRange[id_noRange]))
  dd_ <- data.frame(range=c(rep("range", length(rates_range)), rep("non range", length(rates_noRange))), vals=c(rates_range, rates_noRange), weights=c(weight_range, weight_noRange))

  # The palette with grey:
  cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  p <- ggplot(dd_) + geom_density(aes(log10(vals), weight=weights, color=range, fill=range, alpha=range)) 
  p1<- p + theme_bw(base_size = 25) + 
    guides(fill = guide_legend(ncol=1, byrow=T)) +
    theme(legend.position = "top", 
          # legend.justification = c("right", "top"), 
          legend.box.just = "left", 
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          legend.box.spacing = margin(6),
          legend.title = element_blank(),
          legend.text=element_text(size=19),
          legend.spacing.y = unit(0.5, "lines")) + 
    scale_color_manual(values=cbp2) +
    scale_fill_manual(values=cbp2)+
    scale_alpha_manual(values=c("range"=0.3, "non range"=1))+
    labs( x = "log10(clock rate)",
          y = "density")
  ggsave(figurePath, p1, width=width, height=height)
  
}


###### SET THRESHOLDS #######
sa_percent <- 5
burnin <- 30

################# Combine logs ############################

# combine and read logs for FBD run
if (!file.exists(paste0(sRangesPath, "SA_Georgia_L2_singleSample.combined.log"))){
  combine_logs(sRangesPath, "SA_Georgia_L2_singleSample", "log", burnin)
}
SA_relaxed_all <- read.table(paste0(sRangesPath, "SA_Georgia_L2_singleSample.combined.log"), sep = "\t", header = TRUE)

# combine and read logs for SRFBD run
if (!file.exists(paste0(sRangesPath, "Georgia_L2_allserial.combined.log"))){
  combine_logs(sRangesPath, "Georgia_L2_allserial", "log", burnin, seed=T)
}
sRange_relaxed_all <- read.table(paste0(sRangesPath, "Georgia_L2_allserial.combined.log"), sep = "\t", header = TRUE)


################# Combine trees ############################

# combine and read logs for FBD run
if (!file.exists(paste0(sRangesPath, "SA_Georgia_L2_singleSample.combined.trees"))){
  combine_logs(sRangesPath, "SA_Georgia_L2_singleSample", "trees", burnin)
}

# combine and read logs for SRFBD run
if (!file.exists(paste0(sRangesPath, "Georgia_L2_allserial.combined.trees"))){
  combine_logs(sRangesPath, "Georgia_L2_allserial", "trees", burnin, seed=T)
}


cmd <- paste0("java -jar ", branchRateAnalyserPath, " -tree ", sRangesPath,
              "Georgia_L2_allserial.combined.trees -out ", sRangesPath,
                          "Georgia_L2_allserial.combined.branchRate.log", " -burnin 0 -normalise True")
if (!file.exists(paste0(sRangesPath,
                        "Georgia_L2_allserial.combined.branchRate.log"))){
  system(cmd)
}

plot_clock_ranges_difference(paste0(sRangesPath,"Georgia_L2_allserial.combined.branchRate.log"), 
                             paste0(figure_path, "srfbd_clock_rangeDiff", figure_extension),
                             width=6, height=8)

cmd <- paste0("java -jar ", branchRateAnalyserPath, " -tree ", sRangesPriorPath,
              "Georgia_L2_allserial.233457685745.trees -out ", sRangesPriorPath,
              "Georgia_L2_allserial.prior.branchRate.log", " -burnin 10 -normalise True")
if (!file.exists(paste0(sRangesPriorPath,
                        "Georgia_L2_allserial.prior.branchRate.log"))){
  system(cmd)
}

plot_clock_ranges_difference(paste0(sRangesPriorPath,"Georgia_L2_allserial.prior.branchRate.log"), 
                             paste0(figure_path, "prior_srfbd_clock_rangeDiff", figure_extension),
                             width=6, height=8)

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

# The palette with grey:
cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

sa_r_l <- nrow(SA_relaxed_all)
sR_r_l <- nrow(sRange_relaxed_all)

model <- c(rep("FBD", sa_r_l), rep("SRFBD", sR_r_l))

R0_r_all <- data.frame(r0_1=c(SA_relaxed_all$R0SVEpi.i0, sRange_relaxed_all$R0.1),
                       r0_2=c(SA_relaxed_all$R0SVEpi.i1, sRange_relaxed_all$R0.2),
                       b_1 = c(SA_relaxed_all$becomeUninfectiousRateSVEpi.i0, sRange_relaxed_all$becomeUninfectiousRate.1),
                       b_2 = c(SA_relaxed_all$becomeUninfectiousRateSVEpi.i1, sRange_relaxed_all$becomeUninfectiousRate.2),
                       s = c(SA_relaxed_all$samplingProportionSVEpi.i1, sRange_relaxed_all$samplingProportion.2),
                       orig = c(SA_relaxed_all$originBDMMPrime, sRange_relaxed_all$origin),
                       treeHeight = c(SA_relaxed_all$TreeHeight, sRange_relaxed_all$TreeHeight),
                       clock = c(SA_relaxed_all$ORCRatesStat.mean, sRange_relaxed_all$ORCRatesStat.c.Georgia_L2.mean),
                       clock_v = c(SA_relaxed_all$ORCRatesStat.variance, sRange_relaxed_all$ORCRatesStat.c.Georgia_L2.variance),
                       clock_coefV = c(SA_relaxed_all$ORCRatesStat.coefficientOfVariation, sRange_relaxed_all$ORCRatesStat.c.Georgia_L2.coefficientOfVariation),
                       model=model)

plot_params <-function(data, x_="model", y_, xlabel, labels, figName, title=F, width=10, height=10, prior=F, fbd=F){
  p<-ggplot(data, aes(x=.data[[y_]], fill=.data[[x_]], color=.data[[x_]])) +
    ylab("Density") +
    xlab(xlabel) +#aes(x=.data[[x_]], y=.data[[y_]],fill=.data[[x_]])) + ylab(label) + xlab("") +
    theme_bw(base_size = 25) + 
    guides(fill = guide_legend(ncol=1, byrow=T))+
    theme(legend.position="top",
      #legend.position = c(0.95, 0.95),
          # legend.justification = c("right", "top"), 
          legend.box.just = "left", 
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          legend.box.spacing = margin(6),
          legend.title = element_blank(),
          legend.text=element_text(size=19),
          legend.spacing.y = unit(0.5, "lines"),
      # legend.spacing.x = unit(6, "lines")
          # legend.key.size = unit(1, "cm")
          # axis.title.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.x=element_blank()
          )+
    geom_density(alpha=0.7)+
    # geom_violin(trim=T, alpha=0.7, color="grey",draw_quantiles = c(0.25, 0.5, 0.75))+
    # scale_fill_manual(values = cbp2[c(6,1)], labels = labels)+
    # scale_color_manual(values = cbp2[c(6,1)], labels = labels)
    scale_fill_manual(values = cbp2[c(6,7)], labels = labels)+
    scale_color_manual(values = cbp2[c(6,7)], labels = labels)
  if (prior & fbd){
    p <- p +  scale_fill_manual(values = cbp2[c(6,1)], labels = labels)+
       scale_color_manual(values = cbp2[c(6,1)], labels = labels)
  }
  if (prior & !fbd){
    p <- p +  scale_fill_manual(values = cbp2[c(7,1)], labels = labels)+
      scale_color_manual(values = cbp2[c(7,1)], labels = labels)
  }
  if (title){
    p <- p + facet_grid(. ~ title)
  }
  ggsave(figName, p, width=width, height=height)
}

labels <- c("Reproductive number",
            "Reproductive number",
            "Becoming uninfectious rate",
            "Becoming uninfectious rate",
            "Smpling proportion",
            "Origin",
            "Tree Height",
            "log10(clock rate)",
            "Clock rate variance",
            "Clock rate coeff. of variation")
titles <- c(T,T,T,T,F,F,F,F,F,F)
title_val <- c("Before sampling start", "After sampling start", "Before sampling start", "After sampling start", "","","","","","") 
figNames <- c("r1","r2","b1","b2","s","origin","height","clock_mean","clock_var","clock_cv")
n <- ncol(R0_r_all)
for (i in 1:(n-1)){
  figName <- paste0(figure_path,figNames[i],figure_extension)
  hpd1<- HPDinterval(as.mcmc(R0_r_all[which(R0_r_all$model=="SRFBD"), i]))
  l1 <- paste0("SRFBD, ", round(median(R0_r_all[which(R0_r_all$model=="SRFBD"), i]),2),
               " [",round(hpd1[1],2), ",",round(hpd1[2],2),"]")
  hpd2<- HPDinterval(as.mcmc(R0_r_all[which(R0_r_all$model=="FBD"), i]))
  l2 <- paste0("FBD, ", round(median(R0_r_all[which(R0_r_all$model=="FBD"), i]),2),
               " [",round(hpd2[1],2), ",",round(hpd2[2],2),"]")
  if (i==8){
    l1 <- paste0("SRFBD, ", formatC(median(R0_r_all[which(R0_r_all$model=="SRFBD"), i]), format = "G", digits = 2),
                 " [",formatC(hpd1[1], format = "G", digits = 2), ",",formatC(hpd1[2], format = "G", digits = 2),"]")
    l2 <- paste0("FBD, ", formatC(median(R0_r_all[which(R0_r_all$model=="FBD"), i]), format = "G", digits = 2),
                 " [",formatC(hpd2[1], format = "G", digits = 2), ",",formatC(hpd2[2], format = "G", digits = 2),"]")
    R0_r_all$clock <- log10(R0_r_all$clock)
  }
  if (titles[i]){
    d <- R0_r_all
    d$title <- title_val[i]
    plot_params(d, y_=colnames(R0_r_all)[i], xlabel=labels[i], labels=c(l2, l1),
                figName=figName, title=T, width=6, height=8)
  } else {
    plot_params(R0_r_all, y_=colnames(R0_r_all)[i], xlabel=labels[i], labels=c(l2, l1),
                figName=figName, title=F, width=6, height=8)
  }
}

sRange_relaxed_all_prior <- read.table(paste0(sRangesPriorPath, "Georgia_L2_allserial.233457685745.log"), sep = "\t", header = TRUE)
sRange_relaxed_all_prior <- tail(sRange_relaxed_all_prior, -(round(nrow(sRange_relaxed_all_prior)*0.1)))

sR_prior_r_l <- nrow(sRange_relaxed_all_prior)
sR_r_l <- nrow(sRange_relaxed_all)

model <- c(rep("SRFBD Prior", sR_prior_r_l), rep("SRFBD", sR_r_l))

R0_r_all_prior <- data.frame(r0_1=c(sRange_relaxed_all_prior$R0.1, sRange_relaxed_all$R0.1),
                             r0_2=c(sRange_relaxed_all_prior$R0.2, sRange_relaxed_all$R0.2),
                             b_1 = c(sRange_relaxed_all_prior$becomeUninfectiousRate.1, sRange_relaxed_all$becomeUninfectiousRate.1),
                             b_2 = c(sRange_relaxed_all_prior$becomeUninfectiousRate.2, sRange_relaxed_all$becomeUninfectiousRate.2),
                             s = c(sRange_relaxed_all_prior$samplingProportion.2, sRange_relaxed_all$samplingProportion.2),
                             orig = c(sRange_relaxed_all_prior$origin, sRange_relaxed_all$origin),
                             treeHeight = c(sRange_relaxed_all_prior$TreeHeight, sRange_relaxed_all$TreeHeight),
                             clock = c(sRange_relaxed_all_prior$ORCRatesStat.c.Georgia_L2.mean, sRange_relaxed_all$ORCRatesStat.c.Georgia_L2.mean),
                             clock_v = c(sRange_relaxed_all_prior$ORCRatesStat.c.Georgia_L2.variance, sRange_relaxed_all$ORCRatesStat.c.Georgia_L2.variance),
                             clock_coefV = c(sRange_relaxed_all_prior$ORCRatesStat.c.Georgia_L2.coefficientOfVariation, sRange_relaxed_all$ORCRatesStat.c.Georgia_L2.coefficientOfVariation),
                             model=model)
n <- ncol(R0_r_all_prior)
dir.create(paste0(figure_path,"prior"))

for (i in 1:(n-1)){
  figName <-  paste0(figure_path,"prior/",figNames[i],figure_extension)
  hpd1<- HPDinterval(as.mcmc(R0_r_all_prior[which(R0_r_all_prior$model=="SRFBD"), i]))
  l2 <- paste0("SRFBD, ", round(median(R0_r_all_prior[which(R0_r_all_prior$model=="SRFBD"), i]),2),
               " [",round(hpd1[1],2), ",",round(hpd1[2],2),"]")
  hpd2<- HPDinterval(as.mcmc(R0_r_all_prior[which(R0_r_all_prior$model=="SRFBD Prior"), i]))
  l1 <- paste0("SRFBD Prior, ", round(median(R0_r_all_prior[which(R0_r_all_prior$model=="SRFBD Prior"), i]),2),
               " [",round(hpd2[1],2), ",",round(hpd2[2],2),"]")
  if (i==8){
    l2 <- paste0("SRFBD, ", formatC(median(R0_r_all_prior[which(R0_r_all_prior$model=="SRFBD"), i]), format = "G", digits = 2),
                 " [",formatC(hpd1[1], format = "G", digits = 2), ",",formatC(hpd1[2], format = "G", digits = 2),"]")
    l1 <- paste0("SRFBD Prior, ", formatC(median(R0_r_all_prior[which(R0_r_all_prior$model=="SRFBD Prior"), i]), format = "G", digits = 2),
                 " [",formatC(hpd2[1], format = "G", digits = 2), ",",formatC(hpd2[2], format = "G", digits = 2),"]")
    R0_r_all_prior$clock <- log10(R0_r_all_prior$clock)
  }
  if (titles[i]){
    d <- R0_r_all_prior
    d$title <- title_val[i]
    plot_params(d, y_=colnames(R0_r_all_prior)[i], xlabel=labels[i], labels=c(l2, l1),
                figName=figName, title=T, width=6, height=8, prior=T)
  } else {
    plot_params(R0_r_all_prior, y_=colnames(R0_r_all_prior)[i], xlabel=labels[i], labels=c(l2, l1),
                figName=figName, title=F, width=6, height=8, prior=T)
  }
}



FBD_relaxed_all_prior <- read.table(paste0(sRangesPriorPath, "SA_Georgia_L2_singleSample.log"), sep = "\t", header = TRUE)
FBD_relaxed_all_prior <- tail(FBD_relaxed_all_prior, -(round(nrow(FBD_relaxed_all_prior)*0.1)))

FBD_prior_r_l <- nrow(FBD_relaxed_all_prior)
FBD_r_l <- nrow(SA_relaxed_all)

model <- c(rep("FBD Prior", FBD_prior_r_l), rep("FBD", FBD_r_l))

R0_r_all_fbd_prior <- data.frame(r0_1=c(FBD_relaxed_all_prior$R0SVEpi.i0, SA_relaxed_all$R0SVEpi.i0),
                       r0_2=c(FBD_relaxed_all_prior$R0SVEpi.i1, SA_relaxed_all$R0SVEpi.i1),
                       b_1 = c(FBD_relaxed_all_prior$becomeUninfectiousRateSVEpi.i0, SA_relaxed_all$becomeUninfectiousRateSVEpi.i0),
                       b_2 = c(FBD_relaxed_all_prior$becomeUninfectiousRateSVEpi.i1, SA_relaxed_all$becomeUninfectiousRateSVEpi.i1),
                       s = c(FBD_relaxed_all_prior$samplingProportionSVEpi.i1, SA_relaxed_all$samplingProportionSVEpi.i1),
                       orig = c(FBD_relaxed_all_prior$originBDMMPrime, SA_relaxed_all$originBDMMPrime),
                       treeHeight = c(FBD_relaxed_all_prior$TreeHeight, SA_relaxed_all$TreeHeight),
                       clock = c(FBD_relaxed_all_prior$ORCRatesStat.mean, SA_relaxed_all$ORCRatesStat.mean),
                       clock_v = c(FBD_relaxed_all_prior$ORCRatesStat.variance, SA_relaxed_all$ORCRatesStat.variance),
                       clock_coefV = c(FBD_relaxed_all_prior$ORCRatesStat.coefficientOfVariation, SA_relaxed_all$ORCRatesStat.coefficientOfVariation),
                       model=model)
n <- ncol(R0_r_all_fbd_prior)
dir.create(paste0(figure_path,"prior"))

for (i in 1:(n-1)){
  figName <-  paste0(figure_path,"prior/fbd_",figNames[i],figure_extension)
  hpd1<- HPDinterval(as.mcmc(R0_r_all_fbd_prior[which(R0_r_all_fbd_prior$model=="FBD"), i]))
  l2 <- paste0("FBD, ", round(median(R0_r_all_fbd_prior[which(R0_r_all_fbd_prior$model=="FBD"), i]),2),
               " [",round(hpd1[1],2), ",",round(hpd1[2],2),"]")
  hpd2<- HPDinterval(as.mcmc(R0_r_all_fbd_prior[which(R0_r_all_fbd_prior$model=="FBD Prior"), i]))
  l1 <- paste0("FBD Prior, ", round(median(R0_r_all_fbd_prior[which(R0_r_all_fbd_prior$model=="FBD Prior"), i]),2),
               " [",round(hpd2[1],2), ",",round(hpd2[2],2),"]")
  if (i==8){
    l2 <- paste0("FBD, ", formatC(median(R0_r_all_fbd_prior[which(R0_r_all_fbd_prior$model=="FBD"), i]), format = "G", digits = 2),
                 " [",formatC(hpd1[1], format = "G", digits = 2), ",",formatC(hpd1[2], format = "G", digits = 2),"]")
    l1 <- paste0("FBD Prior, ", formatC(median(R0_r_all_fbd_prior[which(R0_r_all_fbd_prior$model=="FBD Prior"), i]), format = "G", digits = 2),
                 " [",formatC(hpd2[1], format = "G", digits = 2), ",",formatC(hpd2[2], format = "G", digits = 2),"]")
    R0_r_all_fbd_prior$clock<-log10(R0_r_all_fbd_prior$clock)
  }
  if (titles[i]){
    d <- R0_r_all_fbd_prior
    d$title <- title_val[i]
    plot_params(d, y_=colnames(R0_r_all_fbd_prior)[i], xlabel=labels[i], labels=c(l2, l1),
                figName=figName, title=T, width=6, height=8, prior=T, fbd=T)
  } else {
    plot_params(R0_r_all_fbd_prior, y_=colnames(R0_r_all_fbd_prior)[i], xlabel=labels[i], labels=c(l2, l1),
                figName=figName, title=F, width=6, height=8, prior=T, fbd=T)
  }
}


ii






# for (i in 1:(n-1)){
#   figName <- paste0(figure_path,"prior/",figNames[i],figure_extension)
#   if (titles[i]){
#     d <- R0_r_all_prior
#     d$title <- title_val[i]
#     plot_params(d, y_=colnames(R0_r_all)[i], label=labels[i], figName=figName,
#                 title=T, width=6, height=8)
#   } else {
#     plot_params(R0_r_all_prior, y_=colnames(R0_r_all)[i], label=labels[i],
#                 figName=figName, title=F, width=6, height=8)
#   }
# }

library(coda)

log <-  read.table(paste0(sRangesPath, "Georgia_L2_allserial.combined.transmission_2.log"),
                   header=TRUE, sep="\t")
# burnIn <- 0.1
# log <- tail(log, -round(nrow(log)*burnIn)) # take out 10% burnin
log<-log[,1:(ncol(log)-1)]
columnsWithValue <- colnames(log)[colSums(log != "0.0,0.0") != 0]
# we don't need "Sample" column
columnsWithValue <- columnsWithValue[2:length(columnsWithValue)]
log <- log[, columnsWithValue]

n <- length(columnsWithValue)
d <- data.frame(from=character(n), to=character(n), prob=numeric(n), length_median=numeric(n),
                length_hpd_l=numeric(n), length_hpd_h=numeric(n), nodes_median=numeric(n),
                nodes_hpd_l=numeric(n), nodes_hpd_h=numeric(n))
# Split each string by ":"
split_strings <- strsplit(columnsWithValue, "\\.")
# Extract the first and second elements
d$from <- sapply(split_strings, function(x) x[1])
d$to <- sapply(split_strings, function(x) x[2])

for (i in 1:n){
  split_strings <- strsplit(log[,columnsWithValue[i]], ",")
  length <- sapply(split_strings, function(x) as.double(x[1]))
  nodes <- sapply(split_strings, function(x) as.double(x[2]))
  
  d$prob[i] <- length(length)
  idx <- which(length != 0)
  length <- length[idx] # no length means no speciation between the two species recorded.
  nodes <- nodes[idx] # 0 intermediate nodes are legit
  d$prob[i] <- length(length)/d$prob[i]
  
  if (length(idx)>1){
    d$length_median[i] <- median(length)
    length_hpd <- HPDinterval(as.mcmc(length))
    d$length_hpd_l[i] <- length_hpd[1]
    d$length_hpd_h[i] <- length_hpd[2]
    
    d$nodes_median[i] <- median(nodes)
    nodes_hpd <- HPDinterval(as.mcmc(nodes))
    d$nodes_hpd_l[i] <- nodes_hpd[1]
    d$nodes_hpd_h[i] <- nodes_hpd[2]
  }
}

d <- d[order(d$prob,decreasing = TRUE),]
head(d, 30)

d_20 <- d[which(d$prob>0.30),]
d_20$prob <- round(d_20$prob, 3)
write.csv(d_20, paste0(figure_path, "directed_transmission_2.csv"))
