rm(list = ls())
gc()

library(stringr)
library(deeptime)
library(coda)
library(ggridges)
library(ggplot2)

addToDistr <- function(log){
  if ((max(l)-min(l))>=0.5){
    
  }
}


wd <- "~/Documents/Source/beast2.7/sRanges-transmission-material/hiv/"
setwd(wd)
log <- read.table("inference/hiv_relaxedStart_env.tipAge.log", header = T)
burnIn <- 0.1
log <- tail(log, -round(nrow(log)*burnIn)) 
colnames <- colnames(log)[2:dim(log)[2]]
ranges<-sapply(strsplit(colnames,"_"), `[`, 1)
ranges_unique <- unique(ranges)

sp <- c("A_first", "B_first", "C_first", "D_first", "E_first", "F_first", "G_first", "H_first",
        "I_first", "K_first", "L_first", 
        "A_last", "B_last", "C_last", "D_last", "E_last", "F_last", "G_last", "H_last", "I_last", "K_last")

start <- c(30,30,14.2724161533196,12.7954634337114,7.54379562043796,12.913179507235,
         5.95698083691826,12.5080172076652,14.9130732375086,3.5823754789272,14.2724161533196,
         1.99270072992701,11.3794483049063,4.0456204379562,4.77372262773723,1.85583941605839,
         2.28336755646817,5.71059835745014,9.60378297660528,8.23706542567753,3.2512315270936)

end <- c(17.6465753424658,17.6657534246575,13.9835616438356,12.5123287671233,6.91780821917808,
         12.4149393820884,5.71059835745014,12.4395776300352,14.6776180698152,3.29228243021346,1.80656934306569,
         0,10.9979466119097,0,0,0,0,0,0,0,0)

priors <- data.frame(species=sp, start=start, end=end)


all_string_listwise <- unlist(lapply(ranges, unique))
ranges_single <- names(which(table(all_string_listwise)==1))
ranges <- names(which(table(all_string_listwise)>1))
start_r <-c()
start_r_low <- c()
start_r_high <- c()
start_r_min <- c()
start_r_max <- c()
start_r_iqr <- c()
start_r_iqr_prior <- c()
end_r <-c()
end_r_low <- c()
end_r_high <- c()
end_r_min <- c()
end_r_max <- c()
end_r_iqr <- c()
end_r_iqr_prior <- c()


distr <- c()
endpoint <- c()
range <- c()
range_s <- c()
range_e <- c()
median_s <- c()
median_s_prior<-c()
median_e <- c()
median_e_prior<-c()

# for (r in ranges_single){
#   if ((max(log[,paste0(r,"_first")])-min(log[,paste0(r,"_first")]))>0.1){
#     start_l <- length(log[,paste0(r,"_first")])
#     distr <- c(distr, log[,paste0(r,"_first")])
#     endpoint <- c(endpoint, rep("start", start_l))
#     range <- c(range, rep(r, start_l))
#     
#     id <- which(priors$species==paste0(r,"_first"))
#     un <- runif(5000, min=priors$end[id], max=priors$start[id])
#     start_r_iqr <- c(start_r_iqr, IQR(log[,paste0(r,"_first")]))
#     start_r_iqr_prior <- c(start_r_iqr_prior, IQR(un))
#     
#     # end_r_iqr <- c(end_r_iqr, NA)
#     # end_r_iqr_prior <- c(end_r_iqr_prior, NA)
#     range_s <- c(range_s, r)
#     median_s<- c(median_s, median(log[,paste0(r,"_first")]))
#     median_s_prior<- c(median_s_prior, median(un))
#     
#   }
#   start_r <- c(start_r, median(log[,paste0(r,"_first")])) 
#   end_r <- c(end_r, NA)
#   hpd_f<- HPDinterval(as.mcmc(log[,paste0(r,"_first")]))
#   start_r_low <- c(start_r_low, hpd_f[1])
#   start_r_high <- c(start_r_high, hpd_f[2])
#   start_r_min <- c(start_r_min, min(log[,paste0(r,"_first")]))
#   start_r_max <- c(start_r_max, max(log[,paste0(r,"_first")]))
#   
#   
#   
#   # hpd_l<- HPDinterval(as.mcmc(log[,paste0(r,"_last")]))
#   end_r_low<- c(end_r_low, NA)
#   end_r_high<- c(end_r_high, NA)
#   end_r_min <- c(end_r_min, NA)
#   end_r_max <- c(end_r_max, NA)
#   
# }
  # start_l <- length(log[,paste0(r,"_first")])
  # end_l <- length(log[,paste0(r,"_last")])
  # distr <- c(distr, log[,paste0(r,"_first")], log[,paste0(r,"_last")])
  # endpoint <- c(endpoint, rep("start", start_l), rep("end", end_l))
  # range <- c(range, rep(r, (start_l+end_l)))
  
for (r in ranges){
  if ((max(log[,paste0(r,"_first")])-min(log[,paste0(r,"_first")]))>0.01){
    start_l <- length(log[,paste0(r,"_first")])
    distr <- c(distr, log[,paste0(r,"_first")])
    endpoint <- c(endpoint, rep("start", start_l))
    range <- c(range, rep(r, start_l))
    
    id <- which(priors$species==paste0(r,"_first"))
    un <- runif(5000, min=priors$end[id], max=priors$start[id])
    start_r_iqr <- c(start_r_iqr, IQR(log[,paste0(r,"_first")]))
    start_r_iqr_prior <- c(start_r_iqr_prior, IQR(un))
    range_s <- c(range_s, r)
    median_s<- c(median_s, median(log[,paste0(r,"_first")]))
    median_s_prior<- c(median_s_prior, median(un))
    
  }
  if ((max(log[,paste0(r,"_last")])-min(log[,paste0(r,"_last")]))>0.01){
    end_l <- length(log[,paste0(r,"_last")])
    distr <- c(distr, log[,paste0(r,"_last")])
    endpoint <- c(endpoint, rep("end", end_l))
    range <- c(range, rep(r, end_l))
    
    id <- which(priors$species==paste0(r,"_last"))
    un <- runif(5000, min=priors$end[id], max=priors$start[id])
    end_r_iqr <- c(end_r_iqr, IQR(log[,paste0(r,"_last")]))
    end_r_iqr_prior <- c(end_r_iqr_prior, IQR(un))
    range_e <- c(range_e, r)
    median_e<- c(median_e, median(log[,paste0(r,"_last")]))
    median_e_prior<- c(median_e_prior, median(un))
  }
  
  start_r <- c(start_r, median(log[,paste0(r,"_first")]))
  end_r <- c(end_r, median(log[,paste0(r,"_last")]))
  
  hpd_f<- HPDinterval(as.mcmc(log[,paste0(r,"_first")]))
  start_r_low <- c(start_r_low, hpd_f[1])
  start_r_high <- c(start_r_high, hpd_f[2])
  start_r_min <- c(start_r_min, min(log[,paste0(r,"_first")]))
  start_r_max <- c(start_r_max, max(log[,paste0(r,"_first")]))

  
  hpd_l<- HPDinterval(as.mcmc(log[,paste0(r,"_last")]))
  end_r_low<- c(end_r_low, hpd_l[1])
  end_r_high<- c(end_r_high, hpd_l[2])
  end_r_min <- c(end_r_min, min(log[,paste0(r,"_last")]))
  end_r_max <- c(end_r_max, max(log[,paste0(r,"_last")]))

}



ranges_full_s <- data.frame(range=gsub("_", " ", range_s),
                          start_r_iqr=start_r_iqr,
                          start_r_iqr_prior=start_r_iqr_prior,
                          median_s_prior=median_s_prior,
                          median_s=median_s)
ranges_full_e <- data.frame(range=gsub("_", " ", range_e),
                            end_r_iqr=end_r_iqr,
                            end_r_iqr_prior=end_r_iqr_prior,
                            median_e_prior=median_e_prior,
                            median_e=median_e)



ranges_full <- data.frame(distr=distr, endpoint=endpoint, range=range)


length(end_r)<-length(start_r)
ranges_stat <- data.frame(range=c(ranges_single,ranges),
                          median_start=start_r,
                          median_end=end_r,
                          start_r_min=start_r_min,
                          start_r_low=start_r_low,
                          start_r_max=start_r_max,
                          start_r_high=start_r_high,
                          end_r_min=end_r_min,
                          end_r_low=end_r_low,
                          end_r_max=end_r_max,
                          end_r_high=end_r_high)




p<-ggplot(ranges_full_s, aes(x = range, y = start_r_iqr_prior/start_r_iqr)) +
  geom_point(aes(color = "First"), size = 5, alpha=0.7) + 
  geom_point(data=ranges_full_e, 
             aes(x = range, y = end_r_iqr_prior/end_r_iqr, color = "Last"), 
             size = 5, alpha=0.7) +
  # geom_point(data=ranges_full_s, 
  #            aes(x = range, y = median_s/median_s_prior, color = "Median Start"), 
  #            size = 3, alpha=0.7)+
  # geom_point(data=ranges_full_e, 
  #            aes(x = range, y = median_e/median_e_prior, color = "Median Last"), 
  #            size = 3, alpha=0.7) +
  geom_hline(yintercept = 1)+ 
  ylab("IQR ratio (prior/posterior)") + 
  xlab("Species") + 
  scale_color_manual(values = c("First" = "#FF7F50", "Last" = "#43AA8B")) +
  theme_minimal()  + labs(color="Sample") +
  theme(legend.position = c(0.80, 0.85),
        legend.justification = c("right", "bottom"),text=element_text(size=21)) + coord_flip()
      
ggsave("figures/IQR_.pdf", width=8, height=11)



priors$endpoint<- str_split_fixed(priors$species, "_", 2)[,2]
priors$species<- str_split_fixed(priors$species, "_", 2)[,1]
ranges_prior_full <- data.frame(distr=numeric(), endpoint=character(), range=character())
for (i in 1:nrow(priors)){
  # e<-"end"
  # if (priors$endpoint[i]=="first"){
  #   e <- "start"
  # }
  ddd<- data.frame(distr=runif(100000, min=priors$end[i], max=priors$start[i]),
                   endpoint=rep(priors$endpoint[i], 100000), range=rep(priors$species[i], 100000))
  ranges_prior_full <- rbind(ranges_prior_full,ddd)
}


ranges_full$endpoint<-str_replace_all(ranges_full$endpoint, "start", "first")
ranges_full$endpoint<-str_replace_all(ranges_full$endpoint, "end", "last")

ranges_full$mode<-"posterior"
ranges_prior_full$mode<-"prior"
d4<- rbind(ranges_full,ranges_prior_full)
d4$distr<- 2006.226-d4$distr
p<-ggplot(d4) + theme_bw(base_size = 19)+
  geom_density(aes(distr, fill=mode), alpha=0.5, color=NA)+ 
  facet_wrap(range ~ endpoint, scales = "free", ncol=6)  +xlab("") +
  theme(axis.text.x = element_text(angle=-60, vjust=0.5), legend.position = c(0.7, 0.1),
        legend.title = element_blank())


ggsave("figures/endpoint_post_prior.pdf", width=15, height=12)


#  ranges_stat <- ranges_stat[-which(ranges_stat$median_start<0.001),]
# # 
# # ranges_stat_wm<- ranges_full[which(ranges_full$range=="Pygoscelis_grandis"),]
# # 
# # 
# p<-ggplot(ranges_stat) + #geom_point(aes(x=median_start, y=range), size=0.6) + geom_point(aes(x=median_end, y=range), size=0.6)+
#   # geom_segment(aes(x = median_start, xend=median_end, y = range, yend=range)) +
#   geom_density_ridges(data=ranges_full,aes(x = distr, y = as.factor(range), fill = endpoint),
#                        scale = 2, rel_min_height = 0.01, color = "white", alpha = .7, bandwidth = 0.5)+
#   scale_x_reverse() + theme_minimal()
# 
# p+geom_density_ridges(data=ranges_prior_full,aes(x = distr, y = as.factor(range)),
#                       scale = 2, rel_min_height = 0.01, color = "white", alpha = .7, bandwidth = 0.5)
# 
# p<-ggplot(ranges_stat) + #geom_point(aes(x=median_start, y=range), size=0.6) + geom_point(aes(x=median_end, y=range), size=0.6)+
#   # geom_segment(aes(x = median_start, xend=median_end, y = range, yend=range)) +
#   geom_density_ridges(data=ranges_full,aes(x = distr, y = as.factor(range), fill = endpoint),
#                       color = "white", alpha = .7, bandwidth = 0.5)+
#   scale_x_reverse() + theme_minimal()
# 
# p+geom_density_ridges(data=ranges_prior_full,aes(x = distr, y = as.factor(range)),
#                       color = "white", alpha = .7, bandwidth = 0.5) + facet_grid(range ~ endpoint)
# 
# ranges_full$mode<-"posterior"
# ranges_prior_full$mode<-"prior"
# p<-ggplot(rbind(ranges_full,ranges_prior_full)) + 
#   geom_density(aes(distr, fill=mode), alpha=0.5, color=NA)+ 
#   facet_wrap(range ~ endpoint, scales = "free", ncol=6) + theme_bw() +xlab("") +
#   theme(axis.text.x = element_text(angle=-60, vjust=0.5))
#  
# idx<- which(ranges_full$distr==0 & ranges_full$range %in% ranges_single)
# 
# 
# ggplot(ranges_stat_wm, aes(x = distr, y = range)) +
#   stat_density_ridges(geom = "density_ridges_gradient", quantile_lines = TRUE) +
#   scale_y_discrete(expand = expand_scale(mult = c(0.01, .7)))+
#   theme_classic() + ylab("") +guides(y="none")
# # 
# dt <- data.frame(x=c(1:200),y=rnorm(200))
# dens <- density(dt$y)
# df <- data.frame(x=dens$x, y=dens$y)
# probs <- c(0, 0.25, 0.5, 0.75, 1)
# quantiles <- quantile(dt$y, prob=probs)
# df$quant <- factor(findInterval(df$x,quantiles))
# ggplot(df, aes(x,y)) + geom_line() + geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) + scale_x_continuous(breaks=quantiles)
