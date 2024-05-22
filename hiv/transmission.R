library(coda)
wd <- "~/Documents/Source/beast2.7/sRanges-transmission-material/hiv/"
setwd(wd)

# this tool is from https://github.com/jugne/TnT
transmission_analyser_path <- "~/Documents/Source/TnT/out/artifacts/TransmissionAnalyser_jar/TransmissionAnalyser.jar"
system(paste("java -jar",transmission_analyser_path,"-tree inference/hiv_relaxedStart_env.45512.trees",
             "-out inference/hiv_relaxedStart_env.transmission.log"))

log <-  read.table("inference/hiv_relaxedStart_env.transmission.log", header=TRUE, sep="\t")
burnIn <- 0.1
log <- tail(log, -round(nrow(log)*burnIn)) # take out 10% burnin

system(paste("java -jar",transmission_analyser_path,"-tree inference/hiv_env.45512.trees",
             "-out inference/hiv_env.transmission.log"))

log2 <-  read.table("inference/hiv_env.transmission.log", header=TRUE, sep="\t")
log2 <- tail(log2, -round(nrow(log2)*burnIn)) # take out 10% burnin

data1<- round(sort(colMeans(log)[colMeans(log)!=0], decreasing = T), 3)
d <- data.frame(relaxed=data1)
rownames(d)<-names(data1)
write.csv(d, "figures/relaxed_direct_transmission_prob.csv")

data2<- round(sort(colMeans(log2)[colMeans(log2)!=0], decreasing = T), 3)
d <- data.frame(sampled=data2)
rownames(d)<-names(data2)
write.csv(d, "figures/direct_transmission_prob.csv")


columnsWithValue <- colnames(log)[colSums(log != "0.0,0.0") != 0]
# we don't need "Sample" column
columnsWithValue <- columnsWithValue[2:length(columnsWithValue)]
log <- log[, columnsWithValue]
columnsWithValue <- colnames(log)[colSums(log != "0.0,-1.0") != 0]
log <- log[, columnsWithValue]
n <- length(columnsWithValue)
d <- data.frame(from=character(n), to=character(n), prob=numeric(n), prob_direct=numeric(n), length_median=numeric(n),
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
  
  total <- length(length)
  idx <- which(length != 0)
  length <- length[idx] # no length means no speciation between the two species recorded.
  nodes <- nodes[idx] # 0 intermediate nodes are legit
  d$prob[i] <- length(length)/total
  direct_ids <- which(nodes==0)
  d$prob_direct[i] <- length(direct_ids)/total
  
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
  
  dd<-d[which(d$prob_direct!=0),]
  dd_ <- dd[, c(1,2,4)]

write.csv(dd_, "figures/direct_transmission_prob.csv")
  
}

