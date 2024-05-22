rm(list = ls())

library(seqinr)
library(ape)
setwd("~/Documents/Source/beast2.7/sRanges-transmission-material/hiv/data/")

pol<- read.alignment("set7_pol_Vrancken.fasta", format = "fasta")
pol_matrix <- as.matrix(pol)

consensus_pol<-matrix("",24,1449)
# consensus_pol<-rep(NA, length(unique(substr(pol$nam,1,6))))
i<-1
for (name in unique(substr(pol$nam,1,6))){
  consensus_pol[i,] <- con(pol_matrix[which(substr(pol$nam, 1,6)==name),], method="majority")
  i<-i+1
}

for (j in 1:24){
  label <- unique(substr(pol$nam,1,6))[j]
  if (substr(label,6,6)=="c"){
    label<- substr(label,1,4)
  }
  write.fasta(consensus_pol[j, ], label, "pol_consensus.fasta", open = "a", nbchar = 60, as.string = FALSE)
}


env<-read.alignment("set7_env_Vrancken.fasta", format = "fasta")
env_matrix <- as.matrix(env)
consensus_env<-matrix("",25,1557)
i<-1
for (name in unique(substr(env$nam,1,5))){
  if (name != "F05cl"){
    consensus_env[i,] <- con(env_matrix[which(substr(env$nam, 1,5)==name),], method="majority")
  }else {
    consensus_env[i,] <- env_matrix[which(substr(env$nam, 1,5)==name),]
  }
  
  label <- unique(substr(env$nam,1,5))[i]
  if (substr(label,5,5)=="l"){
    label<- substr(label,1,3)
  }
  write.fasta(consensus_env[i, ], label, "env_consensus.fasta", open = "a", nbchar = 60, as.string = FALSE)
  
  
  i<-i+1
}











