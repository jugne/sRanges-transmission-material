library(weights)
library(ggplot2)

wd <- "/Users/jugne/Documents/Source/beast2.7/sRanges-transmission-material/inference/bulkseq_allserial/no_removal_all_change/"
# Define the filename
filename <- paste0(wd, "relaxed_clock/new_sampling_prior/Georgia_L2_allserial.233457685745_rates.txt")

# Open the connection to the file
connection <- file(filename, open = "r")

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

id_range <- which(rates_range>1e-5)
id_noRange <- which(rates_noRange>1e-5)
wtd.hist(rates_noRange, breaks= 100000,weight=weight_noRange)
wtd.hist(rates_range, breaks= 50000,weight=weight_range)
wtd.hist(rates_noRange[id_noRange], breaks= 1000,weight=weight_noRange[id_noRange], xlim=c(0, 2e-3), density = 30,
         col = rgb(1, 0.8, 0.8, 0.5), add=T)
wtd.hist(rates_range[id_range], breaks= 1000,weight=weight_range[id_range], xlim=c(0, 2e-3), density = 30,
         col = rgb(0.6, 0.8, 1, 0.5))

plot(density(rates_range[id_range]))# Close the connection to the file
plot(density(rates_noRange[id_noRange]))


dd <- data.frame(range=c(rep("yes", length(rates_range)), rep("no", length(rates_noRange))), vals=c(rates_range, rates_noRange))
dd <- data.frame(range=c(rep("yes", length(id_range)), rep("no", length(id_noRange))), vals=c(rates_range[id_range], rates_noRange[id_noRange]))

id_range <- which(rates_range<1e-5)
id_noRange <- which(rates_noRange<1e-5)
wtd.hist(rates_noRange[id_noRange], breaks= 1000,weight=weight_noRange[id_noRange], xlim=c(0, 1e-5), density = 30,
         col = rgb(1, 0.8, 0.8, 0.5))
wtd.hist(rates_range[id_range], breaks= 1000,weight=weight_range[id_range], xlim=c(0, 1e-5), density = 30,
         col = rgb(0.6, 0.8, 1, 0.5), add=T)


wd <- "/Users/jugne/Documents/Source/beast2.7/sRanges-transmission-material/inference/bulkseq_allserial/no_removal_all_change/"
# Define the filename
filename <- paste0(wd, "relaxed_clock/new_sampling_prior/Georgia_L2_allserial.233457685745_rates_norm.txt")

# Open the connection to the file
connection <- file(filename, open = "r")

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

id_range <- which(rates_range>1e-5)
id_noRange <- which(rates_noRange>1e-5)
 dd <- data.frame(range=c(rep("range", length(id_range)), rep("non range", length(id_noRange))), vals=c(rates_range[id_range], rates_noRange[id_noRange]), weights=c(weight_range[id_range], weight_noRange[id_noRange]))
 dd_ <- data.frame(range=c(rep("range", length(rates_range)), rep("non range", length(rates_noRange))), vals=c(rates_range, rates_noRange), weights=c(weight_range, weight_noRange))
 p<- ggplot(dd, aes(range, vals)) +geom_density(aes(weight=weights))
 ggplot(dd_, aes(range, log10(vals), weight=weights)) + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75))
 ggplot(dd, aes(range, log10(vals), weight=weights)) + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75))


wtd.hist(log10(rates_noRange), breaks=500,weight=weight_noRange, density = 30,
         col = rgb(1, 0.8, 0.8, 0.5), add=T)
wtd.hist(log10(rates_range), breaks= 500,weight=weight_range, density = 30,
         col = rgb(0.6, 0.8, 1, 0.9))

wtd.hist(rates_noRange, breaks=50,weight=weight_noRange, ylim=c(0,45), density = 30,
         col = rgb(1, 0.8, 0.8, 0.9), add=T)
wtd.hist(rates_range, breaks= 100,weight=weight_range, ylim=c(0,45), density = 30,
         col = rgb(0.6, 0.8, 1, 0.9))

wtd.hist(rates_noRange, breaks=100,weight=weight_noRange, density = 30,
         col = rgb(1, 0.8, 0.8, 0.9))
wtd.hist(rates_range, breaks=100,weight=weight_range, density = 30,
         col = rgb(0.6, 0.8, 1, 0.9))

wtd.hist(rates_noRange[id_noRange], breaks= 100,weight=weight_noRange[id_noRange], xlim=c(0, 2e-3), density = 30,
         col = rgb(1, 0.8, 0.8, 0.5), add=T)
wtd.hist(rates_range[id_range], breaks= 100,weight=weight_range[id_range], density = 30,
         col = rgb(0.6, 0.8, 1, 0.5))

plot(density(log10(rates_range)))# Close the connection to the file
plot(density(log10(rates_noRange)))



id_range <- which(rates_range<1e-5)
id_noRange <- which(rates_noRange<1e-5)
wtd.hist(rates_noRange[id_noRange], breaks= 100,weight=weight_noRange[id_noRange], xlim=c(0, 1e-5), density = 15,
         col = rgb(1, 0.8, 0.8, 0.9), add=T)
wtd.hist(rates_range[id_range], breaks= 100,weight=weight_range[id_range], xlim=c(0, 1e-5), density = 15,
         col = rgb(0.6, 0.8, 1, 0.9))

