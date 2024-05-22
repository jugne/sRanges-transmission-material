## Read Tb fasta files and change their names. 
## Also add empty sequences for samples marking the start of the range.

library(seqinr)
library(lubridate)
wd <- "/Users/jugne/Documents/Source/beast2.7/sRanges-transmission-material/"
setwd(wd)
data_path <- paste0(wd, "data/bulkseq_withinhost/Georgia_L2_M3bulk_newid.fasta")

f <- read.fasta(data_path)

f_hosts <- unique(unlist(lapply(names(f), function(x) {
  split_items <- unlist(strsplit(x, "/"))
  split_items[3]
})))

grouped_strings <- vector("list", length(f_hosts))
i = 1
for (h in f_hosts){
  matching_lists <- names(f)[grep(paste0("/",h), names(f))]
  grouped_strings_ <- unlist(lapply(matching_lists, paste, collapse = ", "))
  grouped_strings[[i]] <- grouped_strings_
  i=i+1
}

taxon_names <- vector("list", length(f_hosts))
taxon_dates <- vector("list", length(f_hosts))
taxon_dates_dec <- vector ("list", length(f_hosts))
taxon_seqs <- vector ("list", length(f_hosts))
l_seq<-length(getSequence(f[1])[[1]])

i = 1
for (h in grouped_strings){
  l <- length(h)
  taxon_names[[i]] <- append(taxon_names[[i]],
                             paste0(unlist(str_split(h[1], "/"))[3], "_first"))
  year <- strtoi(str_split(unlist(str_split(h[1], "/"))[2], "-")[[1]][1])-1
  date <- paste0(year,substr(unlist(str_split(h[1], "/"))[2], 5, 13))
  taxon_dates[[i]] <- append(taxon_dates[[i]], 
                             date)
  taxon_dates_dec[[i]] <- append(taxon_dates_dec[[i]], decimal_date(as.Date(date)))
  taxon_seqs[[i]] <- append(taxon_seqs[[i]], str_c(strrep("?", l_seq), collapse=""))
  if (l>1){
    h_dates_tmp <- unlist(lapply(h, function(x) {
      split_items <- unlist(strsplit(x, "/"))
      split_items[2]
    }))
    h_sorted <- h[order(h_dates_tmp)]
    for (j in 1:(length(h)-1)){
      taxon_names[[i]] <- append(taxon_names[[i]], 
                                 paste0(unlist(str_split(h_sorted[j], "/"))[3], "_", j))
      taxon_dates[[i]] <- append(taxon_dates[[i]], 
                                 unlist(str_split(h_sorted[j], "/"))[2])
      taxon_dates_dec[[i]] <- append(taxon_dates_dec[[i]], decimal_date(as.Date(unlist(str_split(h_sorted[j], "/"))[2])))
      taxon_seqs[[i]] <- append(taxon_seqs[[i]], str_c(getSequence(f[h_sorted[j]])[[1]], collapse=""))
    }
    j<- length(h)
    taxon_names[[i]] <- append(taxon_names[[i]],
                               paste0(unlist(str_split(h_sorted[j], "/"))[3], "_last"))
    taxon_dates[[i]] <- append(taxon_dates[[i]], 
                               unlist(str_split(h_sorted[j], "/"))[2])
    taxon_dates_dec[[i]] <- append(taxon_dates_dec[[i]], decimal_date(as.Date(unlist(str_split(h_sorted[j], "/"))[2])))
    taxon_seqs[[i]] <- append(taxon_seqs[[i]], str_c(getSequence(f[h_sorted[j]])[[1]], collapse=""))
  } else{
    taxon_names[[i]] <- append(taxon_names[[i]],
                               paste0(unlist(str_split(h[1], "/"))[3], "_last"))
    taxon_dates[[i]] <- append(taxon_dates[[i]], 
                               unlist(str_split(h[1], "/"))[2])
    taxon_dates_dec[[i]] <- append(taxon_dates_dec[[i]], decimal_date(as.Date(unlist(str_split(h[1], "/"))[2])))
    taxon_seqs[[i]] <- append(taxon_seqs[[i]], str_c(getSequence(f[h[1]])[[1]], collapse=""))
  }
 
  i=i+1
}

firsts <- c()
sampling_dates <- c()
max_sample_time<- max(unlist(taxon_dates_dec))
for (m in 1:length(taxon_dates_dec)){
  firsts <- c(firsts, max_sample_time-taxon_dates_dec[[m]][1])
  
    
}

sequence_strings <- c()
taxon_strings <- c()
taxon_dates_xml <- c()
ranges <- c()
range_refs<- c()
range_taxSet <- c()
firsts_tax <- c()
k<-1

for (i in 1:length(taxon_names)){
  tax<- taxon_names[[i]]
  range <- paste0('<stratigraphicRange id="r',i,'" spec="StratigraphicRange" ')
  range_refs <- c(range_refs, paste0('<stratigraphicRange idref="r',i,'"/>'))
  
  for (j in 1:length(tax)){
    if(j==1){
      
      range <- paste0(range, 'firstOccurrence="@',tax[j],'" ')
      taxSet <- paste0('<patientTaxonSets spec="TaxonSet" id="taxset_',i,'"> <taxon idref="',tax[j],'"/> ')
      firsts_tax <- c(firsts_tax, paste0('<taxon idref="',tax[j],'"/>'))
      sampling_dates <- c(sampling_dates, paste0('<samplingDates spec="sa.evolution.tree.SamplingDate" taxon="@',tax[j],'" lower="',firsts[i]-1,'" upper="',firsts[i],'"/>'))
    }else if(j==length(tax)){
      range <- paste0(range, 'lastOccurrence="@',tax[j],'"/>')
      taxSet <- paste0(taxSet, '<taxon idref="',tax[j],'"/> </patientTaxonSets>')
      
    } else {
      if (j==2){
        range <- paste0(range, 'occurenceTaxonSet="@taxset_',i,'" ')
      }
      taxSet <- paste0(taxSet, '<taxon idref="',tax[j],'"/> ')
    }
    
    seq <- taxon_seqs[[i]][j]
    sequence_strings<- c(sequence_strings, paste0('<sequence id="seq_',i,'_',j,'" spec="Sequence" taxon="',tax[j],'" totalcount="4" value="',seq,'"/>'))
    taxon_strings <- c(taxon_strings, paste0('<taxon id="',tax[j],'" spec="Taxon"/>'))
    taxon_dates_xml <- c(taxon_dates_xml, paste0(tax[j], '=', taxon_dates[[i]][j]))
  }
  ranges <- c(ranges, range)
  range_taxSet <- c(range_taxSet, taxSet)
}
i1 <- sort(unlist(taxon_dates), index.return=TRUE)$ix
seq_str<- str_c(sequence_strings[i1], collapse = '\n\t\t')
tax_str <- str_c(taxon_strings[i1], collapse = '\n\t\t')
date_str <- paste0('value="',str_c(taxon_dates_xml, collapse=','), '">')
ranges_str <- str_c(ranges, collapse="\n\t\t\t\t")
ranges_refs_str <- str_c(range_refs, collapse="\n\t\t")
taxSet_str <- str_c(range_taxSet, collapse='\n\t\t')
firsts_tax_str <- str_c(firsts_tax, collapse='\n\t\t\t')
sampling_dates_str <- str_c(sampling_dates, collapse="\n\t\t\t")

# Make xml from template

inf <- readLines(paste0(wd, "Georgia_L2_bulk_relaxedStart_template.xml"))

inf  <- gsub(pattern = "<insertTaxonSet/>",
               replace = tax_str, x=inf)
inf  <- gsub(pattern = "<insertAlignment/>",
             replace = seq_str, x=inf)
inf  <- gsub(pattern = "<insertTaxonSets/>",
             replace = taxSet_str, x = inf)
inf  <- gsub(pattern = "<insertRangesRefs/>",
             replace = ranges_refs_str, x = inf)
inf  <- gsub(pattern = "<insertStratigraphicRanges/>",
             replace = ranges_str, x = inf)
inf  <- gsub(pattern = "<insertDates/>",
             replace = date_str, x = inf)
inf  <- gsub(pattern = "<insertRandomWalkerTaxonSet/>",
             replace = firsts_tax_str, x = inf)
inf  <- gsub(pattern = "<insertSamplingDates/>",
             replace = sampling_dates_str, x = inf)

writeLines(inf, con='Georgia_L2_bulk_relaxedStart.xml')



