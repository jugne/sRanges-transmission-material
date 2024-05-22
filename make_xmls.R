## Read Tb fasta files and change their names. 
## Also add empty sequences for samples marking the start of the range.

library(seqinr)
library(stringr)
library(lubridate)
wd <- "/Users/jugne/Documents/Source/beast2.7/sRanges-transmission-material/"
setwd(wd)
data_path <- paste0(wd, "data/bulkseq_allserial/Georgia_L2_allserialbulk_filtered_newid.fasta")
outfile_name <- 'inference/bulkseq_allserial/no_removal_all_change/relaxed_clock/Georgia_L2_allserial.xml'
template_name <- "inference/bulkseq_allserial/no_removal_all_change/relaxed_clock/Georgia_L2_allserial_template.xml"
# data_path <- paste0(wd, "data/bulkseq_withinhost/Georgia_L2_M3bulk_newid.fasta")
# outfile_name <- 'inference/bulkseq_withinhostGeorgia_L2_bulk.xml'
# template_name <- "inference/bulkseq_withinhost/Georgia_L2_bulk_template.xml"
# data_path <- paste0(wd, "data/bulkseq_allserial_nonserial/Georgia_L2_allserialbulk_100nonserial_newid.fasta")
# outfile_name <- 'inference/bulkseq_allserial_nonserial/no_removal_all_change/relaxed_clock/Georgia_L2_allserial_nonserial.xml'
# template_name <- "inference/bulkseq_allserial_nonserial/no_removal_all_change/relaxed_clock/Georgia_L2_allserial_nonserial_template.xml"

f <- read.fasta(data_path)

f_hosts <- unique(unlist(lapply(names(f), function(x) {
  split_items <- unlist(strsplit(x, "/"))
  split_items[3]
})))

grouped_strings <- vector("list", length(f_hosts))
i = 1
for (h in f_hosts){
  matching_lists <- names(f)[grep(paste0("/",h, "$"), names(f))]
  grouped_strings_ <- unlist(lapply(matching_lists, paste, collapse = ", "))
  grouped_strings[[i]] <- grouped_strings_
  i=i+1
}

taxon_names <- vector("list", length(f_hosts))
taxon_dates <- vector("list", length(f_hosts))
taxon_dates_dec <- vector ("list", length(f_hosts))
taxon_seqs <- vector ("list", length(f_hosts))

i = 1
for (h in grouped_strings){
  l <- length(h)
  
    h_dates_tmp <- unlist(lapply(h, function(x) {
      split_items <- unlist(strsplit(x, "/"))
      split_items[2]
    }))
    h_sorted <- h[order(h_dates_tmp)]
 
  taxon_names[[i]] <- append(taxon_names[[i]],
                             paste0(unlist(str_split(h_sorted[1], "/"))[3], "_first"))
  taxon_dates[[i]] <- append(taxon_dates[[i]],
                             unlist(str_split(h_sorted[1], "/"))[2])
  taxon_dates_dec[[i]] <- append(taxon_dates_dec[[i]], decimal_date(as.Date(unlist(str_split(h_sorted[1], "/"))[2])))
  taxon_seqs[[i]] <- append(taxon_seqs[[i]], str_c(getSequence(f[h_sorted[1]])[[1]], collapse=""))
  if (l>1){
    if (l>2){
    for (j in 2:(length(h)-1)){
      taxon_names[[i]] <- append(taxon_names[[i]], 
                                 paste0(unlist(str_split(h_sorted[j], "/"))[3], "_", j))
      taxon_dates[[i]] <- append(taxon_dates[[i]], 
                                 unlist(str_split(h_sorted[j], "/"))[2])
      taxon_dates_dec[[i]] <- append(taxon_dates_dec[[i]], decimal_date(as.Date(unlist(str_split(h_sorted[j], "/"))[2])))
      taxon_seqs[[i]] <- append(taxon_seqs[[i]], str_c(getSequence(f[h_sorted[j]])[[1]], collapse=""))
    }}
    j<- length(h)
    taxon_names[[i]] <- append(taxon_names[[i]],
                               paste0(unlist(str_split(h_sorted[j], "/"))[3], "_last"))
    taxon_dates[[i]] <- append(taxon_dates[[i]], 
                               unlist(str_split(h_sorted[j], "/"))[2])
    taxon_dates_dec[[i]] <- append(taxon_dates_dec[[i]], decimal_date(as.Date(unlist(str_split(h_sorted[j], "/"))[2])))
    taxon_seqs[[i]] <- append(taxon_seqs[[i]], str_c(getSequence(f[h_sorted[j]])[[1]], collapse=""))
  }
  # } else {
  #   taxon_names[[i]] <- append(taxon_names[[i]],
  #                              paste0(unlist(str_split(h[1], "/"))[3], "_last"))
  #   taxon_dates[[i]] <- append(taxon_dates[[i]], 
  #                              unlist(str_split(h[1], "/"))[2])
  #   taxon_dates_dec[[i]] <- append(taxon_dates_dec[[i]], decimal_date(as.Date(unlist(str_split(h[1], "/"))[2])))
  # }
  
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
l_seq<-length(getSequence(f[1])[[1]])
for (i in 1:length(taxon_names)){
  tax<- taxon_names[[i]]
  range <- paste0('<stratigraphicRange id="r',i,'" spec="StratigraphicRange" ')
  range_refs <- c(range_refs, paste0('<stratigraphicRange idref="r',i,'"/>'))
  
  for (j in 1:length(tax)){
    if(j==1){
      range <- paste0(range, 'firstOccurrence="@',tax[j],'" ')
      firsts_tax <- c(firsts_tax, paste0('<taxon idref="',tax[j],'"/>'))
      if (length(tax)>1){
        taxSet <- paste0('<patientTaxonSets spec="TaxonSet" id="taxset_',i,'"> <taxon idref="',tax[j],'"/> ')
      }

    }
    if(j==length(tax)){
      range <- paste0(range, 'lastOccurrence="@',tax[j],'"/>')
      if (length(tax)>1){
        taxSet <- paste0(taxSet, '<taxon idref="',tax[j],'"/> </patientTaxonSets>')
      }
      

    } else if (length(tax)>1){
      if (j==2){
        range <- paste0(range, 'occurenceTaxonSet="@taxset_',i,'" ')
      }
      if (j!=1){
        taxSet <- paste0(taxSet, '<taxon idref="',tax[j],'"/> ')
      }
    }
    
    seq <- taxon_seqs[[i]][j]
    sequence_strings<- c(sequence_strings, paste0('<sequence id="seq_',i,'_',j,'" spec="Sequence" taxon="',tax[j],'" totalcount="4" value="',seq,'"/>'))
    taxon_strings <- c(taxon_strings, paste0('<taxon id="',tax[j],'" spec="Taxon"/>'))
    taxon_dates_xml <- c(taxon_dates_xml, paste0(tax[j], '=', taxon_dates[[i]][j]))
  }
  ranges <- c(ranges, range)
  if (length(tax)>1){
    range_taxSet <- c(range_taxSet, taxSet)
  }
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

inf <- readLines(paste0(wd, template_name))

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

writeLines(inf, con=outfile_name)


# # Given data
# data_str <- paste0("PT0000185_first=2014-04-02,PT0000185_last=2015-05-04,PT0000310_first=2016-05-25,PT0000310_last=2016-09-08,PT0000368_first=2012-05-08,PT0000368_last=2014-01-13,PT0000503_first=2014-01-30,PT0000503_2=2014-11-12,PT0000503_last=2016-06-08,PT0000920_first=2008-10-28,PT0000920_2=2008-12-24,PT0000920_3=2009-05-07,PT0000920_last=2009-06-10,PT0001013_first=2013-02-15,PT0001013_last=2014-02-27,PT0001140_first=2015-12-11,PT0001140_last=2016-04-15,PT0001227_first=2008-07-28,PT0001227_2=2009-02-25,PT0001227_3=2009-05-07,PT0001227_last=2009-05-18,PT0001228_first=2010-07-08,PT0001228_2=2010-08-12,PT0001228_last=2011-10-12,PT0001268_first=2011-01-14,PT0001268_2=2012-01-25,PT0001268_3=2012-05-01,PT0001268_last=2012-05-25,PT0001346_first=2009-03-19,PT0001346_last=2011-10-06,PT0001399_first=2008-09-09,PT0001399_2=2008-12-17,PT0001399_last=2009-05-06,PT0001515_first=2016-01-20,PT0001515_last=2016-03-30,PT0001619_first=2011-09-30,PT0001619_2=2012-07-10,PT0001619_3=2012-10-16,PT0001619_last=2013-03-19,PT0001705_first=2008-06-03,PT0001705_2=2008-08-29,PT0001705_last=2008-12-29,PT0001706_first=2015-05-05,PT0001706_last=2015-08-21,PT0001728_first=2013-05-16,PT0001728_last=2016-06-06,PT0001815_first=2008-10-03,PT0001815_2=2008-10-06,PT0001815_3=2008-12-24,PT0001815_last=2009-01-27,PT0001819_first=2015-12-30,PT0001819_last=2016-04-07,PT0002014_first=2010-11-09,PT0002014_2=2011-02-08,PT0002014_last=2012-02-09,PT0002027_first=2010-07-15,PT0002027_2=2010-08-24,PT0002027_3=2010-09-07,PT0002027_last=2011-01-17,PT0002095_first=2011-11-02,PT0002095_2=2011-11-22,PT0002095_3=2012-06-08,PT0002095_last=2012-06-19,PT0002122_first=2011-01-18,PT0002122_2=2012-06-06,PT0002122_3=2012-09-20,PT0002122_last=2014-06-20,PT0002140_first=2009-10-22,PT0002140_2=2010-03-18,PT0002140_last=2010-07-13,PT0002159_first=2010-08-12,PT0002159_last=2010-10-15,PT0002192_first=2009-03-16,PT0002192_2=2009-04-27,PT0002192_3=2009-05-22,PT0002192_last=2009-05-25,PT0002445_first=2008-08-04,PT0002445_2=2009-01-09,PT0002445_last=2009-05-08,PT0002622_first=2009-04-24,PT0002622_last=2009-05-22,PT0002720_first=2015-11-30,PT0002720_last=2016-03-10,PT0002756_first=2010-05-28,PT0002756_last=2010-12-30,PT0002791_first=2016-07-05,PT0002791_2=2016-10-04,PT0002791_last=2017-03-09,PT0002890_first=2016-10-24,PT0002890_last=2016-12-06,PT0002893_first=2019-06-02,PT0002893_last=2019-07-12,PT0002917_first=2017-01-20,PT0002917_last=2018-01-25,PT0002977_first=2013-06-26,PT0002977_last=2015-11-12,PT0003027_first=2009-07-13,PT0003027_2=2010-08-10,PT0003027_3=2010-11-05,PT0003027_last=2011-02-25,PT0003199_first=2010-10-25,PT0003199_2=2012-07-02,PT0003199_last=2012-12-25,PT0003420_first=2013-03-22,PT0003420_2=2013-04-25,PT0003420_last=2013-05-24,PT0003617_first=2015-08-18,PT0003617_last=2015-12-15,PT0003898_first=2009-05-13,PT0003898_2=2009-05-27,PT0003898_3=2009-06-04,PT0003898_4=2011-09-16,PT0003898_last=2012-10-30,PT0004070_first=2016-10-28,PT0004070_last=2017-02-08,PT0004163_first=2016-11-03,PT0004163_2=2017-02-27,PT0004163_3=2017-05-31,PT0004163_last=2018-01-18,PT0004348_first=2016-02-22,PT0004348_last=2016-06-28,PT0004375_first=2012-01-23,PT0004375_last=2012-08-03,PT0004458_first=2018-07-11,PT0004458_last=2019-06-24,PT0004684_first=2012-04-18,PT0004684_2=2012-07-24,PT0004684_3=2012-09-05,PT0004684_last=2012-12-04,PT0004798_first=2009-03-13,PT0004798_last=2009-05-15,PT0004830_first=2009-08-25,PT0004830_2=2009-10-23,PT0004830_3=2009-11-13,PT0004830_last=2010-01-18,PT0004947_first=2009-02-27,PT0004947_2=2009-06-01,PT0004947_last=2009-06-15,PT0004964_first=2010-01-25,PT0004964_last=2013-02-19,PT0005001_first=2017-02-22,PT0005001_last=2017-04-19,PT0005062_first=2012-09-10,PT0005062_2=2012-11-08,PT0005062_3=2012-12-11,PT0005062_4=2013-01-11,PT0005062_5=2013-02-07,PT0005062_6=2013-06-13,PT0005062_last=2014-06-11,PT0005314_first=2008-11-21,PT0005314_2=2009-01-21,PT0005314_last=2009-03-30,PT0006080_first=2009-11-03,PT0006080_2=2011-09-30,PT0006080_3=2012-03-28,PT0006080_last=2012-05-10,PT0006547_first=2010-07-09,PT0006547_last=2011-01-27,PT0006597_first=2008-09-10,PT0006597_last=2008-10-22,PT0006682_first=2011-08-17,PT0006682_last=2012-05-31,PT0006683_first=2009-08-19,PT0006683_2=2009-10-16,PT0006683_3=2009-11-20,PT0006683_4=2009-12-21,PT0006683_last=2010-01-27,PT0006793_first=2016-11-17,PT0006793_2=2017-02-28,PT0006793_last=2017-07-26,PT0006954_first=2010-11-22,PT0006954_2=2010-11-26,PT0006954_3=2010-12-24,PT0006954_last=2011-02-24,PT0007055_first=2008-11-12,PT0007055_last=2009-05-07,PT0007246_first=2016-05-13,PT0007246_last=2016-08-15,PT0007348_first=2009-02-19,PT0007348_2=2009-04-30,PT0007348_3=2009-08-19,PT0007348_4=2009-12-04,PT0007348_last=2010-04-15,PT0007876_first=2012-07-20,PT0007876_last=2013-03-15,PT0008072_first=2010-10-05,PT0008072_2=2010-10-19,PT0008072_3=2011-02-18,PT0008072_last=2011-08-24,PT0008440_first=2015-03-02,PT0008440_2=2017-01-23,PT0008440_last=2017-05-25,PT0008576_first=2012-07-03,PT0008576_last=2012-12-24,PT0008933_first=2010-03-12,PT0008933_last=2010-06-25,PT0009125_first=2013-12-26,PT0009125_last=2017-02-15,PT0009245_first=2011-08-26,PT0009245_2=2011-09-28,PT0009245_3=2012-02-17,PT0009245_4=2012-05-25,PT0009245_last=2012-09-03,PT0009255_first=2020-03-16,PT0009255_last=2021-02-11,PT0009310_first=2011-02-02,PT0009310_last=2012-06-21,PT0009515_first=2014-04-04,PT0009515_last=2016-03-14,PT0009652_first=2008-06-19,PT0009652_2=2009-02-23,PT0009652_last=2009-04-30,PT0010043_first=2009-08-18,PT0010043_2=2010-01-05,PT0010043_last=2010-08-17,PT0010156_first=2009-03-24,PT0010156_2=2009-12-07,PT0010156_3=2011-12-02,PT0010156_last=2012-01-27,PT0010246_first=2017-04-06,PT0010246_last=2017-06-26,PT0010658_first=2008-10-23,PT0010658_2=2008-12-12,PT0010658_3=2009-02-18,PT0010658_last=2009-07-09,PT0010703_first=2014-04-14,PT0010703_2=2016-11-24,PT0010703_last=2017-01-31,PT0010709_first=2017-02-13,PT0010709_last=2017-05-22,PT0010750_first=2010-06-18,PT0010750_2=2010-10-13,PT0010750_3=2010-11-04,PT0010750_4=2010-11-11,PT0010750_5=2011-09-30,PT0010750_last=2012-06-05,PT0010850_first=2012-03-22,PT0010850_2=2012-06-17,PT0010850_3=2012-09-24,PT0010850_last=2013-05-24,PT0010851_first=2012-11-16,PT0010851_last=2015-05-28,PT0011184_first=2009-08-26,PT0011184_2=2009-09-23,PT0011184_3=2009-12-01,PT0011184_4=2011-01-14,PT0011184_5=2011-02-17,PT0011184_last=2012-09-17,PT0011277_first=2009-07-03,PT0011277_2=2009-09-04,PT0011277_3=2010-03-24,PT0011277_4=2010-10-08,PT0011277_last=2010-11-16,PT0011457_first=2008-10-07,PT0011457_2=2009-03-16,PT0011457_last=2009-03-18,PT0011723_first=2017-01-20,PT0011723_last=2018-10-03,PT0011769_first=2010-08-20,PT0011769_2=2010-11-29,PT0011769_3=2011-08-23,PT0011769_last=2012-12-04,PT0011823_first=2016-12-06,PT0011823_2=2017-03-16,PT0011823_last=2017-05-15,PT0011865_first=2013-01-11,PT0011865_last=2016-03-28,PT0012325_first=2009-09-14,PT0012325_last=2012-05-18,PT0012611_first=2015-06-09,PT0012611_last=2015-12-22,PT0012652_first=2008-11-18,PT0012652_2=2009-02-18,PT0012652_3=2009-03-18,PT0012652_4=2009-03-20,PT0012652_last=2009-06-11,PT0012685_first=2009-09-04,PT0012685_2=2011-09-01,PT0012685_last=2014-05-01,PT0013537_first=2008-07-10,PT0013537_2=2008-10-15,PT0013537_3=2008-11-10,PT0013537_4=2008-12-25,PT0013537_last=2009-03-27,PT0013620_first=2016-09-15,PT0013620_last=2016-12-08,PT0013682_first=2011-11-01,PT0013682_2=2012-04-20,PT0013682_last=2012-05-21,PT0013793_first=2016-06-30,PT0013793_last=2018-03-01,PT0013996_first=2010-10-07,PT0013996_2=2010-11-04,PT0013996_last=2011-03-10,PT0013998_first=2010-12-23,PT0013998_2=2011-02-15,PT0013998_3=2013-02-11,PT0013998_4=2013-04-26,PT0013998_last=2013-05-24,PT0014046_first=2008-11-12,PT0014046_2=2009-03-16,PT0014046_last=2009-06-10,PT0014140_first=2010-03-09,PT0014140_2=2011-08-25,PT0014140_last=2012-05-01,PT0014214_first=2008-11-05,PT0014214_2=2009-02-04,PT0014214_3=2009-04-30,PT0014214_last=2009-05-29,PT0014332_first=2008-05-14,PT0014332_2=2008-09-29,PT0014332_3=2008-11-12,PT0014332_4=2009-01-22,PT0014332_last=2009-03-16,PT0014349_first=2008-11-06,PT0014349_last=2009-02-06,PT0015022_first=2021-04-13,PT0015022_2=2021-05-11,PT0015022_last=2021-09-15,PT0015022b_first=2021-10-20,PT0015022b_last=2021-11-19,PT0015231_first=2009-07-17,PT0015231_last=2012-09-19,PT0015235_first=2017-01-27,PT0015235_last=2017-04-28,PT0015273_first=2012-02-29,PT0015273_2=2012-03-02,PT0015273_3=2012-09-14,PT0015273_last=2013-01-21,PT0015309_first=2017-01-12,PT0015309_2=2017-05-29,PT0015309_last=2017-06-21,PT0015335_first=2010-05-13,PT0015335_last=2010-09-03,PT0015386_first=2013-06-07,PT0015386_last=2014-04-01,PT0015912_first=2008-10-31,PT0015912_2=2009-01-20,PT0015912_last=2009-03-23,PT0016117_first=2008-06-20,PT0016117_2=2008-07-17,PT0016117_3=2008-09-19,PT0016117_4=2008-11-20,PT0016117_last=2009-02-10,PT0016700_first=2021-08-04,PT0016700_2=2021-08-06,PT0016700_last=2021-09-08,PT0016719_first=2022-03-10,PT0016719_2=2022-03-22,PT0016719_3=2022-03-29,PT0016719_4=2022-04-04,PT0016719_5=2022-04-15,PT0016719_last=2022-06-15,PT0016722_first=2021-11-03,PT0016722_2=2021-11-04,PT0016722_3=2021-12-03,PT0016722_last=2022-02-07,PT0017314_first=2022-03-25,PT0017314_2=2022-04-01,PT0017314_last=2022-04-26,PT0017345_first=2021-03-17,PT0017345_last=2021-04-20,PT0017365_first=2021-11-11,PT0017365_2=2021-11-12,PT0017365_last=2021-12-14,PT0017427_first=2021-02-25,PT0017427_last=2021-03-26,PT0017539_first=2021-09-17,PT0017539_2=2021-09-20,PT0017539_last=2021-10-20,PT0018280_first=2021-06-09,PT0018280_last=2021-07-09")
# 
# data_str <- "PT0000185_first=2014-04-02,PT0000185_last=2015-05-04,PT0000310_first=2016-05-25,PT0000310_last=2016-09-08,PT0000368_first=2012-05-08,PT0000368_last=2014-01-13,PT0000503_first=2014-01-30,PT0000503_2=2014-11-12,PT0000503_last=2016-06-08,PT0000920_first=2008-10-28,PT0000920_2=2008-12-24,PT0000920_3=2009-05-07,PT0000920_last=2009-06-10,PT0001013_first=2013-02-15,PT0001013_last=2014-02-27,PT0001140_first=2015-12-11,PT0001140_last=2016-04-15,PT0001227_first=2008-07-28,PT0001227_2=2009-02-25,PT0001227_3=2009-05-07,PT0001227_last=2009-05-18,PT0001228_first=2010-07-08,PT0001228_2=2010-08-12,PT0001228_last=2011-10-12,PT0001268_first=2011-01-14,PT0001268_2=2012-01-25,PT0001268_3=2012-05-01,PT0001268_last=2012-05-25,PT0001346_first=2009-03-19,PT0001346_last=2011-10-06,PT0001399_first=2008-09-09,PT0001399_2=2008-12-17,PT0001399_last=2009-05-06,PT0001515_first=2016-01-20,PT0001515_last=2016-03-30,PT0001619_first=2011-09-30,PT0001619_2=2012-07-10,PT0001619_3=2012-10-16,PT0001619_last=2013-03-19,PT0001705_first=2008-06-03,PT0001705_2=2008-08-29,PT0001705_last=2008-12-29,PT0001706_first=2015-05-05,PT0001706_last=2015-08-21,PT0001728_first=2013-05-16,PT0001728_last=2016-06-06,PT0001815_first=2008-10-03,PT0001815_2=2008-10-06,PT0001815_3=2008-12-24,PT0001815_last=2009-01-27,PT0001819_first=2015-12-30,PT0001819_last=2016-04-07,PT0002014_first=2010-11-09,PT0002014_2=2011-02-08,PT0002014_last=2012-02-09,PT0002027_first=2010-07-15,PT0002027_2=2010-08-24,PT0002027_3=2010-09-07,PT0002027_last=2011-01-17,PT0002095_first=2011-11-02,PT0002095_2=2011-11-22"
# data_str<-"PT0013537_2=2008-10-15,PT0013537_3=2008-11-10,PT0013537_4=2008-12-25,PT0013537_last=2009-03-27,PT0013620_first=2016-09-15,PT0013620_last=2016-12-08,PT0013682_first=2011-11-01,PT0013682_2=2012-04-20,PT0013682_last=2012-05-21,PT0013793_first=2016-06-30,PT0013793_last=2018-03-01,PT0013996_first=2010-10-07,PT0013996_2=2010-11-04,PT0013996_last=2011-03-10,PT0013998_first=2010-12-23,PT0013998_2=2011-02-15,PT0013998_3=2013-02-11,PT0013998_4=2013-04-26,PT0013998_last=2013-05-24,PT0014046_first=2008-11-12,PT0014046_2=2009-03-16,PT0014046_last=2009-06-10,PT0014140_first=2010-03-09,PT0014140_2=2011-08-25,PT0014140_last=2012-05-01,PT0014214_first=2008-11-05,PT0014214_2=2009-02-04,PT0014214_3=2009-04-30,PT0014214_last=2009-05-29,PT0014332_first=2008-05-14,PT0014332_2=2008-09-29,PT0014332_3=2008-11-12,PT0014332_4=2009-01-22,PT0014332_last=2009-03-16,PT0014349_first=2008-11-06,PT0014349_last=2009-02-06,PT0015022_first=2021-04-13,PT0015022_2=2021-05-11,PT0015022_last=2021-09-15,PT0015022b_first=2021-10-20,PT0015022b_last=2021-11-19,PT0015231_first=2009-07-17,PT0015231_last=2012-09-19,PT0015235_first=2017-01-27,PT0015235_last=2017-04-28,PT0015273_first=2012-02-29,PT0015273_2=2012-03-02,PT0015273_3=2012-09-14,PT0015273_last=2013-01-21,PT0015309_first=2017-01-12,PT0015309_2=2017-05-29,PT0015309_last=2017-06-21,PT0015335_first=2010-05-13,PT0015335_last=2010-09-03,PT0015386_first=2013-06-07,PT0015386_last=2014-04-01,PT0015912_first=2008-10-31,PT0015912_2=2009-01-20,PT0015912_last=2009-03-23,PT0016117_first=2008-06-20,PT0016117_2=2008-07-17,PT0016117_3=2008-09-19,PT0016117_4=2008-11-20,PT0016117_last=2009-02-10,PT0016700_first=2021-08-04,PT0016700_2=2021-08-06,PT0016700_last=2021-09-08,PT0016719_first=2022-03-10,PT0016719_2=2022-03-22,PT0016719_3=2022-03-29,PT0016719_4=2022-04-04,PT0016719_5=2022-04-15,PT0016719_last=2022-06-15,PT0016722_first=2021-11-03,PT0016722_2=2021-11-04,PT0016722_3=2021-12-03,PT0016722_last=2022-02-07,PT0017314_first=2022-03-25,PT0017314_2=2022-04-01,PT0017314_last=2022-04-26,PT0017345_first=2021-03-17,PT0017345_last=2021-04-20,PT0017365_first=2021-11-11,PT0017365_2=2021-11-12,PT0017365_last=2021-12-14,PT0017427_first=2021-02-25,PT0017427_last=2021-03-26,PT0017539_first=2021-09-17,PT0017539_2=2021-09-20,PT0017539_last=2021-10-20,PT0018280_first=2021-06-09,PT0018280_last=2021-07-09"
# # Convert data string into a data frame
# data_list <- strsplit(data_str, ",")[[1]]
# data_frame <- as.data.frame(do.call(rbind, strsplit(data_list, "=")), stringsAsFactors = FALSE)
# colnames(data_frame) <- c("ID", "Date")
# 
# # Remove rows where the ID doesn't end with "_first"
# filtered_data <- subset(data_frame, grepl("_first$", ID))
# 
# print(filtered_data)
# 
# result_str <- paste(filtered_data$ID, filtered_data$Date, sep="=", collapse=",")
# cat(result_str)

