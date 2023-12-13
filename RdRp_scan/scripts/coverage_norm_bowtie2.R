#!/usr/bin/env Rscript --vanilla

######  setwd("/Volumes/@geva/Artem/NGS_mosquito_metagenomics_method_comparison/mapping_summary")

library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
cat(paste0("\nArguments passed to the script: ", paste(args, collapse = " "), "\n-----"))

dir<- args[1]
sample_id <- args[2]
ref_id <- args[3]
bowtie2_mapinfo_file_path <- args[4]
cov_file_path <- args[5]
vcf_file_path <- args[6]
summary_folder <- args[7]

cat(paste0("\nDirectory : ", dir))
cat(paste0("\nSample : ", sample_id))
cat(paste0("\nReference : ", ref_id))
cat(paste0("\nMapping info file path : ", bowtie2_mapinfo_file_path))
cat(paste0("\nCoverage file path : ", cov_file_path))
cat(paste0("\nSNV file path : ", vcf_file_path))

####### parse bowtie2 specifically for read numbers rather than through the percentage of aligned which is rounded

### These are paired reads, however, bowtie2 counts each pair as 1 read, and then splits those into mates.
### This number is multiplied by 2 to count the total paired mates.
pe_mates <- as.numeric(str_extract(read_lines(
  bowtie2_mapinfo_file_path, skip = 1, n_max = 1),
  pattern = "\\d+(?=\\s+\\()")) * 2


### This line provides the number of paired mates that did not map concordantly or discordantly
unal_pe_mates_concord_or_discord <- as.numeric(str_extract(
  read_lines(bowtie2_mapinfo_file_path, skip = 11, n_max = 1),
  pattern = "\\d+(?=\\s+\\()"))

### Now account for unpaired reads that were provided to bowtie2
unp_reads <- as.numeric(str_extract(read_lines(
  bowtie2_mapinfo_file_path, skip = 14, n_max = 1),
  pattern = "\\d+(?=\\s+\\()"))

### This is the number of unpaired reads that did not map at all.
unal_unp_reads <- as.numeric(str_extract(read_lines(
  bowtie2_mapinfo_file_path, skip = 15, n_max = 1),
  pattern = "\\d+(?=\\s+\\()"))

total_reads <- pe_mates + unp_reads
mapped_reads <- total_reads - unal_pe_mates_concord_or_discord - unal_unp_reads

cat(paste0("\n\n total_reads ", total_reads, "\n\n"))
cat(paste0("\n\n mapped_reads ", mapped_reads, "\n\n"))

####
cov_df <- read.table(cov_file_path, header=FALSE)
colnames(cov_df) = c("ref_id","position","depth")
norm_cov_df <- cov_df
norm_cov_df$depth_norm <- (norm_cov_df$depth/(as.numeric(total_reads)))*10^6

norm_cov_df$sample_id <- sample_id
cov_df$depth <- as.numeric(cov_df$depth)
cov_df$mapped_pos <- 0
if (sum(cov_df$depth) > 0){
  cov_df[which(!cov_df$depth == 0),]$mapped_pos <- 1
}

percentage_covered <- (sum(cov_df$mapped_pos)/max(as.numeric(cov_df$position)))*100
cat(paste0("\n percentage_covered ", percentage_covered))


####
average_cov <- sum(norm_cov_df$depth)/length(norm_cov_df$depth)
cat(paste0("\n average_cov  ", average_cov ))

norm_average_cov <- (average_cov/(as.numeric(total_reads)))*10^6
cat(paste0("\n norm_average_cov  ", norm_average_cov ))

####
summary_df <- tibble(sample_id, ref_id, total_reads, mapped_reads, average_cov, norm_average_cov, percentage_covered)
####
write.table(x = norm_cov_df, file = gzfile(paste0(summary_folder, "/", ref_id, "_", sample_id, "_normdf.cov.gz")),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = summary_df, file = gzfile(paste0(summary_folder, "/", ref_id, "_", sample_id, "_mapsum.txt.gz")),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##### SNV data
#### Quick verification of vcf format
vcf_lines <- read_lines(vcf_file_path)
for (i in 1:length(vcf_lines)){
  if (str_detect(vcf_lines[i], "^#CHROM")){
    print(paste0("VCF cols: ", vcf_lines[i]))
    
    lines_to_skip <- i
    print(paste0("Skip ", lines_to_skip, " lines"))
    
    break
  }
}

vcf <- read_delim(vcf_file_path, "\t", escape_double = FALSE, 
                  col_names = FALSE, trim_ws = TRUE, skip = lines_to_skip, col_types = cols())


if(length(vcf) > 1){
  colnames(vcf) = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  
  library(plyr)
  vcf <- plyr::ddply(vcf, .(POS), mutate, DP = word(word(INFO,1,sep=";"),2, sep="="),
                     AF = as.numeric(word(word(INFO,2,sep=";"),2, sep="=")),
                     SB = as.numeric(word(word(INFO,3,sep=";"),2, sep="=")),
                     DP4_refFOR = as.numeric(word(word(word(INFO,4,sep=";"),2, sep="="),1,sep=",")),
                     DP4_refREV = as.numeric(word(word(word(INFO,4,sep=";"),2, sep="="),2,sep=",")),
                     DP4_altFOR = as.numeric(word(word(word(INFO,4,sep=";"),2, sep="="),3,sep=",")),
                     DP4_altREV = as.numeric(word(word(word(INFO,4,sep=";"),2, sep="="),4,sep=",")))
}else{
  vcf <- tibble(CHROM  = NA, POS  = NA, ID  = NA, REF  = NA, ALT  = NA, QUAL  = NA, FILTER  = NA, INFO  = NA,
  DP  = NA, AF  = NA, SB  = NA, DP4_refFOR  = NA, DP4_refREV  = NA, DP4_altFOR  = NA, DP4_altREV = NA)
}

vcf$INFO <- NULL
vcf1 <- vcf
vcf1$REF[vcf1$REF == "T"] <- "U"
vcf1$ALT[vcf1$ALT == "T"] <- "U"
vcf1$SNP <- paste0("g.", vcf1$POS,"_", vcf1$REF, ">", vcf1$ALT)
vcf1$sample_id <- sample_id
vcf1$ref_id <- ref_id



write.table(x = vcf1, file = gzfile(paste0(summary_folder, "/", ref_id, "_", sample_id, "_snv.vcf.gz")),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

cat(paste0("\n\nNormalized coverage and mapping summary files written : \n", 
           paste0(summary_folder, "/", ref_id, "_", sample_id, "_normdf.cov.gz\n"),
           paste0(summary_folder, "/", ref_id, "_", sample_id, "_mapsum.txt.gz\n"),
           paste0(summary_folder, "/", ref_id, "_", sample_id, "_snv.vcf.gz")))



