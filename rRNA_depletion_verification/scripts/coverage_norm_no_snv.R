#!/usr/bin/env Rscript --vanilla

library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
cat(paste0("\nArguments passed to the script: ", paste(args, collapse = " "), "\n-----"))

dir<- args[1]
sample_id <- args[2]
ref_id <- args[3]
mapinfo_file_path <- args[4]
cov_file_path <- args[5]
summary_folder <- args[6]


cat(paste0("\nDirectory : ", dir))
cat(paste0("\nSample : ", sample_id))
cat(paste0("\nReference : ", ref_id))
cat(paste0("\nMapping info file path : ", mapinfo_file_path))
cat(paste0("\nCoverage file path : ", cov_file_path))


####
mapinfo_reads <- read_lines(mapinfo_file_path, skip = 16, n_max = 1)
total_reads <- str_extract(mapinfo_reads, pattern = "\\d+")
mapinfo_mapped_reads <- read_lines(mapinfo_file_path, skip = 18, n_max = 1)
mapped_reads <- str_extract(mapinfo_mapped_reads, pattern = "\\d+")

mapinfo_lines <- read_lines(mapinfo_file_path)

for (l in 1:length(mapinfo_lines)){
  if (str_detect(mapinfo_lines[l], "^Coverage info:")){
    l_to_skip = l+3
    break
  }
}

percentage_uncovered <- read_lines(mapinfo_file_path, skip = l_to_skip, n_max = 1)
percentage_uncovered_number <- as.numeric(str_remove_all(str_extract(percentage_uncovered, 
                                                          pattern = " [[:digit:]]+\\.[[:digit:]]+ "),
                                              " "))
percentage_covered <- as.character(100 - percentage_uncovered_number)


####
cov_df <- read.table(cov_file_path, header=FALSE)
colnames(cov_df) = c("ref_id","position","depth")
norm_cov_df <- cov_df
norm_cov_df$depth_norm <- (norm_cov_df$depth/(as.numeric(total_reads)))*10^6
norm_cov_df$sample_id <- sample_id

####
average_cov <- sum(norm_cov_df$depth)/length(norm_cov_df$depth)
norm_average_cov <- (average_cov/(as.numeric(total_reads)))*10^6

####
#summary_df <- tibble(sample_id, ref_id, total_reads, mapped_reads, average_cov, norm_average_cov)
summary_df <- tibble(sample_id, ref_id, total_reads, mapped_reads, average_cov, norm_average_cov, percentage_covered)
####
write.table(x = norm_cov_df, file = gzfile(paste0(summary_folder, "/", ref_id, "_", sample_id, "_normdf.cov.gz")),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = summary_df, file = gzfile(paste0(summary_folder, "/", ref_id, "_", sample_id, "_mapsum.txt.gz")),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


# ##### SNV data
# 
# vcf <- read_delim(vcf_file_path, "\t", escape_double = FALSE, 
#                   col_names = FALSE, trim_ws = TRUE, skip = 18, col_types = cols())
# 
# if(length(vcf) > 1){
#   colnames(vcf) = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
#   
#   library(plyr)
#   vcf <- plyr::ddply(vcf, .(POS), mutate, DP = word(word(INFO,1,sep=";"),2, sep="="),
#                      AF = as.numeric(word(word(INFO,2,sep=";"),2, sep="=")),
#                      SB = as.numeric(word(word(INFO,3,sep=";"),2, sep="=")),
#                      DP4_refFOR = as.numeric(word(word(word(INFO,4,sep=";"),2, sep="="),1,sep=",")),
#                      DP4_refREV = as.numeric(word(word(word(INFO,4,sep=";"),2, sep="="),2,sep=",")),
#                      DP4_altFOR = as.numeric(word(word(word(INFO,4,sep=";"),2, sep="="),3,sep=",")),
#                      DP4_altREV = as.numeric(word(word(word(INFO,4,sep=";"),2, sep="="),4,sep=",")))
# }else{
#   vcf <- tibble(CHROM  = NA, POS  = NA, ID  = NA, REF  = NA, ALT  = NA, QUAL  = NA, FILTER  = NA, INFO  = NA,
#   DP  = NA, AF  = NA, SB  = NA, DP4_refFOR  = NA, DP4_refREV  = NA, DP4_altFOR  = NA, DP4_altREV = NA)
# }
# 
# vcf$INFO <- NULL
# vcf1 <- vcf
# vcf1$REF[vcf1$REF == "T"] <- "U"
# vcf1$ALT[vcf1$ALT == "T"] <- "U"
# vcf1$SNP <- paste0("g.", vcf1$POS,"_", vcf1$REF, ">", vcf1$ALT)
# vcf1$sample_id <- sample_id
# vcf1$ref_id <- sample_id
# 
# 
# 
# write.table(x = vcf1, file = gzfile(paste0(summary_folder, "/", ref_id, "_", sample_id, "_.vcf.gz")),
#             quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

cat(paste0("\n\nNormalized coverage and mapping summary files written : \n", 
           paste0(summary_folder, "/", ref_id, "_", sample_id, "_normdf.cov.gz\n"),
           paste0(summary_folder, "/", ref_id, "_", sample_id, "_mapsum.txt.gz\n")
           ))



