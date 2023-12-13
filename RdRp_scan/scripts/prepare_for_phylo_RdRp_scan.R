library(tidyverse)
library(seqinr)

R_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/R_orf_df_curated_with_chimeric_contig_id.txt", delim = "\t")

orf_ids <- sort(unique(R_df$orf_id))

orfs <-  read.fasta("/full_path_to/wd/RdRp_scan/analysis/R_orfs_curated.fasta", 
           seqtype = "AA", strip.desc = T)

R_orfs <- orfs[orf_ids]

write.fasta(sequences = R_orfs, names = names(R_orfs), 
            file.out = "/full_path_to/wd/RdRp_scan/analysis/R_232orfs_aa.fasta")
