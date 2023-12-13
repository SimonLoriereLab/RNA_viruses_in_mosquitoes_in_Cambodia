#!/usr/bin/env Rscript --vanilla

library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
cat(paste0("\nArguments passed to the script: ", paste(args, collapse = " "), "\n-----"))

dir<- args[1]
sample_id <- args[2]
map_ref_tab <- args[3]
ref_dir <- args[4]

if (map_ref_tab != "all_to_all"){
  print("specific reference X sample table provided")
  cat(paste0("\nDirectory : ", dir))
  cat(paste0("\nSample : ", sample_id))
  cat(paste0("\nReference table : ", map_ref_tab))
  
  df <- read.delim(map_ref_tab)
  ref_str <- as.character(df[which(df$sample == sample_id),]$ref)
  ref_vec <- str_split(ref_str, pattern = ",")[[1]]
  
  write.table(x = ref_vec, file = paste0("reflist_sample_", sample_id, ".txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  cat(paste0("\n\nReference list written : ", 
             paste0("reflist_sample_", sample_id, ".txt\n\n")))
  
} else {
  print("all refs X all samples")
  cat(paste0("\nDirectory : ", dir))
  cat(paste0("\nSample : ", sample_id))
  cat(paste0("\nReference table : ", map_ref_tab))
  cat(paste0("\nReference directory : ", ref_dir))
  
  
  ref_vec <- str_remove_all(list.files(path = ref_dir, ".fasta$"), ".fasta$")
  
  cat(ref_vec)
  
  write.table(x = ref_vec, file = paste0("reflist_sample_", sample_id, ".txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  cat(paste0("\n\nReference list written : ", 
             paste0("reflist_sample_", sample_id, ".txt\n\n")))
  
}
