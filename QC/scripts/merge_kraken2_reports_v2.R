
###### Run on MAESTRO HPC

library(tidyverse)

path_to_output <- "/full_path_to/wd/QC/analysis/Kraken2/"
dir.create(path_to_output)

#### prep list of samples
samples_to_scrape <- read_lines("/full_path_to/wd/QC/metadata/sample_list1.tsv", num_threads = 95)

kraken2_path_p1 <- "/full_path_to/wd/auto_metaviromics/analyzed_samples/"
kraken2_path_p2 <- "/Kraken2/"
kraken2_path_p3 <- "_ntDB_kraken2_report.txt" ### nt DB


not_found <- c("")
count = 0
for (s in samples_to_scrape){
  path_to_kraken2_report <- paste0(kraken2_path_p1, s, kraken2_path_p2, s, kraken2_path_p3)
  print(path_to_kraken2_report)
  if (file.exists(path_to_kraken2_report)){
    kraken2_report <- read_delim(path_to_kraken2_report, 
                                 delim = "\t", 
                                 col_names = c("percent_frags", "num_frags", 
                                               "num_frags_directly_to_taxon", 
                                               "num_minimizers_for_taxon",
                                               "num_distinct_minimizers_for_taxon", 
                                               "rank_code", "taxid", "taxon_name"), num_threads = 95)
    kraken2_report$sample_id <- s
      if (count == 0){
        kraken2_df <- kraken2_report
      }else{
        kraken2_df <- bind_rows(kraken2_df, kraken2_report)
      }
      count = count + 1
      print(count)
  }else{
    not_found <- c(not_found, s)
  }
  
}

date <- str_remove_all(Sys.Date(), "-")

write.table(x = kraken2_df, file = gzfile(paste0(path_to_output, "kraken2_reports_CMVM_", date, ".txt.gz")),
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)

write.table(x = not_found[!not_found == ""], file = paste0(path_to_output, "kraken2_not_found_CMVM_", date, ".txt"),
            append = F, quote = F, sep = "\t", row.names = F, col.names = F)

