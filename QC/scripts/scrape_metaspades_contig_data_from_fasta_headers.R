library(tidyverse)

path_to_output <- "/full_path_to/wd/QC/analysis/metaspades/"
dir.create(path_to_output)

#### prep list of samples
mdf <- read_delim("/full_path_to/wd/QC/metadata/CMVM_full_metadata_updated_20211021.tsv", delim = "\t") %>% 
  select(Sample, Demultiplexed_and_trimmed)


samples_to_scrape <- mdf[which(mdf$Demultiplexed_and_trimmed == 1),]$Sample

metaspades_path_p1 <- "/full_path_to/wd/auto_metaviromics/analyzed_samples/"
metaspades_path_p2 <- "/metaSPAdes/metaSPAdes_output_"
metaspades_path_p3 <- "/contigs.fasta"


not_found <- c("")
count = 0
for (s in samples_to_scrape){
  path_to_contig_fa <- paste0(metaspades_path_p1, s, metaspades_path_p2, s, metaspades_path_p3)
  print(path_to_contig_fa)
  if (file.exists(path_to_contig_fa)){
    cont_fa <- read_lines(file = path_to_contig_fa)

    if (length(cont_fa) > 0){
      headers <- str_remove(cont_fa[str_detect(cont_fa, ">")], ">")
      contig_len <- as.numeric(str_remove(str_remove(headers, "NODE_[[:digit:]]+_length_"), "_cov_[[:print:]]+"))
      contig_kmer_cov <- as.numeric(str_remove(str_remove(headers, "NODE_[[:digit:]]+_length_[[:digit:]]+_cov_"), "_[[:print:]]+"))
      s_temp_df <- tibble(sample_id = s, contig_id = headers, contig_len, contig_kmer_cov)
      if (count == 0){
        metaspades_df <- s_temp_df
      }else{
        metaspades_df <- bind_rows(metaspades_df, s_temp_df)
      }
      count = count + 1
      print(count)
    }
  }else{
    not_found <- c(not_found, s)
  }
  
}

date <- str_remove_all(Sys.Date(), "-")

write.table(x = metaspades_df, file = gzfile(paste0(path_to_output, "metaspades_scraped_stats_CMVM_", date, ".txt.gz")),
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)

write.table(x = not_found[!not_found == ""], file = paste0(path_to_output, "metaspades_not_found_CMVM_", date, ".txt"),
            append = F, quote = F, sep = "\t", row.names = F, col.names = F)
