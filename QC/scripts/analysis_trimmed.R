library(tidyverse)
library(ggrepel)

path_to_inout <- "/full_path_to/wd/QC/analysis/Trimmomatic/"
date <- str_remove_all(Sys.Date(), "-")

#### preprocessing
mdf <- read_delim("/full_path_to/wd/QC/metadata/CMVM_full_metadata_updated_20211021.tsv", delim = "\t") %>% 
  select(sample_id = Sample, Batch, LibPool, NTC, LibPool_prep_by = LibPool_prep_by.y, Run_CD, Run_PF, Run_EY, Location, genus_species)

# df_trim <- read_delim(file = paste0(path_to_inout, "trimmomatic_scraped_stats_CMVM_20211005.txt"), delim = "\t")
# df_trim_more <- read_delim(file = paste0(path_to_inout, "trimmomatic_scraped_stats_CMVM_20211122.txt"), delim = "\t")
# df_trim <- bind_rows(df_trim, df_trim_more)


df_trim <- read_delim(file = paste0(path_to_inout, "trimmomatic_scraped_stats_CMVM_20220908.txt"), delim = "\t")
df_trim_more <- read_delim(file = paste0(path_to_inout, "trimmomatic_scraped_stats_CMVM_20220908_1.txt"), delim = "\t")
df_trim <- bind_rows(df_trim, df_trim_more)

class(df_trim$err_fmtime)

df_trim$bis_seq <- 0
count_d <- 1
for (s in unique(df_trim$sample_id)){
  temp_df_trim <- df_trim[which(df_trim$sample_id ==s),]
  
  if (nrow(temp_df_trim) > 1){
    print("duplicates, take later file")
    temp_df_trim <- temp_df_trim %>% arrange(desc(err_fmtime))
    print(temp_df_trim$err_fmtime)
    temp_df_trim[1,]
  }
  
  if (str_detect(s, "bis")){
    print(s)
    temp_df_trim$sample_id <- str_remove(temp_df_trim$sample_id, "bis")
    temp_df_trim$bis_seq <- 1
  }
  
  if(count_d == 1){
    df_trim_dedup <- temp_df_trim
  }else{
    df_trim_dedup <- bind_rows(df_trim_dedup, temp_df_trim)
  }
  count_d = count_d + 1
  
  
}

df_trim <- df_trim_dedup
df_trim$num_qualtrim_kept_reads <- df_trim$num_qualtrim_kept_PE_reads + df_trim$for_only_surv + df_trim$rev_only_surv


unique(df_trim$sample_id)
unique(df_trim$err_fname)
nrow(df_trim)


sum(df_trim$total_num_of_PE_reads)/10^9
sum(df_trim$num_qualtrim_kept_reads)/10^9
sum(df_trim$num_qualtrim_kept_PE_reads)/10^9

#cat(paste0(colnames(df_trim), collapse = ", "))

df_trim <- left_join(df_trim, mdf, by = "sample_id")


df_trim$Batch <-  factor(as.character(df_trim$Batch), levels = as.character(sort(unique(df_trim$Batch))))
class(df_trim$Batch)
levels(df_trim$Batch)


df_trim$LibPool <-  factor(as.character(df_trim$LibPool), levels = as.character(sort(unique(df_trim$LibPool))))
class(df_trim$LibPool)
levels(df_trim$LibPool)

df_trim$NTC <- as.factor(as.character(df_trim$NTC))


#### plots


p0 <- ggplot()
p0 <- p0 + geom_density(data = df_trim, mapping = aes(x = log10(total_num_of_PE_reads)), 
                        color = "grey", fill = "grey", alpha = 1)
p0 <- p0 + geom_density(data = df_trim, mapping = aes(x = log10(num_qualtrim_kept_reads)), 
                        color = "firebrick3", fill = "firebrick3", alpha = 0.1)
p0 <- p0 + xlab("Reads/sample, log10")
p0 <- p0 + theme_classic(base_size = 8)
p0

ggsave(filename = paste0(path_to_inout, "dens_total_PE_reads_total_vs_kept_", date, ".pdf"),
       plot = p0, device = "pdf", width = 7, height = 4, units = "cm")



df_trim$color_bad_libpools <- "OK"
df_trim[which(df_trim$LibPool %in% c("9", "15")),]$color_bad_libpools <- "Illumina reagent error"
df_trim[which(df_trim$LibPool %in% c("9", "15") & df_trim$bis_seq == 1),]$color_bad_libpools <- "re-sequenced, OK"

p0 <- ggplot()
p0 <- p0 + geom_point(data = df_trim, mapping = aes(x = log10(num_qualtrim_kept_reads),
                                                    y = log10(total_num_of_PE_reads),
                                                    color = color_bad_libpools), shape = 1, size = 2)
p0 <- p0 + ylab("Total number of PE reads, log10")
p0 <- p0 + xlab("Total number of PE reads post trimming, log10")

p0 <- p0 + theme_minimal(base_size = 11)
p0 <- p0 + theme(legend.position = "bottom", legend.title = element_blank())
p0

ggsave(filename = paste0(path_to_inout, "total_vs_kept_", date, ".pdf"),
       plot = p0, device = "pdf", width = 15, height = 15, units = "cm")


# #### flag the poorly sequenced samples
# df_trim$flag_poor_read_qual <- 0
# 
# df_trim[which(df_trim$prop_qualtrim_kept_PE_reads < 0.8),]$flag_poor_read_qual <- 1
# 
# df_trim[which(df_trim$prop_qualtrim_kept_PE_reads < 0.8),]$sample_id
# 
# df_trim$flag_poor_seq_depth <- 0
# 
# df_trim[which(df_trim$total_num_of_PE_reads < 1000000),]$flag_poor_seq_depth <- 1
# 
# df_trim[which(df_trim$total_num_of_PE_reads < 1000000),]$sample_id
# 
# df_trim_short_sum <- df_trim %>% select(sample_id, total_num_of_PE_reads, num_qualtrim_kept_reads, prop_qualtrim_kept_PE_reads, flag_poor_seq_depth, flag_poor_read_qual)
# cat(paste0(colnames(df_trim_short_sum), collapse = ", "))
# write.table(x = df_trim_short_sum, file = paste0(path_to_inout, "trimmomatic_summary_CMVM_", date, ".txt"),
#             append = F, quote = F, sep = "\t", row.names = F, col.names = T)





