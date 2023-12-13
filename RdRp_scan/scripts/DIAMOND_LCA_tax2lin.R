#!/usr/bin/env Rscript

### works with:
### /full_path_to/db/ncbi_lineages_2022-08-23.csv

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
df_diamond_LCA <- read.delim(args[1], header = FALSE, sep = "\t")
colnames(df_diamond_LCA) = c("qseqid","tax_id", "evalue")


df_lineage_csv <- read.csv(args[2])

df_joined_lineage <- left_join(df_diamond_LCA, df_lineage_csv, by = "tax_id")


### split string of arg[1] added due to format of input file with full path 
s <- strsplit(args[1],"/")
s_end <- sapply(s,tail,1)


write.table(x = df_joined_lineage, file = gzfile(file.path(args[3], paste0("lineage_complete.", s_end, ".gz"))),  
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)

#######
cols_to_collapse <- colnames(df_joined_lineage)[4:length(colnames(df_joined_lineage))]
#######

df_joined_lineage_collapsed <- df_joined_lineage %>% 
  mutate_at(cols_to_collapse, as.character) %>%
  unite("Collapsed_lineage", superkingdom:varietas, 
        sep = "|", na.rm = TRUE, remove = FALSE) %>% 
  rename(realm = clade) %>%
  select(c(1:4,'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 
           'kingdom', 'realm', 'clade1', 'clade2'))

write.table(x = df_joined_lineage_collapsed, file = gzfile(file.path(args[3],paste0("lineage_collapsed.",s_end,".gz"))), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)



