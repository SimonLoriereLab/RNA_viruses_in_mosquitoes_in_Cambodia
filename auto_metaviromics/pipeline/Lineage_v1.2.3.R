#!/usr/bin/env Rscript

##### v1.2.3 current lineage.csv does not have 'species1', column numbers corrected
## "/full_path_to/ncbi_lineages_2021-07-23.csv"

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
df_diamond_raw <- read.delim(args[1], header = FALSE, sep = "\t")


colnames(df_diamond_raw) = c("qseqid", "sseqid", "pident", "length", "mismatch",
                             "gapopen", "qstart", "qend", "sstart", "send", "slen", "evalue",
                             "bitscore", "stitle", "tax_id")

df_diamond_raw_f6 <- df_diamond_raw %>% select(qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)
### split string of arg[1] added due to format of input file with full path 
s <- strsplit(args[1],"/")
s_end <- sapply(s,tail,1)
write.table(x = df_diamond_raw_f6, file = gzfile(file.path(args[3], paste0("f6.",s_end, ".gz"))),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = FALSE)

print(df_diamond_raw[1,])

### remove the extra taxids, pick first
df_diamond_raw$tax_id <- as.numeric(str_remove(as.character(df_diamond_raw$tax_id), ";[[:print:]]*"))

df_lineage_csv <- read.csv(args[2])

df_joined_lineage <- left_join(df_diamond_raw, df_lineage_csv, by = "tax_id")

write.table(x = df_joined_lineage, file = gzfile(file.path(args[3], paste0("lineage_complete.", s_end, ".gz"))),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)

#######
cols_to_collapse <- colnames(df_joined_lineage)[16:length(colnames(df_joined_lineage))]
#######

df_joined_lineage_collapsed <- df_joined_lineage %>% 
  mutate_at(cols_to_collapse, as.character) %>%
  unite("Collapsed_lineage", superkingdom:varietas, 
        sep = "|", na.rm = TRUE, remove = FALSE) %>% 
  rename(realm = clade) %>%
  select(c(1:16,'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 
           'kingdom', 'realm', 'clade1', 'clade2'))

write.table(x = df_joined_lineage_collapsed, file = gzfile(file.path(args[3],paste0("lineage_collapsed.",s_end,".gz"))),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)


