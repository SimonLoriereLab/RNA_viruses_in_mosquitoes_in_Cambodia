library(tidyverse)
library(seqinr)

####### Removing already annotated known contaminants from set F

### Contig set F df
F_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/contig_set_F_hmmscan_df.txt", delim = "\t")
length(unique(F_df$contig_id))

######## Set T (containing RT)

fnsT <- list.files("/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/not_analyzed_further/with_RT",
                   full.names = T)
count = 0
for (fn in fnsT){
  df <- read_delim(fn, delim = "\t", guess_max = 10000) %>% arrange(desc(bitscore))
  if (count == 0){
    df_RT <- df
  }else{
    df_RT <- bind_rows(df_RT,df)
  }
  count = count + 1
}

count1 = 0
for (q in unique(df_RT$qseqid)){
  dfq <- df_RT[which(df_RT$qseqid == q),]
  dfq1 <- dfq[1,]
  if (count1 == 0){
    df_RT1 <- dfq1
  }else{
    df_RT1 <- bind_rows(df_RT1,dfq1)
  }
  count1 = count1 + 1
}

## check if in F
RT_in_F <- intersect(F_df$contig_id, unique(df_RT1$qseqid))
df_RT2 <- df_RT1[which(df_RT1$qseqid %in% RT_in_F),]
#### removing these RT-containing viruses, mostly Metaviridae which are retrotransposons
F_df_upd1 <- F_df[which(!F_df$contig_id %in% RT_in_F),]

######### Set C


contaminants <- c("/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/not_analyzed_further/contaminants/Togaviridae_family_rdrp_polyprot_hittable.txt",
  "/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/not_analyzed_further/contaminants/Coronaviridae/Coronaviridae_unambiguous_best_bitscore_hits.txt",
  "/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/not_analyzed_further/contaminants/Coronaviridae/Coronaviridae_family_ambiguous_hits_check_later.txt",
  "/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/not_analyzed_further/contaminants/Paramyxoviridae/Paramyxoviridae_family_ambiguous_best_bitscore_hits.txt")

count = 0
for (f in contaminants){
  temp_df <-read_delim(f, delim = "\t") %>% select(contig_id = qseqid, family, species, stitle)
  if (count == 0){
    C_df <- temp_df
  }else{
    C_df <-  bind_rows(C_df, temp_df) %>% distinct()
  }
  count = count + 1
}

FintC_df <- F_df_upd1[which(F_df_upd1$contig_id %in% unique(C_df$contig_id)),]
FnotC_df <- F_df_upd1[which(!F_df_upd1$contig_id %in% unique(C_df$contig_id)),]
sort(unique(C_df$contig_id))

##### No contaminants found in F because of the minsize for orfs 300 aa
### This one is long enough "PBS1c_NODE_11_length_1288_cov_132.524736" 
C_df[which(C_df$contig_id == "PBS1c_NODE_11_length_1288_cov_132.524736"),]
orf_nr_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/all_emboss_getorfs_nr.txt", delim = "\t")
orf_nr_df[which(orf_nr_df$contig_id == "PBS1c_NODE_11_length_1288_cov_132.524736"),]
### has a long enough ORF but it's not RdRP
write.table(x = sort(unique(C_df$contig_id)), 
            file = "/full_path_to/wd/RdRp_scan/analysis/contig_set_C.txt", 
            append = F, quote = F,sep = "\t",row.names = F, col.names = F)


########## Check which of the F contigs were deemed too short during diamond result curation and check them again.
dmnd_ts_dir <- "/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/not_analyzed_further/too_short"
dmnd_ts_files <- list.files(path = dmnd_ts_dir, full.names = T, pattern = "hittable.txt")

count1 = 0
for (f in dmnd_ts_files){
  temp_df <-read_delim(f, delim = "\t") %>% select(contig_id = qseqid, family, species, stitle, length, pident, bitscore)
  if (count1 == 0){
    D_df <- temp_df
  }else{
    D_df <-  bind_rows(D_df, temp_df) %>% distinct()
  }
  count1 = count1 + 1
}

FintD_df <- F_df_upd1[which(F_df_upd1$contig_id %in% unique(D_df$contig_id)),]
FnotD_df <- F_df_upd1[which(!F_df_upd1$contig_id %in% unique(D_df$contig_id)),]


FintD_df$contig_id
DintF_df_upd1 <- D_df[which(D_df$contig_id %in% FintD_df$contig_id),]

###check if contigs are longer than 50% of expected length of the best hits
## < 50% of expected length for best hits or family
temp_df <- DintF_df_upd1[which(!DintF_df_upd1$contig_id %in% c(
                                                   "F1362_NODE_102_length_3099_cov_4.507556")),]
## >= 50% of expected length for best hits or family
temp_df1 <- DintF_df_upd1[which(DintF_df_upd1$contig_id %in% c(
                                                               "F1362_NODE_102_length_3099_cov_4.507556")),]

### Almost all contigs considered too short are indeed too short, 
### but "F1362_NODE_102_length_3099_cov_4.507556" looks like it's near full-size segment,

### create contig set H = F - T - C - D + "F1362_NODE_102_length_3099_cov_4.507556"
H_df <- bind_rows(FnotD_df, F_df_upd1[which(F_df_upd1$contig_id == "F1362_NODE_102_length_3099_cov_4.507556"),])

write.table(x = H_df, 
            file = "/full_path_to/wd/RdRp_scan/analysis/contig_set_H_hmmscan_df.txt", 
            append = F, quote = F,sep = "\t",row.names = F, col.names = T)

write.table(x = H_df$contig_id, 
            file = "/full_path_to/wd/RdRp_scan/analysis/contig_set_H.txt", 
            append = F, quote = F,sep = "\t",row.names = F, col.names = T)




