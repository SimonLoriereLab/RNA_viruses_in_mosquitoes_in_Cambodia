library(tidyverse)
library(seqinr)

########## Check which of the H contigs were discarded during diamond result curation for 
########## unclassified species (Set E), genera (set J), and families (set J) 
########## check them again.


### Contig set H df
H_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/contig_set_H_hmmscan_df.txt", delim = "\t")
length(unique(H_df$contig_id))

####################### Contig set E
## Contaminants
fns1 <- list.files("/full_path_to/wd/RdRp_search/analysis/DIAMOND_species_level/contaminants", full.names = T)
count = 0
for (fn in fns1){
  df <- read_delim(fn, delim = "\t") %>% arrange(desc(bitscore))
  if (count == 0){
    df_contaminants <- df
  }else{
    df_contaminants <- bind_rows(df_contaminants,df)
  }
  count = count + 1
}

count1 = 0
for (q in unique(df_contaminants$qseqid)){
  dfq <- df_contaminants[which(df_contaminants$qseqid == q),]
  dfq1 <- dfq[1,]
  if (count1 == 0){
    df_contaminants1 <- dfq1
  }else{
    df_contaminants1 <- bind_rows(df_contaminants1,dfq1)
  }
  count1 = count1 + 1
}

## check if in H
intersect(H_df$contig_id, unique(df_contaminants1$qseqid))

## RT
fns2 <- list.files("/full_path_to/wd/RdRp_search/analysis/DIAMOND_species_level/RT", full.names = T)
count = 0
for (fn in fns2){
  df <- read_delim(fn, delim = "\t") %>% arrange(desc(bitscore))
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

## check if in H
prev_discarded_RT_in_H <- intersect(H_df$contig_id, unique(df_RT1$qseqid))
df_RT1_in_H <- df_RT1[which(df_RT1$qseqid %in% prev_discarded_RT_in_H),]
H_df_in_RT1 <- H_df[which(H_df$contig_id %in% prev_discarded_RT_in_H),]

unique(H_df_in_RT1$contig_id)
length(unique(H_df_in_RT1$contig_id))
##### Checking why they were considered as RT:
##### Both Aedes aegypti To virus 1 and 2 belong to Metaviridae
##### that could be verified by blastp search and from the publication involving one of them.

### Removing these contigs from the set H
H_df_upd1 <- H_df[which(!H_df$contig_id %in% prev_discarded_RT_in_H),]

## initially considered as too short
fns3 <- list.files("/full_path_to/wd/RdRp_search/analysis/DIAMOND_species_level/too_short", full.names = T)
count = 0
for (fn in fns3){
  df <- read_delim(fn, delim = "\t") %>% arrange(desc(bitscore))
  if (count == 0){
    df_too_short <- df
  }else{
    df_too_short <- bind_rows(df_too_short,df)
  }
  count = count + 1
}

count1 = 0
for (q in unique(df_too_short$qseqid)){
  dfq <- df_too_short[which(df_too_short$qseqid == q),]
  dfq1 <- dfq[1,]
  if (count1 == 0){
    df_too_short1 <- dfq1
  }else{
    df_too_short1 <- bind_rows(df_too_short1,dfq1)
  }
  count1 = count1 + 1
}

## check if in H
intersect(H_df_upd1$contig_id, unique(df_too_short1$qseqid))


################################################################################

### Contig set J.
### Finding which of the H set contigs were discarded as too short for a given virus group during DIAMOND results analysis.


###### !!!!
###### This part could only be run in the original folder with the original results, because some of the hittable files are very large
# dirs <- c(list.dirs(path = "/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level", 
#             full.names = T, recursive = F),
#   list.dirs(path = "/full_path_to/wd/RdRp_search/analysis/DIAMOND_genus_level", 
#             full.names = T, recursive = F))
# dirs <- dirs[which(str_detect(str_remove_all(dirs, "[[:alnum:]\\_]*\\/"), "vir"))]
# 
# 
# 
# count_all = 0
# for (d in dirs){
#   print(d)
# 
#   vir_group <- str_extract(d, "(?<=\\/)[[:alpha:]]+?$")
#   if(str_detect(vir_group, "virus")){
#     all_contigs_hittable <- paste0(vir_group, "_genus_rdrp_polyprot_hittable.txt")
#   }else{
#     all_contigs_hittable <- paste0(vir_group, "_family_rdrp_polyprot_hittable.txt")
#   }
# 
#   
#   
#   d_upstream <- str_remove(d, vir_group)
#   
#   
#   hit_file <- paste0(d_upstream, all_contigs_hittable)
#   n_of_lines <- as.numeric(system(paste0("cat ", hit_file , " | wc -l"), intern = T))
#   df_hits <- read_delim(file = hit_file,
#                                      delim = "\t",
#                                      num_threads = 2,
#                                      guess_max = n_of_lines)
# 
#   df_hits_in_H <- df_hits[which(df_hits$qseqid %in% H_df_upd1$contig_id),]
#   
#   if (nrow(df_hits_in_H) > 0){
#     df_hits_in_H$vir_group <- vir_group
#     count_in_H = 0
#     for (q in unique(df_hits_in_H$qseqid)){
#       df_hits_in_H_q <- df_hits_in_H[which(df_hits_in_H$qseqid == q),] %>% arrange(desc(bitscore))
#       df_hits_in_H_q_bb <- df_hits_in_H_q[1,]
#       if (count_in_H == 0){
#         df_hits_in_H_bb <- df_hits_in_H_q_bb
#       }else{
#         df_hits_in_H_bb <- bind_rows(df_hits_in_H_bb, df_hits_in_H_q_bb)
#       }
#       count_in_H = count_in_H + 1
#     }
#     if (count_all == 0){
#       df_hits_in_H_bb_all <- df_hits_in_H_bb
#     }else{
#       df_hits_in_H_bb_all <- bind_rows(df_hits_in_H_bb_all, df_hits_in_H_bb)
#     }
#     count_all = count_all + 1
#   }
# }
# 
# ### Make sure again none of these are in G or M
# G_set <- read_delim("/full_path_to/wd/RdRp_scan/analysis/contig_set_G.txt", delim = "\t")
# length(unique(G_set$contig_id))
# 
# M_set <- read_delim("/full_path_to/wd/RdRp_scan/analysis/contig_set_M.txt", delim = "\t")
# intersect(df_hits_in_H_bb_all$qseqid, G_set$contig_id)
# intersect(df_hits_in_H_bb_all$qseqid, M_set$contig_id)

#### Look into each one separately, decide which ones to keep.

# write.table(x = df_hits_in_H_bb_all, 
#             file = "/full_path_to/wd/RdRp_scan/analysis/J_intersect_H_diamond_best_bitscore_hits.txt", 
#             append = F, quote = F,sep = "\t", row.names = F, col.names = T)
df_hits_in_H_bb_all <- read_delim("/full_path_to/wd/RdRp_scan/analysis/J_intersect_H_diamond_best_bitscore_hits.txt", 
                                          delim = "\t")
### check genome/segment sizes and discard contigs < 50% of the expected size

df_hits_in_H_bb_all_checked <- read_delim("/full_path_to/wd/RdRp_scan/analysis/CMVM_bioinformatics - J_intersect_H_check.tsv", 
           delim = "\t")

intersect(unique(df_hits_in_H_bb_all_checked$qseqid), unique(df_hits_in_H_bb_all$qseqid))
setdiff(unique(df_hits_in_H_bb_all_checked$qseqid), unique(df_hits_in_H_bb_all$qseqid))
setdiff(unique(df_hits_in_H_bb_all$qseqid), unique(df_hits_in_H_bb_all_checked$qseqid))

contigs_to_discard_JintH <- df_hits_in_H_bb_all_checked[which(df_hits_in_H_bb_all_checked$keep ==  0),]$qseqid

length(contigs_to_discard_JintH)

N_df <- H_df_upd1[which(!H_df_upd1$contig_id %in% contigs_to_discard_JintH),] %>% distinct()

write.table(x = N_df, 
            file = "/full_path_to/wd/RdRp_scan/analysis/N_df.txt", 
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
N_df1 <- N_df %>% select(contig_id) %>% distinct()
write.table(x = N_df1, 
            file = "/full_path_to/wd/RdRp_scan/analysis/contig_set_N.txt", 
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)







