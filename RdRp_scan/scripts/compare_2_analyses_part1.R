library(tidyverse)
library(seqinr)


##### Import RdRp_search (DIAMOND) results
### Contig set A
rdrp_dmnd_results_classified <- "/full_path_to/wd/RdRp_search/analysis/all_families_and_genus_classified/phylo_classified_contigs_and_consensus_w_metadata_correctedAdableaps.txt"
df_dmnd_1 <- read_delim(rdrp_dmnd_results_classified, delim = "\t") %>% select(contig_id = contig_ID, vir_group, phylo_virus_species) %>% distinct()
### Contig set B
rdrp_dmnd_results1_unclassified_contigs <- "/full_path_to/wd/RdRp_search/analysis/DIAMOND_species_level/vir_group_curated_unclassified_viruses.txt"
df_dmnd_2 <- read_delim(rdrp_dmnd_results1_unclassified_contigs, delim = "\t")
df_dmnd_2 <- df_dmnd_2[which(df_dmnd_2$keep == 1), ] %>% select(contig_id, vir_group) %>% distinct()
df_dmnd <- bind_rows(df_dmnd_1, df_dmnd_2)

##### Import RdRp_scan (HMMSCAN) results
df_rdrpscan_0 <- read_delim("/full_path_to/wd/RdRp_scan/2022_07_21_rdrp_results/code_1/hmmscan_tblout_best_AB.txt", delim = "\t")
df_rdrpscan_1 <- df_rdrpscan_0 %>% 
  mutate_all(.funs = str_replace_all, ",", ".") %>%
  mutate_at(c("evalue_full_seq", "score_full_seq", "bias_full_seq", 
              "evalue_best_dom", "score_best_dom", "bias_best_dom"),
            .funs = as.numeric)

df_rdrpscan_1$contig_id <- str_remove(df_rdrpscan_1$rdrp_scan_orf_id, "\\_[[:digit:]]+$")
length(unique(df_rdrpscan_1$rdrp_scan_orf_id))
length(unique(df_rdrpscan_1$contig_id))

####### Filter by e-value
df_rdrpscan <- df_rdrpscan_1[which(df_rdrpscan_1$evalue_full_seq < 0.01),]
df_rdrpscan$taxon <-  str_extract(df_rdrpscan$target, "[[:alpha:]]+vir[[:alpha:]]+") 
length(unique(df_rdrpscan$rdrp_scan_orf_id))
length(unique(df_rdrpscan$contig_id))

####### Remove viruses of prokaryotes
df_rdrpscan <- df_rdrpscan[which(!df_rdrpscan$taxon %in% c("Levivirales", 
                                                      "Mindivirales") |
                                 is.na(df_rdrpscan$taxon)),]

length(unique(df_rdrpscan$rdrp_scan_orf_id))
length(unique(df_rdrpscan$contig_id))



####### Filter by ORF size
### run this once
# orf_nr_df <- tibble(contig_id = str_remove(str_extract(names(orfs_nr), 
#                                                     "^[[:alnum:]\\.\\_]+"),
#                                         "\\_[[:digit:]]+$"), 
#                  emboss_getorf_id = str_extract(names(orfs_nr),
#                                                 "^[[:alnum:]\\.\\_]+"), 
#                  orf_start = str_extract(names(orfs_nr),
#                                          "(?<=\\[)[[:digit:]]+"),
#                  orf_end = str_extract(names(orfs_nr),
#                                        "[[:digit:]]+(?=\\])"),
#                  orf_sense = str_extract(names(orfs_nr),
#                                          "REVERSE"))
# orf_nr_df[which(is.na(orf_nr_df$orf_sense)),]$orf_sense <- "FORWARD"
# 
# write.table(x = orf_nr_df, 
#             file = "/full_path_to/wd/RdRp_scan/analysis/all_emboss_getorfs_nr.txt", 
#             append = F, quote = F, sep = "\t", row.names = F, col.names = T)

orf_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/all_emboss_getorfs.txt", delim = "\t")
length(unique(orf_df$emboss_getorf_id))
length(unique(orf_df$contig_id))

orf_nr_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/all_emboss_getorfs_nr.txt", delim = "\t")

length(unique(orf_nr_df$emboss_getorf_id))
length(unique(orf_nr_df$contig_id))





df_rdrpscan <- left_join(df_rdrpscan, orf_nr_df, by = c("rdrp_scan_orf_id" = "emboss_getorf_id", "contig_id"))


df_rdrpscan$orf_len_nt <- abs(as.numeric(df_rdrpscan$orf_end) - as.numeric(df_rdrpscan$orf_start))


df_rdrpscan <- df_rdrpscan[which(df_rdrpscan$orf_len_nt >= 900),]


cat(paste0(sort(unique(df_rdrpscan$target)), collapse = "\n"))


length(sort(unique(df_rdrpscan$contig_id)))
length(sort(unique(df_rdrpscan$rdrp_scan_orf_id)))

write.table(x = sort(unique(df_rdrpscan$contig_id)), 
            file = "/full_path_to/wd/RdRp_scan/analysis/all_hmmscan_filtered_contig_ids.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = F)














##### Start comparisons
##### Find intersection of two results - set G
df_rdrpscan_match <- df_rdrpscan[which(df_rdrpscan$contig_id %in% df_dmnd$contig_id),]
confirmed_with_2_analyses <- sort(unique(df_rdrpscan_match$contig_id))

### Detected with diamond but not hmmscan, why?
dmnd_not_hmm <- setdiff(df_dmnd$contig_id, confirmed_with_2_analyses)
### Chimeric contigs
chimeric_contigs <- c("F1416_NODEs_39_8_2_finalCxFVlike", 
                      "F1576_NODEs_2_25945_25414_25622_8_26117_15_chimera", 
                      "F1052_NODE_13_663_chimera", 
                      "F1506_NODEs_12_58_205")
chimeric_contig_nodes <- c("F1416_NODE_39", "F1416_NODE_8", "F1416_NODE_2",
                           "F1576_NODE_2", "F1576_NODE_25945", "F1576_NODE_25414", "F1576_NODE_25622", "F1576_NODE_8", "F1576_NODE_26117", "F1576_NODE_15",
                           "F1052_NODE_13", "F1052_NODE_663",
                           "F1506_NODE_12", "F1506_NODE_58", "F1506_NODE_205")
df_rdrpscan$contig_id_node <- str_extract(df_rdrpscan$contig_id, "^[[:alnum:]]+\\_NODE_[[:digit:]]+")
df_rdrpscan_match_chimeric <- df_rdrpscan[which(df_rdrpscan$contig_id_node %in% chimeric_contig_nodes),]
df_rdrpscan_match_chimeric$contig_id_node
chimeric_contigs_of_set_G_nodes <- tibble(contig_id_node = c("F1052_NODE_13", 
                                                       "F1416_NODE_2", 
                                                       "F1506_NODE_58", 
                                                       "F1576_NODE_2"), 
                                    contig_id_set_AuB = c("F1052_NODE_13_663_chimera", 
                                                          "F1416_NODEs_39_8_2_finalCxFVlike",
                                                          "F1506_NODEs_12_58_205",
                                                          "F1576_NODEs_2_25945_25414_25622_8_26117_15_chimera"))
chimeric_contigs_of_set_G <- left_join(chimeric_contigs_of_set_G_nodes, df_rdrpscan_match_chimeric) %>% select(contig_id, contig_id_set_AuB)
### Update set G
confirmed_with_2_analyses_1 <- tibble(contig_id = confirmed_with_2_analyses, contig_id_set_AuB = confirmed_with_2_analyses)
set_G <- bind_rows(confirmed_with_2_analyses_1, chimeric_contigs_of_set_G)




### Detected with diamond but not hmmscan, not chimeric contigs, why?
dmnd_not_hmm_1 <- setdiff(df_dmnd$contig_id, set_G$contig_id_set_AuB)

#### Run extract_dmnd_hits_for_non-hmmscan_hits.R on cluster, to subset all the diamond hits for all HMM-SCAN hits
#### Import an intersection of all diamond with hmmscan hits
dmnd_hits_not_hmm_file <- "/full_path_to/wd/RdRp_scan/analysis/lineage_complete.DIAMOND_hits_non-hmmscan.txt"
n_of_lines <- as.numeric(system(paste0("cat ", dmnd_hits_not_hmm_file, " | wc -l"), intern = T))
dmnd_hits_not_hmm_df <- read_delim(file = dmnd_hits_not_hmm_file,
                  delim = "\t",
                  num_threads = 2,
                  guess_max = n_of_lines)

unique(dmnd_hits_not_hmm_df$qseqid)
length(unique(dmnd_hits_not_hmm_df$qseqid))

#### Investigating why these hits are missing
### Reason 1:
### In DIAMOND search I used only PB2 of Orthomyxoviruses and 
### it's possible in the RdRP scan propfiles PA and PB1 are more prominent subunits.
### There is no easy way to verify this, therefore, for now,
### I can assume that whatever RdRP segment HMM-scan picked up, 
### it will appear in the same samples as diamond-detected PB2. 
### This means, I will have to update Orthomyxo phylogenies with PB1 segment.
### Unclear reasons for other hits.

### Set aside set M

set_M <- tibble(contig_id = sort(unique(dmnd_hits_not_hmm_df$qseqid)),
                contig_id_set_AuB = sort(unique(dmnd_hits_not_hmm_df$qseqid)))
write.table(x = set_M, 
            file = "/full_path_to/wd/RdRp_scan/analysis/contig_set_M.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)

intersect(set_G$contig_id, set_M$contig_id) #### check


df_dmnd_1[which(df_dmnd_1$contig_id %in% set_M$contig_id),]


write.table(x = set_G, 
            file = "/full_path_to/wd/RdRp_scan/analysis/contig_set_G.txt", 
            append = F, quote = F,sep = "\t",row.names = F, col.names = T)


######################################################
######################################################



###### Set F
df_set_F <- df_rdrpscan[which(!df_rdrpscan$contig_id %in% set_G$contig_id),]
new_rdrpscan_contigs <- sort(unique(df_set_F$contig_id))

write.table(x = df_set_F, 
            file = "/full_path_to/wd/RdRp_scan/analysis/contig_set_F_hmmscan_df.txt", 
            append = F, quote = F,sep = "\t",row.names = F, col.names = T)

write.table(x = df_set_F$contig_id, 
            file = "/full_path_to/wd/RdRp_scan/analysis/contig_set_F.txt", 
            append = F, quote = F,sep = "\t",row.names = F, col.names = T)











