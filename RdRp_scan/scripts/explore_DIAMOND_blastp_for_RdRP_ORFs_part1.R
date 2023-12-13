library(tidyverse)
library(seqinr)

##### Import DIAMOND blastp results
LCA_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/lineage_collapsed.DIAMOND_blastp_hits_nrcontigs_rdrp_LCA.txt.gz",
                     delim = "\t")

LCA_df$orf_id <- LCA_df$qseqid
LCA_df$contig_id <- str_remove(LCA_df$qseqid, "_[[:digit:]]+$")

unique(LCA_df$contig_id)

LCA_df_non_cellular <- LCA_df[which(LCA_df$superkingdom == "Viruses" |
                                      is.na(LCA_df$superkingdom) |
                                      LCA_df$tax_id == 0 &
                                      !str_detect(LCA_df$Collapsed_lineage, "cellular")), ]

LCA_df_non_cellular <- LCA_df_non_cellular[which(!str_detect(LCA_df_non_cellular$Collapsed_lineage, "cellular") |
                                                   is.na(LCA_df_non_cellular$Collapsed_lineage)), ]


N_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/contig_set_N.txt", delim = "\t")
G_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/contig_set_G.txt", delim = "\t")

#### Verify if all G set contigs are detected as viral
sum(sort(unique(G_df$contig_id)) == sort(intersect(unique(LCA_df_non_cellular$contig_id), unique(G_df$contig_id)))) == length(unique(G_df$contig_id))

#### Verify if DIAMOND blastp results are available for all orfs
orfs <- read.fasta("/full_path_to/wd/RdRp_scan/analysis/all_samples_renamed_contigs_redundancy_removed_code1_aa_nr_rdrp.fasta")
setdiff(sort(unique(LCA_df$qseqid)), sort(unique(names(orfs))))
setdiff(sort(unique(names(orfs))), sort(unique(LCA_df$qseqid)))

### Check if some of the contigs have both cellular and non-cellular contigs
LCA_df_ambiguous_cellular <- LCA_df[which(str_detect(LCA_df$Collapsed_lineage, "cellular") &
                                                        LCA_df$contig_id %in% LCA_df_non_cellular$contig_id), ]
LCA_df_ambiguous_non_cellular <-  LCA_df_non_cellular[which(LCA_df_non_cellular$contig_id %in% LCA_df_ambiguous_cellular$contig_id),]
LCA_df_non_cellular1 <- bind_rows(LCA_df_non_cellular, LCA_df_ambiguous_cellular)
LCA_df_non_cellular1$possibly_cellular_contig <- 0
LCA_df_non_cellular1[which(LCA_df_non_cellular1$contig_id %in% LCA_df_ambiguous_cellular$contig_id),]$possibly_cellular_contig <- 1





### Import RdRP scan taxonomy
GN_hmmscan_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/GN_hmmscan_df.txt",
                     delim = "\t") %>% select(-accession1, -accession2, -description_of_target)


### L set is the selected non-cellular by DIAMOND blastp GN set orfs
L_set_df <- left_join(LCA_df_non_cellular1, GN_hmmscan_df, by = c("orf_id"="rdrp_scan_orf_id", "contig_id"))



############### Add DIAMOND blastp best hit taxonomy
dmnd_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/lineage_collapsed.DIAMOND_blastp_hits_nrcontigs_rdrp.txt.gz",
                     delim = "\t", guess_max = 5301)


###### Best bitscore
count = 0
for (q in unique(dmnd_df$qseqid)){
  dmnd_df_q <- dmnd_df[which(dmnd_df$qseqid == q),] %>% arrange(desc(bitscore))
  dmnd_df_q <- dmnd_df_q[1,]
  if (count == 0){
    dmnd_df_bb <- dmnd_df_q
  }else{
    dmnd_df_bb <- bind_rows(dmnd_df_bb, dmnd_df_q)
  }
  count = count + 1
}

write.table(x = dmnd_df_bb,
            file = "/full_path_to/wd/RdRp_scan/analysis/lineage_collapsed.DIAMOND_blastp_hits_nrcontigs_rdrp_BestBitscore.txt",
            sep = "\t", append = F, row.names = F, col.names = T, quote = F)

dmnd_df_bb <- read_delim("/full_path_to/wd/RdRp_scan/analysis/lineage_collapsed.DIAMOND_blastp_hits_nrcontigs_rdrp_BestBitscore.txt",
                      delim = "\t")

colnames(dmnd_df_bb) <- paste0(colnames(dmnd_df_bb), "_bb")

### L set is the selected non-cellular by DIAMOND blastp GN set orfs
L_set_df1_bb <- left_join(L_set_df, dmnd_df_bb, by = c("orf_id"="qseqid_bb"))
L_set_df1_bb$rdrp_scan_vir_group <- str_extract(L_set_df1_bb$target, "[[:alpha:]]+vir[[ia]][[:alpha:]]+")
L_set_df1_bb[which(is.na(L_set_df1_bb$rdrp_scan_vir_group)),]$rdrp_scan_vir_group <- str_extract(L_set_df1_bb[which(is.na(L_set_df1_bb$rdrp_scan_vir_group)),]$target,
                                                                                                 "clstr_[[:digit:]]+")

unique(L_set_df1_bb$target)

temp <- L_set_df1_bb %>% select(rdrp_scan_vir_group, target) %>% distinct() #### check if word extraction was correct

### Removing RdRP-scan hits with no DIAMOND blastp hits
L_set_df2_bb <- L_set_df1_bb[which(!is.na(L_set_df1_bb$pident_bb)),]

write.table(x = L_set_df2_bb,
            file = "/full_path_to/wd/RdRp_scan/analysis/L_df_hmmscan.LCA.BestBitscore_hit.txt",
            sep = "\t", append = F, row.names = F, col.names = T, quote = F)


###### Best evalue

count = 0
for (q in unique(dmnd_df$qseqid)){
  dmnd_df_q <- dmnd_df[which(dmnd_df$qseqid == q),] %>% arrange(evalue)
  dmnd_df_q <- dmnd_df_q[1,]
  if (count == 0){
    dmnd_df_be <- dmnd_df_q
  }else{
    dmnd_df_be <- bind_rows(dmnd_df_be, dmnd_df_q)
  }
  count = count + 1
}

write.table(x = dmnd_df_be,
            file = "/full_path_to/wd/RdRp_scan/analysis/lineage_collapsed.DIAMOND_blastp_hits_nrcontigs_rdrp_BestEvalue.txt",
            sep = "\t", append = F, row.names = F, col.names = T, quote = F)

dmnd_df_be <- read_delim("/full_path_to/wd/RdRp_scan/analysis/lineage_collapsed.DIAMOND_blastp_hits_nrcontigs_rdrp_BestEvalue.txt",
                         delim = "\t")

colnames(dmnd_df_be) <- paste0(colnames(dmnd_df_be), "_be")

### L set is the selected non-cellular by DIAMOND blastp GN set orfs
L_set_df1_be <- left_join(L_set_df, dmnd_df_be, by = c("orf_id"="qseqid_be"))
L_set_df1_be$rdrp_scan_vir_group <- str_extract(L_set_df1_be$target, "[[:alpha:]]+vir[[ia]][[:alpha:]]+")
L_set_df1_be[which(is.na(L_set_df1_be$rdrp_scan_vir_group)),]$rdrp_scan_vir_group <- str_extract(L_set_df1_be[which(is.na(L_set_df1_be$rdrp_scan_vir_group)),]$target,
                                                                                                 "clstr_[[:digit:]]+")

temp <- L_set_df1_be %>% select(rdrp_scan_vir_group, target) %>% distinct() #### check if word extraction was correct


##### Check contigs with no matches
unclassified_df<- dmnd_df[which(dmnd_df$qseqid %in% L_set_df1_be[which(L_set_df1_be$tax_id == 0),]$qseqid),]
unclassified_df_be<- dmnd_df_be[which(dmnd_df_be$qseqid_be %in% L_set_df1_be[which(L_set_df1_be$tax_id == 0),]$qseqid),]
##### Some just don't have taxids but others didn't have hits at all and they were marked as taxid 0 in LCA analysis


### Removing RdRP-scan hits with no DIAMOND blastp hits
L_set_df2_be <- L_set_df1_be[which(!is.na(L_set_df1_be$pident_be)),]


length(unique(L_set_df2_be$qseqid))

write.table(x = L_set_df2_be,
            file = "/full_path_to/wd/RdRp_scan/analysis/L_df_hmmscan.LCA.BestEvalue_hit.txt",
            sep = "\t", append = F, row.names = F, col.names = T, quote = F)
