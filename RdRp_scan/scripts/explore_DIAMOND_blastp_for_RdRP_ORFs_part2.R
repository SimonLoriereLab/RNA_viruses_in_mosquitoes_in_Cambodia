library(tidyverse)
library(seqinr)


L_set_df2_be <- read_delim("/full_path_to/wd/RdRp_scan/analysis/L_df_hmmscan.LCA.BestEvalue_hit.txt", 
                         delim = "\t")



L_set_df2_be$RdRP_scan_match_with_DMND_LCA <- 0
L_set_df2_be[which(str_detect(L_set_df2_be$Collapsed_lineage, L_set_df2_be$rdrp_scan_vir_group)),]$RdRP_scan_match_with_DMND_LCA <- 1
sum(L_set_df2_be$RdRP_scan_match_with_DMND_LCA)


##### Visualize the results
p1 <- ggplot()
p1 <- p1 + geom_point(data = L_set_df2_be, 
                      mapping = aes(x = score_full_seq, 
                                    y = bitscore_be, alpha = RdRP_scan_match_with_DMND_LCA))
print(p1)


p2 <- ggplot()
p2 <- p2 + geom_point(data = L_set_df2_be, 
                      mapping = aes(x = pident_be, 
                                    y = log10(orf_len), size = length_be, 
                                    alpha = RdRP_scan_match_with_DMND_LCA))
p2 <- p2 + facet_wrap(facets = "rdrp_scan_vir_group", ncol = 4)
p2 <- p2 + theme_minimal()
print(p2)



########## Run filtering by size script

R_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/R_df_hmmscan.LCA.BestEvalue_hit.txt", 
                   delim = "\t")

R_orfs <- read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/R_ORFs.fasta", seqtype = "AA")

setdiff(names(R_orfs), R_df$qseqid)


R_df$RdRP_scan_match_with_DMND <- "no match"
R_df[which(str_detect(R_df$Collapsed_lineage, R_df$rdrp_scan_vir_group)),]$RdRP_scan_match_with_DMND <- "LCA"
R_df[which(str_detect(R_df$Collapsed_lineage_be, R_df$rdrp_scan_vir_group)),]$RdRP_scan_match_with_DMND <- "BestHit"

R_df$RdRP_scan_match_with_DMND <- factor(R_df$RdRP_scan_match_with_DMND, 
                                         levels = c("no match", "BestHit", "LCA"), ordered = T)


check_no_match <- R_df[which(R_df$RdRP_scan_match_with_DMND == "no match"),]


unique(R_df$RdRP_scan_match_with_DMND) 

colnames(R_df)

temp <- R_df[which(R_df$RdRP_scan_match_with_DMND == "no match"),] %>% select(qseqid, Collapsed_lineage, Collapsed_lineage_be, rdrp_scan_vir_group)



##### Visualize the results
p1 <- ggplot()
p1 <- p1 + geom_point(data = R_df, 
                      mapping = aes(x = score_full_seq, 
                                    y = bitscore_be, color = RdRP_scan_match_with_DMND), alpha = 0.5)
print(p1)



dir.create(path = "/full_path_to/wd/RdRp_scan/analysis/plots")
p2 <- ggplot()
p2 <- p2 + geom_point(data = R_df, 
                      mapping = aes(x = pident_be, 
                                    y = orf_len/3, 
                                    #size = length_be/3, 
                                    color = RdRP_scan_match_with_DMND), shape = 1)
#p2 <- p2 + ylim(c(min(R_df$orf_len)/3, max(R_df$orf_len)/3))
p2 <- p2 + scale_y_continuous(breaks = c(2000, 4000, 6000))
p2 <- p2 + scale_x_continuous(breaks = c(0, 50, 100))
p2 <- p2 + scale_color_manual(values = c("no match" = "black", "BestHit" = "firebrick3"))

p2 <- p2 + facet_wrap(facets = "rdrp_scan_vir_group", ncol = 6)
p2 <- p2 + theme_minimal(base_size = 8)
p2 <- p2 + theme(legend.position = "bottom")
print(p2)

ggsave(filename = "/full_path_to/wd/RdRp_scan/analysis/plots/selected_orfs_both_methods.pdf", 
       plot = p2, device = "pdf",width = 12, height = 10, units = "cm")




R_df$orf_id_original <- R_df$orf_id


length(unique(R_df$contig_id))  ####### check contigs with multiple orfs

R_df_orfs <- R_df %>% distinct(contig_id, .keep_all = T) %>% select(orf_id)
dup_contigs <- str_remove(setdiff(R_df$orf_id, R_df_orfs$orf_id), "\\_[[:digit:]]+$")
R_df_dup_contigs <- R_df[which(R_df$contig_id %in% dup_contigs),]


####### Transform domain positions according to their ORF positions

R_df_dup_contigs$qstart_be_nt <- NA
R_df_dup_contigs$qend_be_nt <- NA
for (i in 1:nrow(R_df_dup_contigs)){
  if (R_df_dup_contigs[i,]$orf_sense == "FORWARD"){
    R_df_dup_contigs[i,]$qstart_be_nt <- R_df_dup_contigs[i,]$orf_start + 3*(R_df_dup_contigs[i,]$qstart_be - 1)
    R_df_dup_contigs[i,]$qend_be_nt <- R_df_dup_contigs[i,]$orf_start + 3*(R_df_dup_contigs[i,]$qend_be) - 1
  } else {
    R_df_dup_contigs[i,]$qstart_be_nt <- R_df_dup_contigs[i,]$orf_start - 3*(R_df_dup_contigs[i,]$qstart_be - 1)
    R_df_dup_contigs[i,]$qend_be_nt <- R_df_dup_contigs[i,]$orf_start - 3*(R_df_dup_contigs[i,]$qend_be) + 1
  }
}




dir.create("/full_path_to/wd/RdRp_scan/analysis/plots/contigs_with_several_orfs/")

for (cont in sort(unique(R_df_dup_contigs$contig_id))){
  R_df_dup_contig <- R_df_dup_contigs[which(R_df_dup_contigs$contig_id == cont),]
  R_df_dup_contig$orf_id_short <- str_extract(R_df_dup_contig$orf_id, "[[:digit:]]+$")

  ### Plot
  p3 <- ggplot()
  p3 <- p3 + geom_segment(data = R_df_dup_contig, 
                          mapping = aes(x = orf_start, xend = orf_end, 
                                        y = orf_id_short, yend = orf_id_short, color = orf_sense),
                          size = 1)
  p3 <- p3 + geom_segment(data = R_df_dup_contig, 
                          mapping = aes(x = qstart_be_nt, xend = qend_be_nt, 
                                        y = orf_id_short, yend = orf_id_short, alpha = pident_be), 
                          size = 2, color = "firebrick3")
  p3 <- p3 + geom_text(data = R_df_dup_contig, 
                          mapping = aes(x = (qstart_be_nt + qend_be_nt)/2, 
                                        y = orf_id_short, 
                                        label = rdrp_scan_vir_group), nudge_y = 0.5)
  p3 <- p3 + geom_text(data = R_df_dup_contig, 
                       mapping = aes(x = (qstart_be_nt + qend_be_nt)/2, 
                                     y = orf_id_short, 
                                     label = sseqid_be), nudge_y = -0.5)
  p3 <- p3 + scale_color_manual(values = c("FORWARD" = "#f7941d", "REVERSE" = "#2e3192"))
  
  p3 <- p3 + scale_alpha_continuous(range = c(0.3,1))
  p3 <- p3 + ggtitle(label = cont)
  p3 <- p3 + theme_minimal(base_size = 8)
  p3 <- p3 + theme(legend.position = "bottom")
  
  ggsave(filename = paste0("/full_path_to/wd/RdRp_scan/analysis/plots/contigs_with_several_orfs/",cont,".pdf"), 
         plot = p3, device = "pdf",width = 15, height = 10, units = "cm")
  
}

#################################################
######### Dealing with these contigs one by one.
#################################################

### Contig F1621_NODE_24_length_3179_cov_5.991997
R_df_dup_contigs[which(R_df_dup_contigs$contig_id == "F1621_NODE_24_length_3179_cov_5.991997"),]$contig_set
### This contig has not been read mapped yet, it's from contig set N
R_df_dup_contigs[which(R_df_dup_contigs$contig_id == "F1621_NODE_24_length_3179_cov_5.991997"),]$orf_start
(3177-1942)%/%3
(3177-1942)%%3
(3177-1942)/3
### frameshift +2/-1, possibly an assembly mistake
R_df_dup_contigs[which(R_df_dup_contigs$contig_id == "F1621_NODE_24_length_3179_cov_5.991997"),]$orf_end
### checking the best hit and alignment in Geneious
#### the orfs are overlapping by ~34 aa. 
#### Concatenated using the alignment to the reference.

R_df_dup_contigs_curated <- R_df_dup_contigs

R_df_dup_contigs_curated[which(R_df_dup_contigs_curated$contig_id == "F1621_NODE_24_length_3179_cov_5.991997"),]$orf_id <- "F1621_NODE_24_length_3179_cov_5.991997_1concat2"



### F1358_NODE_2309_length_1852_cov_1.580412 
R_df_dup_contigs[which(R_df_dup_contigs$contig_id == "F1358_NODE_2309_length_1852_cov_1.580412"),]$contig_set
R_df_dup_contigs[which(R_df_dup_contigs$contig_id == "F1358_NODE_2309_length_1852_cov_1.580412"),]$orf_start
(1852-956)%/%3
(1852-956)%%3
(1852-956)/3
### frameshift +2/-1, possibly an assembly mistake
### checking the best hit and alignment in Geneious
#### the orfs are overlapping.
#### Concatenated using the alignment to the reference. Removed overlapping regions and an overhang in the beginning of orf 1.
R_df_dup_contigs_curated[which(R_df_dup_contigs_curated$contig_id == "F1358_NODE_2309_length_1852_cov_1.580412"),]$orf_id <- "F1358_NODE_2309_length_1852_cov_1.580412_1concat2"


### F1015_NODE_34_length_4438_cov_44.128451
### Very simple case, where one of the orfs corresponds to a hypothetical protein that is somehow similar with clstr 378 of RdRP scan.
### I will not keep it, because there is another orf that corresponds to RdRp.
R_df_dup_contigs_curated1 <- R_df_dup_contigs_curated[which(!R_df_dup_contigs_curated$orf_id == "F1015_NODE_34_length_4438_cov_44.128451_1"),]

### F1194_NODE_83_length_4344_cov_522.683143
### the same case as above
R_df_dup_contigs_curated2 <- R_df_dup_contigs_curated1[which(!R_df_dup_contigs_curated1$orf_id == "F1194_NODE_83_length_4344_cov_522.683143_1"),]

### F1061_NODE_159_length_4490_cov_42.580383
### the same case as above, but a different best hit, so RdRp is annotated as hypothetical protein 2, which is 94% aa identical to Aedes binegev RdRp
R_df_dup_contigs_curated3 <- R_df_dup_contigs_curated2[which(!R_df_dup_contigs_curated2$orf_id == "F1061_NODE_159_length_4490_cov_42.580383_1"),]



R_df_nondup_contigs <- R_df[which(!R_df$contig_id %in% unique(R_df_dup_contigs_curated3$contig_id)),]


R_df_curated <- bind_rows(R_df_nondup_contigs, R_df_dup_contigs_curated3)




###### Remove possibly non-cellular by LCA or best hit
check3 <- R_df_curated[which(str_detect(R_df_curated$Collapsed_lineage, regex("bact", ignore_case = T))),]
check4 <- R_df_curated[which(str_detect(R_df_curated$Collapsed_lineage, regex("euka", ignore_case = T))),]
check5 <- R_df_curated[which(str_detect(R_df_curated$Collapsed_lineage, regex("archaea", ignore_case = T))),]
check31 <- R_df_curated[which(str_detect(R_df_curated$Collapsed_lineage_be, regex("bact", ignore_case = T))),]
check41 <- R_df_curated[which(str_detect(R_df_curated$Collapsed_lineage_be, regex("euka", ignore_case = T))),]
check51 <- R_df_curated[which(str_detect(R_df_curated$Collapsed_lineage_be, regex("archaea", ignore_case = T))),]

R_df_curated <- R_df_curated[which(!R_df_curated$contig_id == unique(check31$contig_id)),]
R_df_curated <- R_df_curated[which(!R_df_curated$contig_id == unique(check41$contig_id)),]


R_df_curated <- R_df_curated[which(!is.na(R_df_curated$tax_id_be)),]

















p4 <- ggplot()
p4 <- p4 + geom_point(data = R_df_curated, 
                      mapping = aes(x = pident_be, 
                                    y = orf_len/3, 
                                    #size = length_be/3, 
                                    color = RdRP_scan_match_with_DMND), shape = 1)
p4 <- p4 + scale_y_continuous(breaks = c(2000, 4000, 6000))
p4 <- p4 + scale_x_continuous(breaks = c(0, 50, 100))
p4 <- p4 + scale_color_manual(values = c("no match" = "black", "BestHit" = "firebrick3"))

p4 <- p4 + facet_wrap(facets = "rdrp_scan_vir_group", ncol = 5)
p4 <- p4 + theme_minimal(base_size = 8)
p4 <- p4 + theme(legend.position = "bottom")
print(p4)

ggsave(filename = "/full_path_to/wd/RdRp_scan/analysis/plots/selected_orfs_both_methods_curated_multiorf.pdf", 
       plot = p4, device = "pdf",width = 12, height = 10, units = "cm")






length(unique(R_df_curated$orf_id))
length(unique(R_df_curated$contig_id))



p5 <- ggplot()
p5 <- p5 + geom_histogram(data = R_df_curated, 
                        mapping = aes(x = (orf_len/3)), binwidth = 300, 
                        color = "grey", fill = "grey", alpha = 1)
p5 <- p5 + geom_histogram(data = R_df_curated[which(R_df_curated$RdRP_scan_match_with_DMND == "BestHit"),], 
                        mapping = aes(x = (orf_len/3)),  binwidth = 300, 
                        color = "black", fill = "black", alpha = 1)
p5 <- p5 + xlab("Protein sequence length, aa")
p5 <- p5 + theme_classic(base_size = 8)
p5

ggsave(filename = "/full_path_to/wd/RdRp_scan/analysis/plots/selected_orfs_both_methods_curated_multiorf_hist_protlen.pdf",
       plot = p5, device = "pdf", width = 5, height = 5, units = "cm")

p6 <- ggplot()
p6 <- p6 + geom_histogram(data = R_df_curated, 
                          mapping = aes(x = pident_be), binwidth = 5, 
                          color = "grey", fill = "grey", alpha = 1)
p6 <- p6 + geom_histogram(data = R_df_curated[which(R_df_curated$RdRP_scan_match_with_DMND == "BestHit"),], 
                          mapping = aes(x = pident_be),  binwidth = 5, 
                          color = "black", fill = "black", alpha = 1)
p6 <- p6 + xlab("Identity of the best DIAMOND BLASTp hit, % aa")
p6 <- p6 + theme_classic(base_size = 8)
p6

ggsave(filename = "/full_path_to/wd/RdRp_scan/analysis/plots/selected_orfs_both_methods_curated_multiorf_hist_protpident.pdf",
       plot = p6, device = "pdf", width = 5, height = 5, units = "cm")


########## Matching taxonomies more precisely


R_df_curated$RdRP_scan_match_with_DMND <- "no match"
R_df_curated[which(str_detect(R_df_curated$Collapsed_lineage, R_df_curated$rdrp_scan_vir_group)),]$RdRP_scan_match_with_DMND <- "LCA"
R_df_curated[which(str_detect(R_df_curated$Collapsed_lineage_be, R_df_curated$rdrp_scan_vir_group)),]$RdRP_scan_match_with_DMND <- "BestHit"

#### verifying a mismatch versus absence of taxonomic annotation
check_no_match <- R_df_curated[which(R_df_curated$RdRP_scan_match_with_DMND == "no match"),] %>% 
  select(Collapsed_lineage_be, rdrp_scan_vir_group) %>% distinct()

for (i in 1:nrow(check_no_match)){
  temp_check_no_match <- check_no_match[i,]
  rdrp_taxon <- str_extract(temp_check_no_match$rdrp_scan_vir_group, "vir[[:alpha:]]+$")
  print(rdrp_taxon)
  corresponding_be_taxon <- str_extract(temp_check_no_match$Collapsed_lineage_be, paste0("[[:alpha:]]+", rdrp_taxon))
  if (is.na(corresponding_be_taxon)){
    R_df_curated[which(R_df_curated$Collapsed_lineage_be == temp_check_no_match$Collapsed_lineage_be &
                         R_df_curated$rdrp_scan_vir_group == temp_check_no_match$rdrp_scan_vir_group),]$RdRP_scan_match_with_DMND <- "BestHit unclassified species or genus"
    print(R_df_curated[which(R_df_curated$Collapsed_lineage_be == temp_check_no_match$Collapsed_lineage_be &
                         R_df_curated$rdrp_scan_vir_group == temp_check_no_match$rdrp_scan_vir_group),]$RdRP_scan_match_with_DMND)
  } else {
    print(corresponding_be_taxon)
    
    R_df_curated[which(R_df_curated$Collapsed_lineage_be == temp_check_no_match$Collapsed_lineage_be &
                         R_df_curated$rdrp_scan_vir_group == temp_check_no_match$rdrp_scan_vir_group),]$RdRP_scan_match_with_DMND <- "actual mismatch"
    
    print(R_df_curated[which(R_df_curated$Collapsed_lineage_be == temp_check_no_match$Collapsed_lineage_be &
                               R_df_curated$rdrp_scan_vir_group == temp_check_no_match$rdrp_scan_vir_group),]$RdRP_scan_match_with_DMND)
    
  }
  
}


unique(R_df_curated$RdRP_scan_match_with_DMND)

R_df_curated$RdRP_scan_match_with_DMND <- factor(R_df_curated$RdRP_scan_match_with_DMND, 
                                         levels = c("BestHit unclassified species or genus","actual mismatch", "BestHit"), ordered = T)













p51 <- ggplot()
p51 <- p51 + geom_histogram(data = R_df_curated, 
                          mapping = aes(x = (orf_len/3), fill = RdRP_scan_match_with_DMND), 
                          binwidth = 300, position = "stack")
p51 <- p51 + scale_fill_manual(values = c("BestHit unclassified species or genus" = "#fee0d2", 
                                          "actual mismatch" = "#fc9272", 
                                          "BestHit" = "#de2d26"))
p51 <- p51 + xlab("Protein sequence length, aa")
p51 <- p51 + theme_classic(base_size = 8)
p51 <- p51 + theme(legend.position = "none")


ggsave(filename = "/full_path_to/wd/RdRp_scan/analysis/plots/selected_orfs_both_methods_curated_multiorf_hist_protlen_v2.pdf",
       plot = p51, device = "pdf", width = 5, height = 5, units = "cm")





p61 <- ggplot()
p61 <- p61 + geom_histogram(data = R_df_curated, 
                          mapping = aes(x = pident_be, fill = RdRP_scan_match_with_DMND), 
                          binwidth = 5, position = "stack")
p61 <- p61 +scale_fill_manual(values = c("BestHit unclassified species or genus" = "#fee0d2", 
                                          "actual mismatch" = "#fc9272", 
                                          "BestHit" = "#de2d26"))
p61 <- p61 + xlab("Identity of the best DIAMOND BLASTp hit, % aa")
p61 <- p61 + theme_classic(base_size = 8)
p61 <- p61 + theme(legend.position = "none")

ggsave(filename = "/full_path_to/wd/RdRp_scan/analysis/plots/selected_orfs_both_methods_curated_multiorf_hist_protpident_v2.pdf",
       plot = p61, device = "pdf", width = 5, height = 5, units = "cm")

















temp <- R_df_curated[which(R_df_curated$RdRP_scan_match_with_DMND == "BestHit"),]

length(unique(temp$orf_id))
length(unique(temp$contig_id))


temp <- R_df_curated[which(R_df_curated$RdRP_scan_match_with_DMND == "BestHit unclassified species or genus"),]

length(unique(temp$orf_id))
length(unique(temp$contig_id))


temp <- R_df_curated[which(R_df_curated$RdRP_scan_match_with_DMND == "actual mismatch"),] %>% 
  select(orf_id, contig_id, pident_be, Collapsed_lineage_be, rdrp_scan_vir_group) %>% distinct()

length(unique(temp$orf_id))
length(unique(temp$contig_id))










unique(R_df_curated$contig_id)

write.table(x = R_df_curated, file = "/full_path_to/wd/RdRp_scan/analysis/R_orf_df_curated.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)

#### Extract R contigs updated
### I ran extract_R_contigs.R on cluster first.

R_fasta <- read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/R_contigs.fasta",seqtype = "DNA")
R_fasta_non_cellular <- R_fasta[unique(R_df_curated$contig_id)]
write.fasta(sequences = R_fasta_non_cellular, names = names(R_fasta_non_cellular), 
            file.out = "/full_path_to/wd/RdRp_scan/analysis/R_contigs_non_cellular.fasta")

