library(tidyverse)
library(seqinr)
library(ggtree)
library(treeio)


load("/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree_tips_clean.ggtree_midpoint_by_figtree_rectangular_v3_data.rdata")


tip_df_ord <- ggt6$data[which(ggt6$data$isTip == T),]



R_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/R_orf_df_curated_with_chimeric_contig_id.txt", 
                   delim = "\t")
##### get rid of duplicate hits
temp <- R_df[which(str_detect(R_df$orf_id, "F1621_NODE_24_length_3179_cov_5.991997")),]
nrow(temp) == 1
R_df <- R_df[which(!(R_df$orf_id == "F1621_NODE_24_length_3179_cov_5.991997_1concat2" & R_df$evalue_be == max(temp$evalue_be))),] ### concatenated orfs
temp <- R_df[which(str_detect(R_df$orf_id, "F1621_NODE_24_length_3179_cov_5.991997")),]
nrow(temp) == 1

temp <- R_df[which(str_detect(R_df$orf_id, "F1358_NODE_2309_length_1852_cov_1.580412")),]
nrow(temp) == 1
R_df <- R_df[which(!(R_df$orf_id == "F1358_NODE_2309_length_1852_cov_1.580412_1concat2" & is.na(R_df$species))),] ### concatenated orfs
temp <- R_df[which(str_detect(R_df$orf_id, "F1358_NODE_2309_length_1852_cov_1.580412")),]
nrow(temp) == 1

R_df_reduced <- R_df %>% select(orf_id, 
                                contig_id, 
                                contig_set, 
                                target, 
                                orf_start, orf_end, orf_sense, orf_len, 
                                sseqid_be, slen_be, pident_be, length_be)

R_BE_tax_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/read_mapping/best_hit_taxonomy.txt", 
                          delim = "\t")

R_BE_tax_df_reduced <- R_BE_tax_df %>% select(contig_id,
                                              rdrp_scan_vir_group, 
                                              BestHit_phylum = phylum, 
                                              BestHit_class = class, 
                                              BestHit_order = order,
                                              BestHit_family = family,
                                              BestHit_genus = genus,
                                              BestHit_species = species,
                                              BestHit_pident = pident,
                                              suggested_vir_sp)

R_df_reduced1 <- left_join(R_df_reduced, R_BE_tax_df_reduced, by = "contig_id")



##### Import information from phylogenies after RdRp search (DIAMOND blastx) first analysis
df_previous_analysis <- read_delim("/full_path_to/wd/RdRp_search/analysis/all_families_and_genus_classified/phylo_classified_contigs_and_consensus_w_metadata_correctedAdableaps.txt",
           delim = "\t")


all_df <- left_join(R_df_reduced1, df_previous_analysis, by = c("contig_id" = "contig_ID"),)

write.table(x = all_df, file = "/full_path_to/wd/RdRp_scan/analysis/R_set_merged_hit_data_from_2_analyses.txt", 
            append = F,quote = F,sep = "\t", row.names = F, col.names = T)

### Flaviviridae
Flaviviridae_df <- all_df[which(all_df$BestHit_family == "Flaviviridae"), ]
## only keep original contigs, not consensuses from mapping
Flaviviridae_df <- Flaviviridae_df[which(Flaviviridae_df$original_sample_ID == Flaviviridae_df$read_map_consensus_sample_ID | 
                                           is.na(Flaviviridae_df$read_map_consensus_sample_ID)),]


################################### cISF_clade
## There is a clade of cISFs that I have already put in the phylogeny before, so let's extract these.
cISF_clade_df <- Flaviviridae_df[which(!is.na(Flaviviridae_df$read_map_consensus_sample_ID)),]
cISF_clade_df$followup_phylo_group <- "cISF_clade"
## let's annotate it as a follow-up phylogeny clade
# followup_phylo_df <- cISF_clade_df

dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/")
write.table(x = cISF_clade_df, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/cISF_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)




################################### third_ISF_clade
## Compared to previous analysis I found one more Flaviviridae, potentially Flavivirus
third_ISF_clade_df <- Flaviviridae_df[which(is.na(Flaviviridae_df$read_map_consensus_sample_ID)),]
third_ISF_clade_df$contig_id
#### it is very distant - 26% aa identity based on blastp
third_ISF_clade_df$BestHit_pident
third_ISF_clade_df$pident_be
third_ISF_clade_df$length_be
(third_ISF_clade_df$orf_len+1)/3
third_ISF_clade_df$slen_be
third_ISF_clade_df$BestHit_species

Flaviviridae_df$comment_v3 <- NA

third_ISF_clade_df$comment_v3 <- "New very distant Flaviviridae. Assigned a new name. On RdRp scan phylogeny close to a clade of large insect flaviviruses, but not previously described cISFs and dISFs. Mosquito and tick flaviviruses are also in that clade."
third_ISF_clade_df$suggested_vir_sp <- "Mummution_virus"
third_ISF_clade_df$followup_phylo_group <- "third_ISF_clade"


# third_ISF_clade <- tip_df_ord[tip_df_ord[which(tip_df_ord$merged_vir_sp == "Lampyris_noctiluca_flavivirus_1"),]$node : 
#              tip_df_ord[which(tip_df_ord$merged_vir_sp == "Tacheng_tick_virus_8"),]$node, ]



third_ISF_clade <- tip_df_ord[tip_df_ord[which(tip_df_ord$merged_vir_sp == "Lampyris_noctiluca_flavivirus_1"),]$node : 
                                tip_df_ord[which(tip_df_ord$merged_vir_sp == "Xinzhou_spider_virus_3"),]$node, ]


### modifying the reference sequence accession to the newer sequence with full poluyprotein
third_ISF_clade[which(third_ISF_clade$tip_labels == "Kitrinoviricota_Amarillovirales_BAM99939.1_hypothetical_polyprotein"),]$tip_labels <- "BAM78287"


third_ISF_clade_df$ref_suggestions <- paste0(third_ISF_clade[which(third_ISF_clade$status == "reference"),]$tip_labels, collapse = ", ")



# followup_phylo_df <- bind_rows(followup_phylo_df, third_ISF_clade_df)
# 
# followup_phylo_df %>% select(pident_be, suggested_vir_sp, contig_id)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/third_ISF_clade/")
write.table(x = third_ISF_clade_df, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/third_ISF_clade/third_ISF_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)


### Phenuiviridae
Phenuiviridae_df <- all_df[which(all_df$BestHit_family == "Phenuiviridae"), ]
## only keep original contigs, not consensuses from mapping
Phenuiviridae_df <- Phenuiviridae_df[which(Phenuiviridae_df$original_sample_ID == Phenuiviridae_df$read_map_consensus_sample_ID | 
                                           is.na(Phenuiviridae_df$read_map_consensus_sample_ID)),]


################################### Phasi_clade
## There is a clade of Phasis that I have already put in the phylogeny before, so let's extract these.
Phasi_clade_df <- Phenuiviridae_df[which(!is.na(Phenuiviridae_df$read_map_consensus_sample_ID)),]
Phasi_clade_df$followup_phylo_group <- "Phasi_clade"
## let's annotate it as a follow-up phylogeny clade
# followup_phylo_df <- Phasi_clade_df

dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Phasi_clade/")
write.table(x = Phasi_clade_df, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasi_clade/Phasi_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)









### Orthomyxoviridae
Orthomyxoviridae_df1 <- all_df[which(all_df$BestHit_family == "Orthomyxoviridae"), ]
Orthomyxoviridae_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "F1019_NODE_753_length_2446_cov_43.823923_1")),]$node:
                                      tip_df_ord[which(str_detect(tip_df_ord$label, "NOT_ASSIGNED_Henan_2_orthomyxo_like_virus")),]$node, ]

Orthomyxoviridae_clade_df <- all_df[which(all_df$orf_id %in% 
                                            Orthomyxoviridae_clade[which(Orthomyxoviridae_clade$status == "novel"),]$tip_labels),] %>% 
  discard(~all(is.na(.)))

setdiff(Orthomyxoviridae_clade_df$orf_id, Orthomyxoviridae_df1$orf_id)
setdiff(Orthomyxoviridae_df1$orf_id, Orthomyxoviridae_clade_df$orf_id)
union(Orthomyxoviridae_df1$orf_id, Orthomyxoviridae_clade_df$orf_id)


Orthomyxoviridae_df <- bind_rows(Orthomyxoviridae_clade_df, Orthomyxoviridae_df1) %>% distinct()

## only keep original contigs, not consensuses from mapping
Orthomyxoviridae_df <- Orthomyxoviridae_df[which(Orthomyxoviridae_df$original_sample_ID == Orthomyxoviridae_df$read_map_consensus_sample_ID | 
                                             is.na(Orthomyxoviridae_df$read_map_consensus_sample_ID)),]

################################### Quaranja_clade
Quaranja_clade_df <- Orthomyxoviridae_df[which(!Orthomyxoviridae_df$contig_id == "F1508_NODE_144_length_2236_cov_314.370931"),] %>% discard(~all(is.na(.)))
Quaranja_clade_df$followup_phylo_group <- "Quaranja_clade"

### Orthomyxoviridae sp is actually 89% aa identical to Wuhan Mosquito Virus 4, discovered in 2016, 
### However, to avoid further confusion (Orthomyxoviridae sp might be not unique, Wuhan Mosquito Virus 4 - not the closest)
### Give the virus a new unique name
### even though it should probably be the same species with the two viruses mentioned above 
Quaranja_clade_df[which(Quaranja_clade_df$contig_id == "F1413_NODE_102_length_2466_cov_759.949813"),]$suggested_vir_sp <- "Pageromars_virus"


### Even though F1398_NODE_269_length_2478_cov_48.295089 is 85% (aa) close to either Guadeloupe mosquito quaranja-like virus 1 or Palmetto orthomyxo-like virus
### It seems to form a separate phylogenetic group and to avoid confusion, give it a novel unique name.
Quaranja_clade_df[which(Quaranja_clade_df$contig_id == "F1398_NODE_269_length_2478_cov_48.295089"),]$suggested_vir_sp <- "Panessemat_virus"

Quaranja_clade_df$suggested_vir_sp
Quaranja_clade_df[which(Quaranja_clade_df$contig_id == "F1019_NODE_753_length_2446_cov_43.823923"),]$suggested_vir_sp <- "Spoisillett_virus"
Quaranja_clade_df$suggested_vir_sp


Quaranja_clade <- tip_df_ord[tip_df_ord[which(tip_df_ord$merged_vir_sp == "Henan_2_orthomyxo_like_virus"),]$node, ]
Quaranja_clade_df$ref_suggestions <- paste0(Quaranja_clade[which(Quaranja_clade$status == "reference"),]$tip_labels, collapse = ", ")






dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/")
write.table(x = Quaranja_clade_df, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/Quaranja_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)









################################### Thogoto_clade
Thogoto_clade_df <- Orthomyxoviridae_df[which(Orthomyxoviridae_df$contig_id == "F1508_NODE_144_length_2236_cov_314.370931"),] %>% discard(~all(is.na(.)))
Thogoto_clade_df$followup_phylo_group <- "Thogoto_clade"
### Naming the on more divergent sequence
Thogoto_clade_df[which(Thogoto_clade_df$contig_id == "F1508_NODE_144_length_2236_cov_314.370931"),]$suggested_vir_sp <- "Norayellet_virus"

Thogoto_clade <- tip_df_ord[(tip_df_ord[which(tip_df_ord$merged_vir_sp == "F1508_NODE_144_length_2236_cov_314.370931"),]$node-1), ]
Thogoto_clade_df$ref_suggestions <- paste0(Thogoto_clade[which(Thogoto_clade$status == "reference"),]$tip_labels, collapse = ", ")





dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Thogoto_clade/")
write.table(x = Thogoto_clade_df, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Thogoto_clade/Thogoto_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)














### Sobelivirales
Sobelivirales_df1 <- all_df[which(all_df$BestHit_order == "Sobelivirales"), ]
Sobelivirales_df2 <- all_df[which(all_df$suggested_vir_sp == "Culex_associated_Luteo_like_virus"), ] ## clusters within one of the Sobeli clades
Sobelivirales_df <- bind_rows(Sobelivirales_df1, Sobelivirales_df2)
## only keep original contigs, not consensuses from mapping
Sobelivirales_df <- Sobelivirales_df[which(Sobelivirales_df$original_sample_ID == Sobelivirales_df$read_map_consensus_sample_ID | 
                                                   is.na(Sobelivirales_df$read_map_consensus_sample_ID)),]


################################### Sobeli_clade
Sobeli_clade_df <- Sobelivirales_df[which(!Sobelivirales_df$contig_id == "F1508_NODE_144_length_2236_cov_314.370931"),] %>% discard(~all(is.na(.)))
Sobeli_clade_df$followup_phylo_group <- "Sobeli_clade"

Sobeli_clade_df$suggested_vir_sp


Sobeli_clade_df[which(Sobeli_clade_df$contig_id == "F1127_NODE_300_length_3117_cov_993.590790"),]$suggested_vir_sp <- "Peninsere_virus"
Sobeli_clade_df[which(Sobeli_clade_df$contig_id == "F1502_NODE_85_length_2916_cov_12.682279"),]$suggested_vir_sp <- "Perkeld_virus"
Sobeli_clade_df[which(Sobeli_clade_df$contig_id == "F1032_NODE_17_length_3181_cov_1736.585413"),]$suggested_vir_sp <- "Pirounced_virus"
Sobeli_clade_df[which(Sobeli_clade_df$contig_id == "F1276_NODE_48_length_2574_cov_4.979754"),]$suggested_vir_sp <- "Poterfers_virus"






Sobeli_clade <- tip_df_ord[ tip_df_ord[which(tip_df_ord$merged_vir_sp == "Hubei_sobemo_like_virus_41_subsp2"),]$node : 
                              tip_df_ord[which(tip_df_ord$merged_vir_sp == "Hubei_myriapoda_virus_9"),]$node, ]
Sobeli_clade_df$ref_suggestions <- paste0(Sobeli_clade[which(Sobeli_clade$status == "reference"),]$tip_labels, collapse = ", ")



length(Sobeli_clade[which(Sobeli_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/")
write.table(x = Sobeli_clade_df, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/Sobeli_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)









################################### Monjiviricetes_clade
Monjiviricetes_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "F1097_NODE_1_length_12217_cov_2063.855123_4")),]$node:
                                      tip_df_ord[which(str_detect(tip_df_ord$label, "6UEB\\|Chain_A\\|Large_structural_protein\\|Rabies_virus_")),]$node, ]



Monjiviricetes_clade_df <- all_df[which(all_df$orf_id %in% 
                                                      Monjiviricetes_clade[which(Monjiviricetes_clade$status == "novel"),]$tip_labels),] %>% 
  discard(~all(is.na(.)))


## only keep original contigs, not consensuses from mapping
Monjiviricetes_clade_df <- Monjiviricetes_clade_df[which(Monjiviricetes_clade_df$original_sample_ID == Monjiviricetes_clade_df$read_map_consensus_sample_ID | 
                                                 is.na(Monjiviricetes_clade_df$read_map_consensus_sample_ID)),]



sort(Monjiviricetes_clade[which(Monjiviricetes_clade$status == "novel"),]$tip_labels)


Monjiviricetes_clade_df$ref_suggestions <- paste0(Monjiviricetes_clade[which(Monjiviricetes_clade$status == "reference"),]$tip_labels, collapse = ", ")
Monjiviricetes_clade_df$followup_phylo_group <- "Monjiviricetes_clade"

Monjiviricetes_clade_df$suggested_vir_sp
Monjiviricetes_clade_df[which(Monjiviricetes_clade_df$contig_id == "F1543_NODE_2_length_10836_cov_97.382989"),]$suggested_vir_sp <- "Poterfers_virus"
Monjiviricetes_clade_df[which(Monjiviricetes_clade_df$contig_id == "F1097_NODE_1_length_12217_cov_2063.855123"),]$suggested_vir_sp <- "Preculath_virus"
Monjiviricetes_clade_df[which(Monjiviricetes_clade_df$contig_id == "F1032_NODE_1_length_13213_cov_14.139915"),]$suggested_vir_sp <- "Prepisiers_virus"
Monjiviricetes_clade_df[which(Monjiviricetes_clade_df$contig_id == "F1107_NODE_4_length_4693_cov_19.984260"),]$suggested_vir_sp <- "Raistely_virus"
Monjiviricetes_clade_df[which(Monjiviricetes_clade_df$contig_id == "F1187_NODE_3_length_11590_cov_460.255743"),]$suggested_vir_sp <- "Regreagly_virus"
Monjiviricetes_clade_df[which(Monjiviricetes_clade_df$contig_id == "F1535_NODE_4_length_12152_cov_28.768951"),]$suggested_vir_sp <- "Reloonsia_virus"
Monjiviricetes_clade_df[which(Monjiviricetes_clade_df$contig_id == "F1011_NODE_1_length_12910_cov_159.228082"),]$suggested_vir_sp <- "Retriended_virus"
sort(Monjiviricetes_clade_df$suggested_vir_sp)

### further discarding contigs that are too short, partial, or potential EVEs
Monjiviricetes_drop_contigs <- c("F1632_NODE_476_length_2529_cov_2.204931", "F1632_NODE_737_length_2202_cov_1.634839", 
                                    "F1469_NODE_1254_length_1943_cov_3.096398", "F1198_NODE_97_length_1341_cov_56.915241")
Monjiviricetes_clade_df_drop <- Monjiviricetes_clade_df[which(Monjiviricetes_clade_df$contig_id %in% Monjiviricetes_drop_contigs),]
Monjiviricetes_clade_df_keep <- Monjiviricetes_clade_df[which(!Monjiviricetes_clade_df$contig_id %in% Monjiviricetes_drop_contigs),]




length(Monjiviricetes_clade[which(Monjiviricetes_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/")
write.table(x = Monjiviricetes_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/Monjiviricetes_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Monjiviricetes_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/Monjiviricetes_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)










################################### Peribunyaviridae_clades
##### This includes phenui clade that has already been separately processed, extract other clades


Peribunyaviridae_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "Schistocephalus_solidus_bunya_like_virus")),]$node:
                                      tip_df_ord[which(str_detect(tip_df_ord$label, "Hubei_diptera_virus_3")),]$node, ]





Peribunyaviridae_clade_df <- all_df[which(all_df$orf_id %in% 
                                          Peribunyaviridae_clade[which(Peribunyaviridae_clade$status == "novel"),]$tip_labels),] %>% 
  discard(~all(is.na(.)))


intersect(Phasi_clade_df$orf_id, Peribunyaviridae_clade_df$orf_id)


## only keep original contigs, not consensuses from mapping
Peribunyaviridae_clade_df <- Peribunyaviridae_clade_df[which(Peribunyaviridae_clade_df$original_sample_ID == Peribunyaviridae_clade_df$read_map_consensus_sample_ID | 
                                                           is.na(Peribunyaviridae_clade_df$read_map_consensus_sample_ID)),]



sort(Peribunyaviridae_clade[which(Peribunyaviridae_clade$status == "novel"),]$tip_labels)


Peribunyaviridae_clade_df$ref_suggestions <- paste0(Peribunyaviridae_clade[which(Peribunyaviridae_clade$status == "reference"),]$tip_labels, collapse = ", ")
Peribunyaviridae_clade_df$followup_phylo_group <- "Peribunyaviridae_clade"

Peribunyaviridae_clade_df$suggested_vir_sp
Peribunyaviridae_clade_df[which(Peribunyaviridae_clade_df$contig_id == "F1543_NODE_27_length_5990_cov_6.343724"),]$suggested_vir_sp <- "Snarittled_virus"
Peribunyaviridae_clade_df[which(Peribunyaviridae_clade_df$contig_id == "F1011_NODE_2_length_5406_cov_31.927864"),]$suggested_vir_sp <- "Snowerces_virus"
Peribunyaviridae_clade_df[which(Peribunyaviridae_clade_df$contig_id == "F1422_NODE_7_length_6690_cov_40.209495"),]$suggested_vir_sp <- "Snownnecut_virus"
Peribunyaviridae_clade_df[which(Peribunyaviridae_clade_df$contig_id == "F1486_NODE_8_length_6793_cov_9.764619"),]$suggested_vir_sp <- "Stevained_virus"
Peribunyaviridae_clade_df[which(Peribunyaviridae_clade_df$contig_id == "F1628_NODE_3_length_6523_cov_75.503556"),]$suggested_vir_sp <- "Struded_virus"
Peribunyaviridae_clade_df[which(Peribunyaviridae_clade_df$contig_id == "F1472_NODE_1_length_7201_cov_146.442065"),]$suggested_vir_sp <- "Substaily_virus"
Peribunyaviridae_clade_df[which(Peribunyaviridae_clade_df$contig_id == "F1600_NODE_2_length_5756_cov_13.689177"),]$suggested_vir_sp <- "Sweadeo_virus"
Peribunyaviridae_clade_df[which(Peribunyaviridae_clade_df$contig_id == "F1628_NODE_1_length_7594_cov_95.775700"),]$suggested_vir_sp <- "Throwbor_virus"

sort(Peribunyaviridae_clade_df$suggested_vir_sp)

### further discarding contigs that are too short, partial, or potential EVEs
Peribunyaviridae_drop_contigs <- c("F1172_NODE_144_length_2542_cov_1.231604", 
                                   "F1621_NODE_1440_length_951_cov_4.206473", 
                                   "F1621_NODE_24_length_3179_cov_5.991997")

Peribunyaviridae_clade_df_drop <- Peribunyaviridae_clade_df[which(Peribunyaviridae_clade_df$contig_id %in% Peribunyaviridae_drop_contigs),]
Peribunyaviridae_clade_df_keep <- Peribunyaviridae_clade_df[which(!Peribunyaviridae_clade_df$contig_id %in% Peribunyaviridae_drop_contigs),]

sort(Peribunyaviridae_clade_df_keep$suggested_vir_sp)


length(Peribunyaviridae_clade[which(Peribunyaviridae_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Peribunyaviridae_clade/")
write.table(x = Peribunyaviridae_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Peribunyaviridae_clade/Peribunyaviridae_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Peribunyaviridae_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Peribunyaviridae_clade/Peribunyaviridae_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)




####################################################################################

Phasmaviridae_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "Pythium_polare_bunya_like_RNA_virus_1")),]$node:
                                     tip_df_ord[which(str_detect(tip_df_ord$label, "Negarnaviricota_Ellioviricetes_AIN37024.1_polymerase__Mojui_dos_Campos_virus_")),]$node, ]




Phasmaviridae_clade_df <- all_df[which(all_df$orf_id %in% 
                                            Phasmaviridae_clade[which(Phasmaviridae_clade$status == "novel"),]$tip_labels),] %>% 
  discard(~all(is.na(.)))


## only keep original contigs, not consensuses from mapping
Phasmaviridae_clade_df <- Phasmaviridae_clade_df[which(Phasmaviridae_clade_df$original_sample_ID == Phasmaviridae_clade_df$read_map_consensus_sample_ID | 
                                                               is.na(Phasmaviridae_clade_df$read_map_consensus_sample_ID)),]



sort(Phasmaviridae_clade[which(Phasmaviridae_clade$status == "novel"),]$tip_labels)


Phasmaviridae_clade_df$ref_suggestions <- paste0(Phasmaviridae_clade[which(Phasmaviridae_clade$status == "reference"),]$tip_labels, collapse = ", ")
Phasmaviridae_clade_df$followup_phylo_group <- "Phasmaviridae_clade"

Phasmaviridae_clade_df$suggested_vir_sp
# Phasmaviridae_clade_df[which(Phasmaviridae_clade_df$contig_id == ""),]$suggested_vir_sp <- "_virus"


sort(Phasmaviridae_clade_df$suggested_vir_sp)

### further discarding contigs that are too short, partial, or potential EVEs
Phasmaviridae_drop_contigs <- c("")

Phasmaviridae_clade_df_drop <- Phasmaviridae_clade_df[which(Phasmaviridae_clade_df$contig_id %in% Phasmaviridae_drop_contigs),]
Phasmaviridae_clade_df_keep <- Phasmaviridae_clade_df[which(!Phasmaviridae_clade_df$contig_id %in% Phasmaviridae_drop_contigs),]

sort(Phasmaviridae_clade_df_keep$suggested_vir_sp)


length(Phasmaviridae_clade[which(Phasmaviridae_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/")
write.table(x = Phasmaviridae_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/Phasmaviridae_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Phasmaviridae_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/Phasmaviridae_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)








################################### Iflaviridae_clade
Iflaviridae_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "Pisuviricota_Picornavirales_AVK80197.1_polyprotein__Sogatella_furcifera_honeydew_virus_")),]$node:
                                   tip_df_ord[which(str_detect(tip_df_ord$label, "Pisuviricota_Picornavirales_AWK77886.1_polyprotein__Victoria_bee_virus_1_")),]$node, ]




Iflaviridae_clade_df <- all_df[which(all_df$orf_id %in% 
                                            Iflaviridae_clade[which(Iflaviridae_clade$status == "novel"),]$tip_labels),] %>% 
  discard(~all(is.na(.)))



## only keep original contigs, not consensuses from mapping
Iflaviridae_clade_df <- Iflaviridae_clade_df[which(Iflaviridae_clade_df$original_sample_ID == Iflaviridae_clade_df$read_map_consensus_sample_ID | 
                                                               is.na(Iflaviridae_clade_df$read_map_consensus_sample_ID)),]



sort(Iflaviridae_clade[which(Iflaviridae_clade$status == "novel"),]$tip_labels)


Iflaviridae_clade_df$ref_suggestions <- paste0(Iflaviridae_clade[which(Iflaviridae_clade$status == "reference"),]$tip_labels, collapse = ", ")
Iflaviridae_clade_df$followup_phylo_group <- "Iflaviridae_clade"

Iflaviridae_clade_df$suggested_vir_sp
# Iflaviridae_clade_df[which(Iflaviridae_clade_df$contig_id == ""),]$suggested_vir_sp <- ""


sort(Iflaviridae_clade_df$suggested_vir_sp)

### further discarding contigs that are too short, partial, or potential EVEs
Iflaviridae_drop_contigs <- c()

Iflaviridae_clade_df_drop <- Iflaviridae_clade_df[which(Iflaviridae_clade_df$contig_id %in% Iflaviridae_drop_contigs),]
Iflaviridae_clade_df_keep <- Iflaviridae_clade_df[which(!Iflaviridae_clade_df$contig_id %in% Iflaviridae_drop_contigs),]

sort(Iflaviridae_clade_df_keep$suggested_vir_sp)


length(Iflaviridae_clade[which(Iflaviridae_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Iflaviridae_clade/")
write.table(x = Iflaviridae_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Iflaviridae_clade/Iflaviridae_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Iflaviridae_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Iflaviridae_clade/Iflaviridae_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)




################################### Dicistroviridae_clade
Dicistroviridae_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "Pisuviricota_Picornavirales_ya20_JAAOEH010000223_1_JAAOEH010000223.1_43275733_OV.13_NODE_222_truseq_orf.1333")),]$node:
                                   tip_df_ord[which(str_detect(tip_df_ord$label, "Pisuviricota_Picornavirales_ya20_JAAOEH010003388_1_JAAOEH010003388.1_21123149_OV.60_NODE_3374_truseq_orf.20325")),]$node, ]




Dicistroviridae_clade_df <- all_df[which(all_df$orf_id %in% 
                                       Dicistroviridae_clade[which(Dicistroviridae_clade$status == "novel"),]$tip_labels),] %>% 
  discard(~all(is.na(.)))


## only keep original contigs, not consensuses from mapping
Dicistroviridae_clade_df <- Dicistroviridae_clade_df[which(Dicistroviridae_clade_df$original_sample_ID == Dicistroviridae_clade_df$read_map_consensus_sample_ID | 
                                                     is.na(Dicistroviridae_clade_df$read_map_consensus_sample_ID)),]



sort(Dicistroviridae_clade[which(Dicistroviridae_clade$status == "novel"),]$tip_labels)


Dicistroviridae_clade_df$ref_suggestions <- paste0(Dicistroviridae_clade[which(Dicistroviridae_clade$status == "reference"),]$tip_labels, collapse = ", ")
Dicistroviridae_clade_df$followup_phylo_group <- "Dicistroviridae_clade"

Dicistroviridae_clade_df$suggested_vir_sp
Dicistroviridae_clade_df[which(Dicistroviridae_clade_df$contig_id == "F1429_NODE_4_length_9909_cov_248.781510"),]$suggested_vir_sp <- "Unictinger_virus"
Dicistroviridae_clade_df[which(Dicistroviridae_clade_df$contig_id == "F1539_NODE_1_length_9238_cov_54.876511"),]$suggested_vir_sp <- "Untradvers_virus"

sort(Dicistroviridae_clade_df$suggested_vir_sp)

### further discarding contigs that are too short, partial, or potential EVEs
Dicistroviridae_drop_contigs <- c()

Dicistroviridae_clade_df_drop <- Dicistroviridae_clade_df[which(Dicistroviridae_clade_df$contig_id %in% Dicistroviridae_drop_contigs),]
Dicistroviridae_clade_df_keep <- Dicistroviridae_clade_df[which(!Dicistroviridae_clade_df$contig_id %in% Dicistroviridae_drop_contigs),]

sort(Dicistroviridae_clade_df_keep$suggested_vir_sp)


length(Dicistroviridae_clade[which(Dicistroviridae_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Dicistroviridae_clade/")
write.table(x = Dicistroviridae_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Dicistroviridae_clade/Dicistroviridae_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Dicistroviridae_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Dicistroviridae_clade/Dicistroviridae_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)






################################### Caliciviridae_clade
Caliciviridae_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "Pisuviricota_Picornavirales_AVM87220.1_polyprotein__Wenling_yellow_goosefish_calicivirus_")),]$node:
                                       tip_df_ord[which(str_detect(tip_df_ord$label, "Pisuviricota_Picornavirales_Xinjiang_6_picorna_like_virus_12")),]$node, ]




Caliciviridae_clade_df <- all_df[which(all_df$orf_id %in% 
                                           Caliciviridae_clade[which(Caliciviridae_clade$status == "novel"),]$tip_labels),] 



## only keep original contigs, not consensuses from mapping
Caliciviridae_clade_df1 <- Caliciviridae_clade_df[which(Caliciviridae_clade_df$original_sample_ID == Caliciviridae_clade_df$read_map_consensus_sample_ID | 
                                                             is.na(Caliciviridae_clade_df$read_map_consensus_sample_ID)),] %>% 
  discard(~all(is.na(.)))




sort(Caliciviridae_clade[which(Caliciviridae_clade$status == "novel"),]$tip_labels)


Caliciviridae_clade_df$ref_suggestions <- paste0(Caliciviridae_clade[which(Caliciviridae_clade$status == "reference"),]$tip_labels, collapse = ", ")
Caliciviridae_clade_df$followup_phylo_group <- "Caliciviridae_clade"

Caliciviridae_clade_df$suggested_vir_sp
# Caliciviridae_clade_df[which(Caliciviridae_clade_df$contig_id == ""),]$suggested_vir_sp <- "_virus"


sort(Caliciviridae_clade_df$suggested_vir_sp)

### further discarding contigs that are too short, partial, or potential EVEs
Caliciviridae_drop_contigs <- c()

Caliciviridae_clade_df_drop <- Caliciviridae_clade_df[which(Caliciviridae_clade_df$contig_id %in% Caliciviridae_drop_contigs),]
Caliciviridae_clade_df_keep <- Caliciviridae_clade_df[which(!Caliciviridae_clade_df$contig_id %in% Caliciviridae_drop_contigs),]

sort(Caliciviridae_clade_df_keep$suggested_vir_sp)


length(Caliciviridae_clade[which(Caliciviridae_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Caliciviridae_clade/")
write.table(x = Caliciviridae_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Caliciviridae_clade/Caliciviridae_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Caliciviridae_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Caliciviridae_clade/Caliciviridae_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)










################################### Picorna_like_clade
Picorna_like_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "Pisuviricota_Picornavirales_YP_009345020.1_hypothetical_protein__Zhejiang_mosquito_virus_1_")),]$node:
                                     tip_df_ord[which(str_detect(tip_df_ord$label, "Pisuviricota_Picornavirales_YP_009337062.1_hypothetical_protein__Hubei_picornalike_virus_61_")),]$node, ]




Picorna_like_clade_df <- all_df[which(all_df$orf_id %in% 
                                         Picorna_like_clade[which(Picorna_like_clade$status == "novel"),]$tip_labels),] 


## only keep original contigs, not consensuses from mapping
Picorna_like_clade_df1 <- Picorna_like_clade_df[which(Picorna_like_clade_df$original_sample_ID == Picorna_like_clade_df$read_map_consensus_sample_ID | 
                                                          is.na(Picorna_like_clade_df$read_map_consensus_sample_ID)),] %>% 
  discard(~all(is.na(.)))




sort(Picorna_like_clade[which(Picorna_like_clade$status == "novel"),]$tip_labels)


Picorna_like_clade_df$ref_suggestions <- paste0(Picorna_like_clade[which(Picorna_like_clade$status == "reference"),]$tip_labels, collapse = ", ")
Picorna_like_clade_df$followup_phylo_group <- "Picorna_like_clade"

Picorna_like_clade_df$suggested_vir_sp
# Picorna_like_clade_df[which(Picorna_like_clade_df$contig_id == ""),]$suggested_vir_sp <- "_virus"


sort(Picorna_like_clade_df$suggested_vir_sp)

### further discarding contigs that are too short, partial, or potential EVEs
Picorna_like_drop_contigs <- c()

Picorna_like_clade_df_drop <- Picorna_like_clade_df[which(Picorna_like_clade_df$contig_id %in% Picorna_like_drop_contigs),]
Picorna_like_clade_df_keep <- Picorna_like_clade_df[which(!Picorna_like_clade_df$contig_id %in% Picorna_like_drop_contigs),]

sort(Picorna_like_clade_df_keep$suggested_vir_sp)


length(Picorna_like_clade[which(Picorna_like_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Picorna_like_clade/")
write.table(x = Picorna_like_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Picorna_like_clade/Picorna_like_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Picorna_like_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Picorna_like_clade/Picorna_like_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)



################################### Mesoniviridae_clade
Mesoniviridae_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "Pisuviricota_Nidovirales_ATV90887.1_putative_replicase_polyprotein__Yichang_virus_")),]$node:
                                    tip_df_ord[which(str_detect(tip_df_ord$label, "Pisuviricota_Nidovirales_QNM37795.1_putative_replicase__Frankliniella_occidentalis_associated_mesonivirus_1_")),]$node, ]




Mesoniviridae_clade_df <- all_df[which(all_df$orf_id %in% 
                                        Mesoniviridae_clade[which(Mesoniviridae_clade$status == "novel"),]$tip_labels),] 



## only keep original contigs, not consensuses from mapping
Mesoniviridae_clade_df1 <- Mesoniviridae_clade_df[which(Mesoniviridae_clade_df$original_sample_ID == Mesoniviridae_clade_df$read_map_consensus_sample_ID | 
                                                        is.na(Mesoniviridae_clade_df$read_map_consensus_sample_ID)),] %>% 
  discard(~all(is.na(.)))




sort(Mesoniviridae_clade[which(Mesoniviridae_clade$status == "novel"),]$tip_labels)


Mesoniviridae_clade_df$ref_suggestions <- paste0(Mesoniviridae_clade[which(Mesoniviridae_clade$status == "reference"),]$tip_labels, collapse = ", ")
Mesoniviridae_clade_df$followup_phylo_group <- "Mesoniviridae_clade"

Mesoniviridae_clade_df$suggested_vir_sp
# Mesoniviridae_clade_df[which(Mesoniviridae_clade_df$contig_id == ""),]$suggested_vir_sp <- "_virus"


sort(Mesoniviridae_clade_df$suggested_vir_sp)

### further discarding contigs that are too short, partial, or potential EVEs
Mesoniviridae_drop_contigs <- c()

Mesoniviridae_clade_df_drop <- Mesoniviridae_clade_df[which(Mesoniviridae_clade_df$contig_id %in% Mesoniviridae_drop_contigs),]
Mesoniviridae_clade_df_keep <- Mesoniviridae_clade_df[which(!Mesoniviridae_clade_df$contig_id %in% Mesoniviridae_drop_contigs),]

sort(Mesoniviridae_clade_df_keep$suggested_vir_sp)


length(Mesoniviridae_clade[which(Mesoniviridae_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Mesoniviridae_clade/")
write.table(x = Mesoniviridae_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Mesoniviridae_clade/Mesoniviridae_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Mesoniviridae_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Mesoniviridae_clade/Mesoniviridae_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)






################################### Reovirales_clade
Reovirales_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "Duplornaviricota_Reovirales_AEQ75466.1_VP1__Scylla_serrata_reovirus_SZ_2007_")),]$node:
                                     tip_df_ord[which(str_detect(tip_df_ord$label, "Duplornaviricota_Reovirales_APG79093.1_RNA_dependent_RNA_polymerase__Hubei_lepidoptera_virus_3_")),]$node, ]




Reovirales_clade_df <- all_df[which(all_df$orf_id %in% 
                                         Reovirales_clade[which(Reovirales_clade$status == "novel"),]$tip_labels),] 



## only keep original contigs, not consensuses from mapping
Reovirales_clade_df1 <- Reovirales_clade_df[which(Reovirales_clade_df$original_sample_ID == Reovirales_clade_df$read_map_consensus_sample_ID | 
                                                          is.na(Reovirales_clade_df$read_map_consensus_sample_ID)),] %>% 
  discard(~all(is.na(.)))




sort(Reovirales_clade[which(Reovirales_clade$status == "novel"),]$tip_labels)


Reovirales_clade_df$ref_suggestions <- paste0(Reovirales_clade[which(Reovirales_clade$status == "reference"),]$tip_labels, collapse = ", ")
Reovirales_clade_df$followup_phylo_group <- "Reovirales_clade"

Reovirales_clade_df$suggested_vir_sp
Reovirales_clade_df[which(Reovirales_clade_df$contig_id == "F1168_NODE_3_length_4269_cov_1634.693878"),]$suggested_vir_sp <- "Usnessius_virus"


sort(Reovirales_clade_df$suggested_vir_sp)

### further discarding contigs that are too short, partial, or potential EVEs
Reovirales_drop_contigs <- c()

Reovirales_clade_df_drop <- Reovirales_clade_df[which(Reovirales_clade_df$contig_id %in% Reovirales_drop_contigs),]
Reovirales_clade_df_keep <- Reovirales_clade_df[which(!Reovirales_clade_df$contig_id %in% Reovirales_drop_contigs),]

sort(Reovirales_clade_df_keep$suggested_vir_sp)


length(Reovirales_clade[which(Reovirales_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Reovirales_clade/")
write.table(x = Reovirales_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Reovirales_clade/Reovirales_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Reovirales_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Reovirales_clade/Reovirales_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)






################################### Ghabrivirales_clade
Ghabrivirales_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "NOT_ASSIGNED_APG77587.1_hypothetical_protein")),]$node:
                                  tip_df_ord[which(str_detect(tip_df_ord$label, "F1186_NODE_10_length_7325_cov_98.964924_2")),]$node, ]




Ghabrivirales_clade_df <- all_df[which(all_df$orf_id %in% 
                                      Ghabrivirales_clade[which(Ghabrivirales_clade$status == "novel"),]$tip_labels),] 



## only keep original contigs, not consensuses from mapping
Ghabrivirales_clade_df1 <- Ghabrivirales_clade_df[which(Ghabrivirales_clade_df$original_sample_ID == Ghabrivirales_clade_df$read_map_consensus_sample_ID | 
                                                    is.na(Ghabrivirales_clade_df$read_map_consensus_sample_ID)),] %>% 
  discard(~all(is.na(.)))




sort(Ghabrivirales_clade[which(Ghabrivirales_clade$status == "novel"),]$tip_labels)


Ghabrivirales_clade_df$ref_suggestions <- paste0(Ghabrivirales_clade[which(Ghabrivirales_clade$status == "reference"),]$tip_labels, collapse = ", ")
Ghabrivirales_clade_df$followup_phylo_group <- "Ghabrivirales_clade"

Ghabrivirales_clade_df$suggested_vir_sp
Ghabrivirales_clade_df[which(Ghabrivirales_clade_df$contig_id == "F1409_NODE_1832_length_1127_cov_2.578358"),]$suggested_vir_sp <- "Vistafter_virus"
Ghabrivirales_clade_df[which(Ghabrivirales_clade_df$contig_id == "F1469_NODE_738_length_2404_cov_1.864198"),]$suggested_vir_sp <- "Wrunlize_virus"
Ghabrivirales_clade_df[which(Ghabrivirales_clade_df$contig_id == "F1183_NODE_2_length_5307_cov_3.704303"),]$suggested_vir_sp <- "Adlissaust_virus"
Ghabrivirales_clade_df[which(Ghabrivirales_clade_df$contig_id == "F1470_NODE_2_length_7157_cov_428.287525"),]$suggested_vir_sp <- "Aftionautor_virus"
Ghabrivirales_clade_df[which(Ghabrivirales_clade_df$contig_id == "F1490_NODE_6_length_7134_cov_7.580732"),]$suggested_vir_sp <- "Belfgalbs_virus"
Ghabrivirales_clade_df[which(Ghabrivirales_clade_df$contig_id == "F1183_NODE_376_length_1516_cov_4.340862"),]$suggested_vir_sp <- "Coundenality_virus"
Ghabrivirales_clade_df[which(Ghabrivirales_clade_df$contig_id == "F1497_NODE_7_length_5809_cov_200.726625"),]$suggested_vir_sp <- "Demadles_virus"
Ghabrivirales_clade_df[which(Ghabrivirales_clade_df$contig_id == "F1497_NODE_15_length_4526_cov_148.730485"),]$suggested_vir_sp <- "Fanomonover_virus"
Ghabrivirales_clade_df[which(Ghabrivirales_clade_df$contig_id == "F1187_NODE_10_length_7026_cov_14.692727"),]$suggested_vir_sp <- "Flatesals_virus"
Ghabrivirales_clade_df[which(Ghabrivirales_clade_df$contig_id == "F1367_NODE_33_length_3638_cov_33.225230"),]$suggested_vir_sp <- "Gainquerated_virus"



sort(Ghabrivirales_clade_df$suggested_vir_sp)

### further discarding contigs that are too short, partial, or potential EVEs
Ghabrivirales_drop_contigs <- c()

Ghabrivirales_clade_df_drop <- Ghabrivirales_clade_df[which(Ghabrivirales_clade_df$contig_id %in% Ghabrivirales_drop_contigs),]
Ghabrivirales_clade_df_keep <- Ghabrivirales_clade_df[which(!Ghabrivirales_clade_df$contig_id %in% Ghabrivirales_drop_contigs),]

sort(Ghabrivirales_clade_df_keep$suggested_vir_sp)


length(Ghabrivirales_clade[which(Ghabrivirales_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/")
write.table(x = Ghabrivirales_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/Ghabrivirales_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Ghabrivirales_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/Ghabrivirales_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)








################################### Negevirus_clade
Negevirus_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "Kitrinoviricota_Martellivirales_YP_009345002.1_RdRp__Wuhan_insect_virus_9_")),]$node:
                                     tip_df_ord[which(str_detect(tip_df_ord$label, "NOT_ASSIGNED_YP_009553328.1_hypothetical_protein__Osedax_japonicus_RNA_virus_1_")),]$node, ]




Negevirus_clade_df <- all_df[which(all_df$orf_id %in% 
                                         Negevirus_clade[which(Negevirus_clade$status == "novel"),]$tip_labels),] 



## only keep original contigs, not consensuses from mapping
Negevirus_clade_df1 <- Negevirus_clade_df[which(Negevirus_clade_df$original_sample_ID == Negevirus_clade_df$read_map_consensus_sample_ID | 
                                                          is.na(Negevirus_clade_df$read_map_consensus_sample_ID)),] %>% 
  discard(~all(is.na(.)))




sort(Negevirus_clade[which(Negevirus_clade$status == "novel"),]$tip_labels)


Negevirus_clade_df$ref_suggestions <- paste0(Negevirus_clade[which(Negevirus_clade$status == "reference"),]$tip_labels, collapse = ", ")
Negevirus_clade_df$followup_phylo_group <- "Negevirus_clade"

Negevirus_clade_df$suggested_vir_sp
Negevirus_clade_df[which(Negevirus_clade_df$contig_id == "F1018_NODE_1_length_14998_cov_54.605501"),]$suggested_vir_sp <- "Saljeroized_virus"
Negevirus_clade_df[which(Negevirus_clade_df$contig_id == "F1038_NODE_1_length_11101_cov_60.073149"),]$suggested_vir_sp <- "Simimissemian_virus"
Negevirus_clade_df[which(Negevirus_clade_df$contig_id == "F1147_NODE_1_length_11305_cov_1652.532000"),]$suggested_vir_sp <- "Snubberrored_virus"
Negevirus_clade_df[which(Negevirus_clade_df$contig_id == "F1359_NODE_1_length_14961_cov_86.217563"),]$suggested_vir_sp <- "Undsholming_virus"
Negevirus_clade_df[which(Negevirus_clade_df$contig_id == "F1359_NODE_2_length_14873_cov_82.619989"),]$suggested_vir_sp <- "Verstrife_virus"
Negevirus_clade_df[which(Negevirus_clade_df$contig_id == "F1371_NODE_2_length_10261_cov_47.214776"),]$suggested_vir_sp <- "Exquased_virus"






sort(unique(Negevirus_clade_df$suggested_vir_sp))

### further discarding contigs that are too short, partial, or potential EVEs
Negevirus_drop_contigs <- c()

Negevirus_clade_df_drop <- Negevirus_clade_df[which(Negevirus_clade_df$contig_id %in% Negevirus_drop_contigs),]
Negevirus_clade_df_keep <- Negevirus_clade_df[which(!Negevirus_clade_df$contig_id %in% Negevirus_drop_contigs),]

sort(Negevirus_clade_df_keep$suggested_vir_sp)


length(Negevirus_clade[which(Negevirus_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Negevirus_clade/")
write.table(x = Negevirus_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Negevirus_clade/Negevirus_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Negevirus_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Negevirus_clade/Negevirus_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)




################################### Negevirus_Tanay_clade
Negevirus_Tanay_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "NOT_ASSIGNED_AFI24675.1_hypothetical_protein_1__Santana_virus_")),]$node:
                                 tip_df_ord[which(str_detect(tip_df_ord$label, "F1530_NODE_2_length_9569_cov_814.715997_1")),]$node, ]




Negevirus_Tanay_clade_df <- all_df[which(all_df$orf_id %in% 
                                     Negevirus_Tanay_clade[which(Negevirus_Tanay_clade$status == "novel"),]$tip_labels),] 



## only keep original contigs, not consensuses from mapping
Negevirus_Tanay_clade_df1 <- Negevirus_Tanay_clade_df[which(Negevirus_Tanay_clade_df$original_sample_ID == Negevirus_Tanay_clade_df$read_map_consensus_sample_ID | 
                                                  is.na(Negevirus_Tanay_clade_df$read_map_consensus_sample_ID)),] %>% 
  discard(~all(is.na(.)))




sort(Negevirus_Tanay_clade[which(Negevirus_Tanay_clade$status == "novel"),]$tip_labels)


Negevirus_Tanay_clade_df$ref_suggestions <- paste0(Negevirus_Tanay_clade[which(Negevirus_Tanay_clade$status == "reference"),]$tip_labels, collapse = ", ")
Negevirus_Tanay_clade_df$followup_phylo_group <- "Negevirus_Tanay_clade"

Negevirus_Tanay_clade_df$suggested_vir_sp

sort(unique(Negevirus_Tanay_clade_df$suggested_vir_sp))

### further discarding contigs that are too short, partial, or potential EVEs
Negevirus_drop_contigs <- c()

Negevirus_Tanay_clade_df_drop <- Negevirus_Tanay_clade_df[which(Negevirus_Tanay_clade_df$contig_id %in% Negevirus_drop_contigs),]
Negevirus_Tanay_clade_df_keep <- Negevirus_Tanay_clade_df[which(!Negevirus_Tanay_clade_df$contig_id %in% Negevirus_drop_contigs),]

sort(Negevirus_Tanay_clade_df_keep$suggested_vir_sp)


length(Negevirus_Tanay_clade[which(Negevirus_Tanay_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Negevirus_Tanay_clade/")
write.table(x = Negevirus_Tanay_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Negevirus_Tanay_clade/Negevirus_Tanay_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Negevirus_Tanay_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Negevirus_Tanay_clade/Negevirus_Tanay_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)











################################### Tymovirales_clade
Tymovirales_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "Kitrinoviricota_Tymovirales_QCC30253.1_polymerase__Peach_marafivirus_D_")),]$node:
                                 tip_df_ord[which(str_detect(tip_df_ord$label, "Kitrinoviricota_Tymovirales_AMN92730.1_replication_associated_polyprotein__Fusarium_graminearum_mycotymovirus_1_")),]$node, ]




Tymovirales_clade_df <- all_df[which(all_df$orf_id %in% 
                                     Tymovirales_clade[which(Tymovirales_clade$status == "novel"),]$tip_labels),] 



## only keep original contigs, not consensuses from mapping
Tymovirales_clade_df1 <- Tymovirales_clade_df[which(Tymovirales_clade_df$original_sample_ID == Tymovirales_clade_df$read_map_consensus_sample_ID | 
                                                  is.na(Tymovirales_clade_df$read_map_consensus_sample_ID)),] %>% 
  discard(~all(is.na(.)))




sort(Tymovirales_clade[which(Tymovirales_clade$status == "novel"),]$tip_labels)


Tymovirales_clade_df$ref_suggestions <- paste0(Tymovirales_clade[which(Tymovirales_clade$status == "reference"),]$tip_labels, collapse = ", ")
Tymovirales_clade_df$followup_phylo_group <- "Tymovirales_clade"

Tymovirales_clade_df$suggested_vir_sp
# Tymovirales_clade_df[which(Tymovirales_clade_df$contig_id == ""),]$suggested_vir_sp <- ""








### further discarding contigs that are too short, partial, or potential EVEs
Tymovirales_drop_contigs <- c()

Tymovirales_clade_df_drop <- Tymovirales_clade_df[which(Tymovirales_clade_df$contig_id %in% Tymovirales_drop_contigs),]
Tymovirales_clade_df_keep <- Tymovirales_clade_df[which(!Tymovirales_clade_df$contig_id %in% Tymovirales_drop_contigs),]

sort(Tymovirales_clade_df_keep$suggested_vir_sp)


length(Tymovirales_clade[which(Tymovirales_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Tymovirales_clade/")
write.table(x = Tymovirales_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Tymovirales_clade/Tymovirales_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Tymovirales_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Tymovirales_clade/Tymovirales_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)






################################### Permutotetraviridae_clade
Permutotetraviridae_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "NOT_ASSIGNED_Jiangsu_2_permutotetra_like_virus_2")),]$node:
                                   tip_df_ord[which(str_detect(tip_df_ord$label, "Permutotetraviridae_Permutotetraviridae_guangxi_permutotetra_like_virus_1")),]$node, ]




Permutotetraviridae_clade_df <- all_df[which(all_df$orf_id %in% 
                                       Permutotetraviridae_clade[which(Permutotetraviridae_clade$status == "novel"),]$tip_labels),] 



## only keep original contigs, not consensuses from mapping
Permutotetraviridae_clade_df1 <- Permutotetraviridae_clade_df[which(Permutotetraviridae_clade_df$original_sample_ID == Permutotetraviridae_clade_df$read_map_consensus_sample_ID | 
                                                      is.na(Permutotetraviridae_clade_df$read_map_consensus_sample_ID)),] %>% 
  discard(~all(is.na(.)))




sort(Permutotetraviridae_clade[which(Permutotetraviridae_clade$status == "novel"),]$tip_labels)


Permutotetraviridae_clade_df$ref_suggestions <- paste0(Permutotetraviridae_clade[which(Permutotetraviridae_clade$status == "reference"),]$tip_labels, collapse = ", ")
Permutotetraviridae_clade_df$followup_phylo_group <- "Permutotetraviridae_clade"

Permutotetraviridae_clade_df$suggested_vir_sp
Permutotetraviridae_clade_df[which(Permutotetraviridae_clade_df$contig_id == "F1052_NODE_9_length_4423_cov_7.729625"),]$suggested_vir_sp <- "Gretteness_virus"








### further discarding contigs that are too short, partial, or potential EVEs
Permutotetraviridae_drop_contigs <- c()

Permutotetraviridae_clade_df_drop <- Permutotetraviridae_clade_df[which(Permutotetraviridae_clade_df$contig_id %in% Permutotetraviridae_drop_contigs),]
Permutotetraviridae_clade_df_keep <- Permutotetraviridae_clade_df[which(!Permutotetraviridae_clade_df$contig_id %in% Permutotetraviridae_drop_contigs),]

sort(Permutotetraviridae_clade_df_keep$suggested_vir_sp)


length(Permutotetraviridae_clade[which(Permutotetraviridae_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Permutotetraviridae_clade/")
write.table(x = Permutotetraviridae_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Permutotetraviridae_clade/Permutotetraviridae_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Permutotetraviridae_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Permutotetraviridae_clade/Permutotetraviridae_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)







################################### Nodamuvirales_clade
Nodamuvirales_clade <- tip_df_ord[ tip_df_ord[which(str_detect(tip_df_ord$label, "Kitrinoviricota_Nodamuvirales_ya20_JAAOEH010001993_1_JAAOEH010001993.1_26281228_OV.3_NODE_1983_truseq_orf.11957")),]$node:
                                           tip_df_ord[which(str_detect(tip_df_ord$label, "Kitrinoviricota_Nodamuvirales_ADW54431.1_RNA_dependent_RNA_polymerase__Santeuil_nodavirus_")),]$node, ]




Nodamuvirales_clade_df <- all_df[which(all_df$orf_id %in% 
                                               Nodamuvirales_clade[which(Nodamuvirales_clade$status == "novel"),]$tip_labels),] 



## only keep original contigs, not consensuses from mapping
Nodamuvirales_clade_df1 <- Nodamuvirales_clade_df[which(Nodamuvirales_clade_df$original_sample_ID == Nodamuvirales_clade_df$read_map_consensus_sample_ID | 
                                                                      is.na(Nodamuvirales_clade_df$read_map_consensus_sample_ID)),] %>% 
  discard(~all(is.na(.)))




sort(Nodamuvirales_clade[which(Nodamuvirales_clade$status == "novel"),]$tip_labels)


Nodamuvirales_clade_df$ref_suggestions <- paste0(Nodamuvirales_clade[which(Nodamuvirales_clade$status == "reference"),]$tip_labels, collapse = ", ")
Nodamuvirales_clade_df$followup_phylo_group <- "Nodamuvirales_clade"

Nodamuvirales_clade_df$suggested_vir_sp
Nodamuvirales_clade_df[which(Nodamuvirales_clade_df$contig_id == "F1531_NODE_222_length_3068_cov_2435.778626"),]$suggested_vir_sp <- "Murmaxied_virus"








### further discarding contigs that are too short, partial, or potential EVEs
Nodamuvirales_drop_contigs <- c()

Nodamuvirales_clade_df_drop <- Nodamuvirales_clade_df[which(Nodamuvirales_clade_df$contig_id %in% Nodamuvirales_drop_contigs),]
Nodamuvirales_clade_df_keep <- Nodamuvirales_clade_df[which(!Nodamuvirales_clade_df$contig_id %in% Nodamuvirales_drop_contigs),]

sort(Nodamuvirales_clade_df_keep$suggested_vir_sp)


length(Nodamuvirales_clade[which(Nodamuvirales_clade$status == "reference"),]$tip_labels)


dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/")
write.table(x = Nodamuvirales_clade_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/Nodamuvirales_clade_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = Nodamuvirales_clade_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/Nodamuvirales_clade_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)











################################### other_clades

phylo_dir_paths <- list.dirs(path = "/full_path_to/wd/RdRp_scan/analysis/phylo", full.names = T, recursive = F)
clade_dir_paths <- phylo_dir_paths[which(str_detect(phylo_dir_paths, "\\_clade"))]
clade_dir_paths <- clade_dir_paths[which(!str_detect(clade_dir_paths, "other_clade"))] ## make sure to exclude those not considered for followup phylogenies

count = 0
for (cdp in clade_dir_paths){
  cons_df <- read_delim(file = paste0(cdp, "/consensus_df.txt"))
  if (count == 0){
    concat_cons_df <- cons_df
  } else{
    concat_cons_df <- bind_rows(concat_cons_df, cons_df)
  }
  count = count + 1
}

unique(concat_cons_df$suggested_vir_sp)
unique(concat_cons_df$contig_id)

write.table(x = concat_cons_df, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/consensus_df_follow_up_clades.txt", 
            append = F, quote = F, sep = "\t",row.names = F, col.names = T)



other_clades <- tip_df_ord[which(!tip_df_ord$tip_labels %in% unique(concat_cons_df$orf_id) & str_detect(tip_df_ord$label, "F1[[:digit:]]{3}\\_NODE")),] 

other_clades_df <- all_df[which(all_df$orf_id %in% 
                                  other_clades[which(other_clades$status == "novel"),]$tip_labels),] 


## only keep original contigs, not consensuses from mapping
other_clades_df1 <- other_clades_df[which(other_clades_df$original_sample_ID == other_clades_df$read_map_consensus_sample_ID | 
                                                          is.na(other_clades_df$read_map_consensus_sample_ID)),] %>% 
  discard(~all(is.na(.)))

sum(unique(other_clades_df$contig_id) == unique(other_clades_df1$contig_id)) == length(unique(other_clades_df1$contig_id)) #### check




##### Check if some viruses were missed in the follow-up clades
colnames(other_clades_df1)


other_clades_df1_minimal <- other_clades_df1 %>% select(orf_id, contig_id, 
                                                        rdrp_scan_vir_group, BestHit_class, BestHit_order, 
                                                        BestHit_family, BestHit_genus, BestHit_species, BestHit_pident)


other_clades_df1_minimal_to_check <- other_clades_df1_minimal[which(other_clades_df1_minimal$BestHit_class %in% 
                                                                      c("Monjiviricetes", "Insthoviricetes", "Ellioviricetes") | 
                                 other_clades_df1_minimal$BestHit_species %in% 
                                   c("Tanay virus") ), ]
##### see if any of these were dropped for some reason
phylo_dir_paths <- list.dirs(path = "/full_path_to/wd/RdRp_scan/analysis/phylo", full.names = T, recursive = F)
clade_dir_paths <- phylo_dir_paths[which(str_detect(phylo_dir_paths, "\\_clade"))]
clade_dir_paths <- clade_dir_paths[which(!str_detect(clade_dir_paths, "other_clade"))] ## make sure to exclude those not considered for followup phylogenies
count1 = 0
for (cdp in clade_dir_paths){
  skip_to_next <- FALSE
  clade_name <- str_extract(cdp, "[[:alpha:]]+_clade$")
  tryCatch(
    expr = {
      dropped_df <- read_delim(file = paste0(cdp, "/", clade_name, "_novel_sequence_data_dropped.txt")) %>% 
        select(orf_id, contig_id)
    },
    error = function(e){ 
      skip_to_next <<- TRUE
    }
  )
  if(skip_to_next) { next } 
  
  dropped_df$clade_name <- clade_name

  if (count1 == 0){
    concat_dropped_df <- dropped_df
  } else{
    concat_dropped_df <- bind_rows(concat_dropped_df, dropped_df)
  }
  count1 = count1 + 1
}
other_clades_df1_minimal_to_check_dropped_from_phylo <- other_clades_df1_minimal_to_check[which(other_clades_df1_minimal_to_check$orf_id %in% concat_dropped_df$orf_id),]
#### These are either partial genomes/segments or EVEs, because of tehir lengths. We can completely drop them from the analysis now

dir.create("/full_path_to/wd/RdRp_scan/analysis/phylo/other_clades")

write.table(x = other_clades_df1_minimal_to_check_dropped_from_phylo, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/other_clades/contigs_dropped_from_analysis_and_submission.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)



other_clades_df1_minimal_to_check <- other_clades_df1_minimal_to_check[which(!other_clades_df1_minimal_to_check$orf_id %in% other_clades_df1_minimal_to_check_dropped_from_phylo$orf_id),]

write.table(x = other_clades_df1_minimal_to_check, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/other_clades/other_clades_contigs_to_check.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)


other_clades_df <- other_clades_df[which(!other_clades_df$orf_id %in% other_clades_df1_minimal_to_check_dropped_from_phylo$orf_id),]


sort(other_clades[which(other_clades$status == "novel"),]$tip_labels)

other_clades_df$ref_suggestions <- paste0(other_clades[which(other_clades$status == "reference"),]$tip_labels, collapse = ", ")
other_clades_df$followup_phylo_group <- "other_clades"

other_clades_df$suggested_vir_sp






### further discarding contigs that are too short, partial, or potential EVEs
other_clades_drop_contigs <- c()

other_clades_df_drop <- other_clades_df[which(other_clades_df$contig_id %in% other_clades_drop_contigs),]
other_clades_df_keep <- other_clades_df[which(!other_clades_df$contig_id %in% other_clades_drop_contigs),]

sort(other_clades_df_keep$suggested_vir_sp)


length(other_clades[which(other_clades$status == "reference"),]$tip_labels)

write.table(x = other_clades_df_keep, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/other_clades/other_clades_novel_sequence_data.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)
write.table(x = other_clades_df_drop, 
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/other_clades/other_clades_novel_sequence_data_dropped.txt",
            append = F, quote = F,sep = "\t", row.names = F, col.names = T)




















