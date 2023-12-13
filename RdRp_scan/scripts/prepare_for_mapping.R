library(tidyverse)
library(seqinr)

### Ran extract_R_contigs.R on cluster first
### Here I am subsetting only the updated contigs and putting chimeric contigs where needed.
R_fasta <- read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/R_contigs_non_cellular.fasta",seqtype = "DNA") ## all R contigs, without chimeric (scaffolds).


### Import R df
R_df_curated <- read_delim("/full_path_to/wd/RdRp_scan/analysis/R_orf_df_curated.txt", delim = "\t")

length(sort(unique(R_df_curated$contig_id)))

colnames(R_df_curated)[str_detect(colnames(R_df_curated), "contig")]

##### Import set G where we'll see chimeric contig IDs corresponding to RdRp only original contig IDs that are in R set.
G_df <- read_delim(file = "/full_path_to/wd/RdRp_scan/analysis/contig_set_G.txt", delim = "\t")


### Import previously mapped and curated contigs and get chimeras from there.
df4 <- read_delim("/full_path_to/wd/RdRp_search/analysis/all_families_and_genus_classified/phylo_classified_contigs_and_consensus_w_metadata_correctedAdableaps.txt")
length(sort(unique(df4$contig_ID)))


##### Verify set difference
dmnd_curated_not_in_R <- setdiff(sort(unique(df4$contig_ID)), sort(unique(R_df_curated$contig_id)))
setdiff(sort(unique(R_df_curated$contig_id)), sort(unique(df4$contig_ID)))

###
chimeric_R <- G_df[which(G_df$contig_id_set_AuB %in% dmnd_curated_not_in_R),]

### Substitute for chimeric in R

R_set <- sort(unique(R_df_curated$contig_id))
R_set[which(R_set %in% chimeric_R$contig_id)] <- chimeric_R$contig_id_set_AuB 

R_set

write.table(x = R_set, file = "/full_path_to/wd/RdRp_scan/analysis/contig_set_R_final_with_chimeras.txt",
            sep = "\t",append = F, quote = F,row.names = F,col.names = F)


R_df_curated[which(R_df_curated$contig_id %in% chimeric_R$contig_id),]$contig_id <- chimeric_R$contig_id_set_AuB

write.table(x = R_df_curated, file = "/full_path_to/wd/RdRp_scan/analysis/R_orf_df_curated_with_chimeric_contig_id.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)


##### Get chimeric fastas
dir.create("/full_path_to/wd/RdRp_scan/analysis/R_contigs")


# df4[which(df4$contig_ID %in% chimeric_R$contig_id_set_AuB),] %>% select(vir_group, contig_ID)
# 
# s1 <- read.fasta(file = "/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/Flaviviridae/read_mapping/map_refs/F1416_NODEs_39_8_2_finalCxFVlike.fasta", 
#                  seqtype = "DNA")
# s2 <- read.fasta(file = "/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/Permutotetraviridae/read_mapping/sense_corrected_contigs/F1052_NODE_13_663_chimera.fasta", 
#                  seqtype = "DNA")
# s3 <- read.fasta(file = "/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/Mesoniviridae/read_mapping/sense_corrected_contigs/F1576_NODEs_2_25945_25414_25622_8_26117_15_chimera.fasta", 
#                  seqtype = "DNA")
# s4 <- read.fasta(file = "/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/Totiviridae/read_mapping/selected_cons3/F1506_NODEs_12_58_205_F1506_MapVerif_F1506_detect_F1506_round2_corrected.fa", 
#                  seqtype = "DNA")
# 
# chimeric_contigs <- c(s1, s2, s3, s4)
# 
# names(chimeric_contigs)
# 
# names_correct <- c("F1416_NODEs_39_8_2_finalCxFVlike","F1052_NODE_13_663_chimera", 
#                    "F1576_NODEs_2_25945_25414_25622_8_26117_15_chimera", "F1506_NODEs_12_58_205")
# 
# write.fasta(sequences = chimeric_contigs, names = names_correct, file.out = "/full_path_to/wd/RdRp_scan/analysis/R_contigs/chimeric_contigs.fasta")
chimeric_contigs <- read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/R_contigs/chimeric_contigs.fasta", 
                             seqtype = "DNA")

names()

R_set[which(!R_set %in% names(chimeric_contigs))]

R_fasta_non_chimeric <- R_fasta[R_set[which(!R_set %in% names(chimeric_contigs))]]

R_fasta_updated <- c(R_fasta_non_chimeric, chimeric_contigs)


write.fasta(sequences = R_fasta_updated, names = names(R_fasta_updated), 
            file.out = "/full_path_to/wd/RdRp_scan/analysis/R_contigs/R_contigs_with_chimeric.fasta")


##### Correct strand sense according to RdRp scan RdRp domain sense
dir.create("/full_path_to/wd/RdRp_scan/analysis/read_mapping/contigs_correct_sense", recursive = T)
R_df_curated$orf_sense


for (i in unique(R_df_curated$contig_id)){
  cdf <- R_df_curated[which(R_df_curated$contig_id == i),]
  
  if (nrow(cdf) == 0){
    next
  }
  
  if (cdf$orf_sense == "REVERSE"){
    temp <- R_fasta_updated[i]
    correct_fasta <- list(as.SeqFastadna(object = rev(comp(temp[[1]])), name = names(R_fasta_updated[i]), Annot = names(R_fasta_updated[i])))
    names(correct_fasta) <- names(R_fasta_updated[i])
  }else(
    correct_fasta <- R_fasta_updated[i]
  )
  write.fasta(sequences = correct_fasta, names = names(correct_fasta), 
              file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/read_mapping/contigs_correct_sense/", names(correct_fasta), ".fasta"))
}


for (i in unique(chimeric_R$contig_id_set_AuB)){
  
  correct_fasta <- R_fasta_updated[i]
  write.fasta(sequences = correct_fasta, names = names(correct_fasta), 
              file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/read_mapping/contigs_correct_sense/", names(correct_fasta), ".fasta"))
}

### Import sample list
samples <- read_lines("/full_path_to/wd/RdRp_scan/analysis/read_mapping/sample_list_RdRpScan_Round1.txt")
samples_on_geva_nodes <- str_extract(R_df_curated$contig_id, "^[[:alnum:]]+(?=\\_)")
samples_on_common_nodes <- setdiff(samples, samples_on_geva_nodes)

write_lines(x = samples_on_geva_nodes, 
            file = "/full_path_to/wd/RdRp_scan/analysis/read_mapping/sample_list_RdRpScan_Round1_geva_nodes.txt", 
            append = F)
write_lines(x = samples_on_common_nodes, 
            file = "/full_path_to/wd/RdRp_scan/analysis/read_mapping/sample_list_RdRpScan_Round1_common_nodes.txt", 
            append = F)


R_df_curated1 <- R_df_curated
R_df_curated1[which(R_df_curated1$contig_set == "G"),]$contig_id
R_df_curated1$contig_id_set_AuB <- R_df_curated1$contig_id
R_df_curated1[which(R_df_curated1$contig_id %in% chimeric_R$contig_id),]$contig_id_set_AuB <- chimeric_R$contig_id_set_AuB

N_set_priority_1 <-  R_df_curated1[which(R_df_curated1$contig_set == "N"),]$contig_id
G_set_priority_2 <- R_df_curated1[which(R_df_curated1$contig_set == "G"),]$contig_id_set_AuB


### Copy into folders by mapping priority (priority simply refers to the order in which I launched the mappings on HPC on different nodes: dedicated or common)
ref_folder_priority1 <- "/full_path_to/wd/RdRp_scan/analysis/read_mapping/contigs_correct_sense_priority1"
ref_folder_priority2 <- "/full_path_to/wd/RdRp_scan/analysis/read_mapping/contigs_correct_sense_priority2"
dir.create(ref_folder_priority1)
dir.create(ref_folder_priority2)

for (cont in N_set_priority_1){
  f <- paste0("/full_path_to/wd/RdRp_scan/analysis/read_mapping/contigs_correct_sense/", 
              cont, 
              ".fasta")
  f_new <- paste0(ref_folder_priority1, "/", 
              cont, 
              ".fasta")
  file.copy(from = f, to = f_new, overwrite = T,copy.mode = T, copy.date = T )
}

for (cont in G_set_priority_2){
  f <- paste0("/full_path_to/wd/RdRp_scan/analysis/read_mapping/contigs_correct_sense/", 
              cont, 
              ".fasta")
  f_new <- paste0(ref_folder_priority2, "/", 
                  cont, 
                  ".fasta")
  file.copy(from = f, to = f_new, overwrite = T,copy.mode = T, copy.date = T )
}









