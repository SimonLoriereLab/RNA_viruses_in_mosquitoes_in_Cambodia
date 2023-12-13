library(tidyverse)
library(seqinr)
library(ggtree)

df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/R_set_merged_hit_data_from_2_analyses.txt", delim = "\t")
unique(df$contig_id)
unique(df$consensus_cds_ID)

other_clades_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/other_clades/other_clades_novel_sequence_data.txt")


viruses_not_in_followup_clades <- sort(unique(other_clades_df$suggested_vir_sp))
yet_to_name <- viruses_not_in_followup_clades[which(str_detect(string = viruses_not_in_followup_clades, pattern = "_NODE"))]

nonsense_words_verified <- str_split("Abillumpty
Abradivilines
Afficallowely
Albuilitusly
Allescrokeness
Ammanver
Anarindottina
Appingrabludly
Aurovirths
Autioncer
Bakeuprily
Balkhesia
Cahewisemili
Cioughtles
Clainforians
Comblicompons
Constranars
Corringskyot
Crearetracons
Dialanscovere
Distredfologra
Dranizinef
Enoctampy
Erminizes
Eveninate
Fluousser
Grannonnaho
Idembrear
Inaturn
Indubiout
Irestatersion
Irricalied
Judissly
Knocate
Mcelliess
Moundetelly
Mucklity
Murpreadaperic
Nogreepyth
Opagrainglors
Panthumpness
Pectimpers
Proptituckey
Prosveged
Proxidires
Rapticial
Recoroporters
Reflattratine
Rejoyoven
Replairmators
Reromeracy
Robeancess
Scarelveness
Shoonismitches
Stratereighted
Thogencentays
Troveremish
Unnamentia
Unprectory
Vaneathly", "\\s")[[1]]


naming_df <- tibble(contig_id = yet_to_name, suggested_vir_sp = paste0(nonsense_words_verified[1:length(yet_to_name)], "_virus"))



other_clades_df_named_earlier <- other_clades_df[which(str_detect(other_clades_df$suggested_vir_sp, pattern = "_NODE", negate = T)),]

other_clades_df_yet_to_name <- other_clades_df[which(str_detect(other_clades_df$suggested_vir_sp, pattern = "_NODE", negate = F)),]

other_clades_df_yet_to_name$suggested_vir_sp <- NULL

other_clades_df_freshly_named <- left_join(other_clades_df_yet_to_name, naming_df, by = "contig_id")

other_clades_df_all_named <- bind_rows(other_clades_df_named_earlier, other_clades_df_freshly_named) 

### select consensus sequences based on mapping coverage
### first round of read mapping results
read_map_df <-  read_delim("/full_path_to/wd/RdRp_scan/analysis/read_mapping/mapping_summary_bowtie2_RdRpScan_Round1_analysis/results_mapped_only.txt", 
                           delim = "\t")


read_map_df_other_clades <- read_map_df[which(read_map_df$ref_id %in% other_clades_df_all_named$contig_id & read_map_df$selected == "yes"),]



read_map_df_other_clades$RdRp_scan_map_consensus_id <- paste0(read_map_df_other_clades$ref_id, "_",
                                                                     read_map_df_other_clades$sample_id, "_RdRpScan_Round1")

read_map_df_other_clades$contig_id <- read_map_df_other_clades$ref_id

consensus_df <- left_join(read_map_df_other_clades, other_clades_df_all_named, by  = "contig_id")



### Parse consensus IDs

consensus_df$original_sample_ID <- str_extract(consensus_df$RdRp_scan_map_consensus_id, "^F[[:digit:]]{4}(?=\\_)")

consensus_df$original_metaspades_contig_node_IDs <- str_extract(consensus_df$RdRp_scan_map_consensus_id, "(?<=^F[[:digit:]]{4}\\_NODE)[[s]]*\\_[[:digit:][\\_]]+(?=[[:alpha:]]+)")
consensus_df$original_metaspades_contig_node_IDs <- str_remove(consensus_df$original_metaspades_contig_node_IDs, "s")
consensus_df$original_metaspades_contig_node_IDs <- str_remove(consensus_df$original_metaspades_contig_node_IDs, "_$")
consensus_df$original_metaspades_contig_node_IDs <- str_replace_all(consensus_df$original_metaspades_contig_node_IDs, "_", "C")


consensus_df$read_map_consensus_sample_ID_ <- str_extract_all(consensus_df$RdRp_scan_map_consensus_id, "F[[:digit:]]{4}\\_")

consensus_df$read_map_consensus_sample_ID <- NA

take_last <- function(x){
  x <- x[[1]][[1]]
  y <- str_remove(x[length(x)], "_")
  return(y)
}

for (i in 1:nrow(consensus_df)){
  consensus_df[i,'read_map_consensus_sample_ID'] <- take_last(x = consensus_df[i,'read_map_consensus_sample_ID_'])
}

consensus_df$read_map_consensus_sample_ID_ <- NULL

consensus_df$read_map_consensus_sample_ID


### New sequence IDs for further analyses

consensus_df$vir_sp_seq_ID_v3 <- paste0(consensus_df$suggested_vir_sp, "_",
                                        consensus_df$original_sample_ID, 
                                        consensus_df$original_metaspades_contig_node_IDs,
                                        consensus_df$read_map_consensus_sample_ID)

##### For these viruses, only keep the consensus from the original sample, as there is no phylogeny follow up later.
consensus_df <- consensus_df[which(consensus_df$original_sample_ID == consensus_df$read_map_consensus_sample_ID),] %>% 
  select(vir_sp_seq_ID_v3, 
         original_metaspades_contig_node_IDs,
         original_sample_ID, read_map_consensus_sample_ID, 
         orf_id, contig_id, RdRp_scan_map_consensus_id, orf_start, orf_end, orf_sense) %>%
  distinct()



dir_cons1 <- "/full_path_to/wd/RdRp_scan/analysis/read_mapping/cons1_renamed/"
extracted_consensus_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/other_clades/extracted_cons1_renamed/"
dir.create(extracted_consensus_directory)
extracted_orf_nt_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/other_clades/extracted_cons1_orf_nt/"
dir.create(extracted_orf_nt_directory)
extracted_orf_aa_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/other_clades/extracted_cons1_orf_aa/"
dir.create(extracted_orf_aa_directory)
for (i in 1:nrow(consensus_df)){
  consensus_fasta <- seqinr::read.fasta(file = paste0(dir_cons1, consensus_df[i,]$RdRp_scan_map_consensus_id, ".fa"), seqtype = "DNA")
  extract_fn <- paste0(extracted_consensus_directory, consensus_df[i,]$vir_sp_seq_ID_v3, ".fa")
  seqinr::write.fasta(sequences = consensus_fasta, names = consensus_df[i,]$vir_sp_seq_ID_v3, file.out = extract_fn)
  ###### translate
  
  if(consensus_df[i,]$orf_sense == "FORWARD"){
    nt_seq <- as.character(consensus_fasta[[1]])
    orf_nt <- as.SeqFastadna(object = nt_seq[consensus_df[i,]$orf_start:consensus_df[i,]$orf_end], 
                             name = consensus_df[i,]$vir_sp_seq_ID_v3, 
                             Annot = paste0(">", consensus_df[i,]$vir_sp_seq_ID_v3))
  }else{
    nt_seq <- rev(comp(as.character(consensus_fasta[[1]]), ambiguous = TRUE))
    orf_nt <- rev(comp(as.SeqFastadna(object = nt_seq[consensus_df[i,]$orf_end:consensus_df[i,]$orf_start], 
                                      name = consensus_df[i,]$vir_sp_seq_ID_v3, 
                                      Annot = paste0(">", consensus_df[i,]$vir_sp_seq_ID_v3)), 
                       ambiguous = TRUE))
    
  }
  extract_fn2 <- paste0(extracted_orf_nt_directory, consensus_df[i,]$vir_sp_seq_ID_v3, ".fa")
  seqinr::write.fasta(sequences = orf_nt, names = consensus_df[i,]$vir_sp_seq_ID_v3, file.out = extract_fn2)
  orf_aa <- translate(seq = as.character(orf_nt),
                      frame = 0, sens = "F", NAstring = "X")
  extract_fn3 <- paste0(extracted_orf_aa_directory, consensus_df[i,]$vir_sp_seq_ID_v3, ".fa")
  seqinr::write.fasta(sequences = orf_aa, names = consensus_df[i,]$vir_sp_seq_ID_v3, file.out = extract_fn3)
}

write.table(x = consensus_df, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/other_clades/consensus_df.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)



#### Prepare consensus metadata
consensus_df$verif_mosq_sp


### import sample metadata
CMVM_sample_mdf <- read_delim("/full_path_to/wd/RdRp_scan/metadata/CMVM_full_metadata_updated_20211021.tsv",
                              delim = "\t") %>% select(sample_id = Sample, Location)

CMVM_sample_mdf$Country <- "Cambodia"
CMVM_sample_mdf$Location <- paste0(CMVM_sample_mdf$Country, ":", CMVM_sample_mdf$Location)
CMVM_sample_mdf$Year <- "2019"

colnames(CMVM_sample_mdf)

### import mosquito species verification data
mosq_mdf <- read_delim("/full_path_to/wd/mosquito_species_verification/analysis/SKA/CMVM_mosquito_species_verification_gradual_curation - Sheet1.tsv",
                       delim = "\t")


mdf <- left_join(CMVM_sample_mdf, mosq_mdf)  
colnames(consensus_df)
novel_seq_mdf <- consensus_df %>% 
  select(vir_sp_seq_ID_v3,
         original_metaspades_contig_node_IDs,
         original_sample_ID, read_map_consensus_sample_ID, 
         orf_id, contig_id)

novel_seq_mdf1 <- left_join(novel_seq_mdf, mdf, by  = c("read_map_consensus_sample_ID"  = "sample_id"))


novel_seq_mdf1$host_taxon <- novel_seq_mdf1$verif_mosq_sp
novel_seq_mdf1[which(novel_seq_mdf1$SKA_SNP_dist_dend == "could not be analyzed"),]$host_taxon <- paste0(novel_seq_mdf1[which(novel_seq_mdf1$SKA_SNP_dist_dend == "could not be analyzed"),]$morph_mosq_sp, "_morph" )

novel_seq_mdf1$host_genus <- str_extract(novel_seq_mdf1$host_taxon, "^[[:alpha:]]+")

novel_seq_mdf1[which(novel_seq_mdf1$SKA_SNP_dist_dend == "could not be analyzed"),]$host_genus<- paste0(novel_seq_mdf1[which(novel_seq_mdf1$SKA_SNP_dist_dend == "could not be analyzed"),]$host_genus, "_morph" )

write.table(x = novel_seq_mdf1, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/other_clades/updated_consensus_df.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)








