library(tidyverse)
library(seqinr)
library(treeio)

Quaranja_clade_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/Quaranja_clade_novel_sequence_data.txt")

### select consensus sequences based on mapping coverage
### first round of read mapping results
read_map_df <-  read_delim("/full_path_to/wd/RdRp_scan/analysis/read_mapping/mapping_summary_bowtie2_RdRpScan_Round1_analysis/results_mapped_only.txt", 
           delim = "\t")


read_map_df_Quaranja_clade <- read_map_df[which(read_map_df$ref_id %in% Quaranja_clade_df$contig_id & read_map_df$selected == "yes"),]



read_map_df_Quaranja_clade$RdRp_scan_map_consensus_id <- paste0(read_map_df_Quaranja_clade$ref_id, "_",
                      read_map_df_Quaranja_clade$sample_id, "_RdRpScan_Round1")

read_map_df_Quaranja_clade$contig_id <- read_map_df_Quaranja_clade$ref_id

consensus_df <- left_join(read_map_df_Quaranja_clade, Quaranja_clade_df, by  = "contig_id")



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

dir_cons1 <- "/full_path_to/wd/RdRp_scan/analysis/read_mapping/cons1_renamed/"
extracted_consensus_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/extracted_cons1_renamed/"
dir.create(extracted_consensus_directory)
extracted_orf_nt_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/extracted_cons1_orf_nt/"
dir.create(extracted_orf_nt_directory)
extracted_orf_aa_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/extracted_cons1_orf_aa/"
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

write.table(x = consensus_df, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/consensus_df.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)



#### Prepare consensus metadata in the format used for tree annotation later
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
  select(nt_seq_ID = vir_sp_seq_ID_v3, 
         tip_labels = vir_sp_seq_ID_v3,
         original_metaspades_contig_node_IDs,
         original_sample_ID, read_map_consensus_sample_ID, 
         orf_id, contig_id)

novel_seq_mdf1 <- left_join(novel_seq_mdf, mdf, by  = c("read_map_consensus_sample_ID"  = "sample_id"))


novel_seq_mdf1$host_taxon <- novel_seq_mdf1$verif_mosq_sp
novel_seq_mdf1[which(novel_seq_mdf1$SKA_SNP_dist_dend == "could not be analyzed"),]$host_taxon <- paste0(novel_seq_mdf1[which(novel_seq_mdf1$SKA_SNP_dist_dend == "could not be analyzed"),]$morph_mosq_sp, "_morph" )

novel_seq_mdf1$host_genus <- str_extract(novel_seq_mdf1$host_taxon, "^[[:alpha:]]+")

novel_seq_mdf1[which(novel_seq_mdf1$SKA_SNP_dist_dend == "could not be analyzed"),]$host_genus<- paste0(novel_seq_mdf1[which(novel_seq_mdf1$SKA_SNP_dist_dend == "could not be analyzed"),]$host_genus, "_morph" )

write.table(x = novel_seq_mdf1, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/updated_consensus_df.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)










#### REFERENCE SEQUENCES

### Use DIAMOND blastp results to get the hits
blastp_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/lineage_collapsed.DIAMOND_blastp_hits_nrcontigs_rdrp.txt.gz", 
                        delim = "\t") %>% 
  select(orf_id = qseqid, sseqid, pident, length, stitle, tax_id, superkingdom)
blastp_df1 <- blastp_df[which(blastp_df$superkingdom == "Viruses"),]

blastp_df2 <- blastp_df1[which(blastp_df1$orf_id %in% sort(unique(consensus_df$orf_id))),]


### Add potentially more distant DIAMOND blastp results from the web blastp
web_blastp_colnames <- c("qseqid", "sseqid", "pident", "length", "mismatches", 
                         "gapopens", "qstart", "qend", "sstart", "send", 
                         "evalue", "bitscore", "perc_positives")


web_blastp_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/1PN6MD51016-Alignment.txt", 
                        delim = "\t", col_names = web_blastp_colnames, comment = "#") %>% 
  select(qseqid, sseqid, pident, length)






#### Let's also add RdRp_scan sequences from the same clade
rdrp_scan_clade_accns <- str_remove(str_remove(str_split(consensus_df$ref_suggestions, ",")[[1]], "[[:print:]]+NOT_ASSIGNED_"), "\\.[[:digit:]]{1}\\_[[:print:]]+$")

refs <- sort(unique(c(blastp_df2$sseqid, rdrp_scan_clade_accns, web_blastp_df$sseqid)))

##### add sequences from publications on quaranjaviruses
################## Tree from Dudas & Batson 2022
tree_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/Dudas_2022_PB1_increasing_midpoint.nwk"
tree <- read.newick(tree_file)
tip_labels <- tree$tip.label

start <- F
tl_select <- c()
for (tl in tip_labels){
  if (start == F){
    if (tl == "'QRD99912.1|Aedes_detritus_orthomyxo-like_virus'"){
      start <- T
      tl_select <- c(tl_select, tl)
    }
  }else{
    if(tl == "'UMO75734.1|Hainan_orthomyxo-like_virus_2'"){
      start <- F
      tl_select <- c(tl_select, tl)
    }else{
      tl_select <- c(tl_select, tl)
    }
  }
}

tl_select1 <- str_extract(tl_select, "[[:alnum:]\\.\\_]+(?=\\|)") 

refs1 <- sort(unique(c(refs, tl_select1)))
refs1 <- refs1[!refs1=="NOT_ASSIGNED_Henan_2_orthomyxo_like_virus"]

##### Check on NCBI Viruses
paste0(refs1, collapse = ", ")


#### Reimport metadata
NCBI_Virus_mdf <- read_csv("/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/Quaranja_metadata_public.csv")

NCBI_Virus_mdf$nt_seq_ID <- paste0(NCBI_Virus_mdf$Accession, "_", str_replace_all(NCBI_Virus_mdf$Species, " ", "_"))
NCBI_Virus_mdf$nt_seq_ID <- paste0(NCBI_Virus_mdf$Accession, "_", str_replace_all(NCBI_Virus_mdf$Species, " ", "_"))
NCBI_Virus_mdf$host_genus <- str_extract(NCBI_Virus_mdf$Host, "^[[:alpha:]]+")
NCBI_Virus_mdf$Year <- str_extract(NCBI_Virus_mdf$Collection_Date, "^[[:digit:]]{4}")
selected_mdf <- NCBI_Virus_mdf %>% 
  select(nt_accession = Accession,
         nt_seq_ID = nt_seq_ID, 
         tip_labels = nt_seq_ID, 
         host_taxon = Host, 
         host_genus, 
         Country, 
         Location = Geo_Location, 
         Year, 
         nt_slen = Length, 
         nt_title = GenBank_Title)

selected_mdf_upd_novel <- bind_rows(selected_mdf, novel_seq_mdf1)



temp_check <- selected_mdf_upd_novel[which(is.na(selected_mdf_upd_novel$host_genus)),] ## check if indeed these are all sediment/soil
selected_mdf_upd_novel[which(is.na(selected_mdf_upd_novel$host_taxon)),]$host_taxon <- "sediment/soil"
selected_mdf_upd_novel[which(is.na(selected_mdf_upd_novel$host_genus)),]$host_genus <- "sediment/soil"


write.table(x = selected_mdf_upd_novel, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/updated_metadata_with_novel_seq.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)



######################## Dealing with CDS
####### Extract CDS from NCBI nucleotide, not NCBI Virus, because the headers are more informative
paste0(sort(selected_mdf$nt_accession), collapse = ", ")
### using this accession list to NCBI and extract CDS sequences and import them here.
#### Put the file here!
cds_fasta_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/Quaranja_clade_refs_cds_nt.txt"


### There might be more automatic solutions for CDS extractions with EMBOSS or CL implementation of entrez. But web NCBI is the fastest for now.

############################################################################################
############################################################################################
############################################################################################
cds_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/ref_cds_nt"
dir.create(cds_clade_dir)
aa_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/ref_cds_aa"
dir.create(aa_clade_dir)


cds_fasta <- seqinr::read.fasta(file = cds_fasta_file, seqtype = "DNA", whole.header = T)
names(cds_fasta)
length(selected_mdf$nt_accession)


cds_names_df <- tibble(original_cds_headers = names(cds_fasta),
                       nt_accession = str_remove(str_extract(names(cds_fasta), "(?<=^lcl\\|)[[:print:]]+(?=_cds)"), "\\.[[:digit:]]*"),
                       cds_product = str_replace_all(str_extract(names(cds_fasta), "(?<=\\[protein\\=)[[:print:]]+?(?=\\])"), " ", "_")
)


ref_df <- selected_mdf
setdiff(sort(cds_names_df$nt_accession), sort(ref_df$nt_accession))

cds_names_df1 <- left_join(cds_names_df, ref_df,  by = c( "nt_accession"= "nt_accession"))
cds_names_df1$new_short_name_cds <- cds_names_df1$nt_seq_ID
cds_names_df1$cds_status <- "OK"


genome_fraction_threshold_for_cds <- 0.7

for (cds_name in names(cds_fasta)){
  cds <- cds_fasta[cds_name]
  ## print(cds)
  cds_seq <- as.character(cds[[1]])
  
  
  ### CDS verification
  
  genome_fraction <- length(cds_seq)/as.numeric(cds_names_df1[which(cds_names_df1$original_cds_headers == cds_name), ]$nt_slen)
  
  if (genome_fraction < genome_fraction_threshold_for_cds){
    cds_names_df1[which(cds_names_df1$original_cds_headers == cds_name), ]$cds_status <- paste0("cds_below_threshold_of_genome_fraction_", 
                                                                                                genome_fraction_threshold_for_cds)
    print("cds too short")
    print(length(cds_seq))
    print(cds_names_df1[which(cds_names_df1$original_cds_headers == cds_name), ]$nt_slen)
    print(cds_names_df1[which(cds_names_df1$original_cds_headers == cds_name), ]$nt_title)
    print(cds_names_df1[which(cds_names_df1$original_cds_headers == cds_name), ]$cds_product)
    next
  }
  last_codon <- translate(seq = as.character(cds_seq[(length(cds_seq)-2):(length(cds_seq))]),
                          frame = 0, sens = "F", NAstring = "X")
  if (last_codon == "*"){
    ### remove last stop codon from cds
    cds_seq <- as.character(cds_seq[1:(length(cds_seq)-3)])
  }
  ### Translate without last codon
  cds_trans <- translate(seq = cds_seq,
                         frame = 0, sens = "F", NAstring = "X", ambiguous = T)
  
  if (length(cds_trans[cds_trans=="*"])>0){
    cds_names_df1[which(cds_names_df1$original_cds_headers == cds_name), ]$cds_status <- "stop_codon_detected"
  }
  
  cds_names_df1[which(cds_names_df1$original_cds_headers == cds_name), ]$new_short_name_cds
  
  seqinr::write.fasta(sequences = cds_seq, names = cds_names_df1[which(cds_names_df1$original_cds_headers == cds_name), ]$new_short_name_cds,
              file.out  = paste0(cds_clade_dir,"/", cds_names_df1[which(cds_names_df1$original_cds_headers == cds_name), ]$new_short_name_cds, ".fasta"))
  
  seqinr::write.fasta(sequences = cds_trans, names = cds_names_df1[which(cds_names_df1$original_cds_headers == cds_name), ]$new_short_name_cds,
              file.out  = paste0(aa_clade_dir,"/", cds_names_df1[which(cds_names_df1$original_cds_headers == cds_name), ]$new_short_name_cds, ".fasta"))
  
}


####### After this proceed to concatenating sequence datasets in Geneious, then alignment in mafft and tree in IQ-TREE



# #### After checking v2 tree, remove the arthropod-associated clade that doesn't have any mosquito associated viruses
# arthro_clade_str <- "MT153531_Coleopteran_orthomyxo-related_virus_OKIAV184
# MG972993_Photinus_pyralis_orthomyxo-like_virus_2
# MT153375_Hemipteran_orthomyxo-related_virus_OKIAV183
# MT153551_Hemipteran_orthomyxo-related_virus_OKIAV182
# MT153452_Coleopteran_orthomyxo-related_virus_OKIAV186
# KX883873_Hubei_orthomyxo-like_virus_1
# MT153423_Lepidopteran_orthomyxo-related_virus_OKIAV1731
# MT153484_Lepidopteran_orthomyxo-related_virus_OKIAV178
# MT153534_Coleopteran_orthomyxo-related_virus_OKIAV179
# MT153541_Blattodean_orthomyxo-related_virus_OKIAV181
# MW288202_Raphidiopteran_orthomyxo-related_virus_OKIAV180
# MG972988_Photinus_pyralis_orthomyxo-like_virus_1
# MT153548_Dermapteran_orthomyxo-related_virus_OKIAV170
# MW288219_Hymenopteran_orthomyxo-related_virus_OKIAV175
# MW288228_Hymenopteran_orthomyxo-related_virus_OKIAV174
# MW288220_Hymenopteran_orthomyxo-related_virus_OKIAV171
# MT153463_Phasmatodean_orthomyxo-related_virus_OKIAV172
# MT153296_Hymenopteran_orthomyxo-related_virus_OKIAV173
# KM817627_Wuhan_Mothfly_Virus
# MK026596_Old_quarry_swamp_virus
# MW033635_Soybean_thrips_quaranja-like_virus_1
# KX883882_Hubei_earwig_virus_1
# KX883858_Hubei_earwig_virus_1
# MN167480_Lestrade_virus
# MT153496_Siphonapteran_orthomyxo-related_virus_OKIAV157
# MT153412_Coleopteran_orthomyxo-related_virus_OKIAV158
# OM405131_Quaranjavirus_sp.
# MN830237_Uumaja_virus
# NC_052931_Quaranjavirus_johnstonense
# KX883844_Beihai_orthomyxo-like_virus_1
# MW310365_Bactrocera_correcta_orthomyxo-like_virus_isolate_Bl
# MW310389_Bactrocera_correcta_orthomyxo-like_virus_isolate_Bz
# MT153421_Phasmatodean_orthomyxo-related_virus_OKIAV167
# MT153485_Hemipteran_orthomyxo-related_virus_OKIAV188
# ON872578_Bat_faecal_associated_orthomyxo-like_virus_1
# MW784037_Hainan_orthomyxo-like_virus_2
# MW784037_Hainan_orthomyxo-like_virus_2_0001
# MT153419_Hemipteran_orthomyxo-related_virus_OKIAV191
# OM650310_Halyomorpha_halys_orthomyxo-like_virus_1
# MW288163_Neuropteran_orthomyxo-related_virus_OKIAV190
# "
# arthro_clade <- str_split(arthro_clade_str, "\\s")[[1]]
# arthro_clade <- arthro_clade[which(!arthro_clade == "")]
# sort(unique(selected_mdf_upd_novel[which(selected_mdf_upd_novel$nt_seq_ID %in% arthro_clade),]$host_genus))
# selected_mdf_upd_novel[which(selected_mdf_upd_novel$nt_seq_ID %in% arthro_clade & selected_mdf_upd_novel$host_taxon == "sediment/soil"),]
# #### remove from sequence list
# keep_for_v3 <- selected_mdf_upd_novel[which(!selected_mdf_upd_novel$nt_seq_ID %in% arthro_clade),]$nt_seq_ID
# 
# ##remove two more sequences because they are short
# 
# keep_for_v3 <- keep_for_v3[which(!keep_for_v3 %in% c("MW784007_Beijing_sediment_orthomyxo-like_virus", 
#                                      "MW784011_Xinjiang_orthomyxo-like_virus_3"))]
# 
# v2_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/Quaranja_clade_aa_v2.fasta",
#                                 seqtype = "AA")
# 
# names(v2_aa_seq)
# 
# v3_aa_seq <- v2_aa_seq[keep_for_v3]
# 
# names(v3_aa_seq)
# 
# seqinr::write.fasta(sequences = v3_aa_seq , names =names(v3_aa_seq), file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/Quaranja_clade_aa_v3.fasta")
# ##keep the closest sequence as an outgroup
# new_outg <- "MW288163_Neuropteran_orthomyxo-related_virus_OKIAV190"
# new_outg_seq <- v2_aa_seq[new_outg]
# seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/MW288163_Neuropteran_orthomyxo-related_virus_OKIAV190.fasta")





# ####### Zooming in on different clades
# 
# #### v31 tree
# clade_v31_str <- "OL700087_Wuhan_Mosquito_Virus_4
# KM817623_Wuhan_Mosquito_Virus_4
# OL700088_Wuhan_Mosquito_Virus_4
# MW520421_Culex_modestus_orthomyxo-like_virus
# MN513372_Culex_orthomyxo-like_virus
# MK440637_Wuhan_Mosquito_Virus_4
# Pageromars_virus_F1413C102F1245
# Pageromars_virus_F1413C102F1413
# MW317143_Orthomyxoviridae_sp."
# clade_v31 <- str_split(clade_v31_str, "\\s")[[1]]
# clade_v31 <- clade_v31[which(!clade_v31 == "")]
# keep_for_v31 <- clade_v31
# v3_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/Quaranja_clade_aa_v3.fasta",
#                                 seqtype = "AA")
# v31_aa_seq <- v3_aa_seq[keep_for_v31]
# names(v31_aa_seq)
# seqinr::write.fasta(sequences = v31_aa_seq , names =names(v31_aa_seq), file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/Quaranja_clade_aa_v31.fasta")
# ##keep the closest sequence as an outgroup
# new_outg <- "MN053836_Guadeloupe_mosquito_quaranja-like_virus_3"
# new_outg_seq <- v2_aa_seq[new_outg]
# seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/MN053836_Guadeloupe_mosquito_quaranja-like_virus_3.fasta")
# 
# 
# #### v32 tree
# clade_v32_str <- "MW434422_Wuhan_Mosquito_Virus_6
# MW434426_Wuhan_Mosquito_Virus_6
# OM817558_Wuhan_Mosquito_Virus_6
# MW434432_Wuhan_Mosquito_Virus_6
# MW434423_Wuhan_Mosquito_Virus_6
# MW434430_Wuhan_Mosquito_Virus_6
# MW434425_Wuhan_Mosquito_Virus_6
# MW434421_Wuhan_Mosquito_Virus_6
# MW434431_Wuhan_Mosquito_Virus_6
# MW434428_Wuhan_Mosquito_Virus_6
# MW434429_Wuhan_Mosquito_Virus_6
# MW434424_Wuhan_Mosquito_Virus_6
# MW520424_Culex_pipiens_orthomyxo-like_virus
# KM817625_Wuhan_Mosquito_Virus_6
# Wuhan_Mosquito_Virus_6_F1115C44F1120
# MW452277_Wuhan_Mosquito_Virus_6
# MF176271_Wuhan_Mosquito_Virus_6
# MF176360_Wuhan_Mosquito_Virus_6
# MF176320_Wuhan_Mosquito_Virus_6
# MF176249_Wuhan_Mosquito_Virus_6
# MK440648_Wuhan_Mosquito_Virus_6
# Wuhan_Mosquito_Virus_6_F1115C44F1115
# Wuhan_Mosquito_Virus_6_F1115C44F1576
# Wuhan_Mosquito_Virus_6_F1115C44F1422
# OL700090_Wuhan_Mosquito_Virus_6"
# clade_v32 <- str_split(clade_v32_str, "\\s")[[1]]
# clade_v32 <- clade_v32[which(!clade_v32 == "")]
# keep_for_v32 <- clade_v32
# v3_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/Quaranja_clade_aa_v3.fasta",
#                                 seqtype = "AA")
# v32_aa_seq <- v3_aa_seq[keep_for_v32]
# seqinr::write.fasta(sequences = v32_aa_seq , names =names(v32_aa_seq), file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/Quaranja_clade_aa_v32.fasta")
# ##keep the closest sequence as an outgroup
# new_outg <- "OL700154_XiangYun_orthomyxo-like_virus_2"
# new_outg_seq <- v2_aa_seq[new_outg]
# seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/OL700154_XiangYun_orthomyxo-like_virus_2.fasta")
# 
# 
# #### v33 tree
# clade_v33_str <- "Guadeloupe_mosquito_quaranja_like_virus_1_F1058C389F1059
# Guadeloupe_mosquito_quaranja_like_virus_1_F1058C389F1062
# Guadeloupe_mosquito_quaranja_like_virus_1_F1058C389F1068
# Guadeloupe_mosquito_quaranja_like_virus_1_F1058C389F1058
# Guadeloupe_mosquito_quaranja_like_virus_1_F1058C389F1065
# Guadeloupe_mosquito_quaranja_like_virus_1_F1058C389F1067
# OL343761_Palmetto_orthomyxo-like_virus
# MW434249_Guadeloupe_mosquito_quaranja-like_virus_1
# MN053838_Guadeloupe_mosquito_quaranja-like_virus_1
# Panessemat_virus_F1398C269F1230
# Panessemat_virus_F1398C269F1540
# Panessemat_virus_F1398C269F1398
# Panessemat_virus_F1398C269F1387
# Panessemat_virus_F1398C269F1396
# Panessemat_virus_F1398C269F1551
# MW520417_Aedes_detritus_orthomyxo-like_virus
# KX898491_Whidbey_virus"
# clade_v33 <- str_split(clade_v33_str, "\\s")[[1]]
# clade_v33 <- clade_v33[which(!clade_v33 == "")]
# keep_for_v33 <- clade_v33
# v3_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/Quaranja_clade_aa_v3.fasta",
#                                 seqtype = "AA")
# v33_aa_seq <- v3_aa_seq[keep_for_v33]
# seqinr::write.fasta(sequences = v33_aa_seq , names =names(v33_aa_seq), file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/Quaranja_clade_aa_v33.fasta")
# ##keep the closest sequence as an outgroup
# new_outg <- "BK059432_Aedes_orthomyxo-like_virus_2"
# new_outg_seq <- v2_aa_seq[new_outg]
# seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/BK059432_Aedes_orthomyxo-like_virus_2.fasta")
# 




####### Update after including previously missing novel sequences


#### After checking v4 tree, remove the arthropod-associated clade that doesn't have any mosquito associated viruses
select_clade_str <- "MF176249_Wuhan_Mosquito_Virus_6
OL700153_XiangYun_orthomyxo-like_virus_1
MW784029_Heilongjiang_sediment_orthomyxo-like_virus
MW784026_Xinjiang_sediment_orthomyxo-like_virus_2
MW784025_Henan_sediment_orthomyxo-like_virus
MW784032_Yunnan_sediment_orthomyxo-like_virus
OL700156_XiangYun_orthomyxo-like_virus_4
MW784015_Xinjiang_sediment_orthomyxo-like_virus_4
OL700157_XiangYun_orthomyxo-like_virus_5
MW434314_Usinis_virus
Panessemat_virus_F1398C269F1540
Panessemat_virus_F1398C269F1230
Panessemat_virus_F1398C269F1398
Panessemat_virus_F1398C269F1387
Panessemat_virus_F1398C269F1396
MF176271_Wuhan_Mosquito_Virus_6
OL700087_Wuhan_Mosquito_Virus_4
KM817623_Wuhan_Mosquito_Virus_4
Guadeloupe_mosquito_quaranja_like_virus_1_F1058C389F1059
Guadeloupe_mosquito_quaranja_like_virus_1_F1058C389F1062
Guadeloupe_mosquito_quaranja_like_virus_1_F1058C389F1068
Guadeloupe_mosquito_quaranja_like_virus_1_F1058C389F1058
Guadeloupe_mosquito_quaranja_like_virus_1_F1058C389F1067
Guadeloupe_mosquito_quaranja_like_virus_1_F1058C389F1065
MW434249_Guadeloupe_mosquito_quaranja-like_virus_1
OL343761_Palmetto_orthomyxo-like_virus
MN053838_Guadeloupe_mosquito_quaranja-like_virus_1
Panessemat_virus_F1398C269F1551
MW520417_Aedes_detritus_orthomyxo-like_virus
KX898491_Whidbey_virus
BK059432_Aedes_orthomyxo-like_virus_2
MN053837_Guadeloupe_mosquito_quaranja-like_virus_2
MF176325_Aedes_alboannulatus_orthomyxo-like_virus
MW452276_Wuhan_Mosquito_Virus_5
KM817624_Wuhan_Mosquito_Virus_5
OL700089_Wuhan_Mosquito_Virus_5
KM817622_Wuhan_Mosquito_Virus_3
OL700086_Wuhan_Mosquito_Virus_3
MW434228_Astopletus_virus
MW434229_Astopletus_virus
ON059763_Byreldi_virus
MF176320_Wuhan_Mosquito_Virus_6
MF176360_Wuhan_Mosquito_Virus_6
MK440648_Wuhan_Mosquito_Virus_6
Wuhan_Mosquito_Virus_6_F1115C44F1120
MT129718_Splett_orthomyxo-like_virus
MT129719_Splett_orthomyxo-like_virus
KM817615_Jingshan_Fly_Virus_1
MT153448_Dipteran_orthomyxo-related_virus_OKIAV195
Wuhan_Mosquito_Virus_6_F1115C44F1115
Wuhan_Mosquito_Virus_6_F1115C44F1576
Wuhan_Mosquito_Virus_6_F1115C44F1422
OL700090_Wuhan_Mosquito_Virus_6
MW520424_Culex_pipiens_orthomyxo-like_virus
MW452277_Wuhan_Mosquito_Virus_6
KM817625_Wuhan_Mosquito_Virus_6
MW434424_Wuhan_Mosquito_Virus_6
MW434422_Wuhan_Mosquito_Virus_6
MW434426_Wuhan_Mosquito_Virus_6
MW434431_Wuhan_Mosquito_Virus_6
OM817558_Wuhan_Mosquito_Virus_6
MW434423_Wuhan_Mosquito_Virus_6
MW434430_Wuhan_Mosquito_Virus_6
MW434425_Wuhan_Mosquito_Virus_6
MW434421_Wuhan_Mosquito_Virus_6
MW434428_Wuhan_Mosquito_Virus_6
MW434432_Wuhan_Mosquito_Virus_6
MW434429_Wuhan_Mosquito_Virus_6
OL700154_XiangYun_orthomyxo-like_virus_2
MW434320_Usinis_virus
Usinis_virus_F1083C108F1515
Usinis_virus_F1083C108F1083
Spoisillett_virus_F1019C753F1019
KM817620_Wuhan_Louse_Fly_Virus_3
KM817621_Wuhan_Louse_Fly_Virus_4
MN053836_Guadeloupe_mosquito_quaranja-like_virus_3
KM817619_Shuangao_Insect_Virus_4
MW199253_Orthomyxoviridae_sp.
MT153481_Dipteran_orthomyxo-related_virus_OKIAV193
MN513372_Culex_orthomyxo-like_virus
ON860461_Tierp_virus
MW520421_Culex_modestus_orthomyxo-like_virus
MK440637_Wuhan_Mosquito_Virus_4
OL700088_Wuhan_Mosquito_Virus_4
Pageromars_virus_F1413C102F1245
Pageromars_virus_F1413C102F1413
MW317143_Orthomyxoviridae_sp.
ON860462_Husby_virus
ON059769_Byreska_virus
KM817617_Sanxia_Water_Strider_Virus_3
MW784019_Xinjiang_sediment_orthomyxo-like_virus_5
"
select_clade <- str_split(select_clade_str, "\\s")[[1]]
select_clade <- select_clade[which(!select_clade == "")]
keep_for_v5 <- selected_mdf_upd_novel[which(selected_mdf_upd_novel$nt_seq_ID %in% select_clade),]$nt_seq_ID

##remove two more sequences because they are short

keep_for_v5 <- keep_for_v5[which(!keep_for_v5 %in% c("MW784007_Beijing_sediment_orthomyxo-like_virus", 
                                                     "MW784011_Xinjiang_orthomyxo-like_virus_3"))]

v4_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/Quaranja_clade_aa_v4.fasta",
                                seqtype = "AA")

names(v4_aa_seq)

v5_aa_seq <- v4_aa_seq[keep_for_v5]

names(v5_aa_seq)

seqinr::write.fasta(sequences = v5_aa_seq , names =names(v5_aa_seq), file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/Quaranja_clade_aa_v5.fasta")
##keep the closest sequence as an outgroup
new_outg <- "MW784037_Hainan_orthomyxo-like_virus_2"
new_outg_seq <- v4_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/aln/aa_seq/MW784037_Hainan_orthomyxo-like_virus_2.fasta")











