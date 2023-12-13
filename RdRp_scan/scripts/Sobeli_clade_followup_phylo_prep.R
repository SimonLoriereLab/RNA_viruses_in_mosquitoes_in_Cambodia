library(tidyverse)
library(seqinr)
library(treeio)

Sobeli_clade_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/Sobeli_clade_novel_sequence_data.txt")

### select consensus sequences based on mapping coverage
### first round of read mapping results
read_map_df <-  read_delim("/full_path_to/wd/RdRp_scan/analysis/read_mapping/mapping_summary_bowtie2_RdRpScan_Round1_analysis/results_mapped_only.txt", 
           delim = "\t")


read_map_df_Sobeli_clade <- read_map_df[which(read_map_df$ref_id %in% Sobeli_clade_df$contig_id & read_map_df$selected == "yes"),]



read_map_df_Sobeli_clade$RdRp_scan_map_consensus_id <- paste0(read_map_df_Sobeli_clade$ref_id, "_",
                      read_map_df_Sobeli_clade$sample_id, "_RdRpScan_Round1")

read_map_df_Sobeli_clade$contig_id <- read_map_df_Sobeli_clade$ref_id

consensus_df <- left_join(read_map_df_Sobeli_clade, Sobeli_clade_df, by  = "contig_id")



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
extracted_consensus_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/extracted_cons1_renamed/"
dir.create(extracted_consensus_directory)
extracted_orf_nt_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/extracted_cons1_orf_nt/"
dir.create(extracted_orf_nt_directory)
extracted_orf_aa_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/extracted_cons1_orf_aa/"
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

write.table(x = consensus_df, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/consensus_df.txt", 
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




write.table(x = novel_seq_mdf1, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/updated_consensus_df.txt", 
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


web_blastp_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/TGK1B6PU01N-Alignment.txt", 
                        delim = "\t", col_names = web_blastp_colnames, comment = "#") %>% 
  select(qseqid, sseqid, pident, length)


#### Let's also add RdRp_scan sequences from the same clade
rdrp_scan_clade_accns <- str_remove(str_remove(str_split(consensus_df$ref_suggestions, ",")[[1]], "[[:print:]]+Sobelivirales_"), "\\.[[:digit:]]{1}\\_[[:print:]]+$")

refs <- sort(unique(c(blastp_df2$sseqid, rdrp_scan_clade_accns, web_blastp_df$sseqid)))
refs <- str_remove(refs, "\\.[[:digit:]]$")





##### Check on NCBI Viruses
paste0(refs, collapse = ", ")


#### Reimport metadata
NCBI_Virus_mdf <- read_csv("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/Sobeli_metadata_public.csv")


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


write.table(x = selected_mdf_upd_novel, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/updated_metadata_with_novel_seq.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)



######################## Dealing with CDS
####### Extract CDS from NCBI nucleotide, not NCBI Virus, because the headers are more informative
paste0(sort(selected_mdf$nt_accession), collapse = ",")

####### Here, because there could be multiple CDS, use protein accessions to extract CDS
####### Correspondence between the prot and nt accessions is encoded in standard CDS names when extraction from NCBI Protein same as from NCBI Nucleotide

paste0(sort(refs), collapse = ",")


### using this accession list to NCBI and extract CDS sequences and import them here.
#### Put the file here!
cds_fasta_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/Sobeli_clade_refs_cds_nt.txt"


### There might be more automatic solutions for CDS extractions with EMBOSS or CL implementation of entrez. But web NCBI is the fastest for now.

############################################################################################
############################################################################################
############################################################################################
cds_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/ref_cds_nt"
dir.create(cds_clade_dir)
aa_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/ref_cds_aa"
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

### remove this endogenous element "MN272068"
cds_names_df <- cds_names_df[which(!cds_names_df$nt_accession == "MN272068"),]


### remove capsid proteins
sort(unique(cds_names_df$cds_product))

cds_names_df <- cds_names_df[which(!cds_names_df$cds_product == "capsid_protein"),]

cds_names_df1 <- left_join(cds_names_df, ref_df,  by = c( "nt_accession"= "nt_accession"))
cds_names_df1$new_short_name_cds <- cds_names_df1$nt_seq_ID
cds_names_df1$cds_status <- "OK"


cds_fasta <- cds_fasta[cds_names_df1$original_cds_headers]


genome_fraction_threshold_for_cds <- 0 ### no threshold, because it's not a polyprotein, cds can be quite small

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

##### In Geneious correct several Humaita Tubiacanga virus ORFs.


####### Zooming in on different clades
#### v21 tree
clade_v21_str <- "Humaita_Tubiacanga_virus_F1058C249F1058
Humaita_Tubiacanga_virus_F1058C249F1336
MW199222_Virus_sp.
MZ396017_Solemoviridae_sp.
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1021
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1381
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1531
KX882872_Hubei_sobemo-like_virus_41
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1532
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1525
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1222
Humaita_Tubiacanga_virus_F1058C249F1002
MN053809_Humaita-Tubiacanga_virus
MN053819_Humaita-Tubiacanga_virus
Humaita_Tubiacanga_virus_F1058C249F1086
Humaita_Tubiacanga_virus_F1058C249F1063
ON125137_Freshwater_macrophyte_associated_sobeli-like_virus_2
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1219
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1383
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1019
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1527
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1536
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1533
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1537
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1094
NC_032253_Hubei_sobemo-like_virus_41
MW509554_Laem_Singh_virus
NC_033300_Wenzhou_shrimp_virus_9
MN053817_Humaita-Tubiacanga_virus
MT025121_Khabarov_virus
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1528
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1223
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1215
Hubei_sobemo_like_virus_41_subsp1_F1020C17F1020
Culex_luteo_like_virus_F1039C258F1039
Culex_luteo_like_virus_F1039C258F1563
MF176307_Culex_luteo-like_virus
MW520388_Culex_impudicus_luteo-like_virus
MT568531_Pyongtaek_Culex_luteo-like_virus
MK440657_Berrek_virus
MW434790_Culex_mosquito_virus_6
MW434802_Culex_mosquito_virus_6
MH188030_Culex_mosquito_virus_6
MW434798_Culex_mosquito_virus_6
MW434800_Culex_mosquito_virus_6
OM817536_Culex_mosquito_virus_6
OM743945_Pine_Lake_virus
Broome_luteo_like_virus_1_F1261C1F1261
MT498822_Broome_luteo-like_virus_1
MN661095_Atrato_Sobemo-like_virus_4
MT568532_Pyongtaek_Culex_Solemovirus
MN661091_Atrato_Sobemo-like_virus_3
MN661093_Atrato_Sobemo-like_virus_3
Humaita_Tubiacanga_virus_F1058C249F1062
MN053823_Humaita-Tubiacanga_virus
MW650676_Humaita-Tubiacanga_virus
MW773212_Humaita-Tubiacanga_virus
Humaita_Tubiacanga_virus_F1058C249F1060
Humaita_Tubiacanga_virus_F1058C249F1001
Humaita_Tubiacanga_virus_F1058C249F1335
Humaita_Tubiacanga_virus_F1058C249F1067
Humaita_Tubiacanga_virus_F1058C249F1503
Humaita_Tubiacanga_virus_F1058C249F1065
Humaita_Tubiacanga_virus_F1058C249F1059
Humaita_Tubiacanga_virus_F1058C249F1064
Humaita_Tubiacanga_virus_F1058C249F1339
MW199221_Virus_sp.
MW199223_Virus_sp.
MW239244_Riboviria_sp.
MZ443605_Leuven_Luteo-like_virus_1
MW239226_Riboviria_sp.
MN831442_Cycas_revoluta_sobemo-like_virus
MN831444_Teucrium_fruticans_sobemo-like_virus
MW239188_Riboviria_sp.
MZ484471_Solemoviridae_sp.
NC_033452_Wuchan_romanomermis_nematode_virus_3
NC_032865_Beihai_sobemo-like_virus_27
MW239454_Riboviria_sp.
MW239461_Riboviria_sp.
MW239404_Riboviria_sp.
MW239435_Riboviria_sp.
NC_032251_Hubei_sobemo-like_virus_40
KX882864_Hubei_sobemo-like_virus_40
NC_033089_Shuangao_sobemo-like_virus_2
MW239382_Riboviria_sp.
MF189984_Barns_Ness_beadlet_anemone_sobemo-like_virus_1
"
clade_v21 <- str_split(clade_v21_str, "\\s")[[1]]
clade_v21 <- clade_v21[which(!clade_v21 == "")]
keep_for_v21 <- clade_v21
v2_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v2.fasta",
                                seqtype = "AA")
v21_aa_seq <- v2_aa_seq[keep_for_v21]
names(v21_aa_seq)
seqinr::write.fasta(sequences = v21_aa_seq , names =names(v21_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v21.fasta")
##keep the closest sequence as an outgroup
new_outg <- "NC_032193_Beihai_sobemo-like_virus_26"
new_outg_seq <- v2_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))




#### v211 tree
clade_v211_str <- "Hubei_sobemo_like_virus_41_subsp2_F1381C255F1019
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1536
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1533
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1094
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1527
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1537
NC_032253_Hubei_sobemo-like_virus_41
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1528
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1223
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1021
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1531
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1532
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1219
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1525
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1383
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1215
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1381
Hubei_sobemo_like_virus_41_subsp2_F1381C255F1222
KX882872_Hubei_sobemo-like_virus_41
Hubei_sobemo_like_virus_41_subsp1_F1020C17F1020
Culex_luteo_like_virus_F1039C258F1039
Culex_luteo_like_virus_F1039C258F1563
MF176307_Culex_luteo-like_virus
MT568531_Pyongtaek_Culex_luteo-like_virus
MW520388_Culex_impudicus_luteo-like_virus
MK440657_Berrek_virus
MW434802_Culex_mosquito_virus_6
MH188030_Culex_mosquito_virus_6
MW434800_Culex_mosquito_virus_6
MW434790_Culex_mosquito_virus_6
MW434798_Culex_mosquito_virus_6
OM817536_Culex_mosquito_virus_6
OM743945_Pine_Lake_virus
Broome_luteo_like_virus_1_F1261C1F1261
MT498822_Broome_luteo-like_virus_1
MN661095_Atrato_Sobemo-like_virus_4
MT568532_Pyongtaek_Culex_Solemovirus
MN661091_Atrato_Sobemo-like_virus_3
MN661093_Atrato_Sobemo-like_virus_3
MN053809_Humaita-Tubiacanga_virus
Humaita_Tubiacanga_virus_F1058C249F1063
Humaita_Tubiacanga_virus_F1058C249F1062
MN053819_Humaita-Tubiacanga_virus
MN053817_Humaita-Tubiacanga_virus
Humaita_Tubiacanga_virus_F1058C249F1086
Humaita_Tubiacanga_virus_F1058C249F1002
MN053823_Humaita-Tubiacanga_virus
Humaita_Tubiacanga_virus_F1058C249F1058
Humaita_Tubiacanga_virus_F1058C249F1059
Humaita_Tubiacanga_virus_F1058C249F1064
MW650676_Humaita-Tubiacanga_virus
Humaita_Tubiacanga_virus_F1058C249F1001
Humaita_Tubiacanga_virus_F1058C249F1339
Humaita_Tubiacanga_virus_F1058C249F1067
Humaita_Tubiacanga_virus_F1058C249F1335
Humaita_Tubiacanga_virus_F1058C249F1336
Humaita_Tubiacanga_virus_F1058C249F1503
Humaita_Tubiacanga_virus_F1058C249F1060
MW773212_Humaita-Tubiacanga_virus
Humaita_Tubiacanga_virus_F1058C249F1065
"
clade_v211 <- str_split(clade_v211_str, "\\s")[[1]]
clade_v211 <- clade_v211[which(!clade_v211 == "")]
keep_for_v211 <- clade_v211
v2_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v2.fasta",
                                seqtype = "AA")
v211_aa_seq <- v2_aa_seq[keep_for_v211]
names(v211_aa_seq)
seqinr::write.fasta(sequences = v211_aa_seq , names =names(v211_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v211.fasta")
##keep the closest sequence as an outgroup
new_outg <- "MW199223_Virus_sp."
new_outg_seq <- v2_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))





####### Zooming in on different clades
#### v22 tree
clade_v22_str <- "MZ396021_Solemoviridae_sp.
NC_033713_Wuhan_insect_virus_34
MW199225_Virus_sp.
MZ395979_Solemoviridae_sp.
MW199220_Virus_sp.
MZ396026_Solemoviridae_sp.
MW199224_Virus_sp.
NC_033489_Wuhan_house_centipede_virus_5
MW239390_Riboviria_sp.
MZ443595_Leuven_Sobemo-like_virus_2
MT138167_Solemoviridae_sp.
MT138168_Solemoviridae_sp.
NC_032233_Hubei_sobemo-like_virus_49
MW826514_Solemoviridae_sp.
MN627439_Erysiphe_necator_associated_sobemo-like_virus_3
MW239343_Riboviria_sp.
MZ375169_Riboviria_sp.
NC_032232_Hubei_sobemo-like_virus_48
KU754509_La_Tardoire_virus
MF893251_Medway_virus
MW239498_Riboviria_sp.
MN661085_Atrato_Sobemo-like_virus_1
MN661087_Atrato_Sobemo-like_virus_1
MW520380_Evros_sobemo-like_virus
MN609861_Leveillula_taurica_associated_sobemo-like_virus_1
MW239506_Riboviria_sp.
OL700122_XiangYun_luteo-sobemo-like_virus_1
OL700123_XiangYun_luteo-sobemo-like_virus_1
MW520373_Thassos_sobemo-like_virus
OL700124_XiangYun_luteo-sobemo-like_virus_1
MZ771233_Solemoviridae_sp.
MW239189_Riboviria_sp.
NC_033472_Wuhan_insect_virus_17
NC_032191_Hubei_sobemo-like_virus_47
NC_032195_Hubei_sobemo-like_virus_46
Poterfers_virus_F1276C48F1276
MT138166_Solemoviridae_sp.
MW239503_Riboviria_sp.
MW239215_Riboviria_sp.
MT138164_Solemoviridae_sp.
NC_032238_Hubei_sobemo-like_virus_45
NC_032213_Hubei_sobemo-like_virus_44
MW239281_Riboviria_sp.
MN661089_Atrato_Sobemo-like_virus_2
MW239471_Riboviria_sp.
MN725049_Frankliniella_occidentalis_associated_sobemo-like_virus_1"
clade_v22 <- str_split(clade_v22_str, "\\s")[[1]]
clade_v22 <- clade_v22[which(!clade_v22 == "")]
keep_for_v22 <- clade_v22
v2_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v2.fasta",
                                seqtype = "AA")
v22_aa_seq <- v2_aa_seq[keep_for_v22]
names(v22_aa_seq)
seqinr::write.fasta(sequences = v22_aa_seq , names =names(v22_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v22.fasta")
##keep the closest sequence as an outgroup
new_outg <- "MW239542_Riboviria_sp."
new_outg_seq <- v2_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))







####### Zooming in on different clades
####v23 tree


clade_v23_str <- "MW434812_Guadeloupe_mosquito_virus
MW434813_Guadeloupe_mosquito_virus
MK285337_Renna_virus
MT361065_Kwale_mosquito_virus
MN053789_Guadeloupe_mosquito_virus
MN053791_Guadeloupe_mosquito_virus
MN053787_Guadeloupe_mosquito_virus
Colokingly_virus_F1527C217F1019
Colokingly_virus_F1527C217F1100
Colokingly_virus_F1527C217F1382
Colokingly_virus_F1527C217F1527
Colokingly_virus_F1527C217F1534
Colokingly_virus_F1527C217F1094
Colokingly_virus_F1527C217F1021
Colokingly_virus_F1527C217F1211
Colokingly_virus_F1527C217F1381
Colokingly_virus_F1527C217F1096
MW030488_Hubei_mosquito_virus_2
MW434902_Marma_virus
Peninsere_virus_F1127C300F1056
Pirounced_virus_F1032C17F1316
Pirounced_virus_F1032C17F1056
Pirounced_virus_F1032C17F1585
Pirounced_virus_F1032C17F1271
Pirounced_virus_F1032C17F1305
Pirounced_virus_F1032C17F1032
Peninsere_virus_F1127C300F1032
Pirounced_virus_F1032C17F1266
Pirounced_virus_F1032C17F1034
Peninsere_virus_F1127C300F1461
Pirounced_virus_F1032C17F1607
Pirounced_virus_F1032C17F1604
Pirounced_virus_F1032C17F1461
Pirounced_virus_F1032C17F1325
NC_033305_Hubei_mosquito_virus_2
Comfornsto_virus_F1216C3F1183
MW030472_Hubei_mosquito_virus_2
MW239353_Riboviria_sp.
Peninsere_virus_F1127C300F1127
Pirounced_virus_F1032C17F1127
Peninsere_virus_F1127C300F1156
Peninsere_virus_F1127C300F1132
Peninsere_virus_F1127C300F1450
MW030476_Hubei_mosquito_virus_2
MW030482_Hubei_mosquito_virus_2
KX882764_Hubei_mosquito_virus_2
MW030470_Hubei_mosquito_virus_2
MW030484_Hubei_mosquito_virus_2
MW030486_Hubei_mosquito_virus_2
MW030478_Hubei_mosquito_virus_2
MW030480_Hubei_mosquito_virus_2
LC513833_Culex_inatomii_luteo-like_virus
MW434875_Marma_virus
MW434863_Marma_virus
MW434872_Marma_virus
MK440646_Marma_virus
MT361063_Kisumu_mosquito_virus
Guangzhou_sobemo_like_virus_F1517C31F1016
Guangzhou_sobemo_like_virus_F1517C31F1517
Guangzhou_sobemo_like_virus_F1517C31F1083
Guangzhou_sobemo_like_virus_F1517C31F1351
Guangzhou_sobemo_like_virus_F1517C31F1078
Guangzhou_sobemo_like_virus_F1517C31F1339
Guangzhou_sobemo_like_virus_F1517C31F1071
MT361061_Guangzhou_sobemo-like_virus
Guangzhou_sobemo_like_virus_F1517C31F1171
Guangzhou_sobemo_like_virus_F1517C31F1073
Guangzhou_sobemo_like_virus_F1517C31F1352
Guangzhou_sobemo_like_virus_F1517C31F1070
MT361053_Guangzhou_sobemo-like_virus
MT361055_Guangzhou_sobemo-like_virus
MT096519_Wenzhou_sobemo-like_virus_4
NC_033138_Wenzhou_sobemo-like_virus_4
MZ556263_Sichuan_mosquito_sobemo-like_virus
MT361059_Guangzhou_sobemo-like_virus
Guangzhou_sobemo_like_virus_F1517C31F1017
MN053797_Guadeloupe_mosquito_virus
MN053803_Guadeloupe_mosquito_virus
MH703049_Yongsan_sobemo-like_virus_1
MW310331_Zeugodacus_cucurbitae_sobemovirus_isolate_Bc
MW310351_Zeugodacus_cucurbitae_sobemovirus
MW199227_Virus_sp.
NC_032184_Hubei_sobemo-like_virus_9
MW239337_Riboviria_sp.
MW039355_Soybean_thrips_sobemo-like_virus_9
NC_032199_Hubei_sobemo-like_virus_10
MW239384_Riboviria_sp.
NC_032236_Hubei_sobemo-like_virus_8
MN661101_Atrato_Sobemo-like_virus_6
MW251326_Tartas_insect-associated_virus
MN551127_Plasmopara_viticola_lesion_associated_sobemo-like_1
Comfornsto_virus_F1216C3F1216
MN661099_Atrato_Sobemo-like_virus_5
MN661107_Atrato_Sobemo-like_virus_5
MN661097_Atrato_Sobemo-like_virus_5
MT896210_Aedes_sobemo-like_virus
MT896211_Aedes_sobemo-like_virus
MW434821_Kellev_virus
MZ375201_Riboviria_sp.
MW239512_Riboviria_sp.
MZ375164_Riboviria_sp.
MW239414_Riboviria_sp.
MZ443571_Vespula_vulgaris_Sobemo-like_virus_1
MW239377_Riboviria_sp.
MT138162_Solemoviridae_sp.
MK026570_Blue_fish_point_virus
MZ443630_Nelson_Sobemo-like_virus_1
MW239270_Riboviria_sp.
"
clade_v23 <- str_split(clade_v23_str, "\\s")[[1]]
clade_v23 <- clade_v23[which(!clade_v23 == "")]
keep_for_v23 <- clade_v23



v2_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v2.fasta",
                                seqtype = "AA")
v23_aa_seq <- v2_aa_seq[keep_for_v23]
names(v23_aa_seq)
seqinr::write.fasta(sequences = v23_aa_seq , names =names(v23_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v23.fasta")
##keep the closest sequence as an outgroup
new_outg <- "NC_033279_Sanxia_water_strider_virus_10"
new_outg_seq <- v2_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))


####### Zooming in on different clades
#### v231 tree
clade_v231_str <- "Peninsere_virus_F1127C300F1056
Pirounced_virus_F1032C17F1056
Pirounced_virus_F1032C17F1032
Pirounced_virus_F1032C17F1271
Pirounced_virus_F1032C17F1305
Pirounced_virus_F1032C17F1585
Pirounced_virus_F1032C17F1316
Peninsere_virus_F1127C300F1032
Pirounced_virus_F1032C17F1266
Pirounced_virus_F1032C17F1034
Pirounced_virus_F1032C17F1607
Peninsere_virus_F1127C300F1461
Pirounced_virus_F1032C17F1461
Pirounced_virus_F1032C17F1604
Pirounced_virus_F1032C17F1325
LC513833_Culex_inatomii_luteo-like_virus
Peninsere_virus_F1127C300F1127
Peninsere_virus_F1127C300F1132
Peninsere_virus_F1127C300F1156
Pirounced_virus_F1032C17F1127
Peninsere_virus_F1127C300F1450
KX882764_Hubei_mosquito_virus_2
MW030476_Hubei_mosquito_virus_2
MW030482_Hubei_mosquito_virus_2
MW030472_Hubei_mosquito_virus_2
MW030484_Hubei_mosquito_virus_2
MW030478_Hubei_mosquito_virus_2
MW030488_Hubei_mosquito_virus_2
MW030486_Hubei_mosquito_virus_2
MW030480_Hubei_mosquito_virus_2
MW030470_Hubei_mosquito_virus_2
MW434902_Marma_virus
MW434875_Marma_virus
MW434863_Marma_virus
MW434872_Marma_virus
MK440646_Marma_virus
NC_033305_Hubei_mosquito_virus_2
MT361063_Kisumu_mosquito_virus
MN661099_Atrato_Sobemo-like_virus_5
MN661097_Atrato_Sobemo-like_virus_5
MN661107_Atrato_Sobemo-like_virus_5
MT896210_Aedes_sobemo-like_virus
MT896211_Aedes_sobemo-like_virus
MH703049_Yongsan_sobemo-like_virus_1
Guangzhou_sobemo_like_virus_F1517C31F1016
Guangzhou_sobemo_like_virus_F1517C31F1071
Guangzhou_sobemo_like_virus_F1517C31F1083
Guangzhou_sobemo_like_virus_F1517C31F1078
Guangzhou_sobemo_like_virus_F1517C31F1351
Guangzhou_sobemo_like_virus_F1517C31F1339
Guangzhou_sobemo_like_virus_F1517C31F1517
MT361061_Guangzhou_sobemo-like_virus
Guangzhou_sobemo_like_virus_F1517C31F1171
Guangzhou_sobemo_like_virus_F1517C31F1070
Guangzhou_sobemo_like_virus_F1517C31F1073
Guangzhou_sobemo_like_virus_F1517C31F1352
MT361059_Guangzhou_sobemo-like_virus
MT361055_Guangzhou_sobemo-like_virus
MT096519_Wenzhou_sobemo-like_virus_4
MZ556263_Sichuan_mosquito_sobemo-like_virus
NC_033138_Wenzhou_sobemo-like_virus_4
MT361053_Guangzhou_sobemo-like_virus
Guangzhou_sobemo_like_virus_F1517C31F1017
MW434812_Guadeloupe_mosquito_virus
MK285337_Renna_virus
MW434813_Guadeloupe_mosquito_virus
MT361065_Kwale_mosquito_virus
MN053797_Guadeloupe_mosquito_virus
MN053803_Guadeloupe_mosquito_virus
MN053791_Guadeloupe_mosquito_virus
MN053789_Guadeloupe_mosquito_virus
MN053787_Guadeloupe_mosquito_virus
Comfornsto_virus_F1216C3F1183
Comfornsto_virus_F1216C3F1216
MW434821_Kellev_virus
Colokingly_virus_F1527C217F1100
Colokingly_virus_F1527C217F1094
Colokingly_virus_F1527C217F1534
Colokingly_virus_F1527C217F1021
Colokingly_virus_F1527C217F1382
Colokingly_virus_F1527C217F1211
Colokingly_virus_F1527C217F1381
Colokingly_virus_F1527C217F1527
Colokingly_virus_F1527C217F1096
Colokingly_virus_F1527C217F1019
MW310331_Zeugodacus_cucurbitae_sobemovirus_isolate_Bc
MW310351_Zeugodacus_cucurbitae_sobemovirus
MW199227_Virus_sp.
"
clade_v231 <- str_split(clade_v231_str, "\\s")[[1]]
clade_v231 <- clade_v231[which(!clade_v231 == "")]
keep_for_v231 <- clade_v231
v2_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v2.fasta",
                                seqtype = "AA")
v231_aa_seq <- v2_aa_seq[keep_for_v231]
names(v231_aa_seq)
seqinr::write.fasta(sequences = v231_aa_seq , names =names(v231_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v231.fasta")
##keep the closest sequence as an outgroup
new_outg <- "MZ443571_Vespula_vulgaris_Sobemo-like_virus_1"
new_outg_seq <- v2_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))



####### Zooming in on different clades
#### v2311 tree
##remove some of the mapping duplicates because viruses are too close
novelseq_to_remove <- c("Pirounced_virus_F1032C17F1127", 
                        "Peninsere_virus_F1127C300F1461", 
                        "Peninsere_virus_F1127C300F1032",
                        "Peninsere_virus_F1127C300F1056")

clade_v2311_str <- "Peninsere_virus_F1127C300F1056
Pirounced_virus_F1032C17F1305
Pirounced_virus_F1032C17F1585
Pirounced_virus_F1032C17F1271
Pirounced_virus_F1032C17F1316
Pirounced_virus_F1032C17F1032
Pirounced_virus_F1032C17F1266
Pirounced_virus_F1032C17F1034
Pirounced_virus_F1032C17F1056
Peninsere_virus_F1127C300F1032
Pirounced_virus_F1032C17F1607
Peninsere_virus_F1127C300F1461
Pirounced_virus_F1032C17F1604
Pirounced_virus_F1032C17F1461
Pirounced_virus_F1032C17F1325
LC513833_Culex_inatomii_luteo-like_virus
Peninsere_virus_F1127C300F1127
Peninsere_virus_F1127C300F1156
Pirounced_virus_F1032C17F1127
Peninsere_virus_F1127C300F1132
Peninsere_virus_F1127C300F1450
MW030476_Hubei_mosquito_virus_2
KX882764_Hubei_mosquito_virus_2
MW030482_Hubei_mosquito_virus_2
MW030484_Hubei_mosquito_virus_2
MW030472_Hubei_mosquito_virus_2
MW030480_Hubei_mosquito_virus_2
MW030488_Hubei_mosquito_virus_2
MW030478_Hubei_mosquito_virus_2
MW030486_Hubei_mosquito_virus_2
MW030470_Hubei_mosquito_virus_2
MW434902_Marma_virus
MW434875_Marma_virus
MW434863_Marma_virus
MW434872_Marma_virus
MK440646_Marma_virus
NC_033305_Hubei_mosquito_virus_2
MT361063_Kisumu_mosquito_virus
MH703049_Yongsan_sobemo-like_virus_1
MN661099_Atrato_Sobemo-like_virus_5
MN661097_Atrato_Sobemo-like_virus_5
MN661107_Atrato_Sobemo-like_virus_5
MT896210_Aedes_sobemo-like_virus
MT896211_Aedes_sobemo-like_virus
Guangzhou_sobemo_like_virus_F1517C31F1016
Guangzhou_sobemo_like_virus_F1517C31F1083
Guangzhou_sobemo_like_virus_F1517C31F1078
Guangzhou_sobemo_like_virus_F1517C31F1517
Guangzhou_sobemo_like_virus_F1517C31F1351
Guangzhou_sobemo_like_virus_F1517C31F1339
Guangzhou_sobemo_like_virus_F1517C31F1071
MT361061_Guangzhou_sobemo-like_virus
Guangzhou_sobemo_like_virus_F1517C31F1070
Guangzhou_sobemo_like_virus_F1517C31F1073
Guangzhou_sobemo_like_virus_F1517C31F1352
Guangzhou_sobemo_like_virus_F1517C31F1171
MT361059_Guangzhou_sobemo-like_virus
MT361055_Guangzhou_sobemo-like_virus
MT361053_Guangzhou_sobemo-like_virus
Guangzhou_sobemo_like_virus_F1517C31F1017
MT096519_Wenzhou_sobemo-like_virus_4
MZ556263_Sichuan_mosquito_sobemo-like_virus
NC_033138_Wenzhou_sobemo-like_virus_4
MW434812_Guadeloupe_mosquito_virus
MK285337_Renna_virus
MW434813_Guadeloupe_mosquito_virus
MT361065_Kwale_mosquito_virus
MN053797_Guadeloupe_mosquito_virus
MN053803_Guadeloupe_mosquito_virus
MN053791_Guadeloupe_mosquito_virus
MN053789_Guadeloupe_mosquito_virus
MN053787_Guadeloupe_mosquito_virus
Comfornsto_virus_F1216C3F1183
Comfornsto_virus_F1216C3F1216
Colokingly_virus_F1527C217F1100
Colokingly_virus_F1527C217F1382
Colokingly_virus_F1527C217F1534
Colokingly_virus_F1527C217F1021
Colokingly_virus_F1527C217F1094
Colokingly_virus_F1527C217F1211
Colokingly_virus_F1527C217F1381
Colokingly_virus_F1527C217F1527
Colokingly_virus_F1527C217F1096
Colokingly_virus_F1527C217F1019
"
clade_v2311 <- str_split(clade_v2311_str, "\\s")[[1]]
clade_v2311 <- clade_v2311[which(!clade_v2311 == "")]
clade_v2311 <- clade_v2311[which(!clade_v2311 %in% novelseq_to_remove)]
keep_for_v2311 <- clade_v2311
v2_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v2.fasta",
                                seqtype = "AA")
v2311_aa_seq <- v2_aa_seq[keep_for_v2311]
names(v2311_aa_seq)
seqinr::write.fasta(sequences = v2311_aa_seq , names =names(v2311_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v2311_dedup.fasta")
##keep the closest sequence as an outgroup
new_outg <- "MW434821_Kellev_virus"
new_outg_seq <- v2_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))


####### Zooming in on different clades
#### v23111 tree
clade_v23111_str <- "Peninsere_virus_F1127C300F1056
Pirounced_virus_F1032C17F1305
Pirounced_virus_F1032C17F1585
Pirounced_virus_F1032C17F1271
Pirounced_virus_F1032C17F1316
Pirounced_virus_F1032C17F1032
Pirounced_virus_F1032C17F1266
Pirounced_virus_F1032C17F1034
Pirounced_virus_F1032C17F1056
Peninsere_virus_F1127C300F1032
Pirounced_virus_F1032C17F1607
Peninsere_virus_F1127C300F1461
Pirounced_virus_F1032C17F1604
Pirounced_virus_F1032C17F1461
Pirounced_virus_F1032C17F1325
LC513833_Culex_inatomii_luteo-like_virus
Peninsere_virus_F1127C300F1127
Peninsere_virus_F1127C300F1156
Pirounced_virus_F1032C17F1127
Peninsere_virus_F1127C300F1132
Peninsere_virus_F1127C300F1450
MW030476_Hubei_mosquito_virus_2
KX882764_Hubei_mosquito_virus_2
MW030482_Hubei_mosquito_virus_2
MW030484_Hubei_mosquito_virus_2
MW030472_Hubei_mosquito_virus_2
MW030480_Hubei_mosquito_virus_2
MW030488_Hubei_mosquito_virus_2
MW030478_Hubei_mosquito_virus_2
MW030486_Hubei_mosquito_virus_2
MW030470_Hubei_mosquito_virus_2
MW434902_Marma_virus
MW434875_Marma_virus
MW434863_Marma_virus
MW434872_Marma_virus
MK440646_Marma_virus
NC_033305_Hubei_mosquito_virus_2
"
clade_v23111 <- str_split(clade_v23111_str, "\\s")[[1]]
clade_v23111 <- clade_v23111[which(!clade_v23111 == "")]
keep_for_v23111 <- clade_v23111
v2_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v2.fasta",
                                seqtype = "AA")
v23111_aa_seq <- v2_aa_seq[keep_for_v23111]
names(v23111_aa_seq)
seqinr::write.fasta(sequences = v23111_aa_seq , names =names(v23111_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v23111.fasta")
##keep the closest sequence as an outgroup
new_outg <- "MT361063_Kisumu_mosquito_virus"
new_outg_seq <- v2_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))























####### Zooming in on different clades
#### v24 tree
clade_v24_str <- "MT568530_Pyongtaek_Culex_sobemo-like_virus
NC_033128_Wenzhou_sobemo-like_virus_3
OL700081_Wenzhou_sobemo-like_virus_3
OL700082_Wenzhou_sobemo-like_virus_3
MT568542_Pyongtaek_Culex_sobemo-like_virus
LC512855_Wenzhou_sobemo-like_virus_3
LC512854_Wenzhou_sobemo-like_virus_3
LC512856_Wenzhou_sobemo-like_virus_3
LC512857_Wenzhou_sobemo-like_virus_3
MK440647_Culex_associated_luteo_like_virus
MW434133_Culex-associated_Luteo-like_virus
MW434131_Culex-associated_Luteo-like_virus
MW434130_Culex-associated_Luteo-like_virus
MH188016_Culex-associated_Luteo-like_virus
MW434132_Culex-associated_Luteo-like_virus
Culex_associated_Luteo_like_virus_F1042C4F1042
MW317137_Luteoviridae_sp.
MW434129_Culex-associated_Luteo-like_virus
"
clade_v24 <- str_split(clade_v24_str, "\\s")[[1]]
clade_v24 <- clade_v24[which(!clade_v24 == "")]
keep_for_v24 <- clade_v24
v2_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v2.fasta",
                                seqtype = "AA")
v24_aa_seq <- v2_aa_seq[keep_for_v24]
names(v24_aa_seq)
seqinr::write.fasta(sequences = v24_aa_seq , names =names(v24_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/Sobeli_clade_aa_v24.fasta")
##keep the closest sequence as an outgroup
new_outg <- "MZ443576_Vespula_vulgaris_Luteo-like_virus_2"
new_outg_seq <- v2_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))



####### Zooming in on different clades
#### v241 tree
clade_v241_str <- "MW434130_Culex-associated_Luteo-like_virus
MW434132_Culex-associated_Luteo-like_virus
MH188016_Culex-associated_Luteo-like_virus
MW434133_Culex-associated_Luteo-like_virus
MW434131_Culex-associated_Luteo-like_virus
MW434129_Culex-associated_Luteo-like_virus
MW317137_Luteoviridae_sp.
Culex_associated_Luteo_like_virus_F1042C4F1042
MK440647_Culex_associated_luteo_like_virus
"
clade_v241 <- str_split(clade_v241_str, "\\s")[[1]]
clade_v241 <- clade_v241[which(!clade_v241 == "")]
keep_for_v241 <- clade_v241
v2_nt_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/nt_seq/Sobeli_clade_nt.fasta",
                                seqtype = "AA")
v241_nt_seq <- v2_nt_seq[keep_for_v241]
names(v241_nt_seq)
seqinr::write.fasta(sequences = v241_nt_seq , names =names(v241_nt_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/nt_seq/Sobeli_clade_nt_v241.fasta")
##keep the closest sequence as an outgroup
new_outg <- "LC512857_Wenzhou_sobemo-like_virus_3"
new_outg_seq <- v2_nt_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/aln/nt_seq/", 
                                      new_outg, ".fasta"))




