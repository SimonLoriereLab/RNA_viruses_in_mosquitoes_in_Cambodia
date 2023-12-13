library(tidyverse)
library(seqinr)
library(treeio)

Phasmaviridae_clade_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/Phasmaviridae_clade_novel_sequence_data.txt")

### select consensus sequences based on mapping coverage
### first round of read mapping results
read_map_df <-  read_delim("/full_path_to/wd/RdRp_scan/analysis/read_mapping/mapping_summary_bowtie2_RdRpScan_Round1_analysis/results_mapped_only.txt", 
           delim = "\t")


read_map_df_Phasmaviridae_clade <- read_map_df[which(read_map_df$ref_id %in% Phasmaviridae_clade_df$contig_id & read_map_df$selected == "yes"),]



read_map_df_Phasmaviridae_clade$RdRp_scan_map_consensus_id <- paste0(read_map_df_Phasmaviridae_clade$ref_id, "_",
                      read_map_df_Phasmaviridae_clade$sample_id, "_RdRpScan_Round1")

read_map_df_Phasmaviridae_clade$contig_id <- read_map_df_Phasmaviridae_clade$ref_id

consensus_df <- left_join(read_map_df_Phasmaviridae_clade, Phasmaviridae_clade_df, by  = "contig_id")



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
extracted_consensus_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/extracted_cons1_renamed/"
dir.create(extracted_consensus_directory)
extracted_orf_nt_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/extracted_cons1_orf_nt/"
dir.create(extracted_orf_nt_directory)
extracted_orf_aa_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/extracted_cons1_orf_aa/"
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

write.table(x = consensus_df, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/consensus_df.txt", 
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




write.table(x = novel_seq_mdf1, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/updated_consensus_df.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)


# ##### Manually correct frames for some of the sequences. Reimport and trim parts beyond stop codons
orf_aa_w_stops <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/Phasmaviridae_orf_aa_manually_corrected.fasta",
                                     seqtype = "AA")
extracted_orf_aa_directory1 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/extracted_cons1_orf_aa_corrected/"
dir.create(extracted_orf_aa_directory1)

for (n in names(orf_aa_w_stops)){
  s <- str_split(paste0(as.character(orf_aa_w_stops[n][[1]]), collapse = ""), "\\*")
  t <- tibble(seq = s[[1]], len  = str_length(s[[1]]))
  s1 <- as.character(t[which(t$len == max(t$len)),][1,'seq'])
  s2 <- as.SeqFastaAA(s1, name = n, Annot = n)
  seqinr::write.fasta(sequences = s2, names = n,
                      file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/extracted_cons1_orf_aa_corrected/",
                                        n, ".fasta"))
}





#### REFERENCE SEQUENCES

### Use DIAMOND blastp results to get the hits
blastp_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/lineage_collapsed.DIAMOND_blastp_hits_nrcontigs_rdrp.txt.gz", 
                        delim = "\t") %>% 
  select(orf_id = qseqid, sseqid, pident, length, stitle, tax_id, superkingdom)
blastp_df1 <- blastp_df[which(blastp_df$superkingdom == "Viruses"),]

blastp_df2 <- blastp_df1[which(blastp_df1$orf_id %in% sort(unique(consensus_df$orf_id))),]


blastp_df2$sseqid


### Add potentially more distant DIAMOND blastp results from the web blastp
web_blastp_colnames <- c("qseqid", "sseqid", "pident", "length", "mismatches", 
                         "gapopens", "qstart", "qend", "sstart", "send", 
                         "evalue", "bitscore", "perc_positives")


web_blastp_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/VBE3B3HY013-Alignment.txt", 
                        delim = "\t", col_names = web_blastp_colnames, comment = "#") %>% 
  select(qseqid, sseqid, pident, length)

sort(unique(web_blastp_df$sseqid))

#### Let's also add RdRp_scan sequences from the same clade
rdrp_scan_clade_accns <- str_remove(str_remove(str_split(consensus_df$ref_suggestions, ",")[[1]], "[[:print:]]+viricetes_"), "\\.[[:digit:]]{1}\\_[[:print:]]+$")
rdrp_scan_clade_accns <- str_remove(rdrp_scan_clade_accns, "[[:print:]]*NOT_ASSIGNED_")
rdrp_scan_clade_accns[which(rdrp_scan_clade_accns == "Henan_sediment_articulavirus")] <- "UMO75703"
rdrp_scan_clade_accns[which(rdrp_scan_clade_accns == "Hainan_phenui_like_virus_3")] <- "QYF49548"
rdrp_scan_clade_accns[which(rdrp_scan_clade_accns == "Sichuan_10_bunya_like_virus_3")] <- "QYF49513"
rdrp_scan_clade_accns[which(rdrp_scan_clade_accns == "Sichuan_bunya_like_virus_2")] <- "QYF49512"

refs <- sort(unique(c(blastp_df2$sseqid, rdrp_scan_clade_accns, web_blastp_df$sseqid)))
refs <- str_remove(refs, "\\.[[:digit:]]$")
refs[which(refs == "Henan_orthophasma_like_virus_1")] <- "QYF49524"
length(sort(unique(refs)))



##### Check on NCBI Viruses
paste0(sort(unique(refs)), collapse = ", ")


#### Reimport metadata
NCBI_Virus_mdf <- read_csv("/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/Phasmaviridae_metadata_public.csv")


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
selected_mdf_upd_novel$nt_seq_ID <- str_replace_all(selected_mdf_upd_novel$nt_seq_ID, "[\\.\\/]", "_")


# selected_mdf_upd_novel[which(str_detect(selected_mdf_upd_novel$nt_seq_ID,"OP589941")),]
# selected_mdf[which(str_detect(selected_mdf$nt_seq_ID,"OP589941")),]

temp_check <- selected_mdf_upd_novel[which(is.na(selected_mdf_upd_novel$host_genus)),] ## check 


write.table(x = selected_mdf_upd_novel, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/updated_metadata_with_novel_seq.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)



######################## Dealing with CDS
####### Extract CDS from NCBI nucleotide, not NCBI Virus, because the headers are more informative
paste0(sort(selected_mdf$nt_accession), collapse = ",")

####### Here, because there could be multiple CDS, use protein accessions to extract CDS
####### Correspondence between the prot and nt accessions is encoded in standard CDS names when extraction from NCBI Protein same as from NCBI Nucleotide

paste0(unique(sort(refs)), collapse = ",")

length(unique(sort(refs)))


### using this accession list to NCBI and extract CDS sequences and import them here.
#### Put the file here!
cds_fasta_file1 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/Phasmaviridae_clade_refs_cds_nt.txt"
# cds_fasta_file2 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/Phasmaviridae_clade_refs_cds_nt2.txt"

### There might be more automatic solutions for CDS extractions with EMBOSS or CL implementation of entrez. But web NCBI is the fastest for now.

############################################################################################
############################################################################################
############################################################################################
cds_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/ref_cds_nt"
dir.create(cds_clade_dir)
aa_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/ref_cds_aa"
dir.create(aa_clade_dir)


cds_fasta1 <- seqinr::read.fasta(file = cds_fasta_file1, seqtype = "DNA", whole.header = T)
# cds_fasta2 <- seqinr::read.fasta(file = cds_fasta_file2, seqtype = "DNA", whole.header = T)

names(cds_fasta1)
length(selected_mdf$nt_accession)
length(names(cds_fasta1))
# length(names(cds_fasta2))

# length(names(cds_fasta1)) + length(names(cds_fasta2))

cds_fasta <- cds_fasta1
# cds_fasta <- c(cds_fasta1, cds_fasta2)



cds_names_df <- tibble(original_cds_headers = names(cds_fasta),
                       nt_accession = str_remove(str_extract(names(cds_fasta), "(?<=^lcl\\|)[[:print:]]+(?=_cds)"), "\\.[[:digit:]]*"),
                       cds_product = str_replace_all(str_extract(names(cds_fasta), "(?<=\\[protein\\=)[[:print:]]+?(?=\\])"), " ", "_"),
                       prot_accession = str_remove(str_extract(names(cds_fasta), "(?<=_cds_)[[:alnum:]\\_\\.]+(?= )"), "\\.[[:digit:]]*\\_[[:digit:]]*")
)





ref_df <- selected_mdf
# paste0(setdiff(sort(cds_names_df$nt_accession), sort(ref_df$nt_accession)), collapse = ",") ### check





######### Checking if any CDS are missing, troubleshooting
# paste0(setdiff(unique(refs), sort(cds_names_df$prot_accession)),  collapse = ",") ## retry to extract from NCBI prot again
# prot_accns <- setdiff(unique(refs), sort(cds_names_df$prot_accession))
# intersect(cds_names_df$prot_accession, prot_accns)
# intersect(refs, prot_accns)
# intersect(str_remove(sort(unique(web_blastp_df$sseqid)),"\\.[[:digit:]]"), prot_accns)
# intersect(str_remove(blastp_df2$sseqid, "\\.[[:digit:]]"), prot_accns)
# intersect(str_remove(rdrp_scan_clade_accns, "\\.[[:digit:]]"), prot_accns)






















### Reimport genbank files and get metadata from them. They are not in NCBI Virus perhaps because they are MAG: ?
# gb_parse_metadata_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/Phasmaviridae_to_parse_metadata.gb"
# count = 0
# entry_ends = FALSE
# for(i in 1:length(read_lines(gb_parse_metadata_file))){
#   entry_ends = FALSE
#   l <-  read_lines(gb_parse_metadata_file)[i]
#   # print(l)
#   if (str_detect(l, "LOCUS") == TRUE){
#     entry_ends = FALSE
#     nt_accession <- str_remove(str_extract(l, "LOCUS\\s+[[:alnum:]\\_]+"), "LOCUS\\s*")
#     nt_slen <- str_remove(str_extract(l, "(?<=\\s)[[:digit:]]+(?=\\s+bp)"), "\\s*")
#     nt_title <- str_remove(read_lines(gb_parse_metadata_file)[i+1], "DEFINITION\\s*")
#   }
#   
#   if (str_detect(l, "\\/host\\=") == TRUE){
#     host_taxon <- str_remove_all(str_remove(l, "\\s*\\/host\\="), '"')
#   }
#   
#   if (str_detect(l, "\\/organism\\=") == TRUE){
#     Species <- str_remove_all(str_remove(l, "\\s*\\/organism\\="), '"')
#   }
#   
#   if (str_detect(l, "\\/country\\=") == TRUE){
#     Location <- str_remove_all(str_remove(l, "\\s*\\/country\\="), '"')
#   }
#   
#   if (str_detect(l, "\\/collection_date\\=") == TRUE){
#     Year <- str_remove_all(str_remove(l, "\\s*\\/collection_date\\="), '"')
#   }
# 
#   if (str_detect(l, "^\\/\\/") == TRUE){
#     entry_ends = TRUE
#   }
#   
#   if (entry_ends == TRUE){
#     ## save row and add count
#     count = count + 1
#     metadata_row <- tibble(nt_accession, nt_slen, nt_title, host_taxon, Location, Year, Species)
#     
#     if (count == 1){
#       metadata <- metadata_row
#     }else{
#       metadata <- bind_rows(metadata, metadata_row)
#     }
#   }
# }
# 
# metadata$nt_seq_ID <- paste0(metadata$nt_accession, "_", str_replace_all(metadata$Species, " ", "_"))
# metadata$nt_slen <- as.numeric(as.character(metadata$nt_slen))
# ref_df <- bind_rows(ref_df, metadata)



unique(na.omit(str_extract(ref_df$nt_seq_ID, "[\\.\\+\\*\\?\\^\\$\\(\\)\\[\\]\\{\\}\\|\\\\\\/]")))
ref_df$nt_seq_ID <- str_replace_all(ref_df$nt_seq_ID, "[\\.\\/]", "_")
selected_mdf_upd_novel1 <- selected_mdf_upd_novel
# selected_mdf_upd_novel1 <- bind_rows(selected_mdf_upd_novel, metadata)
selected_mdf_upd_novel1$nt_seq_ID <- str_replace_all(selected_mdf_upd_novel1$nt_seq_ID, "[\\.\\/]", "_")


write.table(x = selected_mdf_upd_novel1, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/updated_metadata_with_novel_seq1.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)





paste0(setdiff(sort(ref_df$nt_accession), sort(cds_names_df$nt_accession)), collapse = ",")

paste0(setdiff(sort(refs), sort(cds_names_df$prot_accession)), collapse = ",")


cds_names_df1 <- inner_join(cds_names_df, ref_df, by = "nt_accession")

cds_names_df1$new_short_name_cds <- cds_names_df1$nt_seq_ID
cds_names_df1$cds_status <- "OK"


sort(unique(cds_names_df1$cds_product))



cds_fasta <- cds_fasta[cds_names_df1$original_cds_headers]


genome_fraction_threshold_for_cds <- 0 ### no threshold, because it's not a polyprotein, cds can be quite small

for (cds_name in names(cds_fasta)){
  cds <- cds_fasta[cds_name]
  print(cds_name)
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
  
  
  
  #### special case of alternative genetic code
  if (str_detect(cds_name, "NC_033284") == TRUE){
    cds_trans <- translate(seq = cds_seq,
                           frame = 0, sens = "F", NAstring = "X", ambiguous = T, numcode = 6)
  }else{
    cds_trans <- translate(seq = cds_seq,
                           frame = 0, sens = "F", NAstring = "X", ambiguous = T)
  }
  

  
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


####### Zooming in on different clades
#### Phasmaviridae_clade1 tree
Phasmaviridae_clade1_str <- "KP710260_Ferak_feravirus
KP710250_Ferak_feravirus
KP710256_Ferak_feravirus
KP710259_Ferak_feravirus
KP710258_Ferak_feravirus
KP710257_Ferak_feravirus
KP710263_Ferak_feravirus
KP710254_Ferak_feravirus
KP710262_Ferak_feravirus
KP710247_Ferak_feravirus
KP710255_Ferak_feravirus
NC_043031_Ferak_feravirus
KP710253_Ferak_feravirus
KP710252_Ferak_feravirus
KP710248_Ferak_feravirus
KP710261_Ferak_feravirus
KP710249_Ferak_feravirus
KP710251_Ferak_feravirus
MN661015_Guagua_virus
MW039262_Neuropteran_feravirus
MT153364_Hemipteran_feravirus
NC_055192_Sanxia_sawastrivirus
MF893255_Notori_virus
MT153462_Hymenopteran_phasma-related_virus_OKIAV244
KX884782_Hubei_bunya-like_virus_10
NC_030752_Insect_wuhivirus
OL752429_Aphid_bunyavirus_1
MZ209865_Sanya_sesamia_inferens_phasmavirus_1
MZ209951_Sanya_sesamia_inferens_phasmavirus_1
NC_040462_Saesbyeol_virus
MZ209660_Hangzhou_frankliniella_intonsa_phasmavirus_1
MT153437_Lepidopteran_phasma-related_virus_OKIAV246
Cavaretrus_virus_F1611C2F1608
Cavaretrus_virus_F1611C2F1611
Cavaretrus_virus_F1611C2F1049
KT966494_Terena_virus
MW288188_Hemipteran_phasma-related_virus_OKIAV245
MZ209649_Hangzhou_phasmavirus_1
KP710239_Jonchet_jonvirus
NC_038706_Jonchet_jonvirus
KP710236_Jonchet_jonvirus
KP710237_Jonchet_jonvirus
KP710235_Jonchet_jonvirus
KP710234_Jonchet_jonvirus
KP710233_Jonchet_jonvirus
KP710238_Jonchet_jonvirus
MZ202272_Mikado_virus
MZ202275_Mikado_virus
MZ202269_Spilikins_virus
Cloizable_virus_F1036C2F1036
"
Phasmaviridae_clade1 <- str_split(Phasmaviridae_clade1_str, "\\s")[[1]]
Phasmaviridae_clade1 <- Phasmaviridae_clade1[which(!Phasmaviridae_clade1 == "")]

Phasmaviridae_clade_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/aln/aa_seq/Phasmaviridae_clade_aa.fasta",
                                seqtype = "AA")
Phasmaviridae_clade1_aa_seq <- Phasmaviridae_clade_aa_seq[Phasmaviridae_clade1]
names(Phasmaviridae_clade1_aa_seq)
seqinr::write.fasta(sequences = Phasmaviridae_clade1_aa_seq , names =names(Phasmaviridae_clade1_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/aln/aa_seq/Phasmaviridae_clade1_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "MT153380_Collembolan_phasma-related_virus_OKIAV223"
new_outg_seq <- Phasmaviridae_clade_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))


#### Phasmaviridae_clade2 tree
Phasmaviridae_clade2_str <- "MT361031_Aedes_phasmavirus
MT361034_Aedes_phasmavirus
MT361037_Aedes_phasmavirus
MT361028_Aedes_phasmavirus
MT361040_Aedes_phasmavirus
MT361043_Aedes_phasmavirus
MT361049_Aedes_phasmavirus
Aedes_phasmavirus_F1083C4F1083
MW434650_Orthophasmavirus_bastukasense
MW434651_Orthophasmavirus_barstukasense
MW434652_Orthophasmavirus_barstukasense
MW434654_Orthophasmavirus_barstukasense
MW434660_Orthophasmavirus_barstukasense
Aedes_phasmavirus_F1083C4F1520
ON955170_Kuusamo_orthophasmavirus_4
MN513376_Flen_bunya-like_virus
NC_040759_Yongsan_bunyavirus_1
MW434693_Orthophasmavirus_miglotasense
MW434694_Orthophasmavirus_moglotasense
MW434698_Orthophasmavirus_miglotasense
MW434699_Orthophasmavirus_moglotasense
MW434704_Orthophasmavirus_moglotasense
MW434703_Orthophasmavirus_moglotasense
MW434700_Orthophasmavirus_miglotasense
MW434692_Orthophasmavirus_moglotasense
MW434697_Orthophasmavirus_miglotasense
MW434696_Orthophasmavirus_miglotasense
MW434695_Orthophasmavirus_moglotasense
MW434705_Orthophasmavirus_miglotasense
MN661018_Culex_orthophasmavirus
MF176293_Culex_orthophasmavirus
MF176313_Culex_orthophasmavirus
NC_055215_Culex_orthophasmavirus
MK440628_Anjon_virus
OL700053_Culex_orthophasmavirus
MW452284_Wuhan_mosquito_orthophasmavirus_2
NC_031312_Wuhan_mosquito_orthophasmavirus_2
OL700113_XiangYun_bunya-arena-like_virus_12
ON955163_Hameenlinna_orthophasmavirus_1
ON955161_Hameenlinna_orthophasmavirus_1
ON955162_Hameenlinna_orthophasmavirus_1
ON955171_Lestijarvi_orthophasmavirus_1
ON955167_Kuusamo_orthophasmavirus_2
MN661021_Orthophasmavirus_coredoense
ON059782_Niukluk_phantom_orthophasmavirus
ON059784_Niukluk_phantom_orthophasmavirus
MN168168_Niukluk_phantom_orthophasmavirus
ON059783_Niukluk_phantom_orthophasmavirus
NC_034462_Kigluaik_phantom_orthophasmavirus
KX884805_Hubei_diptera_virus_6
OK491514_Hubei_diptera_virus_6
MT153491_Dipteran_phasma-related_virus_OKIAV224
MT153492_Dipteran_phasma-related_virus_OKIAV225
MW310369_Bactrocera_latifrons_orthophasmavirus
MW452283_Wuhan_mosquito_orthophasmavirus_1
OL700085_Wuhan_mosquito_orthophasmavirus_1
NC_031307_Wuhan_mosquito_orthophasmavirus_1
NC_055394_Anopheles_orthophasmavirus
MT129676_Peralta_bunya-like_virus
MT129677_Peralta_bunya-like_virus
KY354237_Apis_bunyavirus_2
NC_038745_Nome_phantom_orthophasmavirus
KX884807_Hubei_diptera_virus_7
"
Phasmaviridae_clade2 <- str_split(Phasmaviridae_clade2_str, "\\s")[[1]]
Phasmaviridae_clade2 <- Phasmaviridae_clade2[which(!Phasmaviridae_clade2 == "")]

Phasmaviridae_clade_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/aln/aa_seq/Phasmaviridae_clade_aa.fasta",
                                                    seqtype = "AA")
Phasmaviridae_clade2_aa_seq <- Phasmaviridae_clade_aa_seq[Phasmaviridae_clade2]
names(Phasmaviridae_clade2_aa_seq)
seqinr::write.fasta(sequences = Phasmaviridae_clade2_aa_seq , names =names(Phasmaviridae_clade2_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/aln/aa_seq/Phasmaviridae_clade2_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "MN164619_Pink_bollworm_virus_2"
new_outg_seq <- Phasmaviridae_clade_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Phasmaviridae_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))

