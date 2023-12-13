library(tidyverse)
library(seqinr)
library(treeio)

Nodamuvirales_clade_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/Nodamuvirales_clade_novel_sequence_data.txt")

### select consensus sequences based on mapping coverage
### first round of read mapping results
read_map_df <-  read_delim("/full_path_to/wd/RdRp_scan/analysis/read_mapping/mapping_summary_bowtie2_RdRpScan_Round1_analysis/results_mapped_only.txt", 
           delim = "\t")


read_map_df_Nodamuvirales_clade <- read_map_df[which(read_map_df$ref_id %in% Nodamuvirales_clade_df$contig_id & read_map_df$selected == "yes"),]



read_map_df_Nodamuvirales_clade$RdRp_scan_map_consensus_id <- paste0(read_map_df_Nodamuvirales_clade$ref_id, "_",
                      read_map_df_Nodamuvirales_clade$sample_id, "_RdRpScan_Round1")

read_map_df_Nodamuvirales_clade$contig_id <- read_map_df_Nodamuvirales_clade$ref_id

consensus_df <- left_join(read_map_df_Nodamuvirales_clade, Nodamuvirales_clade_df, by  = "contig_id")



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
extracted_consensus_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/extracted_cons1_renamed/"
dir.create(extracted_consensus_directory)
extracted_orf_nt_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/extracted_cons1_orf_nt/"
dir.create(extracted_orf_nt_directory)
extracted_orf_aa_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/extracted_cons1_orf_aa/"
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

write.table(x = consensus_df, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/consensus_df.txt", 
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




write.table(x = novel_seq_mdf1, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/updated_consensus_df.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)




##### Manually correct frames for some of the sequences. Reimport and trim parts beyond stop codons
# orf_aa_w_stops <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/Nodamuvirales_orf_aa_manually_corrected.fasta",
#                                      seqtype = "AA")
# extracted_orf_aa_directory1 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/extracted_cons1_orf_aa_corrected/"
# dir.create(extracted_orf_aa_directory1)
# 
# for (n in names(orf_aa_w_stops)){
#   s <- str_split(paste0(as.character(orf_aa_w_stops[n][[1]]), collapse = ""), "\\*")
#   t <- tibble(seq = s[[1]], len  = str_length(s[[1]]))
#   s1 <- as.character(t[which(t$len == max(t$len)),][1,'seq'])
#   s2 <- as.SeqFastaAA(s1, name = n, Annot = n)
#   seqinr::write.fasta(sequences = s2, names = n,
#                       file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/extracted_cons1_orf_aa_corrected/",
#                                         n, ".fasta"))
# }




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


web_blastp_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/ZYUSZ2HF016-Alignment.txt", 
                        delim = "\t", col_names = web_blastp_colnames, comment = "#") %>% 
  select(qseqid, sseqid, pident, length)

sort(unique(web_blastp_df$sseqid))

#### Let's also add RdRp_scan sequences from the same clade
rdrp_scan_clade_accns <- str_remove(str_remove(str_split(consensus_df$ref_suggestions, ",")[[1]], "[[:print:]]+Nodamuvirales_"), "\\.[[:digit:]]{1}\\_[[:print:]]+$")
rdrp_scan_clade_accns <- str_remove(rdrp_scan_clade_accns, "[[:print:]]*NOT_ASSIGNED_")

# rdrp_scan_clade_accns <- str_remove(rdrp_scan_clade_accns, "^[[:alnum:]]*_")
rdrp_scan_clade_accns <- rdrp_scan_clade_accns[which(!str_detect(rdrp_scan_clade_accns, "\\|"))] ### remove ones for which accessions not found
rdrp_scan_clade_accns[which(rdrp_scan_clade_accns == "")] <- ""

rdrp_scan_clade_accns <- rdrp_scan_clade_accns[which(!rdrp_scan_clade_accns %in% c())]

refs <- sort(unique(c(blastp_df2$sseqid, rdrp_scan_clade_accns, web_blastp_df$sseqid)))
refs <- str_remove(refs, "\\.[[:digit:]]$")
# refs[which(refs == "")] <- ""
length(sort(unique(refs)))



##### Check on NCBI Viruses
paste0(sort(unique(refs)), collapse = ", ")


#### Reimport metadata
NCBI_Virus_mdf <- read_csv("/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/Nodamuvirales_metadata_public.csv")


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


write.table(x = selected_mdf_upd_novel, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/updated_metadata_with_novel_seq.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)



######################## Dealing with CDS
####### Extract CDS from NCBI nucleotide, not NCBI Virus, because the headers are more informative
paste0(sort(selected_mdf$nt_accession), collapse = ",")

####### Here, because there could be multiple CDS, use protein accessions to extract CDS
####### Correspondence between the prot and nt accessions is encoded in standard CDS names when extraction from NCBI Protein same as from NCBI Nucleotide


### Only extract prot refs that are viruses based on presence in NCBI Virus, because a lot of accessions fail to be processed in batch extraction in MCBI protein
# refs <- read_lines("/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/Nodamuvirales_accessions_viruses.txt")

paste0(unique(sort(refs)), collapse = ",")

length(unique(sort(refs)))






### using this accession list to NCBI and extract CDS sequences and import them here.
#### Put the file here!
cds_fasta_file1 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/Nodamuvirales_clade_refs_cds_nt1.txt"
cds_fasta_file2 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/Nodamuvirales_clade_refs_cds_nt2.txt"
cds_fasta_file3 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/Nodamuvirales_clade_refs_cds_nt3.txt"
cds_fasta_file4 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/Nodamuvirales_clade_refs_cds_nt4.txt"



### There might be more automatic solutions for CDS extractions with EMBOSS or CL implementation of entrez. But web NCBI is the fastest for now.

############################################################################################
############################################################################################
############################################################################################
cds_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/ref_cds_nt"
dir.create(cds_clade_dir)
aa_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/ref_cds_aa"
dir.create(aa_clade_dir)


cds_fasta1 <- seqinr::read.fasta(file = cds_fasta_file1, seqtype = "DNA", whole.header = T)
cds_fasta2 <- seqinr::read.fasta(file = cds_fasta_file2, seqtype = "DNA", whole.header = T)
cds_fasta3 <- seqinr::read.fasta(file = cds_fasta_file3, seqtype = "DNA", whole.header = T)
cds_fasta4 <- seqinr::read.fasta(file = cds_fasta_file4, seqtype = "DNA", whole.header = T)

names(cds_fasta1)
length(selected_mdf$nt_accession)
length(names(cds_fasta1))
# length(names(cds_fasta2))

# length(names(cds_fasta1)) + length(names(cds_fasta2))

# cds_fasta <- cds_fasta1
cds_fasta <- c(cds_fasta1, cds_fasta2, cds_fasta3, cds_fasta4)


# 
# 
# cds_fasta<- seqinr::read.fasta(file = cds_fasta_file, seqtype = "DNA", whole.header = T)
# 

names(cds_fasta)
length(selected_mdf$nt_accession)
length(names(cds_fasta))


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
### ignoring these 






### Reimport genbank files and get metadata from them. They are not in NCBI Virus perhaps because they are MAG: ?
# gb_parse_metadata_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/Nodamuvirales_to_parse_metadata.gb"
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


write.table(x = selected_mdf_upd_novel1, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/updated_metadata_with_novel_seq1.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)


length(unique(refs))


paste0(setdiff(sort(ref_df$nt_accession), sort(cds_names_df$nt_accession)), collapse = ",")
paste0(setdiff(sort(cds_names_df$prot_accession), sort(refs)), collapse = ",")
paste0(setdiff(sort(refs), sort(cds_names_df$prot_accession)), collapse = ",")



### These are viral proteins that don't have link to protein/nuccore because they are from PDB or UniProt.
### Ignoring these.


cds_names_df1 <- inner_join(cds_names_df, ref_df, by = "nt_accession")

cds_names_df1$new_short_name_cds <- cds_names_df1$nt_seq_ID
cds_names_df1$cds_status <- "OK"


sort(unique(cds_names_df1$cds_product))

### remove at least capsid proteins
# cds_names_df1 <- cds_names_df1[which(!cds_names_df1$cds_product %in% c("capsid")),]

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
#### Nodamuvirales_clade1 tree
Nodamuvirales_clade1_str <- "MN862357_Emilia_sonchifolia_noda-like_virus
MW239267_Riboviria_sp_
MT482496_PNG_bee_virus_14
MT317176_Rice_Noda-like_virus
OL472288_Koper_noda-like_virus_2
ON260786_rice-associated_noda-like_virus_1
MW239460_Riboviria_sp_
ON049872_Nodaviridae_sp_
KX883009_Hubei_noda-like_virus_10
ON872547_Bat_faecal_associated_nodavirus_3
MW239448_Riboviria_sp_
MW239466_Riboviria_sp_
MW826495_Nodaviridae_sp_
MW897079_Guangxi_forest_noda-like_virus
MW897128_Henan_sediment_noda-like_virus_4
MW314615_Lasius_neglectus_noda-like_virus_1
OP884074_Army_ant_associated_Nodavirus_3
KJ632942_Mosinovirus
NC_033311_Hubei_orthoptera_virus_4
OM954255_Flumine_nodavirus_27
MW239459_Riboviria_sp_
NC_032732_Shuangao_insect_virus_10
MZ209797_Guiyang_nodavirus_1
MW897119_Zhejiang_farmland_noda-like_virus
MZ209764_Guiyang_nodavirus_1
OM954257_Flumine_nodavirus_3
MW897091_Sichuan_forest_noda-like_virus_3
AY962576_Pieris_rapae_virus
KX883132_Wuhan_house_centipede_virus_6
KX883198_Wuhan_house_centipede_virus_6
NC_033482_Wuhan_house_centipede_virus_6
MH719096_Mayapan_virus
MH794142_Culannivirus_DW-2019a
KX883080_Hubei_noda-like_virus_5
MT753152_Mayapan_virus
OL700048_XiangYun_tombus-noda-like_virus_11
MW897098_Henan_sediment_noda-like_virus_2
NC_003691_Pariacoto_virus
ON049869_Nodaviridae_sp_
MW897103_Sichuan_forest_noda-like_virus_5
KX883167_Hubei_noda-like_virus_6
Murmaxied_virus_F1531C222F1533
Murmaxied_virus_F1531C222F1534
Murmaxied_virus_F1531C222F1537
Murmaxied_virus_F1531C222F1527
Murmaxied_virus_F1531C222F1531
Murmaxied_virus_F1531C222F1532
Murmaxied_virus_F1531C222F1217
Murmaxied_virus_F1531C222F1093
Murmaxied_virus_F1531C222F1094
Murmaxied_virus_F1531C222F1096
NC_033139_Sanxia_water_strider_virus_16
MZ443627_Nelson_wasp-associated_virus_3
OL700186_XiangYun_tombus-noda-like_virus_7
KP263548_Lunovirus
KY214439_Caninovirus_sp_
OL700180_XiangYun_tombus-noda-like_virus_1
ON049855_Nodaviridae_sp_
KX883225_Hubei_noda-like_virus_7
ON049805_Nodaviridae_sp_
MW897073_Xinjiang_sediment_noda-like_virus_2
MW826493_Nodaviridae_sp_
KX883130_Hubei_noda-like_virus_9
NC_033307_Hubei_noda-like_virus_9
MK533155_Tetranychus_urticae-associated_nodavirus_A
ON746452_Tianjin_Nodav_tick_virus_1
OL472289_Koper_noda-like_virus_3
OL700185_XiangYun_tombus-noda-like_virus_6
KX883237_Hubei_noda-like_virus_8
NC_033309_Hubei_noda-like_virus_8
MK533156_Tetranychus_urticae-associated_nodavirus_B
OL472290_Koper_noda-like_virus_4
MW239321_Riboviria_sp_
MZ375165_Riboviria_sp_
OM954246_Flumine_nodavirus_11
OM953982_Flumine_nodavirus_10
MW239203_Riboviria_sp_
OL700188_XiangYun_tombus-noda-like_virus_9
KX883125_Hubei_noda-like_virus_12
MW897072_Beijing_sediment_noda-like_virus_4
MN513380_Hammarskog_noda-like_virus
NC_033077_Sanxia_water_strider_virus_17
KX883010_Hubei_noda-like_virus_11
MT096526_Hubei_noda-like_virus_11
MW897075_Beijing_sediment_noda-like_virus_6
GU976287_Alphanodavirus_HB-2007_CHN
MT138113_Nodaviridae_sp_
MN918672_Picornavirales_sp_
NC_032760_Beihai_sphaeromadae_virus_2
MG874757_Redspotted_grouper_nervous_necrosis_virus
MZ527260_Redspotted_grouper_nervous_necrosis_virus
KP455643_Redspotted_grouper_nervous_necrosis_virus
KU705814_Mouse_grouper_Nervous_Necrosis_Virus
KF668184_Redspotted_grouper_nervous_necrosis_virus
MN309751_Hulong_grouper_nervous_necrosis_virus
HQ859895_Tiger_grouper_Nervous_Necrosis_Virus
AY721616_Redspotted_grouper_nervous_necrosis_virus
FJ803915_Sea_bass_Iberian_betanodavirus
KY354681_Redspotted_grouper_nervous_necrosis_virus
NC_008040_Redspotted_grouper_nervous_necrosis_virus
HQ859908_Tiger_grouper_Nervous_Necrosis_Virus
HQ859911_Asian_seabass_Nervous_Necrosis_Virus
FJ789783_Redspotted_grouper_nervous_necrosis_virus
AB373028_Redspotted_grouper_nervous_necrosis_virus
FJ748760_Redspotted_grouper_nervous_necrosis_virus
JN189865_Striped_jack_nervous_necrosis_virus
OM513982_Betanodavirus_sp_
MN097800_Betanodavirus_sp_
MN097798_Betanodavirus_sp_
KT390713_Grouper_betanodavirus
MZ054261_Redspotted_grouper_nervous_necrosis_virus
HQ859913_Tiger_grouper_Nervous_Necrosis_Virus
GQ402010_Redspotted_grouper_nervous_necrosis_virus
AY369136_Redspotted_grouper_nervous_necrosis_virus
KY315688_Redspotted_grouper_nervous_necrosis_virus
HQ859900_Golden_pompano_nervous_necrosis_virus
MZ449210_Redspotted_grouper_nervous_necrosis_virus
GQ402012_Redspotted_grouper_nervous_necrosis_virus
NC_024492_Senegalese_sole_Iberian_betanodavirus
JN189909_Striped_jack_nervous_necrosis_virus
JN189831_Redspotted_grouper_nervous_necrosis_virus
KY354688_Redspotted_grouper_nervous_necrosis_virus
OL955902_Striped_jack_nervous_necrosis_virus
JN189833_Striped_jack_nervous_necrosis_virus
KY354693_Redspotted_grouper_nervous_necrosis_virus
NC_013460_Tiger_puffer_nervous_necrosis_virus
OL955905_Solea_solea_betanodavirus
KX883026_Beihai_noda-like_virus_17
KX883058_Beihai_noda-like_virus_17
NC_032172_Beihai_noda-like_virus_18
KX883056_Beihai_mantis_shrimp_virus_6
KX883019_Beihai_noda-like_virus_26
KP969992_Craigies_Hill_virus
KP969988_Craigies_Hill_virus
KP969990_Craigies_Hill_virus
KP969989_Craigies_Hill_virus
KP969987_Craigies_Hill_virus
KP969981_Craigies_Hill_virus
KP969985_Craigies_Hill_virus
MH384328_Craigies_Hill_virus
MH384347_Craigies_Hill_virus
KP969980_Craigies_Hill_virus
KP969979_Craigies_Hill_virus
KP969982_Craigies_Hill_virus
KP969975_Craigies_Hill_virus
KP969978_Craigies_Hill_virus
KP969971_Craigies_Hill_virus
KP969973_Craigies_Hill_virus
KP969972_Craigies_Hill_virus
KP714084_Craigies_Hill_virus
KU754526_Craigmillar_Park_virus
MW310387_Zeugodacus_tau_noda-like_virus_isolate_Bz
MZ443596_Leuven_Nodavirus-like_1
NC_032872_Hubei_diptera_virus_16
KX883139_Hubei_noda-like_virus_13
KX883162_Hubei_noda-like_virus_13
KX883140_Hubei_noda-like_virus_15
Easynal_virus_F1383C255F1384
Easynal_virus_F1383C255F1383
Easynal_virus_F1383C255F1375
Easynal_virus_F1383C255F1385
Easynal_virus_F1383C255F1380
Easynal_virus_F1383C255F1374
OL472284_Leveillula_taurica_associated_noda-like_virus_1
OL472285_Leveillula_taurica_associated_noda-like_virus_1
OL472283_Leveillula_taurica_associated_noda-like_virus_1
OL472286_Leveillula_taurica_associated_noda-like_virus_1
OL472287_Leveillula_taurica_associated_noda-like_virus_1
HM228873_Bat_guano_associated_nodavirus_GF-4n
MN617794_Erysiphe_necator_associated_nodavirus_4
MN617791_Erysiphe_necator_associated_nodavirus_1
KX774634_Lone_star_tick_nodavirus
OM954249_Flumine_nodavirus_17
OM953985_Flumine_nodavirus_16
MT085291_Riboviria_sp_
MT085295_Riboviria_sp_
MZ209970_Sanya_nodavirus_1
NC_032647_Beihai_shrimp_virus_6
MF190050_Barns_Ness_serrated_wrack_noda-like_virus_2
MF190049_Barns_Ness_serrated_wrack_noda-like_virus_1
KF510033_Betegovirus_SF
AB083060_Sclerophthora_macrospora_virus_A
KX883034_Beihai_noda-like_virus_4
NC_033166_Wenzhou_bivalvia_virus_3
MW897065_Sichuan_sediment_noda-like_virus_4
MW897125_Sichuan_sediment_noda-like_virus_13
KX883059_Beihai_noda-like_virus_5
"
Nodamuvirales_clade1 <- str_split(Nodamuvirales_clade1_str, "\\s")[[1]]
Nodamuvirales_clade1 <- Nodamuvirales_clade1[which(!Nodamuvirales_clade1 == "")]

Nodamuvirales_clade_aa1_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/aln/aa_seq/Nodamuvirales_clade_aa.fasta",
                                seqtype = "AA")
Nodamuvirales_clade1_aa_seq <- Nodamuvirales_clade_aa1_seq[Nodamuvirales_clade1]
names(Nodamuvirales_clade1_aa_seq)
seqinr::write.fasta(sequences = Nodamuvirales_clade1_aa_seq , names =names(Nodamuvirales_clade1_aa_seq),
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/aln/aa_seq/Nodamuvirales_clade1_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "BK059724_Opistofel_virus"
new_outg_seq <- Nodamuvirales_clade_aa1_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq),
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Nodamuvirales_clade/aln/aa_seq/",
                                      new_outg, ".fasta"))

