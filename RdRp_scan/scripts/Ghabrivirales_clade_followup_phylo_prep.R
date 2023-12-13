library(tidyverse)
library(seqinr)
library(treeio)

Ghabrivirales_clade_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/Ghabrivirales_clade_novel_sequence_data.txt")

### select consensus sequences based on mapping coverage
### first round of read mapping results
read_map_df <-  read_delim("/full_path_to/wd/RdRp_scan/analysis/read_mapping/mapping_summary_bowtie2_RdRpScan_Round1_analysis/results_mapped_only.txt", 
           delim = "\t")


read_map_df_Ghabrivirales_clade <- read_map_df[which(read_map_df$ref_id %in% Ghabrivirales_clade_df$contig_id & read_map_df$selected == "yes"),]



read_map_df_Ghabrivirales_clade$RdRp_scan_map_consensus_id <- paste0(read_map_df_Ghabrivirales_clade$ref_id, "_",
                      read_map_df_Ghabrivirales_clade$sample_id, "_RdRpScan_Round1")

read_map_df_Ghabrivirales_clade$contig_id <- read_map_df_Ghabrivirales_clade$ref_id

consensus_df <- left_join(read_map_df_Ghabrivirales_clade, Ghabrivirales_clade_df, by  = "contig_id")



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
extracted_consensus_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/extracted_cons1_renamed/"
dir.create(extracted_consensus_directory)
extracted_orf_nt_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/extracted_cons1_orf_nt/"
dir.create(extracted_orf_nt_directory)
extracted_orf_aa_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/extracted_cons1_orf_aa/"
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

write.table(x = consensus_df, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/consensus_df.txt", 
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




write.table(x = novel_seq_mdf1, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/updated_consensus_df.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)




##### Manually correct frames for some of the sequences. Reimport and trim parts beyond stop codons
orf_aa_w_stops <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/Ghabrivirales_orf_aa_manually_corrected.fasta",
                                     seqtype = "AA")
extracted_orf_aa_directory1 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/extracted_cons1_orf_aa_corrected/"
dir.create(extracted_orf_aa_directory1)

for (n in names(orf_aa_w_stops)){
  s <- str_split(paste0(as.character(orf_aa_w_stops[n][[1]]), collapse = ""), "\\*")
  t <- tibble(seq = s[[1]], len  = str_length(s[[1]]))
  s1 <- as.character(t[which(t$len == max(t$len)),][1,'seq'])
  s2 <- as.SeqFastaAA(s1, name = n, Annot = n)
  seqinr::write.fasta(sequences = s2, names = n,
                      file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/extracted_cons1_orf_aa_corrected/",
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


web_blastp_df1 <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/WD71K7U7013-Alignment.txt", 
                        delim = "\t", col_names = web_blastp_colnames, comment = "#") %>% 
  select(qseqid, sseqid, pident, length)

web_blastp_df2 <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/WD72CCGW013-Alignment.txt", 
                            delim = "\t", col_names = web_blastp_colnames, comment = "#") %>% 
  select(qseqid, sseqid, pident, length)

web_blastp_df <- bind_rows(web_blastp_df1, web_blastp_df2)


sort(unique(web_blastp_df$sseqid))

#### Let's also add RdRp_scan sequences from the same clade
rdrp_scan_clade_accns <- str_remove(str_remove(str_split(consensus_df$ref_suggestions, ",")[[1]], "[[:print:]]+Ghabrivirales_"), "\\.[[:digit:]]{1}\\_[[:print:]]+$")
rdrp_scan_clade_accns <- str_remove(rdrp_scan_clade_accns, "[[:print:]]*NOT_ASSIGNED_")
rdrp_scan_clade_accns <- str_remove(rdrp_scan_clade_accns, "[[:print:]]*Durnavirales_")
rdrp_scan_clade_accns <- str_remove(rdrp_scan_clade_accns, "[[:print:]]*Patatavirales_")
rdrp_scan_clade_accns <- str_remove(rdrp_scan_clade_accns, "[[:print:]]*Botybirnavirus_")
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
NCBI_Virus_mdf <- read_csv("/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/Ghabrivirales_metadata_public.csv")


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


write.table(x = selected_mdf_upd_novel, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/updated_metadata_with_novel_seq.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)



######################## Dealing with CDS
####### Extract CDS from NCBI nucleotide, not NCBI Virus, because the headers are more informative
paste0(sort(selected_mdf$nt_accession), collapse = ",")

####### Here, because there could be multiple CDS, use protein accessions to extract CDS
####### Correspondence between the prot and nt accessions is encoded in standard CDS names when extraction from NCBI Protein same as from NCBI Nucleotide


### Only extract prot refs that are viruses based on presence in NCBI Virus, because a lot of accessions fail to be processed in batch extraction in MCBI protein
refs <- read_lines("/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/Ghabrivirales_accessions_viruses.txt")

paste0(unique(sort(refs)), collapse = ",")

length(unique(sort(refs)))






### using this accession list to NCBI and extract CDS sequences and import them here.
#### Put the file here!
cds_fasta_file1 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/Ghabrivirales_clade_refs_cds_nt.txt"
cds_fasta_file2 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/Ghabrivirales_clade_refs_cds_nt2.txt"
cds_fasta_file3 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/Ghabrivirales_clade_refs_cds_nt3.txt"
cds_fasta_file4 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/Ghabrivirales_clade_refs_cds_nt4.txt"



### There might be more automatic solutions for CDS extractions with EMBOSS or CL implementation of entrez. But web NCBI is the fastest for now.

############################################################################################
############################################################################################
############################################################################################
cds_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/ref_cds_nt"
dir.create(cds_clade_dir)
aa_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/ref_cds_aa"
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
# gb_parse_metadata_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/Ghabrivirales_to_parse_metadata.gb"
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


write.table(x = selected_mdf_upd_novel1, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/updated_metadata_with_novel_seq1.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)


length(unique(refs))


paste0(setdiff(sort(ref_df$nt_accession), sort(cds_names_df$nt_accession)), collapse = ",")
paste0(setdiff(sort(cds_names_df$prot_accession), sort(refs)), collapse = ",")
paste0(setdiff(sort(refs), sort(cds_names_df$prot_accession)), collapse = ",")



### These are viral proteins that don't have link to nuccore because they are from PDB or UniProt.
### Ignoring these.
### 6K32_A,6TY8_A,Q02119


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
#### Ghabrivirales_clade1 tree
Ghabrivirales_clade1_str <- "OP020071_Hattula_totivirus_1
OP019901_Hattula_totivirus_1
OP019870_Hattula_totivirus_1
OP019868_Hattula_totivirus_1
OP019869_Hattula_totivirus_1
OP019871_Hattula_totivirus_1
OP020070_Hattula_totivirus_1
OP020072_Hattula_totivirus_1
MN661080_Pisingos_virus
MN661078_Pisingos_virus
MF344587_Murici_virus
MN661079_Pisingos_virus
MN661081_Pisingos_virus
MN784057_Australian_Anopheles_totivirus
NC_035674_Australian_Anopheles_totivirus
MN784056_Australian_Anopheles_totivirus
MF073201_Australian_Anopheles_totivirus
MW699041_Culex_inatomii_totivirus
MW699042_Culex_inatomii_totivirus
MW520387_Culex_modestus_totivirus
LC514398_Culex_inatomii_totivirus
Culex_inatomii_totivirus_F1570C24F1570
Culex_inatomii_totivirus_F1570C24F1582
OL700195_XiangYun_toti-like_virus_6
OP264907_Culex_modestus_totivirus
MT498830_Fitzroy_Crossing_toti-like_virus_1
MK440653_Lindangsbacken_virus
OM743943_Green_Valley_Lake_virus
MH188048_Totivirus-like_Culex_mosquito_virus_1
MW434958_Totivirus-like_Culex_mosquito_virus_1
OM817549_Totivirus-like_Culex_mosquito_virus_1
MW251336_Zyryana_toti-like_virus
MZ209742_Hangzhou_totivirus_9
MN745083_Bactrocera_dorsalis_toti-like_virus_2
MN053725_Aedes_aegypti_totivirus
MN053722_Aedes_aegypti_totivirus
MN053726_Aedes_aegypti_totivirus
MN053734_Aedes_aegypti_totivirus
MN053732_Aedes_aegypti_totivirus
MN053728_Aedes_aegypti_totivirus
LC496074_Aedes_aegypti_totivirus
MN053723_Aedes_aegypti_totivirus
MN053733_Aedes_aegypti_totivirus
MN053730_Aedes_aegypti_totivirus
MN053727_Aedes_aegypti_totivirus
MN053729_Aedes_aegypti_totivirus
MW434940_Aedes_aegypti_totivirus
MT435499_Aedes_aegypti_totivirus_2
MZ556103_Lactea_totivirus
ON860472_Koversta_virus
MW520394_Drama_totivirus
MN661076_Embera_virus
"
Ghabrivirales_clade1 <- str_split(Ghabrivirales_clade1_str, "\\s")[[1]]
Ghabrivirales_clade1 <- Ghabrivirales_clade1[which(!Ghabrivirales_clade1 == "")]

Ghabrivirales_clade_aa1_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/Ghabrivirales_clade_aa1.fasta",
                                seqtype = "AA")
Ghabrivirales_clade1_aa_seq <- Ghabrivirales_clade_aa1_seq[Ghabrivirales_clade1]
names(Ghabrivirales_clade1_aa_seq)
seqinr::write.fasta(sequences = Ghabrivirales_clade1_aa_seq , names =names(Ghabrivirales_clade1_aa_seq),
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/Ghabrivirales_clade1_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "MZ209749_Fushun_totivirus_5"
new_outg_seq <- Ghabrivirales_clade_aa1_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq),
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/",
                                      new_outg, ".fasta"))



####### Zooming in on different clades
#### Ghabrivirales_clade2 tree
Ghabrivirales_clade2_str <- "Andiett_virus_F1132C2F1132
Andbire_virus_F1046C18F1046
MZ852358_Drosophila-associated_totivirus_5
MT129768_Keenan_toti-like_virus
MZ852363_Drosophila-associated_totivirus_4
OP020048_Hanko_toti-like_virus_1
BK059880_Aloadae_toti-like_virus_2
"
Ghabrivirales_clade2 <- str_split(Ghabrivirales_clade2_str, "\\s")[[1]]
Ghabrivirales_clade2 <- Ghabrivirales_clade2[which(!Ghabrivirales_clade2 == "")]

Ghabrivirales_clade_aa1_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/Ghabrivirales_clade_aa1.fasta",
                                              seqtype = "AA")
Ghabrivirales_clade2_aa_seq <- Ghabrivirales_clade_aa1_seq[Ghabrivirales_clade2]
names(Ghabrivirales_clade2_aa_seq)
seqinr::write.fasta(sequences = Ghabrivirales_clade2_aa_seq , names =names(Ghabrivirales_clade2_aa_seq),
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/Ghabrivirales_clade2_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "ON746544_Jiamusi_Totiv_tick_virus_1"
new_outg_seq <- Ghabrivirales_clade_aa1_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq),
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/",
                                      new_outg, ".fasta"))





####### Zooming in on different clades
#### Ghabrivirales_clade3 tree
Ghabrivirales_clade3_str <- "MF176261_Hubei_chryso-like_virus_1
MF176280_Hubei_chryso-like_virus_1
MF176388_Hubei_chryso-like_virus_1
MH085092_Hubei_chryso-like_virus_1
MW452288_Chrysoviridae_sp_
MW434089_Hubei_chryso-like_virus_1
OL700058_Hubei_chryso-like_virus_1
KX882962_Hubei_chryso-like_virus_1
OK413943_Hubei_chryso-like_virus_1
Hubei_chryso_like_virus_1_F1427C12F1427
Hubei_chryso_like_virus_1_F1427C12F1431
Hubei_chryso_like_virus_1_F1427C12F1425
LC514396_Hubei_chryso-like_virus_1
MW520399_Soufli_chryso-like_virus
OP264884_Eskilstorp_virus
MK440655_Eskilstorp_virus
MW520400_Xanthi_chryso-like_virus
MT498826_Broome_chryso-like_virus_1
KX882959_Alphachrysovirus_shuangaoense
NC_055224_Alphachrysovirus_shuangaoense
MH745159_Alphachrysovirus_shuangaoense
KX882964_Alphachrysovirus_shuangaoense
MW520414_Evros_chryso-like_virus
OP019910_Lestijarvi_alphachrysovirus
OP019845_Lestijarvi_alphachrysovirus
OP019911_Lestijarvi_alphachrysovirus
OP019837_Enontekio_alphachrysovirus
OP019841_Hanko_alphachrysovirus
MW434094_Keturi_virus
MN661047_Alphachrysovirus_saladoense
KX882973_Hubei_chryso-like_virus_2
"
Ghabrivirales_clade3 <- str_split(Ghabrivirales_clade3_str, "\\s")[[1]]
Ghabrivirales_clade3 <- Ghabrivirales_clade3[which(!Ghabrivirales_clade3 == "")]

Ghabrivirales_clade_aa1_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/Ghabrivirales_clade_aa1.fasta",
                                                 seqtype = "AA")
Ghabrivirales_clade3_aa_seq <- Ghabrivirales_clade_aa1_seq[Ghabrivirales_clade3]
names(Ghabrivirales_clade3_aa_seq)
seqinr::write.fasta(sequences = Ghabrivirales_clade3_aa_seq , names =names(Ghabrivirales_clade3_aa_seq),
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/Ghabrivirales_clade3_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "LC516853_Barley_aphid_RNA_virus_8"
new_outg_seq <- Ghabrivirales_clade_aa1_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq),
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/",
                                      new_outg, ".fasta"))



####### Zooming in on different clades
#### Ghabrivirales_clade4 tree
Ghabrivirales_clade4_str <- "ON812656_Totiviridae_sp_
ON812661_Totiviridae_sp_
ON746559_Dali_Totiv_tick_virus_1
MN617029_Erysiphe_necator_associated_totivirus_8
ON812658_Erysiphe_necator_associated_totivirus_8
ON812740_Totiviridae_sp_
MN812428_Malassezia_sympodialis_mycovirus
MN603497_Malassezia_restricta_virus_MrV40L
NC_038697_Scheffersomyces_segobiensis_virus_L
MZ218643_Totiviridae_sp_
Fanomonover_virus_F1497C15F1497
MN617027_Erysiphe_necator_associated_totivirus_6
MN617028_Erysiphe_necator_associated_totivirus_7
ON812681_Totiviridae_sp_
ON812741_Totiviridae_sp_
ON812606_Totiviridae_sp_
ON812795_Totiviridae_sp_
OP020052_Hanko_totivirus_2
ON746555_Xian_Totiv_tick_virus_1
ON812703_Totiviridae_sp_
OP080565_Tuatara_cloaca-assoiciated_totivrus-3
MN628284_Erysiphales_associated_totivirus_13
MN628264_Erysiphales_associated_mitovirus_1
MN628272_Erysiphales_associated_totivirus_1
MZ218629_Totiviridae_sp_
MZ218553_Totiviridae_sp_
MZ218642_Totiviridae_sp_
KT455458_Delisea_pulchra_totivirus_IndA
KT455457_Delisea_pulchra_totivirus_IndA
KT455451_Delisea_pulchra_totivirus_IndA
MZ600515_Conidiobolus_lamprauges_totivirus_2
MZ218554_Totiviridae_sp_
ON862169_Geotrichum_candidum_totivirus_3a
ON862170_Geotrichum_candidum_totivirus_3b
MZ218667_Totiviridae_sp_
ON746553_Luoyang_Totiv_tick_virus_2
KT784813_Saccharomyces_cerevisiae_virus_L-BC-lus
KX906605_Saccharomyces_cerevisiae_virus_LBCLa
NC_001641_Saccharomyces_cerevisiae_virus_LBCLa
ON746560_Luoyang_Totiv_tick_virus_3
ON862167_Geotrichum_candidum_totivirus_1
ON812734_Totiviridae_sp_
MZ218646_Totiviridae_sp_
MZ218664_Totiviridae_sp_
OW529226_Rhodochaete_parvula_toti-like_virus_1
MZ218584_Totiviridae_sp_
ON746556_Xian_Totiv_tick_virus_2
NC_028481_Red_clover_powdery_mildew-associated_totivirus_2
"
Ghabrivirales_clade4 <- str_split(Ghabrivirales_clade4_str, "\\s")[[1]]
Ghabrivirales_clade4 <- Ghabrivirales_clade4[which(!Ghabrivirales_clade4 == "")]

Ghabrivirales_clade_aa1_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/Ghabrivirales_clade_aa1.fasta",
                                                 seqtype = "AA")
Ghabrivirales_clade4_aa_seq <- Ghabrivirales_clade_aa1_seq[Ghabrivirales_clade4]
names(Ghabrivirales_clade4_aa_seq)
seqinr::write.fasta(sequences = Ghabrivirales_clade4_aa_seq , names =names(Ghabrivirales_clade4_aa_seq),
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/Ghabrivirales_clade4_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "OW529222_Porphyridium_purpureum_toti-like_virus_1"
new_outg_seq <- Ghabrivirales_clade_aa1_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq),
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/",
                                      new_outg, ".fasta"))


####### Zooming in on different clades
#### Ghabrivirales_clade5 tree
Ghabrivirales_clade5_str <- "ON812600_Totiviridae_sp_
ON812599_Totiviridae_sp_
ON746554_Daqing_Totiv_tick_virus_1
MZ218677_Totiviridae_sp_
MZ209960_Sanya_ochthera_mantis_totivirus_1
OL700199_XiangYun_toti-like_virus_10
Ansforinte_virus_F1367C12F1360
Ansforinte_virus_F1367C12F1367
NC_032465_Beihai_blue_swimmer_crab_virus_3
"
Ghabrivirales_clade5 <- str_split(Ghabrivirales_clade5_str, "\\s")[[1]]
Ghabrivirales_clade5 <- Ghabrivirales_clade5[which(!Ghabrivirales_clade5 == "")]

Ghabrivirales_clade_aa1_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/Ghabrivirales_clade_aa1.fasta",
                                                 seqtype = "AA")
Ghabrivirales_clade5_aa_seq <- Ghabrivirales_clade_aa1_seq[Ghabrivirales_clade5]
names(Ghabrivirales_clade5_aa_seq)
seqinr::write.fasta(sequences = Ghabrivirales_clade5_aa_seq , names =names(Ghabrivirales_clade5_aa_seq),
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/Ghabrivirales_clade5_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "NC_032465_Beihai_blue_swimmer_crab_virus_3"
new_outg_seq <- Ghabrivirales_clade_aa1_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq),
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/",
                                      new_outg, ".fasta"))




####### Zooming in on different clades
#### Ghabrivirales_clade6 tree
Ghabrivirales_clade6_str <- "MW588192_Mute_swan_feces_associated_toti-like_virus_1
MW588193_Mute_swan_feces_associated_toti-like_virus_1
MW588191_Mute_swan_feces_associated_toti-like_virus_1
MW588194_Mute_swan_feces_associated_toti-like_virus_2
MW588195_Mute_swan_feces_associated_toti-like_virus_3
MW588197_Mute_swan_feces_associated_toti-like_virus_3
MW588196_Mute_swan_feces_associated_toti-like_virus_3
Vistafter_virus_F1409C1832F1409
MN661074_Atrato_virus
MN661077_Murri_virus
MN803436_Schistocephalus_solidus_toti-like_virus_2
MW826407_Totiviridae_sp_
OP020076_Hattula_totivirus_2
OP020073_Hattula_totivirus_2
OP019873_Hattula_totivirus_2
OP020077_Hattula_totivirus_2
OP020079_Hattula_totivirus_2
OP020078_Hattula_totivirus_2
OP019872_Hattula_totivirus_2
OP019874_Hattula_totivirus_2
OP020080_Hattula_totivirus_2
OP020075_Hattula_totivirus_2
OP019875_Hattula_totivirus_2
OP019876_Hattula_totivirus_2
KP642123_dsRNA_virus_environmental_sample
OP019850_Enontekio_totivirus_2
OP019990_Enontekio_totivirus_2
OP019889_Joutseno_totivirus
OP019864_Hanko_totivirus_8
OP019865_Hanko_totivirus_8
OP020129_Vaasa_totivirus
KP642121_dsRNA_virus_environmental_sample
KP642124_dsRNA_virus_environmental_sample
Wrunlize_virus_F1469C738F1469
Wrunlize_virus_F1469C738F1615
Aftionautor_virus_F1470C2F1470
Aftionautor_virus_F1470C2F1471
Coundenality_virus_F1183C376F1183
Adsentalan_virus_F1186C10F1177
Adsentalan_virus_F1186C10F1186
Adsentalan_virus_F1186C10F1187
Flatesals_virus_F1187C10F1178
Flatesals_virus_F1187C10F1187
Flatesals_virus_F1187C10F1177
Adjoured_virus_F1191C7F1191
Adjoured_virus_F1191C7F1194
Aedes_aegypti_toti_like_virus_F1063C6F1002
Aedes_aegypti_toti_like_virus_F1063C6F1063
KP642128_dsRNA_virus_environmental_sample
MN053720_Aedes_aegypti_toti-like_virus
MN053721_Aedes_aegypti_toti-like_virus
OL343758_Palmetto_toti-like_virus
Adlissaust_virus_F1183C2F1183
Belfgalbs_virus_F1490C6F1490
NC_032733_Hubei_toti-like_virus_10
Agrichlet_virus_F1434C1F1434
Agrichlet_virus_F1434C1F1456
Agrichlet_virus_F1434C1F1407
Agrichlet_virus_F1434C1F1445
LC567881_Culex_vishnui_subgroup_totivirus
LC514295_Culex_vishnui_subgroup_totivirus
Culex_vishnui_subgroup_totivirus_F1143C11F1143
MK440659_Salja_virus
MK440645_Osta_virus
MW434960_Tzifr_virus
MW434962_Tzifr_virus
MW434963_Tzifr_virus
MW434961_Tzifr_virus
MW434964_Tzifr_virus
MW434959_Tzifr_virus
MW434955_Snelk_virus
MW434957_Stinn_virus
Aeriandims_virus_subsp1_F1156C2F1045
Aeriandims_virus_subsp1_F1156C2F1142
Aeriandims_virus_subsp1_F1156C2F1156
Aeriandims_virus_subsp2_F1583C1F1583
OL700194_XiangYun_toti-like_virus_5
MT498828_Fitzroy_Crossing_toti-like_virus_2
OL700198_XiangYun_toti-like_virus_9
MW434949_Lotchka_virus
MW434951_Lotchka_virus
MW434950_Lotchka_virus
MW434945_Gouley_virus
MW434952_Mika_virus
OP020045_Hameenlinna_totivirus_3
OP020042_Hameenlinna_totivirus_3
OP020035_Hameenlinna_totivirus_3
OP020032_Hameenlinna_totivirus_3
OP019859_Hameenlinna_totivirus_3
OP020037_Hameenlinna_totivirus_3
OP020040_Hameenlinna_totivirus_3
OP020069_Hanko_totivirus_9
OP019900_Hanko_totivirus_9
OP019866_Hanko_totivirus_9
OP019867_Hanko_totivirus_9
OP019991_Enontekio_totivirus_3
KP642122_dsRNA_virus_environmental_sample
KP642125_dsRNA_virus_environmental_sample
OP019852_Hameenlinna_totivirus_1
OP019854_Hameenlinna_totivirus_1
OP019855_Hameenlinna_totivirus_1
OP019856_Hameenlinna_totivirus_1
OP020003_Hameenlinna_totivirus_1
OP020014_Hameenlinna_totivirus_1
OP020018_Hameenlinna_totivirus_1
OP020026_Hameenlinna_totivirus_1
OP020002_Hameenlinna_totivirus_1
OP019853_Hameenlinna_totivirus_1
OP019857_Hameenlinna_totivirus_1
OP020020_Hameenlinna_totivirus_1
OP020006_Hameenlinna_totivirus_1
OP020007_Hameenlinna_totivirus_1
OP020016_Hameenlinna_totivirus_1
OP019858_Hameenlinna_totivirus_1
OP020011_Hameenlinna_totivirus_1
OP020015_Hameenlinna_totivirus_1
OP020010_Hameenlinna_totivirus_1
OP020005_Hameenlinna_totivirus_1
OP020021_Hameenlinna_totivirus_1
OP020027_Hameenlinna_totivirus_1
OP019885_Inari_totivirus_1
OP019884_Ilomantsi_totivirus_1
OP020092_Ilomantsi_totivirus_1
OP019886_Inari_totivirus_2
OP019887_Inari_totivirus_2
OP019888_Inari_totivirus_2
OP020067_Hanko_totivirus_6
OP019892_Lestijarvi_totivirus
OP019893_Lestijarvi_totivirus
OP019894_Lestijarvi_totivirus
OP019861_Hanko_totivirus_3
OP019862_Hanko_totivirus_3
OP019902_Hanko_totivirus_3
OP019903_Hanko_totivirus_3
OP019904_Hanko_totivirus_3
OP019909_Hanko_totivirus_3
OP020054_Hanko_totivirus_3
OP019905_Hanko_totivirus_3
OP020053_Hanko_totivirus_3
KF298277_uncultured_virus
OP020060_Hanko_totivirus_5
OP020058_Hanko_totivirus_5
OP020066_Hanko_totivirus_5
OP020059_Hanko_totivirus_5
OP020064_Hanko_totivirus_5
OP020061_Hanko_totivirus_5
OP020063_Hanko_totivirus_5
OP019863_Hanko_totivirus_5
OP020065_Hanko_totivirus_5
OP019907_Palkane_totivirus
OP019908_Palkane_totivirus
OP020114_Kuusamo_totivirus_1
OP020030_Hameenlinna_totivirus_2
OP020031_Hameenlinna_totivirus_2
OP020029_Hameenlinna_totivirus_2
OP020028_Hameenlinna_totivirus_2
OP020068_Hanko_totivirus_7
OP020095_Ilomantsi_totivirus_2
NC_035130_Aedes_alboannulatus_toti-like_virus_1
NC_035131_Aedes_camptorhynchus_toti-like_virus_1
"
Ghabrivirales_clade6 <- str_split(Ghabrivirales_clade6_str, "\\s")[[1]]
Ghabrivirales_clade6 <- Ghabrivirales_clade6[which(!Ghabrivirales_clade6 == "")]

Ghabrivirales_clade_aa1_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/Ghabrivirales_clade_aa1.fasta",
                                                 seqtype = "AA")
Ghabrivirales_clade6_aa_seq <- Ghabrivirales_clade_aa1_seq[Ghabrivirales_clade6]
names(Ghabrivirales_clade6_aa_seq)
seqinr::write.fasta(sequences = Ghabrivirales_clade6_aa_seq , names =names(Ghabrivirales_clade6_aa_seq),
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/Ghabrivirales_clade6_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "MZ218558_Totiviridae_sp_"
new_outg_seq <- Ghabrivirales_clade_aa1_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq),
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Ghabrivirales_clade/aln/aa_seq/",
                                      new_outg, ".fasta"))




