library(tidyverse)
library(seqinr)
library(treeio)

Monjiviricetes_clade_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/Monjiviricetes_clade_novel_sequence_data.txt")

### select consensus sequences based on mapping coverage
### first round of read mapping results
read_map_df <-  read_delim("/full_path_to/wd/RdRp_scan/analysis/read_mapping/mapping_summary_bowtie2_RdRpScan_Round1_analysis/results_mapped_only.txt", 
           delim = "\t")


read_map_df_Monjiviricetes_clade <- read_map_df[which(read_map_df$ref_id %in% Monjiviricetes_clade_df$contig_id & read_map_df$selected == "yes"),]



read_map_df_Monjiviricetes_clade$RdRp_scan_map_consensus_id <- paste0(read_map_df_Monjiviricetes_clade$ref_id, "_",
                      read_map_df_Monjiviricetes_clade$sample_id, "_RdRpScan_Round1")

read_map_df_Monjiviricetes_clade$contig_id <- read_map_df_Monjiviricetes_clade$ref_id

consensus_df <- left_join(read_map_df_Monjiviricetes_clade, Monjiviricetes_clade_df, by  = "contig_id")



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
extracted_consensus_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/extracted_cons1_renamed/"
dir.create(extracted_consensus_directory)
extracted_orf_nt_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/extracted_cons1_orf_nt/"
dir.create(extracted_orf_nt_directory)
extracted_orf_aa_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/extracted_cons1_orf_aa/"
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

write.table(x = consensus_df, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/consensus_df.txt", 
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




write.table(x = novel_seq_mdf1, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/updated_consensus_df.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)


##### Manually correct frames for some of the sequences. Reimport and trim parts beyond stop codons
orf_aa_w_stops <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/Monjiviriceres_orf_aa_manually_corrected.fasta",
                                     seqtype = "AA")
extracted_orf_aa_directory1 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/extracted_cons1_orf_aa_corrected/"
dir.create(extracted_orf_aa_directory1)

for (n in names(orf_aa_w_stops)){
  s <- str_split(paste0(as.character(orf_aa_w_stops[n][[1]]), collapse = ""), "\\*")
  t <- tibble(seq = s[[1]], len  = str_length(s[[1]]))
  s1 <- as.character(t[which(t$len == max(t$len)),][1,'seq']) 
  s2 <- as.SeqFastaAA(s1, name = n, Annot = n)
  seqinr::write.fasta(sequences = s2, names = n, 
                      file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/extracted_cons1_orf_aa_corrected/",
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


web_blastp_df1 <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/U0S5CDGC013-Alignment.txt", 
                        delim = "\t", col_names = web_blastp_colnames, comment = "#") %>% 
  select(qseqid, sseqid, pident, length)

web_blastp_df2 <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/U0S5VUAF013-Alignment.txt", 
                            delim = "\t", col_names = web_blastp_colnames, comment = "#") %>% 
  select(qseqid, sseqid, pident, length)

web_blastp_df <- bind_rows(web_blastp_df1, web_blastp_df2)

sort(unique(web_blastp_df$sseqid))

#### Let's also add RdRp_scan sequences from the same clade
rdrp_scan_clade_accns <- str_remove(str_remove(str_split(consensus_df$ref_suggestions, ",")[[1]], "[[:print:]]+Monjiviricetes_"), "\\.[[:digit:]]{1}\\_[[:print:]]+$")
## remove one strange sequence
rdrp_scan_clade_accns <- rdrp_scan_clade_accns[which(!rdrp_scan_clade_accns == " Lenarviricota_Wolframvirales_QDH88671")]
## Keep not assigned as well as other classes Milneviricetes, Chunqiuviricetes, Yunchangviricetes, which will be removed in later trees
rdrp_scan_clade_accns <- str_remove(str_remove(rdrp_scan_clade_accns, "[[:print:]]+NOT_ASSIGNED_"), "\\.[[:digit:]]{1}\\_[[:print:]]+$")
rdrp_scan_clade_accns <- str_remove(str_remove(rdrp_scan_clade_accns, "[[:print:]]+Milneviricetes_"), "\\.[[:digit:]]{1}\\_[[:print:]]+$")
rdrp_scan_clade_accns <- str_remove(str_remove(rdrp_scan_clade_accns, "[[:print:]]+Chunqiuviricetes_"), "\\.[[:digit:]]{1}\\_[[:print:]]+$")
rdrp_scan_clade_accns <- str_remove(str_remove(rdrp_scan_clade_accns, "[[:print:]]+Yunchangviricetes_"), "\\.[[:digit:]]{1}\\_[[:print:]]+$")

rdrp_scan_clade_accns[which(str_detect(rdrp_scan_clade_accns, "\\|"))] ### these are from pdb, find equivalents in nr prot
rdrp_scan_clade_accns <- rdrp_scan_clade_accns[which(!str_detect(rdrp_scan_clade_accns, "\\|"))]
rdrp_scan_clade_accns <- rdrp_scan_clade_accns[which(!rdrp_scan_clade_accns == "Xinjiang_varicosa_like_virus")] ## substitute this with the corresponding accession QYF49869
rdrp_scan_clade_accns <- sort(c(rdrp_scan_clade_accns, "QYW14734", "UYO79445", "NP_056797", "QYF49869"))


refs <- sort(unique(c(blastp_df2$sseqid, rdrp_scan_clade_accns, web_blastp_df$sseqid)))
refs <- str_remove(refs, "\\.[[:digit:]]$")

length(sort(unique(refs)))



##### Check on NCBI Viruses
paste0(sort(unique(refs)), collapse = ", ")


#### Reimport metadata
NCBI_Virus_mdf <- read_csv("/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/Monjiviricetes_metadata_public.csv")


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


write.table(x = selected_mdf_upd_novel, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/updated_metadata_with_novel_seq.txt", 
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
cds_fasta_file1 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/Monjiviricetes_clade_refs_cds_nt.txt"
cds_fasta_file2 <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/Monjiviricetes_clade_refs_cds_nt2.txt"

### There might be more automatic solutions for CDS extractions with EMBOSS or CL implementation of entrez. But web NCBI is the fastest for now.

############################################################################################
############################################################################################
############################################################################################
cds_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/ref_cds_nt"
dir.create(cds_clade_dir)
aa_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/ref_cds_aa"
dir.create(aa_clade_dir)


cds_fasta1 <- seqinr::read.fasta(file = cds_fasta_file1, seqtype = "DNA", whole.header = T)
cds_fasta2 <- seqinr::read.fasta(file = cds_fasta_file2, seqtype = "DNA", whole.header = T)

names(cds_fasta1)
length(selected_mdf$nt_accession)
length(names(cds_fasta1))
length(names(cds_fasta2))

length(names(cds_fasta1)) + length(names(cds_fasta2))


cds_fasta <- c(cds_fasta1, cds_fasta2)



cds_names_df <- tibble(original_cds_headers = names(cds_fasta),
                       nt_accession = str_remove(str_extract(names(cds_fasta), "(?<=^lcl\\|)[[:print:]]+(?=_cds)"), "\\.[[:digit:]]*"),
                       cds_product = str_replace_all(str_extract(names(cds_fasta), "(?<=\\[protein\\=)[[:print:]]+?(?=\\])"), " ", "_"),
                       prot_accession = str_remove(str_extract(names(cds_fasta), "(?<=_cds_)[[:alnum:]\\_\\.]+(?= )"), "\\.[[:digit:]]*\\_[[:digit:]]*")
)


ref_df <- selected_mdf
paste0(setdiff(sort(cds_names_df$nt_accession), sort(ref_df$nt_accession)), collapse = ",") ### these seem to be endogenous elements

### Reimport genbank files and get metadata from them. They are not in NCBI Virus perhaps because they are MAG: ?
gb_parse_metadata_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/Monjiviricetes_to_parse_metadata.gb"
count = 0
entry_ends = FALSE
for(i in 1:length(read_lines(gb_parse_metadata_file))){
  entry_ends = FALSE
  l <-  read_lines(gb_parse_metadata_file)[i]
  # print(l)
  if (str_detect(l, "LOCUS") == TRUE){
    entry_ends = FALSE
    nt_accession <- str_remove(str_extract(l, "LOCUS\\s+[[:alnum:]\\_]+"), "LOCUS\\s*")
    nt_slen <- str_remove(str_extract(l, "(?<=\\s)[[:digit:]]+(?=\\s+bp)"), "\\s*")
    nt_title <- str_remove(read_lines(gb_parse_metadata_file)[i+1], "DEFINITION\\s*")
  }
  
  if (str_detect(l, "\\/host\\=") == TRUE){
    host_taxon <- str_remove_all(str_remove(l, "\\s*\\/host\\="), '"')
  }
  
  if (str_detect(l, "\\/organism\\=") == TRUE){
    Species <- str_remove_all(str_remove(l, "\\s*\\/organism\\="), '"')
  }
  
  if (str_detect(l, "\\/country\\=") == TRUE){
    Location <- str_remove_all(str_remove(l, "\\s*\\/country\\="), '"')
  }
  
  if (str_detect(l, "\\/collection_date\\=") == TRUE){
    Year <- str_remove_all(str_remove(l, "\\s*\\/collection_date\\="), '"')
  }

  if (str_detect(l, "^\\/\\/") == TRUE){
    entry_ends = TRUE
  }
  
  if (entry_ends == TRUE){
    ## save row and add count
    count = count + 1
    metadata_row <- tibble(nt_accession, nt_slen, nt_title, host_taxon, Location, Year, Species)
    
    if (count == 1){
      metadata <- metadata_row
    }else{
      metadata <- bind_rows(metadata, metadata_row)
    }
  }
}

metadata$nt_seq_ID <- paste0(metadata$nt_accession, "_", str_replace_all(metadata$Species, " ", "_"))
metadata$nt_slen <- as.numeric(as.character(metadata$nt_slen))
ref_df <- bind_rows(ref_df, metadata)



unique(na.omit(str_extract(ref_df$nt_seq_ID, "[\\.\\+\\*\\?\\^\\$\\(\\)\\[\\]\\{\\}\\|\\\\\\/]")))
ref_df$nt_seq_ID <- str_replace_all(ref_df$nt_seq_ID, "[\\.\\/]", "_")

selected_mdf_upd_novel1 <- bind_rows(selected_mdf_upd_novel, metadata)
selected_mdf_upd_novel1$nt_seq_ID <- str_replace_all(selected_mdf_upd_novel1$nt_seq_ID, "[\\.\\/]", "_")


write.table(x = selected_mdf_upd_novel1, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/updated_metadata_with_novel_seq1.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)





paste0(setdiff(sort(ref_df$nt_accession), sort(cds_names_df$nt_accession)), collapse = ",")

paste0(setdiff(sort(refs), sort(cds_names_df$prot_accession)), collapse = ",")


### these are the proteins for which I couldn't extract CDS from NCBI Protein in a batch
## "P13179,P16379,Q8B0H0,Q8B0H5,Q98776"
### but these are mostly VSV and Chandipura virus, and I have many other VSV and Chandipura virus sequences in the ref data already
ref_df[which(str_detect(ref_df$nt_title, "Chandipura virus")),] ## check 
ref_df[which(str_detect(ref_df$nt_title, "Vesicular stomatitis")),] ## check 





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
#### Monjiviricetes_clade1 tree
Monjiviricetes_clade1_str <- "MW434102_Culicidavirus_quitotaense
MW434103_Culicidavirus_quitotaense
MH188031_Culicidavirus_quitotaense
MW434101_Culicidavirus_quitotaense
MW434104_Culicidavirus_quitotaense
MW434100_Culicidavirus_quitotaense
MK440640_Culicidavirus_quitotaense
MH188036_Culicidavirus_culicis
MW452291_Culicidavirus_culicidae
NC_028265_Culicidavirus_culicidae
Imjin_River_virus_1_F1051C1F1051
NC_028482_Culicidavirus_imjinense
ON955149_Hattula_chuvirus
ON955153_Hattula_chuvirus
ON955154_Hattula_chuvirus
ON955148_Hattula_chuvirus
ON955147_Hattula_chuvirus
ON955152_Hattula_chuvirus
MN661024_Atrato_Chu-like_virus_1
MN661026_Atrato_Chu-like_virus_1
KX924630_Doliuvirus_culisetae
MZ078292_Turkana_Chu-like_virus
"
Monjiviricetes_clade1 <- str_split(Monjiviricetes_clade1_str, "\\s")[[1]]
Monjiviricetes_clade1 <- Monjiviricetes_clade1[which(!Monjiviricetes_clade1 == "")]

Monjiviricetes_clade_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/Monjiviricetes_clade_aa.fasta",
                                seqtype = "AA")
Monjiviricetes_clade1_aa_seq <- Monjiviricetes_clade_aa_seq[Monjiviricetes_clade1]
names(Monjiviricetes_clade1_aa_seq)
seqinr::write.fasta(sequences = Monjiviricetes_clade1_aa_seq , names =names(Monjiviricetes_clade1_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/Monjiviricetes_clade1_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "KM817612_Shuangao_Fly_Virus_1"
new_outg_seq <- Monjiviricetes_clade_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))




####### Zooming in on different clades
#### Monjiviricetes_clade2 tree
Monjiviricetes_clade2_str <- "MN053740_Aedes_aegypti_anphevirus
MN053741_Aedes_aegypti_anphevirus
MN053739_Aedes_aegypti_anphevirus
MN053738_Aedes_aegypti_anphevirus
MN053744_Aedes_aegypti_anphevirus
MH037149_Gylbovirus_aagae
MH430663_Gylbovirus_aagae
MN053737_Aedes_aegypti_anphevirus
MN053742_Aedes_aegypti_anphevirus
MN053743_Aedes_aegypti_anphevirus
MH237595_Aedes_aegypti_anphevirus
MG012486_Aedes_aegypti_anphevirus
MH430650_Gylbovirus_aagae
MH430648_Gylbovirus_aagae
MW435012_Gylbovirus_aagae
MW435013_Gylbovirus_aagae
MH430658_Gylbovirus_aagae
Aedes_anphevirus_F1058C3F1058
MH430655_Gylbovirus_aagae
MH430666_Gylbovirus_aagae
Aedes_anphevirus_F1058C3F1059
Aedes_anphevirus_F1058C3F1060
MH430653_Gylbovirus_aagae
MH430651_Gylbovirus_aagae
Aedes_anphevirus_F1058C3F1067
Aedes_anphevirus_F1058C3F1068
Aedes_anphevirus_F1058C3F1062
Aedes_anphevirus_F1058C3F1061
MH430654_Gylbovirus_aagae
MH430652_Gylbovirus_aagae
MH430656_Gylbovirus_aagae
MH430657_Gylbovirus_aagae
MH430665_Gylbovirus_aagae
ON955256_Enontekio_anphevirus_2
MF176245_Culex_mononega-like_virus_1
MF176316_Culex_mononega-like_virus_1
MF176296_Culex_mononega-like_virus_1
MF176356_Culex_mononega-like_virus_1
MF176375_Culex_mononega-like_virus_1
MK440643_Culex_mononega-like_virus_1
OL700050_Culex_mononega-like_virus_1
MW435017_Gordis_virus
MW435018_Gordis_virus
MW435019_Gordis_virus
MW435020_Gordis_virus
MW435014_Gordis_virus
MW435015_Gordis_virus
MW435016_Gordis_virus
MW452303_Anphevirus_sp_
OL700094_Anphevirus_xinchengense
NC_031244_Anphevirus_xinchengense
Xincheng_Mosquito_Virus_F1356C1F1356
Xincheng_Mosquito_Virus_F1356C1F1367
KX148551_Bolahun_virus_variant_1
NC_055112_Gambievirus_bolahunense
KX148553_Gambievirus_senegalense
MH822963_Madalivirus_amazonaense
MH822964_Madalivirus_amazonaense
MH822965_Madalivirus_amapaense
LC514054_Triniovirus_yonagoense
LC514055_Triniovirus_yonagoense
OL700126_XiangYun_mono-chu-like_virus_1
ON860459_Malby_virus
ON955255_Enontekio_anphevirus_1
ON955258_Joensuu_anphevirus
ON955260_Joensuu_anphevirus
ON955261_Joensuu_anphevirus
ON955263_Joensuu_anphevirus
MN053735_Guadeloupe_mosquito_mononega-like_virus
Aedes_albopictus_anphevirus_F1016C1F1016
MW147277_Aedes_albopictus_anphevirus
MT822181_Serbia_mononega-like_virus_1
ON955257_Hanko_anphevirus
MF176332_Doupovirus_australiaense
NC_035133_Doupovirus_australiaense
MF176318_Doupovirus_australiaense
MF176268_Doupovirus_australiaense
MK440623_Doupovirus_australiaense
MN513369_Doupovirus_australiaense
MF176247_Doupovirus_australiaense
OL700051_Doupovirus_australiaense
Prepisiers_virus_F1032C1F1032
MW452297_Mononegavirales_sp_
OL700130_XiangYun_mono-chu-like_virus_5
OL700133_XiangYun_mono-chu-like_virus_8
MN053736_Guadeloupe_mosquito_mononega-like_virus
ON059793_Tolviot_virus
ON059795_Tolviot_virus
ON059794_Tolviot_virus
ON059796_Tolviot_virus
"
Monjiviricetes_clade2 <- str_split(Monjiviricetes_clade2_str, "\\s")[[1]]
Monjiviricetes_clade2 <- Monjiviricetes_clade2[which(!Monjiviricetes_clade2 == "")]

Monjiviricetes_clade_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/Monjiviricetes_clade_aa.fasta",
                                                  seqtype = "AA")
Monjiviricetes_clade2_aa_seq <- Monjiviricetes_clade_aa_seq[Monjiviricetes_clade2]
names(Monjiviricetes_clade2_aa_seq)
seqinr::write.fasta(sequences = Monjiviricetes_clade2_aa_seq , names =names(Monjiviricetes_clade2_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/Monjiviricetes_clade2_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "ON872562_Bat_faecal_associated_anphe-like_virus_1"
new_outg_seq <- Monjiviricetes_clade_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))



####### Zooming in on different clades
#### Monjiviricetes_clade3 tree
Monjiviricetes_clade3_str <- "MN013386_Guadeloupe_Culex_rhabdovirus
MN013387_Guadeloupe_Culex_rhabdovirus
MN013390_Guadeloupe_Culex_rhabdovirus
MN013393_Guadeloupe_Culex_rhabdovirus
MN013389_Guadeloupe_Culex_rhabdovirus
MN013388_Guadeloupe_Culex_rhabdovirus
MN013392_Guadeloupe_Culex_rhabdovirus
MG385079_Grenada_mosquito_rhabdovirus_1
MH707450_Grenada_mosquito_rhabdovirus_1
Grenada_mosquito_rhabdovirus_1_F1429C1F1422
Grenada_mosquito_rhabdovirus_1_F1429C1F1426
Grenada_mosquito_rhabdovirus_1_F1429C1F1420
Grenada_mosquito_rhabdovirus_1_F1429C1F1423
Grenada_mosquito_rhabdovirus_1_F1429C1F1428
Grenada_mosquito_rhabdovirus_1_F1429C1F1042
Grenada_mosquito_rhabdovirus_1_F1429C1F1572
Grenada_mosquito_rhabdovirus_1_F1429C1F1579
Grenada_mosquito_rhabdovirus_1_F1429C1F1040
Grenada_mosquito_rhabdovirus_1_F1429C1F1432
Grenada_mosquito_rhabdovirus_1_F1429C1F1425
Grenada_mosquito_rhabdovirus_1_F1429C1F1575
Grenada_mosquito_rhabdovirus_1_F1429C1F1573
Grenada_mosquito_rhabdovirus_1_F1429C1F1578
Grenada_mosquito_rhabdovirus_1_F1429C1F1431
Grenada_mosquito_rhabdovirus_1_F1429C1F1576
Grenada_mosquito_rhabdovirus_1_F1429C1F1570
Grenada_mosquito_rhabdovirus_1_F1429C1F1582
Grenada_mosquito_rhabdovirus_1_F1429C1F1427
Grenada_mosquito_rhabdovirus_1_F1429C1F1430
Grenada_mosquito_rhabdovirus_1_F1429C1F1571
Grenada_mosquito_rhabdovirus_1_F1429C1F1419
Grenada_mosquito_rhabdovirus_1_F1429C1F1429
KY435948_Cururu_virus
Grenada_mosquito_rhabdovirus_1_F1429C1F1115
Grenada_mosquito_rhabdovirus_1_F1429C1F1124
Grenada_mosquito_rhabdovirus_1_F1429C1F1433
Grenada_mosquito_rhabdovirus_1_F1429C1F1119
Grenada_mosquito_rhabdovirus_1_F1429C1F1121
Grenada_mosquito_rhabdovirus_1_F1429C1F1118
Grenada_mosquito_rhabdovirus_1_F1429C1F1113
Grenada_mosquito_rhabdovirus_1_F1429C1F1114
Grenada_mosquito_rhabdovirus_1_F1429C1F1112
Grenada_mosquito_rhabdovirus_1_F1429C1F1120
Grenada_mosquito_rhabdovirus_1_F1429C1F1123
Grenada_mosquito_rhabdovirus_1_F1429C1F1125
Grenada_mosquito_rhabdovirus_1_F1429C1F1126
MN661109_Atrato_Rhabdo-like_virus_2
Ardicrand_virus_F1568C1F1277
Ardicrand_virus_F1568C1F1568
Antlested_virus_F1034C1F1034
MW452302_Wuhan_Mosquito_Virus_9
NC_031303_Wuhan_Mosquito_Virus_9
OL700091_Wuhan_Mosquito_Virus_9
MW434774_Stang_virus
MW434775_Stang_virus
MZ202303_Cimo_rhabdovirus_III
MW434767_Elisy_virus
MW434768_Elisy_virus
OL700136_XiangYun_mono-chu-like_virus_11
"
Monjiviricetes_clade3 <- str_split(Monjiviricetes_clade3_str, "\\s")[[1]]
Monjiviricetes_clade3 <- Monjiviricetes_clade3[which(!Monjiviricetes_clade3 == "")]

Monjiviricetes_clade_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/Monjiviricetes_clade_aa.fasta",
                                                  seqtype = "AA")
Monjiviricetes_clade3_aa_seq <- Monjiviricetes_clade_aa_seq[Monjiviricetes_clade3]
names(Monjiviricetes_clade3_aa_seq)
seqinr::write.fasta(sequences = Monjiviricetes_clade3_aa_seq , names =names(Monjiviricetes_clade3_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/Monjiviricetes_clade3_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "MZ771217_Rhabdoviridae_sp_"
new_outg_seq <- Monjiviricetes_clade_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))



####### Zooming in on different clades
#### Monjiviricetes_clade5 tree
Monjiviricetes_clade5_str <- "Preculath_virus_F1097C1F1526
Preculath_virus_F1097C1F1537
Preculath_virus_F1097C1F1533
Preculath_virus_F1097C1F1527
Preculath_virus_F1097C1F1525
Preculath_virus_F1097C1F1532
Preculath_virus_F1097C1F1536
Preculath_virus_F1097C1F1531
Preculath_virus_F1097C1F1098
Preculath_virus_F1097C1F1099
Preculath_virus_F1097C1F1021
Preculath_virus_F1097C1F1089
Preculath_virus_F1097C1F1019
Preculath_virus_F1097C1F1097
Preculath_virus_F1097C1F1090
Preculath_virus_F1097C1F1534
Preculath_virus_F1097C1F1371
Preculath_virus_F1097C1F1380
Reloonsia_virus_F1535C4F1020
Reloonsia_virus_F1535C4F1535
Reloonsia_virus_F1535C4F1524
BK059423_San_Gabriel_mononegavirus
San_Gabriel_mononegavirus_F1076C1F1350
San_Gabriel_mononegavirus_F1076C1F1168
San_Gabriel_mononegavirus_F1076C1F1518
San_Gabriel_mononegavirus_F1076C1F1076
San_Gabriel_mononegavirus_F1076C1F1016
San_Gabriel_mononegavirus_F1076C1F1197
San_Gabriel_mononegavirus_F1076C1F1345
Retriended_virus_F1011C1F1011
MN567480_Primus_virus
OL700129_XiangYun_mono-chu-like_virus_4
MN661034_Atrato_Rhabdo-like_virus_3
ON955251_Joutseno_rhabdovirus_1
"
Monjiviricetes_clade5 <- str_split(Monjiviricetes_clade5_str, "\\s")[[1]]
Monjiviricetes_clade5 <- Monjiviricetes_clade5[which(!Monjiviricetes_clade5 == "")]

Monjiviricetes_clade_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/Monjiviricetes_clade_aa.fasta",
                                                  seqtype = "AA")
Monjiviricetes_clade5_aa_seq <- Monjiviricetes_clade_aa_seq[Monjiviricetes_clade5]
names(Monjiviricetes_clade5_aa_seq)
seqinr::write.fasta(sequences = Monjiviricetes_clade5_aa_seq , names =names(Monjiviricetes_clade5_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/Monjiviricetes_clade5_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "NC_031276_Wuhan_Ant_Virus"
new_outg_seq <- Monjiviricetes_clade_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))



####### Zooming in on different clades
#### Monjiviricetes_clade5 tree
Monjiviricetes_clade5_str <- "MT240825_Almendravirus_menghai
MT240827_Almendravirus_menghai
MT240824_Almendravirus_menghai
MT240828_Almendravirus_menghai
NC_040602_Almendravirus_menghai
LC270812_Almendravirus_menghai
OP264923_Almendravirus_menghai
Poterfers_virus_F1543C2F1543
NC_031957_Almendravirus_cootbay
MW890016_Shanxi_Armigeres_subalbatus_rhabdovirus
OK491499_Xiangshan_rhabdo-like_virus_1
NC_031958_Almendravirus_chico
NC_039200_Almendravirus_balsa
MW890015_Shanxi_Arboretum_virus
OK491516_Almendravirus_arboretum
NC_025393_Almendravirus_arboretum
NC_025395_Almendravirus_almendras
"
Monjiviricetes_clade5 <- str_split(Monjiviricetes_clade5_str, "\\s")[[1]]
Monjiviricetes_clade5 <- Monjiviricetes_clade5[which(!Monjiviricetes_clade5 == "")]

Monjiviricetes_clade_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/Monjiviricetes_clade_aa.fasta",
                                                  seqtype = "AA")
Monjiviricetes_clade5_aa_seq <- Monjiviricetes_clade_aa_seq[Monjiviricetes_clade5]
names(Monjiviricetes_clade5_aa_seq)
seqinr::write.fasta(sequences = Monjiviricetes_clade5_aa_seq , names =names(Monjiviricetes_clade5_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/Monjiviricetes_clade5_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "OP589941_Rhabdoviridae_sp_"
new_outg_seq <- Monjiviricetes_clade_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))







####### Zooming in on different clades
#### Monjiviricetes_clade6 tree
Monjiviricetes_clade6_str <- "MW434771_Merhavirus_merida
MW434772_Merhavirus_merida
MW434770_Merhavirus_merida
NC_040599_Merhavirus_merida
MW434769_Merhavirus_merida
MW434773_Merhavirus_merida
OM817546_Merhavirus_merida
MH188000_Culex_rhabdovirus
MH310083_Merhavirus_merida
OL700074_Merhavirus_merida
MT577803_Culex_rhabdovirus
MF882997_Merida-like_virus_Turkey
NC_040532_Merida-like_virus_KE-2017a
MK440621_Merhavirus_merida
Begrized_virus_F1588C1F1462
Begrized_virus_F1588C1F1588
Ausist_virus_F1341C1F1079
Ausist_virus_F1341C1F1344
Ausist_virus_F1341C1F1167
Ausist_virus_F1341C1F1337
Ausist_virus_F1341C1F1341
Ausist_virus_F1341C1F1003
NC_025384_Merhavirus_tritaeniorhynchus
OL700055_Merhavirus_tritaeniorhynchus
MW452298_Merhavirus_tritaeniorhynchus
LC514403_Merhavirus_tritaeniorhynchus
LC026102_Merhavirus_tritaeniorhynchus
Arlicasen_virus_F1276C1F1256
Arlicasen_virus_F1276C1F1276
ON955141_Enontekio_merhavirus
KF310911_Beaumont_virus
MW520372_Evros_rhabdovirus_2
Arindly_virus_F1098C3F1095
Arindly_virus_F1098C3F1219
Arindly_virus_F1098C3F1090
Arindly_virus_F1098C3F1099
Arindly_virus_F1098C3F1087
Arindly_virus_F1098C3F1098
BK059424_Formosus_virus
ON955142_Hattula_rhabdovirus
ON955246_Hattula_rhabdovirus
ON955143_Inari_rhabdovirus
MF176269_Ohlsrhavirus_culex
NC_035132_Ohlsrhavirus_culex
MF176358_Ohlsrhavirus_culex
NC_055292_Ohlsrhavirus_northcreek
NC_028484_Ohlsrhavirus_tongilchon
LC514056_Ohlsrhavirus_pseudovishnui
LC514057_Ohlsrhavirus_pseudovishnui
OL700054_Ohlsrhavirus_pseudovishnui
MH188003_Ohlsrhavirus_culex
KU248086_Riverside_virus_1
KU248087_Riverside_virus_1
NC_040669_Riverside_virus_1
Regreagly_virus_F1187C3F1177
Regreagly_virus_F1187C3F1178
Regreagly_virus_F1187C3F1187
KY768857_Ohlsrhavirus_ohlsdorf
NC_055477_Ohlsrhavirus_ohlsdorf
KY768859_Ohlsrhavirus_ohlsdorf
KY768860_Ohlsrhavirus_ohlsdorf
KY768861_Ohlsrhavirus_ohlsdorf
ON955253_Ohlsrhavirus_ohlsdorf
ON955254_Ohlsrhavirus_ohlsdorf
MF344595_Ohlsrhavirus_lobeira
MW826481_Rhabdoviridae_sp_
MZ209762_Hangzhou_tipula_scripta_rhabdovirus_1
MZ771226_Rhabdoviridae_sp_
ON059778_Lantra_virus
ON059779_Lantra_virus
"
Monjiviricetes_clade6 <- str_split(Monjiviricetes_clade6_str, "\\s")[[1]]
Monjiviricetes_clade6 <- Monjiviricetes_clade6[which(!Monjiviricetes_clade6 == "")]

Monjiviricetes_clade_aa_seq <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/Monjiviricetes_clade_aa.fasta",
                                                  seqtype = "AA")
Monjiviricetes_clade6_aa_seq <- Monjiviricetes_clade_aa_seq[Monjiviricetes_clade6]
names(Monjiviricetes_clade6_aa_seq)
seqinr::write.fasta(sequences = Monjiviricetes_clade6_aa_seq , names =names(Monjiviricetes_clade6_aa_seq), 
                    file.out = "/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/Monjiviricetes_clade6_aa.fasta")
##keep the closest sequence as an outgroup
new_outg <- "OK491503_Xiangshan_rhabdo-like_virus_5"
new_outg_seq <- Monjiviricetes_clade_aa_seq[new_outg]
seqinr::write.fasta(sequences = new_outg_seq , names =names(new_outg_seq), 
                    file.out = paste0("/full_path_to/wd/RdRp_scan/analysis/phylo/Monjiviricetes_clade/aln/aa_seq/", 
                                      new_outg, ".fasta"))





