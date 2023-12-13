library(tidyverse)
library(seqinr)

cISF_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/cISF_novel_sequence_data.txt", delim = "\t")

### select consensus sequences based on mapping coverage
### first round of read mapping results
read_map_df <-  read_delim("/full_path_to/wd/RdRp_scan/analysis/read_mapping/mapping_summary_bowtie2_RdRpScan_Round1_analysis/results_mapped_only.txt", 
           delim = "\t")


read_map_df_cISF <- read_map_df[which(read_map_df$ref_id %in% cISF_df$contig_id & read_map_df$selected == "yes"),]



read_map_df_cISF$RdRp_scan_map_consensus_id <- paste0(read_map_df_cISF$ref_id, "_",
                      read_map_df_cISF$sample_id, "_RdRpScan_Round1")

read_map_df_cISF$contig_id <- read_map_df_cISF$ref_id

consensus_df <- left_join(read_map_df_cISF, cISF_df, by  = "contig_id")



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

consensus_df$vir_sp_seq_ID_v3 <- paste0(consensus_df$phylo_virus_species, "_",
                                         consensus_df$original_sample_ID, 
                                         consensus_df$original_metaspades_contig_node_IDs,
                                         consensus_df$read_map_consensus_sample_ID)

dir_cons1 <- "/full_path_to/wd/RdRp_scan/analysis/read_mapping/cons1_renamed/"
extracted_consensus_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/extracted_cons1_renamed/"
dir.create(extracted_consensus_directory)
extracted_orf_nt_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/extracted_cons1_orf_nt/"
dir.create(extracted_orf_nt_directory)
extracted_orf_aa_directory <- "/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/extracted_cons1_orf_aa/"
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

write.table(x = consensus_df, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/consensus_df.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)

#### Manually correct in Geneious:
### the chimeric contig - Accessaty virus
### Cell_fusing_agent_virus_F1335C10F1065 (extra N at the beginning)







#### REFERENCE SEQUENCES

### Use DIAMOND blastp results to get the hits
blastp_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/lineage_collapsed.DIAMOND_blastp_hits_nrcontigs_rdrp.txt.gz", 
                        delim = "\t") %>% 
  select(orf_id = qseqid, sseqid, pident, length, stitle, tax_id, superkingdom)
blastp_df1 <- blastp_df[which(blastp_df$superkingdom == "Viruses"),]

blastp_df2 <- blastp_df1[which(blastp_df1$orf_id %in% sort(unique(consensus_df$orf_id))),]

##### Sort by release date on NCBI Viruses and select only submitted later than 2020-01-30
paste0(sort(unique(blastp_df2$sseqid)), collapse = ", ")

##### Reimport metadata table
new_cISF_df <- read_csv("/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/metadata_new_public_cISF.csv")


#### Let's extract sequences from the previously curated cISF alignment and update metadata

ns_aln_previous <- seqinr::read.fasta(file = "/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/Flaviviridae/phylo/aa_tree/cISF_full_polyprot_aa_aln_v2.fasta.gap_clean.dedup.derec.outg_nsORF_669_3266_correct.fasta",
                   seqtype = "AA")


names(ns_aln_previous)


prev_tree_metadata <- read_delim("/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/Flaviviridae/phylo/aa_tree/cISF_full_polyprot_aa_aln_v2.fasta.gap_clean.dedup.derec.outg_sORF_1_668_correct.fasta.contree.tip.metadata.full_novelseq_renamed_v2.tsv",delim = "\t")
prev_tree_metadata_b <- read_delim("/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/Flaviviridae/phylo/aa_tree/cISF_full_polyprot_aa_aln_v2.fasta.gap_clean.dedup.derec.outg_sORF_1_668_correct.fasta.contree.tip.metadata.brief_novelseq_renamed_v2.tsv",delim = "\t")


#### substitute "Accrood_virus_F1616C2F1616" back to "Palm_Creek_virus_F1616C2F1616"
write_file(x = str_replace_all(read_file("/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/Flaviviridae/phylo/aa_tree/cISF_full_polyprot_aa_aln_v2.fasta.gap_clean.dedup.derec.outg_sORF_1_668_correct.fasta.contree.tip.metadata.brief_novelseq_renamed_v2.tsv"),
                               "Accrood_virus_F1616C2F1616", "Palm_Creek_virus_F1616C2F1616"), 
           file = "/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/metadata_previous_analysis_RdRp_search.txt", append = F)

write_file(x = str_replace_all(read_file("/full_path_to/wd/RdRp_search/analysis/DIAMOND_family_level/Flaviviridae/phylo/aa_tree/cISF_full_polyprot_aa_aln_v2.fasta.gap_clean.dedup.derec.outg_sORF_1_668_correct.fasta.contree.tip.metadata.full_novelseq_renamed_v2.tsv"),
                               "Accrood_virus_F1616C2F1616", "Palm_Creek_virus_F1616C2F1616"), 
           file = "/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/metadata_previous_analysis_RdRp_search_full.txt", append = F)


prev_tree_metadata <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/metadata_previous_analysis_RdRp_search_full.txt",
                                 delim = "\t")
prev_tree_metadata_b <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/metadata_previous_analysis_RdRp_search.txt",
                                   delim = "\t")
### Check if names are matching 
prev_mdf_novel <- prev_tree_metadata[which(prev_tree_metadata$novel_seq == 1),]
setdiff(prev_mdf_novel$tip_join_var, consensus_df$vir_sp_seq_ID_v3)

#### check if some of the new references are already in the dataset or have been discarded from it for some reason
not_in_the_dataset <- setdiff(new_cISF_df$Accession, prev_tree_metadata$tip_join_var)
the_newest_cISF <- new_cISF_df[which(new_cISF_df$Accession %in% not_in_the_dataset),]

#### Update cISF phylogeny metadata
cat(paste0(colnames(prev_tree_metadata_b), collapse = ", "))
selected_mdf <- prev_tree_metadata %>% select(tip_join_var, nt_seq_ID, tip_labels, host_taxon, host_genus, Country, Location, Year, novel_seq, nt_slen, nt_title)
# 
# colnames(prev_tree_metadata)




the_newest_cISF$nt_seq_ID <-  paste0(the_newest_cISF$Accession, "_", str_replace_all(str_extract(the_newest_cISF$GenBank_Title, "[[:print:]]+(?= strain)"), " ", "_"))
the_newest_cISF$Year <- str_extract(the_newest_cISF$Collection_Date, "[[:digit:]]{4}")
the_newest_cISF1 <- the_newest_cISF %>% select(tip_join_var = Accession, 
                           nt_seq_ID,
                           host_taxon = Host, 
                           Country = Country, Location = Country,
                           Year, nt_slen = Length, nt_title = GenBank_Title)
the_newest_cISF1$host_taxon <- "Aedes"
the_newest_cISF1[which(the_newest_cISF1$tip_join_var  == "MZ202252"), ]$host_taxon <- "Culex"
the_newest_cISF1$host_genus <-  the_newest_cISF1$host_taxon



selected_mdf$Year <- as.character(selected_mdf$Year)

selected_mdf_upd <- bind_rows(selected_mdf, the_newest_cISF1)

write.table(x = selected_mdf_upd, file = "/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/updated_metadata.txt", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)




######################## Dealing with CDS
####### Extract CDS from NCBI nucleotide, not NCBI Virus, because the headers are more informative
paste0(sort(selected_mdf_upd[which(is.na(selected_mdf_upd$novel_seq)),]$tip_join_var), collapse = ", ")
### using this accession list to NCBI and extract CDS sequences and import them here.
#### Put the file here!
cds_fasta_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/cISF_refs_upd_cds_nt.txt"


### There might be more automatic solutions for CDS extractions with EMBOSS or CL implementation of entrez. But web NCBI is the fastest for now.

############################################################################################
############################################################################################
############################################################################################
cds_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/ref_cds_nt"
dir.create(cds_clade_dir)
aa_clade_dir <- "/full_path_to/wd/RdRp_scan/analysis/phylo/cISF_clade/ref_cds_aa"
dir.create(aa_clade_dir)


cds_fasta <- seqinr::read.fasta(file = cds_fasta_file, seqtype = "DNA", whole.header = T)
names(cds_fasta)
length(selected_mdf_upd[which(is.na(selected_mdf_upd$novel_seq)),]$tip_join_var)


cds_names_df <- tibble(original_cds_headers = names(cds_fasta),
                       nt_accession = str_remove(str_extract(names(cds_fasta), "(?<=^lcl\\|)[[:print:]]+(?=_cds)"), "\\.[[:digit:]]*"),
                       cds_product = str_replace_all(str_extract(names(cds_fasta), "(?<=\\[protein\\=)[[:print:]]+?(?=\\])"), " ", "_")
)

setdiff(sort(cds_names_df$nt_accession), ref_df$tip_join_var)

ref_df <- selected_mdf_upd[which(is.na(selected_mdf_upd$novel_seq)),]

cds_names_df1 <- left_join(cds_names_df, ref_df,  by = c( "nt_accession"= "tip_join_var"))
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







