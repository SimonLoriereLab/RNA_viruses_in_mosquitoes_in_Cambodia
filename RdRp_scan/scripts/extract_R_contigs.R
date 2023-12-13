library(tidyverse)
library(seqinr)

R_contig_set <- "/pasteur/zeus/projets/p02/GEVA/users/Artem/Mosquito_RNA_viruses_Cambodia/RdRp_scan/analysis/R_df_hmmscan.LCA.BestEvalue_hit.txt"

df <- read_delim(R_contig_set, delim = "\t")

contig_headers <- sort(unique(df$contig_id))

all_contigs_file <- "/pasteur/zeus/projets/p02/GEVA/users/Artem/Mosquito_RNA_viruses_Cambodia/redundancy_removal/analysis/extracted_renamed_contigs/all_samples_renamed_contigs_redundancy_removed.fasta"
all_contigs <- read.fasta(file = all_contigs_file, seqtype = "DNA")
contigs <- all_contigs[contig_headers]

write.fasta(sequences = contigs, names = names(contigs), 
            file.out = "/pasteur/zeus/projets/p02/GEVA/users/Artem/Mosquito_RNA_viruses_Cambodia/RdRp_scan/analysis/R_contigs.fasta")