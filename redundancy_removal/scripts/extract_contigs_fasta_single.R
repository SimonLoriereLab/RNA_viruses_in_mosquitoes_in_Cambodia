library(tidyverse)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)
commandArgs()

s <- args[1]

metadata_dir <- "/full_path_to/wd/redundancy_removal/metadata"
output_dir <- "/full_path_to/wd/redundancy_removal/analysis/extracted_renamed_contigs"
dir.create(output_dir)

  
dir_to_sample1 <- "/full_path_to/wd/auto_metaviromics/analyzed_samples"
dir_to_sample2 <- "/metaSPAdes/metaSPAdes_output_"

contig_path <- paste0(dir_to_sample1, "/", s, dir_to_sample2, s, "/contigs.fasta")
contig_fasta <- read.fasta(file = contig_path, seqtype = "DNA")
new_headers <- paste0(s, "_", names(contig_fasta))
print(new_headers[1:3])
write.fasta(sequences = contig_fasta, names = new_headers, file.out = paste0(output_dir, "/", s, "_renamed.fasta"))
