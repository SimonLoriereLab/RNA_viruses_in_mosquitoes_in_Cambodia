#!/bin/bash
#SBATCH -J DMND_x
#SBATCH --mem=500000
#SBATCH -c 95
#SBATCH --qos=geva
#SBATCH -p geva

directory="/full_path_to/wd/RdRp_search/analysis/DIAMOND"
input="/full_path_to/wd/redundancy_removal/analysis/extracted_renamed_contigs/all_samples_renamed_contigs_redundancy_removed.fasta"
output="/full_path_to/wd/RdRp_search/analysis/DIAMOND/DIAMOND_hits_nrcontigs.txt"
diamond_db="/full_path_to/db/nrprot_tax_20210717.dmnd"
lineage="/full_path_to/db/ncbi_lineages_2021-07-23.csv"
lineage_annotation_script="/full_path_to/wd/RdRp_search/scripts/Lineage_v1.2.3.R"

module load diamond/2.0.6 R/4.1.0

cd $directory

diamond blastx -d $diamond_db -q $input --sensitive --top 5 \
-b 2.0 -c 1.0 --threads 95 \
-f 6 qseqid sseqid pident length mismatch gapopen \
qstart qend sstart send slen evalue bitscore stitle staxids -o $output

Rscript --vanilla $lineage_annotation_script $output $lineage $directory

echo "DONE"
