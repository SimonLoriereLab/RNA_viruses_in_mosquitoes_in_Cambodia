#!/bin/bash
#SBATCH -J DMND_p_LCA
#SBATCH --mem=500000
#SBATCH -c 95
#SBATCH --qos=geva
#SBATCH -p geva

directory="/full_path_to/wd/RdRp_scan/analysis"
input="/full_path_to/wd/RdRp_scan/analysis/all_samples_renamed_contigs_redundancy_removed_code1_aa_nr_rdrp.fasta"
output="/full_path_to/wd/RdRp_scan/analysis/DIAMOND_blastp_hits_nrcontigs_rdrp_LCA.txt"
diamond_db="/full_path_to/db/nrprot_tax_20220814.dmnd"
lineage="/full_path_to/db/ncbi_lineages_2022-08-23.csv"
lineage_annotation_script="/full_path_to/wd/RdRp_scan/scripts/DIAMOND_LCA_tax2lin.R"


module load diamond/2.0.4 R/4.1.0

cd $directory

# diamond blastp -d $diamond_db -q $input --ultra-sensitive --top 5 \
# -b 2.0 -c 1.0 --threads ${SLURM_CPUS_PER_TASK} -f 102 -o $output

Rscript --vanilla $lineage_annotation_script $output $lineage $directory


echo "DONE"
