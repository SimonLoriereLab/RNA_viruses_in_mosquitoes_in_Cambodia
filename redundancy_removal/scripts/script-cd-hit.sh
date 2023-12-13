#!/bin/bash
#SBATCH -J redund_remove
#SBATCH --mem=500000
#SBATCH -c 95
#SBATCH -q geva
#SBATCH -p geva

_directory="/full_path_to/wd/redundancy_removal/analysis/extracted_renamed_contigs"

cd $_directory

# cat $_directory/*_renamed.fasta > $_directory/all_samples_renamed_contigs.fasta

#Load the necessary modules
module load blast+/2.10.0 cd-hit/4.8.1

cd-hit-est -i $_directory/all_samples_renamed_contigs.fasta \
-o $_directory/all_samples_renamed_contigs_redundancy_removed.fasta -c 0.9 -n 9 -sf 1 -T 95 -M 500000 -d 0




echo "DONE"
