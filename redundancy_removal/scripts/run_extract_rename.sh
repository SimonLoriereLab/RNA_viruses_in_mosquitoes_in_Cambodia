#!/bin/bash
#SBATCH -J extract_rename
#SBATCH --mem=1000
#SBATCH -c 1
#SBATCH --qos=normal
#SBATCH -p common

_list=`sed -n ${SLURM_ARRAY_TASK_ID}p $1`

#Load the necessary modules
module load R/4.1.0

echo "$_list"

Rscript /full_path_to/wd/redundancy_removal/scripts/extract_contigs_fasta_single.R $_list

echo "DONE"
