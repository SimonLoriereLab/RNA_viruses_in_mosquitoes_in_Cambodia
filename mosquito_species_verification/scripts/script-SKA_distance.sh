#!/bin/bash

#SBATCH -J SKA_dist
#SBATCH -p geva
#SBATCH -q geva
#SBATCH -c 95
#SBATCH --mem=500000

#To be changed:
_directory="/full_path_to/wd"
_output_folder="/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances"
_sk_filenames_file="/full_path_to/wd/mosquito_species_verification/metadata/split_kmer_filenames_1.txt"

if [ ! -d "$_output_folder" ]; then
mkdir -pv $_output_folder
fi

cd $_output_folder
ska distance -i 0.2 -s 1000  -f $_sk_filenames_file -o Culicidae_transcriptome_comparison_all
