#!/bin/bash

#SBATCH -J wrap_k2_extract
#SBATCH -p geva
#SBATCH -q geva
#SBATCH -c 1
#SBATCH --mem=5000

#To be changed:
_directory="/full_path_to/wd"
_output_folder="/full_path_to/wd/mosquito_species_verification/analysis/SKA/split_kmer_files"
_sample_list="/full_path_to/wd/mosquito_species_verification/metadata/sample_list1.tsv" #### prepared already


_SKA_split_kmers_jobserrors="/full_path_to_scratch/temp_SKA_errouts/SKA_split_kmers_jobserrors"
_SKA_split_kmers_jobsoutputs="/full_path_to_scratch/temp_SKA_errouts/SKA_split_kmers_jobsoutputs"

cd $_directory

if [ ! -d "$_SKA_split_kmers_jobserrors" ]; then
mkdir -pv $_SKA_split_kmers_jobserrors
fi

if [ ! -d "$_SKA_split_kmers_jobsoutputs" ]; then
mkdir -pv $_SKA_split_kmers_jobsoutputs
fi

if [ ! -d "$_output_folder" ]; then
mkdir -pv $_output_folder
fi

_nb_jobs=`wc -l < $_sample_list` #Computes number of files to process

echo $_nb_jobs

sbatch --wait --array=1-$_nb_jobs -o $_SKA_split_kmers_jobsoutputs/slurm-%A_%a.out -e \
$_SKA_split_kmers_jobserrors/slurm-%A_%a.err \
$_directory/mosquito_species_verification/scripts/script-SKA_split_kmers.sh \
$_directory $_sample_list $_output_folder || exit 1

exit 0
