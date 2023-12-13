#!/bin/bash

#SBATCH -J wrap_k2_extract
#SBATCH -p geva
#SBATCH -q geva
#SBATCH -c 1
#SBATCH --mem=5000

#To be changed:
_directory="/full_path_to/wd"

_script_path="/full_path_to/wd/QC/scripts/KrakenTools/extract_kraken_reads.py"

_sample_list="QC/metadata/sample_list1.tsv" #### prepared already


_k2extract_jobserrors="/full_path_to_scratch/temp_Kraken2_errouts/k2extract_jobserrors"
_k2extract_jobsoutputs="/full_path_to_scratch/temp_Kraken2_errouts/k2extract_jobsoutputs"

cd $_directory

if [ ! -d "$_k2extract_jobserrors" ]; then
mkdir -pv $_k2extract_jobserrors
fi

if [ ! -d "$_k2extract_jobsoutputs" ]; then
mkdir -pv $_k2extract_jobsoutputs
fi

_nb_jobs=`wc -l < $_sample_list` #Computes number of files to process

#_nb_jobs=1 #Computes number of files to process

echo $_nb_jobs

sbatch --wait --array=1-$_nb_jobs -o $_k2extract_jobsoutputs/slurm-%A_%a.out -e \
$_k2extract_jobserrors/slurm-%A_%a.err \
$_directory/QC/scripts/Kraken2/script-Kraken2_extract_reads.sh \
$_directory $_sample_list $_script_path || exit 1

exit 0
