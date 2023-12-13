#!/bin/bash

#SBATCH -J wrap_kraken2
#SBATCH -p geva
#SBATCH -q geva
#SBATCH -c 1
#SBATCH --mem=8000

###############################################################################
# This script is written in case trimmed data are in auto_metaviromics
# in separate folders for each sample.
###############################################################################
#To be changed:
_directory="/full_path_to/wd"

_DB_path="/full_path_to/nt_2021-10-25/kraken/2.1.1/nt/"

_sample_list="QC/metadata/analyzed_sample_list_by_20211021.tsv" #### prepared already

#### or make one (sample name before _R1.fastq.gz or _R2.fastq.gz )
###############################################################################

_krakentwo_jobserrors="$_directory/QC/analysis/Kraken2/Kraken2_jobserrors"
_krakentwo_jobsoutputs="$_directory/QC/analysis/Kraken2/Kraken2_jobsoutputs"

cd $_directory

if [ ! -d "$_krakentwo_jobserrors" ]; then
mkdir -pv $_krakentwo_jobserrors
fi

if [ ! -d "$_krakentwo_jobsoutputs" ]; then
mkdir -pv $_krakentwo_jobsoutputs
fi

_nb_jobs=`wc -l < $_sample_list` #Computes number of files to process
echo $_nb_jobs

sbatch --wait --array=1-$_nb_jobs -o $_krakentwo_jobsoutputs/slurm-%A_%a.out -e \
$_krakentwo_jobserrors/slurm-%A_%a.err \
$_directory/QC/scripts/Kraken2/script-Kraken2_nt.sh \
$_directory $_sample_list $_DB_path || exit 1

exit 0
