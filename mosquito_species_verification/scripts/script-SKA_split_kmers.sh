#!/bin/bash
#SBATCH -J k2_extract
#SBATCH -p common
#SBATCH -q normal
#SBATCH -c 1
#SBATCH --mem=5000

_directory=$1
_sample_id=`sed -n ${SLURM_ARRAY_TASK_ID}p $2`
_output_folder=$3
_input_folder="${_directory}/auto_metaviromics/analyzed_samples/${_sample_id}/Kraken2"
###############################################################################

cd $_output_folder

ska fastq -o ${_sample_id} ${_input_folder}/${_sample_id}_R1.Culicidae.qual.fastq ${_input_folder}/${_sample_id}_R2.Culicidae.qual.fastq
