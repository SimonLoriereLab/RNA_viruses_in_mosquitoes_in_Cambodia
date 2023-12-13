#!/bin/bash
#SBATCH -J k2_extract
#SBATCH -p common
#SBATCH -q normal
#SBATCH -c 1
#SBATCH --mem=5000

_directory=$1
cd $_directory
_sample_id=`sed -n ${SLURM_ARRAY_TASK_ID}p $2`
_script_path=$3
_output_folder="${_directory}/auto_metaviromics/analyzed_samples/${_sample_id}/Kraken2"
###############################################################################

module load Python/3.8.3

### Culicidae taxonomy ID is 7157

$_script_path -k "${_output_folder}/${_sample_id}_ntDB_kraken2_output.txt" \
-s1 "${_directory}/auto_metaviromics/analyzed_samples/${_sample_id}/Trimmed/${_sample_id}_R1.qual.fastq.gz" \
-s2 "${_directory}/auto_metaviromics/analyzed_samples/${_sample_id}/Trimmed/${_sample_id}_R2.qual.fastq.gz" \
-o "${_output_folder}/${_sample_id}_R1.Culicidae.qual.fastq" \
-o2 "${_output_folder}/${_sample_id}_R2.Culicidae.qual.fastq" \
--fastq-output \
--taxid 7157 \
--report "${_output_folder}/${_sample_id}_ntDB_kraken2_report.txt" \
--include-children
