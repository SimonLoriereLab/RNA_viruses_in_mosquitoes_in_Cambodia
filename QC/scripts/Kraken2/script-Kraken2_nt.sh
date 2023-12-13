#!/bin/bash
#SBATCH -J kraken2
#SBATCH -p geva
#SBATCH -q geva
#SBATCH -c 95
#SBATCH --mem=50000

module load blast+/2.10.1
module load kraken/2.1.1

_directory=$1
cd $_directory
_sample_id=`sed -n ${SLURM_ARRAY_TASK_ID}p $2`
_DB_path=$3
_output_folder="${_directory}/auto_metaviromics/analyzed_samples/${_sample_id}/Kraken2"
###############################################################################

if [ ! -d "$_output_folder" ]; then
mkdir -pv $_output_folder
fi


kraken2 --db $_DB_path --gzip-compressed --threads 10 --memory-mapping --use-names --report-zero-counts --report-minimizer-data \
--output "${_output_folder}/${_sample_id}_ntDB_kraken2_output.txt" \
--report "${_output_folder}/${_sample_id}_ntDB_kraken2_report.txt" \
--paired \
"${_directory}/auto_metaviromics/analyzed_samples/${_sample_id}/Trimmed/${_sample_id}_R1.qual.fastq.gz" \
"${_directory}/auto_metaviromics/analyzed_samples/${_sample_id}/Trimmed/${_sample_id}_R2.qual.fastq.gz"
