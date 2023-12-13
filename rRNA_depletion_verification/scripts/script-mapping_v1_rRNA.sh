#!/bin/bash
#SBATCH -J mapAB

#SBATCH -p clcbio
#SBATCH --qos=clcbio
#SBATCH -c 4
#SBATCH --mem=10000


###############################################################################
_directory=$1
cd $_directory


_sample_id=$2
_ref_id=`sed -n ${SLURM_ARRAY_TASK_ID}p $3`
_reference_folder=$4
_summary_folder=$5
_scripts=$6
_sample_input_folder=$7
_sample_output_folder=$8

###############################################################################

###############################################################################
##### Careful here, if the name of the pipeline changes
##### ( e. g. not auto_metaviromics-master but auto_metaviromics-some_other_version)
_sample_input_path="${_directory}/auto_metaviromics/analyzed_samples/${_sample_id}/${_sample_input_folder}"
_sample_output_path="${_directory}/auto_metaviromics/analyzed_samples/${_sample_id}/${_sample_output_folder}"

if [ ! -d "$_sample_output_path" ]; then
mkdir $_sample_output_path
fi

_reference_path=${_reference_folder}
_consensus_path=${_summary_folder}
###############################################################################

###############################################################################
cd $_sample_output_path

module load clc-assembly-cell/5.1.0

srun -p clcbio --qos=clcbio -c 4 --mem=10000 \
clc_mapper -o $_sample_output_path/${_ref_id}_${_sample_id}.cas \
-a local -t 1 -r random -l 0.9 -s 0.9 -x 2 -g 3 -e 3 \
--noprogress \
-d ${_reference_path}/${_ref_id}.fasta \
-q -p fb ss 100 1500 \
-i ${_sample_input_path}/${_sample_id}_R1.qual.fastq.gz ${_sample_input_path}/${_sample_id}_R2.qual.fastq.gz
#-i ${_sample_input_path}/${_sample_id}_R1.fastq.qual.gz ${_sample_input_path}/${_sample_id}_R2.fastq.qual.gz

echo "\n\n\n--------- alignment done ---------\n\n\n"

srun -p clcbio --qos=clcbio -c 4 --mem=10000 \
clc_mapping_info -c -p  fb  ss  100  1500  -e 10 \
$_sample_output_path/${_ref_id}_${_sample_id}.cas > $_sample_output_path/${_ref_id}_${_sample_id}_mapinfo.txt

echo "\n\n\n--------- info extracted ---------\n\n\n"

srun -p clcbio --qos=clcbio -c 4 --mem=10000 \
clc_cas_to_sam -a $_sample_output_path/${_ref_id}_${_sample_id}.cas \
-o $_sample_output_path/${_ref_id}_${_sample_id}.bam -f 33 -u

echo "\n\n\n--------- converted to bam ---------\n\n\n"

module load samtools/1.10
srun -p common --qos=normal -c 4 --mem=8000 samtools sort \
-o $_sample_output_path/${_ref_id}_${_sample_id}.sorted.bam \
$_sample_output_path/${_ref_id}_${_sample_id}.bam \
-@ 4 -m 4G

echo "\n\n\n--------- bam sorted ---------\n\n\n"

module load bedtools/2.29.2
srun -p common --qos=normal -c 2 --mem=4000 bedtools genomecov \
-d -ibam $_sample_output_path/${_ref_id}_${_sample_id}.sorted.bam > $_sample_output_path/${_ref_id}_${_sample_id}.cov

echo "\n\n\n--------- coverage calculated, cov created ---------\n\n\n"

module load samtools/1.10
srun -p common --qos=normal -c 2 --mem=4000 \
samtools faidx ${_reference_path}/${_ref_id}.fasta

echo "\n\n\n--------- faidx done ---------\n\n\n"

module unload samtools/1.10
module load samtools/1.10
module load lofreq/2.1.4

srun -p common --qos=normal -c 4 --mem=8000 \
lofreq call --force-overwrite -f ${_reference_path}/${_ref_id}.fasta \
-o $_sample_output_path/${_ref_id}_${_sample_id}.vcf \
$_sample_output_path/${_ref_id}_${_sample_id}.sorted.bam

echo "\n\n\n--------- varcall done, vcf created ---------\n\n\n"

cd ${_summary_folder}
Rscript $_scripts/coverage_norm.R ${_directory} ${_sample_id} ${_ref_id} \
$_sample_output_path/${_ref_id}_${_sample_id}_mapinfo.txt \
$_sample_output_path/${_ref_id}_${_sample_id}.cov \
$_sample_output_path/${_ref_id}_${_sample_id}.vcf \
${_summary_folder}


##### This one is if lofreq fails or freezes. Manually rerun again and mute commands above.
# cd ${_summary_folder}
# Rscript $_scripts/coverage_norm_no_snv.R ${_directory} ${_sample_id} ${_ref_id} \
# $_sample_output_path/${_ref_id}_${_sample_id}_mapinfo.txt \
# $_sample_output_path/${_ref_id}_${_sample_id}.cov \
# ${_summary_folder}
