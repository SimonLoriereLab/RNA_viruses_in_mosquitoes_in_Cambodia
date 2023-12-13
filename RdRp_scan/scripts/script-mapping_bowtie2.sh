#!/bin/bash
#SBATCH -J map

#SBATCH -p common
#SBATCH --qos=normal

#SBATCH -c 1
#SBATCH --mem=5000



###############################################################################
_directory=$1
_sample_id=$2
cd $_directory
_ref_id=`sed -n ${SLURM_ARRAY_TASK_ID}p $3`
_reference_folder=$4
_summary_folder=$5
_multimap1st_jobserrors="$6/slurm-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}.err"
echo "$_multimap1st_jobserrors"


###############################################################################

###############################################################################
##### Careful here, if the name of the pipeline changes
##### ( e. g. not auto_metaviromics-master but auto_metaviromics-some_other_version)
_sample_input_path="/full_path_to/wd/auto_metaviromics/analyzed_samples/${_sample_id}/Trimmed"
_sample_output_path="/full_path_to/wd/auto_metaviromics/analyzed_samples/${_sample_id}/Multiref_mapping_RdRpScan_Round1"

if [ ! -d "$_sample_output_path" ]; then
mkdir $_sample_output_path
fi

_reference_path="${_reference_folder}"
_consensus_path="${_summary_folder}"
###############################################################################

###############################################################################
cd $_sample_output_path

module load ivar/1.0.1
module load samtools/1.10
module load bowtie2/2.3.5.1


if [ -f "${_reference_path}/${_ref_id}.fasta" ]; then

  _ref_fasta_ext="fasta"


elif [ -f "${_reference_path}/${_ref_id}.fa" ]; then
  _ref_fasta_ext="fa"
fi

echo ${_ref_fasta_ext}



##### Using very fast local option
###
bowtie2 -p ${SLURM_CPUS_PER_TASK} \
--very-fast-local --no-unal \
-x ${_reference_path}/${_ref_id}-index \
-1 ${_sample_input_path}/${_sample_id}_R1.qual.fastq.gz \
-2 ${_sample_input_path}/${_sample_id}_R2.qual.fastq.gz \
-U ${_sample_input_path}/${_sample_id}_R12.unpair.fastq.gz \
> $_sample_output_path/${_ref_id}_${_sample_id}_bt2.sam 2> $_sample_output_path/${_ref_id}_${_sample_id}_bt2stderr.txt

echo "\n\n\n--------- very fast local bowtie2 mapping done ---------\n\n\n"






grep -v "XS:" $_sample_output_path/${_ref_id}_${_sample_id}_bt2.sam > \
$_sample_output_path/${_ref_id}_${_sample_id}_bt2_unique.sam

echo "\n\n\n--------- grep unique done ---------\n\n\n"


samtools view -bh -o $_sample_output_path/${_ref_id}_${_sample_id}_bt2_unique.bam \
$_sample_output_path/${_ref_id}_${_sample_id}_bt2_unique.sam

echo "\n\n\n--------- sambam done ---------\n\n\n"





#### 1 Thread
samtools sort \
-o $_sample_output_path/${_ref_id}_${_sample_id}_bt2_unique.sorted.bam \
$_sample_output_path/${_ref_id}_${_sample_id}_bt2_unique.bam \
-@ 1 -m 5G






echo "\n\n\n--------- bam sorted ---------\n\n\n"

echo "\n\n\n--------- bam indexing started ---------\n\n\n"

samtools index $_sample_output_path/${_ref_id}_${_sample_id}_bt2_unique.sorted.bam

echo "\n\n\n--------- bam indexing done ---------\n\n\n"



## generate 50% consensus with 3X as min cov if not we put an N
samtools mpileup -d 0 -aa -A -Q 0 $_sample_output_path/${_ref_id}_${_sample_id}_bt2_unique.sorted.bam | ivar consensus -p $_consensus_path/${_ref_id}_${_sample_id}_RdRpScan_Round1 -q 20 -t 0.5 -m 3 -n N




echo "\n\n\n--------- consensus generated 1 ---------\n\n\n"
module load bedtools/2.29.2

bedtools genomecov \
-d -ibam $_sample_output_path/${_ref_id}_${_sample_id}_bt2_unique.sorted.bam > $_sample_output_path/${_ref_id}_${_sample_id}.cov

echo "\n\n\n--------- coverage calculated, cov created ---------\n\n\n"

## samtools faidx $_consensus_path/${_ref_id}_${_sample_id}_RdRpScan_Round1.fa
samtools faidx ${_reference_path}/${_ref_id}.${_ref_fasta_ext}

echo "\n\n\n--------- faidx done ---------\n\n\n"




module load lofreq/2.1.4

# lofreq call -f $_consensus_path/${_ref_id}_${_sample_id}_RdRpScan_Round1.fa \
# -o $_sample_output_path/${_ref_id}_${_sample_id}.vcf \
# $_sample_output_path/${_ref_id}_${_sample_id}_bt2_unique.sorted.bam



#### 1 Thread
lofreq call --force -f ${_reference_path}/${_ref_id}.${_ref_fasta_ext} \
-o $_sample_output_path/${_ref_id}_${_sample_id}.vcf \
$_sample_output_path/${_ref_id}_${_sample_id}_bt2_unique.sorted.bam


# lofreq call-parallel --pp-threads ${SLURM_CPUS_PER_TASK} \
# --force -f ${_reference_path}/${_ref_id}.${_ref_fasta_ext} \
# -o $_sample_output_path/${_ref_id}_${_sample_id}.vcf \
# $_sample_output_path/${_ref_id}_${_sample_id}_bt2_unique.sorted.bam



echo "\n\n\n--------- varcall done, vcf created ---------\n\n\n"




cd ${_summary_folder}
Rscript $_directory/scripts/coverage_norm_bowtie2.R ${_directory} ${_sample_id} ${_ref_id} \
$_sample_output_path/${_ref_id}_${_sample_id}_bt2stderr.txt \
$_sample_output_path/${_ref_id}_${_sample_id}.cov \
$_sample_output_path/${_ref_id}_${_sample_id}.vcf \
${_summary_folder}
