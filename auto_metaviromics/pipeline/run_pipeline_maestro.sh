#!/bin/bash

_directory=$1
_data_folder=$3
_Snakefile=$4
_cores=$5

echo

cd $_directory
_list=`sed -n ${SLURM_ARRAY_TASK_ID}p $2`
_sample_path="${_directory}/${_data_folder}/${_list}"
cd $_directory/auto_metaviromics/pipeline


#### This is for MAESTRO cluster
module load graalvm/ce-java8-20.0.0 Trimmomatic/0.39 bowtie2/2.3.5.1 samtools/1.10 \
bedtools/2.29.2 diamond/2.0.6 Python/3.8.1 htslib/1.10 primer3/2.4.0 ClustalW/2.1 \
golden/3.4.1 tabix/0.2.6 perl/5.30.1 aragorn/1.2.38 blast+/2.10.0	infernal/1.1.3 \
prodigal/2.6.3 ncbitools/20170106 minced/0.4.2 signalp/4.1 R/3.6.2 gnuplot/5.2.8 \
MUMmer/3.23 fasta/3.6 ruby/2.7.0 mafft/7.467 megahit/1.2.9 vcftools/0.1.16 SPAdes/3.15.2 snakemake/6.6.1


echo 'sources and modules are OK'

# snakemake -s $_Snakefile --cores $_cores --configfile config.yaml \
# --config sample=${_list} --keep-going --nolock \
# --reason || exit 1


# snakemake -s $_Snakefile --cores $_cores --configfile config.yaml \
# --config sample=${_list} --keep-going --nolock \
# --forcerun sum_DIAMOND \
# --reason --dry-run || exit 1

snakemake -s $_Snakefile --cores $_cores --configfile config.yaml \
--config sample=${_list} --keep-going --nolock \
--reason --forcerun filtering_low_qual_paired merge_unpaired assembling_metaSPAdes || exit 1


echo "DONE"
