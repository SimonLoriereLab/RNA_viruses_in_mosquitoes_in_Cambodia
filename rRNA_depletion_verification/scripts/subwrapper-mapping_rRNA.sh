#!/bin/bash
#SBATCH -J subwrap_mapAB
#SBATCH -p geva
#SBATCH -q geva
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000

module load clc-assembly-cell/5.1.0

_directory=$1

_sample_id=`sed -n ${SLURM_ARRAY_TASK_ID}p $2`

_map_ref_tab=$3
_reference_folder=$4
_summary_folder=$5

_multimap_jobsoutputs_folder=$6
_multimap_jobserrors_folder=$7

_scripts=$8

_temp_ref_folder=$9

_analysis=${10}

_sample_input_folder=${11}
_sample_output_folder=${12}

###############################################################################

cd $_temp_ref_folder
Rscript $_scripts/get_ref_list.R ${_directory} ${_sample_id} ${_map_ref_tab} ${_reference_folder}

cd $_analysis
###############################################################################
_multimap_internal__jobsoutputs="$_multimap_jobsoutputs_folder/${_sample_id}__internal__outputs"
_multimap_internal__jobserrors="$_multimap_jobserrors_folder/${_sample_id}__internal__errors"


if [ ! -d "$_multimap_internal__jobsoutputs" ]; then
mkdir $_multimap_internal__jobsoutputs
fi

if [ ! -d "$_multimap_internal__jobserrors" ]; then
mkdir $_multimap_internal__jobserrors
fi
###############################################################################

_nb_jobs=`wc -l < $_temp_ref_folder/reflist_sample_${_sample_id}.txt` #Computes number of files to process
echo $_nb_jobs
echo $_temp_ref_folder/reflist_sample_${_sample_id}.txt

sbatch --wait --array=1-$_nb_jobs -o $_multimap_internal__jobsoutputs/slurm-%A_%a.out -e \
$_multimap_internal__jobserrors/slurm-%A_%a.err \
$_scripts/script-mapping_v1_rRNA.sh \
$_directory $_sample_id $_temp_ref_folder/reflist_sample_${_sample_id}.txt \
$_reference_folder $_summary_folder $_scripts $_sample_input_folder $_sample_output_folder || exit 1

echo "subwrapper done"
