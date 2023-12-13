#!/bin/bash
#SBATCH -J subwrap_map

#SBATCH -p common
#SBATCH --qos=normal


#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000

_directory=$1
_sample_id=`sed -n ${SLURM_ARRAY_TASK_ID}p $2`
_map_ref_tab=$3
_reference_folder=$4
_summary_folder=$5
_multimap_jobserrors=$6
_multimap_jobsoutputs=$7
_reflist_files_dir=$8

###############################################################################

Rscript $_directory/scripts/get_ref_list.R ${_directory} ${_sample_id} $_map_ref_tab $_reference_folder $_reflist_files_dir

###############################################################################

_multimap1st_jobserrors="$_multimap_jobserrors/${_sample_id}_1st_errors"
_multimap1st_jobsoutputs="$_multimap_jobsoutputs/${_sample_id}_1st_outputs"

if [ ! -d "$_multimap1st_jobserrors" ]; then
mkdir $_multimap1st_jobserrors
fi

if [ ! -d "$_multimap1st_jobsoutputs" ]; then
mkdir $_multimap1st_jobsoutputs
fi

_nb_jobs=`wc -l < $_reflist_files_dir/reflist_sample_${_sample_id}.txt` #Computes number of files to process
echo $_nb_jobs
echo $_reflist_files_dir/reflist_sample_${_sample_id}.txt

sbatch --wait --array=1-$_nb_jobs -o $_multimap1st_jobsoutputs/slurm-%A_%a.out -e \
$_multimap1st_jobserrors/slurm-%A_%a.err \
$_directory/scripts/script-mapping_bowtie2.sh \
$_directory $_sample_id $_reflist_files_dir/reflist_sample_${_sample_id}.txt $_reference_folder $_summary_folder \
$_multimap1st_jobserrors || exit 1

echo "subwrapper done"
