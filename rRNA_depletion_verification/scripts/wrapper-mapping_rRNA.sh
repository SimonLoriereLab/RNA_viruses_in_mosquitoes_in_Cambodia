#!/bin/bash

#SBATCH -J wrap_mapAB
#SBATCH -p geva
#SBATCH -q geva
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000


module load R/4.0.2
###############################################################################
# It is assumed that trimmed data are in auto_metaviromics in separate folders for each sample.
###############################################################################
#To be changed:
_directory="/full_path_to/wd"
_scripts="$_directory/rRNA_depletion_verification/scripts"
_analysis="$_directory/rRNA_depletion_verification/analysis/Aedes_Culex_depleted"
_map_ref_tab="$_directory/rRNA_depletion_verification/metadata/table_map_refs_Aedes_Culex_rRNA.txt"
# _map_ref_tab="$_directory/rRNA_depletion_verification/metadata/table_map_refs_Aedes_Culex_rRNA_F1126.txt" ## run just for one sample F1126 without lofreq for which it freezes



#_map_ref_tab="all_to_all" #### this is an alternative when all references in ref folder to all samples. Don't use unless few pairs.



_reference_folder="$_directory/rRNA_depletion_verification/rRNA_reference_sequences"
_summary_folder="$_analysis/mapping_rRNA"  ##### the name of the summary and consensus output folder
_sample_list="$_directory/rRNA_depletion_verification/metadata/sample_list_Aedes_Culex_rRNA_mapping.txt" 
# _sample_list="$_directory/rRNA_depletion_verification/metadata/sample_list_Aedes_Culex_rRNA_mapping_F1126.txt" ## run just for one sample F1126 without lofreq for which it freezes
_temp_ref_folder="$_analysis/temp_ref_rRNA"


#### where the fastqs are located in each sample folder
_sample_input_folder="Trimmed"
#### where mapping results will be located in each sample folder
_sample_output_folder="rRNA_mapping"

#### or make one (sample name before _R1.fastq.gz or _R2.fastq.gz ,
###############################################################################

_multimap_jobserrors="$_analysis/rRNA_multimap_jobserrors_AB_r1"
_multimap_jobsoutputs="$_analysis/rRNA_multimap_jobsoutputs_AB_r1"

cd $_directory

if [ ! -d "${_analysis}" ]; then
mkdir -pv ${_analysis}
fi

if [ ! -d "$_multimap_jobserrors" ]; then
mkdir -pv $_multimap_jobserrors
fi

if [ ! -d "$_multimap_jobsoutputs" ]; then
mkdir -pv $_multimap_jobsoutputs
fi

if [ ! -d "${_summary_folder}" ]; then
mkdir -pv ${_summary_folder}
fi

if [ ! -d "${_temp_ref_folder}" ]; then
mkdir -pv ${_temp_ref_folder}
fi


_nb_jobs=`wc -l < $_sample_list` #Computes number of files to process
echo $_nb_jobs

sbatch --wait --array=1-$_nb_jobs -o $_multimap_jobsoutputs/slurm-%A_%a.out -e \
$_multimap_jobserrors/slurm-%A_%a.err \
$_scripts/subwrapper-mapping_rRNA.sh \
$_directory $_sample_list $_map_ref_tab $_reference_folder $_summary_folder \
$_multimap_jobsoutputs $_multimap_jobserrors $_scripts $_temp_ref_folder $_analysis \
$_sample_input_folder $_sample_output_folder || exit 1

exit 0
