#!/bin/bash

#SBATCH -J wrap_map

#SBATCH -p common
#SBATCH -q normal

##SBATCH -p geva
##SBATCH -q geva

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000


module load R/3.6.2
module load ivar/1.0.1
###############################################################################
# This pipeline is written in case trimmed data are in auto_metaviromics
# in separate folders for each sample.
###############################################################################
#To be changed:
_directory="/full_path_to/wd/RdRp_scan"

#_map_ref_tab="${_directory}/analysis/read_mapping/map_table_.txt" ####  THIS IS TO PREPARE. See EXAMPLE FOR SALOME table_map_refs_RdRpScan_Round1.txt. The refs are putative virus names, comma separated.

_map_ref_tab="all_to_all" #### this is an alternative when all references in ref folder to all samples. Don't use unless few pairs.

_reference_folder="${_directory}/analysis/read_mapping/contigs_correct_sense_priority1"
#_reference_folder="${_directory}/analysis/read_mapping/contigs_correct_sense_priority2"

_summary_folder="${_directory}/analysis/read_mapping/mapping_summary_bowtie2_RdRpScan_Round1"  ##### the name of the summary and consensus output folder

_sample_list="${_directory}/analysis/read_mapping/sample_list_RdRpScan_Round1_common_nodes.txt" #### prepared already

_reflist_files_dir="${_directory}/analysis/read_mapping/reflist_file_dir_bowtie2_RdRpScan_Round1"

#### or make one (sample name before _R1.fastq.gz or _R2.fastq.gz )
###############################################################################

_multimap_jobserrors="${_directory}/analysis/read_mapping/multimap_jobserrors_bowtie2_RdRpScan_Round1"
_multimap_jobsoutputs="${_directory}/analysis/read_mapping/multimap_jobsoutputs_bowtie2_RdRpScan_Round1"


module load bowtie2/2.3.5.1
for ref in $(ls $_reference_folder)
do
  ref_id=$(echo $ref | sed -r 's/.[fast]+$//')
  echo "\n\n"
  echo ${ref_id}
  echo "\n\n"
  bowtie2-build \
  -f ${_reference_folder}/${ref} \
  ${_reference_folder}/${ref_id}-index
done

echo "bowtie2-build done"

cd $_directory

if [ ! -d "$_multimap_jobserrors" ]; then
mkdir $_multimap_jobserrors
fi

if [ ! -d "$_multimap_jobsoutputs" ]; then
mkdir $_multimap_jobsoutputs
fi

if [ ! -d "${_summary_folder}" ]; then
mkdir ${_summary_folder}
fi


_nb_jobs=`wc -l < $_sample_list` #Computes number of files to process
echo $_nb_jobs

sbatch --wait --array=1-$_nb_jobs%9 -o $_multimap_jobsoutputs/slurm-%A_%a.out -e \
$_multimap_jobserrors/slurm-%A_%a.err \
$_directory/scripts/subwrapper-mapping_bowtie2.sh \
$_directory $_sample_list $_map_ref_tab $_reference_folder $_summary_folder \
$_multimap_jobserrors $_multimap_jobsoutputs $_reflist_files_dir || exit 1

exit 0
