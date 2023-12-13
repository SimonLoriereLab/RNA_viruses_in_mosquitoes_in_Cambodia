#!/bin/bash

#SBATCH -q geva
#SBATCH -p geva
#SBATCH -c 1
#SBATCH --mem=1000

#SBATCH -J WRAP_META_PIPELINE
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|,$s\]$s\$|]|" \
        -e ":1;s|^\($s\)\($w\)$s:$s\[$s\(.*\)$s,$s\(.*\)$s\]|\1\2: [\3]\n\1  - \4|;t1" \
        -e "s|^\($s\)\($w\)$s:$s\[$s\(.*\)$s\]|\1\2:\n\1  - \3|;p" $1 | \
   sed -ne "s|,$s}$s\$|}|" \
        -e ":1;s|^\($s\)-$s{$s\(.*\)$s,$s\($w\)$s:$s\(.*\)$s}|\1- {\2}\n\1  \3: \4|;t1" \
        -e    "s|^\($s\)-$s{$s\(.*\)$s}|\1-\n\1  \2|;p" | \
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)-$s[\"']\(.*\)[\"']$s\$|\1$fs$fs\2|p" \
        -e "s|^\($s\)-$s\(.*\)$s\$|\1$fs$fs\2|p" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p" | \
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]; idx[i]=0}}
      if(length($2)== 0){  vname[indent]= ++idx[indent] };
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) { vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, vname[indent], $3);
      }
   }'
}

eval $(parse_yaml config.yaml)

_directory=$outdir
_data_folder=$data ### this has to be specified, otherwise downstream doesn't work
_Snakefile=$Snakefile
_cores=$default_cpus
_mem=$default_mem
_partition=$partition
_qos=$qos

echo $_directory
echo $_data_folder
echo $_Snakefile
echo $_cores
echo $_mem
echo $_partition
echo $_qos

_auto_metaviromics_jobserrors="$_directory/auto_metaviromics_jobserrors"
_auto_metaviromics_jobsoutputs="$_directory/auto_metaviromics_jobsoutputs"

cd $_directory

if [ ! -d "$_auto_metaviromics_jobserrors" ]; then
mkdir $_auto_metaviromics_jobserrors
fi

if [ ! -d "$_auto_metaviromics_jobsoutputs" ]; then
mkdir $_auto_metaviromics_jobsoutputs
fi

# Prepare info for job array :
# - list files to process
# ls -p $_directory/$_data_folder | grep -v / > listname.txt
# awk -F  "_" '{print $1}' listname.txt > listname2.txt # list of files in
# awk '!seen[$0]++' listname2.txt > sample_list.txt #Samples to process, by their names

_nb_jobs=`wc -l < sample_list_lp16.txt` #get number of files to process

echo $_nb_jobs

####   Running on MAESTRO
sbatch -J META_PIPELINE -c $_cores --mem=$_mem -p $_partition -q $_qos --wait --array=1-$_nb_jobs -o $_auto_metaviromics_jobsoutputs/slurm-%A_%a.out -e \
$_auto_metaviromics_jobserrors/slurm-%A_%a.err \
$_directory/auto_metaviromics/pipeline/run_pipeline_maestro.sh \
$_directory sample_list_lp16.txt $_data_folder $_Snakefile $_cores|| exit 1


exit 0
