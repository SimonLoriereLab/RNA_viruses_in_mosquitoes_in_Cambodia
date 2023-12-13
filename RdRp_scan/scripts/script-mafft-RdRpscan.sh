#!/bin/bash

#SBATCH -q geva
#SBATCH -p geva
#SBATCH --cpus-per-task=95
#SBATCH --mem=500000

module load fasta ruby
module load mafft/7.467

# _reference_dataset="/full_path_to/wd/RdRp_scan/RdRp-scan-main/Phylogenies/RdRp-scan.CLUSTALO_0.4.fasta"
# _sequences_to_add="/full_path_to/wd/RdRp_scan/analysis/R_ORFs.fasta"
# _output_aln="/full_path_to/wd/RdRp_scan/analysis/phylo/R_ORFs_aln.fasta"
#
# #########
# mafft --localpair \
# --inputorder --add ${_sequences_to_add} --keeplength \
# --thread ${SLURM_CPUS_PER_TASK} \
# ${_reference_dataset} > ${_output_aln}
#
# echo "done"

#### Updated the list of orfs and rerunning
_reference_dataset="/full_path_to/wd/RdRp_scan/analysis/phylo/RdRp-scan.CLUSTALO_0.4_dups_marked_clean_headers.fasta"
_sequences_to_add="/full_path_to/wd/RdRp_scan/analysis/R_232orfs_aa.fasta"
_output_aln="/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_aln.fasta"

#########
mafft --localpair \
--reorder \
--add ${_sequences_to_add} --keeplength \
--thread ${SLURM_CPUS_PER_TASK} \
${_reference_dataset} > ${_output_aln}

echo "done"
