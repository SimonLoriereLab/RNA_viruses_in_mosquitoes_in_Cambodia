#!/bin/bash

#SBATCH -q geva
#SBATCH -p geva
#SBATCH --cpus-per-task=95
#SBATCH --mem=500000

module load FastTree/2.1.11

# _aln="/full_path_to/wd/RdRp_scan/analysis/phylo/R_ORFs_aln_names_modified.fasta"
# _output_tree="/full_path_to/wd/RdRp_scan/analysis/phylo/R_ORFs_fasttree.tre"
#
# #########
#
# FastTree ${_aln} > ${_output_tree}
#
# echo "done"



########## RdRp scan verification
# _aln="/full_path_to/wd/RdRp_scan/analysis/phylo/RdRp-scan.CLUSTALO_0.4_dups_marked_clean_headers.fasta"
# _output_tree="/full_path_to/wd/RdRp_scan/analysis/phylo/RdRp-scan.CLUSTALO_0.4_fasttree.tre"
#
# #########
#
# FastTree ${_aln} > ${_output_tree}
#
# echo "done"


######### Updated contig list, rerun phylogeny
_aln="/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_aln.fasta"
_output_tree="/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree.tre"

#########

FastTree ${_aln} > ${_output_tree}

echo "done"
