#!/bin/bash

#SBATCH -q geva
#SBATCH -p geva
#SBATCH --cpus-per-task=47
#SBATCH --mem=250000

module load fasta ruby
module load mafft/7.467
module load goalign/0.3.1
module load IQ-TREE/2.0.6
module load FastTree/2.1.11


####### I ran commands below consecutively and intependantly for each clade.
####### If any of them need to be run again uncomment the needed block and launch the script.

####### cISF_clade
#### V3. Added missed sequences, new outgroup.
### Run in 2 rounds
### First align without the outgroup
## consensus sequences and refseq set of polyprotein sequences
# _seq_fname="cISF_full_polyprot_aa_v3.fasta"
# _aln_fname="cISF_full_polyprot_aa_aln_v3.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="cISF_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# Add an outgroup to the alignment
# _aln_fname="cISF_full_polyprot_aa_aln_v3_NSregion.fasta"
# _out_aln_fname="cISF_full_polyprot_aa_aln_v3_NSregion_outg.fasta"

### Then add the outgroup
#
# _outg_fname="LC540441_Tabanus_rufidens_flavivirus.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="cISF_clade"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
#
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"




####### third_ISF_clade
### Align without the outgroup
## consensus sequences and refseq set of polyprotein sequences
# _seq_fname="third_ISF_clade_aa.fasta"
# _aln_fname="third_ISF_clade_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="third_ISF_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}

####### V2: added one missing seq from spiders and substitited Kobu-sho seq for a longer version
## consensus sequences and refseq set of polyprotein sequences
# _seq_fname="third_ISF_clade_aa_v2.fasta"
# _aln_fname="third_ISF_clade_aa_v2_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="third_ISF_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
#
#
# # _aln_fname="third_ISF_clade_aa_v2_aln.fasta"
# # _directory="/full_path_to/wd/RdRp_scan"
# # _vir_group="third_ISF_clade"
# # _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
#
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA
#
# echo "tree is done"


### Phasi_clade
# _seq_fname="Phasi_clade_aa_v3.fasta"
# _aln_fname="Phasi_clade_aa_v3_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Phasi_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
#
#
# # _aln_fname="Phasi_clade_aa_v2_aln.fasta"
# # _directory="/full_path_to/wd/RdRp_scan"
# # _vir_group="Phasi_clade"
# # _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
#
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o NC_031298_Fly_phasivirus
#
# echo "tree is done"










# ### Quaranja_clade v2
# _seq_fname="Quaranja_clade_aa_v2.fasta"
# _aln_fname="Quaranja_clade_aa_v2_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Quaranja_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# # mafft --genafpair --maxiterate 1000 --thread 95 \
# # ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Quaranja_clade_aa_v2_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg_fname="MW784037_Hainan_orthomyxo-like_virus_2.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# # mafft --thread 95 --quiet --inputorder --keeplength --add \
# # ${_outg_folder}/${_outg_fname} \
# # ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap_clean.fasta"
# # goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap_clean.fasta -m LG+F+I+G4 -bb 1000 \
# -nt AUTO -mem 500G -st AA -o MW784037_Hainan_orthomyxo-like_virus_2
#
# echo "tree is done"




# ### Quaranja_clade v3
# _seq_fname="Quaranja_clade_aa_v3.fasta"
# _aln_fname="Quaranja_clade_aa_v3_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Quaranja_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Quaranja_clade_aa_v3_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg_fname="MW288163_Neuropteran_orthomyxo-related_virus_OKIAV190.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o MW288163_Neuropteran_orthomyxo-related_virus_OKIAV190
#
# echo "tree is done"


# ### Quaranja_clade v4
# _seq_fname="Quaranja_clade_aa_v4.fasta"
# _aln_fname="Quaranja_clade_aa_v4_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Quaranja_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# # mafft --genafpair --maxiterate 1000 --thread 95 \
# # ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Quaranja_clade_aa_v4_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg_fname="MW784037_Hainan_orthomyxo-like_virus_2.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# # mafft --thread 95 --quiet --inputorder --keeplength --add \
# # ${_outg_folder}/${_outg_fname} \
# # ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap_clean.fasta"
# # goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
#
# FastTree ${_tre_folder}/${_out_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_out_aln_fname}.gap_clean.fasttree.tre
#
# # iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA -o MW784037_Hainan_orthomyxo-like_virus_2


# echo "tree is done"
#
# ### Quaranja_clade v5
# _seq_fname="Quaranja_clade_aa_v5.fasta"
# _aln_fname="Quaranja_clade_aa_v5_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Quaranja_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Quaranja_clade_aa_v5_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg_fname="MW784037_Hainan_orthomyxo-like_virus_2.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o MW784037_Hainan_orthomyxo-like_virus_2
#
# echo "tree is done"









# ### Quaranja_clade v31
# _seq_fname="Quaranja_clade_aa_v31.fasta"
# _aln_fname="Quaranja_clade_aa_v31_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Quaranja_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Quaranja_clade_aa_v31_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg_fname="MN053836_Guadeloupe_mosquito_quaranja-like_virus_3.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o MN053836_Guadeloupe_mosquito_quaranja-like_virus_3
#
# echo "tree is done"



# ### Quaranja_clade v32
# _seq_fname="Quaranja_clade_aa_v32.fasta"
# _aln_fname="Quaranja_clade_aa_v32_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Quaranja_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Quaranja_clade_aa_v32_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg_fname="OL700154_XiangYun_orthomyxo-like_virus_2.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o OL700154_XiangYun_orthomyxo-like_virus_2
#
# echo "tree is done"


### Quaranja_clade v33
# _seq_fname="Quaranja_clade_aa_v33.fasta"
# _aln_fname="Quaranja_clade_aa_v33_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Quaranja_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Quaranja_clade_aa_v33_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg_fname="BK059432_Aedes_orthomyxo-like_virus_2.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o BK059432_Aedes_orthomyxo-like_virus_2
#
# echo "tree is done"














### Thogoto_clade
# _seq_fname="Thogoto_clade_aa.fasta"
# _aln_fname="Thogoto_clade_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Thogoto_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# # _aln_fname="Thogoto_clade_aa_aln.fasta"
# # _directory="/full_path_to/wd/RdRp_scan"
# # _vir_group="Thogoto_clade"
# # _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA
#
# echo "tree is done"
#
# ### Try without gap trimming, alignment is already quite short.
# ######### tree
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA
#
# echo "tree is done"


















# ## Sobeli_clade
# _seq_fname="Sobeli_clade_aa.fasta"
# _aln_fname="Sobeli_clade_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Sobeli_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# # _aln_fname="Sobeli_clade_aa_aln.fasta"
# # _directory="/full_path_to/wd/RdRp_scan"
# # _vir_group="Sobeli_clade"
# # _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA
#
# echo "tree is done"
#
# ### Try without gap trimming, if alignment is already quite short.
# ######### tree
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA
#
# echo "tree is done"
#
#







# ### Sobeli_clade v21
# _seq_fname="Sobeli_clade_aa_v21.fasta"
# _aln_fname="Sobeli_clade_aa_v21_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Sobeli_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Sobeli_clade_aa_v21_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="NC_032193_Beihai_sobemo-like_virus_26"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"



# ### Sobeli_clade v211
# _seq_fname="Sobeli_clade_aa_v211.fasta"
# _aln_fname="Sobeli_clade_aa_v211_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Sobeli_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Sobeli_clade_aa_v211_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MW199223_Virus_sp."
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"





# ### Sobeli_clade v22
# _seq_fname="Sobeli_clade_aa_v22.fasta"
# _aln_fname="Sobeli_clade_aa_v22_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Sobeli_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Sobeli_clade_aa_v22_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MW239542_Riboviria_sp."
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"





# ### Sobeli_clade v23
# _seq_fname="Sobeli_clade_aa_v23.fasta"
# _aln_fname="Sobeli_clade_aa_v23_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Sobeli_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Sobeli_clade_aa_v23_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="NC_033279_Sanxia_water_strider_virus_10"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"




# ### Sobeli_clade v231
# _seq_fname="Sobeli_clade_aa_v231.fasta"
# _aln_fname="Sobeli_clade_aa_v231_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Sobeli_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Sobeli_clade_aa_v231_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MZ443571_Vespula_vulgaris_Sobemo-like_virus_1"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"






# ### Sobeli_clade v2311
# _seq_fname="Sobeli_clade_aa_v2311.fasta"
# _aln_fname="Sobeli_clade_aa_v2311_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Sobeli_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Sobeli_clade_aa_v2311_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MW434821_Kellev_virus"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
#
# echo "tree is done"


# ### Sobeli_clade v2311_dedup
# _seq_fname="Sobeli_clade_aa_v2311_dedup.fasta"
# _aln_fname="Sobeli_clade_aa_v2311_dedup_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Sobeli_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Sobeli_clade_aa_v2311_dedup_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MW434821_Kellev_virus"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# ### Rerun w/o some duplicates. Model already known.
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m LG+R3 -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"






### Sobeli_clade v23111
# _seq_fname="Sobeli_clade_aa_v23111.fasta"
# _aln_fname="Sobeli_clade_aa_v23111_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Sobeli_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Sobeli_clade_aa_v23111_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MT361063_Kisumu_mosquito_virus"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"







# ### Sobeli_clade v24
# _seq_fname="Sobeli_clade_aa_v24.fasta"
# _aln_fname="Sobeli_clade_aa_v24_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Sobeli_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Sobeli_clade_aa_v24_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MZ443576_Vespula_vulgaris_Luteo-like_virus_2"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"


# ### Sobeli_clade v241
# _aln_fname="Sobeli_clade_nt_v241_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Sobeli_clade"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/nt_aln"
#
# ## Then add the outgroup
# _outg="LC512857_Wenzhou_sobemo-like_virus_3"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/nt_tree"
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_aln_folder}/${_aln_fname} \
# ${_tre_folder}/${_aln_fname}
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname} -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st DNA -o ${_outg}
#
# echo "tree is done"





# ## Monjiviricetes_clade
# _seq_fname="Monjiviricetes_clade_aa.fasta"
# _aln_fname="Monjiviricetes_clade_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Monjiviricetes_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA
#
# echo "tree is done"
#
# ### Try without gap trimming, if alignment is already quite short.
# ######### tree
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA
#
# echo "tree is done"





## Monjiviricetes_clade, fast alignment
# _seq_fname="Monjiviricetes_clade_aa.fasta"
# _aln_fname="Monjiviricetes_clade_aa_aln1.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Monjiviricetes_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA
#
# echo "tree is done"
#
# ### Try without gap trimming, if alignment is already quite short.
# ######### tree
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA
#
# echo "tree is done"



# ## Monjiviricetes_clade1
# _seq_fname="Monjiviricetes_clade1_aa.fasta"
# _aln_fname="Monjiviricetes_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Monjiviricetes_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Monjiviricetes_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="KM817612_Shuangao_Fly_Virus_1"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"





# ## Monjiviricetes_clade2
# _seq_fname="Monjiviricetes_clade2_aa.fasta"
# _aln_fname="Monjiviricetes_clade2_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Monjiviricetes_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Monjiviricetes_clade2_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="ON872562_Bat_faecal_associated_anphe-like_virus_1"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"




# ## Monjiviricetes_clade3
# _seq_fname="Monjiviricetes_clade3_aa.fasta"
# _aln_fname="Monjiviricetes_clade3_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Monjiviricetes_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Monjiviricetes_clade3_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MZ771217_Rhabdoviridae_sp_"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"



# ## Monjiviricetes_clade4
# _seq_fname="Monjiviricetes_clade4_aa.fasta"
# _aln_fname="Monjiviricetes_clade4_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Monjiviricetes_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Monjiviricetes_clade4_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="NC_031276_Wuhan_Ant_Virus"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"



# ## Monjiviricetes_clade5
# _seq_fname="Monjiviricetes_clade5_aa.fasta"
# _aln_fname="Monjiviricetes_clade5_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Monjiviricetes_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Monjiviricetes_clade5_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="OP589941_Rhabdoviridae_sp_"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"






# ## Monjiviricetes_clade6
# _seq_fname="Monjiviricetes_clade6_aa.fasta"
# _aln_fname="Monjiviricetes_clade6_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Monjiviricetes_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Monjiviricetes_clade6_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="OK491503_Xiangshan_rhabdo-like_virus_5"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"


#
# # Peribunyaviridae_clade, fast alignment
# _seq_fname="Peribunyaviridae_clade_aa.fasta"
# _aln_fname="Peribunyaviridae_clade_aa_aln1.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Peribunyaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# # mafft --retree 1 --thread 95 \
# # ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# # goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"
#
# # ### Try without gap trimming, if alignment is already quite short.
# # ######### tree
# # cp ${_out_aln_file} \
# # ${_tre_folder}/${_aln_fname}.fasta
# #
# # cd ${_tre_folder}
# #
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
# #
# # echo "tree is done"



# ## Peribunyaviridae_clade1
# _seq_fname="Peribunyaviridae_clade1_aa.fasta"
# _aln_fname="Peribunyaviridae_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Peribunyaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Peribunyaviridae_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MG995845_Yunnan_manyleaf_rhizome_virus"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"


# ## Peribunyaviridae_clade2
# _seq_fname="Peribunyaviridae_clade2_aa.fasta"
# _aln_fname="Peribunyaviridae_clade2_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Peribunyaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Peribunyaviridae_clade2_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MF190051_Barns_Ness_serrated_wrack_bunya_phlebo-like_virus_1"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"



## Peribunyaviridae_clade3. No outgroup.
# _seq_fname="Peribunyaviridae_clade3_aa.fasta"
# _aln_fname="Peribunyaviridae_clade3_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Peribunyaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA
#
# echo "tree is done"




# ## Peribunyaviridae_clade4
# _seq_fname="Peribunyaviridae_clade4_aa.fasta"
# _aln_fname="Peribunyaviridae_clade4_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Peribunyaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Peribunyaviridae_clade4_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="ON872551_Bat_faecal_associated_bunyavirus_8"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"


# ## Peribunyaviridae_clade5
# _seq_fname="Peribunyaviridae_clade5_aa.fasta"
# _aln_fname="Peribunyaviridae_clade5_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Peribunyaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Peribunyaviridae_clade5_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MG967342_Blechomonas_maslovi_leishbunyavirus_1"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"



# ## Peribunyaviridae_clade6
# _seq_fname="Peribunyaviridae_clade6_aa.fasta"
# _aln_fname="Peribunyaviridae_clade6_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Peribunyaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Peribunyaviridae_clade6_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MT153358_Dipteran_phenui-related_virus_OKIAV281"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"





# Phasmaviridae_clade, fast alignment
# _seq_fname="Phasmaviridae_clade_aa.fasta"
# _aln_fname="Phasmaviridae_clade_aa_aln1.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Phasmaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"
#
# # ### Try without gap trimming, if alignment is already quite short.
# # ######### tree
# # cp ${_out_aln_file} \
# # ${_tre_folder}/${_aln_fname}.fasta
# #
# # cd ${_tre_folder}
# #
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
# #
# # echo "tree is done"








# ## Phasmaviridae_clade1
# _seq_fname="Phasmaviridae_clade1_aa.fasta"
# _aln_fname="Phasmaviridae_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Phasmaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Phasmaviridae_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MT153380_Collembolan_phasma-related_virus_OKIAV223"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"


# ## Phasmaviridae_clade2
# _seq_fname="Phasmaviridae_clade2_aa.fasta"
# _aln_fname="Phasmaviridae_clade2_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Phasmaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Phasmaviridae_clade2_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MN164619_Pink_bollworm_virus_2"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"





# ### Iflaviridae_clade, fast alignment
# _seq_fname="Iflaviridae_clade_aa.fasta"
# _aln_fname="Iflaviridae_clade_aa_aln1.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Iflaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"
#
# # ### Try without gap trimming, if alignment is already quite short.
# # ######### tree
# # cp ${_out_aln_file} \
# # ${_tre_folder}/${_aln_fname}.fasta
# #
# # cd ${_tre_folder}
# #
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
# #
# # echo "tree is done"







# ## Iflaviridae_clade1
# _seq_fname="Iflaviridae_clade1_aa.fasta"
# _aln_fname="Iflaviridae_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Iflaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Iflaviridae_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="KX779452_Rolda_virus"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"


# ## Iflaviridae_clade2
# _seq_fname="Iflaviridae_clade2_aa.fasta"
# _aln_fname="Iflaviridae_clade2_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Iflaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Iflaviridae_clade2_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="OM622379_Picornavirales_sp_"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"




# ## Iflaviridae_clade3
# _seq_fname="Iflaviridae_clade3_aa.fasta"
# _aln_fname="Iflaviridae_clade3_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Iflaviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Iflaviridae_clade3_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="NC_033455_Wuhan_insect_virus_13"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"





# ### Dicistroviridae_clade, fast alignment
# _seq_fname="Dicistroviridae_clade_aa.fasta"
# _aln_fname="Dicistroviridae_clade_aa_aln1.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Dicistroviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"
#
# # ### Try without gap trimming, if alignment is already quite short.
# # ######### tree
# # cp ${_out_aln_file} \
# # ${_tre_folder}/${_aln_fname}.fasta
# #
# # cd ${_tre_folder}
# #
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
# #
# # echo "tree is done"





# ## Dicistroviridae_clade1
# _seq_fname="Dicistroviridae_clade1_aa.fasta"
# _aln_fname="Dicistroviridae_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Dicistroviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Dicistroviridae_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="JF423195_Big_Sioux_River_virus"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"






# ### Caliciviridae_clade, fast alignment
# _seq_fname="Caliciviridae_clade_aa.fasta"
# _aln_fname="Caliciviridae_clade_aa_aln1.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Caliciviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"
#
# # ### Try without gap trimming, if alignment is already quite short.
# # ######### tree
# # cp ${_out_aln_file} \
# # ${_tre_folder}/${_aln_fname}.fasta
# #
# # cd ${_tre_folder}
# #
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
# #
# # echo "tree is done"




# ## Caliciviridae_clade1
# _seq_fname="Caliciviridae_clade1_aa.fasta"
# _aln_fname="Caliciviridae_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Caliciviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Caliciviridae_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MW826415_Polycipiviridae_sp_"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"




# ### Picorna_like_clade, fast alignment
# _seq_fname="Picorna_like_clade_aa.fasta"
# _aln_fname="Picorna_like_clade_aa_aln1.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Picorna_like_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"
#
# # ### Try without gap trimming, if alignment is already quite short.
# # ######### tree
# # cp ${_out_aln_file} \
# # ${_tre_folder}/${_aln_fname}.fasta
# #
# # cd ${_tre_folder}
# #
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
# #
# # echo "tree is done"







# ## Picorna_like_clade1
# _seq_fname="Picorna_like_clade1_aa.fasta"
# _aln_fname="Picorna_like_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Picorna_like_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Picorna_like_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MZ209811_Guiyang_Solinvi-like_virus_3"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"






# ### Mesoniviridae_clade, fast alignment
# _seq_fname="Mesoniviridae_clade_aa.fasta"
# _aln_fname="Mesoniviridae_clade_aa_aln1.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Mesoniviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"
#
# # ### Try without gap trimming, if alignment is already quite short.
# # ######### tree
# # cp ${_out_aln_file} \
# # ${_tre_folder}/${_aln_fname}.fasta
# #
# # cd ${_tre_folder}
# #
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
# #
# # echo "tree is done"





# ## Mesoniviridae_clade1
# _seq_fname="Mesoniviridae_clade1_aa.fasta"
# _aln_fname="Mesoniviridae_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Mesoniviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Mesoniviridae_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MN609866_Leveillula_taurica_associated_alphamesonivirus_1"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"




# ## Mesoniviridae_clade11
# _seq_fname="Mesoniviridae_clade11_aa.fasta"
# _aln_fname="Mesoniviridae_clade11_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Mesoniviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Mesoniviridae_clade11_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="NC_036586_Alphamesonivirus_10"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"




# ### Reovirales_clade, fast alignment
# _seq_fname="Reovirales_clade_aa.fasta"
# _aln_fname="Reovirales_clade_aa_aln1.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Reovirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"
#
# # ### Try without gap trimming, if alignment is already quite short.
# # ######### tree
# # cp ${_out_aln_file} \
# # ${_tre_folder}/${_aln_fname}.fasta
# #
# # cd ${_tre_folder}
# #
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
# #
# # echo "tree is done"





# ## Reovirales_clade1
# _seq_fname="Reovirales_clade1_aa.fasta"
# _aln_fname="Reovirales_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Reovirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Reovirales_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="NC_027567_Lutzomyia_reovirus_1"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"


# ## Reovirales_clade2
# _seq_fname="Reovirales_clade2_aa.fasta"
# _aln_fname="Reovirales_clade2_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Reovirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Reovirales_clade2_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="KX884702_Hubei_reo-like_virus_14"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"



# ### Ghabrivirales_clade, fast alignment
# _seq_fname="Ghabrivirales_clade_aa.fasta"
# _aln_fname="Ghabrivirales_clade_aa_aln1.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Ghabrivirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"
#
# # ### Try without gap trimming, if alignment is already quite short.
# # ######### tree
# # cp ${_out_aln_file} \
# # ${_tre_folder}/${_aln_fname}.fasta
# #
# # cd ${_tre_folder}
# #
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
# #
# # echo "tree is done"




# ### Ghabrivirales_clade, fast alignment. Redo after replacement of CP (ORF1) sequences with the corresponding RdRP (ORF2) sequences.
# ### Also sequence Aedes_aegypti_totivirus_F1506C12C58C205F1506 was manually corrected in Geneious to remove the stop codon in the ORF2.
# _seq_fname="Ghabrivirales_clade_aa1.fasta"
# _aln_fname="Ghabrivirales_clade_aa1_aln1.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Ghabrivirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"
#
# # ### Try without gap trimming, if alignment is already quite short.
# # ######### tree
# # cp ${_out_aln_file} \
# # ${_tre_folder}/${_aln_fname}.fasta
# #
# # cd ${_tre_folder}
# #
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
# #
# # echo "tree is done"








# ## Ghabrivirales_clade1
# _seq_fname="Ghabrivirales_clade1_aa.fasta"
# _aln_fname="Ghabrivirales_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Ghabrivirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Ghabrivirales_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MZ209749_Fushun_totivirus_5"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"




# ## Ghabrivirales_clade2
# _seq_fname="Ghabrivirales_clade2_aa.fasta"
# _aln_fname="Ghabrivirales_clade2_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Ghabrivirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Ghabrivirales_clade2_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="ON746544_Jiamusi_Totiv_tick_virus_1"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"


# ## Ghabrivirales_clade3
# _seq_fname="Ghabrivirales_clade3_aa.fasta"
# _aln_fname="Ghabrivirales_clade3_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Ghabrivirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Ghabrivirales_clade3_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="LC516853_Barley_aphid_RNA_virus_8"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"


# ## Ghabrivirales_clade4
# _seq_fname="Ghabrivirales_clade4_aa.fasta"
# _aln_fname="Ghabrivirales_clade4_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Ghabrivirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Ghabrivirales_clade4_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="OW529222_Porphyridium_purpureum_toti-like_virus_1"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"


# ## Ghabrivirales_clade5
# _seq_fname="Ghabrivirales_clade5_aa.fasta"
# _aln_fname="Ghabrivirales_clade5_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Ghabrivirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Ghabrivirales_clade5_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="NC_032465_Beihai_blue_swimmer_crab_virus_3"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"



# ## Ghabrivirales_clade6
# _seq_fname="Ghabrivirales_clade6_aa.fasta"
# _aln_fname="Ghabrivirales_clade6_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Ghabrivirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Ghabrivirales_clade6_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MZ218558_Totiviridae_sp_"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"




# ### Negevirus_clade, fast alignment
# # _seq_fname="Negevirus_clade_aa.fasta"
# _aln_fname="Negevirus_clade_aa_aln2.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Negevirus_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# # mafft --retree 1 --thread 95 \
# # ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# # goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"
#
# # ### Try without gap trimming, if alignment is already quite short.
# # ######### tree
# # cp ${_out_aln_file} \
# # ${_tre_folder}/${_aln_fname}.fasta
# #
# # cd ${_tre_folder}
# #
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
# #
# # echo "tree is done"




# ## Negevirus_clade1
# _seq_fname="Negevirus_clade1_aa.fasta"
# _aln_fname="Negevirus_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Negevirus_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Negevirus_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="OL700120_XiangYun_hepe-virga-like_virus_7"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"




# ## Negevirus_clade2
# _seq_fname="Negevirus_clade2_aa.fasta"
# _aln_fname="Negevirus_clade2_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Negevirus_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Negevirus_clade2_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="NC_032437_Beihai_anemone_virus_1"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"


# ## Negevirus_clade3
# _seq_fname="Negevirus_clade3_aa.fasta"
# _aln_fname="Negevirus_clade3_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Negevirus_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Negevirus_clade3_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="KU754539_Boutonnet_virus"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"


# ## Negevirus_clade4
# _seq_fname="Negevirus_clade4_aa.fasta"
# _aln_fname="Negevirus_clade4_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Negevirus_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Negevirus_clade4_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="NC_032824_Hubei_virga-like_virus_16"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"


# ## Negevirus_clade5
# _seq_fname="Negevirus_clade5_aa.fasta"
# _aln_fname="Negevirus_clade5_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Negevirus_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Negevirus_clade5_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MH614308_Bombus-associated_virus_Vir1"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"



# ## Negevirus_clade6
# _seq_fname="Negevirus_clade6_aa.fasta"
# _aln_fname="Negevirus_clade6_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Negevirus_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Negevirus_clade6_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MT227315_Bemisia_tabaci_virga-like_virus_1"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"





# ### Negevirus_Tanay_clade, fast alignment
# _seq_fname="Negevirus_Tanay_clade_aa.fasta"
# _aln_fname="Negevirus_Tanay_clade_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Negevirus_Tanay_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"
#
# # ### Try without gap trimming, if alignment is already quite short.
# # ######### tree
# # cp ${_out_aln_file} \
# # ${_tre_folder}/${_aln_fname}.fasta
# #
# # cd ${_tre_folder}
# #
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
# #
# # echo "tree is done"

# ## Negevirus_Tanay_clade1
# _seq_fname="Negevirus_Tanay_clade1_aa.fasta"
# _aln_fname="Negevirus_Tanay_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Negevirus_Tanay_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Negevirus_Tanay_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="OK491483_Xiangshan_martelli-like_virus_1"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"













# ### Tymovirales_clade, fast alignment
# _seq_fname="Tymovirales_clade_aa.fasta"
# _aln_fname="Tymovirales_clade_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Tymovirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"

# ### Try without gap trimming, if alignment is already quite short.
# ######### tree
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA
#
# echo "tree is done"




# ## Tymovirales_clade1
# _seq_fname="Tymovirales_clade1_aa.fasta"
# _aln_fname="Tymovirales_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Tymovirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Tymovirales_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="MZ556269_Sichuan_mosquito_tymo-like_virus"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"





### Permutotetraviridae_clade, fast alignment
# _seq_fname="Permutotetraviridae_clade_aa.fasta"
# _aln_fname="Permutotetraviridae_clade_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Permutotetraviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"
#
# # ### Try without gap trimming, if alignment is already quite short.
# # ######### tree
# # cp ${_out_aln_file} \
# # ${_tre_folder}/${_aln_fname}.fasta
# #
# # cd ${_tre_folder}
# #
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
# #
# # echo "tree is done"


# ## Permutotetraviridae_clade1
# _seq_fname="Permutotetraviridae_clade1_aa.fasta"
# _aln_fname="Permutotetraviridae_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Permutotetraviridae_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Permutotetraviridae_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="NC_030845_Egaro_virus"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"




# ## Nodamuvirales_clade, fast alignment
# _seq_fname="Nodamuvirales_clade_aa.fasta"
# _aln_fname="Nodamuvirales_clade_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Nodamuvirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
#
# #### fast alignment
# mafft --retree 1 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
# ### Gap trimming
# _in_aln_file="${_aln_folder}/${_aln_fname}"
# _out_aln_file="${_aln_folder}/${_aln_fname}.gap_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.2 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.gap_clean.fasta
#
# cd ${_tre_folder}
#
# # iqtree2 -s ${_tre_folder}/${_aln_fname}.gap_clean.fasta -m MFP -bb 1000 \
# # -nt AUTO -mem 500G -st AA
#
# FastTree ${_tre_folder}/${_aln_fname}.gap_clean.fasta > ${_tre_folder}/${_aln_fname}.gap_clean.fasttree.tre
#
#
#
#
# echo "tree is done"

# ### Try without gap trimming, if alignment is already quite short.
# ######### tree
# cp ${_out_aln_file} \
# ${_tre_folder}/${_aln_fname}.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_aln_fname}.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA
#
# echo "tree is done"


## Nodamuvirales_clade1
# _seq_fname="Nodamuvirales_clade1_aa.fasta"
# _aln_fname="Nodamuvirales_clade1_aa_aln.fasta"
#
# ##########
# _directory="/full_path_to/wd/RdRp_scan"
# _vir_group="Nodamuvirales_clade"
# _seq_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
# if [ ! -d "$_aln_folder" ]; then
# mkdir -pv $_aln_folder
# fi
# ########
# mafft --genafpair --maxiterate 1000 --thread 95 \
# ${_seq_folder}/${_seq_fname} > ${_aln_folder}/${_aln_fname}
#
# ## Add an outgroup to the alignment
# _out_aln_fname="Nodamuvirales_clade1_aa_aln_outg.fasta"
#
# ## Then add the outgroup
# _outg="BK059724_Opistofel_virus"
# _outg_fname="${_outg}.fasta"
# _outg_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_seq"
# _aln_folder="${_directory}/analysis/phylo/${_vir_group}/aln/aa_aln"
#
# ########
# mafft --thread 95 --quiet --inputorder --keeplength --add \
# ${_outg_folder}/${_outg_fname} \
# ${_aln_folder}/${_aln_fname} > ${_aln_folder}/${_out_aln_fname}
#
# echo "alignment is done"
#
# _tre_folder="${_directory}/analysis/phylo/${_vir_group}/tree/aa_tree"
#
#
# ### Gap trimming only for more than 50% missing
# _in_aln_file="${_aln_folder}/${_out_aln_fname}"
# _out_aln_file="${_aln_folder}/${_out_aln_fname}.gap50_clean.fasta"
# goalign clean sites -i ${_in_aln_file} -c 0.5 > ${_out_aln_file}
# echo "gap trimming is done"
#
# ######### tree
#
# if [ ! -d "$_tre_folder" ]; then
# mkdir -pv $_tre_folder
# fi
#
# cp ${_out_aln_file} \
# ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta
#
# cd ${_tre_folder}
#
# iqtree2 -s ${_tre_folder}/${_out_aln_fname}.gap50_clean.fasta -m MFP -bb 1000 \
# -nt AUTO -mem 500G -st AA -o ${_outg}
#
# echo "tree is done"
