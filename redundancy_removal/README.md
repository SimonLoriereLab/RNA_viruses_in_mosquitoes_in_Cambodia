1. Extract contigs from metaSPAdes assembly files and prepend sample IDs to each contig ID.

Launch [extract_rename_wrapper.sh](scripts/extract_rename_wrapper.sh) that will launch [run_extract_rename.sh](scripts/run_extract_rename.sh) then [extract_contigs_fasta_single.R](scripts/extract_contigs_fasta_single.R) for each sample in the [sample_list1.tsv](metadata/sample_list1.tsv)

2. Concatenate (I used cat) all the newly generated fasta files into one (here named all_samples_renamed_contigs.fasta but not uploaded because too large).

3. Remove redundancy with CD-HIT-EST, output in all_samples_renamed_contigs_redundancy_removed.fasta (not uploaded because too large).

Launch [script-cd-hit.sh](scripts/script-cd-hit.sh)
