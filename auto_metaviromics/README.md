# auto_metaviromics (simplified version)

Generic pipeline to trim, assemble and run sequence similarity search on multiple samples in parallel. In this project it was used only for trimming and de novo assembly, because the sequence similarity search was then run separately on non-redundant contigs from all samples combined.

### If running the pipeline fully with DIAMOND blast and tax lineage annotation there are a few things to prepare.

1. Place fastq.gz files into a data folder, in the working directory.

2. Prepare diamond DB annotated with NCBI taxonomy ID in any preferred directory.

    i. Download most recent [taxdmp.zip](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip) and unpack.
    
    ii. Download most recent [pdb.accession2taxid.gz](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz).
    
    iii. Download most recent database in fasta format (e.g. nr-prot).
    
    iv. Run job_diamond_db.sh, make sure the paths to all the arguments are correct.

3. Prepare NCBI taxonomic lineage database csv file.

    i. Install [**ncbitax2lin**](https://github.com/zyxue/ncbitax2lin), needs Python 3.7 (worked on 3.8.3 for me).
    
    ii. Follow **ncbitax2lin** to run it with **taxdmp**.

4. Prepare the pipeline directory in the working directory

```console
mkdir -pv auto_metaviromics/pipeline
```

5. Check the following scripts in *auto_metaviromics/pipeline* and provide all the needed info, mask/unmask lines depending on the task, set snakemake flags ...:

config.yaml

wrapper-jobarray_v1.sh

run_pipeline_maestro.sh

Snakefile_v326

Lineage_v1.2.3.R

6. Give execution permissions to the scripts.

```console
chmod -R -x auto_metaviromics/pipeline/
```

7. Before running better try a dry run adding --dry-run inside run_pipeline_maestro.sh as one of the snakemake arguments. Check the outputs in *working_directory/auto_metaviromics_jobsoutputs* to make sure the correct rules are running as specified in the config.

```console
sbatch wrapper-jobarray_v1.sh
```

8. Remove --dry-run and run the pipeline.

```console
sbatch wrapper-jobarray_v1.sh
```

Acknowledgements. The pipeline was initially prepared by Sohta Ishikawa and Etienne Simon-Loriere ([auto_metaviromics](https://github.com/saishikawa/auto_metaviromics)). It was then simplified with help from Salome Steinke, some functions removed and put into Snakemake.