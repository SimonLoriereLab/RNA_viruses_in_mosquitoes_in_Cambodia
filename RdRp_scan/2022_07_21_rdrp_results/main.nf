nextflow.enable.dsl=2

params.fasta = 'data/sequences.fasta'
params.results = 'results/'
params.profilepath = 'RdRp-scan/Profile_db_and_alignments/'
params.profilename = 'RdRp_HMM_profile_CLUSTALO.db'

process getORF {
    label 'emboss'

    publishDir "${params.results}/code_${code}/", mode: 'copy'

    input:
    path fasta
    each code

    output:
    tuple val(code), path("*_aa.fasta")

    script:
    """
    getorf -minsize 600 -table $code -find 0 -sequence ${fasta} -outseq ${fasta.baseName}_code${code}_aa.fasta
    """
}

process removeRedundency {
    label 'cdhit'

    input:
    tuple val(code), path(prot)

    output:
    tuple val(code), path("${prot.baseName}_nr.fasta")

    script:
    """
    cd-hit -i $prot -o ${prot.baseName}_nr.fasta -c 0.98 -T ${task.cpus} -M ${task.memory.toMega()}
    """
}

process scanProfiles{
    label 'hmmer'

    input:
    path profilepath
    val profilename
    tuple val(code), path(prot)

    output:
    tuple val(code), path("hmmscan_out.txt")
    tuple val(code), path("hmmscan_tblout.txt")
    tuple val(code), path("hmmscan_tblout_best.txt")
    tuple val(code), path("hmmscan_domtblout.txt")
    
    script:
    """
    hmmscan --cpu ${task.cpus} -o hmmscan_out.txt --tblout hmmscan_tblout.txt --domtblout hmmscan_domtblout.txt $profilepath/$profilename $prot
    awk '!x[\$3]++' hmmscan_tblout.txt > hmmscan_tblout_best.txt
    """
}

workflow{
   gencode = [1, 3, 4, 5, 6, 11, 16]
   fasta=file(params.fasta)
   profilepath=file(params.profilepath)
   prot = getORF(fasta,gencode)
   nrprot = removeRedundency(prot)
   nrprotsplit=nrprot.splitFasta(by: 1000, file:true)
   scanout = scanProfiles(profilepath, params.profilename, nrprotsplit)

   hmmout = scanout[0]
   hmmtblout = scanout[1]
   hmmbest = scanout[2]
   hmmdom = scanout[3]

   hmmout.collectFile(storeDir: "${params.results}/"){ n,f -> ["code_${n}/hmmscan_out.txt", f]}
   hmmtblout.collectFile(storeDir: "${params.results}/"){ n,f -> ["code_${n}/hmmscan_tblout.txt", f]}
   hmmbest.collectFile(storeDir: "${params.results}/"){ n,f -> ["code_${n}/hmmscan_tblout_best.txt", f]}
   hmmdom.collectFile(storeDir: "${params.results}/"){ n,f -> ["code_${n}/hmmscan_domtblout.txt", f]}
}
