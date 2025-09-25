#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// process the forward and reverse reads from sanger sequencing
process sanger_reads {
    publishDir "${params.out_dir}/sanger_reads/", mode: "copy"
    label "medium"
	containerOptions "-v ${params.blastdb_path}:${params.blastdb_path}"
    input:
    tuple val(SampleName), path(SamplePath)
	path (blastdb_path)
	val(blastdb_name)

    output:
    tuple val(SampleName), path("${SampleName}_sangerreads.txt"), emit: sangerreads
	path("${SampleName}_reads_blast.csv"), emit:blast_reads

    script:

    """
	# Check if the SamplePath contains any .seq files (NRC-sanger sequencing facility) or *_premix.txt (macrogen - Sanger sequencing) files
    if ls ${SamplePath}/*.seq 1> /dev/null 2>&1; then
        for f in ${SamplePath}/*.seq; do
            echo ">\$f"
            cat "\$f"
            echo
        done > ${SampleName}_sangerreads.txt
	fi
    if ls ${SamplePath}/*_premix.txt 1> /dev/null 2>&1; then
        for f in ${SamplePath}/*_premix.txt; do
	         cat "\$f"
            echo
        done > ${SampleName}_sangerreads.txt
    
    fi
	# Run BLAST on forward and reverse reads
	blastn -db ${blastdb_path}/${blastdb_name} -query "${SampleName}_sangerreads.txt" -out ${SampleName}_reads_blast.csv -outfmt "7 qseqid sseqid length qcovs pident evalue staxids ssciname scomnames stitle" -max_target_seqs 5

    """
}
