#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process mafft {
    publishDir "${params.out_dir}/mafft/", mode: "copy"
    label "low"

    input:
	path (csv)
    path (consensus)
    path (reference_sequences)

    output:
    path "*_msa.fasta"

    script:
    """
	mafft.sh ${csv} ${reference_sequences}
    """
}
