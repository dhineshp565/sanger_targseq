#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process seq_length {
	label 'low'
	publishDir "${params.out_dir}/seqlengths/", mode: "copy"
	input:
	tuple val(SampleName), path(consensus), path(orf)

	output:
	path("*.tsv")

	script:

	"""
	if grep -q "noconsensus" ${consensus}; then
	 	sed -i 's/noconsensus//g' ${consensus}
		cat ${consensus} ${SampleName}_ORF.fasta | seqkit fx2tab -n -l > ${SampleName}_seqlengths.tsv
	else 
		cat ${consensus} ${SampleName}_ORF.fasta | seqkit fx2tab -n -l > ${SampleName}_seqlengths.tsv
	fi

	"""

}
