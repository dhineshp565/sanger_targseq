#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process orfipy {
	label "low"
	publishDir "${params.out_dir}/orf",mode:"copy"
	input:
	tuple val(SampleName),path(fasta)
	output:
	tuple val(SampleName),path ("${SampleName}_ORF.fasta"),emit:orf
	path("${SampleName}_ORF.fasta"),emit:orf_only
	script:
	"""
	orfy.sh ${SampleName} ${fasta}

	"""

}
