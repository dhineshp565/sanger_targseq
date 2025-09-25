#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Process to get the reverse complement of the consensus sequence becuase tracy outputs consensus in reverse direction by default
process reverse_complement {
	label "low"
	publishDir "${params.out_dir}/reverse_complement/", mode: "copy"
	input:
	tuple val(SampleName), path(consensus)
	output:
	tuple val(SampleName),path("${SampleName}_consensus.fasta"), emit:consensus
	path("${SampleName}_consensus.fasta"), emit:cons_only

	script:
	"""
	# check if consensus exist and create empty file if not
	
	if grep -q "noconsensus" ${consensus}; then
		cat ${consensus} > ${SampleName}_consensus.fasta
	
	else
		# Use seqtk to get the reverse complement of the consensus sequence
		seqtk seq -r ${consensus} > ${SampleName}_consensus.fasta
	fi
	
	"""

}
