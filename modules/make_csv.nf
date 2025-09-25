#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// make csv file with headers from the given input

process make_csv {
	publishDir "${params.out_dir}"
	label "low"
	input:
	path(fastq_input)
	output:
	path("samplelist.csv")
	
	script:
	"""
	makecsv.sh ${fastq_input}

	"""

}
