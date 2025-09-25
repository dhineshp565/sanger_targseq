#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Process to run IQ-TREE for phylogenetic analysis
// This process takes multiple MSA files and runs IQ-TREE on each of them
process iqtree {
	publishDir "${params.out_dir}/iqtree/", mode: "copy"
	label "low"

	input:
	path (msa)

	output:
	path "*.treefile"

	script:
	"""
	#run iqtree2 for each msa file
    for file in ${msa};do

		# get the prefix of the file name to check if it is a no_msa file
		prefix=\$(basename "\${file}" .fasta)

		# check if the prefix is "no_msa" and skip iqtree analysis if it is
		if [ "\$prefix" = "no_msa" ]; then
			echo "no virus sequences found for this sample, skipping IQ-TREE analysis" > "\${prefix}.treefile"

		else
			# run iqtree2 with MFP model selection and 1000 bootstrap replicates
			iqtree2 -s "\${file}" -m MFP -bb 1000 -nt AUTO -pre "\${prefix}"
		fi
	done
    """

}
