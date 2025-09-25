#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process trimal {
	label 'low'
	publishDir "${params.out_dir}/trimal/", mode: "copy"
  input:
    path msa

  output:
  path "*"

  script:
    """
   #
    for file in ${msa};do

		# get the prefix of the file name to check if it is a no_msa file
		prefix=\$(basename "\${file}" .fasta)

		# check if the prefix is "no_msa" and skip iqtree analysis if it is
		if [ "\$prefix" = "no_msa" ]; then
			echo "no virus sequences found for this sample, skipping trimal" > "\${prefix}_trimmed.fasta"

		else
			trimal -in "\${file}" -out "\${prefix}_trimmed.fasta" -automated1
		fi
	done
    """
}
