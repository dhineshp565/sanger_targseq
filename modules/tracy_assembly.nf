#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//Assembling the consensus sequence using tracy
process tracy_assembly {
	publishDir "${params.out_dir}/tracy_assembly/",mode:"copy"
	label "low"
	
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path("${SampleName}_consR.fasta"), emit:consR
	path ("${SampleName}_alignment.fasta"), emit:alignment
	path("${SampleName}_log.txt"), emit:tracy_log
	
	
	script:
	"""
	set +e

	tracy assemble ${SamplePath}/*.ab1 > ${SampleName}_log.txt 2>&1
        exit_code=\$?

	# Check if tracy assemble failed and create empty files if so
        if [ \$exit_code -eq 255 ]; then
            echo "tracy assemble failed with exit code 255" >> ${SampleName}_log.txt
            echo ">${SampleName}_consensus" > ${SampleName}_consR.fasta
            echo "noconsensus" >> ${SampleName}_consR.fasta
			echo ">${SampleName}_Noalignment" > ${SampleName}_alignment.fasta
			
	# rename files and sequence headers if tracy assembly was successful
        else
            mv *cons.fa ${SampleName}_consR.fasta
            sed -i 's/Consensus/${SampleName}_consensus/g' ${SampleName}_consR.fasta
            mv *align.fa ${SampleName}_alignment.fasta
			
        fi


	exit 0


	"""

}
