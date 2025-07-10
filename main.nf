#!/usr/bin/env nextflow
nextflow.enable.dsl=2



// make csv file with headers from the given input

process make_csv {
	publishDir "${params.out_dir}"
	input:
	path(fastq_input)
	output:
	path("samplelist.csv")
	
	script:
	"""
	makecsv.sh ${fastq_input}

	"""

}

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

        if [ \$exit_code -eq 255 ]; then
            echo "tracy assemble failed with exit code 255" >> ${SampleName}_log.txt
            echo ">${SampleName}_consensus" > ${SampleName}_consR.fasta
            echo "noconsensus" >> ${SampleName}_consR.fasta
			echo ">${SampleName}_Noalignment" > ${SampleName}_alignment.fasta
			
	
        else
            mv *cons.fa ${SampleName}_consR.fasta
            sed -i 's/Consensus/${SampleName}_consensus/g' ${SampleName}_consR.fasta
            mv *align.fa ${SampleName}_alignment.fasta
			
        fi


	exit 0


	"""
}
// Process to get the reverse complement of the consensus sequence
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
	if grep -q "noconsensus" ${consensus}; then
		cat ${consensus} > ${SampleName}_consensus.fasta
	else
		# Use seqtk to get the reverse complement of the consensus sequence
		seqtk seq -r ${consensus} > ${SampleName}_consensus.fasta
	fi
	
	"""

}
// process the forward and reverse reads from sanger sequencing
process sanger_reads {
    publishDir "${params.out_dir}/sanger_reads/", mode: "copy"
    label "low"
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
	# Run BLAST
	blastn -db ${blastdb_path}/${blastdb_name} -query "${SampleName}_sangerreads.txt" -out ${SampleName}_reads_blast.csv -outfmt "7 qseqid sseqid length qcovs pident evalue staxids ssciname scomnames stitle" -max_target_seqs 5

    """
}


//make html report with rmarkdown
process make_report {
	publishDir "${params.out_dir}/",mode:"copy"
	
	label "low"
	input:
	path(csv)
	path(cons_only)
	path(abricate)
	path(blast_formatted)
	path (blast_reads)
	path(png)
	path(rmdfile)
	output:
	path("sanger_targseq*.html")
	script:
	"""
	
	cp ${csv} samples.csv
	cp ${rmdfile} report.Rmd
	cp ${png} tree.png
	
	Rscript -e 'rmarkdown::render(input="report.Rmd", params=list(csv="samples.csv",png="tree.png"), output_file = paste0("sanger_targseq_results_report_", Sys.Date(), "_", format(Sys.time(), "%H-%M-%S"), ".html"))'

	"""
}
// performs blast of the consensus sequences
process blast_cons {
	publishDir "${params.out_dir}/blast/",mode:"copy"
	containerOptions "-v ${params.blastdb_path}:${params.blastdb_path}"

	label "high"
	input:
	tuple val(SampleName),path(consensus)
	//tuple val(SampleNAme),path(sanger_reads)
	path (blastdb_path)
	val(blastdb_name)

	output:
	path("${SampleName}_report_blast.csv"), emit:blast_formatted
	path("${SampleName}_blast.csv")
	
	script:
	"""
	

	# Run BLAST
	blastn -db ${blastdb_path}/${blastdb_name} -query ${consensus} -out ${SampleName}_blast.csv -outfmt "7 qseqid sseqid length qcovs pident evalue staxids ssciname scomnames stitle" -max_target_seqs 5

	# Check for "# 0 hits found"
	if grep -q "# 0 hits found" "${SampleName}_blast.csv"; then
		# Write header and 'none' row to report
		echo -e "queryid\tsubject_id\talignment length\tquery_coverage\t%identity\tevalue\tstaxids\tsscinames\tscomnames\tstitle" > "${SampleName}_report_blast.csv"
		echo -e "none\tnone\tnone\tnone\tnone\tnone\tnone\tnone\tnone\tnone" >> "${SampleName}_report_blast.csv"
	else
		# remove lines starting with "#" and sort
		grep -v "#" "${SampleName}_blast.csv" | sort | uniq > "${SampleName}_report_blast.csv"
		# add header to the report
		sed -i '1i queryid\tsubject_id\talignment length\tquery_coverage\t%identity\tevalue\tstaxids\tsscinames\tscomnames\tstitle' "${SampleName}_report_blast.csv"
	fi

	




	"""

}
//kraken2 for classification
process kraken2 {
	publishDir "${params.out_dir}/kraken2/",mode:"copy"
	label "high"
	input:
	tuple val(SampleName),path (SamplePath)
	path(db_path)
	
	output:
	path ("${SampleName}_kraken.csv")
	path ("${SampleName}_kraken_report.csv"),emit:(kraken2_raw)
	
	script:
	"""
	kraken2 --db $db_path --output ${SampleName}_kraken.csv --report ${SampleName}_kraken_report.csv --threads 1 ${SamplePath}
	"""
}




process orfipy {
	label "low"
	publishDir "${params.out_dir}/orf",mode:"copy"
	input:
	tuple val(SampleName),path("${SampleName}_consensus.fasta")
	output:
	path ("${SampleName}_ORF.fasta")
	script:
	"""
	orfy.sh ${SampleName} ${SampleName}_consensus.fasta

	"""

}

process abricate{
	publishDir "${params.out_dir}/abricate/",mode:"copy"
	label "low"
	input:
	tuple val(SampleName),path(consensus)
	path(dbdir)
	output:
	path("${SampleName}_abricate.csv"),emit:abricate
	path("${SampleName}_withseq.csv"),emit:withseq
	script:
	"""
	typing.sh ${SampleName} ${consensus} ${dbdir}
	"""
	
}

process make_LIMSfile {
	publishDir "${params.out_dir}/LIMS/",mode:"copy"
	label "low"
	input:
	path (withseq)
	path (software_version_file)
	output:
	path("LIMSfile_*.tsv")
	
	script:
	"""
	date=\$(date '+%Y-%m-%d_%H-%M-%S')
	
	awk 'FNR==1 && NR!=1 { while (/^#F/) getline; } 1 {print}' *_withseq.csv > LIMSfile.tsv

	cat ${software_version_file} LIMSfile.tsv > LIMSfile_\${date}.tsv

	"""
}

process mafft {
    publishDir "${params.out_dir}/mafft/", mode: "copy"
    label "low"

    input:
    path (consensus)
    path (reference_sequences)

    output:
    path "*_msa.fasta"

    script:
    """
	mafft.sh ${reference_sequences}
    """
}

process iqtree {
	publishDir "${params.out_dir}/iqtree/", mode: "copy"
	label "low"

	input:
	path (msa)

	output:
	path "*.treefile"

	script:
	"""
    for file in ${msa};do

		prefix=\$(basename "\${file}" .fasta)

		if [ "\$prefix" = "no_msa" ]; then
			echo "no virus sequences found for this sample, skipping IQ-TREE analysis" > "\${prefix}.treefile"
		else
			iqtree2 -s "\${file}" -m MFP -bb 1000 -nt AUTO -pre "\${prefix}_iqtree"
		fi
	done
    """

}

process ggtree {
    publishDir "${params.out_dir}/ggtree/", mode: "copy"
    label "low"

    input:
    path (treefiles)

    output:
    path "*iqtree.png",emit:png
	

	when:

    script:
    """
	for treefile in ${treefiles}; do
    filename=\$(basename "\$treefile")
		if [[ "\$filename" != *no_msa* ]]; then
			plot_tree.R "\$treefile"
		else 
			echo "no virus sequences found for this sample, skipping ggtree analysis" > "\$filename.png"
		fi
	done
    """
}




workflow {
	data=Channel
	.fromPath(params.input)
	tracy_assembly(make_csv(data).splitCsv(header:true).map{row-> tuple(row.SampleName,row.SamplePath)})
	software_version_file=file("${baseDir}/software_version.tsv")
	sanger_reads(make_csv.out.splitCsv(header:true).map{row-> tuple(row.SampleName,row.SamplePath)},params.blastdb_path,params.blastdb_name)
	dbdir=file("${baseDir}/targseq")
	
	reverse_complement(tracy_assembly.out.consR)
	
	abricate(reverse_complement.out.consensus,dbdir)
	make_LIMSfile(abricate.out.withseq.collect(),software_version_file)
	//taxon=("${baseDir}/taxdb")
	blast_cons(reverse_complement.out.consensus,params.blastdb_path,params.blastdb_name)
	//orfipy(medaka.out.consensus)
	//generate 
	if (params.kraken_db){
		kraken2(reverse_complement.out.consensus,params.kraken_db)
	}

	rmd_file=file("${baseDir}/targseq.Rmd")

	
	refdir="${baseDir}/reference_sequences"
	mafft(reverse_complement.out.cons_only.collect(),refdir)
	iqtree(mafft.out.collect())
	ggtree(iqtree.out.collect())
	if (params.blastdb_name) {
		make_report(make_csv.out,reverse_complement.out.cons_only.collect(),abricate.out.abricate.collect(),blast_cons.out.blast_formatted.collect(),sanger_reads.out.blast_reads.collect(),ggtree.out.png,rmd_file)
	}
}

