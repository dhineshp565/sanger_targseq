#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { make_csv } from './modules/make_csv.nf'
include { tracy_assembly } from './modules/tracy_assembly.nf'
include { reverse_complement } from './modules/reverse_complement.nf'
include { sanger_reads } from './modules/sanger_reads.nf'
include { make_report } from './modules/make_report.nf'
include { blast_cons } from './modules/blast_cons.nf'
include { kraken2 } from './modules/kraken2.nf'
include { orfipy } from './modules/orfipy.nf'
include { abricate } from './modules/abricate.nf'
include { make_LIMSfile } from './modules/make_LIMSfile.nf'
include { mafft } from './modules/mafft.nf'
include { iqtree } from './modules/iqtree.nf'
include { ggtree } from './modules/ggtree.nf'
include { seq_length } from './modules/seq_length.nf'
include { trimal } from './modules/trimal.nf'



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
	
	blast_cons(reverse_complement.out.consensus,params.blastdb_path,params.blastdb_name)
	orfipy(reverse_complement.out.consensus)
	//generate 
	//if (params.kraken_db){
		//kraken2(reverse_complement.out.consensus,params.kraken_db)
	//}

		// Pair consensus and ORF files by sample name
	paired_consensus_orf = reverse_complement.out.consensus.join(orfipy.out.orf)
    .map { sample, consensus, orf -> tuple(sample, consensus, orf) }
	seq_length(paired_consensus_orf)
	

	rmd_file=file("${baseDir}/sanger_targseq.Rmd")
	rmd_file_case=file("${baseDir}/sanger_targseq_case.Rmd")

	
	refdir="${baseDir}/reference_sequences"
	mafft(make_csv.out,reverse_complement.out.cons_only.collect(),refdir)
	//trimal(mafft.out.collect())
	iqtree(mafft.out.collect())
	ggtree(iqtree.out.collect())
	make_report(make_csv.out,reverse_complement.out.cons_only.collect(),abricate.out.abricate.collect(),blast_cons.out.blast_formatted.collect(),ggtree.out.png,rmd_file,rmd_file_case,seq_length.out.collect(),orfipy.out.orf_only.collect())
	
}

