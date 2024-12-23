#!/usr/bin/env nextflow

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

params.R1 = "Read1_file"
params.R2 = "Read2_file"

params.out_fname = "output"

//params.k2_DB = "$projectDir/k2_DB/"
params.pathogen_DB = "$projectDir/CIWARS_Pathogen_DB/" // Bacterial Pathogen DB
params.non_prokaryote_DB = "$projectDir/non-prokaryote-DB/" // Human + RefSeq complete Fungi + RefSeq complete Protozoa + NCBI UniVec. To learn more, visit Kraken2 manual

params.ARG_DB = "$projectDir/DB/DeepARG-DB"
params.ARG_DB_LEN = "$projectDir/DB/DeepARG_DB.len"

params.rpoB_DB = "$projectDir/DB/rpoB"
params.rpoB_DB_LEN = "$projectDir/DB/rpoB.len"

process QC{

    //publishDir "$projectDir", mode: "copy"

    input:
    path R1
    path R2

    output:
    path "${params.out_fname}_qc_R1.fastq.gz", emit: qc_R1
    path "${params.out_fname}_qc_R2.fastq.gz", emit: qc_R2
	
    """
  	$projectDir/fastp -i $R1 -I $R2 -o ${params.out_fname}_qc_R1.fastq.gz -O ${params.out_fname}_qc_R2.fastq.gz --detect_adapter_for_pe --trim_poly_g --trim_poly_x --low_complexity_filter --average_qual 10 --thread 8
    """

}

process filter_non_prokaryote_reads {

    //publishDir "$projectDir", mode: "copy"

    input:
    path qc_R1
    path qc_R2

    output:
    path "${params.out_fname}_ucseqs_1.fastq.gz", emit: filtered_R1
    path "${params.out_fname}_ucseqs_2.fastq.gz", emit: filtered_R2

    """
    kraken2 --threads 16 -db ${params.non_prokaryote_DB} --gzip-compressed --output ${params.out_fname}_k2_out.txt --unclassified-out ${params.out_fname}_ucseqs#.fastq --paired $qc_R1 $qc_R2
    gzip ${params.out_fname}_ucseqs_1.fastq ${params.out_fname}_ucseqs_2.fastq
    rm ${params.out_fname}_k2_out.txt
    """

}


process kraken2 {

    publishDir "$projectDir", mode: "copy"

    input:
    path qc_R1
    path qc_R2

    output:
    path "${params.out_fname}.k2report", emit: k2report

    """
    kraken2 --use-names -db ${params.pathogen_DB} --threads 16 --report ${params.out_fname}.k2report --gzip-compressed --paired $qc_R1 $qc_R2 > ${params.out_fname}.k2report
    """

}

process merge_reads {

    //publishDir "$projectDir", mode: "copy"

    input:
    path R1
    path R2

    output:
    path "${params.out_fname}_final_read_merged.fasta", emit: merged_reads

    """
    vsearch --fastq_mergepairs $R1 --reverse $R2 --fastaout ${params.out_fname}_merged.fasta --fastaout_notmerged_fwd ${params.out_fname}_R1_unmerged.fasta --fastaout_notmerged_rev ${params.out_fname}_R2_unmerged.fasta --threads 4
    cat ${params.out_fname}_merged.fasta ${params.out_fname}_R1_unmerged.fasta ${params.out_fname}_R2_unmerged.fasta > ${params.out_fname}_final_read_merged.fasta
    """
}

process align_reads {

    //publishDir "$projectDir", mode: "copy"

    input:
    path R

    output:
    path "${params.out_fname}_aln_ARG.tsv", emit: aln_ARG
    path "${params.out_fname}_aln_rpoB.tsv", emit: aln_rpoB

    """
    $projectDir/diamond blastx -e 1e-10 --id 80 -k 1 --threads 32 -d ${params.ARG_DB} -q $R -o ${params.out_fname}_aln_ARG.tsv --outfmt 6
    $projectDir/diamond blastx -e 1e-10 --id 40 -k 1 --threads 32 -d ${params.rpoB_DB} -q $R -o ${params.out_fname}_aln_rpoB.tsv --outfmt 6
    """
}


process rpoB_normalization_ARG {

    publishDir "$projectDir", mode: "copy"

    input:
    path aln_ARG
    path aln_rpoB

    output:
    path "${params.out_fname}_rpoB_ARG_norm.tsv", emit: rpoB_ARG_norm

    """
    python $projectDir/rpoB_abund.py -a $aln_ARG -r $aln_rpoB -la ${params.ARG_DB_LEN} -lr ${params.rpoB_DB_LEN} -o ${params.out_fname}_rpoB_ARG_norm.tsv --db deeparg
    """

}

process rpoB_normalization_ARG_drug_wise {

    publishDir "$projectDir", mode: "copy"

    input:
    path rpoB_ARG_norm

    output:
    path "${params.out_fname}_drug_wise_rpoB_norm.tsv"

    """
    python $projectDir/drug_class_wise_norm.py -n $rpoB_ARG_norm -o ${params.out_fname}
    """

}


workflow {

    fw_file_ch = Channel.from(params.R1)
    rev_file_ch = Channel.from(params.R2)
    
    qc_ch = QC(fw_file_ch, rev_file_ch)

    filter_ch = filter_non_prokaryote_reads(qc_ch.qc_R1, qc_ch.qc_R2)

    k2_ch = kraken2(filter_ch.filtered_R1, filter_ch.filtered_R2)

    merged_reads_ch = merge_reads(qc_ch.qc_R1, qc_ch.qc_R2)

    aln_ch = align_reads(merged_reads_ch.merged_reads)

    rpoB_ARG_norm_ch = rpoB_normalization_ARG(aln_ch.aln_ARG, aln_ch.aln_rpoB)

    rpoB_normalization_ARG_drug_wise(rpoB_ARG_norm_ch.rpoB_ARG_norm)

}
