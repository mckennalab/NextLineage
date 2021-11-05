#!/usr/bin/env nextflow

/*

========================================================================================
                         mckenna_lab/paired_end_lineage
========================================================================================
 Converts paired-end DNA sequencing of recorders sequences into lineage calls
----------------------------------------------------------------------------------------

git remote -v*/

/*
 * Read in the sample table and convert entries into a channel of sample information 
 */

// sample  umi     reference       targets fastq1  fastq2  barcode.fastq1  barcode.fastq2  barcode1        barcode2
Channel.fromPath( params.samplesheet )
        .splitCsv(header: true, sep: '\t')
        .map{ tuple(it.sample, it.umi, it.umiLength, it.reference, it.targets, it.primer5, it.primer3, it.fastq1, it.fastq2)}
        .into{sample_table_cutsites; sample_table_alignment}

results_path = "results"
primers_extension = 5

/*
 * Create the cutSites location file and validate their primers
 */
process CreateCutSitesFile {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/cutSites"

    input:
    set sampleId, umi, umiLength, reference, targets, fwd_primer, rev_primer, fastq1, fastq2 from sample_table_cutsites
    
    output:
    tuple sampleId, "${ref_name}.cutSites", "${ref_name}", "${ref_name}.primers" into reference_pack
    
    script:
    ref_name = new File(reference).getName()
    
    """
    cp ${reference} ${ref_name}
    scala /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/fasta_to_cutsites.scala \
        ${ref_name} \
        ${targets} \
        ${fwd_primer} \
        ${rev_primer} \
        ${primers_extension}
    """
}

/*
 * Filter input reads, requiring that at least one of the paired reads align to the lineage recorder sequence
 */
process SetupFilterReads {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/filter_by_alignment"

    input:
    set sampleId, umi, umiLength, reference, targets, fwd_primer, rev_primer, fastq1, fastq2 from sample_table_alignment
    
    output:
    tuple sampleId, "${sampleId}.r1.fq.gz", "${sampleId}.r2.fq.gz", "${sampleId}.flagstat_pre.txt", "${sampleId}.post.reads.txt" into filtered_reads
    
    when:
    params.genome_reference?.trim()

    script:
    ref_name = new File(reference).getName()
    
    """
    bwa mem ${params.genome_reference} ${fastq1} ${fastq2} | samtools sort -o output.bam -
    samtools flagstat output.bam > ${sampleId}.flagstat_pre.txt
    samtools index output.bam
    samtools view -H output.bam > header.sam
    samtools view output.bam ${params.bait_seq_name} | cut -f1 > IDs.txt
    LC_ALL=C grep -w -F -f IDs.txt <(samtools view output.bam) | cat header.sam - | samtools collate -u -O - | samtools fastq -1 ${sampleId}.r1.fq.gz -2 ${sampleId}.r2.fq.gz -0 /dev/null -s /dev/null -n -
    zcat ${sampleId}.r1.fq.gz | wc > ${sampleId}.post.reads.txt
    """
}


/*
 * Create the cutSites location file and validate their primers
 
process ExtractUmis {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/umis"


    input:
    set sampleId, umi, reference, targets, fwd_primer, rev_primer, fastq1, fastq2, indexFastq1, indexFastq2, index1, index2 from sample_table

    when:
    umi == "TRUE"

    output:
    path "${ref_name}.cutSites" into cutSites
    path "${ref_name}" into ref
    path "${ref_name}.primers" into primers

    script:
    ref_name = new File(reference).getName()
    
    """
    cp ${reference} ${ref_name}
    scala /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/fasta_to_cutsites.scala \
        ${ref_name} \
        ${targets} \
        ${fwd_primer} \
        ${rev_primer} \
        ${primers_extension}
    """
}
*/