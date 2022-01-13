#!/usr/bin/env nextflow

/*
========================================================================================
                         mckenna_lab/cas12a_amplicon
========================================================================================
 Converts paired-end DNA sequencing of Cas12a lineage recorders into event calls
----------------------------------------------------------------------------------------
*/


/*
 * Read in the sample table and convert entries into a channel of sample information 
 */


Channel.fromPath( params.samplesheet )
        .splitCsv(header: true, sep: '\t')
        .map{ tuple(it.sample, it.reference, it.targets, it.read1, it.read2)}
        .into{sample_table_cutsites; sample_table_fastqc; sample_table_print}

results_path = "results"
primers_extension = 5

/*
 * Create the cutSites location file and validate their primers
 */
process CreateCutSitesFile {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/01_cutSites"

    input:
    set sampleId, reference, targets, read1, read2 from sample_table_cutsites

    output:
    set sampleId, "${sampleId}.fa", targets, read1, read2, "${sampleId}.fa.cutSites", "${sampleId}.fa.primers" into reference_pack
    
    script:
    
    """
    cp ${reference} ${sampleId}.fa

    scala -J-Xmx4g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/fasta_to_cutsites.py \
        --reference ${sampleId}.fa \
        --sites ${targets} \
        --forward_primer ${params.primer5} \
        --reverse_primer ${params.primer3} \
        --CRISPR_type ${params.crispr}
    """
}
/*
 * Fastqc
 */
process FastqcReport {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/02_fastqc"

    input:
    set sampleId, reference, targets, read1, read2 from sample_table_fastqc
    
    output:
    set sampleId, reference, targets, read1, read2, "${sampleId}_fastqc" into fastqc_results
    
    script:
    
    """
    mkdir ${sampleId}_fastqc
    /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/tools/FastQC/fastqc -o ${sampleId}_fastqc ${read1} ${read2}
    """
}

/*
 * Filter input reads, requiring that at least one of the paired reads align to the lineage recorder sequence
 */
process FilterGenomicReads {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/03_filter_genomic_reads"

    input:
    set sampleId, reference, targets, read1, read2, cutsites, primers from reference_pack
    
    output:
    set sampleId, reference, targets, "${sampleId}.read.fq.gz", "${sampleId}.barcode.fq.gz", cutsites, primers into genome_filtered_reads
    file "${sampleId}_pre_post_reads.txt" into filtered_reads_analysis
    
    when:
    !params.filter_against_genome
    
    script:

    """
    zcat ${read1}| wc > ${sampleId}.pre.reads.txt

    python /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/filter_by_alignment.py \
    --reference_genome ${params.organism_reference} \
    --reference_gestalt_name ${params.bait_seq_name} \
    --read1 ${read1} \
    --barcode ${read2} \
    --output_sample_prefix ${sampleId}

    zcat ${sampleId}.read.fq.gz | wc > ${sampleId}.post.reads.txt
    cat ${sampleId}.pre.reads.txt ${sampleId}.post.reads.txt > ${sampleId}_pre_post_reads.txt
    """
}

process MergeUmis {
    memory '20 GB'
    beforeScript 'chmod o+rw .'
    //errorStrategy 'ignore'
    publishDir "$results_path/05_umi_collapse"

    when:
    params.umiLength

    input:
    set sampleId, reference, targets, read1, read2, cutsites, primers from genome_filtered_reads

    output:
    set sampleId, reference, targets, "${sampleId}_read1.fq.gz", "${sampleId}_read2.fq.gz", cutsites, primers into extracted_umis
    tuple sampleId, "${sampleId}_read1.fq.gz", "${sampleId}_read2.fq.gz", "${sampleId}.umiCounts" into extracted_umis_stats

    script:

    """
    
    java -jar -Xmx16g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/bin/MergeAndCall.jar \
    UMIMerge \
    -inputReads1=${read1} \
    -inputReads2=${read2} \
    -umiStart=${params.umiStart} \
    -minimumUMIReads=${params.minUMIReadCount} \
    -umiLength=${params.umiLength} \
    -umiStatsFile=${sampleId}.umiCounts \
    -primers=${primers} \
    -primersToCheck=${params.primersToCheck} \
    -primerMismatches=${params.primerMismatchCount} \
    -outputReads1=${sampleId}_read1.fq.gz
    -outputReads2=${sampleId}_read2.fq.gz

    """
}

if (params.umiLength) {
    filtered_reads = extracted_umis
} else {
    extracted_umis = genome_filtered_reads
}

/*
 * merge paired reads on overlapping sequences using 
 */
process MergeReads {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/04_ng_trim"

    input:
    set sampleId, reference, targets, read1, read2, cutsites, primers from filtered_reads
    
    output:
    set sampleId, reference, targets, "${sampleId}_merged.fastq.gz", "${sampleId}_single_1.fastq.gz", "${sampleId}_single_2.fastq.gz", cutsites, primers into merged_reads
    
    script:
    """
    /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/NGmerge/NGmerge \
    -1 ${read1} -2 ${read2} -o ${sampleId}_merged.fastq -f ${sampleId}_single
    """
}

/*
 * Zip together read pairs from the unmerged reads into a single, interleaved file
 */
process InterleaveUnmergedReads {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/05_interleaved"

    input:
    set sampleId, reference, targets, mergedReads, read1, read2, cutsites, primers from merged_reads
    
    output:
    set sampleId, reference, targets, mergedReads, "${sampleId}_interleaved.fq.gz", cutsites, primers into interleaved_reads

    script:
    """
    scala -J-Xmx4g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/zip_two_read_files.scala \
    ${read1} ${read2} ${sampleId}_interleaved.fq _ _
    gzip ${sampleId}_interleaved.fq
    """
}

/*
 * Zip together read-pairs from the unmerged reads into a single, interleaved file
 */
process AlignReads {
    memory '32 GB'
    cpus 20
    time '12h'
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/06_alignment"
    
    input:
    set sampleId, reference, targets, mergedReads, interleavedReads, cutsites, primers from interleaved_reads
    
    output:
    set sampleId, reference, targets, "${sampleId}.merged.fasta.gz", "${sampleId}.interleavedReads.fasta.gz", cutsites, primers into aligned_reads
    
    script:

    """
    cp ${mergedReads} merged.fastq.gz
    gunzip merged.fastq.gz

    scala -J-Xmx4g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/align_merged_reads.scala \
    /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/resources/EDNAFULL.Ns_are_zero \
    ${params.aligner} \
    merged.fastq \
    ${reference} \
    ${sampleId}.merged.fasta \
    10 \
    0.5 \
    FALSE 

    cp ${interleavedReads} interleavedReads.fastq.gz
    gunzip interleavedReads.fastq.gz

    scala -J-Xmx4g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/align_merged_reads.scala \
    /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/resources/EDNAFULL.Ns_are_zero \
    ${params.aligner} \
    interleavedReads.fastq \
    ${reference} \
    ${sampleId}.interleavedReads.fasta \
    10 \
    0.5 \
    TRUE

    gzip ${sampleId}.merged.fasta
    gzip ${sampleId}.interleavedReads.fasta
    """
}


/*
 * Given the aligned reads, call the events over the recorders
 */
process CallEvents {
    memory '12 GB'
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/07_event_calls"

    input:
    set sampleId, reference, targets, merged_aligned_reads, interleaved_aligned_reads, cutsites, primers from aligned_reads
    
    
    output:
    set sampleId, reference, targets, "${sampleId}.stats.gz", cutsites, primers into stats_output
    path "${sampleId}.stats.gz" into stats_file_eval

    script:

    """
    cp ${merged_aligned_reads} merged.fasta.gz
    gunzip merged.fasta.gz
    cp ${interleaved_aligned_reads} interleaved.fasta.gz
    gunzip interleaved.fasta.gz

    java -jar -Xmx11g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/bin/MergeAndCall.jar \
    DeepSeq \
    -inputFileUnmerged=interleaved.fasta \
    -inputMerged=merged.fasta \
    -cutSites=${cutsites} \
    -outputStats=${sampleId}.stats \
    -primersEachEnd=${primers} \
    -sample=${sampleId} \
    -primerMismatches=${params.primerMismatchCount} \
    -primersToCheck=FORWARD \
    -requiredMatchingProp=${params.alignmentThreshold} \
    -requiredRemainingBases=${params.minimumReadLength} \
    -cutsiteWindow=5

    rm merged.fasta
    rm interleaved.fasta

    gzip ${sampleId}.stats
    """
}

/*
 * Generate a number of derivative files from the stats 
 */
process BreakdownFiles {
    beforeScript 'chmod o+rw .'
    errorStrategy 'ignore'
    publishDir "$results_path/08_breakdown_files"

    input:
    set sampleId, reference, targets, stats, cutsites, primers from stats_output
    
    output:
    tuple sampleId, reference, targets, "${sampleId}.perBase", "${sampleId}.topReadEvents", "${sampleId}.topReadEventsNew", "${sampleId}.topReadCounts", "${sampleId}.allReadCounts" into breakdown_files
    
    script:

    """
    cp ${stats} ./sample.stats.gz
    gunzip ./sample.stats.gz

    scala -J-Xmx4g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/stats_to_javascript_tables.scala \
    sample.stats \
    ${sampleId}.perBase \
    ${sampleId}.topReadEvents \
    ${sampleId}.topReadCounts \
    ${sampleId}.allReadCounts \
    ${cutsites} \
    ${sampleId}.topReadEventsNew \
    ${params.convert_unknown_to_none} \
    ${reference}

    rm sample.stats
    """
}

process ReadFilterQC {
    beforeScript 'chmod o+rw .'

    publishDir "$results_path/09_read_filter_assessment"

    input:
    file pre_post from filtered_reads_analysis.toList()
    
    output:
    path "genome_filter_assessment_mqc.tsv" into read_filter_assessment
    
    script:

    """
    python /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/read_qc.py \
    --pre_post_stats ${pre_post.collect().join(",")} \
    --output genome_filter_assessment_mqc.tsv
    """
}

process StatsAssessment {
    errorStrategy 'ignore'
    publishDir "$results_path/10_stats_assessment"

    input:
    file stats from stats_file_eval.toList()

    output:
    tuple "stat_meta_qc_assessment_mqc.tsv", "target_assessment_mqc.tsv" into stats_assessment
    
    script:

    """
    python /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/stats_qc.py \
    --stats ${stats.collect().join(",")} \
    --output_edits target_assessment_mqc.tsv \
    --output_stats stat_meta_qc_assessment_mqc.tsv
    """
}
