#!/usr/bin/env nextflow

/*

========================================================================================
                         mckenna_lab/RNA_10X_single_cell_lineage
========================================================================================
 Converts paired-end DNA sequencing of recorders sequences into lineage calls
----------------------------------------------------------------------------------------
params.organism_reference
params.bait_seq_name
params.reference
params.targets
params.primer5
params.primer3
params.primers_extension
params.aligner
params.primerMismatchCount
params.alignmentThreshold
params.minimumReadLength
params.convert_unknown_to_none
params.trim_read_bases
*/


/*
 * Read in the sample table and convert entries into a channel of sample information 
 */


Channel.fromPath( params.samplesheet )
        .splitCsv(header: true, sep: '\t')
        .map{ tuple(it.sample, it.reference, it.sites, it.read, it.barcode, it.umi)}
        .into{sample_table_cutsites; sample_table_fastqc}

results_path = "results"
primers_extension = 5

/*
 * Create the cutSites location file and validate their primers
 */
process CreateCutSitesFile {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/01_cutSites"

    input:
    set sampleId, reference, targets, read, barcode, umi from sample_table_cutsites

    output:
    set sampleId, "${sampleId}.fa", targets, read, barcode, umi, "${sampleId}.fa.cutSites", "${sampleId}.fa.primers" into reference_pack
    
    script:
    
    """
    cp ${reference} ${sampleId}.fa

    scala -J-Xmx4g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/fasta_to_cutsites.scala \
        ${sampleId}.fa \
        ${targets} \
        ${params.primer5} \
        ${params.primer3} \
        ${params.primers_extension}
    """
}
/*
 * Fastqc
 */
process FastqcReport {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/02_fastqc"

    input:
    set sampleId, reference, targets, read, barcode, umi from sample_table_fastqc
    
    output:
    set sampleId, reference, targets, read, barcode, umi, "${sampleId}_fastqc" into fastqc_results
    
    script:
    
    """
    mkdir ${sampleId}_fastqc
    fastqc -o ${sampleId}_fastqc ${read} ${barcode} ${umi}
    """
}

/*
 * Filter input reads, requiring that at least one of the paired reads align to the lineage recorder sequence
 */
process FilterGenomicReads {
    errorStrategy 'ignore'
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/03_filter_genomic_reads"

    input:
    set sampleId, reference, targets, read, barcode, umi, cutsites, primers from reference_pack
    
    output:
    set sampleId, reference, targets, "${sampleId}.read.fq.gz", "${sampleId}.barcode.fq.gz", "${sampleId}.umi.fq.gz", cutsites, primers into filtered_reads
    file "${sampleId}_pre_post_reads.txt" into filtered_reads_analysis
    
    script:
    
    umi_string = "--umi " + umi
    fl = file(umi)
    if (!fl.exists()) {
        umi_string = ""
    }

    """
    zcat ${read}| wc > ${sampleId}.pre.reads.txt

    python /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/filter_by_alignment.py \
    --reference_genome ${params.organism_reference} \
    --reference_gestalt_name ${params.bait_seq_name} \
    --read1 ${read} \
    --barcode ${barcode} \
    --output_sample_prefix ${sampleId} \
    ${umi_string}

    zcat ${sampleId}.read.fq.gz | wc > ${sampleId}.post.reads.txt
    cat ${sampleId}.pre.reads.txt ${sampleId}.post.reads.txt > ${sampleId}_pre_post_reads.txt
    """
}

/*
 * move the second read fastq file to the front of the first read
 */
process MoveFirstRead {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/04_collapse_to_sequence"

    input:
    set sampleId, reference, targets, read, barcode, umi, cutsites, primers from filtered_reads
    
    output:
    set sampleId, reference, targets, "${sampleId}.umiMerged.fastq.gz", cutsites, primers into combined_reads

    script:
    
    umi_string = "--umi " + umi
    fl = file(umi)
    if (!fl.exists()) {
        umi_string = ""
    }

    """
    python /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/attach_read_UMIs.py --read1 ${read} --barcode ${barcode} --trimbases ${params.trim_read_bases} --outputfastq ${sampleId}.umiMerged.fastq.gz ${umi_string}
    """
}

process ExtractUmis {
    beforeScript 'chmod o+rw .'
    //errorStrategy 'ignore'
    publishDir "$results_path/05_umi_collapse"

    input:
    set sampleId, reference, targets, combined_reads, cutsites, primers from combined_reads

    output:
    set sampleId, reference, targets, "${sampleId}.fasta", cutsites, primers into extracted_umis
    tuple sampleId, "${sampleId}.fasta", "${sampleId}.umiCounts" into extracted_umis_stats

    script:

    """
    
    java -jar -Xmx16g /dartfs/rc/lab/M/McKennaLab/resources/GESTALT/SingleCellLineage/UMIMerge/target/scala-2.12/MergeAndCall.jar \
    UMIMerge \
    -inputReads1=${combined_reads} \
    -umiStart=${params.umiStart} \
    -minimumUMIReads=${params.minUMIReadCount} \
    -umiLength=${params.umiLength} \
    -umiStatsFile=${sampleId}.umiCounts \
    -primers=${primers} \
    -primersToCheck=${params.primersToCheck} \
    -primerMismatches=${params.primerMismatchCount} \
    -outputReads1=${sampleId}.fasta

    """
}


/*
 * Zip together read-pairs from the unmerged reads into a single, interleaved file
 */
process AlignReads {
    beforeScript 'chmod o+rw .'
    errorStrategy 'ignore'
    publishDir "$results_path/06_alignment"

    input:
    set sampleId, reference, targets, umis_fasta, cutsites, primers from extracted_umis
    
    output:
    set sampleId, reference, targets, "${sampleId}.merged.fasta.gz", cutsites, primers into aligned_reads
    
    script:

    """
    scala -J-Xmx16g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/align_merged_reads.scala \
    /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/resources/EDNAFULL.Ns_are_zero \
    ${params.aligner} \
    ${umis_fasta} \
    ${reference} \
    ${sampleId}.merged.fasta \
    10 \
    0.5

    gzip ${sampleId}.merged.fasta

    """
}


/*
 * Given the aligned reads, call the events over the recorders
 */
process CallEvents {
    beforeScript 'chmod o+rw .'
    errorStrategy 'ignore'
    publishDir "$results_path/07_event_calls"

    input:
    set sampleId, reference, targets, aligned_fasta, cutsites, primers from aligned_reads
    
    output:
    set sampleId, reference, targets, "${sampleId}.stats.gz", cutsites, primers into stats_output
    path "${sampleId}.stats.gz" into stats_file_eval

    script:

    """
    cp ${aligned_fasta} merged.fasta.gz
    gunzip merged.fasta.gz


    java -jar -Xmx16g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/single_cell/bin/MergeAndCall.jar \
    DeepSeq \
    -inputMerged=merged.fasta \
    -cutSites=${cutsites} \
    -outputStats=${sampleId}.stats \
    -primersEachEnd=${primers} \
    -sample=${sampleId} \
    -primerMismatches=${params.primerMismatchCount} \
    -primersToCheck=FORWARD \
    -requiredMatchingProp=${params.alignmentThreshold} \
    -requiredRemainingBases=${params.minimumReadLength}

    rm merged.fasta

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
