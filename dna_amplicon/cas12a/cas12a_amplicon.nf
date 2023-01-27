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

    python /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/fasta_to_cutsites.py \
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

process MergeUmis {
    memory '20 GB'
    beforeScript 'chmod o+rw .'
    //errorStrategy 'ignore'
    publishDir "$results_path/03_umi_collapse"

    when:
    params.umiLength

    input:
    set sampleId, reference, targets, read1, read2, cutsites, primers from reference_pack

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
    -outputReads1=${sampleId}_read1.fq.gz \
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
    set sampleId, reference, targets, "out.extendedFrags.fastq.gz", "out.notCombined_1.fastq", "out.notCombined_2.fastq", cutsites, primers into merged_reads
    
    script:
    """
    /dartfs/rc/lab/M/McKennaLab/resources/tools/bin/flash --min-overlap 10 --max-mismatch-density 0.02 ${read1} ${read2}
    gzip out.extendedFrags.fastq
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
    memory '12 GB'
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

    scala -J-Xmx8g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/align_merged_reads.scala \
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

    scala -J-Xmx8g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/align_merged_reads.scala \
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
    tuple sampleId, reference, targets, "${sampleId}.perBase", "${sampleId}.topReadEvents", "${sampleId}.topReadEventsNew", "${sampleId}.topReadCounts", "${sampleId}.allReadCounts", "${sampleId}_event_length_histrogram.txt", cutsites into breakdown_files
    
    script:

    """
    cp ${stats} ./sample.stats.gz
    gunzip ./sample.stats.gz

    /home/Kiewit/f003w5r/.local/share/coursier/bin/scala -J-Xmx4g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/stats_to_javascript_tables.scala \
    sample.stats \
    ${sampleId}.perBase \
    ${sampleId}.topReadEvents \
    ${sampleId}.topReadCounts \
    ${sampleId}.allReadCounts \
    ${cutsites} \
    ${sampleId}.topReadEventsNew \
    ${params.convert_unknown_to_none} \
    ${reference}

    python /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/create_length_histogram.py \
    --compressed_stats_file ${stats} \
    --output ${sampleId}_event_length_histrogram.txt

    rm sample.stats
    """
}

process StatsAssessment {
    publishDir "$results_path/09_stats_assessment"

    input:
    file stats from stats_file_eval.toList()

    output:
    path "stat_meta_qc_assessment_mqc.tsv" into stats_assessment_qc
    path "target_assessment_mqc.tsv" into stats_assessment_targets
    
    script:

    """
    python /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/stats_qc.py \
    --stats ${stats.collect().join(",")} \
    --output_edits target_assessment_mqc.tsv \
    --output_stats stat_meta_qc_assessment_mqc.tsv
    """
}

process CreateSamplePlot {
    publishDir "$results_path/10_sample_plot"

    input:
    tuple sample, reference, targets, perBase, topReadEvents, topReadEventsNew, topReadCounts, allReadCounts, event_lengths, cutsites from breakdown_files

    output:
    path "${sample}/read_editing_mutlihistogram.html" into html_file
    path "${sample}/read_editing_mutlihistogram.js" into js_file
    path "${sample}/JS_files.js" into js_vars
    tuple "${sample}/${sample}.perBase", "${sample}/${sample}.topReadEvents", "${sample}/${sample}.topReadEventsNew", "${sample}/${sample}.topReadCounts", "${sample}/${sample}.allReadCounts", "${sample}/${sample}.cutSites" into per_base_info
    tuple sample, "${sample}/" into sample_dir
    
    script:
    """

    mkdir ${sample}
    cp /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/resources/plots/cas12a/read_editing_mutlihistogram.html ./${sample}/
    cp /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/resources/plots/cas12a/read_editing_mutlihistogram.js ./${sample}/
    cp /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/resources/plots/cas12a/edit_lengths.html ./${sample}/
    cp /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/resources/plots/cas12a/draw_figure_1_histograms.js ./${sample}/

    cp ${perBase} ${sample}/${sample}.perBase
    cp ${topReadEvents} ${sample}/${sample}.topReadEvents
    cp ${topReadEventsNew} ${sample}/${sample}.topReadEventsNew
    cp ${topReadCounts} ${sample}/${sample}.topReadCounts
    cp ${allReadCounts} ${sample}/${sample}.allReadCounts
    cp ${cutsites} ${sample}/${sample}.cutSites
    cp ${event_lengths} ${sample}/${sample}_event_lengths.txt

    echo var occurance_file = \\"${sample}.topReadCounts\\" >> ${sample}/JS_files.js
    echo var top_read_melted_to_base = \\"${sample}.topReadEventsNew\\" >> ${sample}/JS_files.js
    echo var per_base_histogram_data = \\"${sample}.perBase\\" >> ${sample}/JS_files.js
    echo var cut_site_file = \\"${sample}.cutSites\\" >> ${sample}/JS_files.js
    echo var event_file = \\"${sample}_event_lengths.txt\\" >> ${sample}/JS_files.js

    python /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/copy_to_web.py --input_dir ./${sample}/ --sample ${sample} --project ${params.project} --webdir ${params.webdir}

    """
}

