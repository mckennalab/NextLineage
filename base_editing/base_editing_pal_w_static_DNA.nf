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
        .map{ tuple(it.sample, it.reference, it.targets, it.reads)}
        .into{sample_table_cutsites; sample_table_seq; sample_table_print}

results_path = "results"
primers_extension = 5

/*
 * Create the cutSites location file and validate their primers
 */
process CreateCutSitesFile {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/01_cutSites"

    input:
    set sampleId, reference, targets, reads from sample_table_cutsites

    output:
    set sampleId, "${sampleId}.fa", targets, reads, "${sampleId}.fa.cutSites", "${sampleId}.fa.primers" into reference_pack
    
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
 * Zip together read-pairs from the unmerged reads into a single, interleaved file
 */
process AlignReads {
    memory '12 GB'
    time '12h'
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/02_alignment"
    
    input:
    set sampleId, reference, targets, mergedReads, cutsites, primers from reference_pack
    
    output:
    set sampleId, reference, targets, "${sampleId}.merged.fasta.gz", cutsites, primers into aligned_reads
    
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


    gzip ${sampleId}.merged.fasta
    """
}

/*
 * Zip together read-pairs from the unmerged reads into a single, interleaved file
 
process AlignReads {
    memory '12 GB'
    time '12h'
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/02_alignment"
    
    input:
    set sampleId, reference, targets, reads, cutsites, primers from reference_pack
    
    output:
    set sampleId, reference, targets, "${sampleId}.merged.fasta.gz", cutsites, primers into aligned_reads
    
    script:

    """
    cp ${reads} merged.fastq.gz
    gunzip merged.fastq.gz

    /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/simd_test/Complete-Striped-Smith-Waterman-Library-1.1/src/ssw_test \
     -c -r -m 10 -x 8 -o 20 -e 1 \
     ${reference} \
     merged.fastq > ${sampleId}.merged.fasta

    gzip ${sampleId}.merged.fasta
    """
}*/


/*
 * Given the aligned reads, call the events over the recorders
 */
process CallEvents {
    memory '12 GB'
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/03_event_calls"

    input:
    set sampleId, reference, targets, merged_aligned_reads, cutsites, primers from aligned_reads
    
    
    output:
    set sampleId, reference, targets, "${sampleId}.stats.gz", cutsites, primers into(stats_output,stats_for_base_calls)
    path "${sampleId}.stats.gz" into stats_file_eval

    script:

    """
    cp ${merged_aligned_reads} merged.fasta.gz
    gunzip merged.fasta.gz

    java -jar -Xmx11g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/bin/MergeAndCall.jar \
    DeepSeq \
    -inputMerged=merged.fasta \
    -cutSites=${cutsites} \
    -outputStats=${sampleId}.stats \
    -primersEachEnd=${primers} \
    -primerMismatches=${params.primerMismatchCount} \
    -primersToCheck=${params.primersToCheck} \
    -requiredMatchingProp=${params.alignmentThreshold} \
    -requiredRemainingBases=${params.minimumReadLength} \
    -cutsiteWindow=5

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
    publishDir "$results_path/04_breakdown_files"

    input:
    set sampleId, reference, targets, stats, cutsites, primers from stats_output
    
    output:
    tuple sampleId, reference, targets, "${sampleId}.perBase", "${sampleId}.topReadEvents", "${sampleId}.topReadEventsNew", "${sampleId}.topReadCounts", "${sampleId}.allReadCounts", cutsites into breakdown_files
    
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

process StatsAssessment {
    publishDir "$results_path/05_stats_assessment"

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
    publishDir "$results_path/06_sample_plot"

    input:
    tuple sample, reference, targets, perBase, topReadEvents, topReadEventsNew, topReadCounts, allReadCounts, cutsites from breakdown_files

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

    cp ${perBase} ${sample}/${sample}.perBase
    cp ${topReadEvents} ${sample}/${sample}.topReadEvents
    cp ${topReadEventsNew} ${sample}/${sample}.topReadEventsNew
    cp ${topReadCounts} ${sample}/${sample}.topReadCounts
    cp ${allReadCounts} ${sample}/${sample}.allReadCounts
    cp ${cutsites} ${sample}/${sample}.cutSites

    echo var occurance_file = \\"${sample}.topReadCounts\\" >> ${sample}/JS_files.js
    echo var top_read_melted_to_base = \\"${sample}.topReadEventsNew\\" >> ${sample}/JS_files.js
    echo var per_base_histogram_data = \\"${sample}.perBase\\" >> ${sample}/JS_files.js
    echo var cut_site_file = \\"${sample}.cutSites\\" >> ${sample}/JS_files.js
    python /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/copy_to_web.py --input_dir ./${sample}/ --sample ${sample} --project ${params.project} --webdir ${params.webdir}

    """
}


/*
 * Generate a number of derivative files from the stats 
 */
process CallBaseEdits {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/07_base_editing"

    input:
    tuple sampleId, reference, targets, stats, cutsites, primers from stats_for_base_calls
    
    output:
    tuple sampleId, "${sampleId}.baseEditCalls" into basecalls
    
    script:

    """
    cp ${stats} ./sample.stats.gz
    gunzip ./sample.stats.gz

    bash /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/base_editing_multi_target_analysis.sh \
    sample.stats \
    ${sampleId} \
    ${sampleId}.baseEditCalls \
    ${cutsites}

    rm sample.stats
    """
}

/*
 * Generate a number of derivative files from the stats 
 */
process BaseEditingSummary {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/08_base_summary"

    input:
    tuple sampleId, baseEditing from basecalls
    
    output:
    tuple sampleId, "${sampleId}_editing_summary.txt", "${sampleId}_editing_positions.txt" into baseEditingSummary
    
    script:

    """
    Rscript /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/base_editing_summary.R \
    ${sampleId} \
    ${baseEditing}
    """
}
