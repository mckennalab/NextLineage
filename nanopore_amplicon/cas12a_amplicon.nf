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
    set sampleId, reference, targets, reads from sample_table_cutsites

    output:
    set sampleId, path("${sampleId}.fa"), targets, reads, path("${sampleId}.fa.cutSites"),path( "${sampleId}.fa.primers") into reference_pack
    
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
process AlignSequences {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/02_aligned_sequences"

    input:
    set sampleId, path(reference), targets, reads, path(cutsites), path(primers) from reference_pack
    
    output:
    set sampleId, path(reference), targets, path("${sampleId}_aligned.fasta"), path(cutsites), path(primers) into aligned_reads
    
    script:
    
    """
    
    /dartfs-hpc/rc/lab/M/McKennaLab/projects/sequential_targets_evan/2022_12_09_EW_mega_25X25/NextLineage/nanopore_amplicon/clique/target/release/clique \
    --reference ${reference} \
    --output-base base \
    --output ${sampleId}_aligned.fasta \
    --read1 ${reads} \
    --read-template test \
    --threads ${params.threads} \
    --use-capture-sequences \
    --only-output-captured-ref align
    """
}

/*
 * Given the aligned reads, call the events over the recorders
 */
process CallEvents {
    memory '12 GB'
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/03_event_calls"

    input:
    set sampleId, path(reference), targets, path(aligned_fastq), path(cutsites), path(primers) from aligned_reads
    
    
    output:
    set sampleId, path(reference), targets, path("${sampleId}.stats.gz"), path(cutsites), path(primers) into stats_output
    
    script:

    """

    java -jar -Xmx11g /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/bin/MergeAndCall.jar \
    DeepSeq \
    -inputMerged=${aligned_fastq} \
    -cutSites=${cutsites} \
    -outputStats=${sampleId}.stats \
    -primersEachEnd=${primers} \
    -primerMismatches=3 \
    -primersToCheck=${params.primersToCheck} \
    -requiredMatchingProp=${params.alignmentThreshold} \
    -requiredRemainingBases=${params.minimumReadLength} \
    -cutsiteWindow=5

    gzip ${sampleId}.stats
    """
}
/*
 * Generate a number of derivative files from the stats 
 */
process BreakdownFiles {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/04_breakdown_files"

    input:
    set sampleId, path(reference), path(targets), path(stats), path(cutsites), path(primers) from stats_output
    
    output:
    tuple sampleId, path(reference), path(targets), path("${sampleId}.perBase"), path("${sampleId}.topReadEvents"), path("${sampleId}.topReadEventsNew"), path("${sampleId}.topReadCounts"), path("${sampleId}.allReadCounts"), path("${sampleId}_event_length_histrogram.txt"), path(cutsites) into breakdown_files
    
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

process CreateSamplePlot {
    publishDir "$results_path/05_sample_plot"

    input:
    tuple sample, path(reference), path(targets), path(perBase), path(topReadEvents), path(topReadEventsNew), path(topReadCounts), path(allReadCounts), path(event_lengths), path(cutsites) from breakdown_files

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

