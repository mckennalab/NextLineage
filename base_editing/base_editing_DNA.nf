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


Channel.fromPath( params.samplesheet )
        .splitCsv(header: true, sep: '\t')
        .map{ tuple(it.sample, it.umi, it.umiLength, it.reference, it.targets, it.primer5, it.primer3, it.fastq1, it.fastq2)}
        .into{sample_table_cutsites; sample_table_filter; sample_table_alignment; sample_table_stats; sample_table_final_tables}

results_path = "results"
primers_extension = 5


/*
 * Create the cutSites location file and validate their primers
 */
process CreateCutSitesFile {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/01_cutSites"

    input:
    set sampleId, umi, umiLength, reference, targets, fwd_primer, rev_primer, fastq1, fastq2 from sample_table_cutsites

    output:
    set sampleId, "${sampleId}.fa.cutSites", "${sampleId}.fa", "${sampleId}.fa.primers" into reference_pack, reference_pack2, reference_pack3
    
    script:
    
    """
    cp ${reference} ${sampleId}.fa
    python /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/fasta_to_cutsites.py \
        --reference ${sampleId}.fa \
        --sites ${targets} \
        --forward_primer ${fwd_primer} \
        --reverse_primer ${rev_primer} \
        --CRISPR_type Cas9
    """
}

/*
 * merge paired reads on overlapping sequences using 
 */
process MergeReads {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/03_ng_trim"

    input:
    set sampleId, umi, umiLength, reference, targets, fwd, rev, read1, read2 from sample_table_filter
    
    output:
    set sampleId, "${sampleId}_merged.fastq.gz", "${sampleId}_single_1.fastq.gz", "${sampleId}_single_2.fastq.gz" into merged_reads
    
    script:
    
    """
    /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/tools/NGmerge-0.3/NGmerge \
    -1 ${read1} -2 ${read2} -o ${sampleId}_merged.fastq -f ${sampleId}_single
    """
}

/*
 * Zip together read pairs from the unmerged reads into a single, interleaved file
 */
process InterleaveUnmergedReads {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/04_interleaved"

    input:
    set sampleId, merged, read1, read2 from merged_reads
    
    output:
    set sampleId, merged, "${sampleId}_interleaved.fq.gz" into merged_and_interleaved_reads
    
    script:
    
    """
    scala /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/zip_two_read_files.scala \
    ${read1} ${read2} ${sampleId}_interleaved.fq _ _
    gzip ${sampleId}_interleaved.fq
    """
}

phased_ref_reads_for_alignment = sample_table_alignment.phase(merged_and_interleaved_reads)


/*
 * Zip together read-pairs from the unmerged reads into a single, interleaved file
 */
process AlignReads {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/05_alignment"

    input:
    val tuple_pack from phased_ref_reads_for_alignment
    
    output:
    set sampleId, "${sampleId}.merged.fasta.gz", "${sampleId}.interleaved.fasta.gz" into aligned_reads
    
    script:
    sampleId = tuple_pack.get(0).get(0)
    reference = tuple_pack.get(0).get(3)
    merged_reads = tuple_pack.get(1).get(1)
    interleaved_reads = tuple_pack.get(1).get(1)

    """
    cp ${merged_reads}  ./merged_reads.fastq.gz
    gunzip ./merged_reads.fastq.gz

    cp ${interleaved_reads}  ./interleaved_reads.fastq.gz
    gunzip ./interleaved_reads.fastq.gz

    scala /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/align_merged_unmerged_reads.scala \
    /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/resources/EDNAFULL.Ns_are_zero \
    ${params.aligner} \
    ./merged_reads.fastq \
    ./interleaved_reads.fastq \
    ${reference} \
    ${sampleId}.merged.fasta \
    ${sampleId}.interleaved.fasta \
    10 \
    0.5

    rm ./merged_reads.fastq
    rm ./interleaved_reads.fastq

    gzip ${sampleId}.merged.fasta
    gzip ${sampleId}.interleaved.fasta 

    """
}

aligned_reads_phased = reference_pack.phase(aligned_reads)

/*
 * Given the aligned reads, call the events over the recorders
 */
process CallEvents {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/06_event_calls"

    input:
    val tuple_pack from aligned_reads_phased
    
    output:
    set sampleId, "${sampleId}.stats.gz" into stats_file
    
    script:
    sampleId = tuple_pack.get(0).get(0)
    reference = tuple_pack.get(0).get(2)
    cutSites = tuple_pack.get(0).get(1)
    primers = tuple_pack.get(0).get(3)
    merged_alignments = tuple_pack.get(1).get(1)
    interleaved_alignments = tuple_pack.get(1).get(1)

    """
    cp ${merged_alignments} merged.fasta.gz
    gunzip merged.fasta.gz
    cp ${interleaved_alignments} interleaved.fasta.gz
    gunzip interleaved.fasta.gz

    java -jar /dartfs/rc/lab/M/McKennaLab/resources/GESTALT/SingleCellLineage/UMIMerge/target/scala-2.12/MergeAndCall.jar \
    DeepSeq \
    -inputFileUnmerged=interleaved.fasta \
    -inputMerged=merged.fasta \
    -cutSites=${cutSites} \
    -outputStats=${sampleId}.stats \
    -primersEachEnd=${primers} \
    -sample=${sampleId} \
    -primerMismatches=3 \
    -primersToCheck=FORWARD \
    -requiredMatchingProp=0.85 \
    -requiredRemainingBases=75

    rm interleaved.fasta
    rm merged.fasta

    gzip ${sampleId}.stats
    """
}

stats_phased = reference_pack2.phase(stats_file)
stats_phased.into{stats_file_for_base_edits; stats_file_phased}

/*
 * Generate a number of derivative files from the stats 
 */
process BreakdownFiles {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/07_breakdown_files"

    input:
    val tuple_pack from stats_file_phased
    
    output:
    set sampleId, "${sampleId}.perBase", "${sampleId}.topReadEvents", "${sampleId}.topReadEventsNew", "${sampleId}.topReadCounts", "${sampleId}.allReadCounts" into breakdown_files
    
    script:
    sampleId = tuple_pack.get(0).get(0)
    reference = tuple_pack.get(0).get(2)
    cutSites = tuple_pack.get(0).get(1)
    primers = tuple_pack.get(0).get(3)
    stats = tuple_pack.get(1).get(1)

    """
    cp ${stats} ./sample.stats.gz
    gunzip ./sample.stats.gz

    scala /dartfs/rc/lab/M/McKennaLab/projects/nextflow_lineage/src/stats_to_javascript_tables.scala \
    sample.stats \
    ${sampleId}.perBase \
    ${sampleId}.topReadEvents \
    ${sampleId}.topReadCounts \
    ${sampleId}.allReadCounts \
    ${cutSites} \
    ${sampleId}.topReadEventsNew \
    ${params.convert_unknown_to_none} \
    ${reference}

    rm sample.stats
    """
}

/*
 * Generate a number of derivative files from the stats 
 */
process CallBaseEdits {
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/08_base_editing"

    input:
    val tuple_pack from stats_file_for_base_edits
    
    output:
    set sampleId, "${sampleId}.allele_counts" into basecalls
    
    script:
    sampleId = tuple_pack.get(0).get(0)
    reference = tuple_pack.get(0).get(2)
    cutSites = tuple_pack.get(0).get(1)
    primers = tuple_pack.get(0).get(3)
    stats = tuple_pack.get(1).get(1)

    """
    cut -f 24 <(zcat ${stats}) | sort | uniq -c | sort -r -n > ${sampleId}.allele_counts 
    """
}
