# DNA amplicon processing pipelines

This folder contains Nextflow pipelines for processing DNA amplicon lineage recorder data. Pipelines are split by their CRISPR type. 
There are a couple of configuration files that have a sharred format across all pipelines in this directory:

- a samplesheet. It __must__ have columns with names sample, reference, targets, read1, read2. This format is explained more below
- nextflow.config file. This file contains configuration details that are shared across samples within the run.