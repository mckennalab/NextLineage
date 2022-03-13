# Pipeline for processing xxCas12a molecular recording data

This pipeline processes short-read molecular recording sequences, generating common output files and plots. 

## Setup 

You'll need to have Nextflow setup. Direction are on [Nexflow website](https://www.nextflow.io/). We assume you'll have a nextflow file in the directory you're running from, though you can change the run script to match where it's located on your setup.

Running the pipeline is a simple command. In the example directory the [shell script](https://github.com/mckennalab/NextLineage/blob/main/dna_amplicon/cas12a/example_setup/run_pipeline.sh) looks something like:

```
./nextflow <path_to_github_checkout>/dna_amplicon/cas12a/cas12a_amplicon.nf --samplesheet FEX_CPF1_tearsheet.txt -resume
```

Here we're calling nextflow on our pipeline file (cas12a_amplicon.nf) using a specified samplesheet. We're also asking Nextflow to resume any previous run; if we've already run the pipeline it should pickup where we left off. This is all hopefully straightfoward. There are two main configuration files, the samplesheet listed above, as well as the nextflow.config file. 

### nextflow.config

The nextflow configuration file contains parameters that are common to the run as a whole. The example file looks something like:

```
params.primer5 = "CTAGCAGATCACCGTAAGGA" 
params.primer3 = "CACTGTACTCTACGCGACTC"
params.aligner = "needle"
params.primerMismatchCount = 3
params.alignmentThreshold = 0.85
params.minimumReadLength = 50
params.convert_unknown_to_none = "FALSE"
params.umiStart = -1
params.umiLength = 12
params.minUMIReadCount = 1
params.primersToCheck = "BOTH"
params.crispr = 'cas12a'
process.executor = 'slurm'
```

We can go through each parameter below. All parameters need to be prefixed with ```params.```, which tells Nextflow to make them visible to all pipeline steps. 

- ```primer5```: This is the 5' primer that's expected to start the first read. Multiple steps in the pipeline will assess reads for this sequence, given the mismatch parameters, and remove reads missing this sequence. Why is this used? Occasionally we'll accidentally pull out off-target locations in the genome, and the two primers are used to filter out these sequences.
- ```primer3```: the same as primer5, this sequence is expected to start read 2, but this is **oriented 5' to 3' on the global sequence**. So your second read sequence will start with the reverse compliment of this sequence. This was done to make cut/paste a little easier, but it does lead to some confusion. A diagram:

```
setup:
-------------------- reference ----------------
----> actual primer1    /      primer2 <-------
-------------> read1      read2 <--------------

----> primer5 above / primer3 above      ----->
```

- ```aligner```: two options, either ```needle```, a standard Needleman/Wunsch affine gap alignment, or the ```convex``` aligner that uses a convex gap function to assign longer deletions more accurately. For amplicon sequencing ```needle``` should be your default, as convex is significantly more expensive computationally. 

- ```primerMismatchCount:```: how many mismatches we allow when matching ```primer5``` and ```primer3``` to the reads. 

- ```alignmentThreshold```: after we've aligned reads, how many non-gap sequences must be matches vs mismatches. Sometimes PCR will pickup on endogenous sequences, converting their ends to the primer sequences (and won't be caught with our primer5/primer3 filter above). These sequences will often be poor alignments to our amplicon reference, and this check throws them out if less than XX % of bases match the reference. This check ignore indels (for obvious reasons). 

- ```minimumReadLength```: after all our cleaning steps, how many bases of the combined reads has to be left. This avoids aligning primer-dimer fragments or other problematic fragments.

- ```convert_unknown_to_none```: We try to merge the forward and reverse reads to fully cover the target sites. If we can't, we continue with the paired reads. Sometimes these reads won't cover a target site (both ends are too short, either due to quality filtering or experimental design). In this 'uncovered' case the site will be called _UNKNOWN_. Recorders containing _UNKNOWN_ are excluded from many of the visualizations. This parameter changes our perspective, converting these _UNKNOWN_ to _NONE_ calls, which indicates you'll assume they're wild-type target sequences. Use with care, this is really only intended for diagnostics.

- ```umiStart```: where the UMI sequence (if included) starts within the reads. Read 1 starts at position 0 and offsets are positive numbers. If the UMI is on read 2, you use negitive numbers. Read 2 offsets also start at -1, so if the UMI is at the 9th base in read 1 its position would be 8, for read 2 this would be -9. 

- ```umiLength```: the length of the UMI sequence. **Setting this parameter means all samples in the run have UMIs of this length** (the samefor ```umiStart above```). For now if you have multiple UMI lengths per run, split it into multiple runs. 

- ```minUMIReadCount```: when using UMIs, how many reads are required to call the UMI successfully captured (and for it to be used in downstream analysis). 

- ```primersToCheck```: for primer5 and primer3, do we want to check for both (_BOTH_), not require either (_NONE_), or just require the primer5 (_FORWARD_) or just the primer3 (_REVERSE_). Generally we set this to _BOTH_, though in some cases (like scRNA-seq) we can only rely on one of the primers being present (_FORWARD_). 

- ```crispr```: which CRISPR enzyme we're using. In this case, cas12a. 

We also have one non ```params.``` setting here, ```process.executor```. This setting can be left out, but if you have a local compute farm this allows Nextflow to run on a distributed computing system. For instance at Dartmouth, if you login to [Discovery](https://rc.dartmouth.edu/index.php/discovery-overview/), you can set this parameter to ```slurm``` and Nextflow will run each job in parallel, which immensely speeds things up.

### samplesheet

Our sample sheet describes the files associated with each of our samples. This is a tab-separated file with the following header and one sample per line:

```sample	reference	targets	read1	read2```

Columns:

- ```sample```: the sample name. **No spaces**, but try to make your life easier and keep these informative, most files will have this prepended to their name

- ```reference```: the reference file of just the amplicon region in fasta format. The reference can include sequencing adapters or be larger segments from a plasmid, etc, but the longer the sequence the more effort the aligner will have to put in. 

- ```targets```: a list of the CRISPR target sequences, with PAMs, one per line. This file should have a header with the name ```sites```. 

- ```read1``` and ```read2```: the read one (forward reads), and read 2 (reverse reads) as compressed fastq files (like something.fq.gz or something.fastq.gz). **Be careful here**: we assume read 1 should align to the reference you provide in the forward orientation, and read 2 the reverse. Depending on your sequencing configuration and reference file, this could be flipped. The pipeline doesn't check names, so if you need to put what Illumina sequenced as the second read into the ```read1``` column (and the first illumina read into ```read2``` column) it's fine, as long as they're oriented so that ```read1``` will be seen before ```read2``` on the reference.
