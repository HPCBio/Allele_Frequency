# Description
Allele frequency of codons in specific targeted regions of the mouse genome.


# Dependencies

This program expects the following tools/languages to be installed as modules and be available in your path:

- Nextflow        tested with version 20.01.0   ( download page https://github.com/nextflow-io/nextflow/releases) 
- MultiQC         tested with version 1.7  ( https://multiqc.info/ ) 
- fastp           tested with version 0.19.5 ( download page https://github.com/OpenGene/fastp )
- SAMtools        tested with version 1.9 ( download page http://www.htslib.org/ )
- picard          tested with version 2.10.1 ( download page https://broadinstitute.github.io/picard/ )
- BWA             tested with version 0.7.17 ( download page http://bio-bwa.sourceforge.net/ )
- seqtk           tested with version 1.2 ( download page https://github.com/lh3/seqtk )
- bam-readcount   tested with version 0.8.0 ( download page https://github.com/genome/bam-readcount )
- RStudio         ( download page https://rstudio.com/products/rstudio/download/ )


# Installation instructions

- Install all dependencies first. You may need to have root access to install some of these tools/languages on a cluster.
- Do not forget to launch the 'hello world' nextflow pipeline (as per https://www.nextflow.io/) to make sure it works fine.
- Clone this repo to a local unix-based cumputing platform like a cluster.

# Data preparation

## Input files

The input files should be placed in the folder <b>data/raw-seq/</b> and arranged as described here.

There should be one subfolder for each amplicon. The subfolder names should match exactly the amplicon names of the four regions as follows:

- BRAF584KF_BRAF584KR
- EGFRF_EGFRR
- BRAF637VF_BRAF637VR
- HRAS61QF_HRAS61QR

This program expects input files to come from  a Illumina platform using a fluidigm library preparation protocol. See this link for more information:  https://www.fluidigm.com/search?query=DNA%20sequencing

These files should be paired-ended fastq files; that is, there should be a forward read with extension R1.fastq and a reverse read with extension R2.fastq

We know that there are four targeted regions; moreover, there could be several samples per region. Your files must be demultiplexed by amplicon and by sample; moreover the filenames should conform to this pattern:

<pre>
 AmpliconName - SampleName _ * _R[1|2].fastq
</pre>


### Examples of filenames:

- Amplicon name: BRAF584KF_BRAF584KR
- Sample name: GIKO_N811
- Filenames of the two paired-end fastq files:

<pre>
BRAF584KF_BRAF584KR-GIKO_N811_CAGTCTACAT_R1.fastq
BRAF584KF_BRAF584KR-GIKO_N811_CAGTCTACAT_R2.fastq
</pre>

- Amplicon name: EGFRF_EGFRR
- Sample name: GIKO_N811
- Filenames of the two paired-end fastq files:

<pre>
EGFRF_EGFRR-GIKO_N811_CAGTCTACAT_R1.fastq
EGFRF_EGFRR-GIKO_N811_CAGTCTACAT_R2.fastq
</pre>

Finally,  all demultiplexed files for a particular amplicon should be placed inside the subfolder with the corresponding amplicon name.  In other words, all demultiplexed files starting with <b> EGFRF_EGFRR </b> should go in the subfolder named <b> EGFRF_EGFRR </b>. Misplaced files will simply be skipped/ignored by the pipeline.

## The genome files

The genome files should be placed in the folder <b>data/genome/</b>

These are the instructions for downloading and decompressing the genome files:

<pre>

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz
gunzip GRCm38.primary_assembly.genome.fa.gz
gunzip gencode.vM23.primary_assembly.annotation.gtf.gz

</pre>

Also copy or move these files that came with the repo in the data folder to the same location where the genomes files are

- Amplicon_Seqs.fasta
- chr11.bed
- chr6.bed
- chr7.bed



# How to run the pipeline

The pipeline should be run for each amplicon separately. It needs a different configuration file for each amplicon.

The configuration file is just a text file containing one parameter-value pair per line. In other words, varriables set to specific values that the pipeline expects to configure and use during the analysis. Since there are so many of them, it is preferable to put them in a configuration file than to type them at the command line.

We included examples of configuration files in the <b> conf/ </b> folder that came with this repo.

### How to run the pipeline for the EGFRF_EGFRR amplicon

To run the pipeline type this command:

<pre>
nextflow run -c EGFR_codon254F_new.conf  AlleleFrequency_pipeline.nf
</pre>

### Outputs of the pipeline

After the pipeline finishes execution successfully, you should see several folders inside the output path that you specified inside the configuration file.

- <b>tmpResult-fastp_readPrep</b> This folder contains the results of the first step, read preparation with fastp
- <b>BWA-align-n-sort</b> This folder contains the reults of the second step, alignment to genome with BWA
- <b>readcounts</b> This folder contains the results of the third step, calclate allele frequecies with read-count tool
- <b>MultiQC</b> This folder contains the results of the last step, read tracking with MultiQC

In order to generate plots with the data produced by this pipeline, please read instructions in the next section.


# How to run the R script to generate the plots

The input files for this R script are inside the <b> results/ampliconx/readcounts/ </b>.

This folder has many intermediary results; we just need those ending in <b> BAM.summary.readcounts.txt </b> 

- To keep sanity, simply make a new folder on your laptop <b> results/ampliconx/plots/ </b> and copy those file ending in <b> BAM.summary.readcounts.txt </b> from <b> results/ampliconx/readcounts/ </b> to <b> results4ampliconx/plots/ </b>.  You can do it from the command line with the scp command as shown below:

<pre>

mkdir -p mylaptop/results/ampliconx/plots/
scp   results/ampliconx/readcounts/*BAM.summary.readcounts.txt   userid@domain:mylaptop/results/ampliconx/plot/

</pre>

Or you can use your favorite file transfer tool.

- Next, copy the R script that came with this repo from <b> scripts/ </b> to this location

<pre>

scp   scripts/import_and_plot_from_full_counts.R   userid@domain:mylaptop/results/ampliconx/plot/

</pre>

- Start RStudio on your laptop and load the script scripts/import_and_plot_from_full_counts.R
- Edit the script with the correct path to the default folder so that it will find the input files
- Hit the run button to execute the script
- It should generate a PDF file with the plot

That's all

