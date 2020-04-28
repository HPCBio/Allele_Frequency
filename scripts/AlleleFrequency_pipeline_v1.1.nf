#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
               (Fluidigm) Allele Frequency pipeline
========================================================================================
*/

// version
version                     = 1.1

// Credit to Phil Ewels for this segment to generate the help message at the command-line

def helpMessage() {
    log.info"""
    =========================================
     (Fluidigm) Allele Frequency pipeline v${version} 
    =========================================

    Usage:
    The typical command for running the pipeline is as follows:
    
    nextflow run -c <conf>  <AlleleFrequency_pipeline.nf>
    
    where conf is the configuration file for a particular experiment
    
    To override existing values from the command line, please use these parameters:

        Required options
	--genome                Path to fasta file of reference genome to be used for analysis (must be surrounded by quotes)                          
	--outdir                Path to output directory where the results will be saved. (must be surrounded by quotes)
	--reads                 Path to input data (must be surrounded by quotes)
        --ampliconName          Name of the amplicon  (must be surrounded by quotes)
        
	Read preparation options:
	--singleEnd             options: true|false. true for single aka SE reads; false for paired aka PE reads. Default: false.
	--mergeReads            options: true|false. Merge paired reads after trimming. Default: false.
	--guess_adapter         options: true|false. auto-detect adapter from input file. Default: false.
        --trim_lowC             options: true|false. Filter low complexity reads. Default: false.
        --trim_polyG            options: true|false. Trim polyGs from the end of the read. Default: false.
        --forward_adaptor       foward adapter sequence to be clipped off (must be surrounded by quotes). cannot be combined with guess_adapter
	--reverse_adaptor       reverse adapter sequence to be clipped off (must be surrounded by quotes). cannot be combined with guess_adapter
	--min_read_length       minimum length of read to be kept after trimming for downstream analysis. Default: 70
	--max_read_length       maximum length of read to be kept after trimming for downstream analysis. Default: 230
        --read_length           read length before trimming. Default: 250
	--min_BQ                minimum base quality. Default: 20
     
        Alignment and BAM filtering options
        --bwaparams             BWA alignment parameters enclosed in double quotes
        --min_MQ                minimum mapping quality. positive integer. mpileup option. Default: 0 
        --min_coverage          minimum coverage. positive integer. VarScan option. Default: 3
        --p_value               p value. VarScan  option. Default: 0.01
        --min_var_freq          minimum var frequency. VarScan  option. Default: 0.01
        
    """.stripIndent()
}



/*
 * SET UP CONFIGURATION VARIABLES
 */


// Override these values from the command line

// Inputs
params.reads                = false
params.genome               = false
params.outputDir            = false
params.outdir               = false
params.ampliconName         = false

// read prep options
params.singleEnd            = false              /* true for SE read. false for PE reads */
params.guess_adapter        = false              /* if true, fatsp will guess adaptors from the first 1M reads */
params.mergeReads           = false              /* only available for PE reads. Reads must overlap */
params.minOverlap           = '50'               /* minimum overlap for merging PE reads */
params.maxLength            = '121'              /* maximum length after stitching */
params.min_BQ               = '20'               /* minimum base quality */
params.forward_adaptor      = false              /* sequence of forward adaptor surrounded by quotes */
params.reverse_adaptor      = false              /* sequence of reverse adaptor surrounded by quotes */
params.min_read_length      = '70'               /* minimum read length after trimming */
params.max_read_length      = '230'              /* maximum read length after trimming */
params.read_length          = '250'              /* read length before trimming */
params.trim_lowC            = true               /* Filter low complexity reads */
params.trim_polyG           = false              /* Trim polyGs from  read */
params.RGLB                 = "fluidigm-nano"    /* LB field of RG line */
params.RGPL                 = "illumina"         /* PL field of RG line */
params.RGCN                 = "CBC"              /* CN field of RG line */
params.RGPU                 = "${params.ampliconName}"

// alignment mpileup VarScan options
params.bwaOptions           = " "                
params.bwaIdxAlgorithm      = " "                /* Options: -a is | -a bwtsw */
params.min_MQ               = '0'                /* deprecated */
params.p_value              = '0.01'
params.min_var_freq         = '0.01'
params.min_coverage         = '3'
params.TargetRegion         = false
params.segment              = false
params.downsample           = '100000'

// other options
params.email                = false


// cluster-specific variables
params.executor             = 'slurm'
params.Queue                = 'normal'
params.mem                  = '10'
params.cpus                 = '1'
params.trimThreads          = '4'
params.trimMem              = '20'
params.alignThreads         = '10'
params.alignMem             = '250'
params.vcallThreads         = '2'
params.vcallMem             = '50'

// Software stack
params.multiqcMod           = 'MultiQC/1.7-IGB-gcc-4.9.4-Python-3.6.1'
params.fastpMod             = 'fastp/0.19.5-IGB-gcc-4.9.4'
params.samtoolsMod          = 'SAMtools/1.9-IGB-gcc-4.9.4'
params.picardMod            = 'picard/2.10.1-Java-1.8.0_152'
params.bwaMod               = 'BWA/0.7.17-IGB-gcc-4.9.4'
params.VarScanMod           = 'VarScan/2.3.9-Java-1.8.0_152'
params.scriptdir            = "/home/groups/hpcbio/projects/anakk/fluidigm-Oct-2019/src"
params.seqtkMod             = 'seqtk/1.2-IGB-gcc-4.9.4'
params.PEARmod              = 'PEAR/0.9.8-IGB-gcc-4.9.4'
params.bamreadcountMod      = 'bam-readcount/0.8.0-IGB-gcc-4.9.4'
params.perlMod              = 'Perl/5.24.1-IGB-gcc-4.9.4'

// Sanity check
if (!params.ampliconName)                   exit 1, "Must set --ampliconName"
if (!params.outdir)                         exit 1, "Must set --outdir with path to results folder"
if (!params.reads)                          exit 1, "Must set --reads with a valid regexp. SE or PE Illumina short reads are expected."
if (!params.genome)                         exit 1, "Must set --genome with filename of reference genome to be used for analysis"
if (params.mergeReads && params.singleEnd)  exit 1, "--mergeReads and --singleEnd cannot be used together"
if (!params.guess_adapter && !params.singleEnd && !params.forward_adaptor &&  !params.reverse_adaptor)  exit 1, "Must set --forward_adaptor and --reverse_adaptor "
if (!params.guess_adapter &&  params.singleEnd && !params.forward_adaptor)                              exit 1, "Must set --forward_adaptor"
if (!params.guess_adapter && !params.singleEnd && !params.forward_adaptor &&  params.reverse_adaptor)   exit 1, "Must set --forward_adaptor"
if (!params.guess_adapter && !params.singleEnd &&  params.forward_adaptor && !params.reverse_adaptor)   exit 1, "Must set --reverse_adaptor"

// variables for fastp
adapterOptionsPE          = params.guess_adapter ? " --detect_adapter_for_pe " : " --detect_adapter_for_pe --adapter_sequence=${params.forward_adaptor}  --adapter_sequence_r2=${params.reverse_adaptor}  "
adapterOptionsSE          = params.guess_adapter ? " " : " --adapter_sequence=${params.forward_adaptor} "
trimOptions2              = params.trim_lowC     ? ' --low_complexity_filter ' : ' ' 
trimOptions3              = params.trim_polyG    ? ' --trim_poly_g ' : ' --disable_trim_poly_g ' 
params.trimOptions        = " -5 -3 --cut_mean_quality ${params.min_BQ} --length_required ${params.min_read_length} "

// variables for indices
params.genomesDir           = "/home/groups/hpcbio/projects/anakk/fluidigm-Oct-2019/data/nxt-generated"
genomeFile                  = file(params.genome)
genomePrefix                = genomeFile.getBaseName()
genomeStore                 = "${params.genomesDir}/${genomePrefix}"
bedFile                     = file(params.TargetRegion)

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Header log info

runInfo = """

------------------------------------------------------------------------- 
(Fluidigm) Allele Frequency pipeline v${version} 
------------------------------------------------------------------------- 
Pipeline ver               : ${version}

REQUIRED ARGUMENTS
Amplicon                   : ${params.ampliconName}
Reads                      : ${params.reads}
Genome                     : ${params.genome}
Output dir                 : ${params.outdir}

TOOLS
trim tool                  : ${params.fastpMod }
align tool                 : ${params.bwaMod }
samtools                   : ${params.samtoolsMod }
VarScan tool               : ${params.VarScanMod}
bam readcounts tools       : ${params.bamreadcountMod}

READ PREPARATION OPTIONS
read length                : ${params.read_length}
min read len after trim    : ${params.min_read_length}
if(!params.guess_adapter) {
guess.adaptor              : false
forward adaptor            : ${params.forward_adaptor}
reverse adaptor            : ${params.reverse_adaptor}
}
if(params.guess_adapter) {
guess.adaptor              : true
}
if(params.singleEnd) {
Read type                  : SE
merge PE reads             : false
adaptor trimming options   : ${adapterOptionsSE}
}
if(!params.singleEnd) {
Read type                  : PE
Merge PE reads             : ${params.mergeReads}
minimum overlap            : ${params.minOverlap}
adaptor trimming options   : ${adapterOptionsPE}
}
other filter/trimming options  : ${params.trimOptions} ${trimOptions2} ${trimOptions3}

ALIGN MPILEUP VARSCAN OPTIONS
BWA options                : ${params.bwaOptions }
minimum base qual          : ${params.min_MQ}
min coverage               : ${params.min_coverage}
p value                    : ${params.p_value}
min_var_freq               : ${params.min_var_freq}
downsampling               : ${params.downsample}


OTHER OPTIONS
Current home               : $HOME
Current user               : $USER
Current path               : $PWD
Script dir                 : ${params.scriptdir}
Working dir                : $workDir
------------------------------------------------------------------------- 
"""

println(runInfo)

/*
 *
 * Phase 1: Genome preparation
 *
 * Index genome if needed. Create fai, bed, dict files
 *
 */

 
process Prepare_Genome {
    tag                    { gf }
    executor               params.executor
    cpus                   params.cpus
    queue                  params.Queue
    memory                 "${params.mem} GB"
    module                 params.samtoolsMod,params.picardMod 
    storeDir               genomeStore
    
    input:
    file gf from genomeFile

    output:
    file "*.bed" into refBED
    file "*" into ref4VarScan
    
    script:
    """
    samtools faidx $gf
    cat ${gf}.fai | awk -v OFS='\t' '{chr = \$1; len = \$2; print chr, 1, len }' > ${gf}.bed
    java -jar \${EBROOTPICARD}/picard.jar CreateSequenceDictionary R=${gf} O=${gf}.dict 
    """
}


process BWA_Index_Genome {
    tag                    { gf }
    executor               params.executor
    cpus                   params.cpus
    queue                  params.Queue
    memory                 "${params.mem} GB"
    module                 params.bwaMod
    storeDir               genomeStore
    validExitStatus        0

    input:
    file gf from genomeFile

    output:
    file "*.{sa,amb,ann,pac,bwt}" into bwaIndexRef

    script:
    """
    bwa index -p ${gf.getBaseName()} ${params.bwaIdxAlgorithm} ${gf} 
    """
}



/*
 *
 * Phase 2: Read preprocessing/QC
 *
 * fastp  performs QC_PreTrim, Trim+Filter, QC-PostTrim. 
 * Optional step : merge PE reads
 *
 */


Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { reads2Trim1; reads2Trim2 }

if(params.mergeReads && !params.singleEnd) {

        /* Read processing with fastp  and merging PE reads*/     
   	
	process readPrep_fastp_and_merge {
            tag                    "reads: $name"
	    executor               params.executor
	    cpus                   params.trimThreads
	    queue                  params.Queue
	    memory                 "${params.trimMem} GB"
	    module                 params.fastpMod,params.PEARmod 
	    publishDir             "${params.outdir}/tmpResult-fastp_readPrep", mode: 'copy', overwrite: true	    
	    validExitStatus        0,1
	    
	    input:
	    set name, file(reads) from reads2Trim1

	    output:
	    set val(name), file('*.merged_and_R1.trimmed.fq') optional true into trim2bwaResults
	    set val(name), file('*.json') optional true into multiqcI    
	    file '*'
	    	    
	    script:   
	    """
		fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 "${name}.R1.trimmed.fq" --out2 "${name}.R2.trimmed.fq"  \
		${adapterOptionsPE}  ${params.trimOptions} ${trimOptions2} ${trimOptions3} \
		 --thread ${task.cpus} --html "${name}"_fastp.html --json "${name}"_fastp.json
		 
		pear -f "${name}.R1.trimmed.fq" -r "${name}.R2.trimmed.fq" -o ${name} \
		-v ${params.minOverlap} -n ${params.min_read_length} -j ${task.cpus}

		cat "${name}.assembled.fastq" "${name}.unassembled.forward.fastq" > "${name}.merged_and_R1.trimmed.fq"
	    """
    
	}  /* end process */
} /* end if */

if(!params.mergeReads) {

        /* Read processing with fastp  and no merging of PE reads*/     
   	
	process readPrep_fastp_noMerge {
            tag                    "reads: $name"
	    executor               params.executor
	    cpus                   params.trimThreads
	    queue                  params.Queue
	    memory                 "${params.trimMem} GB"
	    module                 params.fastpMod
	    publishDir             "${params.outdir}/tmpResult-fastp_readPrep", mode: 'copy', overwrite: true	    
	    validExitStatus        0,1
	    
	    input:
	    set name, file(reads) from reads2Trim2

	    output:
	    set val(name), file('*.trimmed.fq') optional true into trim2bwaResults
	    set val(name), file('*.json') optional true into  multiqcII    
	    file '*'
	    	    
	    script: 	    
	    if(params.singleEnd){
		"""
		fastp --in1 ${reads[0]} --out1 "${name}.R1.trimmed.fq" \
		${adapterOptionsSE} ${params.trimOptions} ${trimOptions2} ${trimOptions3} \
		--thread ${task.cpus} -w ${task.cpus} --html "${name}"_fastp.html --json "${name}"_fastp.json

		"""
	    } else {
		"""
		fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 "${name}.R1.trimmed.fq" --out2 "${name}.R2.trimmed.fq"  \
		${adapterOptionsPE}  ${params.trimOptions} ${trimOptions2} ${trimOptions3} \
		 --thread ${task.cpus} --html "${name}"_fastp.html --json "${name}"_fastp.json
		"""
	    }  /* end if  singleEnd*/

	}  /* end process */

} /* end if */


/*
 *
 * Phase 3: BAM preparation
 *
 * BWA alignment, deduplication and sort.
 *
 */


process BWA_aln_dedup {
    tag                    { name }
    executor               params.executor
    cpus                   params.alignThreads
    queue                  params.Queue
    memory                 "${params.alignMem} GB"
    errorStrategy          'finish'
    module                 params.bwaMod,params.samtoolsMod,params.picardMod,params.seqtkMod
    publishDir             "${params.outdir}/BWA-align-n-sort", mode: 'copy'

    input:
    set name, file(trimRead) from trim2bwaResults
    file bwaIndex from bwaIndexRef

    output:
    set val(name), file('*.ready.bam')  optional true into ReadyBAM
    file "*.bai" optional true into ReadyBAI
    file "*.stat" optional true into bwaAlnstat    
    file "*"

    script:
    RGLINE="@RG\\tID:${name}\\tSM:${name}\\tLB:${params.RGLB}\\tPL:${params.RGPL}\\tPU:${params.RGPU}\\tCN:${params.RGCN}"
    if(!params.mergeReads && !params.singleEnd){    
    """
    echo we have TWO reads files per sample ${name}
    
    echo downsampling to ${params.downsample}...
    
    seqtk sample -s10017 ${trimRead[0]} ${params.downsample} > ${name}.downsample.R1.fq
    seqtk sample -s10017 ${trimRead[1]} ${params.downsample} > ${name}.downsample.R2.fq

    echo bwa aligning...
    
    bwa mem ${params.bwaOptions} -t ${params.alignThreads} -R \"${RGLINE}\" ${genomePrefix} \
    ${name}.downsample.R1.fq ${name}.downsample.R2.fq | samtools sort -O bam -@  ${params.alignThreads} - > ${name}.ready.bam    

    samtools index ${name}.ready.bam
    samtools flagstat ${name}.ready.bam > ${name}.ready.bam.stat 
    """
    } else {    
    """    
    echo we have ONE read file per sample ${name}
   
    echo downsampling to ${params.downsample}...
    
    seqtk sample -s10017 ${trimRead[0]} ${params.downsample} > ${name}.downsample.R1.fq

    echo bwa aligning...
    
    bwa mem ${params.bwaOptions}  -t ${params.alignThreads} -R \"${rgheader}\"  ${genomePrefix} \
    ${name}.downsample.R1.fq | samtools sort -O bam -@  ${params.alignThreads} - > ${name}.ready.bam    

    samtools index ${name}.ready.bam
    samtools flagstat ${name}.ready.bam > ${name}.ready.bam.stat 
    """
    }
} /* end process */


/*
 *
 * Phase 4: Allele Freq and Variant discovery
 *
 *
 */

process readcounts_and_variants {
    tag { id }
    executor               params.executor
    cpus                   params.vcallThreads
    queue                  params.Queue
    memory                 "${params.vcallMem} GB"
    module                 params.VarScanMod,params.samtoolsMod,params.bamreadcountMod,params.perlMod
    publishDir             "${params.outdir}/readcounts", mode: 'copy', overwrite: true
    validExitStatus        0,1
    errorStrategy          'finish'

    input:
    file refgenome from genomeFile
    file refidx from ref4VarScan
    file region from bedFile
    
    set val(id), file(readybam) from ReadyBAM
    file readybai from ReadyBAI

    output:
    set val(id), file("*.vcf")  optional true into ANNVariantsVCF
    file '*' 

    script:

    """
    samtools mpileup -l ${region} -o ${id}.pileup -B --max-depth 0 -f ${refgenome} ${readybam}
    
    java -jar \${EBROOTVARSCAN}/VarScan.v2.3.9.jar readcounts   ${id}.pileup --output-file  ${id}.varScan.readcounts.txt
    
    java -jar \${EBROOTVARSCAN}/VarScan.v2.3.9.jar mpileup2cns  ${id}.pileup \
    --p-value ${params.p_value} --min-var-freq ${params.min_var_freq} --min-coverage ${params.min_coverage } \
    --variants --output-vcf > ${id}.rawVariantsOnly.vcf
    
    bam-readcount -b ${params.min_BQ} -q ${params.min_MQ} -d 100000 -f ${refgenome} ${readybam} ${params.segment} > ${id}.BAM.readcounts.txt
    perl ${params.scriptdir}/summarize_bamreadcounts.pl -in ${id}.BAM.readcounts.txt -out ${id}.BAM.summary.readcounts.txt
    
    """
}  /* end process */

if(params.mergeReads) {
	process MultiQC_merge {
	    executor               params.executor
	    cpus                   params.cpus
	    queue                  params.Queue
	    memory                 "${params.mem} GB"
	    module                 params.multiqcMod
	    publishDir             "${params.outdir}/MultiQC", mode: 'copy', overwrite: true

	    input:
	    file('/readprep_merge/*')   from multiqcI.collect()
	    file('/BWAaln/*')           from bwaAlnstat.flatten().toList()

	    output:
	    file "*_report.html" into multiqc_report_post
	    file "*_data"

	    script:
	    """
	    multiqc -d -f  .
	    """
	} /* end process */
} else {
	process MultiQC_noMerge {
	    executor               params.executor
	    cpus                   params.cpus
	    queue                  params.Queue
	    memory                 "${params.mem} GB"
	    module                 params.multiqcMod
	    publishDir             "${params.outdir}/MultiQC", mode: 'copy', overwrite: true

	    input:
	    file('/readprep_nomerge/*') from multiqcII.collect()
	    file('/BWAaln/*')           from bwaAlnstat.flatten().toList()

	    output:
	    file "*_report.html" into multiqc_report_post
	    file "*_data"

	    script:
	    """
	    multiqc -d -f  .
	    """
	} /* end process */
}
	
workflow.onComplete {

      def subject = "[Allele Frequency pipeline] Successful: $workflow.runName"
      if(!workflow.success){
          subject = "[Allele Frequency pipeline] FAILED: $workflow.runName"
      }
      
    finalLog = """
-------------------------------------------------------------------------    
[Allele Frequency pipeline] execution summary
-------------------------------------------------------------------------
Completed at : ${workflow.complete}
Duration     : ${workflow.duration}
Success      : ${workflow.success}
workDir      : ${workflow.workDir}
exit status  : ${workflow.exitStatus}
Error report : ${workflow.errorReport ?: '-'}
-------------------------------------------------------------------------

"""

    ['mail', '-s', subject, params.email ].execute() << runInfo
    
    log.info "[Allele Frequency pipeline] COMPLETED. Sent summary e-mail to $params.email"

}

