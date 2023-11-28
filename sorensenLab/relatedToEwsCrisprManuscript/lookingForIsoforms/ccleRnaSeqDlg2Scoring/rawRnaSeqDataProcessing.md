## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a the CCLE database, or DepMap. The data are all located in SRA: PRJNA523380.

### Description

I am interested in the RNAseq data.

I am going to use [sraDownloader](https://github.com/s-andrews/sradownloader) to get at these data.

I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).

I will use [STAR](https://github.com/alexdobin/STAR) for the alignment. I am following the instructions in the STAR [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

I will use [Salmon](https://github.com/COMBINE-lab/salmon) to do read quantification. There is a great info page for this tool [here](https://salmon.readthedocs.io/en/latest/).

I will then use [samtools](http://www.htslib.org/) to prepare the files before using [FeatureCounts](http://subread.sourceforge.net/featureCounts.html) (part of [Subread](http://subread.sourceforge.net/)) for quantification.

This is kind of a useful website for a general pipeline, [here](https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html), also [here](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html).

### Data pipeline

First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211112_ccleSequencingDataReprocessing
touch sraDataProcessingScript.sh
chmod +x sraDataProcessingScript.sh
touch snakefile
```

Edit the contents of the snakefile to include the text below. I use vim for this.

```python
"""
Author: Christopher Hughes
Affiliation: BCCRC
Aim: Workflow for RNA-Seq data
Date: 20211105
"""

###############################
#working directory
BASE_DIR = "/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211112_ccleSequencingDataReprocessing"


###############################
#locations of tools we will use
#BBDUK = "/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
BBDUK = "/projects/ptx_analysis/chughes/softwareTools/bbmap-38.90/bbduk.sh"
#SALMON = "/home/chughes/softwareTools/salmon-1.5.2/bin/salmon"
SALMON = "/projects/ptx_analysis/chughes/softwareTools/salmon-1.5.2/bin/salmon"
#HISAT2 = "/home/chughes/softwareTools/hisat2-2.2.1/hisat2"
#STAR = "/home/chughes/softwareTools/STAR-2.7.9a/bin/Linux_x86_64/STAR"
STAR = "/projects/ptx_analysis/chughes/softwareTools/STAR-2.7.9a/bin/Linux_x86_64/STAR"
#SAMTOOLS="/home/chughes/softwareTools/samtools-1.12/samtools"
SAMTOOLS="/gsc/software/linux-x86_64-centos7/samtools-1.14/bin/samtools"
#SAMBAMBA="/home/chughes/softwareTools/sambamba-0.8.1/sambamba"
#FEATURECOUNTS="/home/chughes/softwareTools/subread-2.0.3/bin/featureCounts"
FEATURECOUNTS="/projects/ptx_analysis/chughes/softwareTools/subread-2.0.3/bin/featureCounts"
#BAMCOVERAGE="/home/chughes/virtualPython368/bin/bamCoverage"
BAMCOVERAGE="/home/chughes/Virtual_Python383/bin/bamCoverage"


###############################
#locations of our index files
#DATABASE_DIR = "/home/chughes/databases/projectEwsDlg2"
DATABASE_DIR = "/projects/ptx_analysis/chughes/databases/projectEwsDlg2"
STARINDEX = DATABASE_DIR + "/starIndex"
SALMONINDEX = DATABASE_DIR + "/salmonIndex/salmon_index"
GTF = DATABASE_DIR + "/baseGenomeFiles/genome.gtf"
FASTA = DATABASE_DIR + "/baseGenomeFiles/genome.fa"


###############################
#get our list of sample files
SAMPLES, = glob_wildcards("raw/{smp}_1.fastq.gz")
print("There are " + str(len(SAMPLES)) + " total samples to be processed.")
for smp in SAMPLES:
  print("Sample " + smp + " will be processed")


###############################
#processing workflow
rule all:
    input:
      expand("results/{smp}.counts.txt", smp = SAMPLES),
      expand("results/{smp}.sorted.bam.bai", smp = SAMPLES),
      expand("results/{smp}.sorted.bw", smp = SAMPLES),
      expand("quants/{smp}/quant.sf", smp = SAMPLES)

rule bbduk:
  input:
      r1 = "raw/{smp}_1.fastq.gz",
      r2 = "raw/{smp}_2.fastq.gz"
  output:
      ro1 = "results/{smp}_1.clean.fastq.gz",
      ro2 = "results/{smp}_2.clean.fastq.gz"
  message:
      "Processing with BBDuk."
  shell:
      "{BBDUK} in1={input.r1} in2={input.r2} ref=adapters out1={output.ro1} out2={output.ro2} ktrim=r k=23 mink=11 hdist=1 tpe tbo"

rule star:
  input:
      r1 = "results/{smp}_1.clean.fastq.gz",
      r2 = "results/{smp}_2.clean.fastq.gz"
  output:
      "results/{smp}.sorted.bam"
  message:
      "Aligning with STAR."
  shell:
      "{STAR} --runThreadN 8 --genomeDir {STARINDEX} --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat --sjdbGTFfile {GTF} --outStd SAM | {SAMTOOLS} sort -o {output}"

rule bam_indexing:
  input:
      "results/{smp}.sorted.bam"
  output:
      "results/{smp}.sorted.bam.bai"
  message:
      "Indexing BAM with samtools."
  shell:
      "{SAMTOOLS} index {input}"

rule bam_coverage:
  input:
      r1 = "results/{smp}.sorted.bam",
      w2 = "results/{smp}.sorted.bam.bai"
  output:
      "results/{smp}.sorted.bw"
  message:
      "Calculating coverage with deeptools."
  shell:
      "{BAMCOVERAGE} -b {input.r1} -o {output} -p 8"

rule featurecounts:
  input:
      r1 = "results/{smp}.sorted.bam",
      w2 = "results/{smp}.sorted.bam.bai",
      w3 = "results/{smp}.sorted.bw"
  output:
      "results/{smp}.counts.txt"
  message:
      "Counting reads with featureCounts."
  shell:
      "{FEATURECOUNTS} -p --countReadPairs -t exon -g gene_id -a {GTF} -o {output} {input.r1}"

rule salmon:
  input:
      r1 = "results/{smp}_1.clean.fastq.gz",
      r2 = "results/{smp}_2.clean.fastq.gz",
      w3 = "results/{smp}.counts.txt"
  output:
      "quants/{smp}/quant.sf"
  params:
      dir = "quants/{smp}"
  message:
      "Quantifying with salmon."
  shell:
      "{SALMON} quant -i {SALMONINDEX} -l A -p 8 --gcBias --validateMappings -o {params.dir} -1 {input.r1} -2 {input.r2}"
```

Below is the shell script I will use to process these data with snakemake.

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
#sraDownloader="/home/chughes/softwareTools/sradownloader-3.8/sradownloader"
sraDownloader="/projects/ptx_analysis/chughes/softwareTools/sradownloader-3.8/sradownloader"
#sraCacheLocation="/mnt/Data/chughes/sratoolsRepository"
sraCacheLocation="/projects/ptx_results/Sequencing/sraCache"
#workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211112_ccleSequencingDataReprocessing"
workingDirectory="/projects/ptx_results/Sequencing/publishedStudies/sequencing20211112_ccleSequencingDataReprocessing"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results
eval mkdir quants

##loop over the accessions
#for i in SRR8616012 SRR8615497 SRR8616213 SRR8616214 SRR8615592 SRR8615679 SRR8615832 SRR8615859 SRR8615273 SRR8615499 SRR8615521
for i in SRR8615811 SRR8615470 SRR8616196 SRR8615602 SRR8615955 SRR8615610 SRR8615372 SRR8616175 SRR8616152 SRR8616003 SRR8616200 SRR8615908 SRR8615341 SRR8615679 SRR8615521 SRR8615614 SRR8615601 SRR8616193 SRR8616109 SRR8615875 SRR8615607 SRR8616102 SRR8616062 SRR8616107 SRR8616151 SRR8615771 SRR8615430 SRR8615384 SRR8615747 SRR8615916 SRR8615631 SRR8615331 SRR8615579 SRR8616169 SRR8616011 SRR8616001 SRR8615522 SRR8615961 SRR8615286 SRR8615616 SRR8615273 SRR8615832 SRR8616058 SRR8615675 SRR8615977 SRR8615426 SRR8616213 SRR8615527 SRR8616214 SRR8616146 SRR8615592 SRR8615998 SRR8616079 SRR8615625 SRR8615329 SRR8616161 SRR8615368 SRR8616189 SRR8615975 SRR8615807 SRR8615810 SRR8615812 SRR8616067 SRR8616078 SRR8615844 SRR8616019 SRR8615549 SRR8615605 SRR8616112 SRR8615727 SRR8615816 SRR8615703 SRR8615378 SRR8615712 SRR8615727 SRR8615423 SRR8615739 SRR8618307 SRR8615971 SRR8615727 SRR8615497 SRR8616096 SRR8615376 SRR8616005 SRR8618301 SRR8615422 SRR8615429 SRR8618305 SRR8615946 SRR8615697 SRR8616012
do
  printf "Downloading files associated with ${i}."
  eval ${sraDownloader} --outdir ${workingDirectory}/raw ${i}
  ##the file gets renamed upon download, but I just want it to have the SRR id and I can annotate it later
  eval mv ${workingDirectory}/raw/${i}*_1.fastq.gz ${workingDirectory}/raw/${i}_1.fastq.gz
  eval mv ${workingDirectory}/raw/${i}*_2.fastq.gz ${workingDirectory}/raw/${i}_2.fastq.gz
  #eval conda activate snakemake
  eval snakemake --cores 8 --latency-wait 300
  #eval conda deactivate
  eval rm ${workingDirectory}/raw/${i}*.fastq.gz
  eval rm ${sraCacheLocation}/sra/${i}*
  eval rm ${workingDirectory}/results/${i}*.clean.fastq.gz
  eval rm ${workingDirectory}/*.out
  eval rm ${workingDirectory}/*.tab
  eval rm -r ${workingDirectory}/*STAR*
  eval rm ${workingDirectory}/results/${i}*.bam
  eval rm ${workingDirectory}/results/${i}*.bai
done
```

Make another script to process through these accessions with megadepth to process the bigwig files.

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToEwsCrisprManuscript/lookingForIsoforms/ccleRnaSeqDlg2Scoring"
eval cd ${workingDirectory}

##loop over the accessions
for i in SRR8615811 SRR8615470 SRR8616196 SRR8615602 SRR8615955 SRR8615610 SRR8615372 SRR8616175 SRR8616152 SRR8616003 SRR8616200 SRR8615908 SRR8615341 SRR8615679 SRR8615521 SRR8615614 SRR8615601 SRR8616193 SRR8616109 SRR8615875 SRR8615607 SRR8616102 SRR8616062 SRR8616107 SRR8616151 SRR8615771 SRR8615430 SRR8615384 SRR8615747 SRR8615916 SRR8615631 SRR8615331 SRR8615579 SRR8616169 SRR8616011 SRR8616001 SRR8615522 SRR8615961 SRR8615286 SRR8615616 SRR8615273 SRR8615832 SRR8616058 SRR8615675 SRR8615977 SRR8615426 SRR8616213 SRR8615527 SRR8616214 SRR8616146 SRR8615592 SRR8615998 SRR8616079 SRR8615625 SRR8615329 SRR8616161 SRR8615368 SRR8616189 SRR8615975 SRR8615807 SRR8615810 SRR8615812 SRR8616067 SRR8616078 SRR8615844 SRR8616019 SRR8615549 SRR8615605 SRR8616112 SRR8615727 SRR8615816 SRR8615703 SRR8615378 SRR8615712 SRR8615727 SRR8615423 SRR8615739 SRR8618307 SRR8615971 SRR8615727 SRR8615497 SRR8616096 SRR8615376 SRR8616005 SRR8618301 SRR8615422 SRR8615429 SRR8618305 SRR8615946 SRR8615697 SRR8616012
do
  printf "Processing files associated with accession ${i}."
  eval /home/chughes/softwareTools/megadepth-1.2.0/megadepth ${workingDirectory}/coverageFiles/${i}.sorted.bw --annotation ${workingDirectory}/datasetOutputs/dlg2LongIsoformExons.bed --coverage > ${workingDirectory}/coverageFiles/${i}_Dlg2Megadepth.tsv
  eval /home/chughes/softwareTools/megadepth-1.2.0/megadepth ${workingDirectory}/coverageFiles/${i}.sorted.bw --annotation ${workingDirectory}/datasetOutputs/dlg2n3IsoformExons.bed --coverage > ${workingDirectory}/coverageFiles/${i}_Dlg2n3Megadepth.tsv
done
```