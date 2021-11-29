## Reprocessing some WGS data

This document describes the reprocessing of some WGS data from a previous publication. Specifically:

"The Ewing Sarcoma Cell Line Atlas (ESCLA)"
GEO: GSE176339, SRA: PRJNA610192, https://www.biorxiv.org/content/10.1101/2021.06.08.447518v1

### Description

I am interested in the WGS data across the cell lines they profiled. I am going to follow what the original authors did pretty closely, the biggest difference being that I will use an updated annotation. From their manuscript:

DNA was extracted from wildtype EwS cell lines using the NucleoSpin Tissue kit (Macherey Nagel) following manufacturer’s protocol and eluted in H2O. For sequencing, 50 µl of 50 ng/µl DNA were used. After initial DNA quality assessment on a bioanalyzer (DNA Integrity Number at least 7.0), DNA was sequenced on Illumina HiSeq Xten (150 bp, paired-end; Illumina, San Diego, USA) and a PCR-free protocol at the Genomics and Proteomics Core Facility of the German Cancer Research Center (DKFZ, Heidelberg, Germany). Raw sequencing data were aligned to the human reference genome (version hg19) with Burrows-Wheeler Aligner (bwa mem)70. PhiX contamination was excluded, Illumina adapters and duplicates were marked with picard (Broad Institute). Base quality scores were recalibrated with GATK (Broad Institute). Quality of the final alignments was controlled with FASTQC71. Coverage was assessed with samtools depth72 and displayed as average of 90 kb bins. De novo motif finding was performed with HOMER66.

### Data pipeline

First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia
touch sraDataProcessingScript.sh
chmod +x sraDataProcessingScript.sh
touch snakefile
```

Edit the contents of the snakefile to include the text below. I use vim for this.

```python
"""
Author: Christopher Hughes
Affiliation: BCCRC
Aim: Workflow for ChIP-Seq data
Date: 20211129
"""

###############################
#working directory
BASE_DIR = "/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211129_grunewaldEwsAtlasWgs"


###############################
#locations of tools we will use
BBDUK = "/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
BWA = "/home/chughes/softwareTools/bwa-mem2-2.2.1/bwa-mem2"
SAMTOOLS="/home/chughes/softwareTools/samtools-1.12/samtools"
SAMBAMBA="/home/chughes/softwareTools/sambamba-0.8.1/sambamba"
BAMCOVERAGE="~/virtualPython368/bin/bamCoverage"

###############################
#locations of our index files
DATABASE_DIR = "/home/chughes/databases/projectEwsDlg2"
INDEX = DATABASE_DIR + "/baseGenomeFiles/genome.fa"
GTF   = DATABASE_DIR + "/baseGenomeFiles/genome.gtf"
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
      expand("results/{smp}.bam.bai", smp = SAMPLES),
      expand("results/{smp}.bw", smp = SAMPLES)

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

rule bwa:
  input:
      r1 = "results/{smp}_1.clean.fastq.gz",
      r2 = "results/{smp}_2.clean.fastq.gz"
  output:
      "results/{smp}.unsorted.sam"
  message:
      "Aligning with bwa-mem2."
  shell:
      "{BWA} mem -t 16 {INDEX} {input.r1} {input.r2} > {output}"

rule bam_conversion:
  input:
      "results/{smp}.unsorted.sam"
  output:
      "results/{smp}.unsorted.bam"
  message:
      "SAM to BAM conversion with samtools."
  shell:
      "{SAMTOOLS} view -h -b -o {output} {input}"

rule bam_sorting:
  input:
      "results/{smp}.unsorted.bam"
  output:
      "results/{smp}.sorted.bam"
  message:
      "BAM sorting with sambamba."
  shell:
      "{SAMBAMBA} sort -t 6 -o {output} {input}"

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
      "results/{smp}.sorted.bam"
  output:
      "results/{smp}.bw"
  message:
      "Calculating coverage with deeptools."
  shell:
      "{BAMCOVERAGE} -b {input} -o {output} -p 8"
```

The code below downloads the raw data with [sraDownloader](https://github.com/s-andrews/sradownloader) and sends it to snakemake.

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
sraDownloader="/home/chughes/softwareTools/sradownloader-3.8/sradownloader"
sraCacheLocation="/mnt/Data/chughes/sratoolsRepository"
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211129_grunewaldEwsAtlasWgs"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results

##loop over the accessions
for i in SRR11235{318..335}
do
  printf "Downloading files associated with ${i}."
  eval ${sraDownloader} --outdir ${workingDirectory}/raw ${i}
  ##the file gets renamed upon download, but I just want it to have the SRR id and I can annotate it later
  eval mv ${workingDirectory}/raw/${i}*_1.fastq.gz ${workingDirectory}/raw/${i}_1.fastq.gz
  eval mv ${workingDirectory}/raw/${i}*_2.fastq.gz ${workingDirectory}/raw/${i}_2.fastq.gz
  #eval conda activate snakemake
  eval snakemake --cores 8
  #eval conda deactivate
  eval rm ${workingDirectory}/raw/${i}*.fastq.gz
  eval rm ${sraCacheLocation}/sra/${i}*.sra.cache
  eval rm ${workingDirectory}/results/${i}*.clean.fastq.gz
  eval rm ${workingDirectory}/results/${i}.unsorted.sam
  eval rm ${workingDirectory}/results/${i}.unsorted.bam
done
```