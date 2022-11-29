## Processing EwS nanopore data

This document describes the processing of data from nanopore cDNA sequencing of A673 cells at day 0 or day 7 of dox treatment (shRNA against EWS-FLI1).

### Description

I am going to use [minimap2](https://github.com/lh3/minimap2) for this. I am going to follow their standard protocol for mRNA (https://github.com/lh3/minimap2#map-long-splice).

### Data pipeline

Move do the data directory and start the alignment. We used two flow cells, PAM37258 and PAM39767, and each has two libraries multiplexed on it:

PAM37258
noDox = CAAGAAAGTTGTCGGTGTCTTTGTGAC
yesDoc = CTCGATTCCGTTTGTAGTCGTCTGTAC

PAM39767
noDox = CAAGAAAGTTGTCGGTGTCTTTGTGAC
yesDox = CTCGATTCCGTTTGTAGTCGTCTGTAC

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore
```

I am following the protocol detailed [here](https://github.com/epi2me-labs/wf-transcriptomes). Process the data with [pychopper](https://github.com/epi2me-labs/pychopper) prior to alignment. NOTE - I ended up not using pychopper because StringTie2 seems to perform equally as well without using it, and it is slow. So I skipped right to the next step with minimap2.

```shell
##
pychopper -r /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_noDox_report.pdf -u /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_noDox_unclassified.fq -w /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_noDox_rescued.fq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM37258/PAM37258_CAAGAAAGTTGTCGGTGTCTTTGTGAC_pass_fastq_concat.fastq.gz /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_noDox_full_length_output.fq

##
pychopper -r /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_yesDox_report.pdf -u /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_yesDox_unclassified.fq -w /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_yesDox_rescued.fq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM37258/PAM37258_CTCGATTCCGTTTGTAGTCGTCTGTAC_pass_fastq_concat.fastq.gz /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_yesDox_full_length_output.fq

##
pychopper -r /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_noDox_report.pdf -u /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_noDox_unclassified.fq -w /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_noDox_rescued.fq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM39767/PAM39767_CAAGAAAGTTGTCGGTGTCTTTGTGAC_pass_fastq_concat.fastq.gz /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_noDox_full_length_output.fq

##
pychopper -r /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_yesDox_report.pdf -u /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_yesDox_unclassified.fq -w /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_yesDox_rescued.fq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM39767/PAM39767_CTCGATTCCGTTTGTAGTCGTCTGTAC_pass_fastq_concat.fastq.gz /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_yesDox_full_length_output.fq
```

Perform the alignment of the pychopper processed reads using [minimap2](https://github.com/lh3/minimap2). 

```shell
##alignment calls
/home/chughes/softwareTools/minimap2-2.24/minimap2-2.24_x64-linux/minimap2 -ax splice /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.fa /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM37258/PAM37258_CAAGAAAGTTGTCGGTGTCTTTGTGAC_pass_fastq_concat.fastq.gz /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM39767/PAM39767_CAAGAAAGTTGTCGGTGTCTTTGTGAC_pass_fastq_concat.fastq.gz > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.aligned.sam

##
/home/chughes/softwareTools/minimap2-2.24/minimap2-2.24_x64-linux/minimap2 -ax splice /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.fa /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM37258/PAM37258_CTCGATTCCGTTTGTAGTCGTCTGTAC_pass_fastq_concat.fastq.gz /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM39767/PAM39767_CTCGATTCCGTTTGTAGTCGTCTGTAC_pass_fastq_concat.fastq.gz > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/yesDox.aligned.sam
```

Sort the output results in preparation for additional processing.

```shell
##
samtools sort /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.aligned.sam -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.aligned.sorted.bam

##
samtools sort /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/yesDox.aligned.sam -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/yesDox.aligned.sorted.bam
```

Process the files using [StringTie2](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#mix). I am going to use the mix mode that combines short and long read data to improve assignments.

```shell

```








First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211210_soleMscToEwsPmid34341072/rnaSeq
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
BASE_DIR = "/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211210_soleMscToEwsPmid34341072/rnaSeq"


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
#workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211210_soleMscToEwsPmid34341072/rnaSeq"
workingDirectory="/projects/ptx_results/Sequencing/publishedStudies/sequencing20211210_soleMscToEwsPmid34341072/rnaSeq"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results
eval mkdir quants

##loop over the accessions
for i in SRR1{{1807630..1807666},{4743608..4743630}}
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
