## Reprocessing some ChIPseq data

his document describes the re-processing of some ChIPseq data from neuroblastoma cells published in the paper:

Zeid R, Lawlor MA, Poon E, Reyes JM et al. Enhancer invasion shapes MYCN-dependent transcriptional amplification in neuroblastoma. Nat Genet 2018 Apr;50(4):515-523. PMID: 29379199, GEO: GSE80151

### Description

I am interested in the ChIPseq data across the cell lines they profiled. There is a metadata file in the server directory detailed below that describes these data.

I am going to use [SRAExplorer](https://sra-explorer.info/#) to get at these data. I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). Because the authors originally used Bowtie for the alignment, I will as well. There is a pretty good walkthrough of the alignment process [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/03_align_and_filtering.html). I will then use [samtools](http://www.htslib.org/) and [sambamba](https://lomereiter.github.io/sambamba/) to sort and filter the sam and bam files that are generated, prior to peak calling. I ended up writing this into a script that uses snakemake. Proceed with the analysis as below.

### Data pipeline

First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20220711_nbCellLineMycNPmid29379199/chipSeq
touch sraDataProcessingScript.sh
chmod +x sraDataProcessingScript.sh
touch snakefile
```

Edit the contents of the snakefile to include the text below. I use vim for this.

```shell
"""
Author: Christopher Hughes
Affiliation: BCCRC
Aim: Workflow for ChIP-Seq data
Date: 20211030
"""

###############################
#locations of tools we will use
#BBDUK = "/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
BBDUK = "/projects/ptx_analysis/chughes/softwareTools/bbmap-38.90/bbduk.sh"
#BOWTIE2 = "/home/chughes/softwareTools/bowtie2-2.4.4/bowtie2"
BOWTIE2 = "/projects/ptx_analysis/chughes/softwareTools/bowtie2-2.4.4/bowtie2"
#SALMON = "/home/chughes/softwareTools/salmon-1.5.2/bin/salmon"
SALMON = "/projects/ptx_analysis/chughes/softwareTools/salmon-1.5.2/bin/salmon"
#HISAT2 = "/home/chughes/softwareTools/hisat2-2.2.1/hisat2"
#STAR = "/home/chughes/softwareTools/STAR-2.7.9a/bin/Linux_x86_64/STAR"
STAR = "/projects/ptx_analysis/chughes/softwareTools/STAR-2.7.9a/bin/Linux_x86_64/STAR"
#SAMTOOLS="/home/chughes/softwareTools/samtools-1.12/samtools"
SAMTOOLS="/gsc/software/linux-x86_64-centos7/samtools-1.14/bin/samtools"
#SAMBAMBA="/home/chughes/softwareTools/sambamba-0.8.1/sambamba"
SAMBAMBA="/projects/ptx_analysis/chughes/softwareTools/sambamba-0.8.1/sambamba"
#FEATURECOUNTS="/home/chughes/softwareTools/subread-2.0.3/bin/featureCounts"
FEATURECOUNTS="/projects/ptx_analysis/chughes/softwareTools/subread-2.0.3/bin/featureCounts"
#BAMCOVERAGE="/home/chughes/virtualPython368/bin/bamCoverage"
BAMCOVERAGE="/home/chughes/Virtual_Python383/bin/bamCoverage"


###############################
#locations of our index files
#DATABASE_DIR = "/home/chughes/databases/projectEwsDlg2"
DATABASE_DIR = "/projects/ptx_analysis/chughes/databases/projectEwsDlg2"
INDEX = DATABASE_DIR + "/bowtie2Index/human_genome"
GTF = DATABASE_DIR + "/baseGenomeFiles/genome.gtf"
FASTA = DATABASE_DIR + "/baseGenomeFiles/genome.fa"


###############################
#get our list of sample files
SAMPLES, = glob_wildcards("raw/{smp}.fastq.gz")
print("There are " + str(len(SAMPLES)) + " total samples to be processed.")
for smp in SAMPLES:
  print("Sample " + smp + " will be processed")


###############################
#processing workflow
rule all:
    input: 
      expand("results/{smp}.filtered.bam.bai", smp = SAMPLES),
      expand("results/{smp}.filtered.bw", smp = SAMPLES)

rule bbduk:
  input:
      "raw/{smp}.fastq.gz"
  output:
      "results/{smp}.clean.fastq.gz"
  message:
      "Processing with BBDuk."
  shell:
      "{BBDUK} in={input} ref=adapters out={output} ktrim=r k=23 mink=11 hdist=1"

rule bowtie:
  input:
      "results/{smp}.clean.fastq.gz"
  output:
      "results/{smp}.unsorted.sam"
  message:
      "Aligning with bowtie2."
  shell:
      "{BOWTIE2} -p 8 -q --local -x {INDEX} -U {input} -S {output}"

rule bam_conversion:
  input:
      "results/{smp}.unsorted.sam"
  output:
      "results/{smp}.unsorted.bam"
  message:
      "SAM to BAM conversion with samtools."
  shell:
      "{SAMTOOLS} view -h -S -b -o {output} {input}"

rule bam_sorting:
  input:
      "results/{smp}.unsorted.bam"
  output:
      "results/{smp}.sorted.bam"
  message:
      "BAM sorting with sambamba."
  shell:
      "{SAMBAMBA} sort -t 6 -o {output} {input}"

rule bam_filtering:
  input:
      "results/{smp}.sorted.bam"
  output:
      "results/{smp}.filtered.bam"
  message:
      "Filtering BAM with sambamba."
  shell:
      """{SAMBAMBA} view -h -t 6 -f bam -F "[XS] == null and not unmapped and not duplicate" {input} > {output}"""

rule bam_indexing:
  input:
      "results/{smp}.filtered.bam"
  output:
      "results/{smp}.filtered.bam.bai"
  message:
      "Indexing BAM with samtools."
  shell:
      "{SAMTOOLS} index {input}"

rule bam_coverage:
  input:
      r1 = "results/{smp}.filtered.bam",
      w2 = "results/{smp}.filtered.bam.bai"
  output:
      "results/{smp}.filtered.bw"
  message:
      "Calculating coverage with deeptools."
  shell:
      "{BAMCOVERAGE} -b {input.r1} -o {output} --binSize 10 --normalizeUsing BPM --smoothLength 30 --extendReads 150 --centerReads -p 6"
```

Below is the shell script I will use to process these data with snakemake.

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
sraDownloader="/projects/ptx_analysis/chughes/softwareTools/sradownloader-3.8/sradownloader"
#sraCacheLocation="/mnt/Data/chughes/sratoolsRepository"
sraCacheLocation="/projects/ptx_results/Sequencing/sraCache"
#workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20220711_nbCellLineMycNPmid29379199/chipSeq"
workingDirectory="/projects/ptx_results/Sequencing/publishedStudies/sequencing20220711_nbCellLineMycNPmid29379199/chipSeq"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results
eval mkdir quants

##loop over the accessions
for i in SRR3356{395..397} SRR3356{400..402} SRR3356{404..408} SRR5441{926..934}
do
  printf "Downloading files associated with ${i}."
  eval ${sraDownloader} --outdir ${workingDirectory}/raw ${i}
  ##the file gets renamed upon download, but I just want it to have the SRR id and I can annotate it later
  eval mv ${workingDirectory}/raw/${i}*.fastq.gz ${workingDirectory}/raw/${i}.fastq.gz
  #eval mv ${workingDirectory}/raw/${i}*_2.fastq.gz ${workingDirectory}/raw/${i}_2.fastq.gz
  #eval conda activate snakemake
  eval snakemake --cores 8 --latency-wait 300
  #eval conda deactivate
  eval rm ${workingDirectory}/raw/${i}*.fastq.gz
  eval rm ${sraCacheLocation}/sra/${i}*
  eval rm ${workingDirectory}/results/${i}*.clean.fastq.gz
  eval rm ${workingDirectory}/results/${i}*.sam
  eval rm ${workingDirectory}/results/${i}*.sorted.bam
done
```

Note - I had to go back and reprocess these data for coverage because I mistakenly had it only calculate for chr11, so this is why there is a 'coverageMapping.sh' script in the results directory. Otherwise, everything was done as above. Now we move on to peak calling using [MACS](https://github.com/macs3-project/MACS). There is a great walkthrough of this process [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html). Also [here](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) for discussion of deepTools, specifically bamCoverage.

The first file here is for narrow peak calling on the transcription factors.

```shell
#!/bin/bash
#annotationLocation="/home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf"
annotationLocation="/projects/ptx_analysis/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf"
#rawDataOutputDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20220711_nbCellLineMycNPmid29379199/chipSeq"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/sequencing20220711_nbCellLineMycNPmid29379199/chipSeq"
for i in SRR5441928
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR5441927.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR5441931
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR5441930.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR5441934
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR5441933.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR3356397
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR3356396.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR3356408
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR3356406.filtered.bam ${rawDataOutputDirectory}/results/SRR3356407.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR3356404
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR3356401.filtered.bam ${rawDataOutputDirectory}/results/SRR3356402.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done
```

This second file is for broad peak calling with the histone marks.

```shell
#!/bin/bash
#annotationLocation="/home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf"
annotationLocation="/projects/ptx_analysis/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf"
#rawDataOutputDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20220711_nbCellLineMycNPmid29379199/chipSeq"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/sequencing20220711_nbCellLineMycNPmid29379199/chipSeq"
for i in SRR5441926
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR5441927.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR5441929
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR5441930.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR5441932
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR5441933.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR3356395
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR3356396.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR3356405
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR3356406.filtered.bam ${rawDataOutputDirectory}/results/SRR3356407.filtered.bam-f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR3356400
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR3356401.filtered.bam ${rawDataOutputDirectory}/results/SRR3356402.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done
```

The rest of the analysis we can perform in R.