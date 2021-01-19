## Reprocessing some ChIPseq data

This document describes the reprocessing of some ChIPseq data from a previous publication. Specifically:

"Transcriptional Programs Define Intratumoral Heterogeneity of Ewing Sarcoma at Single-Cell Resolution"
Cell Reports, 2020, Pubmed ID: 32049009, GEO: GSE130025

### Getting the raw data

I am interested in the ChIPseq data across the timecourse that they ran, as well as the H3K27ac data as well. They do provide bed files of the peak calls, but it is easier to plot if you have the bam files. The data are detailed below:

```
SRX5620946: GSM3701309: ASP14_d0_input; Homo sapiens; ChIP-Seq
SRR8832666

SRX5620947: GSM3701310: ASP14_d0_H3K27ac; Homo sapiens; ChIP-Seq
SRR8832667

SRX5620948: GSM3701311: ASP14_d7_H3K27ac; Homo sapiens; ChIP-Seq
SRR8832668

SRX5620949: GSM3701312: ASP14_d7_FLI1; Homo sapiens; ChIP-Seq
SRR8832669

SRX5620950: GSM3701313: ASP14_d9_FLI1; Homo sapiens; ChIP-Seq
SRR8832670

SRX5620951: GSM3701314: ASP14_d10_FLI1; Homo sapiens; ChIP-Seq
SRR8832671

SRX5620952: GSM3701315: ASP14_d11_FLI1; Homo sapiens; ChIP-Seq
SRR8832672

SRX5620953: GSM3701316: ASP14_d14_FLI1; Homo sapiens; ChIP-Seq
SRR8832673

SRX5620954: GSM3701317: ASP14_d17_FLI1; Homo sapiens; ChIP-Seq
SRR8832674
```

I am going to use [SRAExplorer](https://sra-explorer.info/#) to get at these data. I am going to save them in the directory `/projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq`. I used the provided shell script code and put it into a file named `fileDownloads.sh`. 

```shell
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR883/009/SRR8832669/SRR8832669.fastq.gz -o SRR8832669.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR883/001/SRR8832671/SRR8832671.fastq.gz -o SRR8832671.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR883/008/SRR8832668/SRR8832668.fastq.gz -o SRR8832668.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR883/007/SRR8832667/SRR8832667.fastq.gz -o SRR8832667.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR883/000/SRR8832670/SRR8832670.fastq.gz -o SRR8832670.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR883/006/SRR8832666/SRR8832666.fastq.gz -o SRR8832666.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR883/003/SRR8832673/SRR8832673.fastq.gz -o SRR8832673.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR883/004/SRR8832674/SRR8832674.fastq.gz -o SRR8832674.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR883/002/SRR8832672/SRR8832672.fastq.gz -o SRR8832672.fastq.gz
```

This tool is so much better than SRAToolkit it is unreal. I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).  

I put the pre-processing and alignment in a script called `preProcessAlignment.sh`.

```shell
#!/bin/bash
bbdukLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbduk.sh"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/"
for i in SRR8832666 SRR8832667 SRR8832668 SRR8832669 SRR8832670 SRR8832671 SRR8832672 SRR8832673 SRR8832674
do
  echo $i
  eval "cd ${rawDataOutputDirectory}"
  bbdukCall="$bbdukLocation in=${rawDataOutputDirectory}${i}.fastq.gz ref=adapters out=${rawDataOutputDirectory}${i}.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1"
  eval $bbdukCall
done
```
Because the authors originally used Bowtie for the alignment, I will as well. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the index with the command `refgenie pull hg38/bowtie2_index`. There is a pretty good walkthrough of the alignment process [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/03_align_and_filtering.html). 

```shell
#!/bin/bash
bowtieLocation="/gsc/software/linux-x86_64-centos7/bowtie2-2.3.4.1/bin/bowtie2"
referenceLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/bowtie2_index/default/hg38"
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"
sambambaLocation="/projects/ptx_analysis/chughes/software/sambamba-0.8.0/sambamba-0.8.0"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/"
for i in SRR8832666 SRR8832667 SRR8832668 SRR8832669 SRR8832670 SRR8832671 SRR8832672 SRR8832673 SRR8832674
do
  echo $i
  eval "cd ${rawDataOutputDirectory}"
  bowtieCall="$bowtieLocation -p 8 -q --local -x $referenceLocation -U ${rawDataOutputDirectory}${i}.clean.fastq.gz -S ${rawDataOutputDirectory}${i}.unsorted.sam"
  eval $bowtieCall
  bamConvertCall="$samtoolsLocation view -h -S -b -o ${rawDataOutputDirectory}${i}.unsorted.bam ${rawDataOutputDirectory}${i}.unsorted.sam"
  eval $bamConvertCall
  bamSortCall="$sambambaLocation sort -t 6 -o ${rawDataOutputDirectory}${i}.sorted.bam ${rawDataOutputDirectory}${i}.unsorted.bam"
  eval $bamSortCall
  eval $sambambaLocation view -h -t 6 -f bam -F "[XS] == null and not unmapped and not duplicate" ${rawDataOutputDirectory}${i}.sorted.bam > ${rawDataOutputDirectory}${i}.filtered.bam
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}${i}.filtered.bam"
  eval $bamIndexCall
done
```

This script didn't completely work. The removing of duplicates command wasn't working, I think because of a text error. It just didn't like the quotes in the command. So, I had to go through the files manually just because I didn't really want to get into a big thing about fixing it. 

```shell
/projects/ptx_analysis/chughes/software/sambamba-0.8.0/sambamba-0.8.0 view -h -t 6 -f bam -F '[XS] == null and not unmapped and not duplicate' ./SRR8832667.sorted.bam > ./SRR8832667.filtered.bam
```

Now we move on to peak calling using [MACS](https://github.com/macs3-project/MACS). There is a great walkthrough of this process [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html). Also [here](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) for discussion of deepTools, specifically bamCoverage.

```shell
#!/bin/bash
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/"
for i in SRR8832669 SRR8832670 SRR8832671 SRR8832672 SRR8832673 SRR8832674
do
  echo $i
  eval macs3 callpeak -t ${rawDataOutputDirectory}${i}.filtered.bam -c ${rawDataOutputDirectory}SRR8832666.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 
done
```

For the H3K27ac data, I rand MACS3 using the --broad tag as well as this is how it was done in the original manuscript. Now we are ready for visualization. There is a good walkthrough of visualization [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html). First, I ran [deeptools](https://deeptools.readthedocs.io/en/develop/) on the bam files to create coverage maps:

```shell
 bamCoverage -b ./SRR8832666.filtered.bam -o ./deeptools/SRR8832666.chr11.bw --binSize 20 --region chr11 --normalizeUsing BPM --smoothLength 60 --extendReads 150 --centerReads -p 6

bamCoverage -b ./$i.filtered.bam -o ./deeptools/$i.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --extendReads 150 --centerReads -p 6
```

Now, I want to create a heatmap around the TSS for genes on chromosome 11. First I need to create a bed file with the chr11 gene regions in it. I can do this using the command line:

```shell
cat /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf |  awk 'OFS="\t" {if ($3=="gene" && $1=="chr11") {print $1,$4-1,$5,$10,$16,$7}}' | tr -d '";' > ./chr11Hg38GeneRegions.bed
```

Now I can run computeMatrix from deeptools.

```shell
computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/chr11Hg38GeneRegions.bed -S /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/SRR88326[67][901234]*.bw --skipZeros -o /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/ewsFli1ChipSeqTssMatrixChr11.gz -p 6 --outFileSortedRegions /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/regionsTssChr11.bed
```

We can create a profile plot from these data.

```shell
plotProfile -m /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/ewsFli1ChipSeqTssMatrixChr11.gz --outFileName /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/regionsTssChr11Profile.pdf --perGroup --refPointLabel "TSS"
```

Alternatively, we can show this as a heatmap.

```shell
plotHeatmap -m /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/ewsFli1ChipSeqTssMatrixChr11.gz --outFileName /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/regionsTssChr11Heatmap.pdf --colorMap RdBu
```

I think for the rest of the visualization, we will carry this out in R.

