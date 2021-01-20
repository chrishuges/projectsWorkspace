## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a the CCLE database, or DepMap. The data are all located in SRA: [PRJNA523380](https://www.ncbi.nlm.nih.gov/bioproject?term=PRJNA523380).

### Getting the raw data

I am interested in all of the Ewing sarcoma data now. I had originally done the A673 and SKNMC files, but I will expand it to include the other lines now. the files are:

```
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/001/SRR8615971/SRR8615971_1.fastq.gz -o SRR8615971_RNAseq_of_QGP1_PANCREAS_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/001/SRR8615971/SRR8615971_2.fastq.gz -o SRR8615971_RNAseq_of_QGP1_PANCREAS_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/002/SRR8615812/SRR8615812_1.fastq.gz -o SRR8615812_RNAseq_of_T47D_BREAST_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/002/SRR8615812/SRR8615812_2.fastq.gz -o SRR8615812_RNAseq_of_T47D_BREAST_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/005/SRR8615675/SRR8615675_1.fastq.gz -o SRR8615675_RNAseq_of_SIMA_AUTONOMIC_GANGLIA_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/005/SRR8615675/SRR8615675_2.fastq.gz -o SRR8615675_RNAseq_of_SIMA_AUTONOMIC_GANGLIA_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/007/SRR8615747/SRR8615747_1.fastq.gz -o SRR8615747_RNAseq_of_RL_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/007/SRR8615747/SRR8615747_2.fastq.gz -o SRR8615747_RNAseq_of_RL_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/005/SRR8615955/SRR8615955_1.fastq.gz -o SRR8615955_RNAseq_of_COV644_OVARY_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/005/SRR8615955/SRR8615955_2.fastq.gz -o SRR8615955_RNAseq_of_COV644_OVARY_2.fastq.gz

QGP1_PANCREAS
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/001/SRR8615971/SRR8615971_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/001/SRR8615971/SRR8615971_2.fastq.gz

T47D_BREAST
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/002/SRR8615812/SRR8615812_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/002/SRR8615812/SRR8615812_2.fastq.gz

SIMA_AUTONOMIC_GANGLIA
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/005/SRR8615675/SRR8615675_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/005/SRR8615675/SRR8615675_2.fastq.gz

RL_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/007/SRR8615747/SRR8615747_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/007/SRR8615747/SRR8615747_2.fastq.gz

COV644_OVARY
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/005/SRR8615955/SRR8615955_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/005/SRR8615955/SRR8615955_2.fastq.gz
```

I am going to use [SRAExplorer](https://sra-explorer.info/#) to get at these data. I am going to save them in the directory `/projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq`. I used the provided shell script code and put it into a file named `fileDownloads.sh`. I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). Because I want to look at junction reads, I am going to do an alignment with STAR. There is a good manual [here](https://github.com/alexdobin/STAR) that I am following for the most part. The first thing I want to do is get an index for STAR. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the STAR index with the command `refgenie pull hg38/star_index`. Lastly, I use [deeptools](https://deeptools.readthedocs.io/en/develop/) to get coverage maps for the samples on chromosome 11. I put all of these functions into a script called `allDataProcessing.sh` and run writing a log of the output to a file called `dataProcessingLog.txt`.

```shell
#!/bin/bash
bbdukLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbduk.sh"
starLocation="/projects/ptx_analysis/chughes/software/STAR-2.7.1a/STAR-2.7.1a/bin/Linux_x86_64/STAR"
referenceLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/star_index/default/"
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/otherLineages/"

#download the raw data files
eval "cd ${rawDataOutputDirectory}"
eval "mkdir ${rawDataOutputDirectory}starResults"

eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/001/SRR8615971/SRR8615971_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/001/SRR8615971/SRR8615971_2.fastq.gz

eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/002/SRR8615812/SRR8615812_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/002/SRR8615812/SRR8615812_2.fastq.gz

eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/005/SRR8615675/SRR8615675_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/005/SRR8615675/SRR8615675_2.fastq.gz

eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/007/SRR8615747/SRR8615747_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/007/SRR8615747/SRR8615747_2.fastq.gz

eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/005/SRR8615955/SRR8615955_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/005/SRR8615955/SRR8615955_2.fastq.gz

#process the downloaded files
for i in SRR8615971 SRR8615812 SRR8615675 SRR8615747 SRR8615955
do
  echo $i
  bbdukCall="$bbdukLocation in1=${rawDataOutputDirectory}${i}_1.fastq.gz in2=${rawDataOutputDirectory}${i}_2.fastq.gz ref=adapters out1=${rawDataOutputDirectory}${i}_1.clean.fastq.gz out2=${rawDataOutputDirectory}${i}_2.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo"
  eval $bbdukCall
  starCall="$starLocation --runThreadN 12 --genomeDir ${referenceLocation} --readFilesIn ${rawDataOutputDirectory}${i}_1.clean.fastq.gz ${rawDataOutputDirectory}${i}_2.clean.fastq.gz --readFilesCommand zcat --sjdbGTFfile ${annotationLocation} --outFileNamePrefix ${rawDataOutputDirectory}starResults/${i}_ --outSAMtype BAM Unsorted"
  eval $starCall
  bamSortCall="$samtoolsLocation sort ${rawDataOutputDirectory}/starResults/${i}_Aligned.out.bam -o ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamSortCall
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamIndexCall
  coverageCall="bamCoverage -b ${rawDataOutputDirectory}/starResults/${i}.sorted.bam -o ${rawDataOutputDirectory}/starResults/${i}.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --centerReads -p 6"
  eval $coverageCall
done
```

The run command in a screen session was:

```shell
./allDataProcessing.sh |& tee -a dataProcessingLog.txt
```

Now we can transfer these files to our local machine and look through the bam alignments in IGV. I use these files in the `20210113File5_creatingDlg2CoverageMaps.Rmd` analysis file.
