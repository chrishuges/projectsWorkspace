## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a previous publication. Specifically:

"EWS-FLI1 utilizes divergent chromatin remodeling mechanisms to directly activate or repress enhancer elements in Ewing sarcoma"
Cancer Cell, 2014, Pubmed ID: 25453903, GEO: GSE61953

### Getting the raw data

I am first interested in just the data for the A673 and SKNMC cell lines where it is +/-EWS-FLI1. As far as I can tell, there is just a single replicate run for these samples with a GFP vector used as a control. The data are detailed below:

```
SRX718181: GSM1517630: RNA-seq A673 shGFP 48 hrs; Homo sapiens; RNA-Seq
SRR1594024

SRX718182: GSM1517631: RNA-seq A673 shFLI1 48 hrs; Homo sapiens; RNA-Seq
SRR1594025

SRX718177: GSM1517626: RNA-seq SKNMC shGFP 48 hrs; Homo sapiens; RNA-Seq
SRR1594020

SRX718178: GSM1517627: RNA-seq SKNMC shFLI1 48 hrs; Homo sapiens; RNA-Seq
SRR1594021

SRX718179: GSM1517628: RNA-seq SKNMC shGFP 96 hrs; Homo sapiens; RNA-Seq
SRR1594022

SRX718180: GSM1517629: RNA-seq SKNMC shFLI1 96 hrs; Homo sapiens; RNA-Seq
SRR1594023
```

I am going to use SRATools to get at these data. I am going to save them in the directory `/projects/ptx_results/Sequencing/publishedStudies/201411RiggiCancerCellPmid25453903`. There is a nice little tutorial on using SRATools to get raw data [here](https://www.biostars.org/p/111040/). On the proteomics-svr02 system, SRATools is stored in `/projects/ptx_analysis/chughes/software/sratoolkit.2.9.6-1-centos_linux64/bin`. To download the files, I wrote a small shell script called `sraPreFetch.sh`.

```shell
#!/bin/bash
sraToolsLocation="/projects/ptx_analysis/chughes/software/sratoolkit.2.9.6-1-centos_linux64/bin/prefetch"
for i in SRR1594024 SRR1594025 SRR1594020 SRR1594021 SRR1594022 SRR1594023
do
  echo $i
  eval $sraToolsLocation -v $i
done
```

I executed it to run the prefetch.

```shell
./sraPreFetch.sh
```

For all of these files, I want to get the raw fastq data for each of the reads. I will write another shell script for this called `fastqPairedDump.sh`.

```shell
#!/bin/bash
fastqDumpLocation="/projects/ptx_analysis/chughes/software/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201411RiggiCancerCellPmid25453903/"
for i in SRR1594024 SRR1594025 SRR1594020 SRR1594021 SRR1594022 SRR1594023
do
  echo $i
  fastqCall="$fastqDumpLocation --gzip --outdir $rawDataOutputDirectory --split-files /home/chughes/ncbi/public/sra/${i}.sra"
  eval $fastqCall
done
```

I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). I wrote a shell script called `filterAdaptersBbduk.sh` to call this script on each file.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201411RiggiCancerCellPmid25453903/"
bbdukLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbduk.sh"
for i in SRR1594024 SRR1594025 SRR1594020 SRR1594021 SRR1594022 SRR1594023
do
  echo $i
  bbdukCall="$bbdukLocation in1=${rawDataOutputDirectory}${i}_1.fastq.gz in2=${rawDataOutputDirectory}${i}_2.fastq.gz ref=adapters out1=${rawDataOutputDirectory}${i}_1.clean.fastq.gz out2=${rawDataOutputDirectory}${i}_2.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo"
  eval $bbdukCall
done
```

It looks like basically all of the reads pass the filtering criteria, just a few bases trimmed on the ends. So, that is good. Here is an example of the output for one of the files:

```
Input is being processed as paired
Started output streams: 0.109 seconds.
Processing time:                115.936 seconds.

Input:                          124274650 reads                 4665411805 bases.
KTrimmed:                       119470 reads (0.10%)    1355476 bases (0.03%)
Trimmed by overlap:             0 reads (0.00%)         0 bases (0.00%)
Total Removed:                  4 reads (0.00%)         1355476 bases (0.03%)
Result:                         124274646 reads (100.00%)       4664056329 bases (99.97%)

Time:                           116.181 seconds.
Reads Processed:        124m    1069.67k reads/sec
Bases Processed:       4665m    40.16m bases/sec
```

Because I want to look at junction reads, I am going to do an alignment with STAR. There is a good manual [here](https://github.com/alexdobin/STAR) that I am following for the most part. The first thing I want to do is get an index for STAR. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the Salmon index with the command `refgenie pull hg38/star_index`. You don't need to do this every time, just when you need to update your index. 

We can now run STAR. We will do this using a script called `starAlignment.sh`.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201411RiggiCancerCellPmid25453903/"
starLocation="/projects/ptx_analysis/chughes/software/STAR-2.7.1a/STAR-2.7.1a/bin/Linux_x86_64/STAR"
referenceLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/star_index/default/"
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"
###########################################
for i in SRR1594020 SRR1594021 SRR1594022 SRR1594023 SRR1594024 SRR1594025
do
  echo $i
  starCall="$starLocation --runThreadN 12 --genomeDir ${referenceLocation} --readFilesIn ${rawDataOutputDirectory}${i}_1.clean.fastq.gz ${rawDataOutputDirectory}${i}_2.clean.fastq.gz --readFilesCommand zcat --sjdbGTFfile ${annotationLocation} --outFileNamePrefix ${rawDataOutputDirectory}starResults/${i}_ --outSAMtype BAM Unsorted"
  eval $starCall
  bamSortCall="$samtoolsLocation sort ${rawDataOutputDirectory}/starResults/${i}_Aligned.out.bam -o ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamSortCall
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamIndexCall
done
```
Now we can transfer these files to our local machine and look through the bam alignments in IGV.

