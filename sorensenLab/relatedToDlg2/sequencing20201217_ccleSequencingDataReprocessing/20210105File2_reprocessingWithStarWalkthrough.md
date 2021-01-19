## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a the CCLE database, or DepMap. The data are all located in SRA: [PRJNA523380](https://www.ncbi.nlm.nih.gov/bioproject?term=PRJNA523380).

### Getting the raw data

I am interested in all of the Ewing sarcoma data, to start. Later I would like to compare with brain tissue if I can find data. The data are detailed below:

```
SRX5415038: RNAseq of A673_BONE
SRR8616012

SRX5414256: RNAseq of SKNMC_BONE
SRR8615497
```

I am going to use SRATools to get at these data. I am going to save them in the directory `/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData`. There is a nice little tutorial on using SRATools to get raw data [here](https://www.biostars.org/p/111040/). On the proteomics-svr02 system, SRATools is stored in `/projects/ptx_analysis/chughes/software/sratoolkit.2.9.6-1-centos_linux64/bin`. To download the files, I wrote a small shell script called `sraPreFetch.sh`.

```shell
#!/bin/bash
sraToolsLocation="/projects/ptx_analysis/chughes/software/sratoolkit.2.9.6-1-centos_linux64/bin/prefetch"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/"
for i in SRR8616012 SRR8615497
do
  echo $i
  eval $sraToolsLocation -O ${rawDataOutputDirectory} -v $i
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
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/"
for i in SRR8616012 SRR8615497
do
  echo $i
  fastqCall="$fastqDumpLocation --gzip --outdir $rawDataOutputDirectory --split-files ${rawDataOutputDirectory}${i}.sra"
  eval $fastqCall
done
```

I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). I wrote a shell script called `filterAdaptersBbduk.sh` to call this script on each file.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/"
bbdukLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbduk.sh"
for i in SRR8616012 SRR8615497
do
  echo $i
  bbdukCall="$bbdukLocation in1=${rawDataOutputDirectory}${i}_1.fastq.gz in2=${rawDataOutputDirectory}${i}_2.fastq.gz ref=adapters out1=${rawDataOutputDirectory}${i}_1.clean.fastq.gz out2=${rawDataOutputDirectory}${i}_2.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo"
  eval $bbdukCall
done
```

The results look pretty normal with almost everything passing the filter.

```shell
SRR8616012
java -ea -Xmx25167m -Xms25167m -cp /projects/ptx_analysis/chughes/software/bbmap_v38_87/current/ jgi.BBDuk in1=/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/SRR8616012_1.fastq.gz in2=/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/SRR8616012_2.fastq.gz ref=adapters out1=/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/SRR8616012_1.clean.fastq.gz out2=/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/SRR8616012_2.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo
Executing jgi.BBDuk [in1=/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/SRR8616012_1.fastq.gz, in2=/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/SRR8616012_2.fastq.gz, ref=adapters, out1=/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/SRR8616012_1.clean.fastq.gz, out2=/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/SRR8616012_2.clean.fastq.gz, ktrim=r, k=23, mink=11, hdist=1, tpe, tbo]
Version 38.87

maskMiddle was disabled because useShortKmers=true
0.151 seconds.
Initial:
Memory: max=25291m, total=25291m, free=24367m, used=924m

Added 217135 kmers; time:       0.416 seconds.
Memory: max=25291m, total=25291m, free=23311m, used=1980m

Input is being processed as paired
Started output streams: 0.354 seconds.
Processing time:                497.505 seconds.

Input:                          220756414 reads                 22296397814 bases.
KTrimmed:                       10016050 reads (4.54%)  760071862 bases (3.41%)
Trimmed by overlap:             4431678 reads (2.01%)   22474972 bases (0.10%)
Total Removed:                  6247296 reads (2.83%)   782546834 bases (3.51%)
Result:                         214509118 reads (97.17%)        21513850980 bases (96.49%)

Time:                           498.278 seconds.
Reads Processed:        220m    443.04k reads/sec
Bases Processed:      22296m    44.75m bases/sec
```

Because I want to look at junction reads, I am going to do an alignment with STAR. There is a good manual [here](https://github.com/alexdobin/STAR) that I am following for the most part. The first thing I want to do is get an index for STAR. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the Salmon index with the command `refgenie pull hg38/star_index`. You don't need to do this every time, just when you need to update your index. 

We can now run STAR. We will do this using a script called `starAlignment.sh`.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/"
starLocation="/projects/ptx_analysis/chughes/software/STAR-2.7.1a/STAR-2.7.1a/bin/Linux_x86_64/STAR"
referenceLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/star_index/default/"
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"
###########################################
for i in SRR8616012 SRR8615497
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

I am adding to this after the fact. I want to create bigwig files that I can use for coverage analysis downstream. I can do this directly from the chr11 bam files that I created in the above code.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/"
#########################################
for i in SRR8616012 SRR8615497
do
  echo $i
  coverageCall="bamCoverage -b ${rawDataOutputDirectory}${i}.chr11.bam -o ${rawDataOutputDirectory}${i}.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --centerReads -p 6"
  eval $coverageCall
done
```

I use these files in the `20210113File4_creatingDlg2CoverageMaps.Rmd` analysis file.

