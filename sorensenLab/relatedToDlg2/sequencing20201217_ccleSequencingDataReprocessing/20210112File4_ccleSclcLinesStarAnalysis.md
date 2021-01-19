## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a the CCLE database, or DepMap. The data are all located in SRA: [PRJNA523380](https://www.ncbi.nlm.nih.gov/bioproject?term=PRJNA523380).

### Getting the raw data

I am interested in all of the Ewing sarcoma data now. I had originally done the A673 and SKNMC files, but I will expand it to include the other lines now. the files are:

```
SRX5414898: RNAseq of NCIH69_LUNG
SRR8616152

SRX5414889: RNAseq of NCIH526_LUNG
SRR8616161

SRX5414938: RNAseq of NCIH2227_LUNG
SRR8616112

SRX5414330: RNAseq of CORL279_LUNG
SRR8615423

SRX5414412: RNAseq of NCIH1048_LUNG
SRR8615341

SRX5415047: RNAseq of DMS153_LUNG
SRR8616003
```

I am going to use SRATools to get at these data. I am going to save them in the directory `/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/additionalEwingLines`. There is a nice little tutorial on using SRATools to get raw data [here](https://www.biostars.org/p/111040/). On the proteomics-svr02 system, SRATools is stored in `/projects/ptx_analysis/chughes/software/sratoolkit.2.9.6-1-centos_linux64/bin`. I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). Because I want to look at junction reads, I am going to do an alignment with STAR. There is a good manual [here](https://github.com/alexdobin/STAR) that I am following for the most part. The first thing I want to do is get an index for STAR. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the STAR index with the command `refgenie pull hg38/star_index`. I put all of these functions into a script called `allDataProcessing.sh` and run writing a log of the output to a file called `dataProcessingLog.txt`.

```shell
#!/bin/bash
sraToolsLocation="/projects/ptx_analysis/chughes/software/sratoolkit.2.9.6-1-centos_linux64/bin/prefetch"
fastqDumpLocation="/projects/ptx_analysis/chughes/software/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump"
bbdukLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbduk.sh"
starLocation="/projects/ptx_analysis/chughes/software/STAR-2.7.1a/STAR-2.7.1a/bin/Linux_x86_64/STAR"
referenceLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/star_index/default/"
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sclcLines/"
for i in SRR8616152 SRR8616161 SRR8616112 SRR8615423 SRR8615341 SRR8616003
do
  echo $i
  eval "cd ${rawDataOutputDirectory}"
  eval "mkdir ${rawDataOutputDirectory}starResults"
  eval $sraToolsLocation -O ${rawDataOutputDirectory} -v $i
  fastqCall="$fastqDumpLocation --gzip --outdir $rawDataOutputDirectory --split-files ${rawDataOutputDirectory}${i}.sra"
  eval $fastqCall
  bbdukCall="$bbdukLocation in1=${rawDataOutputDirectory}${i}_1.fastq.gz in2=${rawDataOutputDirectory}${i}_2.fastq.gz ref=adapters out1=${rawDataOutputDirectory}${i}_1.clean.fastq.gz out2=${rawDataOutputDirectory}${i}_2.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo"
  eval $bbdukCall
  starCall="$starLocation --runThreadN 12 --genomeDir ${referenceLocation} --readFilesIn ${rawDataOutputDirectory}${i}_1.clean.fastq.gz ${rawDataOutputDirectory}${i}_2.clean.fastq.gz --readFilesCommand zcat --sjdbGTFfile ${annotationLocation} --outFileNamePrefix ${rawDataOutputDirectory}starResults/${i}_ --outSAMtype BAM Unsorted"
  eval $starCall
  bamSortCall="$samtoolsLocation sort ${rawDataOutputDirectory}/starResults/${i}_Aligned.out.bam -o ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamSortCall
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamIndexCall
done
```

The run command in a screen session was:

```shell
./allDataProcessing.sh |& tee -a dataProcessingLog.txt
```

After the data processing was done I went back and removed the .sra files because they are quite large and not needed anymore. Now we can transfer these files to our local machine and look through the bam alignments in IGV.

I am adding to this after the fact. I want to create bigwig files that I can use for coverage analysis downstream. I can do this directly from the chr11 bam files that I created in the above code.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sclcLines/starResults/"
#########################################
for i in SRR8616152 SRR8616161 SRR8616112 SRR8615423 SRR8615341 SRR8616003
do
  echo $i
  coverageCall="bamCoverage -b ${rawDataOutputDirectory}${i}.chr11.bam -o ${rawDataOutputDirectory}${i}.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --centerReads -p 6"
  eval $coverageCall
done
```

I use these files in the `20210113File4_creatingDlg2CoverageMaps.Rmd` analysis file.
