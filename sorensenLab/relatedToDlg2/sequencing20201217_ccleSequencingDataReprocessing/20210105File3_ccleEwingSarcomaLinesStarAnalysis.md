## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a the CCLE database, or DepMap. The data are all located in SRA: [PRJNA523380](https://www.ncbi.nlm.nih.gov/bioproject?term=PRJNA523380).

### Getting the raw data

I am interested in all of the Ewing sarcoma data now. I had originally done the A673 and SKNMC files, but I will expand it to include the other lines now. the files are:

```
SRX5415189: RNAseq of EWS502_BONE
SRR8616213

SRX5415188: RNAseq of EW8_BONE
SRR8616214

SRX5414808: RNAseq of RDES_BONE
SRR8615592

SRX5414721: RNAseq of SKES1_BONE
SRR8615679

SRX5414568: RNAseq of TC71_BONE
SRR8615832

SRX5414541: RNAseq of SKPNDW_BONE
SRR8615859

SRX5414480: RNAseq of CADOES1_BONE
SRR8615273

SRX5414254: RNAseq of SKNEP1_BONE
SRR8615499

SRX5414232: RNAseq of MHHES1_BONE
SRR8615521
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
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/additionalEwingLines/"
for i in SRR8616213 SRR8616214 SRR8615592 SRR8615679 SRR8615832 SRR8615859 SRR8615273 SRR8615499 SRR8615521
do
  echo $i
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

