## Processing some RNAseq data

This document details analysis of the CHLA cell line RNAseq data with Salmon. There is a great walkthrough [here](http://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html) that I am following for the most part. The first thing I want to do is get an index for Salmon. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the Salmon index with the command `refgenie pull hg38/salmon_partial_sa_index`. You don't need to do this every time, just when you need to update your index. 

The raw data are stored here: `/projects/analysis/analysis32/PX1273/CDP5RANXX_5`. Transfer the data to our own server.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/2019/sequencing20190801_chla9Chla10Hek293RnaSeq/"
rawDataInputDirectory="/projects/analysis/analysis32/PX1273/CDP5RANXX_5/"

######################################
for i in ATCACG CGATGT TTAGGC TGACCA ACAGTG GCCAAT CAGATC ACTTGA GATCAG
do
  echo $i
  eval "scp ${rawDataInputDirectory}PX1273_${i}/75bp/CDP5RANXX_5_1_${i}_75bp.concat_chastity_passed.fastq.gz ${rawDataOutputDirectory}${i}_1.fastq.gz"
  eval "scp ${rawDataInputDirectory}PX1273_${i}/75bp/CDP5RANXX_5_2_${i}_75bp.concat_chastity_passed.fastq.gz ${rawDataOutputDirectory}${i}_2.fastq.gz"
done
```

We can now run Salmon as well as STAR. We will do this using a script called `runSalmonAndStar.sh`.

```shell
#!/bin/bash
bbdukLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbduk.sh"
salmonLocation="/projects/ptx_analysis/chughes/software/salmon-1.4/bin/salmon"
indexLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/salmon_partial_sa_index/default/"
starLocation="/projects/ptx_analysis/chughes/software/STAR-2.7.1a/STAR-2.7.1a/bin/Linux_x86_64/STAR"
referenceLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/star_index/default/"
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/2019/sequencing20190801_chla9Chla10Hek293RnaSeq/"

###########################################
for i in ATCACG CGATGT TTAGGC TGACCA ACAGTG GCCAAT CAGATC ACTTGA GATCAG
do
  echo $i
  bbdukCall="$bbdukLocation in1=${rawDataOutputDirectory}${i}_1.fastq.gz in2=${rawDataOutputDirectory}${i}_2.fastq.gz ref=adapters out1=${rawDataOutputDirectory}${i}_1.clean.fastq.gz out2=${rawDataOutputDirectory}${i}_2.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo"
  eval $bbdukCall
  ##
  salmonCall="$salmonLocation quant -i $indexLocation -l A -1 ${rawDataOutputDirectory}${i}_1.clean.fastq.gz -2 ${rawDataOutputDirectory}${i}_2.clean.fastq.gz --gcBias --validateMappings -o ${rawDataOutputDirectory}${i}_quant"
  eval $salmonCall
  ##
  starCall="$starLocation --runThreadN 12 --genomeDir ${referenceLocation} --readFilesIn ${rawDataOutputDirectory}${i}_1.clean.fastq.gz ${rawDataOutputDirectory}${i}_2.clean.fastq.gz --readFilesCommand zcat --sjdbGTFfile ${annotationLocation} --outFileNamePrefix ${rawDataOutputDirectory}starResults/${i}_ --outSAMtype BAM Unsorted"
  eval $starCall
  ##
  bamSortCall="$samtoolsLocation sort ${rawDataOutputDirectory}/starResults/${i}_Aligned.out.bam -o ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamSortCall
  ##
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamIndexCall
done
```

For downstream analysis, file file layout is:

sampleName | cell | treatment | barcode | experiment
| --- | --- | --- | --- | --- |
chla9_rep1 | chla9 | none | ATCACG | PX1273
chla9_rep2 | chla9 | none | CGATGT | PX1273
chla9_rep3 | chla9 | none | TTAGGC | PX1273
chla10_rep1 | chla10 | none | TGACCA | PX1273
chla10_rep2 | chla10 | none | ACAGTG | PX1273
chla10_rep3 | chla10 | none | GCCAAT | PX1273
hek293a_rep1 | hek293 | none | CAGATC | PX1273
hek293a_rep2 | hek293 | none | ACTTGA | PX1273
hek293a_rep3 | hek293 | none | GATCAG | PX1273

Now that we have these data, process them downstream with R.