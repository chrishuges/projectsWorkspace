## Processing some RNAseq data

This document details analysis of the A673 timecourse RNAseq data with Salmon. There is a great walkthrough [here](http://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html) that I am following for the most part. The first thing I want to do is get an index for Salmon. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the Salmon index with the command `refgenie pull hg38/salmon_partial_sa_index`. You don't need to do this every time, just when you need to update your index.

The raw data are stored here: `/projects/analysis/analysis34/`. Transfer the data to our own server.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/2021/sequencing20210421_a673EwsFli1Timecourse/"
rawDataInputDirectory="/projects/analysis/analysis34/"
bbdukLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbduk.sh"
salmonLocation="/projects/ptx_analysis/chughes/software/salmon-1.4/bin/salmon"
indexLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/alias/hg38/salmon_partial_sa_index/default/"
starLocation="/projects/ptx_analysis/chughes/software/STAR-2.7.1a/STAR-2.7.1a/bin/Linux_x86_64/STAR"
referenceLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/alias/hg38/star_index/default/"
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/alias/hg38/gencode_gtf/default/hg38.gtf.gz"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"

######################################
for i in ACAGTG ACTTGA ATCACG CAGATC CGATGT GCCAAT TGACCA TTAGGC
do
  echo $i
  eval "scp ${rawDataInputDirectory}PX1955/HG3N5CCX2_6/PX1955_${i}/150bp/HG3N5CCX2_6_1_${i}_150bp.concat.fastq.gz ${rawDataOutputDirectory}${i}_setA_1.fastq.gz"
  eval "scp ${rawDataInputDirectory}PX1955/HG3N5CCX2_6/PX1955_${i}/150bp/HG3N5CCX2_6_2_${i}_150bp.concat.fastq.gz ${rawDataOutputDirectory}${i}_setA_2.fastq.gz"
  eval "scp ${rawDataInputDirectory}PX1956/HG3N5CCX2_7/PX1956_${i}/150bp/HG3N5CCX2_7_1_${i}_150bp.concat.fastq.gz ${rawDataOutputDirectory}${i}_setB_1.fastq.gz"
  eval "scp ${rawDataInputDirectory}PX1956/HG3N5CCX2_7/PX1956_${i}/150bp/HG3N5CCX2_7_2_${i}_150bp.concat.fastq.gz ${rawDataOutputDirectory}${i}_setB_2.fastq.gz"
  eval "scp ${rawDataInputDirectory}PX1957/HG3N5CCX2_8/PX1957_${i}/150bp/HG3N5CCX2_8_1_${i}_150bp.concat.fastq.gz ${rawDataOutputDirectory}${i}_setC_1.fastq.gz"
  eval "scp ${rawDataInputDirectory}PX1957/HG3N5CCX2_8/PX1957_${i}/150bp/HG3N5CCX2_8_2_${i}_150bp.concat.fastq.gz ${rawDataOutputDirectory}${i}_setC_2.fastq.gz"
done

######################################
for i in ACAGTG_setA ACTTGA_setA ATCACG_setA CAGATC_setA CGATGT_setA GCCAAT_setA TGACCA_setA TTAGGC_setA ACAGTG_setB ACTTGA_setB ATCACG_setB CAGATC_setB CGATGT_setB GCCAAT_setB TGACCA_setB TTAGGC_setB ACAGTG_setC ACTTGA_setC ATCACG_setC CAGATC_setC CGATGT_setC GCCAAT_setC TGACCA_setC TTAGGC_setC
do
  ###
  umiToolsCall="umitools reformat_fastq -l ${rawDataOutputDirectory}${i}_1.fastq.gz -r ${rawDataOutputDirectory}${i}_2.fastq.gz -L ${rawDataOutputDirectory}${i}_1.umi.fastq.gz -R ${rawDataOutputDirectory}${i}_2.umi.fastq.gz"
  eval $umiToolsCall
  ###
  bbdukCall="$bbdukLocation in1=${rawDataOutputDirectory}${i}_1.umi.fastq.gz in2=${rawDataOutputDirectory}${i}_2.umi.fastq.gz ref=adapters out1=${rawDataOutputDirectory}${i}_1.clean.fastq.gz out2=${rawDataOutputDirectory}${i}_2.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo"
  eval $bbdukCall
  ##
  salmonCall="$salmonLocation quant -i $indexLocation -l A -1 ${rawDataOutputDirectory}${i}_1.clean.fastq.gz -2 ${rawDataOutputDirectory}${i}_2.clean.fastq.gz --gcBias --validateMappings -o ${rawDataOutputDirectory}${i}_quant"
  eval $salmonCall
  ##
  #starCall="$starLocation --runThreadN 12 --genomeDir ${referenceLocation} --readFilesIn ${rawDataOutputDirectory}${i}_1.clean.fastq.gz ${rawDataOutputDirectory}${i}_2.clean.fastq.gz --readFilesCommand zcat --sjdbGTFfile ${annotationLocation} --outFileNamePrefix ${rawDataOutputDirectory}${i}_ --outSAMtype BAM Unsorted"
  #eval $starCall
  ##
  #bamSortCall="$samtoolsLocation sort ${rawDataOutputDirectory}${i}_Aligned.out.bam -o ${rawDataOutputDirectory}${i}.sorted.bam"
  #eval $bamSortCall
  ##
  #bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}${i}.sorted.bam"
  #eval $bamIndexCall
done
```

For downstream analysis, file file layout is:

sampleName | cell | treatment | barcode | experiment
| --- | --- | --- | --- | --- |
day0_rep1 | a673 | none | ATCACG | PX1955
day0_rep2 | a673 | none | ATCACG | PX1956
day0_rep3 | a673 | none | ATCACG | PX1957
day7_rep1 | a673 | none | CGATGT | PX1955
day7_rep2 | a673 | none | CGATGT | PX1956
day7_rep3 | a673 | none | CGATGT | PX1957
day9_rep1 | a673 | none | TTAGGC | PX1955
day9_rep2 | a673 | none | TTAGGC | PX1956
day9_rep3 | a673 | none | TTAGGC | PX1957
day10_rep1 | a673 | none | TGACCA | PX1955
day10_rep2 | a673 | none | TGACCA | PX1956
day10_rep3 | a673 | none | TGACCA | PX1957
day11_rep1 | a673 | none | ACAGTG | PX1955
day11_rep2 | a673 | none | ACAGTG | PX1956
day11_rep3 | a673 | none | ACAGTG | PX1957
day14_rep1 | a673 | none | GCCAAT | PX1955
day14_rep2 | a673 | none | GCCAAT | PX1956
day14_rep3 | a673 | none | GCCAAT | PX1957
day17_rep1 | a673 | none | CAGATC | PX1955
day17_rep2 | a673 | none | CAGATC | PX1956
day17_rep3 | a673 | none | CAGATC | PX1957
day22_rep1 | a673 | none | ACTTGA | PX1955
day22_rep2 | a673 | none | ACTTGA | PX1956
day22_rep3 | a673 | none | ACTTGA | PX1957

Now that we have these data, process them downstream with R.
