## Processing some RNAseq data

This document details analysis of the CHLA cell line RNAseq data with Salmon. There is a great walkthrough [here](http://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html) that I am following for the most part. The first thing I want to do is get an index for Salmon. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the Salmon index with the command `refgenie pull hg38/salmon_partial_sa_index`. You don't need to do this every time, just when you need to update your index.

The raw data are stored here: `/projects/analysis/analysis34/PX1936/HG2HKCCX2_8/`. Transfer the data to our own server.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/2021/sequencing20210312_mg633Ybx1Guc6ApaDeepSeq/"
rawDataInputDirectory="/projects/analysis/analysis34/PX1936/HG2HKCCX2_8/"

######################################
for i in ATCACG CGATGT
do
  echo $i
  eval "scp ${rawDataInputDirectory}PX1936_${i}/150bp/HG2HKCCX2_8_1_${i}_150bp.concat.fastq.gz ${rawDataOutputDirectory}${i}_1.fastq.gz"
  eval "scp ${rawDataInputDirectory}PX1936_${i}/150bp/HG2HKCCX2_8_2_${i}_150bp.concat.fastq.gz ${rawDataOutputDirectory}${i}_2.fastq.gz"
done
```

We can now run Salmon as well as STAR. We will do this using a script called `runSalmonAndStar.sh`. We will also preprocess the data with UMI processing tools found [here](https://github.com/weng-lab/umitools).

```shell
#!/bin/bash
bbdukLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbduk.sh"
salmonLocation="/projects/ptx_analysis/chughes/software/salmon-1.4/bin/salmon"
indexLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/salmon_partial_sa_index/default/"
starLocation="/projects/ptx_analysis/chughes/software/STAR-2.7.1a/STAR-2.7.1a/bin/Linux_x86_64/STAR"
referenceLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/star_index/default/"
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/2021/sequencing20210312_mg633Ybx1Guc6ApaDeepSeq/"

###########################################
for i in ATCACG CGATGT
do
  echo $i
  ##
  umiToolsCall="umitools reformat_fastq -l ${rawDataOutputDirectory}${i}_1.fastq.gz -r ${rawDataOutputDirectory}${i}_2.fastq.gz -L ${rawDataOutputDirectory}${i}_1.umi.fastq.gz -R ${rawDataOutputDirectory}${i}_2.umi.fastq.gz"
  eval $umiToolsCall
  ##
  bbdukCall="$bbdukLocation in1=${rawDataOutputDirectory}${i}_1.umi.fastq.gz in2=${rawDataOutputDirectory}${i}_2.umi.fastq.gz ref=adapters out1=${rawDataOutputDirectory}${i}_1.clean.fastq.gz out2=${rawDataOutputDirectory}${i}_2.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo"
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
mg633Parent | parent | none | CGTGAT (ATCACG) | PX1936
mg633Ybx1Guc6 | guc6 | none | ACATCG (CGATGT) | PX1936

Now that we have these data, process them downstream with R.
