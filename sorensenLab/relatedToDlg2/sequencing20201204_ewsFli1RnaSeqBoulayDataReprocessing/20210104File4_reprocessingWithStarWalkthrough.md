## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a previous publication. Specifically:

"Cancer-Specific Retargeting of BAF Complexes by a Prion-like Domain"
Cell, 2017, Pubmed ID: 28844694, GEO: GSE94278

### Getting the raw data

I am first interested in just the data for the A673 and SKNMC cell lines where it is +/-EWS-FLI1. As far as I can tell, there is just a single replicate run for these samples with a GFP vector used as a control. The data are detailed below:

```
SRX2527813: GSM2472190: RNA-seq in A673 cell line infected with shEWSFLI1; Homo sapiens; RNA-Seq
SRR5217667

SRX2527814: GSM2472191: RNA-seq in A673 cell line infected with shGFP; Homo sapiens; RNA-Seq
SRR5217668

SRX2527815: GSM2472217: RNA-seq in SKNMC cell line infected with shEWSFLI1; Homo sapiens; RNA-Seq
SRR5217669

SRX2527816: GSM2472218: RNA-seq in SKNMC cell line infected with shGFP; Homo sapiens; RNA-Seq
SRR5217670
```

I am going to use SRATools to get at these data. I am going to save them in the directory `/projects/ptx_results/Sequencing/publishedStudies/201411RiggiCancerCellPmid25453903`. There is a nice little tutorial on using SRATools to get raw data [here](https://www.biostars.org/p/111040/). On the proteomics-svr02 system, SRATools is stored in `/projects/ptx_analysis/chughes/software/sratoolkit.2.9.6-1-centos_linux64/bin`. To download the files, I wrote a small shell script called `sraPreFetch.sh`.

```shell
#!/bin/bash
sraToolsLocation="/projects/ptx_analysis/chughes/software/sratoolkit.2.9.6-1-centos_linux64/bin/prefetch"
for i in SRR5217667 SRR5217668 SRR5217669 SRR5217670
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
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201709BoulayCellPmid28844694/"
for i in SRR5217667 SRR5217668 SRR5217669 SRR5217670
do
  echo $i
  fastqCall="$fastqDumpLocation --gzip --outdir $rawDataOutputDirectory --split-files /home/chughes/ncbi/public/sra/${i}.sra"
  eval $fastqCall
done
```

Because I want to look at junction reads, I am going to do an alignment with STAR. There is a good manual [here](https://github.com/alexdobin/STAR) that I am following for the most part. The first thing I want to do is get an index for STAR. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the Salmon index with the command `refgenie pull hg38/star_index`. You don't need to do this every time, just when you need to update your index. 

We can now run STAR. We will do this using a script called `starAlignment.sh`.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201709BoulayCellPmid28844694/"
starLocation="/projects/ptx_analysis/chughes/software/STAR-2.7.1a/STAR-2.7.1a/bin/Linux_x86_64/STAR"
referenceLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/star_index/default/"
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
###########################################
for i in SRR5217667 SRR5217668 SRR5217669 SRR5217670
do
  echo $i
  starCall="$starLocation --runThreadN 12 --genomeDir ${referenceLocation} --readFilesIn ${rawDataOutputDirectory}${i}_1.clean.fastq.gz ${rawDataOutputDirectory}${i}_2.clean.fastq.gz --readFilesCommand zcat --sjdbGTFfile ${annotationLocation} --outFileNamePrefix ${rawDataOutputDirectory}starResults/${i}_ --outSAMtype BAM SortedByCoordinate"
  eval $starCall
done
```

After the alignment is done, I create indexes for the bam files.

```shell
for i in *.bam; do echo $i; /projects/ptx_analysis/chughes/software/samtools-1.9/samtools index $i; done
```

Now we can transfer these files to our local machine and look through the bam alignments in IGV.




