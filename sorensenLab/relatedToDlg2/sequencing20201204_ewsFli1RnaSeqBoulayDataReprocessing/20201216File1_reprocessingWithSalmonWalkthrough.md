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

So now instead of doing an alignment and looking for reads that way, I am going to process the data with Salmon and attempt to look for differential isoform expression. There is a great walkthrough [here](http://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html) that I am following for the most part. The first thing I want to do is get an index for Salmon. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the Salmon index with the command `refgenie pull hg38/salmon_partial_sa_index`. You don't need to do this every time, just when you need to update your index. 

We can now run Salmon. We will do this using a script called `runSalmon.sh`.

```shell
#!/bin/bash
salmonLocation="/projects/ptx_analysis/chughes/software/salmon-1.4/bin/salmon"
indexLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/salmon_partial_sa_index/default/"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201709BoulayCellPmid28844694/"
for i in SRR5217667 SRR5217668 SRR5217669 SRR5217670
do
  echo $i
  salmonCall="$salmonLocation quant -i $indexLocation -l A -1 ${rawDataOutputDirectory}${i}_1.fastq.gz -2 ${rawDataOutputDirectory}${i}_2.fastq.gz --gcBias --validateMappings -o ${rawDataOutputDirectory}${i}_quant"
  eval $salmonCall
done
```






