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

I am going to use SRATools to get at these data. I am going to save them in the directory `/projects/ptx_results/Sequencing/publishedStudies/201709BoulayCellPmid28844694`. There is a nice little tutorial on using SRATools to get raw data [here](https://www.biostars.org/p/111040/). On the proteomics-svr02 system, SRATools is stored in `/projects/ptx_analysis/chughes/software/sratoolkit.2.9.6-1-centos_linux64/bin`. To download the files, I wrote a small shell script called `sraPreFetch.sh`.

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
  fastqCall="$fastqDumpLocation --outdir $rawDataOutputDirectory --split-files /home/chughes/ncbi/public/sra/${i}.sra"
  eval $fastqCall
done
```

I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). I wrote a shell script called `filterAdaptersBbduk.sh` to call this script on each file.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201709BoulayCellPmid28844694/"
bbdukLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbduk.sh"
for i in SRR5217667 SRR5217668 SRR5217669 SRR5217670
do
  echo $i
  bbdukCall="bbdukLocation in1=${rawDataOutputDirectory}${i}_1.fastq.gz in2=${rawDataOutputDirectory}${i}_2.fastq.gz ref=adapters out1=${rawDataOutputDirectory}${i}_1.clean.fastq.gz out2=${rawDataOutputDirectory}${i}_2.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo"
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

So, I think we can safely move on to the alignment. First I need to get my database. I will get it from Ensembl `ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/`. I am going to store it in the directory `/projects/ptx_analysis/chughes/databases/HomoSapiensEnsemblGRCh38_rel102`

```shell
wget ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-102/gff3/homo_sapiens/Homo_sapiens.GRCh38.102.gff3.gz
gunzip *.gz
```

Prepare to do the alignment by creating a shell script to run it for all of the file pairs. This script is called `bbmapAlignmentHuman.sh`.

```shell
#!/bin/bash/
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201709BoulayCellPmid28844694/"
bbmapLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbmap.sh"
referenceLocation="/projects/ptx_analysis/chughes/databases/HomoSapiensEnsemblGRCh38_rel102/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
###########################################
for i in SRR5217667 SRR5217668 SRR5217669 SRR5217670
do
  echo $i
  bbmapCall="$bbmapLocation in1=${rawDataOutputDirectory}${i}_1.clean.fastq.gz in2=${rawDataOutputDirectory}${i}_2.clean.fastq.gz out=${rawDataOutputDirectory}${i}.clean.sam ref=$referenceLocation"
  eval $bbmapCall
done
```

Now we can run [featurecounts](http://bioinf.wehi.edu.au/featureCounts/) on the output. The featurecounts tool is stored in `/projects/ptx_analysis/chughes/software/subread-2.0.1-Linux-x86_64/bin/`. I will execute this using a script called `featureCountsProcessing.sh`. I checked the gtf file and it is in the format required by featurecounts.

```shell
#!/bin/bash/
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201709BoulayCellPmid28844694/"
featureCountsLocation="/projects/ptx_analysis/chughes/software/subread-2.0.1-Linux-x86_64/bin/featureCounts"
gtfLocation="/projects/ptx_analysis/chughes/databases/HomoSapiensEnsemblGRCh38_rel102/Homo_sapiens.GRCh38.102.gtf"
#########################################
for i in SRR5217667 SRR5217668 SRR5217669 SRR5217670
do
  echo $i
  featureCountsCall="-t exon -g gene_id -a $featureCountsLocation -o ${rawDataOutputDirectory}${i}_counts.txt ${rawDataOutputDirectory}${i}.clean.sam"
  eval $featureCountsCall
done
```

That is the end of the re-processing. Now, I will go play around with these data in R in order to see if I can see something of interest for DLG2!

