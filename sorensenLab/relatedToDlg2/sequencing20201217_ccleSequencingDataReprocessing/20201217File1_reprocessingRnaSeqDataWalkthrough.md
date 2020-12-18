## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a the CCLE database, or DepMap. The data are all located in SRA: PRJNA523380.

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
  eval "rm *.sra"
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
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/"
bbmapLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbmap.sh"
referenceLocation="/projects/ptx_analysis/chughes/databases/HomoSapiensEnsemblGRCh38_rel102/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
###########################################
for i in SRR8616012 SRR8615497
do
  echo $i
  bbmapCall="$bbmapLocation in1=${rawDataOutputDirectory}${i}_1.clean.fastq.gz in2=${rawDataOutputDirectory}${i}_2.clean.fastq.gz out=${rawDataOutputDirectory}${i}.clean.sam ref=$referenceLocation trimreaddescriptions=t"
  eval $bbmapCall
done
```

Now we can run [featurecounts](http://bioinf.wehi.edu.au/featureCounts/) on the output. The featurecounts tool is stored in `/projects/ptx_analysis/chughes/software/subread-2.0.1-Linux-x86_64/bin/`. I will execute this using a script called `featureCountsProcessing.sh`. I checked the gtf file and it is in the format required by featurecounts.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/"
featureCountsLocation="/projects/ptx_analysis/chughes/software/subread-2.0.1-Linux-x86_64/bin/featureCounts"
gtfLocation="/projects/ptx_analysis/chughes/databases/HomoSapiensEnsemblGRCh38_rel102/Homo_sapiens.GRCh38.102.gtf"
#########################################
for i in SRR8616012 SRR8615497
do
  echo $i
  featureCountsCall="$featureCountsLocation -t exon -g gene_id -a $gtfLocation -o ${rawDataOutputDirectory}${i}_counts.txt ${rawDataOutputDirectory}${i}.clean.sam"
  eval $featureCountsCall
done
```

I want to create bam files so I can look at the reads aligning with DLG2. To do this, I will use samtools in a script called `bamFileCreation.sh`. 

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"
#########################################
for i in SRR8616012 SRR8615497
do
  echo $i
  bamCreateCall="$samtoolsLocation view -b ${rawDataOutputDirectory}${i}.clean.sam > ${rawDataOutputDirectory}${i}.clean.bam"
  eval $bamCreateCall
  bamSortCall="$samtoolsLocation sort ${rawDataOutputDirectory}${i}.clean.bam -o ${rawDataOutputDirectory}${i}.sorted.bam"
  eval $bamSortCall
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}${i}.sorted.bam"
  eval $bamIndexCall
done
```

These files are pretty large, so I will just extract DLG2 at the end. The location for this in GRCh38 is chr11:83455009-85628535, so we will just do the whole chromosome 11. I will do this using a script called `chrExtraction.sh`. To get the chromosome name format, I used the command `samtools view -H your.bam` to see if they had the 'chr' prefix. They didn't.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"
#########################################
for i in SRR8616012 SRR8615497
do
  echo $i
  chr11ExtractCall="$samtoolsLocation view -b ${rawDataOutputDirectory}${i}.sorted.bam 11 > ${rawDataOutputDirectory}${i}.chr11.bam"
  eval $chr11ExtractCall
  chr11IndexCall="$samtoolsLocation index ${rawDataOutputDirectory}${i}.chr11.bam"
  eval $chr11IndexCall
done
```

The last thing I did was to delete the sam files. They are massive, so I didn't want to keep them. They can be re-generated easily enough using the above code, if necessary. That is the end of the re-processing. Now, I will go play around with these data in R in order to see if I can see something of interest for DLG2!

