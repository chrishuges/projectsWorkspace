## Reprocessing some RIPseq data

This document describes the reprocessing of some RIPseq data from a previous publication. Specifically:

"YB-3 substitutes YB-1 in global mRNA binding"
RNA Biology, 2020, Pubmed ID: 31944153, GEO: GSE130781

### Getting the raw data

I am interested in the RIPseq data for YB-1. The data are detailed below:

```
GSM3753523: YB-1 RIP-Seq wt PolyA (IP_YB1_polyA), run 1; Homo sapiens; RIP-Seq
SRR9019705

GSM3753524: YB-1 RIP-Seq wt PolyA (IP_YB1_polyA), run 2; Homo sapiens; RIP-Seq
SRR9019706

GSM3753525: YB-1 RIP-Seq wt Total (IP_YB1_H), run 1; Homo sapiens; RIP-Seq
SRR9019707

GSM3753526: YB-1 RIP-Seq wt Total (IP_YB1_H), run 2; Homo sapiens; RIP-Seq
SRR9019708

GSM3753540: YB-1 RIP-Seq wt Total (Bethyl1_YB1), run 1; Homo sapiens; RIP-Seq
SRR9019722

GSM3753541: YB-1 RIP-Seq wt Total (Bethyl1_YB1), run 2; Homo sapiens; RIP-Seq
SRR9019723
```

I am going to use [SRAExplorer](https://sra-explorer.info/#) to get at these data. I am going to save them in the directory `/projects/ptx_results/Sequencing/publishedStudies/202004LyabinRnaBiologyPmid31944153/yb1RipSeq`. I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).  We can now run Salmon as well as STAR. 

```shell
#!/usr/bin/env bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR901/008/SRR9019708/SRR9019708.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR901/002/SRR9019722/SRR9019722.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR901/003/SRR9019723/SRR9019723.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR901/005/SRR9019705/SRR9019705.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR901/007/SRR9019707/SRR9019707.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR901/006/SRR9019706/SRR9019706.fastq.gz

###########################################
bbdukLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbduk.sh"
salmonLocation="/projects/ptx_analysis/chughes/software/salmon-1.4/bin/salmon"
indexLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/salmon_partial_sa_index/default/"
starLocation="/projects/ptx_analysis/chughes/software/STAR-2.7.1a/STAR-2.7.1a/bin/Linux_x86_64/STAR"
referenceLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/star_index/default/"
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/202004LyabinRnaBiologyPmid31944153/yb1RipSeq/"

###########################################
for i in SRR9019705 SRR9019706 SRR9019707 SRR9019708 SRR9019722 SRR9019723
do
  echo $i
  ##
  bbdukCall="$bbdukLocation in=${rawDataOutputDirectory}${i}.fastq.gz ref=adapters out=${rawDataOutputDirectory}${i}.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1"
  eval $bbdukCall
  ##
  salmonCall="$salmonLocation quant -i $indexLocation -l A -1 ${rawDataOutputDirectory}${i}.clean.fastq.gz --gcBias --validateMappings -o ${rawDataOutputDirectory}${i}_quant"
  eval $salmonCall
  ##
  starCall="$starLocation --runThreadN 12 --genomeDir ${referenceLocation} --readFilesIn ${rawDataOutputDirectory}${i}.clean.fastq.gz --readFilesCommand zcat --sjdbGTFfile ${annotationLocation} --outFileNamePrefix ${rawDataOutputDirectory}starResults/${i}_ --outSAMtype BAM Unsorted"
  eval $starCall
  ##
  bamSortCall="$samtoolsLocation sort ${rawDataOutputDirectory}/starResults/${i}_Aligned.out.bam -o ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamSortCall
  ##
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamIndexCall
done
```

I think for the rest of the visualization, we will carry this out in R.
