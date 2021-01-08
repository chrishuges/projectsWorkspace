## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from the [Human Protein Atlas](https://www.proteinatlas.org/). I am specifically interested in the brain data. They have deposited the raw files in ENA with the accession [ERP006650](https://www.ebi.ac.uk/ena/browser/view/PRJEB6971).

### Getting the raw data

The first thing I need to do is download the data. I am only interested in brain tissue and there are only a few accessions for this. The main page listed above didn't seem to have all the files, so I ended up going [here](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2836/samples/?full=true&s_page=11&s_pagesize=25&s_sortby=col_31&s_sortorder=ascending) to see them all. The files I am interested in are:

```
brain_3c
ERR315455

brain_3b
ERR315477

brain_a
ERR315432
```

The tissue annotation for all of these appears to be cerebral cortex. Since these files are hosted on ENA, I can just use wget to download them to my machine. I will do this using a shell script called `enaDataDownload.sh`.

```shell
#!/bin/bash
ftpLocation="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR315/"
downloadDirectory="/projects/ptx_results/Sequencing/publishedStudies/humanProteinAtlasRnaSequencingData"
for i in ERR315455 ERR315477 ERR315432
do
  echo $i
  eval "cd ${downloadDirectory}"
  enaDownloadRead1="wget ${ftpLocation}${i}/${i}_1.fastq.gz"
  enaDownloadRead2="wget ${ftpLocation}${i}/${i}_2.fastq.gz"
  eval $enaDownloadRead1
  printf "downloading read 1 complete\n"
  eval $enaDownloadRead2
  printf "downloading read 1 complete\n"
done
```

I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). I wrote a shell script called `filterAdaptersBbduk.sh` to call this script on each file.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/humanProteinAtlasRnaSequencingData/"
bbdukLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbduk.sh"
for i in ERR315455 ERR315477 ERR315432
do
  echo $i
  bbdukCall="$bbdukLocation in1=${rawDataOutputDirectory}${i}_1.fastq.gz in2=${rawDataOutputDirectory}${i}_2.fastq.gz ref=adapters out1=${rawDataOutputDirectory}${i}_1.clean.fastq.gz out2=${rawDataOutputDirectory}${i}_2.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo"
  eval $bbdukCall
done
```

It looks like basically all of the reads pass the filtering criteria, just a few bases trimmed on the ends. So, that is good. Here is an example of the output for one of the files:

```
Input is being processed as paired
Started output streams: 0.100 seconds.
Processing time:                169.503 seconds.

Input:                          56933296 reads          5750262896 bases.
KTrimmed:                       8115170 reads (14.25%)  252277888 bases (4.39%)
Trimmed by overlap:             2746070 reads (4.82%)   16313526 bases (0.28%)
Total Removed:                  441852 reads (0.78%)    268591414 bases (4.67%)
Result:                         56491444 reads (99.22%)         5481671482 bases (95.33%)

Time:                           169.731 seconds.
Reads Processed:      56933k    335.43k reads/sec
Bases Processed:       5750m    33.88m bases/sec
```

Because I want to look at junction reads, I am going to do an alignment with STAR. There is a good manual [here](https://github.com/alexdobin/STAR) that I am following for the most part. The first thing I want to do is get an index for STAR. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the Salmon index with the command `refgenie pull hg38/star_index`. You don't need to do this every time, just when you need to update your index. 

We can now run STAR. We will do this using a script called `starAlignment.sh`.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/humanProteinAtlasRnaSequencingData/"
starLocation="/projects/ptx_analysis/chughes/software/STAR-2.7.1a/STAR-2.7.1a/bin/Linux_x86_64/STAR"
referenceLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/star_index/default/"
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"
###########################################
for i in ERR315455 ERR315477 ERR315432
do
  echo $i
  starCall="$starLocation --runThreadN 12 --genomeDir ${referenceLocation} --readFilesIn ${rawDataOutputDirectory}${i}_1.clean.fastq.gz ${rawDataOutputDirectory}${i}_2.clean.fastq.gz --readFilesCommand zcat --sjdbGTFfile ${annotationLocation} --outFileNamePrefix ${rawDataOutputDirectory}starResults/${i}_ --outSAMtype BAM Unsorted"
  eval $starCall
  bamSortCall="$samtoolsLocation sort ${rawDataOutputDirectory}/starResults/${i}_Aligned.out.bam -o ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamSortCall
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamIndexCall
done
```
Now we can transfer these files to our local machine and look through the bam alignments in IGV.
