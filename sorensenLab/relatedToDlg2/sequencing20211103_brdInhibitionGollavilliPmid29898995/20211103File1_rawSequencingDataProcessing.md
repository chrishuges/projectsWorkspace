## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a previous publication. Specifically:

"BET bromodomain dependency in EWS/ETS driven Ewing Sarcoma"
GEO: 113604, PMID: 29898995

### Description

I am interested in the CHLA10 RNAseq data where they have knocked down EWS-FLI1 expression. I want to reprocess these to get expression data, .

I am going to use [sraDownloader](https://github.com/s-andrews/sradownloader) to get at these data. I am going to save them in the directory `/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/`. I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). Because the authors originally used Bowtie for the alignment, I will as well. For the index, I created one as described in `..\relatedToDlg2\sequencing20211029_databaseSetupAndSnakemake\20211029File1_downloadingDatabasesAndIndexes.md`. There is a pretty good walkthrough of the alignment process [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/03_align_and_filtering.html). I will then use [samtools](http://www.htslib.org/) and [sambamba](https://lomereiter.github.io/sambamba/) to sort and filter the sam and bam files that are generated, prior to peak calling. I ended up writing this into a script that uses snakemake. Proceed with the analysis as below.

### Data pipeline

First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia
touch sraDataProcessingScript.sh
chmod +x sraDataProcessingScript.sh
touch snakefile
```

Edit the contents of the snakefile to include the text below. I use vim for this.