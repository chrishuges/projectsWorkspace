## Running MAJIQ

This document describes the process of running [MAJIQ](https://biociphers.bitbucket.io/majiq-docs-academic/getting-started-guide/quick-overview.html). I am following the instructions at that link. 

### Data pipeline

First we need to make our annotation file.

```shell
majiq build  -c settings_file.ini /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gff3 -o /home/chughes/databases/projectEwsDlg2/majiqFiles  -j 8

```

I actually ended up following [this document](https://biociphers.bitbucket.io/majiq-docs-academic/gallery/heterogen-vignette.html) more closely. I created a settings file as below.

```shell
# global parameters
[info]
# where do we find the input bam files?
bamdirs=/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/results
# this is used by voila to produce links to UCSC genome browser
genome=GRCh38
# if our experiments are stranded, we can get better results by specifying so
# (a little bit more work if we have mixed strandedness)
strandness=reverse

# we divide input experiments (bam or sj) into "build groups"
# we find these files in bamdirs or sjdirs
[experiments]
# EWS-FLI1 high
EWSFL1high=ATCACG_setA.sorted.bam,ATCACG_setB.sorted.bam,ATCACG_setC.sorted.bam

# EWS-FLI1 low
EWSFL1low=CGATGT_setA.sorted.bam,CGATGT_setB.sorted.bam,CGATGT_setC.sorted.bam
```

Then I ran MAJIQ with the command below.

```shell
 majiq build /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gff3 -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build -c /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/settingsNew.ini
```
