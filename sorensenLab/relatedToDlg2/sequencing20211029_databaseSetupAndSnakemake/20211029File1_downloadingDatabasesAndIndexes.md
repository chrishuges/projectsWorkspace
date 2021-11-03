## Preparing indexes

This document describes the process of preparing indexes for the DLG2 project.


## Bowtie2 index

First we will prepare the index for Bowtie2. I am following their instructions as detailed [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer). Version 104 was current at the time of download (Oct. 29, 2021). First I did the download with a shell script as detailed below and wrote a log so I could track the version.

```shell
#!/bin/bash
baseDirectory="/home/chughes/databases/projectEwsDlg2/"

##download the database files
eval cd ${baseDirectory}baseGenomeFiles
eval wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
eval gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
eval mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa
eval wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz  
eval gzip -d Homo_sapiens.GRCh38.104.gtf.gz
eval mv Homo_sapiens.GRCh38.104.gtf genome.gtf
```

Now we can build the bowtie2 index. This was run as a shell script and logged as above.

```shell
#!/bin/bash

##locations and software tools
baseDirectory="/home/chughes/databases/projectEwsDlg2/"
bowtie2Location="/home/chughes/softwareTools/bowtie2-2.4.4/"

##build the index
eval cd ${baseDirectory}bowtie2Index
eval ${bowtie2Location}bowtie2-build ${baseDirectory}baseGenomeFiles/genome.fa human_genome
```

## HISAT2 index

Now we will prepare the index for HISAT2. I am following their instructions as detailed [here](http://daehwankimlab.github.io/hisat2/howto/). Version 104 was current at the time of download (Oct. 29, 2021).

```shell
#!/bin/bash
workingDirectory="/home/chughes/databases/projectEwsDlg2/"
hisat2Location="/home/chughes/softwareTools/hisat2-2.2.1/"

eval cd /home/chughes/databases/projectEwsDlg2/hisat2Index
eval ${hisat2Location}hisat2_extract_splice_sites.py ${workingDirectory}baseGenomeFiles/genome.gtf > genome.ss
eval ${hisat2Location}hisat2_extract_exons.py ${workingDirectory}baseGenomeFiles/genome.gtf > genome.exon
eval ${hisat2Location}hisat2-build -p 16 ${workingDirectory}baseGenomeFiles/genome.fa genome
```

For now, I don't have any other indexes to make, but I may revisit and add here later on.