## Preparing indexes

This document describes the process of preparing indexes for the DLG2 project.


## Bowtie2 index

First we will prepare the index for Bowtie2. I am following their instructions as detailed [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer). Version 38 was current at the time of download (Oct. 29, 2021). First I did the download with a shell script as detailed below and wrote a log so I could track the version.

```shell
#!/bin/bash
baseDirectory="/home/chughes/databases/projectEwsDlg2/"

##download the database files
eval cd ${baseDirectory}baseGenomeFiles
eval wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
eval gzip -d GRCh38.primary_assembly.genome.fa.gz
eval mv GRCh38.primary_assembly.genome.fa genome.fa
eval wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
eval gzip -d gencode.v38.annotation.gtf
eval mv gencode.v38.annotation.gtf genome.gtf
eval wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
eval wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
```

Now we can build the bowtie2 index. This was run as a shell script and logged as above.

```shell
#!/bin/bash

##locations and software tools
baseDirectory="/home/chughes/databases/projectEwsDlg2/"
bowtie2Location="/home/chughes/softwareTools/bowtie2-2.4.4/"

##build the index
eval mkdir ${baseDirectory}bowtie2Index
eval cd ${baseDirectory}bowtie2Index
eval ${bowtie2Location}bowtie2-build ${baseDirectory}baseGenomeFiles/genome.fa human_genome
```

## HISAT2 index

Now we will prepare the index for HISAT2. I am following their instructions as detailed [here](http://daehwankimlab.github.io/hisat2/howto/). Version 104 was current at the time of download (Oct. 29, 2021).

```shell
#!/bin/bash
baseDirectory="/home/chughes/databases/projectEwsDlg2/"
hisat2Location="/home/chughes/softwareTools/hisat2-2.2.1/"

eval mkdir ${baseDirectory}hisat2Index
eval cd ${baseDirectory}hisat2Index
eval ${hisat2Location}hisat2_extract_splice_sites.py ${baseDirectory}baseGenomeFiles/genome.gtf > genome.ss
eval ${hisat2Location}hisat2_extract_exons.py ${baseDirectory}baseGenomeFiles/genome.gtf > genome.exon
eval ${hisat2Location}hisat2-build -p 16 ${baseDirectory}baseGenomeFiles/genome.fa genome
eval ${hisat2Location}hisat2-build -p 16 --exon ${baseDirectory}hisat2Index/genome.exon --ss ${baseDirectory}hisat2Index/genome.ss ${baseDirectory}baseGenomeFiles/genome.fa genome_tran
```

## SALMON index

The main Salmon page is [here](https://github.com/COMBINE-lab/salmon). I want to make a decoy-aware transcriptome using the whole genome as a reference. Details for this are [here](https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/).

```shell
mkdir salmonIndex
cd baseGenomeFiles
grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat gencode.v38.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz

## script for building library
#!/bin/bash
baseDirectory="/home/chughes/databases/projectEwsDlg2/"
salmonLocation="/home/chughes/softwareTools/salmon-1.5.2/bin/"

eval cd ${baseDirectory}salmonIndex
eval ${salmonLocation}salmon index -t ${baseDirectory}baseGenomeFiles/gentrome.fa.gz -d ${baseDirectory}baseGenomeFiles/decoys.txt -p 8 -i salmon_index --gencode
```

For now, I don't have any other indexes to make, but I may revisit and add here later on.