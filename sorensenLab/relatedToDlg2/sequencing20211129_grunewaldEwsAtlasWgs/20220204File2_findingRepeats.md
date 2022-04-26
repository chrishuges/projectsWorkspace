## Finding repeat regions in hg38

This document describes the analysis of the human genome to find GGAA repeats.

### Description

In order to do this I am going to use the tool [Tandem Repeats Finder](https://github.com/Benson-Genomics-Lab/TRF). You can pass it the whole genome, but I find the output format a bit unweildy, especially because the Gencode database we are using has many more 'sequences' than just the normal chromosomes. Really, all I am interested in for this analysis is the primary chromosomes. So, I am going to use [samtools faidx](http://www.htslib.org/doc/samtools-faidx.html) to extract them individually and process them with repeat finder.

We can do this using a simple shell script.

```shell
#!/bin/bash

#set the directory location where you want to work
baseDirectory="/projects/ptx_analysis/chughes/databases/projectEwsDlg2/baseGenomeFiles"
samtools="/gsc/software/linux-x86_64-centos7/samtools-1.14/bin/samtools"
tandemRepeatFinder="/projects/ptx_analysis/chughes/softwareTools/tandemRepeatFinder-4.0.9/tandemRepeatFinder"

#migrate to the base directory
eval cd ${baseDirectory}/tandemRepeats

#create and process the individual chromosome files
for i in chr{1..23} chrX chrY
do
  printf "extracting ${i} from genome"
  eval $samtools faidx -o ${baseDirectory}/individualChromosomes/${i}.fa ${baseDirectory}/genome.fa ${i}
  printf "processing tandem repeats from ${i}"
  eval ${tandemRepeatFinder} ${baseDirectory}/individualChromosomes/${i}.fa 2 7 7 80 10 50 500 -f -h -d
done 
```

Test command

/projects/ptx_analysis/chughes/softwareTools/tandemRepeatFinder-4.0.9/tandemRepeatFinder /projects/ptx_analysis/chughes/databases/projectEwsDlg2/baseGenomeFiles/individualChromosomes/testFile.fa 2 7 7 80 10 50 500 -f -h -d

I processed the output from the tandem repeat analysis in order to make a bed file that I can use with bedtools closest.

```shell
#for stranded processing
/home/chughes/softwareTools/bedtools-2.3.0/bedtools closest -s -D "b" -a /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211129_grunewaldEwsAtlasWgs/dataset_msatRepeats.bed -b /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211129_grunewaldEwsAtlasWgs/dataset_gtfGenesOnly.bed > dataset_hg38MsatClosestGenes.tsv

#for unstranded processing
/home/chughes/softwareTools/bedtools-2.3.0/bedtools closest -D "b" -a /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211129_grunewaldEwsAtlasWgs/dataset_msatRepeats.bed -b /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211129_grunewaldEwsAtlasWgs/dataset_gtfGenesOnly.bed > dataset_hg38MsatClosestGenesUnstranded.tsv
```