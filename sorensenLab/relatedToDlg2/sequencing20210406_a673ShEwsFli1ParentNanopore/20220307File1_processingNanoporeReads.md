## Processing the nanopore read data

This document describes the processing of the nanopore data acquired for A673 cells.

## Workflow

I am going to use [minimap2](https://github.com/lh3/minimap2) to do the alignment. After, I will use [Freddie](https://github.com/vpc-ccg/freddie) to do isoform detection. I am just following the directions as detailed on these websites. First, I will perform the alignment.

```shell
#!/bin/bash

#tool locations
minimapLocation="/home/chughes/softwareTools/minimap2-2.24/minimap2"
fastaLocation="/home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.fa"
samtoolsLocation="/home/chughes/softwareTools/samtools-1.12/samtools"
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore"

#process the data
eval $minimapLocation -a -x splice -t 6 $fastaLocation $workingDirectory/PAG30106_pass_concat.fastq.gz > $workingDirectory/PAG30106.unsorted.sam
eval samtools sort $workingDirectory/PAG30106.unsorted.sam -o $workingDirectory/PAG30106.sorted.bam
eval samtools index $workingDirectory/PAG30106.sorted.bam
```

Now we need to run the Freddie tool. There are different steps for this, and as I mention above, I just follow what they say on the Freddie GitHub page.

```shell
conda env create -f envs/freddie.yml
conda activate freddie

##the fastq needed to be decompressed for this to work
py/freddie_split.py --reads /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/PAG30106.fastq --bam /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/PAG30106.sorted.bam --outdir /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/freddieSplit -t 6

##
py/freddie_segment.py -s /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/freddieSplit --outdir  /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/freddieSegment -t 6

##
export GRB_LICENSE_FILE="/home/chughes/softwareTools/gurobi.lic"
py/freddie_cluster.py --segment-dir /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/freddieSegment --outdir /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/freddieCluster

##
py/freddie_isoforms.py --split-dir /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/freddieSplit --cluster-dir /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/freddieCluster --output /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/freddieIsoforms/PAG30106.isoforms.gtf -t 6
```










