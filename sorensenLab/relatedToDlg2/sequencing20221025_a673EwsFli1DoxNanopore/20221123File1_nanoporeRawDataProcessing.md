## Processing EwS nanopore data

This document describes the processing of data from nanopore cDNA sequencing of A673 cells at day 0 or day 7 of dox treatment (shRNA against EWS-FLI1).

### Description

I am going to use [minimap2](https://github.com/lh3/minimap2) for this. I am going to follow their standard protocol for mRNA (https://github.com/lh3/minimap2#map-long-splice).

### Data pipeline

Move do the data directory and start the alignment. We used two flow cells, PAM37258 and PAM39767, and each has two libraries multiplexed on it:

PAM37258
noDox = CAAGAAAGTTGTCGGTGTCTTTGTGAC
yesDoc = CTCGATTCCGTTTGTAGTCGTCTGTAC

PAM39767
noDox = CAAGAAAGTTGTCGGTGTCTTTGTGAC
yesDox = CTCGATTCCGTTTGTAGTCGTCTGTAC

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore
```

I am following the protocol detailed [here](https://github.com/epi2me-labs/wf-transcriptomes). Process the data with [pychopper](https://github.com/epi2me-labs/pychopper) prior to alignment. NOTE - I ended up not using pychopper because StringTie2 seems to perform equally as well without using it, and it is slow. So I skipped right to the next step with minimap2.

```shell
##
pychopper -r /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_noDox_report.pdf -u /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_noDox_unclassified.fq -w /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_noDox_rescued.fq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM37258/PAM37258_CAAGAAAGTTGTCGGTGTCTTTGTGAC_pass_fastq_concat.fastq.gz /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_noDox_full_length_output.fq

##
pychopper -r /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_yesDox_report.pdf -u /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_yesDox_unclassified.fq -w /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_yesDox_rescued.fq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM37258/PAM37258_CTCGATTCCGTTTGTAGTCGTCTGTAC_pass_fastq_concat.fastq.gz /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM37258_yesDox_full_length_output.fq

##
pychopper -r /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_noDox_report.pdf -u /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_noDox_unclassified.fq -w /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_noDox_rescued.fq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM39767/PAM39767_CAAGAAAGTTGTCGGTGTCTTTGTGAC_pass_fastq_concat.fastq.gz /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_noDox_full_length_output.fq

##
pychopper -r /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_yesDox_report.pdf -u /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_yesDox_unclassified.fq -w /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_yesDox_rescued.fq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM39767/PAM39767_CTCGATTCCGTTTGTAGTCGTCTGTAC_pass_fastq_concat.fastq.gz /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/PAM39767_yesDox_full_length_output.fq
```

Perform the alignment of the pychopper processed reads using [minimap2](https://github.com/lh3/minimap2). 

```shell
##alignment calls
/home/chughes/softwareTools/minimap2-2.24/minimap2-2.24_x64-linux/minimap2 -ax splice /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.fa /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM37258/PAM37258_CAAGAAAGTTGTCGGTGTCTTTGTGAC_pass_fastq_concat.fastq.gz /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM39767/PAM39767_CAAGAAAGTTGTCGGTGTCTTTGTGAC_pass_fastq_concat.fastq.gz > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.aligned.sam

##
/home/chughes/softwareTools/minimap2-2.24/minimap2-2.24_x64-linux/minimap2 -ax splice /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.fa /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM37258/PAM37258_CTCGATTCCGTTTGTAGTCGTCTGTAC_pass_fastq_concat.fastq.gz /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/raw/IX11100/PAM39767/PAM39767_CTCGATTCCGTTTGTAGTCGTCTGTAC_pass_fastq_concat.fastq.gz > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/yesDox.aligned.sam

##
/home/chughes/softwareTools/minimap2-2.24/minimap2-2.24_x64-linux/minimap2 -ax splice /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.fa /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/raw/PAG30106.fastq > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/results/noDoxParental.aligned.sam
```

Sort the output results in preparation for additional processing.

```shell
##
samtools sort /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.aligned.sam -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.aligned.sorted.bam

##
samtools sort /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/yesDox.aligned.sam -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/yesDox.aligned.sorted.bam

##
samtools sort /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/results/noDoxParental.aligned.sam -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/results/noDoxParental.aligned.sorted.bam
```

Process the files using [StringTie2](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#mix). I am going to use the mix mode that combines short and long read data to improve assignments.

```shell
##
/home/chughes/softwareTools/stringtie-2.2.1/stringtie-2.2.1.Linux_x86_64/stringtie --mix -G /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.transcripts.gtf /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20220725_chrisA673EwsFli1ShrnaPolysomes/results/F115829.sorted.bam /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.aligned.sorted.bam

##
/home/chughes/softwareTools/stringtie-2.2.1/stringtie-2.2.1.Linux_x86_64/stringtie --mix -G /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/yesDox.transcripts.gtf /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20220725_chrisA673EwsFli1ShrnaPolysomes/results/TCGCATTG-ACAGCAAG.sorted.bam /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/yesDox.aligned.sorted.bam

##
/home/chughes/softwareTools/stringtie-2.2.1/stringtie-2.2.1.Linux_x86_64/stringtie --mix -G /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/results/noDoxParental.transcripts.gtf /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20220725_chrisA673EwsFli1ShrnaPolysomes/results/F115829.sorted.bam /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/results/noDoxParental.aligned.sorted.bam
```

I want to run [gffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) on these files in order to determine what the accuracy of StringTie was, and if there are novel isoforms being expressed.

```shell
/home/chughes/softwareTools/gffCompare-0.12.6/gffcompare-0.12.6.Linux_x86_64/gffcompare -R -r /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.transcripts.gtf /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210406_a673ShEwsFli1ParentNanopore/results/noDoxParental.transcripts.gtf /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/yesDox.transcripts.gtf
```




## Creating a FASTA database for proteomics processing

One thing I want to do is create a fasta database with these sequences and re-search some proteomics data to see what comes out. I am going to use bedtools to get the fasta sequences for my transcripts. In order for this tool to work properly, I need to get the data into a [bed12 format](https://bedtools.readthedocs.io/en/latest/content/general-usage.html). So, I am going to use [bedops](https://bedops.readthedocs.io/en/latest/content/installation.html#installation-via-bioconda) for this, specifically the conda version. I created an environment called 'ngs' to work in. I followed the conda instructions from [here](https://genomics.sschmeier.com/ngs-tools/index.html). 

```shell
 conda create -n ngs python=3
 conda activate ngs
 conda install bedops

##run gff2bed
gff2bed < /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.transcripts.filtered.gff > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.transcripts.filtered.bed

##run bedtools again
/home/chughes/softwareTools/bedtools-2.3.0/bedtools getfasta -name -s -split -fo /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.transcripts.filtered.fa  -fi /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.fa -bed /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.transcripts.filtered.bed
```

The above code didn't really work, so I ended up creating the bed12 file in R and running bedtools on it.

```shell
/home/chughes/softwareTools/bedtools-2.3.0/bedtools getfasta -name -s -split -fo /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.transcripts.filtered.fa  -fi /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.fa -bed /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20221025_a673EwsFli1DoxNanopore/results/noDox.transcripts.filtered2.bed
```



