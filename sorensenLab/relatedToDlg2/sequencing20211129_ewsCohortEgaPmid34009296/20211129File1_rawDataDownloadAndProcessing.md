## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a previous publication. Specifically:

"ERG transcription factors have a splicing regulatory function involving RBFOX2 that is altered in the EWS-FLI1 oncogenic fusion"
Nucleic Acids Research, 2021, Pubmed ID: 34009296, EGA: EGAS00001003333

### Description

I am interested in the RNAseq data across the entire cohort of EwS patients. These data are controlled access, so I had to request to be allowed to work with the data. These are the information provided for the RNAseq provided in the manuscript:

EWS-FLI1-dependent splicing events were identified by comparing doxycycline-untreated versus doxycycline-treated A673/TR/shEF cells. ERG-, FLI1 or RBFOX2-dependent splicing events were identified by comparing cells treated with siCTRL or with the corresponding siRNA. RNA was isolated as described above and sample integrity was evaluated using a Bioanlayzer instrument (Agilent). Only samples with RNA Integrity Number above 9 were used. Libraries were performed using the TruSeq Stranded mRNA Library Preparation Kit. Equimolar pools of libraries were sequenced on a Illumina HiSeq 2500 machine using paired-end reads and High Output run mode allowing 200 million raw reads per sample. Raw reads were mapped to the human reference genome hg19 using the STAR aligner (v.2.5.0a) (23). PCR-duplicated reads and low mapping quality reads (MQ<20) were removed using Picard tools and SAMtools, respectively. We next used rMATS (v3.0.9) (24), an event-based tool, to identify differentially spliced events using RNA-seq data. Five distinct alternative splicing events were analyzed using rMATS: skipped exons (SE), alternative 3′ splice sites (A3SS), alternative 5′ splice sites (A5SS), mutually exclusive exons (MXE) and retained introns (RI). Briefly, rMATS uses a count-based model, to calculate percent of spliced-in (PSI) value among replicates, using both spliced reads and reads that mapped to the exon body. We used three different thresholds to identify differentially spliced events between two groups: each splicing event has to be (i) supported by at least 15 unique reads, (ii) |ΔPSI| > 10%; (iii) FDR < 0.05. For RT-qPCR validation, we choose a set of alternative splicing events with ΔPSI values spread across a wide range. Events were selected from a list of events with >50 reads supporting the event in at least one condition (to allow detection by RT-PCR). Gene expression analysis was performed as follows: aligned reads were counted using htseq-count v.0.6.1p1 (25) and normalized according to the DESeq size factors method (26). We used fold change ≥2 and FDR <0.05 as the determination of differentially expressed genes (DEG).

I am going to use [sraDownloader](https://github.com/s-andrews/sradownloader) to get at these data.

I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).

I will use [STAR](https://github.com/alexdobin/STAR) for the alignment. I am following the instructions in the STAR [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

I will use [Salmon](https://github.com/COMBINE-lab/salmon) to do read quantification. There is a great info page for this tool [here](https://salmon.readthedocs.io/en/latest/).

For library construction, they don't mention anything about whether the libraries are stranded or not, so we will assume they aren't.

I will then use [samtools](http://www.htslib.org/) to prepare the files before using [FeatureCounts](http://subread.sourceforge.net/featureCounts.html) (part of [Subread](http://subread.sourceforge.net/)) for quantification.

I will then use the bamCoverage function from [deeptools](https://deeptools.readthedocs.io/en/develop/) to get coverage estimates across the transcriptome.

### Data pipeline

First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211129_ewsCohortEgaPmid34009296
touch sraDataProcessingScript.sh
chmod +x sraDataProcessingScript.sh
touch snakefile
```

Edit the contents of the snakefile to include the text below. I use vim for this.

```python
"""
Author: Christopher Hughes
Affiliation: BCCRC
Aim: Workflow for RNA-Seq data
Date: 20211105
"""

###############################
#working directory
BASE_DIR = "/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211210_soleMscToEwsPmid34341072/rnaSeq"


###############################
#locations of tools we will use
BBDUK = "/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
SALMON = "/home/chughes/softwareTools/salmon-1.5.2/bin/salmon"
#HISAT2 = "/home/chughes/softwareTools/hisat2-2.2.1/hisat2"
STAR = "/home/chughes/softwareTools/STAR-2.7.9a/bin/Linux_x86_64/STAR"
SAMTOOLS="/home/chughes/softwareTools/samtools-1.12/samtools"
#SAMBAMBA="/home/chughes/softwareTools/sambamba-0.8.1/sambamba"
FEATURECOUNTS="/home/chughes/softwareTools/subread-2.0.3/bin/featureCounts"
BAMCOVERAGE="/home/chughes/virtualPython368/bin/bamCoverage"


###############################
#locations of our index files
DATABASE_DIR = "/home/chughes/databases/projectEwsDlg2"
STARINDEX = DATABASE_DIR + "/starIndex"
SALMONINDEX = DATABASE_DIR + "/salmonIndex/salmon_index"
GTF = DATABASE_DIR + "/baseGenomeFiles/genome.gtf"
FASTA = DATABASE_DIR + "/baseGenomeFiles/genome.fa"


###############################
#get our list of sample files
SAMPLES, = glob_wildcards("raw/{smp}_1.fastq.gz")
print("There are " + str(len(SAMPLES)) + " total samples to be processed.")
for smp in SAMPLES:
  print("Sample " + smp + " will be processed")


###############################
#processing workflow
rule all:
    input:
      expand("results/{smp}.counts.txt", smp = SAMPLES),
      expand("results/{smp}.sorted.bam.bai", smp = SAMPLES),
      expand("results/{smp}.sorted.bw", smp = SAMPLES),
      expand("quants/{smp}/quant.sf", smp = SAMPLES)

rule bbduk:
  input:
      r1 = "raw/{smp}_1.fastq.gz",
      r2 = "raw/{smp}_2.fastq.gz"
  output:
      ro1 = "results/{smp}_1.clean.fastq.gz",
      ro2 = "results/{smp}_2.clean.fastq.gz"
  message:
      "Processing with BBDuk."
  shell:
      "{BBDUK} in1={input.r1} in2={input.r2} ref=adapters out1={output.ro1} out2={output.ro2} ktrim=r k=23 mink=11 hdist=1 tpe tbo"

rule star:
  input:
      r1 = "results/{smp}_1.clean.fastq.gz",
      r2 = "results/{smp}_2.clean.fastq.gz"
  output:
      "results/{smp}.sorted.bam"
  message:
      "Aligning with STAR."
  shell:
      "{STAR} --runThreadN 8 --genomeDir {STARINDEX} --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat --sjdbGTFfile {GTF} --outStd SAM | {SAMTOOLS} sort -o {output}"

rule bam_indexing:
  input:
      "results/{smp}.sorted.bam"
  output:
      "results/{smp}.sorted.bam.bai"
  message:
      "Indexing BAM with samtools."
  shell:
      "{SAMTOOLS} index {input}"

rule bam_coverage:
  input:
      r1 = "results/{smp}.sorted.bam",
      w2 = "results/{smp}.sorted.bam.bai"
  output:
      "results/{smp}.sorted.bw"
  message:
      "Calculating coverage with deeptools."
  shell:
      "{BAMCOVERAGE} -b {input.r1} -o {output} -p 8"

rule featurecounts:
  input:
      r1 = "results/{smp}.sorted.bam",
      w2 = "results/{smp}.sorted.bam.bai",
      w3 = "results/{smp}.sorted.bw"
  output:
      "results/{smp}.counts.txt"
  message:
      "Counting reads with featureCounts."
  shell:
      "{FEATURECOUNTS} -p --countReadPairs -t exon -g gene_id -a {GTF} -o {output} {input.r1}"

rule salmon:
  input:
      r1 = "results/{smp}_1.clean.fastq.gz",
      r2 = "results/{smp}_2.clean.fastq.gz",
      w3 = "results/{smp}.counts.txt"
  output:
      "quants/{smp}/quant.sf"
  params:
      dir = "quants/{smp}"
  message:
      "Quantifying with salmon."
  shell:
      "{SALMON} quant -i {SALMONINDEX} -l A -p 8 --gcBias --validateMappings -o {params.dir} -1 {input.r1} -2 {input.r2}"
```

Below is the shell script I will use to process these data with snakemake. I had to write a special shell script for this because the download tool is so unreliable (they recently updated it, so maybe it is better now).

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211129_ewsCohortEgaPmid34009296"
pyega3="/home/chughes/softwareTools/pyega3/ega-download-client-master"
eval mkdir raw
eval mkdir results
eval mkdir quants

##loop over the accessions
for j in T{1..50}
do
  eval cd ${workingDirectory}
  printf "Downloading files associated with ${j}."
  ############get the files associated with an accession
  egaFileId=($(awk 'BEGIN {FS="\t"; OFS="\t"} {print $3, $4}' EGAD00001004493/delimited_maps/Sample_File.map | grep "^${j}_" | awk 'BEGIN {FS="\t"; OFS=" "} {print $2}'))
  rnaFileId=($(awk 'BEGIN {FS="\t"; OFS="\t"} {print $3, $4}' EGAD00001004493/delimited_maps/Sample_File.map | grep "^${j}_" | awk 'BEGIN {FS="\t"; OFS=" "} {print $1}'))
  ############download the file
  eval cd ${pyega3}
  for i in $( seq 0 $(( ${#egaFileId[@]} - 1)) )
  do
    if [ -f "${rawDataOutputDirectory}/raw/${rnaFileId[$i]::-4}" ]; then
      printf "Raw file for ${rnaFileId[$i]::-4} already exists, skipping file.\n\n"
    else
      while [ ! -f "${rawDataOutputDirectory}/raw/${rnaFileId[$i]::-4}" ]
      do
        printf "File ${rnaFileId[$i]::-4} doesn't exist, attempting download."
        eval python -m pyega3.pyega3 -c 10 -cf ${rawDataOutputDirectory}credential_file.json fetch ${egaFileId[$i]} --saveto ${rawDataOutputDirectory}/raw/${rnaFileId[$i]::-4}
      done
    fi
  done
  
  #############processing the files associated with an accession
  printf "Combining read 1 files."
  eval cat raw/${j}_*_R1_*.gz > raw/${j}_1.fastq.gz
  printf "Combining read 2 files."
  eval cat raw/${j}_*_R2_*.gz > raw/${j}_2.fastq.gz
  printf "Running snakemake."
  #eval conda activate snakemake
  eval snakemake --cores 8
  #eval conda deactivate
  printf "Cleaning up."
  eval rm ${workingDirectory}/raw/${i}*.fastq.gz
  eval rm ${workingDirectory}/raw/${i}*.md5
  eval rm ${workingDirectory}/results/${i}*.clean.fastq.gz
  eval rm ${workingDirectory}/*.out
  eval rm ${workingDirectory}/*.tab
  eval rm -r ${workingDirectory}/*STAR*
done
```


for j in T{33..57}
do
  ############downloading code
  egaFileId=($(awk 'BEGIN {FS="\t"; OFS="\t"} {print $3, $4}' EGAD00001004493/delimited_maps/Sample_File.map | grep "^${j}_" | awk 'BEGIN {FS="\t"; OFS=" "} {print $2}'))
  rnaFileId=($(awk 'BEGIN {FS="\t"; OFS="\t"} {print $3, $4}' EGAD00001004493/delimited_maps/Sample_File.map | grep "^${j}_" | awk 'BEGIN {FS="\t"; OFS=" "} {print $1}'))
  #fileListLength=$(( ${#egaFileId[@]} - 1))
  for i in $( seq 0 $(( ${#egaFileId[@]} - 1)) )
  do
    if [ -f "${rawDataOutputDirectory}${rnaFileId[$i]::-4}" ]; then
      printf "Raw file for ${rnaFileId[$i]::-4} already exists, skipping file.\n\n"
    else
      while [ ! -f "${rawDataOutputDirectory}${rnaFileId[$i]::-4}" ]
      do
        printf "File ${rnaFileId[$i]::-4} doesn't exist, attempting download."
        eval pyega3 -c 10 -cf ${rawDataOutputDirectory}credential_file.json fetch ${egaFileId[$i]} --saveto ${rawDataOutputDirectory}${rnaFileId[$i]::-4}
      done
    fi
  done

#############processing code
  eval "cat ${j}_*_R1_*.gz > ${j}_R1_combined.fastq.gz"
  eval "cat ${j}_*_R2_*.gz > ${j}_R2_combined.fastq.gz"
  ##
  salmonCall="$salmonLocation quant -i $indexLocation -l A -1 ${rawDataOutputDirectory}${j}_R1_combined.fastq.gz -2 ${j}_R2_combined.fastq.gz --gcBias --validateMappings -o ${rawDataOutputDirectory}${j}_quant"
  eval $salmonCall

  ##
  starCall="$starLocation --runThreadN 12 --genomeDir ${referenceLocation} --readFilesIn ${rawDataOutputDirectory}${j}_R1_combined.fastq.gz ${rawDataOutputDirectory}${j}_R2_combined.fastq.gz --readFilesCommand zcat --sjdbGTFfile ${annotationLocation} --outFileNamePrefix ${rawDataOutputDirectory}${j}_ --outSAMtype BAM Unsorted"
    eval $starCall

  ##
  bamSortCall="$samtoolsLocation sort ${rawDataOutputDirectory}${j}_Aligned.out.bam -o ${rawDataOutputDirectory}${j}.sorted.bam"
  eval $bamSortCall

  ##
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}${j}.sorted.bam"
  eval $bamIndexCall

  ##
  coverageCall="bamCoverage -b ${rawDataOutputDirectory}${j}.sorted.bam -o ${rawDataOutputDirectory}${j}.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --centerReads -p 6"
  eval $coverageCall

  ##
  eval "rm ${j}*.fastq.gz"
  eval "rm ${j}*.md5"
  eval "rm *.out.bam"
  #done
done
