## Snakemake test run

This document describes some snakemake testing to see about using it for the DLG2 project instead of the bash scripts I normally do.

## ChIP-Seq data

I am going to test with a ChIP-Seq pipeline from the Grunewald data as in the directory below. I have already installed Snakemake via Mamba as described on the Snakemake [homepage](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

First, I will download the data files. I do this with a simple bash script. First make a couple of directories and activate snakemake.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia
mkdir raw
mkdir results
conda activate snakemake
conda install -c conda-forge tree
```

Now I am going to run some basic analyses following a tutorial [here](https://eriqande.github.io/eca-bioinf-handbook/managing-workflows-with-snakemake.html), as well as on the Snakemake [homepage](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) itself, and [this page](https://eriqande.github.io/eca-bioinf-handbook/managing-workflows-with-snakemake.html) is also useful, and [this page](http://pedagogix-tagc.univ-mrs.fr/courses/ABD/practical/snakemake/snake_intro.html). Enter the working directory and create a snakefile.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia
touch dataProcessingScript.sh
chmod +x dataProcessingScript.sh
touch snakefile
```

Edit the contents of the snakefile to process your data.

```shell
"""
Author: Christopher Hughes
Affiliation: BCCRC
Aim: Workflow for ChIP-Seq data
Date: 20211030
"""

###############################
#working directory
BASE_DIR = "/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia"


###############################
#locations of tools we will use
BBDUK = "/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
BOWTIE2 = "/home/chughes/softwareTools/bowtie2-2.4.4/bowtie2"
SAMTOOLS="/home/chughes/softwareTools/samtools-1.12/samtools"
SAMBAMBA="/home/chughes/softwareTools/sambamba-0.8.1/sambamba"

###############################
#locations of our index files
DATABASE_DIR = "/home/chughes/databases/projectEwsDlg2"
INDEX = DATABASE_DIR + "/bowtie2Index/human_genome"
GTF   = DATABASE_DIR + "/baseGenomeFiles/genome.gtf"
FASTA = DATABASE_DIR + "/baseGenomeFiles/genome.fa"


###############################
#get our list of sample files
SAMPLES, = glob_wildcards("raw/{smp}.fastq.gz")
print("There are " + str(len(SAMPLES)) + " total samples to be processed.")
for smp in SAMPLES:
  print("Sample " + smp + " will be processed")


###############################
#processing workflow
rule all:
    input: expand("results/{smp}.filtered.bam.bai", smp = SAMPLES)

rule bbduk:
  input:
      "raw/{smp}.fastq.gz"
  output:
      "results/{smp}.clean.fastq.gz"
  message:
      "Processing with BBDuk."
  shell:
      "{BBDUK} in={input} ref=adapters out={output} ktrim=r k=23 mink=11 hdist=1"

rule bowtie:
  input:
      "results/{smp}.clean.fastq.gz"
  output:
      "results/{smp}.unsorted.sam"
  message:
      "Aligning with bowtie2."
  shell:
      "{BOWTIE2} -p 8 -q --local -x {INDEX} -U {input} -S {output}"

rule bam_conversion:
  input:
      "results/{smp}.unsorted.sam"
  output:
      "results/{smp}.unsorted.bam"
  message:
      "Sorting with samtools."
  shell:
      "{SAMTOOLS} view -h -S -b -o {output} {input}"

rule bam_sorting:
  input:
      "results/{smp}.unsorted.bam"
  output:
      "results/{smp}.sorted.bam"
  message:
      "BAM sorting with sambamba."
  shell:
      "{SAMBAMBA} sort -t 6 -o {output} {input}"

rule bam_filtering:
  input:
      "results/{smp}.sorted.bam"
  output:
      "results/{smp}.filtered.bam"
  message:
      "Filtering BAM with sambamba."
  shell:
      """{SAMBAMBA} view -h -t 6 -f bam -F "[XS] == null and not unmapped and not duplicate" {input} > {output}"""

rule bam_indexing:
  input:
      "results/{smp}.filtered.bam"
  output:
      "results/{smp}.filtered.bam.bai"
  message:
      "Indexing BAM with samtools."
  shell:
      "{SAMTOOLS} index {input}"

```

Now we will use a shell script to interface the data download to the processing. I want to do it one at a time to save space, which is why I do it this way. Note, this code below is designed to work with ENA downloads, which for some reason stopped working around the time I was doing this analysis. I wrote another script that uses SRA directly as a temporary fix and it can be found below this.

```shell
#!/bin/bash

workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia/raw"
eval cd ${workingDirectory}
firstAccession=({{097..099},{000..021},{027..087}})
secondAccession=(SRR1476{{0997..1021},{1027..1087}})
for i in {0..85}
do
  downloadFile="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/${firstAccession[$i]}/${secondAccession[$i]}/${secondAccession[$i]}.fastq.gz"
  printf "Download path: $downloadFile\n"
  eval wget $downloadFile
  eval snakemake --cores 8
  eval rm ${secondAccession[$i]}*.fastq.gz
  eval rm ../results/${secondAccession[$i]}.clean.fastq.gz
  eval rm ../results/${secondAccession[$i]}.unsorted.sam
  eval rm ../results/${secondAccession[$i]}.unsorted.bam
done
```

As mentioned above, ENA stopped working for a bit so I needed to do this directly from SRA. I did this using [sraDownloader](https://github.com/s-andrews/sradownloader). This is not ideal, but at least it works.

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
sraDownloader="/home/chughes/softwareTools/sradownloader-3.8/sradownloader"
sraCacheLocation="/mnt/Data/chughes/sratoolsRepository"
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia"
eval cd ${workingDirectory}

##loop over the accessions
for i in SRR1476{{0997..1021},{1027..1087}}
do
  printf "Downloading files associated with ${i}."
  eval ${sraDownloader} --noena --outdir ${workingDirectory}/raw ${i}
  ##the file gets renamed upon download, but I just want it to have the SRR id and I can annotate it later
  eval mv ${workingDirectory}/raw/${i}*.fastq.gz ${workingDirectory}/raw/${i}.fastq.gz
  #eval conda activate snakemake
  eval snakemake --cores 8
  #eval conda deactivate
  eval rm ${workingDirectory}/raw/${i}.fastq.gz
  eval rm ${sraCacheLocation}/sra/${i}.sra.cache
  eval rm ${workingDirectory}/results/${i}.clean.fastq.gz
  eval rm ${workingDirectory}/results/${i}.unsorted.sam
  eval rm ${workingDirectory}/results/${i}.unsorted.bam
done
```


Now we can run the analysis, writing a log so we can track it. I am inside the working directory here, so I can just execute it. If you are not, you need to set the absolute path to the script. I also make the shell script executable in the code below, as well as activating the conda environment. I am a complete snakemake n00b, so I haven't figured out the nuts and bolts of making it work with bash script integration.

```shell
chmod +x sraDataProcessingScript.sh
conda activate snakemake
./sraDataProcessingScript.sh |& tee sraDataProcessingScriptLog.txt
```



