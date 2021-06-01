## Looking at the ENCODE eCLIP data

This document details attempting to run the Yeo lab pipeline for eCLIP data processing. 

### Running the pipelien

I am following the pipeline from the [ENCODE browser](https://www.encodeproject.org/) and selected the eCLIP experiment type, fastq files as the file type. I then clicked the 'Download' button to create a 'files.txt' file that has the batch download links in it. I saved this file in `/projects/ptx_results/Sequencing/publishedStudies/encodeEclipData/fastqFiles`. I then downloaded them to our server using the following command.

```shell
xargs -L 1 curl -O -J -L < files.txt
```

There is a metadata file that also gets downloaded with this and there is some info [here](https://www.encodeproject.org/help/batch-download/) on what it contains. I processed this file in R to make a new manifest file that contains download links just for two paired files (left and right reads) that I will run the pipeline on individually. I am following the ENCODE eCLIP SOP for this.

Program locations:

FastQC: /gsc/software/linux-x86_64-centos7/fastqc-0.11.9/fastqc
Cutadapt: 
STAR: /projects/ptx_analysis/chughes/software/STAR-2.7.6a/STAR-2.7.6a/bin/Linux_x86_64/STAR
Samtools: /projects/ptx_analysis/chughes/software/samtools-1.9/samtools
bedToBigBed: /projects/ptx_analysis/chughes/software/bedtoBigBed/
Bedtools: /gsc/software/linux-x86_64-centos7/bedtools-2.27.1/bin/bedtools
R:
fastq-sort:
umi_tools: 
python:
perl:


I am going to attempt to do this in Anaconda. I am using the conda version at `/gsc/software/linux-x86_64-centos7/Anaconda2-4.2.0/anaconda/bin`. I created a virtual environment to do my work in. I followed the instructions [here](https://www.geeksforgeeks.org/set-up-virtual-environment-for-python-using-anaconda/).

```shell
conda search "^python$"
conda create -n eclip python=2.7.12 anaconda
conda activate eclip
```


conda install -n yourenvname package
conda deactivate
```

