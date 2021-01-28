## Performing ROSE analysis

This document describes the process of assigning super enhancer elements from the ChIPseq data using the [ROSE](https://github.com/stjude/ROSE) tool (more info also [here](http://younglab.wi.mit.edu/super_enhancer_code.html)).

### Setting up ROSE

I am going to try running this on the linux server using the Python 3 version of ROSE. I first downloaded the ZIP of the ROSE repository to `/projects/ptx_analysis/chughes/software/rosePython3`. I then had to run some commands to get it set up:

```shell
PATHTO=/projects/ptx_analysis/chughes/software/rosePython3/ROSE-master/
PYTHONPATH=$PATHTO/lib
export PYTHONPATH
export PATH=$PATH:$PATHTO/bin
ROSE_main.py #this command should work now

#output
ROSE_main.py [options] -g [GENOME] -i [INPUT_REGION_GFF] -r [RANKBY_BAM_FILE] -o [OUTPUT_FOLDER] [OPTIONAL_FLAGS]
```

The next thing I needed to do was prepare my input files. As far as I can tell, ROSE requires:

    * All files : All input files much be in one directory.
    * Annotation file : Annotation file should be in UCSC table track format (https://genome.ucsc.edu/cgi-bin/hgTables). Annotation file should be saved as [GENOME]_refseq.ucsc (example: hg19_refseq.ucsc). Annotation file should be in annotation/ folder in the input files directory.
    * BAM files (of sequencing reads for factor of interest and control) : Files must have chromosome IDs starting with "chr" Files must be sorted and indexed using SAMtools in order for bamToGFF.py to work. (http://samtools.sourceforge.net/samtools.shtml)
    * Peak file of constituent enhancers : File must be in GFF format with the following columns: 1: chromosome (chr#) 2: unique ID for each constituent enhancer region 4: start of constituent 5: end of constituent 7: strand (+,-,.) 9: unique ID for each constituent enhancer region NOTE: if value for column 2 and 9 differ, value in column 2 will be used


I already have the bam files and they are already in a single directory, so we are good there. 