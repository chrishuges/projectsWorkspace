## Performing ROSE analysis

This document describes the process of assigning super enhancer elements from the ChIPseq data using the [ROSE](https://github.com/stjude/ROSE) tool (more info also [here](http://younglab.wi.mit.edu/super_enhancer_code.html)).

#### Setting up ROSE

I am going to try running this on the linux server using the Python 3 version of ROSE. I first downloaded the ZIP of the ROSE repository to `/projects/ptx_analysis/chughes/software/rosePython3`. I then had to run some commands to get it set up:

```bash
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

After looking through this a bit more, the python3 version asks for this 'annotation' file that should be in UCSC table format, and I cannot for the life of me figure out what this is. The python2 version does not require this though, so I am going to try that version. In order to get that going, I will use my virtual python environment. To load it up, use the command `export PATH=~/Virtual_Python2712/bin:$PATH`. There is a version of ROSE on pypi, but it seems to be an old one. The one hosted on their [bitbucket](https://bitbucket.org/young_computation/rose/src/master/) seems to be updated to allow for hg38, which is what we need. I downloaded the repository to `/projects/ptx_analysis/chughes/software/rosePython2`.

So, something I just discovered when I was perusing the python2 rose files is that they provide the 'annotation' folder that contains the UCSC files required by the python3 version. So, I am going to go back to using the python3 version and just use these provided UCSC files.

#### Running ROSE

Run in directory: `/projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq`.

Run on files with the `.sorted.bam` extension. 

```bash
ROSE_main.py -g HG38 -i /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/SRR8832667_peaks.broadPeak -r /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/SRR8832674.filtered.bam -c /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/SRR8832666.filtered.bam -o /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/rose
```

Running this is gives an error that I think is related to the format of the 'GFF' file. The direct MACS output is not completely correct, so I will need to reshape this. The current format of the broad peak file is:

```
chr1    904187  906674  SRR8832667_peak_1       67      .       5.2712  8.66262 6.71981
```

But, it needs to be:

```
1: chromosome (chr#)
2: unique ID for each constituent enhancer region
4: start of constituent
5: end of constituent
7: strand (+,-,.)
9: unique ID for each constituent enhancer region

#example
chr6	MM1S_MED1_DMSO_2_15567		41854892	41856607		.		MM1S_MED1_DMSO_2_15567
```


I think we can do this easily with awk.

```bash
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $4, $2, $3, $6, $4}' SRR8832667_peaks.broadPeak > SRR8832667_peaks.gff
```

Try calling ROSE again.

```bash
ROSE_main.py -g HG38 -i /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/SRR8832667_peaks.gff -r /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/SRR8832674.filtered.bam -c /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/SRR8832666.filtered.bam -o /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/rose
```

It still gives me an error. I think it is still a problem with the GFF file. To make sure of this, I am going to run it with the ROSE provided example files.

```bash
export PATH=~/Virtual_Python2712/bin:$PATH
python ROSE_main.py -g HG18 -i ./ROSE_DATA/data/HG18_MM1S_MED1_1000.gff -r ./ROSE_DATA/data/MM1S_MED1.hg18.bwt.sorted.bam -c ./ROSE_DATA/data/MM1S_WCE.hg18.bwt.sorted.bam -o example -s 12500 -t 2500
```

This works with the python2 version. It took about 5 minutes to complete for that example file that only has 1000 entries. I guess you can subset what you are looking for by using this smaller GFF, which could be useful as we could just provide a subset of DLG2-related regions. It does give some beef with the R code, specifically related to the export of the png image, I will have to fiddle around with this I think.

To fix this error, I edited the R script line `png(filename=plotFileName,height=600,width=600)` to instead be `png(filename=plotFileName,height=600,width=600,type='cairo')`. Re-running this worked fine and all of the output files were created correctly.

Now I want to try with the python3 version just because it seems a bit more flexible in terms of overall use. I edited the R script in the same way in order to avoid an error with creating the png. 

```bash
export PATH=~/Virtual_Python361/bin:$PATH
ROSE_main.py -g HG18 -i ./ROSE_DATA/data/HG18_MM1S_MED1_1000.gff -r ./ROSE_DATA/data/MM1S_MED1.hg18.bwt.sorted.bam -c ./ROSE_DATA/data/MM1S_WCE.hg18.bwt.sorted.bam -o example -s 12500 -t 2500
```

Now I think I am ready to run this on my own data. But I guess I need to figure out this GFF problem and I am really unsure of why it has issues with the format. It seems like it might be because it expects empty columns?

```bash
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $4, '\t', $2, $3, '\t', $6, '\t', $4}' SRR8832667_peaks.broadPeak > SRR8832667_peaks.gff
grep chr11 SRR8832667_peaks.gff > SRR8832667_chr11Peaks.gff
ROSE_main.py -g HG38 -i /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/SRR8832667_chr11Peaks.gff -r /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/SRR8832667.filtered.bam -c /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/SRR8832666.filtered.bam -o /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/rose
```

Alright that GFF worked. Unfortunately, I ran it in a normal session so I might get booted out. I will have to go back and run it in a screen because the program is quite slow. I will also repeat it with the day7 sample from the histone ChIP data. I ended up doing this on a much smaller gff input that just has chromosome 11 in it (about 4500 entries). 

















