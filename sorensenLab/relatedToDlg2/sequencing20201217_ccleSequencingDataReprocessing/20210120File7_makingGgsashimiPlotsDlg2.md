## Reprocessing some RNAseq data

This document describes the creation of some [ggSashimi](https://github.com/guigolab/ggsashimi) plots for Dlg2.

## Workflow

The first thing I did was to follow the installation instructions. This was basically downloading the tool and activating it. I saved the script in the directory `/projects/ptx_analysis/chughes/software/ggSashimi`.

```shell
wget https://raw.githubusercontent.com/guigolab/ggsashimi/master/ggsashimi.py
chmod u+x ggsashimi.py
./ggsashimi.py --help
```

The help page is pretty useful.

```
usage: ggsashimi.py [-h] -b BAM -c COORDINATES [-o OUT_PREFIX] [-S OUT_STRAND]
                    [-M MIN_COVERAGE] [-j JUNCTIONS_BED] [-g GTF] [-s STRAND]
                    [--shrink] [-O OVERLAY] [-A AGGR] [-C COLOR_FACTOR]
                    [--alpha ALPHA] [-P PALETTE] [-L LABELS] [--fix-y-scale]
                    [--height HEIGHT] [--ann-height ANN_HEIGHT]
                    [--width WIDTH] [--base-size BASE_SIZE] [-F OUT_FORMAT]
                    [-R OUT_RESOLUTION] [--debug-info] [--version]

Create sashimi plot for a given genomic region

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     Individual bam file or file with a list of bam files.
                        In the case of a list of files the format is tsv:
                        1col: id for bam file, 2col: path of bam file, 3+col:
                        additional columns
  -c COORDINATES, --coordinates COORDINATES
                        Genomic region. Format: chr:start-end. Remember that
                        bam coordinates are 0-based
  -o OUT_PREFIX, --out-prefix OUT_PREFIX
                        Prefix for plot file name [default=sashimi]
  -S OUT_STRAND, --out-strand OUT_STRAND
                        Only for --strand other than 'NONE'. Choose which
                        signal strand to plot: <both> <plus> <minus>
                        [default=both]
  -M MIN_COVERAGE, --min-coverage MIN_COVERAGE
                        Minimum number of reads supporting a junction to be
                        drawn [default=1]
  -j JUNCTIONS_BED, --junctions-bed JUNCTIONS_BED
                        Junction BED file name [default=no junction file]
  -g GTF, --gtf GTF     Gtf file with annotation (only exons is enough)
  -s STRAND, --strand STRAND
                        Strand specificity: <NONE> <SENSE> <ANTISENSE>
                        <MATE1_SENSE> <MATE2_SENSE> [default=NONE]
  --shrink              Shrink the junctions by a factor for nicer display
                        [default=False]
  -O OVERLAY, --overlay OVERLAY
                        Index of column with overlay levels (1-based)
  -A AGGR, --aggr AGGR  Aggregate function for overlay: <mean> <median>
                        <mean_j> <median_j>. Use mean_j | median_j to keep
                        density overlay but aggregate junction counts
                        [default=no aggregation]
  -C COLOR_FACTOR, --color-factor COLOR_FACTOR
                        Index of column with color levels (1-based)
  --alpha ALPHA         Transparency level for density histogram [default=0.5]
  -P PALETTE, --palette PALETTE
                        Color palette file. tsv file with >=1 columns, where
                        the color is the first column. Both R color names and
                        hexadecimal values are valid
  -L LABELS, --labels LABELS
                        Index of column with labels (1-based) [default=1]
  --fix-y-scale         Fix y-scale across individual signal plots
                        [default=False]
  --height HEIGHT       Height of the individual signal plot in inches
                        [default=2]
  --ann-height ANN_HEIGHT
                        Height of annotation plot in inches [default=1.5]
  --width WIDTH         Width of the plot in inches [default=10]
  --base-size BASE_SIZE
                        Base font size of the plot in pch [default=14]
  -F OUT_FORMAT, --out-format OUT_FORMAT
                        Output file format: <pdf> <svg> <png> <jpeg> <tiff>
                        [default=pdf]
  -R OUT_RESOLUTION, --out-resolution OUT_RESOLUTION
                        Output file resolution in PPI (pixels per inch).
                        Applies only to raster output formats [default=300]
  --debug-info          Show several system information useful for debugging
                        purposes [default=None]
  --version             show program's version number and exit
```

We can first try this on one of our A673 bam files where EWS-FLI1 is expressed. I am interested in a specific exon in the isoform ENST00000426717.6 (short) because I want to show that the 5'UTR of this transcript actually extends beyond the annotation, meaning it is likely not expressed. For this, I will use the region containing the first exon and where the junctions that extend from it end (determined from IGV). I will do this in my A673 data file from CCLE (SRR8616012). The region is `chr11:83,630,396-83,762,454`. 

```shell
/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8616012.sorted.bam -c chr11:83626956-83760000 -g /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf --shrink
```

This didn't work and gave me an error about 'transcript_id'. I seem to remember having this issue before and it had to do with how the annotation was set up in the gtf with not everything having a transcript_id or something. We can try a different GTF, also from refgenie, to see if it helps.

```shell
/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8616012.sorted.bam -c chr11:83626956-83760000 -g /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/ensembl_gtf/default/hg38.gtf --min-coverage 50
```

For some reason it isn't showing the annotation. Maybe it doesn't like our GTF file still? We can try feeding it one with just DLG2 exons.

```shell
grep DLG2 /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf

/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8616012.sorted.bam -c chr11:83626956-83760000 -g /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf --min-coverage 50
```

It gives the error again about the transcript_id. What if we try with the other GTF.

```shell
grep DLG2 /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/ensembl_gtf/default/hg38.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf

/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8616012.sorted.bam -c chr11:83626956-83760000 -g /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf --min-coverage 50 --base-size 9
```

It still doesn't show the gene annotation. What if we feed it only exons. 

```shell
grep exon /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf
awk '!/CDS/ && !/start_codon/' /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly2.gtf

/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8616012.sorted.bam -c chr11:83626956-83760000 -g /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly2.gtf --min-coverage 50 --base-size 9
```

Still not working. Maybe try with a different gene.

```shell
/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8616012.sorted.bam -c chr11:85624866-85638811 -g /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/ensembl_gtf/default/hg38.gtf --min-coverage 50 --base-size 9
```

Ok that also didn't work, making me think it is something to do with our GTF. Lets try again with the gencode GTF that we used during alignment. I have separate code here for different server machines.

```shell
grep DLG2 /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2Only.gtf

grep exon /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2Only.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf

awk '!/CDS/ && !/start_codon/' /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnlyFiltered.gtf

/home/chughes/softwareTools/ggSashimi/ggsashimi.py -b /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/results/SRR7059715.sorted.bam -c chr11:83626956-83760000 -g /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnlyFiltered.gtf --min-coverage 10 --base-size 7 -o chla10Dlg2Upstream

#grep DLG2 /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf

#grep exon /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf

#awk '!/CDS/ && !/start_codon/' /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly2.gtf

#/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8616012.sorted.bam -c chr11:83626956-83760000 -g /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly2.gtf --min-coverage 50 --base-size 9
```

Ok that worked, so it was our GTF. We can filter this a bit more because it is showing many isoforms we don't want. First try this a different way.

```shell
grep DLG2 /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2Only.gtf

awk '!/CDS/ && !/start_codon/ && !/UTR/ && !/stop_codon/ && /exon/' /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2Only.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf

/home/chughes/softwareTools/ggSashimi/ggsashimi.py -b /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/results/SRR7059715.sorted.bam -c chr11:83626956-83760000 -g /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf --min-coverage 10 --base-size 9 -o chla10Dlg2Upstream

#grep DLG2 /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf

#awk '!/CDS/ && !/start_codon/ && !/UTR/ && !/stop_codon/ && /exon/' /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf

#/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8616012.sorted.bam -c chr11:83626956-83760000 -g /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf --min-coverage 50 --base-size 9
```

Ok that works and saves us a line of code. We can try for a specific transcript version.

```shell
grep 'ENST00000426717.6\|ENST00000398309.6' /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2Only.gtf

awk '!/start_codon/ && !/stop_codon/ && !/UTR/ && !/^CDS$/' /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2Only.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf

/home/chughes/softwareTools/ggSashimi/ggsashimi.py -b /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/results/SRR7059715.sorted.bam -c chr11:83626956-83760000 -g /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf --min-coverage 10 --base-size 9 -o chla10Dlg2Upstream

#grep 'ENST00000426717.6\|ENST00000398309.6' /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf

#awk '!/start_codon/ && !/stop_codon/ && !/UTR/ && !/^CDS$/' /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf

#/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8616012.sorted.bam -c chr11:83626956-83760000 -g /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf --min-coverage 50 --base-size 9
```

This works great! The last thing I want to do is pass it multiple BAM files.

```shell
/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/bamFileAccessions.txt -c chr11:83626956-83760000 -g /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf --min-coverage 50 --base-size 9
```

These plots don't seem to play well with Illustrator, so make just one at a time.

```shell
/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8616012.sorted.bam -o a673 -c chr11:83626956-83760000 -g /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf --min-coverage 50 --base-size 9

/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8615497.sorted.bam -o sknmc -c chr11:83626956-83760000 -g /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf --min-coverage 50 --base-size 9

/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/additionalEwingLines/starResults/SRR8615521.sorted.bam -o mhhes1 -c chr11:83626956-83760000 -g /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf --min-coverage 50 --base-size 9
```
