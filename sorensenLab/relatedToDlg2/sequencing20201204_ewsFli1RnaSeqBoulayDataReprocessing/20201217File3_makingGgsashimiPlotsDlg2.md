## Reprocessing some RNAseq data

This document describes the creation of some [ggSashimi](https://github.com/guigolab/ggsashimi) plots for Dlg2.

## Workflow

The first thing I did was to follow the installation instructions. This was basically downloading the tool and activating it. I saved the script in the directory `/projects/ptx_analysis/chughes/software/ggSashimi`.

```shell
wget https://raw.githubusercontent.com/guigolab/ggsashimi/master/sashimi-plot.py
chmod u+x sashimi-plot.py
./sashimi-plot.py --help
```

The help page is pretty useful.

```
usage: sashimi-plot.py [-h] -b BAM -c COORDINATES [-o OUT_PREFIX]
                       [-S OUT_STRAND] [-M MIN_COVERAGE] [-j JUNCTIONS_BED]
                       [-g GTF] [-s STRAND] [--shrink] [-O OVERLAY] [-A AGGR]
                       [-C COLOR_FACTOR] [--alpha ALPHA] [-P PALETTE]
                       [-L LABELS] [--fix-y-scale] [--height HEIGHT]
                       [--ann-height ANN_HEIGHT] [--width WIDTH]
                       [--base-size BASE_SIZE] [-F OUT_FORMAT]
                       [-R OUT_RESOLUTION]

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

```

We can first try this on one of our A673 bam files where EWS-FLI1 is expressed. I took the location of [ENST00000280241.12](http://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000150672;r=11:83455012-85627922;t=ENST00000280241) as my range to plot.

```shell
/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/201709BoulayCellPmid28844694/SRR5217668.chr11.bam -c 11:83455012-84317339 -g /projects/ptx_analysis/chughes/databases/HomoSapiensEnsemblGRCh38_rel102/Homo_sapiens.GRCh38.102.gtf --shrink
```

This didn't actually work for me. It gave an error:

```
Traceback (most recent call last):
  File "/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py", line 632, in <module>
    transcripts, exons = read_gtf(args.gtf, args.coordinates)
  File "/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py", line 289, in read_gtf
    transcript_id = d["transcript_id"]
KeyError: 'transcript_id'
```

I found some discussion online that this was perhaps due to the way the Ensembl annotation file is set up, so I tried one from [GENCODE](https://www.gencodegenes.org/human/).

```
/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/201709BoulayCellPmid28844694/SRR5217668.chr11.bam -c 11:83455012-84317339 -g /projects/ptx_analysis/chughes/databases/HomoSapiensEnsemblGRCh38_rel102/gencode.v36.annotation.gtf --shrink
```

Or, on my desktop running locally.

```
python3 /mnt/c/Users/chughes/Documents/bccrc/softwareRepository/ggSashimi/sashimi-plot.py -b /mnt/c/Users/chughes/Documents/bccrc/projectsRepository/sorensenLab/relatedToDlg2/sequencing20201204_ewsFli1RnaSeqBoulayDataReprocessing/SRR5217668.chr11.bam -c 11:83455012-84317339 -g /mnt/c/Users/chughes/Documents/bccrc/databases/HomoSapiensEnsemblGRCh38_rel102/gencode.v36.annotation.gtf --shrink
```

Note, none of the above actually worked. But, I went back and re-did the alignment with STAR and ran the command below and it works.

```shell
/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b ./SRR5217668_Aligned.sortedByCoord.out.bam -c chr11:83455012-84317339 -g /projects/ptx_analysis/chughes/databases/HomoSapiensEnsemblGRCh38_rel102/Homo_sapiens.GRCh38.102.gtf --shrink
```