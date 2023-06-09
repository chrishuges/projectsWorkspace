## Scanning for motifs in ChIPseq data

This document describes the analysis of our called ChIPseq peaks in the DLG2 gene to determine motif enrichment (specifically, GGAA). 

### Data processing

For this I will use the [FIMO](http://meme-suite.org/tools/fimo) tool, as this is what was used in the original manuscript describing these ChIPseq data. It looks like GSC already has a version of the MEME suite installed at `/gsc/software/linux-x86_64-centos7/meme-4.12.0/bin`. 

```shell
Usage: fimo [options] <motif file> <sequence file>

   Options:
     --alpha <double> (default 1.0)
     --bgfile <background file> (DNA and protein use NRDB by default)
     --max-stored-scores <int> (default=100000)
     --max-strand
     --motif <id> (default=all)
     --motif-pseudo <float> (default=0.1)
     --no-qvalue
     --norc
     --o <output dir> (default=fimo_out)
     --oc <output dir> (default=fimo_out)
     --parse-genomic-coord
     --psp <PSP filename> (default none)
     --prior-dist <PSP distribution filename> (default none)
     --qv-thresh
     --skip-matched-sequence
     --text
     --thresh <float> (default = 1e-4)
     --verbosity [1|2|3|4] (default 2)
     --version (print the version and exit)

   When scanning with a single motif use '-' for <sequence file> to
     read the database from standard input.
   Use '--bgfile --motif--' to read the background from the motif file.
   Use '--bgfile --uniform--' to use a uniform background.
```

The first thing we need is a MEME format file of the motif's we want to search for. We can get these from [JASPAR](http://jaspar.genereg.net/). 

    * [ETS motif](http://jaspar.genereg.net/matrix/MA0475.1/)
    * [GGAA mSats](http://jaspar.genereg.net/matrix/MA0149.1/)

The next thing we need is the sequences of our ChIPseq peaks. We can do this in R. I did it in the file `20210118File4_dlg2ChipSeqPeakSequences.Rmd`. After this I am ready to run the tool.

```shell
/gsc/software/linux-x86_64-centos7/meme-4.12.0/bin/fimo -o fimo_MA0475.1 /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/fimo/MA0475.1.meme /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/fimo/SRR8832674_peakSequences.fasta

/gsc/software/linux-x86_64-centos7/meme-4.12.0/bin/fimo -o fimo_MA0149.1 /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/fimo/MA0149.1.meme /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/fimo/SRR8832674_peakSequences.fasta
```