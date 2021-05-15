## Looking at the ENCODE eCLIP data

This document describes how to look at eCLIP data on a large scale from ENCODE.

### Getting the peak data

To download the data, I went to the [ENCODE browser](https://www.encodeproject.org/) and selected the eCLIP experiment type, GRCh38 for genome, and narrowPeak bed files as the file type. I then clicked the 'Download' button to create a 'files.txt' file that has the batch download links in it. I saved this file in `/projects/ptx_results/Sequencing/publishedStudies/encodeEclipData/bedFiles`. I then downloaded them to our server using the following command.

```shell
xargs -L 1 curl -O -J -L < files.txt
```

There is a metadata file that also gets downloaded with this and there is some info [here](https://www.encodeproject.org/help/batch-download/) on what it contains.

The rest of the data processing I will do in R.
