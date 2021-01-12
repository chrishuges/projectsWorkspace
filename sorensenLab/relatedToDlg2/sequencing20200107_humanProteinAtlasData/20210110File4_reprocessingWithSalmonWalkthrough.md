## Reprocessing some RNAseq data

This document details analysis of the human protein atlas RNAseq data with Salmon. There is a great walkthrough [here](http://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html) that I am following for the most part. The first thing I want to do is get an index for Salmon. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the Salmon index with the command `refgenie pull hg38/salmon_partial_sa_index`. You don't need to do this every time, just when you need to update your index. 

We can now run Salmon. We will do this using a script called `runSalmon.sh`.

```shell
#!/bin/bash
salmonLocation="/projects/ptx_analysis/chughes/software/salmon-1.4/bin/salmon"
indexLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/salmon_partial_sa_index/default/"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/humanProteinAtlasRnaSequencingData/"
for i in ERR315455 ERR315477 ERR315432
do
  echo $i
  salmonCall="$salmonLocation quant -i $indexLocation -l A -1 ${rawDataOutputDirectory}${i}_1.fastq.gz -2 ${rawDataOutputDirectory}${i}_2.fastq.gz --gcBias --validateMappings -o ${rawDataOutputDirectory}${i}_quant"
  eval $salmonCall
done
```

This will allow us to determine which isoform is the highest expressed from the data and we can use this information to for mapping coverage later on. 

