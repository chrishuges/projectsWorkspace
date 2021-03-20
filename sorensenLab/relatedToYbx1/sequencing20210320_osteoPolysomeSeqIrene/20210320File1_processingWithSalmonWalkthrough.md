## Processing some RNAseq data

This document details analysis of the polysomeSeq analysis Irene did on osteosarcoma cells data with Salmon. There is a great walkthrough [here](http://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html) that I am following for the most part. The first thing I want to do is get an index for Salmon. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the Salmon index with the command `refgenie pull hg38/salmon_partial_sa_index`. You don't need to do this every time, just when you need to update your index. 

The raw data are stored here: `/projects/ptx_results/Sequencing/2019/sequencing20190323_osteoPolysomeSeqIrene`. Rename the files.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/2019/sequencing20190323_osteoPolysomeSeqIrene/"

#previous example filename CD26TANXX_2_1_AAACCT_75bp_229093.concat_chastity_passed.fastq.gz
######################################
for i in AAACCT AGATGT ATGGCA CAGCAG CGGCCT CTTTTG GCCATG GCGGAC GCTGTA TACAAG TATCGT TCCGTC
do
  echo $i
  eval "mv ${rawDataInputDirectory}CD26TANXX_2_1_${i}*.fastq.gz ${rawDataOutputDirectory}${i}_1.fastq.gz"
  eval "mv ${rawDataInputDirectory}CD26TANXX_2_2_${i}*.fastq.gz ${rawDataOutputDirectory}${i}_2.fastq.gz"
done
```

We can now run Salmon. We will do this using a script called `runSalmon.sh`.

```shell
#!/bin/bash
bbdukLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbduk.sh"
salmonLocation="/projects/ptx_analysis/chughes/software/salmon-1.4/bin/salmon"
indexLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/alias/hg38/salmon_partial_sa_index/default/"
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/alias/hg38/gencode_gtf/default"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/2019/sequencing20190323_osteoPolysomeSeqIrene/"

###########################################
for i in AAACCT AGATGT ATGGCA CAGCAG CGGCCT CTTTTG GCCATG GCGGAC GCTGTA TACAAG TATCGT TCCGTC
do
  echo $i
  bbdukCall="$bbdukLocation in1=${rawDataOutputDirectory}${i}_1.fastq.gz in2=${rawDataOutputDirectory}${i}_2.fastq.gz ref=adapters out1=${rawDataOutputDirectory}${i}_1.clean.fastq.gz out2=${rawDataOutputDirectory}${i}_2.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo"
  eval $bbdukCall
  ##
  salmonCall="$salmonLocation quant -i $indexLocation -l A -1 ${rawDataOutputDirectory}${i}_1.clean.fastq.gz -2 ${rawDataOutputDirectory}${i}_2.clean.fastq.gz --gcBias --validateMappings -o ${rawDataOutputDirectory}${i}_quant"
  eval $salmonCall
done
```

For downstream analysis, file file layout is:

PROJECT SOW     LIBRARY INDEX   FLOWCELL        LANE    ANTIBODY        UPPER_PROTOCOL  LOWER_PROTOCOL  EXPERIMENTAL_CONDITION  CELL_CONDITION  TAXONOMY_ID     SPECIES SEQ LENGTH      PATIENT_ID      EXTERNAL_IDENTIFIER     ORIGINAL_SOURCE SUBMITTED_SPIKE_IN_SEQUENCE     SUBMITTED_SPIKE_IN_IDENTIFIER   ALERT_CODE      ALERT_NOTES
Poul Sorensen   GSC-1719        B06255  CTTTTG  CD26TANXX       2       None    Ribodepletion 2.1       Illumina Indexing                       9606    Homo sapiens    75      MNNG_GFP        Control_1       PoS_51  None    None
Poul Sorensen   GSC-1719        B06256  GCGGAC  CD26TANXX       2       None    Ribodepletion 2.1       Illumina Indexing                       9606    Homo sapiens    75      MNNG_GFP        Control_2       PoS_52  None    None
Poul Sorensen   GSC-1719        B06257  AAACCT  CD26TANXX       2       None    Ribodepletion 2.1       Illumina Indexing                       9606    Homo sapiens    75      MNNG_GFP        Control_3       PoS_53  None    None
Poul Sorensen   GSC-1719        B06258  GCTGTA  CD26TANXX       2       None    Ribodepletion 2.1       Illumina Indexing                       9606    Homo sapiens    75      MNNG_GFP        MS275_1 PoS_54  None    None
Poul Sorensen   GSC-1719        B06259  CGGCCT  CD26TANXX       2       None    Ribodepletion 2.1       Illumina Indexing                       9606    Homo sapiens    75      MNNG_GFP        MS275_2 PoS_55  None    None
Poul Sorensen   GSC-1719        B06260  TACAAG  CD26TANXX       2       None    Ribodepletion 2.1       Illumina Indexing                       9606    Homo sapiens    75      MNNG_GFP        MS275_3 PoS_56  None    None
Poul Sorensen   GSC-1719        B06261  TATCGT  CD26TANXX       2       None    Ribodepletion 2.1       Illumina Indexing                       9606    Homo sapiens    75      MNNG_GFP        Arsenite_1      PoS_57  None    None
Poul Sorensen   GSC-1719        B06262  GCCATG  CD26TANXX       2       None    Ribodepletion 2.1       Illumina Indexing                       9606    Homo sapiens    75      MNNG_GFP        Arsenite_2      PoS_58  None    None
Poul Sorensen   GSC-1719        B06263  ATGGCA  CD26TANXX       2       None    Ribodepletion 2.1       Illumina Indexing                       9606    Homo sapiens    75      MNNG_GFP        Arsenite_3      PoS_59  None    None
Poul Sorensen   GSC-1719        B06264  TCCGTC  CD26TANXX       2       None    Ribodepletion 2.1       Illumina Indexing                       9606    Homo sapiens    75      MNNG_GFP        MS275_Ars_1     PoS_60  None    None
Poul Sorensen   GSC-1719        B06265  AGATGT  CD26TANXX       2       None    Ribodepletion 2.1       Illumina Indexing                       9606    Homo sapiens    75      MNNG_GFP        MS275_Ars_2     PoS_61  None    None
Poul Sorensen   GSC-1719        B06266  CAGCAG  CD26TANXX       2       None    Ribodepletion 2.1       Illumina Indexing                       9606    Homo sapiens    75      MNNG_GFP        MS275_Ars_3     PoS_62  None    None

Now that we have these data, process them downstream with R.