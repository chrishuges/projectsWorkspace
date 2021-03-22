## Processing some RNAseq data

This document details analysis of the polysomeSeq analysis Renata did on Ewing sarcoma cells data with Salmon. There is a great walkthrough [here](http://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html) that I am following for the most part. The first thing I want to do is get an index for Salmon. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the Salmon index with the command `refgenie pull hg38/salmon_partial_sa_index`. You don't need to do this every time, just when you need to update your index.

The raw data are stored here: `/projects/ptx_results/Sequencing/2020/sequencing20200122_ewingPolysomeSeqRenata`. Rename the files.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/2020/sequencing20200122_ewingPolysomeSeqRenata/"

#previous example filename P14057_1068_S128_R2.fastq.gz
######################################
for i in P14057_1065 P14057_1066 P14057_1067 P14057_1068 P14057_1069 P14057_1070 P14057_1071 P14057_1072 P14057_1073 P14057_1074 P14057_1075 P14057_1076 P14057_1077 P14057_1078 P14057_1079 P14057_1080 P14057_1081 P14057_1082 P14057_1083 P14057_1084 P14057_1085 P14057_1086 P14057_1087 P14057_1088 P14057_1089 P14057_1090 P14057_1091 P14057_1092 P14057_1093 P14057_1094 P14057_1095 P14057_1096 P14057_1097 P14057_1098 P14057_1099 P14057_1100 P14057_1101 P14057_1102
do
  echo $i
  eval "mv ${rawDataInputDirectory}${i}*R1.fastq.gz ${rawDataOutputDirectory}${i}_R1.fastq.gz"
  eval "mv ${rawDataInputDirectory}${i}*R2.fastq.gz ${rawDataOutputDirectory}${i}_R2.fastq.gz"
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
for i in P14057_1065 P14057_1066 P14057_1067 P14057_1068 P14057_1069 P14057_1070 P14057_1071 P14057_1072 P14057_1073 P14057_1074 P14057_1075 P14057_1076 P14057_1077 P14057_1078 P14057_1079 P14057_1080 P14057_1081 P14057_1082 P14057_1083 P14057_1084 P14057_1085 P14057_1086 P14057_1087 P14057_1088 P14057_1089 P14057_1090 P14057_1091 P14057_1092 P14057_1093 P14057_1094 P14057_1095 P14057_1096 P14057_1097 P14057_1098 P14057_1099 P14057_1100 P14057_1101 P14057_1102
do
  echo $i
  bbdukCall="$bbdukLocation in1=${rawDataOutputDirectory}${i}_R1.fastq.gz in2=${rawDataOutputDirectory}${i}_R2.fastq.gz ref=adapters out1=${rawDataOutputDirectory}${i}_R1.clean.fastq.gz out2=${rawDataOutputDirectory}${i}_R2.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo"
  eval $bbdukCall
  ##
  salmonCall="$salmonLocation quant -i $indexLocation -l A -1 ${rawDataOutputDirectory}${i}_R1.clean.fastq.gz -2 ${rawDataOutputDirectory}${i}_R2.clean.fastq.gz --gcBias --validateMappings -o ${rawDataOutputDirectory}${i}_quant"
  eval $salmonCall
done
```

For downstream analysis, file file layout is:

sample.ID	sample	date	condition	RNA	NGI_ID	rep	samplesource	samplesource_simple	combo
P14057_1065	YB1-1	Pca2232 #1	01/22/16	Tx	tot	P14057_1065	1	Pca2232	CHLA10	CHLA10_Tx
P14057_1066	YB1-2	Pca2232 #1	01/22/16	Tx	poly	P14057_1066	1	Pca2232	CHLA10	CHLA10_Tx
P14057_1067	YB1-3	Pca2232 #2	01/22/16	Tx	tot	P14057_1067	2	Pca2232	CHLA10	CHLA10_Tx
P14057_1068	YB1-4	Pca2232 #2	01/22/16	Tx	poly	P14057_1068	2	Pca2232	CHLA10	CHLA10_Tx
P14057_1069	YB1-5	Pca2232 #3	01/22/16	Tx	tot	P14057_1069	3	Pca2232	CHLA10	CHLA10_Tx
P14057_1070	YB1-6	Pca2232 #3	01/22/16	Tx	poly	P14057_1070	3	Pca2232	CHLA10	CHLA10_Tx
P14057_1071	YB1-7	Pca2232 #4	01/22/16	Tx	tot	P14057_1071	4	Pca2232	CHLA10	CHLA10_Tx
P14057_1072	YB1-8	Pca2232 #4	01/22/16	Tx	poly	P14057_1072	4	Pca2232	CHLA10	CHLA10_Tx
P14057_1073	YB1-9	Pca2233 #2	01/22/16	Tx	tot	P14057_1073	2	Pca2233	CHLA10	CHLA10_Tx
P14057_1074	YB1-10	Pca2233 #2	01/22/16	Tx	poly	P14057_1074	2	Pca2233	CHLA10	CHLA10_Tx
P14057_1075	YB1-11	Pca2233 #3	01/22/16	Tx	tot	P14057_1075	3	Pca2233	CHLA10	CHLA10_Tx
P14057_1076	YB1-12	Pca2233 #3	01/22/16	Tx	poly	P14057_1076	3	Pca2233	CHLA10	CHLA10_Tx
P14057_1077	YB1-13	Pca2233 #4	01/22/16	Tx	tot	P14057_1077	4	Pca2233	CHLA10	CHLA10_Tx
P14057_1078	YB1-14	Pca2233 #4	01/22/16	Tx	poly	P14057_1078	4	Pca2233	CHLA10	CHLA10_Tx
P14057_1079	YB1-15	Pca2234 #1	01/22/16	CTL	tot	P14057_1079	1	Pca2234	CHLA10	CHLA10_CTL
P14057_1080	YB1-16	Pca2234 #1	01/22/16	CTL	poly	P14057_1080	1	Pca2234	CHLA10	CHLA10_CTL
P14057_1081	YB1-17	Pca2234 #2	01/22/16	CTL	tot	P14057_1081	2	Pca2234	CHLA10	CHLA10_CTL
P14057_1082	YB1-18	Pca2234 #2	01/22/16	CTL	poly	P14057_1082	2	Pca2234	CHLA10	CHLA10_CTL
P14057_1083	YB1-19	Pca2234 #3	01/22/16	CTL	tot	P14057_1083	3	Pca2234	CHLA10	CHLA10_CTL
P14057_1084	YB1-20	Pca2234 #3	01/22/16	CTL	poly	P14057_1084	3	Pca2234	CHLA10	CHLA10_CTL
P14057_1085	YB1-21	Pca2234 #4	01/22/16	CTL	tot	P14057_1085	4	Pca2234	CHLA10	CHLA10_CTL
P14057_1086	YB1-22	Pca2234 #4	01/22/16	CTL	poly	P14057_1086	4	Pca2234	CHLA10	CHLA10_CTL
P14057_1087	YB1-23	Pca2235 #1	01/22/16	CTL	tot	P14057_1087	1	Pca2235	CHLA10	CHLA10_CTL
P14057_1088	YB1-24	Pca2235 #1	01/22/16	CTL	poly	P14057_1088	1	Pca2235	CHLA10	CHLA10_CTL
P14057_1089	YB1-25	Pca2235 #3	01/22/16	CTL	tot	P14057_1089	3	Pca2235	CHLA10	CHLA10_CTL
P14057_1090	YB1-26	Pca2235 #3	01/22/16	CTL	poly	P14057_1090	3	Pca2235	CHLA10	CHLA10_CTL
P14057_1091	YB1-27	10	03/21/18	CTL	tot	P14057_1091	1	patient	patient	patient_CTL
P14057_1092	YB1-28	10	03/21/18	CTL	poly	P14057_1092	1	patient	patient	patient_CTL
P14057_1093	YB1-29	20	03/21/18	CTL	tot	P14057_1093	2	patient	patient	patient_CTL
P14057_1094	YB1-30	20	03/21/18	CTL	poly	P14057_1094	2	patient	patient	patient_CTL
P14057_1095	YB1-31	30	03/21/18	CTL	tot	P14057_1095	3	patient	patient	patient_CTL
P14057_1096	YB1-32	30	03/21/18	CTL	poly	P14057_1096	3	patient	patient	patient_CTL
P14057_1097	YB1-33	100	03/21/18	Tx	tot	P14057_1097	1	patient	patient	patient_Tx
P14057_1098	YB1-34	100	03/21/18	Tx	poly	P14057_1098	1	patient	patient	patient_Tx
P14057_1099	YB1-35	200	03/21/18	Tx	tot	P14057_1099	2	patient	patient	patient_Tx
P14057_1100	YB1-36	200	03/21/18	Tx	poly	P14057_1100	2	patient	patient	patient_Tx
P14057_1101	YB1-37	300	03/21/18	Tx	tot	P14057_1101	3	patient	patient	patient_Tx
P14057_1102	YB1-38	300	03/21/18	Tx	poly	P14057_1102	3	patient	patient	patient_Tx

Now that we have these data, process them downstream with R.