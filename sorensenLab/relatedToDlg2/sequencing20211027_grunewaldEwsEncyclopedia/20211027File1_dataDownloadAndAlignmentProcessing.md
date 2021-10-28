## Reprocessing some ChIPseq data

This document describes the reprocessing of some ChIPseq data from a previous publication. Specifically:

"The Ewing Sarcoma Cell Line Atlas (ESCLA)"
GEO: GSE176339, https://www.biorxiv.org/content/10.1101/2021.06.08.447518v1

### Getting the raw data

I am interested in the ChIPseq data across the cell lines they profiled. There is a metadata file in the server directory detailed below that describes these data.

I am going to use [SRAExplorer](https://sra-explorer.info/#) to get at these data. I am going to save them in the directory `/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia/`. I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). Because the authors originally used Bowtie for the alignment, I will as well. To this, I am going to use a pre-built index from refgenie following the instructions outlined on [this page](http://refgenie.databio.org/en/latest/install/). To see the list of available indexes, use the command `refgenie listr`. Pull the index with the command `refgenie pull hg38/bowtie2_index`. There is a pretty good walkthrough of the alignment process [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/03_align_and_filtering.html).  

```shell
#!/usr/bin/env bash

###downloading...I had to stop partway through due to a restart, so that is why some are hashed out. I did skip EW24 though because these are paired end for some reason
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/099/SRR14760999/SRR14760999.fastq.gz -o SRR14760999_GSM5363938_A673_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/097/SRR14760997/SRR14760997.fastq.gz -o SRR14760997_GSM5363936_A673_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/098/SRR14760998/SRR14760998.fastq.gz -o SRR14760998_GSM5363937_A673_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/001/SRR14761001/SRR14761001.fastq.gz -o SRR14761001_GSM5363940_A673_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/002/SRR14761002/SRR14761002.fastq.gz -o SRR14761002_GSM5363941_CHLA10_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/000/SRR14761000/SRR14761000.fastq.gz -o SRR14761000_GSM5363939_A673_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/004/SRR14761004/SRR14761004.fastq.gz -o SRR14761004_GSM5363943_CHLA10_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/003/SRR14761003/SRR14761003.fastq.gz -o SRR14761003_GSM5363942_CHLA10_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/006/SRR14761006/SRR14761006.fastq.gz -o SRR14761006_GSM5363945_CHLA10_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/007/SRR14761007/SRR14761007.fastq.gz -o SRR14761007_GSM5363946_CHLA25_ERG_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/008/SRR14761008/SRR14761008.fastq.gz -o SRR14761008_GSM5363947_CHLA25_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/005/SRR14761005/SRR14761005.fastq.gz -o SRR14761005_GSM5363944_CHLA10_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/009/SRR14761009/SRR14761009.fastq.gz -o SRR14761009_GSM5363948_CHLA25_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/010/SRR14761010/SRR14761010.fastq.gz -o SRR14761010_GSM5363949_CHLA25_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/011/SRR14761011/SRR14761011.fastq.gz -o SRR14761011_GSM5363950_CHLA25_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/012/SRR14761012/SRR14761012.fastq.gz -o SRR14761012_GSM5363951_EW1_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/013/SRR14761013/SRR14761013.fastq.gz -o SRR14761013_GSM5363952_EW1_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/014/SRR14761014/SRR14761014.fastq.gz -o SRR14761014_GSM5363953_EW1_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/015/SRR14761015/SRR14761015.fastq.gz -o SRR14761015_GSM5363954_EW1_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/016/SRR14761016/SRR14761016.fastq.gz -o SRR14761016_GSM5363955_EW1_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/018/SRR14761018/SRR14761018.fastq.gz -o SRR14761018_GSM5363957_EW22_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/020/SRR14761020/SRR14761020.fastq.gz -o SRR14761020_GSM5363959_EW22_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/017/SRR14761017/SRR14761017.fastq.gz -o SRR14761017_GSM5363956_EW22_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/019/SRR14761019/SRR14761019.fastq.gz -o SRR14761019_GSM5363958_EW22_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/021/SRR14761021/SRR14761021.fastq.gz -o SRR14761021_GSM5363960_EW22_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/022/SRR14761022/SRR14761022_1.fastq.gz -o SRR14761022_GSM5363961_EW24_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/022/SRR14761022/SRR14761022_2.fastq.gz -o SRR14761022_GSM5363961_EW24_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/023/SRR14761023/SRR14761023_1.fastq.gz -o SRR14761023_GSM5363962_EW24_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/023/SRR14761023/SRR14761023_2.fastq.gz -o SRR14761023_GSM5363962_EW24_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/024/SRR14761024/SRR14761024_1.fastq.gz -o SRR14761024_GSM5363963_EW24_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/024/SRR14761024/SRR14761024_2.fastq.gz -o SRR14761024_GSM5363963_EW24_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/026/SRR14761026/SRR14761026_1.fastq.gz -o SRR14761026_GSM5363965_EW24_input_DNA_Homo_sapiens_ChIP-Seq_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/026/SRR14761026/SRR14761026_2.fastq.gz -o SRR14761026_GSM5363965_EW24_input_DNA_Homo_sapiens_ChIP-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/025/SRR14761025/SRR14761025_1.fastq.gz -o SRR14761025_GSM5363964_EW24_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/025/SRR14761025/SRR14761025_2.fastq.gz -o SRR14761025_GSM5363964_EW24_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/027/SRR14761027/SRR14761027.fastq.gz -o SRR14761027_GSM5363966_EW3_ERG_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/028/SRR14761028/SRR14761028.fastq.gz -o SRR14761028_GSM5363967_EW3_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/029/SRR14761029/SRR14761029.fastq.gz -o SRR14761029_GSM5363968_EW3_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/030/SRR14761030/SRR14761030.fastq.gz -o SRR14761030_GSM5363969_EW3_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/031/SRR14761031/SRR14761031.fastq.gz -o SRR14761031_GSM5363970_EW3_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/032/SRR14761032/SRR14761032.fastq.gz -o SRR14761032_GSM5363971_EW7_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/033/SRR14761033/SRR14761033.fastq.gz -o SRR14761033_GSM5363972_EW7_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/035/SRR14761035/SRR14761035.fastq.gz -o SRR14761035_GSM5363974_EW7_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/034/SRR14761034/SRR14761034.fastq.gz -o SRR14761034_GSM5363973_EW7_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/036/SRR14761036/SRR14761036.fastq.gz -o SRR14761036_GSM5363975_EW7_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/037/SRR14761037/SRR14761037.fastq.gz -o SRR14761037_GSM5363976_MHHES1_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/038/SRR14761038/SRR14761038.fastq.gz -o SRR14761038_GSM5363977_MHHES1_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/041/SRR14761041/SRR14761041.fastq.gz -o SRR14761041_GSM5363980_MHHES1_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/039/SRR14761039/SRR14761039.fastq.gz -o SRR14761039_GSM5363978_MHHES1_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/040/SRR14761040/SRR14761040.fastq.gz -o SRR14761040_GSM5363979_MHHES1_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/042/SRR14761042/SRR14761042.fastq.gz -o SRR14761042_GSM5363981_MIC_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/044/SRR14761044/SRR14761044.fastq.gz -o SRR14761044_GSM5363983_MIC_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/043/SRR14761043/SRR14761043.fastq.gz -o SRR14761043_GSM5363982_MIC_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/045/SRR14761045/SRR14761045.fastq.gz -o SRR14761045_GSM5363984_MIC_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/047/SRR14761047/SRR14761047.fastq.gz -o SRR14761047_GSM5363986_POE_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/046/SRR14761046/SRR14761046.fastq.gz -o SRR14761046_GSM5363985_MIC_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/048/SRR14761048/SRR14761048.fastq.gz -o SRR14761048_GSM5363987_POE_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/049/SRR14761049/SRR14761049.fastq.gz -o SRR14761049_GSM5363988_POE_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/050/SRR14761050/SRR14761050.fastq.gz -o SRR14761050_GSM5363989_POE_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/055/SRR14761055/SRR14761055.fastq.gz -o SRR14761055_GSM5363993_RDES_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/051/SRR14761051/SRR14761051.fastq.gz -o SRR14761051_GSM5363990_POE_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/053/SRR14761053/SRR14761053.fastq.gz -o SRR14761053_GSM5363992_RDES_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/052/SRR14761052/SRR14761052.fastq.gz -o SRR14761052_GSM5363991_RDES_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/057/SRR14761057/SRR14761057.fastq.gz -o SRR14761057_GSM5363995_RDES_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/056/SRR14761056/SRR14761056.fastq.gz -o SRR14761056_GSM5363994_RDES_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/058/SRR14761058/SRR14761058.fastq.gz -o SRR14761058_GSM5363996_RH1_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/060/SRR14761060/SRR14761060.fastq.gz -o SRR14761060_GSM5363998_RH1_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/059/SRR14761059/SRR14761059.fastq.gz -o SRR14761059_GSM5363997_RH1_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/061/SRR14761061/SRR14761061.fastq.gz -o SRR14761061_GSM5363999_RH1_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/063/SRR14761063/SRR14761063.fastq.gz -o SRR14761063_GSM5364001_SKES1_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/062/SRR14761062/SRR14761062.fastq.gz -o SRR14761062_GSM5364000_RH1_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/064/SRR14761064/SRR14761064.fastq.gz -o SRR14761064_GSM5364002_SKES1_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/065/SRR14761065/SRR14761065.fastq.gz -o SRR14761065_GSM5364003_SKES1_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/066/SRR14761066/SRR14761066.fastq.gz -o SRR14761066_GSM5364004_SKES1_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/067/SRR14761067/SRR14761067.fastq.gz -o SRR14761067_GSM5364005_SKES1_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/069/SRR14761069/SRR14761069.fastq.gz -o SRR14761069_GSM5364007_SKNMC_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/068/SRR14761068/SRR14761068.fastq.gz -o SRR14761068_GSM5364006_SKNMC_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/070/SRR14761070/SRR14761070.fastq.gz -o SRR14761070_GSM5364008_SKNMC_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/073/SRR14761073/SRR14761073.fastq.gz -o SRR14761073_GSM5364011_TC32_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/072/SRR14761072/SRR14761072.fastq.gz -o SRR14761072_GSM5364010_SKNMC_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/074/SRR14761074/SRR14761074.fastq.gz -o SRR14761074_GSM5364012_TC32_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/071/SRR14761071/SRR14761071.fastq.gz -o SRR14761071_GSM5364009_SKNMC_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/075/SRR14761075/SRR14761075.fastq.gz -o SRR14761075_GSM5364013_TC32_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/077/SRR14761077/SRR14761077.fastq.gz -o SRR14761077_GSM5364015_TC32_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/076/SRR14761076/SRR14761076.fastq.gz -o SRR14761076_GSM5364014_TC32_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/078/SRR14761078/SRR14761078.fastq.gz -o SRR14761078_GSM5364016_TC71_FLI1_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/080/SRR14761080/SRR14761080.fastq.gz -o SRR14761080_GSM5364018_TC71_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/079/SRR14761079/SRR14761079.fastq.gz -o SRR14761079_GSM5364017_TC71_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/081/SRR14761081/SRR14761081.fastq.gz -o SRR14761081_GSM5364019_TC71_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/083/SRR14761083/SRR14761083.fastq.gz -o SRR14761083_GSM5364021_TC106_ERG_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/084/SRR14761084/SRR14761084.fastq.gz -o SRR14761084_GSM5364022_TC106_H3K4me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/082/SRR14761082/SRR14761082.fastq.gz -o SRR14761082_GSM5364020_TC71_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/085/SRR14761085/SRR14761085.fastq.gz -o SRR14761085_GSM5364023_TC106_H3K27me3_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/086/SRR14761086/SRR14761086.fastq.gz -o SRR14761086_GSM5364024_TC106_H3K27ac_ChIPSeq_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/087/SRR14761087/SRR14761087.fastq.gz -o SRR14761087_GSM5364025_TC106_input_DNA_Homo_sapiens_ChIP-Seq.fastq.gz



######processing of the files...these are the tools we need
bbdukLocation="/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
bowtieLocation="/home/chughes/softwareTools/bowtie2-2.4.4/bowtie2"
referenceLocation="/home/chughes/databases/refgenieManualGenomes/hg38/bowtie2Index_072021/default/"
annotationLocation="/home/chughes/databases/refgenieManualGenomes/hg38/gencodeGtf_072021/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.gtf"
samtoolsLocation="/home/chughes/softwareTools/samtools-1.12/samtools"
sambambaLocation="/home/chughes/softwareTools/sambamba-0.8.1/sambamba"
rawDataOutputDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia/"

##########this is the processing code
for i in SRR1476{0997..1035}
do
  echo $i
  eval "cd ${rawDataOutputDirectory}"
  
  ##bbduk
  bbdukCall="$bbdukLocation in=${rawDataOutputDirectory}${i}*Seq.fastq.gz ref=adapters out=${rawDataOutputDirectory}${i}.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1"
  eval $bbdukCall

  ##bowtie
  bowtieCall="$bowtieLocation -p 8 -q --local -x $referenceLocation -U ${rawDataOutputDirectory}${i}.clean.fastq.gz -S ${rawDataOutputDirectory}${i}.unsorted.sam"
  eval $bowtieCall

  ##bamconvert with samtools
  bamConvertCall="$samtoolsLocation view -h -S -b -o ${rawDataOutputDirectory}${i}.unsorted.bam ${rawDataOutputDirectory}${i}.unsorted.sam"
  eval $bamConvertCall

  ##sort the bam file with sambamba
  bamSortCall="$sambambaLocation sort -t 6 -o ${rawDataOutputDirectory}${i}.sorted.bam ${rawDataOutputDirectory}${i}.unsorted.bam"
  eval $bamSortCall

  ##perform filtering on the files to remove duplicates
  eval $sambambaLocation view -h -t 6 -f bam -F '"[XS] == null and not unmapped and not duplicate"' ${rawDataOutputDirectory}${i}.sorted.bam > ${rawDataOutputDirectory}${i}.filtered.bam

  ##index the bam files
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}${i}.filtered.bam"
  eval $bamIndexCall
done
```




This script didn't completely work. The removing of duplicates command wasn't working, I think because of a text error. It just didn't like the quotes in the command. So, I had to go through the files manually just because I didn't really want to get into a big thing about fixing it. 

```shell
/projects/ptx_analysis/chughes/software/sambamba-0.8.0/sambamba-0.8.0 view -h -t 6 -f bam -F '[XS] == null and not unmapped and not duplicate' ./SRR8832667.sorted.bam > ./SRR8832667.filtered.bam
```

Now we move on to peak calling using [MACS](https://github.com/macs3-project/MACS). There is a great walkthrough of this process [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html). Also [here](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) for discussion of deepTools, specifically bamCoverage.

```shell
#!/bin/bash
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/"
for i in SRR8832669 SRR8832670 SRR8832671 SRR8832672 SRR8832673 SRR8832674
do
  echo $i
  eval macs3 callpeak -t ${rawDataOutputDirectory}${i}.filtered.bam -c ${rawDataOutputDirectory}SRR8832666.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 
done
```

For the H3K27ac data, I rand MACS3 using the --broad tag as well as this is how it was done in the original manuscript. Now we are ready for visualization. There is a good walkthrough of visualization [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html). First, I ran [deeptools](https://deeptools.readthedocs.io/en/develop/) on the bam files to create coverage maps:

```shell
bamCoverage -b ./SRR8832666.filtered.bam -o ./deeptools/SRR8832666.chr11.bw --binSize 20 --region chr11 --normalizeUsing BPM --smoothLength 60 --extendReads 150 --centerReads -p 6

bamCoverage -b ./$i.filtered.bam -o ./deeptools/$i.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --extendReads 150 --centerReads -p 6
```

Now, I want to create a heatmap around the TSS for genes on chromosome 11. First I need to create a bed file with the chr11 gene regions in it. I can do this using the command line:

```shell
cat /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf |  awk 'OFS="\t" {if ($3=="gene" && $1=="chr11") {print $1,$4-1,$5,$10,$16,$7}}' | tr -d '";' > ./chr11Hg38GeneRegions.bed
```

Now I can run computeMatrix from deeptools.

```shell
computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/chr11Hg38GeneRegions.bed -S /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/SRR88326[67][901234]*.bw --skipZeros -o /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/ewsFli1ChipSeqTssMatrixChr11.gz -p 6 --outFileSortedRegions /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/regionsTssChr11.bed
```

We can create a profile plot from these data.

```shell
plotProfile -m /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/ewsFli1ChipSeqTssMatrixChr11.gz --outFileName /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/regionsTssChr11Profile.pdf --perGroup --refPointLabel "TSS"
```

Alternatively, we can show this as a heatmap.

```shell
plotHeatmap -m /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/ewsFli1ChipSeqTssMatrixChr11.gz --outFileName /projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/deeptools/regionsTssChr11Heatmap.pdf --colorMap RdBu
```

I think for the rest of the visualization, we will carry this out in R.

