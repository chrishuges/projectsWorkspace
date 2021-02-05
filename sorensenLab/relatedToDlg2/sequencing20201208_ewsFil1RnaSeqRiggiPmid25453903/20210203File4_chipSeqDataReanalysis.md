## Reprocessing some ChIPseq data

This document describes the reprocessing of some ChIPseq data from a previous publication. Specifically:

"EWS-FLI1 utilizes divergent chromatin remodeling mechanisms to directly activate or repress enhancer elements in Ewing sarcoma"
Cancer Cell, 2014, Pubmed ID: 25453903, GEO: GSE61953

### Getting the raw data

I am interested in all of the Ewing sarcoma data now. I had originally done the A673 and SKNMC files, but I will expand it to include the other lines now. the files are:

```
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/005/SRR1593965/SRR1593965.fastq.gz -o SRR1593965_GSM1517542_WDR5_ChIP-seq_in_SKMNC_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/003/SRR1593963/SRR1593963.fastq.gz -o SRR1593963_GSM1517540_H3K27me3_ChIP-seq_in_SKMNC_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/001/SRR1593961/SRR1593961.fastq.gz -o SRR1593961_GSM1517538_H3K27ac_ChIP-seq_in_SKMNC_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/004/SRR1593964/SRR1593964.fastq.gz -o SRR1593964_GSM1517541_H3K4me3_ChIP-seq_in_SKMNC_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/002/SRR1593962/SRR1593962.fastq.gz -o SRR1593962_GSM1517539_H3K4me1_ChIP-seq_in_SKMNC_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/000/SRR1593960/SRR1593960.fastq.gz -o SRR1593960_GSM1517537_FLI1_ChIP-seq_in_SKMNC_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/007/SRR1593967/SRR1593967.fastq.gz -o SRR1593967_GSM1517544_FLI1_ChIP-seq_in_SKMNC_48hrs_shGFP_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/009/SRR1593969/SRR1593969.fastq.gz -o SRR1593969_GSM1517546_FLI1_ChIP-seq_in_SKMNC_96hrs_shGFP96_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/000/SRR1593970/SRR1593970.fastq.gz -o SRR1593970_GSM1517547_H3K27ac_ChIP-seq_in_SKMNC_96hrs_shGFP_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/006/SRR1593966/SRR1593966.fastq.gz -o SRR1593966_GSM1517543_SKMNC_whole_cell_extract_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/008/SRR1593968/SRR1593968.fastq.gz -o SRR1593968_GSM1517545_H3K27ac_ChIP-seq_in_SKMNC_48hrs_shGFP_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/001/SRR1593971/SRR1593971.fastq.gz -o SRR1593971_GSM1517548_H3K4me1_ChIP-seq_in_SKMNC_96hrs_shGFP_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/004/SRR1593974/SRR1593974.fastq.gz -o SRR1593974_GSM1517551_ELF1_ChIP-seq_in_SKMNC_96_hrs_shGFP_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/002/SRR1593972/SRR1593972.fastq.gz -o SRR1593972_GSM1517549_H3K27me3_ChIP-seq_in_SKMNC_96hrs_shGFP_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/003/SRR1593973/SRR1593973.fastq.gz -o SRR1593973_GSM1517550_GABPA_ChIP-seq_in_SKMNC_96hrs_shGFP_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/005/SRR1593975/SRR1593975.fastq.gz -o SRR1593975_GSM1517552_p300_ChIP-seq_in_SKMNC_96_hrs_shGFP_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/000/SRR1593980/SRR1593980.fastq.gz -o SRR1593980_GSM1517557_H3K4me1_ChIP-seq_in_SKMNC_96hrs_shFLI1_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/001/SRR1593981/SRR1593981.fastq.gz -o SRR1593981_GSM1517558_H3K27me3_ChIP-seq_in_SKMNC_96hrs_shFLI1_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/006/SRR1593976/SRR1593976.fastq.gz -o SRR1593976_GSM1517553_FLI1_ChIP-seq_in_SKMNC_48hrs_shFLI1_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/007/SRR1593977/SRR1593977.fastq.gz -o SRR1593977_GSM1517554_H3K27ac_ChIP-seq_in_SKMNC_48hrs_shFLI1_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/008/SRR1593978/SRR1593978.fastq.gz -o SRR1593978_GSM1517555_FLI1_ChIP-seq_in_SKMNC_96hrs_shFLI1_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/002/SRR1593982/SRR1593982.fastq.gz -o SRR1593982_GSM1517559_GABPA_ChIP-seq_in_SKMNC_96hrs_shFLI1_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/009/SRR1593979/SRR1593979.fastq.gz -o SRR1593979_GSM1517556_H3K27ac_ChIP-seq_in_SKMNC_96hrs_shFLI1_Homo_sapiens_ChIP-Seq.fastq.gz



curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/006/SRR1593986/SRR1593986.fastq.gz -o SRR1593986_GSM1517563_H3K27ac_ChIP-Seq_in_A673_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/008/SRR1593988/SRR1593988.fastq.gz -o SRR1593988_GSM1517565_H3K27me3_ChIP-Seq_in_A673_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/009/SRR1593989/SRR1593989.fastq.gz -o SRR1593989_GSM1517566_H3K4me3_ChIP-Seq_in_A673_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/000/SRR1593990/SRR1593990.fastq.gz -o SRR1593990_GSM1517567_WDR5_ChIP-seq_in_A673_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/001/SRR1593991/SRR1593991.fastq.gz -o SRR1593991_GSM1517568_A673_whole_cell_extract_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/003/SRR1593983/SRR1593983.fastq.gz -o SRR1593983_GSM1517560_ELF1_ChIP-seq_in_SKMNC_96hrs_shFLI1_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/002/SRR1593992/SRR1593992.fastq.gz -o SRR1593992_GSM1517569_FLI1_ChIP-seq_in_A673_48hrs_shGFP_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/005/SRR1593985/SRR1593985.fastq.gz -o SRR1593985_GSM1517562_FLI1_ChIP-Seq_in_A673_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/004/SRR1593984/SRR1593984.fastq.gz -o SRR1593984_GSM1517561_p300_ChIP-seq_in_SKMNC_96hrs_shFLI1_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/007/SRR1593987/SRR1593987.fastq.gz -o SRR1593987_GSM1517564_H3K4me1_ChIP-Seq_in_A673_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/003/SRR1593993/SRR1593993.fastq.gz -o SRR1593993_GSM1517570_H3K27ac_ChIP-seq_in_A673_48hrs_shGFP_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/004/SRR1593994/SRR1593994.fastq.gz -o SRR1593994_GSM1517571_H3K27ac_ChIP-seq_in_A673_96hrs_shGFP_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/006/SRR1593996/SRR1593996.fastq.gz -o SRR1593996_GSM1517573_H3K27ac_ChIP-seq_in_A673_48hrs_shFLI1_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/007/SRR1593997/SRR1593997.fastq.gz -o SRR1593997_GSM1517574_H3K27ac_ChIP-seq_in_A673_96hrs_shFLI1_Homo_sapiens_ChIP-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/005/SRR1593995/SRR1593995.fastq.gz -o SRR1593995_GSM1517572_FLI1_ChIP-seq_in_A673_48hrs_shFLI1_Homo_sapiens_ChIP-Seq.fastq.gz
```

I am going to use [SRAExplorer](https://sra-explorer.info/#) to get at these data. I am going to save them in the directory `/projects/ptx_results/Sequencing/publishedStudies/201411RiggiCancerCellPmid25453903/chipSeq`. I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). 

```shell
#!/bin/bash
bbdukLocation="/projects/ptx_analysis/chughes/software/bbmap_v38_87/bbduk.sh"
bowtieLocation="/gsc/software/linux-x86_64-centos7/bowtie2-2.3.4.1/bin/bowtie2"
referenceLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/bowtie2_index/default/hg38"
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"
sambambaLocation="/projects/ptx_analysis/chughes/software/sambamba-0.8.0/sambamba-0.8.0"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201411RiggiCancerCellPmid25453903/chipSeq/"


#download the raw data files
eval "cd ${rawDataOutputDirectory}"
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/005/SRR1593965/SRR1593965.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/003/SRR1593963/SRR1593963.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/001/SRR1593961/SRR1593961.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/004/SRR1593964/SRR1593964.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/002/SRR1593962/SRR1593962.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/000/SRR1593960/SRR1593960.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/007/SRR1593967/SRR1593967.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/009/SRR1593969/SRR1593969.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/000/SRR1593970/SRR1593970.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/006/SRR1593966/SRR1593966.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/008/SRR1593968/SRR1593968.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/001/SRR1593971/SRR1593971.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/004/SRR1593974/SRR1593974.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/002/SRR1593972/SRR1593972.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/003/SRR1593973/SRR1593973.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/005/SRR1593975/SRR1593975.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/000/SRR1593980/SRR1593980.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/001/SRR1593981/SRR1593981.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/006/SRR1593976/SRR1593976.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/007/SRR1593977/SRR1593977.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/008/SRR1593978/SRR1593978.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/002/SRR1593982/SRR1593982.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/009/SRR1593979/SRR1593979.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/006/SRR1593986/SRR1593986.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/008/SRR1593988/SRR1593988.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/009/SRR1593989/SRR1593989.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/000/SRR1593990/SRR1593990.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/001/SRR1593991/SRR1593991.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/003/SRR1593983/SRR1593983.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/002/SRR1593992/SRR1593992.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/005/SRR1593985/SRR1593985.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/004/SRR1593984/SRR1593984.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/007/SRR1593987/SRR1593987.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/003/SRR1593993/SRR1593993.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/004/SRR1593994/SRR1593994.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/006/SRR1593996/SRR1593996.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/007/SRR1593997/SRR1593997.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/005/SRR1593995/SRR1593995.fastq.gz

#process the downloaded files
for i in SRR1593965 SRR1593963 SRR1593961 SRR1593964 SRR1593962 SRR1593960 SRR1593967 SRR1593969 SRR1593970 SRR1593966 SRR1593968 SRR1593971 SRR1593974 SRR1593972 SRR1593973 SRR1593975 SRR1593980 SRR1593981 SRR1593976 SRR1593977 SRR1593978 SRR1593982 SRR1593979 SRR1593986 SRR1593988 SRR1593989 SRR1593990 SRR1593991 SRR1593983 SRR1593992 SRR1593985 SRR1593984 SRR1593987 SRR1593993 SRR1593994 SRR1593996 SRR1593997 SRR1593995
do
  echo $i
  bbdukCall="$bbdukLocation in=${rawDataOutputDirectory}${i}.fastq.gz ref=adapters out=${rawDataOutputDirectory}${i}.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1"
  eval $bbdukCall
  bowtieCall="$bowtieLocation -p 8 -q --local -x $referenceLocation -U ${rawDataOutputDirectory}${i}.clean.fastq.gz -S ${rawDataOutputDirectory}${i}.unsorted.sam"
  eval $bowtieCall
  bamConvertCall="$samtoolsLocation view -h -S -b -o ${rawDataOutputDirectory}${i}.unsorted.bam ${rawDataOutputDirectory}${i}.unsorted.sam"
  eval $bamConvertCall
  bamSortCall="$sambambaLocation sort -t 6 -o ${rawDataOutputDirectory}${i}.sorted.bam ${rawDataOutputDirectory}${i}.unsorted.bam"
  eval $bamSortCall
  eval $sambambaLocation view -h -t 6 -f bam -F '"[XS] == null and not unmapped and not duplicate"' ${rawDataOutputDirectory}${i}.sorted.bam > ${rawDataOutputDirectory}${i}.filtered.bam
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}${i}.filtered.bam"
  eval $bamIndexCall
done
```

The run command in a screen session was:

```shell
./downloadAndAlignment.sh |& tee -a downloadAndAlignmentLog.txt
```

The next thing we need to do is the peak calling and the coverage mapping. We will do this with MACS as we have done before. This is the SKNMC histone data.

```shell
#!/bin/bash
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201411RiggiCancerCellPmid25453903/chipSeq/"
for i in SRR1593963 SRR1593961 SRR1593964 SRR1593962 SRR1593970 SRR1593968 SRR1593971 SRR1593972 SRR1593980 SRR1593981 SRR1593977 SRR1593979
do
  echo $i
  eval macs3 callpeak -t ${rawDataOutputDirectory}${i}.filtered.bam -c ${rawDataOutputDirectory}SRR1593966.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done
```

These are the A673 histone data.

```shell
#!/bin/bash
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201411RiggiCancerCellPmid25453903/chipSeq/"
for i in SRR1593986 SRR1593988 SRR1593989 SRR1593987 SRR1593993 SRR1593994 SRR1593996 SRR1593997
do
  echo $i
  eval macs3 callpeak -t ${rawDataOutputDirectory}${i}.filtered.bam -c ${rawDataOutputDirectory}SRR1593991.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done
```

This is the SKNMC target protein data.

```shell
#!/bin/bash
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201411RiggiCancerCellPmid25453903/chipSeq/"
for i in SRR1593965 SRR1593960 SRR1593967 SRR1593969 SRR1593974 SRR1593973 SRR1593975 SRR1593976 SRR1593978 SRR1593982
do
  echo $i
  eval macs3 callpeak -t ${rawDataOutputDirectory}${i}.filtered.bam -c ${rawDataOutputDirectory}SRR1593966.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01
done
```

This is the A673 target protein data

```shell
#!/bin/bash
annotationLocation="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201411RiggiCancerCellPmid25453903/chipSeq/"
for i in SRR1593990 SRR1593983 SRR1593992 SRR1593985 SRR1593984 SRR1593995
do
  echo $i
  eval macs3 callpeak -t ${rawDataOutputDirectory}${i}.filtered.bam -c ${rawDataOutputDirectory}SRR1593991.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01
done
```

Script to call them all.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201411RiggiCancerCellPmid25453903/chipSeq/"
eval "cd ${rawDataOutputDirectory}"
eval "./macsPeakCall_sknmcHistones.sh |& tee macsPeakCall_sknmcHistonesLog.txt"
eval "./macsPeakCall_a673Histones.sh |& tee macsPeakCall_a673HistonesLog.txt"
eval "./macsPeakCall_sknmcTargets.sh |& tee macsPeakCall_sknmcTargetsLog.txt"
eval "./macsPeakCall_a673Targets.sh |& tee macsPeakCall_a673TargetsLog.txt"
```

Last thing to do is to run deeptools to get the coverage values.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/201411RiggiCancerCellPmid25453903/chipSeq/"
for i in SRR1593965 SRR1593963 SRR1593961 SRR1593964 SRR1593962 SRR1593960 SRR1593967 SRR1593969 SRR1593970 SRR1593966 SRR1593968 SRR1593971 SRR1593974 SRR1593972 SRR1593973 SRR1593975 SRR1593980 SRR1593981 SRR1593976 SRR1593977 SRR1593978 SRR1593982 SRR1593979 SRR1593986 SRR1593988 SRR1593989 SRR1593990 SRR1593991 SRR1593983 SRR1593992 SRR1593985 SRR1593984 SRR1593987 SRR1593993 SRR1593994 SRR1593996 SRR1593997 SRR1593995
do
  echo $i
  eval bamCoverage -b ./$i.filtered.bam -o ./$i.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --extendReads 150 --centerReads -p 6
done
```

That is finished. I will use these data files in a downstream Rmd analysis file. 





