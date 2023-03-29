## Processing Haifeng CRISPR screen data

This document describes the reprocessing of some CRISPR screen data from Haifeng.

### Description

I am basically following the protocol described [here](https://www.nature.com/articles/s41596-018-0113-7). First I installed the software on sorensenlab server.

```shell
conda create -n mageck-vispr mageck mageck-vispr python=3
conda activate mageck-vispr
```

Run mageck count.

```shell
mageck count -l library.csv -n hzMetabo --sample-label d0_a,a673_r1,a673_r2,a673_r3,d0_r,b20r_r1,b20r_r2,b20r_r3 --fastq Ao_S1_L001_R1_001.fastq.gz A673-1_S3_L001_R1_001.fastq.gz A673-2_S4_L001_R1_001.fastq.gz A673-3_S5_L001_R1_001.fastq.gz Ro_S2_L001_R1_001.fastq.gz B20R-1_S9_L001_R1_001.fastq.gz B20R-2_S10_L001_R1_001.fastq.gz B20R-3_S11_L001_R1_001.fastq.gz


##other library
mageck count -l library.csv -n hzMetabo --sample-label d0_a,a673_r1,a673_r2,a673_r3,d0_r,cb8_r1,cb8_r2,cb8_r3 --fastq Ao_S1_L001_R1_001.fastq.gz A673-1_S3_L001_R1_001.fastq.gz A673-2_S4_L001_R1_001.fastq.gz A673-3_S5_L001_R1_001.fastq.gz Ro_S2_L001_R1_001.fastq.gz CB4-1_S6_L001_R1_001.fastq.gz CB4-2_S7_L001_R1_001.fastq.gz CB4-3_S8_L001_R1_001.fastq.gz
```

Run the mageck test command to get the dependency genes for comparisons.

```shell
mageck test -k hzMetabo.count.txt -t a673_r1,a673_r2,a673_r3 -c d0_a -n hzMetaboA673_rra --remove-zero both --remove-zero-threshold 0
mageck test -k hzMetabo.count.txt -t b20r_r1,b20r_r2,b20r_r3 -c d0_r -n hzMetaboB20r_rra --remove-zero both --remove-zero-threshold 0
mageck test -k hzMetabo.count.txt -t cb8_r1,cb8_r2,cb8_r3 -c d0_r -n hzMetaboCb8_rra --remove-zero both --remove-zero-threshold 0




##the code below here I don't think I ever used
mageck test -k GSC_0131.count.txt -t day23_r1,day23_r2 -c
day0_r1,day0_r2 -n GSC_0131_rra --remove-zero both --removezero-threshold 0



mageck mle --count-table hzMetabo.count.txt --design-matrix designmatrix.txt --norm-method control --control-sgrna nonessential_ctrl_sgrna_
list.txt --output-prefix braf.mle

mageck mle --count-table rawcount.txt --design-matrix designmatrix.
txt --norm-method control --control-sgrna nonessential_ctrl_sgrna_
list.txt --output-prefix braf.mle

```