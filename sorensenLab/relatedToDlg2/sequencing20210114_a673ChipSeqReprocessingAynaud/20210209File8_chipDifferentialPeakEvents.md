## Differential peak calling from ChIPseq data

This document describes analysis of the MACS3 peak data we already have to call for differential binding events. 

### Setting up the analysis

For this work, we are going to use the timecourse ChIPseq data from the Aynaud paper. Because there is a lag in EWS-FLI1 re-expression, we can use the day 7 and 9 files as replicates, and the d14-d17 as well.

```
Control
SRX5620946: GSM3701309: ASP14_d0_input; Homo sapiens; ChIP-Seq
SRR8832666

Treatment 1
SRX5620949: GSM3701312: ASP14_d7_FLI1; Homo sapiens; ChIP-Seq
SRR8832669

SRX5620950: GSM3701313: ASP14_d9_FLI1; Homo sapiens; ChIP-Seq
SRR8832670

Treatment 2
SRX5620953: GSM3701316: ASP14_d14_FLI1; Homo sapiens; ChIP-Seq
SRR8832673

SRX5620954: GSM3701317: ASP14_d17_FLI1; Homo sapiens; ChIP-Seq
SRR8832674
```

I am mostly going to follow a guide on the [MACS GitHub](https://github.com/macs3-project/MACS/wiki/Call-differential-binding-events). From the wiki:

```
Purpose of this step is to use callpeak with -B option to generate bedGraph files for both conditions. There are several things to be remember: 1. --SPMR is not compatible with bdgdiff, so avoid using it; 2. prepare a pen to write down the number of non-redundant reads of both conditions -- you will find such information in runtime message or xls output from callpeak; 3. keep using the same --extsize for both conditions (you can get it from predictd module).
```

Run the commands to get the extsize for our files. It actually looks like we can just use two files, so we don't necessarily need to use the border timepoints as replicates. We can revisit that later.

```shell
#!/bin/bash
processingDir='/projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/'

###
for i in SRR8832669 SRR8832674
do
eval echo $i
eval "macs3 predictd -i ${processingDir}${i}.filtered.bam --outdir ${processingDir}diffPeakCalling"
done
```

From these commands, the fragment size for the first file was 175bp, and for the second, 188bp. The average of these is 181. This is what we will use for extsize.

```shell
#!/bin/bash
processingDir='/projects/ptx_results/Sequencing/publishedStudies/202002AynaudCellReportsPmid32049009/chipSeq/'

###
for i in SRR8832669 SRR8832674
do
eval echo $i
eval "macs3 callpeak -B -t ${processingDir}${i}.filtered.bam -c ${processingDir}SRR8832666.filtered.bam -n ${i} --nomodel --extsize 180 --outdir ${processingDir}diffPeakCalling"
done
```

We need to pay attention to the output of this command to get the tag numbers that will be useful later on, according to the wiki page. These data are:

```
SRR8832669
tags after filtering in treatment: 17948389
tags after filtering in control: 26621844

SRR8832674
tags after filtering in treatment: 17182635
tags after filtering in control: 26621844
```

Or you can call these data from the output using the command `egrep "tags after filtering in treatment|tags after filtering in control" cond1_peaks.xls`. Now we are ready for the differential peak calling analysis.

```shell
macs3 bdgdiff --t1 SRR8832674_treat_pileup.bdg --c1 SRR8832674_control_lambda.bdg --t2 SRR8832669_treat_pileup.bdg --c2 SRR8832669_control_lambda.bdg --d1 26621844 --d2 26621844 -g 60 -l 120 --o-prefix diff_d17_vs_d7
```

The challenge here is that we have the same control for both samples, which is not ideal, but it will have to do for now. Description of the output from the wiki page is:

```
You will get the following three files in working directory:

diff_c1_vs_c2_c3.0_cond1.bed
This file stores regions that are highly enriched in condition 1 comparing to condition 2. The last column in the file represent the log10 likelihood ratio to show how likely the observed signal in condition 1 in this region is from condition 1 comparing to condition 2. Higher the value, bigger the difference.

diff_c1_vs_c2_c3.0_cond2.bed
This file stores regions that are highly enriched in condition 2 comparing to condition 1. The last column in the file represent the log10 likelihood ratio to show how likely the observed signal in condition 2 in this region is from condition 2 comparing to condition 1. Higher the value, bigger the difference.

diff_c1_vs_c2_c3.0_common.bed
This file stores regions that are highly enriched in both condition 1 and condition 2, and the difference between condition 1 and condition 2 is not significant. The last column in the file represent the difference between condition 1 and condition 2 in log10 likelihood ratios.
```

So, I think the first one is the one I am most interested in. We will process these data further in R. 

Actually, I had to come back to this because I need some coverage data for the metagene portion of the analysis now that I have called and filtered the peaks. It looks like deeptools might be the best option for this. If you look in the `/projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/` directory where the A673 and SKNMC data are stored, I used a script named `bigwigCoverageRefseqAll.sh` to process these data into bigwig files. 

It looks like we can use computeMatrix from deeptools to make a metagene plot. It needs a bed file or a GTF as an annotation input though. 


```shell
grep 'DLG2' /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/refseq_gtf/default/hg38.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/dlg2Regions.gtf


computeMatrix scale-regions -R dlg2Regions.gtf -S SRR8616012_refseq.all.bw -b 3000 -a 3000 regionBodyLength 5000 --skipZeros --metagene -o dlg2Metagene.gz --outFileNameMatrix dlg2Metagene



grep 'ENST00000426717.6\|ENST00000398309.6' /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf
awk '!/start_codon/ && !/stop_codon/ && !/UTR/ && !/^CDS$/' /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf
```



deepTools2.0/bin/computeMatrix scale-regions \
  -R genes_chr19_firstHalf.bed genes_chr19_secondHalf.bed \ # separate multiple files with spaces
  -S testFiles/log2ratio_*.bw  \ or use the wild card approach
  -b 3000 -a 3000 \
  --regionBodyLength 5000 \
  --skipZeros -o matrix2_multipleBW_l2r_twoGroups_scaled.gz \
  --outFileNameMatrix matrix2_multipleBW_l2r_twoGroups_scaled.tab \
  --outFileSortedRegions regions2_multipleBW_l2r_twoGroups_genes.bed