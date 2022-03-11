## Running bigwig file comparison

This document describes the running the [multiBigWigSummary](https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html) tool for Dlg2.

## Workflow

The first thing I did was to follow the installation instructions. First I am going to try running the tool on a set of 5 cell lines, the first 5 from CCLE. In order to run the tool, I need a BED file that contains my regions of interest.

```shell
grep "DLG2[^-AS2]" /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211112_ccleSequencingDataReprocessing/multiBigWig/dlg2Only.gtf
```

This gives me all of the DLG2 exons across all isoforms as well as other things, like complete coding regions and UTR's. I am not sure if I want to subset it. Perhaps I can just run the tool and see what happens.

```shell
multiBigwigSummary BED-file --metagene --bwfiles ../results/*.bw --BED dlg2Only.gtf --outFileName dlg2OnlyGtfScores.npz --outRawCounts dlg2OnlyGtfScores.tab
```

That was fast. I think the GTF is not the best way to do this. It seems to not like when som efeatures are repeated, which is understandable. 




```shell
computeMatrix scale-regions -S ../results/SRR8616112.sorted.bw -R dlg2Only.gtf -bs 20 -m 10000 -b 3000 -a 3000 -out testOut.tab.gz --skipZeros --missingDataAsZero

plotHeatmap -m testOut.tab.gz --outFileName testOut.png --colorMap YlGnBu --regionsLabel dlg2 --heatmapHeight 15 --plotTitle 'DLG2' &
```



$ computeMatrix scale-regions \
 -S GCcontent_Mm9_50_5.bw \
 -R RefSeq_genes_uniqNM.bed \
 -bs 50
 -m 10000 -b 3000 -a 3000 \
 -out matrix_GCcont_Mm9_scaledGenes.tab.gz \
 --skipZeros \
 --missingDataAsZero

$ computeMatrix scale-regions \
 -S GCcontent_Dm3_50_5.bw \
 -R Dm530.genes.bed \
 -bs 50
 -m 3000 -b 1000 -a 1000 \
 -out matrix_GCcont_Dm3_scaledGenes.tab.gz \
 --skipZeros --missingDataAsZero

$ plotHeatmap \
 -m matrix_GCcont_Dm3_scaledGenes.tab.gz \
 -out hm_GCcont_Dm3_scaledGenes.png \
 --colorMap YlGnBu \
 --regionsLabel 'fly genes' \
 --heatmapHeight 15 \
 --plotTitle 'GC content fly' &

$ plotHeatmap \
 -m matrix_GCcont_Mm9_scaledGenes.tab.gz \
 -out hm_GCcont_Mm9_scaledGenes.png \
 --colorMap YlGnBu \
 --regionsLabel 'mouse genes' \
 --heatmapHeight 15 \
 --plotTitle 'GC content mouse' &












grep exon /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2Only.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf

awk '!/CDS/ && !/start_codon/' /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnlyFiltered.gtf

/home/chughes/softwareTools/ggSashimi/ggsashimi.py -b /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/results/SRR7059715.sorted.bam -c chr11:83626956-83760000 -g /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnlyFiltered.gtf --min-coverage 10 --base-size 7 -o chla10Dlg2Upstream

#grep DLG2 /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf

#grep exon /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf

#awk '!/CDS/ && !/start_codon/' /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly2.gtf

#/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8616012.sorted.bam -c chr11:83626956-83760000 -g /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly2.gtf --min-coverage 50 --base-size 9
```

Ok that worked, so it was our GTF. We can filter this a bit more because it is showing many isoforms we don't want. First try this a different way.

```shell
grep DLG2 /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2Only.gtf

awk '!/CDS/ && !/start_codon/ && !/UTR/ && !/stop_codon/ && /exon/' /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2Only.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf

/home/chughes/softwareTools/ggSashimi/ggsashimi.py -b /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/results/SRR7059715.sorted.bam -c chr11:83626956-83760000 -g /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf --min-coverage 10 --base-size 9 -o chla10Dlg2Upstream

#grep DLG2 /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf

#awk '!/CDS/ && !/start_codon/ && !/UTR/ && !/stop_codon/ && /exon/' /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf

#/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8616012.sorted.bam -c chr11:83626956-83760000 -g /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf --min-coverage 50 --base-size 9
```

Ok that works and saves us a line of code. We can try for a specific transcript version.

```shell
grep 'ENST00000426717.6\|ENST00000398309.6' /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2Only.gtf

awk '!/start_codon/ && !/stop_codon/ && !/UTR/ && !/^CDS$/' /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2Only.gtf > /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf

/home/chughes/softwareTools/ggSashimi/ggsashimi.py -b /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/results/SRR7059715.sorted.bam -c chr11:83626956-83760000 -g /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf --min-coverage 10 --base-size 9 -o chla10Dlg2Upstream

#grep 'ENST00000426717.6\|ENST00000398309.6' /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf

#awk '!/start_codon/ && !/stop_codon/ && !/UTR/ && !/^CDS$/' /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2Only.gtf > /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf

#/projects/ptx_analysis/chughes/software/ggSashimi/sashimi-plot.py -b /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/starResults/SRR8616012.sorted.bam -c chr11:83626956-83760000 -g /projects/ptx_results/Sequencing/publishedStudies/ccleRnaSequencingData/sashimi/dlg2ExonsOnly.gtf --min-coverage 50 --base-size 9
```

Pass one at a time for the 3 cell lines.

```shell
#CHLA10
/home/chughes/softwareTools/ggSashimi/ggsashimi.py -b /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/results/SRR7059715.sorted.bam -c chr11:83626956-83760000 -g /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf --min-coverage 10 --base-size 9 -o chla10Dlg2Upstream

/home/chughes/softwareTools/ggSashimi/ggsashimi.py -b /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20201204_ewsFli1PrionBoulayPmid28844694/results/SRR5217668.sorted.bam -c chr11:83626956-83760000 -g /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf --min-coverage 10 --base-size 9 -o a673Dlg2Upstream

/home/chughes/softwareTools/ggSashimi/ggsashimi.py -b /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20201204_ewsFli1PrionBoulayPmid28844694/results/SRR5217670.sorted.bam -c chr11:83626956-83760000 -g /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/sashimiPlots/dlg2ExonsOnly.gtf --min-coverage 10 --base-size 9 -o sknmcDlg2Upstream
```
