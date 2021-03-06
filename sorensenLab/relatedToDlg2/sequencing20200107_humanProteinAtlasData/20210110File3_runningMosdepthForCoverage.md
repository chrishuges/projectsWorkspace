## Bedtools coverage analysis

This document describes running [Mosdepth](https://github.com/brentp/mosdepth) to get coverage data for DLG2. 

### Data processing

I already have aligned bam files, so I just need to run mosdepth. I am going to do this on a subset of the alignment data (specifically, chromosome 11) because this is where DLG2 is located it it will speed things up. For this I needed to know the 'chr' annotation and I looked at this in the bam files using the command `/projects/ptx_analysis/chughes/software/samtools-1.9/samtools view -h ./ERR315432.sorted.bam`.

Mosdepth will give coverage for specific regions you pass it, along with per base counts. Note, that for the per-base counts it will combine adjacent counts when they are the same, so sometimes it you won't have a value for every base that you passed it originally. To make the regions file, I used the GTF from the alignment to extract all possible transcripts and regions associated with DLG2 and wrote a bed file called 'targetTranscriptList.bed'. 



```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/humanProteinAtlasRnaSequencingData/starResults/"
samtoolsLocation="/projects/ptx_analysis/chughes/software/samtools-1.9/samtools"
mosdepthLocation="/projects/ptx_analysis/chughes/software/mosdepth-0.3.1/mosdepth"
targetGeneGtf="/projects/ptx_results/Sequencing/publishedStudies/humanProteinAtlasRnaSequencingData/starResults/targetTranscriptList.gtf"
#########################################
for i in ERR315455 ERR315477 ERR315432
do
  echo $i
  chr11ExtractCall="$samtoolsLocation view -b ${rawDataOutputDirectory}${i}.sorted.bam chr11:82775907-86346052 > ${rawDataOutputDirectory}${i}.chr11.bam"
  eval $chr11ExtractCall
  chr11IndexCall="$samtoolsLocation index ${rawDataOutputDirectory}${i}.chr11.bam"
  eval $chr11IndexCall
  eval "cd ${rawDataOutputDirectory}/mosdepthAnalysis"
  eval grep DLG2 /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf > targetTranscriptList.gtf
  eval cat targetTranscriptList.gtf |  awk 'OFS="\t" {if ($3=="exon") {print $1,$4-1,$5,$10,$16,$7}}' | tr -d '";' > targetTranscriptList.bed
  mosdepthCall="$mosdepthLocation --by ${rawDataOutputDirectory}/mosdepthAnalysis/targetTranscriptList.bed ${i} ${rawDataOutputDirectory}${i}.chr11.bam"
  eval $mosdepthCall
done
```

If you wanted to make single bp regions and count that way (this will eliminate the binning from mosdepth), you could use bedtools with the command `/projects/ptx_analysis/chughes/software/bedtools2/bin/bedtools.static.binary makewindows -b targetTranscriptList.bed -w 1 -i src > targetTranscriptsSingleBase.bed`. 

I am adding to this after the fact. I want to create bigwig files that I can use for coverage analysis downstream. I can do this directly from the chr11 bam files that I created in the above code.

```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/humanProteinAtlasRnaSequencingData/starResults/"
targetGeneGtf="/projects/ptx_results/Sequencing/publishedStudies/humanProteinAtlasRnaSequencingData/starResults/targetTranscriptList.gtf"
#########################################
for i in ERR315455 ERR315477 ERR315432
do
  echo $i
  coverageCall="bamCoverage -b ${rawDataOutputDirectory}${i}.sorted.bam -o ${rawDataOutputDirectory}${i}.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --centerReads -p 6"
  eval $coverageCall
done
```

I use these files in the `20210113File6_creatingDlg2CoverageMapsRevisited.Rmd` analysis file.

I went back to this again. I wanted to use the refseq annotation instead. So, I downloaded it from [NCBI FTP](ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/) to `/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/refseq_gtf/default`. I processed it as I did above and ran the bamCoverage analysis in order to get the coverage data.


```shell
grep DLG2 /projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/refseq_gtf/default/GRCh38_latest_genomic_20210202.gtf > targetTranscriptList.gtf

##script
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/Sequencing/publishedStudies/humanProteinAtlasRnaSequencingData/starResults/"
targetGeneGtf="/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/refseq_gtf/default/targetTranscriptList.gtf"
#########################################
for i in ERR315455 ERR315477 ERR315432
do
  echo $i
  coverageCall="bamCoverage -b ${rawDataOutputDirectory}${i}.sorted.bam -o ${rawDataOutputDirectory}${i}_refseq.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --centerReads -p 6"
  eval $coverageCall
done
```

