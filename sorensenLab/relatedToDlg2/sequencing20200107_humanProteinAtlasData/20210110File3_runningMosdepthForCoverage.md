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


