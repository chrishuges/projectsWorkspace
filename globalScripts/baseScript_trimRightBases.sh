#! /bin/bash
bbdukLocation="/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"


##trim the last two bases from R2 reads
for i in *_R2_*.fastq.gz
do
  echo $i
  targetId=$(basename "$i" ".fastq.gz")
  eval $bbdukLocation in=$i out=${targetId}_trim.fastq.gz -ftr=100
done