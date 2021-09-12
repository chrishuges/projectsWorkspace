#! /bin/bash
salmonLocation="/home/chughes/softwareTools/salmon-1.5.1/bin/salmon"
indexLocation="/home/chughes/databases/refgenieManualGenomes/hg38/salmonPartialSaIndex_072021/default/"
starLocation="/home/chughes/softwareTools/STAR-2.7.1a/bin/Linux_x86_64/STAR"
referenceLocation="/home/chughes/databases/refgenieManualGenomes/hg38/starIndex_072021/default/"
annotationLocation="/home/chughes/databases/refgenieManualGenomes/hg38/gencodeGtf_072021/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.gtf"
samtoolsLocation="/home/chughes/softwareTools/samtools-1.12/samtools"
rawDataOutputDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencingEwsCohortDbgapPmid25010205/"


#create the directory to store the files
if [ ! -d $rawDataOutputDirectory ]; then
  printf "Data directory $rawDataOutputDirectory does not exist, creating it.\n"
  eval "mkdir $rawDataOutputDirectory"
else
  printf "Data directory $rawDataOutputDirectory exists, moving on.\n"
fi

eval "cd ${rawDataOutputDirectory}"
eval "mkdir ${rawDataOutputDirectory}starResults"

#download the process the files
for i in SRR5163{671..757}
do
  echo $i
  
  ##prefetch the files from SRA
  eval prefetch --ngc prj_29838.ngc -p ${i}
  eval mv ${i}/${i}*.sra ./${i}.sra
  eval rm -r ${i}
  eval fasterq-dump --ngc prj_29838.ngc -p ${i}.sra
  eval rm ${i}.sra
  eval rm -r fasterq.tmp*
  eval gzip -v ${i}_1.fastq
  eval gzip -v ${i}_2.fastq

  ##salmon analylsis
  salmonCall="$salmonLocation quant -i $indexLocation -l A -1 ${rawDataOutputDirectory}${i}_1.fastq.gz -2 ${rawDataOutputDirectory}${i}_2.fastq.gz --gcBias --validateMappings -o ${rawDataOutputDirectory}${i}_quant"
  eval $salmonCall
  
  ##star alignment
  starCall="$starLocation --runThreadN 12 --genomeDir ${referenceLocation} --readFilesIn ${rawDataOutputDirectory}${i}_1.fastq.gz ${rawDataOutputDirectory}${i}_2.fastq.gz --readFilesCommand zcat --sjdbGTFfile ${annotationLocation} --outFileNamePrefix ${rawDataOutputDirectory}starResults/${i}_ --outSAMtype BAM Unsorted"
  eval $starCall
  
  ##bamsort
  bamSortCall="$samtoolsLocation sort ${rawDataOutputDirectory}/starResults/${i}_Aligned.out.bam -o ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamSortCall
  
  ##index bam files
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamIndexCall
  
  ##calculate coverage on chromosome 11
  coverageCall="bamCoverage -b ${rawDataOutputDirectory}/starResults/${i}.sorted.bam -o ${rawDataOutputDirectory}/starResults/${i}.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --centerReads -p 6"
  eval $coverageCall
done








