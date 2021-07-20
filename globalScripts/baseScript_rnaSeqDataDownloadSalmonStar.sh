#! /bin/bash
bbdukLocation="/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
salmonLocation="/home/chughes/softwareTools/salmon-1.5.1/bin/salmon"
indexLocation="/home/chughes/databases/refgenieManualGenomes/hg38/salmonPartialSaIndex_072021/default/"
starLocation="/home/chughes/softwareTools/STAR-2.7.9a/bin/Linux_x86_64/STAR"
referenceLocation="/home/chughes/databases/refgenieManualGenomes/hg38/starIndex_072021/default/"
annotationLocation="/home/chughes/databases/refgenieManualGenomes/hg38/gencodeGtf_072021/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.gtf"
samtoolsLocation="/home/chughes/softwareTools/samtools-1.12/samtools"
rawDataOutputDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210713_a673DoxTreatmentPmid34009296/"

#download the data files
if [ ! -d $rawDataOutputDirectory ]; then
  printf "Data directory $rawDataOutputDirectory does not exist, creating it.\n"
  eval "mkdir $rawDataOutputDirectory"
else
  printf "Data directory $rawDataOutputDirectory exists, moving on.\n"
fi

eval "cd ${rawDataOutputDirectory}"
eval "mkdir ${rawDataOutputDirectory}starResults"
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/007/SRR8561337/SRR8561337_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/007/SRR8561337/SRR8561337_2.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/004/SRR8561334/SRR8561334_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/004/SRR8561334/SRR8561334_2.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/003/SRR8561333/SRR8561333_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/003/SRR8561333/SRR8561333_2.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/005/SRR8561335/SRR8561335_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/005/SRR8561335/SRR8561335_2.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/007/SRR8561347/SRR8561347_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/007/SRR8561347/SRR8561347_2.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/002/SRR8561332/SRR8561332_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/002/SRR8561332/SRR8561332_2.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/086/SRR11908086/SRR11908086_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/086/SRR11908086/SRR11908086_2.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/085/SRR11908085/SRR11908085_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/085/SRR11908085/SRR11908085_2.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/084/SRR11908084/SRR11908084_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/084/SRR11908084/SRR11908084_2.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/082/SRR11908082/SRR11908082_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/082/SRR11908082/SRR11908082_2.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/083/SRR11908083/SRR11908083_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/083/SRR11908083/SRR11908083_2.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/081/SRR11908081/SRR11908081_1.fastq.gz
eval wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/081/SRR11908081/SRR11908081_2.fastq.gz


#process the downloaded files
for i in SRR8561337 SRR8561334 SRR8561333 SRR8561335 SRR8561347 SRR8561332 SRR11908086 SRR11908085 SRR11908084 SRR11908082 SRR11908083 SRR11908081
do
  echo $i
  ##
  bbdukCall="$bbdukLocation in1=${rawDataOutputDirectory}${i}_1.fastq.gz in2=${rawDataOutputDirectory}${i}_2.fastq.gz ref=adapters out1=${rawDataOutputDirectory}${i}_1.clean.fastq.gz out2=${rawDataOutputDirectory}${i}_2.clean.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo"
  eval $bbdukCall
  
  ##
  salmonCall="$salmonLocation quant -i $indexLocation -l A -1 ${rawDataOutputDirectory}${i}_1.clean.fastq.gz -2 ${rawDataOutputDirectory}${i}_2.clean.fastq.gz --gcBias --validateMappings -o ${rawDataOutputDirectory}${i}_quant"
  eval $salmonCall
  
  ##
  starCall="$starLocation --runThreadN 12 --genomeDir ${referenceLocation} --readFilesIn ${rawDataOutputDirectory}${i}_1.clean.fastq.gz ${rawDataOutputDirectory}${i}_2.clean.fastq.gz --readFilesCommand zcat --sjdbGTFfile ${annotationLocation} --outFileNamePrefix ${rawDataOutputDirectory}starResults/${i}_ --outSAMtype BAM Unsorted"
  eval $starCall
  
  ##
  bamSortCall="$samtoolsLocation sort ${rawDataOutputDirectory}/starResults/${i}_Aligned.out.bam -o ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamSortCall
  
  ##
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamIndexCall
  
  ##
  coverageCall="bamCoverage -b ${rawDataOutputDirectory}/starResults/${i}.sorted.bam -o ${rawDataOutputDirectory}/starResults/${i}.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --centerReads -p 6"
  eval $coverageCall
done





