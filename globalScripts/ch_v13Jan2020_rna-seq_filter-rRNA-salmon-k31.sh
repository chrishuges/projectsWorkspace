#! /bin/bash

##############################################################
##############################################################
#this first part calls the help if the user triggers it
usage(){
printf "
This script will perform basic processing for RiboSeq data.\n
The working directory should only contain the fastq.gz files for the samples you wish to process.\n
Before running, open the shell script and edit the barcodes and working directory for your desired samples.\n
All results will be written to the directory specified in output directory.\n\n"
}

if [[ $1 == -h ]] || [[ $1 == --help ]];
then
  usage
  exit
fi


##############################################################
##############################################################
####input the barcodes

barcodes=('SRR2748141' 'SRR2748142')
working_directory="/projects/ptx_results/OtherDataSets/morrow_2018-NatureMedicine_PMID29334376/rna-seq_mg633"


##############################################################
##############################################################
#locations of software tools

bbduk_location="/projects/ptx_analysis/chughes/software/BBTools/bbmap/bbduk.sh"
salmon_location="/projects/ptx_analysis/chughes/software/salmon-latest_linux_x86_64/bin/salmon"
rRNA_kmers_location="/projects/ptx_analysis/chughes/databases/genbank_rRNA/genbank_aug2019_rRNA.fasta"
fastqc_location="/gsc/software/linux-x86_64-centos7/fastqc-0.11.5/fastqc"
desired_index="/projects/ptx_analysis/chughes/databases/salmon_indexes/salmon-index_human_jan2020/index_k31"


##############################################################
##############################################################
####perform adapter trimming and read filtering

eval "cd $working_directory"
printf "Performing read filtering.\n"
printf "Output will be written to: $PWD\n\n"

for i in "${barcodes[@]}"
do
  #####get the list of files
  file_list=()
  while IFS= read -d $'\0' -r file ; do
    file_list=("${file_list[@]}" "$file")
  done < <(find $PWD -name "*$i*.fastq.gz" -print0 | sort -z)

  #print the input files
  printf "Total number of files in the array: ${#file_list[@]}\n"
  printf "First input file: ${file_list[0]}\n"
  printf "Second input file: ${file_list[1]}\n\n"

  #check if the output file exists, and if not, create it
  if [ -f "${file_list[0]::-9}_filtered.fastq.gz" ]; then
    printf "The trimmed output files for barcode $i already exist. Skipping.\n\n"
  else
    #make the output file names
    outfile_first="${file_list[0]::-9}_filtered.fastq.gz"
    outfile_second="${file_list[1]::-9}_filtered.fastq.gz"
    #call bbduk
    bbduk_call="$bbduk_location in1=${file_list[0]} in2=${file_list[1]} out1=$outfile_first out2=$outfile_second ordered=t"
    #call the command
    eval $bbduk_call
    printf "Finished filtering files associated with barcode: $i.\n\n"
  fi
done



##############################################################
##############################################################
####perform rRNA depletion using the GenBank library

printf "Performing splitting of rRNA reads.\n"
printf "Output will be written to: $PWD\n\n"
for i in "${barcodes[@]}"
do
  #####get the list of files
  file_list=()
  while IFS= read -d $'\0' -r file ; do
    file_list=("${file_list[@]}" "$file")
  done < <(find $PWD -name "*$i*_filtered.fastq.gz" -print0 | sort -z)

  #print the input files
  printf "Total number of files in the array: ${#file_list[@]}\n"
  printf "First input file: ${file_list[0]}\n"
  printf "Second input file: ${file_list[1]}\n\n"

  #check if the output file exists, and if not, create it
  if [ -f "${file_list[0]::-18}_mRNA.fastq.gz" ]; then
    printf "The rRNA depleted output files for barcode $i already exist. Skipping.\n\n"
  else
    #make the output file names
    outfile_mRNA_first="${file_list[0]::-18}_mRNA.fastq.gz"
    outfile_mRNA_second="${file_list[1]::-18}_mRNA.fastq.gz"
    #call bbduk
    bbduk_call="$bbduk_location in1=${file_list[0]} in2=${file_list[1]} out1=$outfile_mRNA_first out2=$outfile_mRNA_second ref=$rRNA_kmers_location k=31 hdist=1 ordered=t stats=${i}_rRNA_stats.txt"
    #call the command
    eval $bbduk_call
    printf "Finished rRNA depletion in files associated with barcode: $i.\n\n"
  fi
done


##############################################################
##############################################################
####peform fastqc on the depleted file

printf "Performing fastqc analysis.\n"
printf "Output will be written to: $PWD\n\n"
for i in "${barcodes[@]}"
do
  #####get the list of files
  file_list=()
  while IFS= read -d $'\0' -r file ; do
    file_list=("${file_list[@]}" "$file")
  done < <(find $PWD -name "*$i*_mRNA.fastq.gz" -print0 | sort -z)

  #print the input files
  printf "Total number of files in the array: ${#file_list[@]}\n"
  printf "First input file: ${file_list[0]}\n"
  printf "Second input file: ${file_list[1]}\n\n"

  #check if the output file exists, and if not, create it
  if [ ! -f "${file_list[0]::-9}_fastqc.zip" ]; then
    #call fastqc
    fastqc_call="$fastqc_location -o $PWD -f fastq ${file_list[0]}"
    #call the command
    eval $fastqc_call
    printf "Finished fastqc analysis for first file associated with barcode: $i.\n\n"
  else
    printf "Already finished fastqc analysis for first file associated with barcode: $i.\n\n"
  fi
  ######
  if [ ! -f "${file_list[1]::-9}_fastqc.zip" ]; then
    #call fastqc
    fastqc_call="$fastqc_location -o $PWD -f fastq ${file_list[1]}"
    #call the command
    eval $fastqc_call
    printf "Finished fastqc analysis for second file associated with barcode: $i.\n\n"
  else
    printf "Already finished fastqc analysis for second file associated with barcode: $i.\n\n"
  fi
done



##############################################################
##############################################################
####perform Salmon analysis

printf "Performing Salmon analysis.\n"
printf "Output will be written to: $PWD\n\n"

#########
for i in "${barcodes[@]}"
do
  #####get the list of files
  file_list=()
  while IFS= read -d $'\0' -r file ; do
    file_list=("${file_list[@]}" "$file")
  done < <(find $PWD -name "*$i*_filtered.fastq.gz" -print0 | sort -z)

  #print the input files
  printf "Total number of files in the array: ${#file_list[@]}\n"
  printf "First input file: ${file_list[0]}\n"
  printf "Second input file: ${file_list[1]}\n\n"

  #check if the output file exists, and if not, create it
  if [ -d "${i}_quant" ]; then
    printf "The quantification output files for barcode $i already exist. Skipping.\n\n"
  else
    #I don't set the library type here, I have Salmon automatically detect it, but usually it is 'IU' for non-stranded paired-end reads
    #call salmon
    salmon_call="$salmon_location quant -i $desired_index -l A -1 ${file_list[0]} -2 ${file_list[1]} --validateMappings -o ${i}_quant --gcBias"
    #call the command
    eval $salmon_call
    printf "Finished Salmon analysis for files associated with barcode: $i.\n\n"
  fi
done


#####final cleanup
eval "scp /projects/ptx_analysis/chughes/shell-scripting/rna-seq/ch_v13Jan2020_rna-seq_filter-rRNA-salmon-k31.sh $PWD"
eval "mv /projects/ptx_analysis/chughes/shell-scripting/rna-seq/data-processing.txt $PWD"
#####delete the large files we no longer need
eval "rm $PWD/*mRNA*.gz"
eval "rm $PWD/*filtered*.gz"
#final print
printf "Processing finished for ${#barcodes[@]} input barcodes."
