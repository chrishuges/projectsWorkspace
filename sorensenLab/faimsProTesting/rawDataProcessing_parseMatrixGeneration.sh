#! /bin/bash

##############################################################
#this first part calls the help if the user triggers it
usage(){
printf "
This script will process proteomic data to generate PeptideShaker Reports.\n
The sample tag and desired base location parameters should be edited before running the script.\n
You sould execute this script in the shell scripting directory.\n
All results will be written to the folder created based on the sample tag and base location.\n\n"
}

if [[ $1 == -h ]] || [[ $1 == --help ]];
then
  usage
  exit
fi


##############################################################
##############################################################
#users must edit the below variables
#this is the text that will identify your files
sampleTag="20201016_HeLaTMT_150um-19um-15cm_1in10_400ng_CV"
#this is the desired location for the output of the data processing process
desiredBaseLocation="/projects/ptx_results/2020/Fusion/10-Oct/Samples/CH/"
#this is the base location where your raw data is stored
dataStorageLocation="/projects/ptx_results/2020/Fusion/10-Oct/Samples/CH/"
#this is text that will be appended to your output folder that is created for this analysis
folderAdapter="dataProcessing_"
#users do not need to edit the statement below
folderToCreate="$desiredBaseLocation$folderAdapter$sampleTag"
#starting directory
startingDirectory="${PWD}"


##############################################################
##############################################################
#locations of software tools
rawtoolsLocation="/projects/ptx_analysis/chughes/software/RawTools/rawtools202/RawTools.exe"


##############################################################
##############################################################
####first we need to make the location for the data storage
makingDirectory="mkdir $folderToCreate"

#check if directory does not exist and create it if it doesn't
if [ ! -d "$folderToCreate" ];
then
  printf "Creating the processing directory.\n"
  eval $makingDirectory
fi

#move into the processing directory
eval "cd $folderToCreate"
printf "\nData processing proceeding in the directory: $PWD\n\n"



##############################################################
##############################################################
####run rawtools to get MGF files
##first find the raw files you want to process
fileList=()
  while IFS= read -d $'\0' -r file ; do
    fileList=("${fileList[@]}" "$file")
  done < <(find $dataStorageLocation -name "*$sampleTag*.raw" -print0 | sort -z)

###now execute rawtools on each of these
for i in "${fileList[@]}"; do
  targetId=$(basename "$i" ".raw")
  if [ -f "$targetId.raw.mgf" ]; then
    printf "MGF output for $targetId already exists, skipping file.\n\n"
  else
    rawtoolsExecution="mono $rawtoolsLocation -f $i -puxmR --chro 1B -o $PWD"
    eval $rawtoolsExecution
  fi
done

eval "mv ${startingDirectory}/dataProcessing.txt $PWD"
printf "\n\nFinished processing a total of ${#fileList[@]} files!\n\n"