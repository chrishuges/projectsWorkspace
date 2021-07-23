#! /bin/bash

#################################
# script for raw file conversion
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/proteomics20210719_ewsCellTotalProteomes/"
thermoRawFileParser="/home/chughes/softwareTools/thermoRawFileParser-1.3.4/ThermoRawFileParser.exe"

#################################
# set the experiment identifiers
experimentList=(chla10IsocloneA sknmcIsocloneA a673ef1)

for j in "${experimentList[@]}"; do
  # find the raw files
  fileList=()
  while IFS= read -d $'\0' -r file ; do
    fileList=("${fileList[@]}" "$file")
  done < <(find $workingDirectory -name "*${j}*.raw" -print0 | sort -z)
  # convert them to mzML
  for i in "${fileList[@]}"; do
    targetId=$(basename "$i" ".raw")
    if [ -f "${targetId}.mzML" ]; then
      printf "mzML output for $targetId already exists, skipping file.\n\n"
    else
      eval mono $thermoRawFileParser -i $i
    fi
  done
  # send files for each experiment to a directory
  if [ -d ${j} ]; then
    printf "directory for $j already exists, skipping file.\n\n"
  else
    eval mkdir ${j}
    eval "mv *${j}*.raw ./${j}"
    eval "mv *${j}*.mzML ./${j}"
  fi
done
