#! /bin/bash

# this script will perform a Philosopher analysis of proteomics data using default settings.
# this script should be executed in the directory where your data is.

workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToSoftwareTesting/proteomics20210719_msfraggerTesting/"

# create the working directory, if necessary
if [ ! -d $workingDirectory ]; then
  printf "Data directory $workingDirectory does not exist, creating it.\n"
  eval "mkdir $workingDirectory"
else
  printf "Data directory $workingDirectory exists, moving on.\n"
fi

# prepare for philosopher analysis
eval philosopher workspace --init
eval philosopher database --reviewed --contam --id UP000005640

# example database name: 2021-07-19-decoys-reviewed-contam-UP000005640.fas

##############################
# convert the raw files to mzML
thermoRawFileParser="/home/chughes/softwareTools/thermoRawFileParser-1.3.4/ThermoRawFileParser.exe"

eval $thermoRawFileParser -i 



###############################
# run msfragger
msfragger="/home/chughes/softwareTools/msfragger-3.2/MSFragger-3.2/MSFragger-3.2.jar"

# example params call: java -jar /home/chughes/softwareTools/msfragger-3.2/MSFragger-3.2/MSFragger-3.2.jar --config closed
# change database name and other parameters of interest

java -Xmx32g -jar $msfragger closed_fragger.params ./controlSp3Temp70/ch_20210704_rapidTrypsinControlSp3_70C_1.mzML

#example: java -Xmx32g -jar /home/chughes/softwareTools/msfragger-3.2/MSFragger-3.2/MSFragger-3.2.jar closed_fragger.params ./controlSp3Temp70/ch_20210704_rapidTrypsinControlSp3_70C_1.mzML



###############################
# perform peptide and protein validation

philosopher peptideprophet --database 2021-07-19-decoys-reviewed-contam-UP000005640.fas --decoy rev_ --ppm --accmass --expectscore --decoyprobs --nonparam ./controlSp3Temp70/ch_20210704_rapidTrypsinControlSp3_70C_1.pepXML

philosopher proteinprophet ./controlSp3Temp70/interact-ch_20210704_rapidTrypsinControlSp3_70C_1.pep.xml

philosopher filter --sequential --razor --picked --tag rev_ --pepxml ./controlSp3Temp70/interact-ch_20210704_rapidTrypsinControlSp3_70C_1.pep.xml --protxml interact.prot.xml

philosopher report








###############################
# using pipeline


#################################
# script for raw file conversion
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToSoftwareTesting/proteomics20210719_msfraggerTesting/"
thermoRawFileParser="/home/chughes/softwareTools/thermoRawFileParser-1.3.4/ThermoRawFileParser.exe"
# find the raw files
fileList=()
  while IFS= read -d $'\0' -r file ; do
    fileList=("${fileList[@]}" "$file")
  done < <(find $workingDirectory -name "*.raw" -print0 | sort -z)
# convert them to mzML
for i in "${fileList[@]}"; do
  targetId=$(basename "$i" ".raw")
  if [ -f "$targetId.mzML" ]; then
    printf "mzML output for $targetId already exists, skipping file.\n\n"
  else
    eval mono $thermoRawFileParser -i $i
  fi
done
# send each to its own directory
for i in "${fileList[@]}"; do
  targetId=$(basename "$i" ".mzML")
  if [ -d $targetId ]; then
    printf "directory for $targetId already exists, skipping file.\n\n"
  else
    eval mkdir $targetId
    eval mv ${targetId}.raw ./$targetId/
    eval mv ${targetId}.mzML ./$targetId/
  fi
done
#################################



## /home/chughes/softwareTools/philosopher-3.4.13/philosopher.yml
eval philosopher workspace --init
eval philosopher database --reviewed --contam --id UP000005640
eval scp /home/chughes/softwareTools/philosopher-3.4.13/philosopher.yml ./

philosopher pipeline --config philosopher.yml CONTROL_1 CONTROL_2 CONTROL_3 HDAC5_1 HDAC5_2 HDAC5_3
