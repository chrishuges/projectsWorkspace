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
  if [ -f "${targetId}.mzML" ]; then
    printf "mzML output for $targetId already exists, skipping file.\n\n"
  else
    eval mono $thermoRawFileParser -i $i
  fi
done
# send each to its own directory
for i in "${fileList[@]}"; do
  targetId=$(basename "$i" ".raw")
  if [ -d ${targetId:12} ]; then
    printf "directory for $targetId already exists, skipping file.\n\n"
  else
    eval mkdir ${targetId:12}
    eval mv ${targetId}.raw ./${targetId:12}/
    eval mv ${targetId}.mzML ./${targetId:12}/
  fi
done
#################################


#################################
# script for raw file conversion that just converts files and does RawTools for TMT
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/proteomics20210930_ewsDtagPtxPmid32948771/"
thermoRawFileParser="/home/chughes/softwareTools/thermoRawFileParser-1.3.4/ThermoRawFileParser.exe"
rawTools="/home/chughes/softwareTools/rawtools2.0.3a/RawTools.exe"
# find the raw files
fileList=()
  while IFS= read -d $'\0' -r file ; do
    fileList=("${fileList[@]}" "$file")
  done < <(find $workingDirectory -name "*.raw" -print0 | sort -z)
# convert them to mzML
for i in "${fileList[@]}"; do
  targetId=$(basename "$i" ".raw")
  if [ -f "${targetId}.mzML" ]; then
    printf "mzML output for $targetId already exists, skipping file.\n\n"
  else
    eval mono $thermoRawFileParser -i $i
    eval mono $rawTools -f $i -q -r TMT11 -uxR -chro 1B 
  fi
done
#################################



## /home/chughes/softwareTools/philosopher-3.4.13/philosopher.yml
philosopher workspace --init
philosopher database --reviewed --contam --id UP000005640
#scp /home/chughes/softwareTools/philosopher-3.4.13/philosopher.yml ./
#scp /home/chughes/softwareTools/philosopher-4.0.0/philosopher.yml ./

/home/chughes/softwareTools/msfragger-3.2/MSFragger-3.2/MSFragger-3.2.jar


philosopher pipeline --config philosopher.yml rapidTrypsinControlSp3_40C_1 rapidTrypsinControlSp3_55C_1 rapidTrypsinControlSp3_70C_1 rapidTrypsinControl_40C_1 rapidTrypsinControl_55C_1 rapidTrypsinControl_70C_1 rapidTrypsinDeoxycholateSp3_40C_1 rapidTrypsinDeoxycholateSp3_55C_1 rapidTrypsinDeoxycholateSp3_70C_1 rapidTrypsinDeoxycholate_40C_1 rapidTrypsinDeoxycholate_55C_1 rapidTrypsinDeoxycholate_70C_1 rapidTrypsinProteaseMaxSp3_40C_1 rapidTrypsinProteaseMaxSp3_55C_1 rapidTrypsinProteaseMaxSp3_70C_1 rapidTrypsinProteaseMax_40C_1 rapidTrypsinProteaseMax_55C_1 rapidTrypsinProteaseMax_70C_1


philosopher pipeline --config philosopher.yml

