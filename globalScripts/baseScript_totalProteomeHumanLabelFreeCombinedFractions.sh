#! /bin/bash

##############################################################
#this first part calls the help if the user triggers it
usage(){
printf "
This script will process proteomic data using RawTools, SearchGUI, and PeptideShaker to generate result reports.\n
The sample tag and desired base location parameters should be edited before running the script.\n
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
sampleTag="20210625_a673ef1_trypsinLyscCombined"
#this is the desired location for the output of the data processing process
desiredBaseLocation="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/proteomics20210625_a673Ef1TotalProteome/"
#this is the base location where your raw data is stored
dataStorageLocation="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/proteomics20210625_a673Ef1TotalProteome/"
#this is text that will be appended to your output folder that is created for this analysis
folderAdapter="dataProcessing_"
#users do not need to edit the statement below
folderToCreate="$desiredBaseLocation$folderAdapter$sampleTag"
#starting directory
startingDirectory="${PWD}"


##############################################################
##############################################################
#locations of software tools
crapDatabase="/home/chughes/databases/proteomes/crap_jan2021.fasta"
rawtoolsLocation="/home/chughes/softwareTools/rawtools2.0.3a/RawTools.exe"
javaLocation="/usr/bin/java"
searchGuiLocation="/home/chughes/softwareTools/SearchGUI-4.0.41/SearchGUI-4.0.41.jar"
peptideshakerLocation="/home/chughes/softwareTools/PeptideShaker-2.0.33/PeptideShaker-2.0.33.jar"



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



##############################################################
##############################################################
####set up the database to use
uniprotVersion=$( date +%b%Y )
databaseExtension="TargetDecoy"

#move into the database directory
eval "cd /home/chughes/databases/proteomes/"
printf "\nDatabase processing proceeding in the directory: $PWD\n\n"

#check if the database exists, and if not, get it
if [ ! -f "/home/chughes/databases/proteomes/uniprotHumanCrap$databaseExtension$uniprotVersion.fasta" ]; then
  printf "Database not found!\n\n"
  eval "wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
  eval "gunzip UP000005640_9606.fasta.gz"
  eval "mv UP000005640_9606.fasta ./uniprotHuman$uniprotVersion.fasta"
  eval "cat ./uniprotHuman$uniprotVersion.fasta $crapDatabase > ./uniprotHumanCrap$uniprotVersion.fasta"
  eval "java -cp $searchGuiLocation eu.isas.searchgui.cmd.FastaCLI -in ./uniprotHumanCrap$uniprotVersion.fasta -decoy"    
  eval "mv ./uniprotHumanCrap${uniprotVersion}_concatenated_target_decoy.fasta ./uniprotHumanCrap$databaseExtension$uniprotVersion.fasta"
  eval "scp ./uniprotHuman$uniprotVersion.fasta $folderToCreate"
  eval "scp ./uniprotHumanCrap$uniprotVersion.fasta $folderToCreate"
  eval "scp ./uniprotHumanCrap$databaseExtension$uniprotVersion.fasta $folderToCreate"
  printf "Database created and ready to use.\n\n"
else
  printf "\nDatabase file already exists: /home/chughes/databases/proteomes/uniprotHumanCrap$databaseExtension$uniprotVersion.fasta\n\n"
  eval "scp /home/chughes/databases/proteomes/uniprotHumanCrap$databaseExtension$uniprotVersion.fasta $folderToCreate"
fi

####process the database to make an annotation set with R
#just in case you are in a screen session, tell the screen session where R is
eval "export PATH=/opt/R/4.1.0/bin:$PATH"

#check if the the annotation index exists, and if not, create it
if [ ! -f "/home/chughes/databases/proteomes/uniprotHuman$uniprotVersion.fasta.annotated.rds" ]; then
  targetId="uniprotHuman$uniprotVersion.fasta"
  printf "Annotation index for $targetId does not exist, so it will be created. This can be slow but will only happen once per index.\n\n"
  rAnnotationIndex="Rscript /mnt/Data/chughes/projectsWorkspace/globalScripts/proteomics/annotateUniprotDatabase.R /home/chughes/databases/proteomes/uniprotHuman$uniprotVersion.fasta"
  eval $rAnnotationIndex
  eval "scp ${PWD}/uniprotHuman$uniprotVersion.fasta.annotated.rds /home/chughes/databases/proteomes/"
else
  printf "\nAnnotation index already exists: /home/chughes/databases/proteomes/uniprotHuman$uniprotVersion.fasta.annotated.rds\n\n"
  eval "scp /home/chughes/databases/proteomes/uniprotHuman$uniprotVersion.fasta.annotated.rds $folderToCreate"
fi


##############################################################
##############################################################
#move into the working directory
eval "cd $folderToCreate"
printf "\nSearch processing proceeding in the directory: $PWD\n\n"

####set up the parameters file
fixedModifications='"Carbamidomethylation of C"'
variableModifications='"Oxidation of M, Acetylation of protein N-term"'
enzyme='Trypsin'
missedCleavages='6'
precTol='1' #in Daltons
fragTol='0.5' #in Daltons
msgfInst='0' #MS-GF+ instrument id option, 0: Low-res LCQ/LTQ, 1: Orbitrap/FTICR, 2: TOF, 3: Q-Exactive (Default).
msgfProt='0' #MS-GF+ protocol id option, 0: Automatic (Default, true), 1: Phosphorylation, 2: iTRAQ, 3: iTRAQPhospho, 4: TMT, 5: Standard.
msgfFrag='3' #MS-GF+ fragmentation id option, 0: Automatic: as written in the spectrum or CID if no info, 1: CID, 2: ETD, 3: HCD (Default). 4: UVPD.

#####call the parameters file creation
callingParameters='java -cp "$searchGuiLocation" eu.isas.searchgui.cmd.IdentificationParametersCLI -out ./databaseSearchParameters.par -prec_ppm 0 -prec_tol "$precTol" -frag_ppm 0 -frag_tol "$fragTol" -fixed_mods "$fixedModifications" -variable_mods "$variableModifications" -enzyme "$enzyme" -mc "$missedCleavages" -msgf_instrument "$msgfInst" -msgf_protocol "$msgfProt" -msgf_fragmentation "$msgfFrag"'

####call the command
if [ ! -f "./databaseSearchParameters.par" ]; then
  printf "\nParameters file does not exist. Creating it.\n\n"
  eval "$callingParameters"
else
  printf "Parameters file exists: databaseSearchParameters.par.\n\n"
fi



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
    rawtoolsExecution="mono $rawtoolsLocation -f $i -p -uxmR -o $PWD"
    eval $rawtoolsExecution
  fi
done


##############################################################
##############################################################
####run the searchGUI analysis
fileList=()
  while IFS= read -d $'\0' -r file ; do
    fileList=("${fileList[@]}" "$file")
  done < <(find $PWD -name "*.mgf" -print0 | sort -z)

###now execute on each of these
for i in "${fileList[@]}"; do
  if [ -f "${i::-4}_searchgui_out.zip" ]; then
    targetId=$(basename "$i" ".zip")
    printf "Search output for $targetId already exists, skipping file.\n\n"
  else
    searchGuiExecution="$javaLocation -Xmx30g -cp $searchGuiLocation eu.isas.searchgui.cmd.SearchCLI -spectrum_files $i -fasta_file ${PWD}/uniprotHumanCrap"$databaseExtension""$uniprotVersion".fasta -output_folder $PWD -id_params ./databaseSearchParameters.par -output_option 1 -xtandem 1 -msgf 1 -comet 1"
    eval $searchGuiExecution
  fi
done


##############################################################
##############################################################
####run peptideshaker and write the reports
####first move the quant files
if [ ! -d "quantFiles" ]; then
  eval "mkdir quantFiles"
  eval "mv *.txt ./quantFiles"
fi

###now execute on each of these
while ! [ -f "${sampleTag}_Default_PSM_Report.txt" ]; do
  printf "\nNo peptide shaker reports output exists, so it will be created.\n"
  #peptide shaker execution
  peptideshakerExecution="$javaLocation -Xmx30g -cp $peptideshakerLocation eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference $sampleTag -identification_files $PWD -spectrum_files $PWD -fasta_file ${PWD}/uniprotHumanCrap"$databaseExtension""$uniprotVersion".fasta -id_params ./databaseSearchParameters.par -out $PWD/$sampleTag.out.psdb"
  eval $peptideshakerExecution
  #reports output
  reportsExecution="$javaLocation -Xmx30g -cp $peptideshakerLocation eu.isas.peptideshaker.cmd.ReportCLI -in $sampleTag.out.psdb -out_reports $PWD -reports 0,3"
  eval $reportsExecution
done
printf "\nFinished all exporting reports for $sampleTag data set.\n\n"


####final output
eval "rm ${PWD}/*.html"
eval "mv ${startingDirectory}/dataProcessing.txt $PWD"
printf "\n\nFinished processing a total of ${#fileList[@]} files!\n\n"