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
sampleTag="20210610_ahmedSilacTest"
#this is the desired location for the output of the data processing process
desiredBaseLocation="/projects/ptx_results/2021/Lumos-VCP/06-Jun/Samples/CH/"
#this is the base location where your raw data is stored
dataStorageLocation="/projects/ptx_results/2021/Lumos-VCP/06-Jun/Samples/CH/"
#this is text that will be appended to your output folder that is created for this analysis
folderAdapter="dataProcessing_"
#users do not need to edit the statement below
folderToCreate="$desiredBaseLocation$folderAdapter$sampleTag"
#starting directory
startingDirectory="${PWD}"


##############################################################
##############################################################
#locations of software tools
crapDatabase="/projects/ptx_analysis/chughes/databases/crap_jan2019.fasta"
rawtoolsLocation="/projects/ptx_analysis/chughes/software/RawTools/rawtools203/RawTools.exe"
javaLocation="/gsc/software/linux-x86_64-centos7/jdk-14.0.2/bin/java"
searchGuiLocation="/projects/ptx_analysis/chughes/software/SearchGUI-4.0.37/SearchGUI-4.0.37.jar"
peptideshakerLocation="/projects/ptx_analysis/chughes/software/PeptideShaker-2.0.31/PeptideShaker-2.0.31.jar"



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
####set up the database to use
uniprotVersion=$( date +%b%Y )
databaseExtension="TargetDecoy"

#check if the database exists, and if not, get it
if [ ! -f "./uniprotHumanCrap$databaseExtension$uniprotVersion.fasta" ]; then
    printf "Database not found!\n\n"
    eval "wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640_9606.fasta.gz"
    eval "gunzip UP000005640_9606.fasta.gz"
    eval "mv UP000005640_9606.fasta ./uniprotHuman$uniprotVersion.fasta"
    eval "cat ./uniprotHuman$uniprotVersion.fasta $crapDatabase > ./uniprotHumanCrap$uniprotVersion.fasta"
    eval "java -cp $searchGuiLocation eu.isas.searchgui.cmd.FastaCLI -in ./uniprotHumanCrap$uniprotVersion.fasta -decoy"    eval "mv ./uniprotHumanCrap${uniprotVersion}_concatenated_target_decoy.fasta ./uniprotHumanCrap$databaseExtension$uniprotVersion.fasta"
    eval "scp ./uniprotHuman$uniprotVersion.fasta /projects/ptx_analysis/chughes/databases/uniprot/"
    eval "scp ./uniprotHumanCrap$uniprotVersion.fasta /projects/ptx_analysis/chughes/databases/uniprot/"
    eval "scp ./uniprotHumanCrap$databaseExtension$uniprotVersion.fasta /projects/ptx_analysis/chughes/databases/uniprot/"
    printf "Database created and ready to use.\n\n"
else
  printf "\nDatabase file already exists: ./uniprotHumanCrap$databaseExtension$uniprotVersion.fasta\n\n"
fi

####process the database to make an annotation set with R
#just in case you are in a screen session, tell the screen session where R is
eval "export PATH=/gsc/software/linux-x86_64-centos7/R-4.0.2/bin:$PATH"

#check if the the annotation index exists, and if not, create it
if [ ! -f "uniprotHuman$uniprotVersion.fasta.annotated.rds" ]; then
  targetId="uniprotHuman$uniprotVersion.fasta"
  printf "Annotation index for $targetId does not exist, so it will be created. This can be slow but will only happen once per index.\n\n"
  rAnnotationIndex="Rscript /projects/ptx_analysis/chughes/projectsWorkspace/globalScripts/annotateUniprotDatabase.R ${PWD}/uniprotHuman$uniprotVersion.fasta"
  eval $rAnnotationIndex
  eval "scp ${PWD}/uniprotHuman$uniprotVersion.fasta.annotated.rds /projects/ptx_analysis/chughes/databases/uniprot/"
fi


##############################################################
##############################################################
####set up the parameters file
fixedModifications=''
variableModifications='"Oxidation of M"'

#####call the parameters file creation
callingParameters='java -cp "$searchGuiLocation" eu.isas.searchgui.cmd.IdentificationParametersCLI -out ./databaseSearchParameters.par -prec_ppm 1 -prec_tol 10 -frag_ppm 0 -frag_tol 0.05 -fixed_mods "$fixedModifications" -variable_mods "$variableModifications" -msgf_instrument 3 -msgf_protocol 5 -msgf_fragmentation 3'

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
    rawtoolsExecution="mono $rawtoolsLocation -f $i -q -uxmR -o $PWD"
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
    searchGuiExecution="$javaLocation -Xmx120g -cp $searchGuiLocation eu.isas.searchgui.cmd.SearchCLI -spectrum_files $i -fasta_file ${PWD}/uniprotHumanCrap"$databaseExtension""$uniprotVersion".fasta -output_folder $PWD -id_params ./databaseSearchParameters.par -output_option 1 -xtandem 1 -msgf 1 -comet 1"
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

##get the list of files to process
fileList=()
  while IFS= read -d $'\0' -r file ; do
    fileList=("${fileList[@]}" "$file")
  done < <(find $PWD -name "*.zip" -print0 | sort -z)

###now execute on each of these
for i in "${fileList[@]}"; do
  targetId=$(basename "$i" ".rawsearchgui_out.zip")
  printf "\n $targetId"
  while ! [ -f "${targetId}_Default_PSM_Report.txt" ]; do
    printf "\nNo peptide shaker reports output exists, so it will be created.\n"
    #peptide shaker execution
    peptideshakerExecution="$javaLocation -Xmx120g -cp $peptideshakerLocation eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference $targetId -identification_files $i -spectrum_files $PWD -fasta_file ${PWD}/uniprotHumanCrap"$databaseExtension""$uniprotVersion".fasta -id_params ./databaseSearchParameters.par -out $PWD/$targetId.psdb"
    eval $peptideshakerExecution
    #reports output
    reportsExecution="$javaLocation -Xmx120g -cp $peptideshakerLocation eu.isas.peptideshaker.cmd.ReportCLI -in $targetId.psdb -out_reports $PWD -reports 0,3"
    eval $reportsExecution
  done
done
printf "\nFinished all exporting reports for $sampleTag data set.\n\n"

####final output
eval "rm ${PWD}/*.html"
eval "mv ${startingDirectory}/dataProcessing.txt $PWD"
eval "scp ${startingDirectory}/rawDataProcessing_m16tmt11CidItHumanTrypsinCombined.sh $PWD"
printf "\n\nFinished processing a total of ${#fileList[@]} files!\n\n"