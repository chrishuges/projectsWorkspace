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
sampleTag="10Oct2016_ADROIT_PrepCompare_TMT10_hph"
#this is the desired location for the output of the data processing process
desiredBaseLocation="/projects/ptx_results/2016/FusionData/10-Oct/Samples/CH/"
#this is the base location where your raw data is stored
dataStorageLocation="/projects/ptx_results/2016/FusionData/10-Oct/Samples/CH/"
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
rawtoolsLocation="/projects/ptx_analysis/chughes/software/RawTools/rawtools202/RawTools.exe"
searchGuiLocation="/projects/ptx_analysis/chughes/software/SearchGUI-4.0.0/SearchGUI-4.0.0.jar"
peptideshakerLocation="/projects/ptx_analysis/chughes/software/PeptideShaker-2.0.0/PeptideShaker-2.0.0.jar"



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

HpH_9.rawsearchgui_out.zip

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
    eval "java -cp $searchGuiLocation eu.isas.searchgui.cmd.FastaCLI -in ./uniprotHumanCrap$uniprotVersion.fasta -decoy"
    eval "mv ./uniprotHumanCrap${uniprotVersion}_concatenated_target_decoy.fasta ./uniprotHumanCrap$databaseExtension$uniprotVersion.fasta"
    eval "scp ./uniprotHuman$uniprotVersion.fasta /projects/ptx_analysis/chughes/databases/uniprot/"
    eval "scp ./uniprotHumanCrap$uniprotVersion.fasta /projects/ptx_analysis/chughes/databases/uniprot/"
    eval "scp ./uniprotHumanCrap$databaseExtension$uniprotVersion.fasta /projects/ptx_analysis/chughes/databases/uniprot/"
    printf "Database created and ready to use.\n\n"
else
  printf "\nDatabase file already exists: ./uniprotHumanCrap$databaseExtension$uniprotVersion.fasta\n\n"
fi

####process the database to make an annotation set with R
#just in case you are in a screen session, tell the screen session where R is
eval "export PATH=/gsc/software/linux-x86_64-centos7/R-3.6.3/bin:$PATH"

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
fixedModifications='"TMTpro of peptide N-term"'
variableModifications='"Oxidation of M, TMTpro of K"'

#####call the parameters file creation
callingParameters='java -cp "$searchGuiLocation" eu.isas.searchgui.cmd.IdentificationParametersCLI -out ./databaseSearchParameters.par -prec_tol 20 -frag_tol 0.5 -fixed_mods "$fixedModifications" -variable_mods "$variableModifications" -msgf_instrument 0 -msgf_protocol 4 -msgf_fragmentation 1 -meta_morpheus_dissociation_type CID'

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
    rawtoolsExecution="mono $rawtoolsLocation -f $i -q -r TMT11 -uxmR -o $PWD"
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
    searchGuiExecution="java -Xmx120g -cp $searchGuiLocation eu.isas.searchgui.cmd.SearchCLI -spectrum_files $i -fasta_file ./uniprotHumanCrap"$databaseExtension""$uniprotVersion".fasta -output_folder $PWD -id_params ./databaseSearchParameters.par -output_option 1 -xtandem 1 -msgf 1 -comet 1 -meta_morpheus 1"
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
while ! [ -f "n_${sampleTag}_1_Default_PSM_Report.txt" ]; do
  printf "\nNo peptide shaker reports output exists, so it will be created.\n"
  #peptide shaker execution
  peptideshakerExecution="java -Xmx120g -cp $peptideshakerLocation eu.isas.peptideshaker.cmd.PeptideShakerCLI -identification_files $PWD -spectrum_files $PWD -fasta_file ./uniprotHumanCrap"$databaseExtension""$uniprotVersion".fasta -id_params ./databaseSearchParameters.par -out $PWD/$sampleTag.out.cpsx"
  eval $peptideshakerExecution
  #reports output
  reportsExecution="java -Xmx120g -cp $peptideshakerLocation eu.isas.peptideshaker.cmd.ReportCLI -in $sampleTag.out.cpsx -out_reports $PWD -reports 0,3 -documentation 0,3"
  eval $reportsExecution
done
printf "\nFinished all exporting reports for $sampleTag data set.\n\n"


####final output
eval "rm *.html"
eval "mv ${startingDirectory}/dataProcessing.txt $PWD"
printf "\n\nFinished processing a total of ${#fileList[@]} files!\n\n"