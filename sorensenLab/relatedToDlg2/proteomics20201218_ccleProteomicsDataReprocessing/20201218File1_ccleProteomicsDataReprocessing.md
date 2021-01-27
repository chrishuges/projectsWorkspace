## Reprocessing CCLE proteomics data

This document details reprocessing of the raw data for the CCLE proteomics performed by the group of Gygi at Harvard.

Quantitative Proteomics of the Cancer Cell Line Encyclopedia
Nusinow et al., Cell, PMID: 31978347, MASSIVE: MSV000085836

### Data download the preparation

The data are hosted on the MASSIVE ftp and I just downloaded them from there into the directory `/projects/ptx_results/OtherDataSets/dataset20201217_ccleProteomicsPmid31978347`. I had to make some changes to some of the raw file names in order to make it a cohesive set in terms of naming convention, but these changes were minor. 

We already have some good scripts for processing the actual data, but there are a total of 38 TMT batches here with 12 fractions in each, so I want it to loop through all of these batches individually. 

Batch names:

```
Prot_01
Prot_02
Prot_03
Prot_04
Prot_05
Prot_06
Prot_07
Prot_08
Prot_09
Prot_10
Prot_11
Prot_12
Prot_13
Prot_14
Prot_15
Prot_16
Prot_17
Prot_18
Prot_19
Prot_20
Prot_21
Prot_22
Prot_23
Prot_24
Prot_25
Prot_26
Prot_27
Prot_28
Prot_29
Prot_30
Prot_31
Prot_32
Prot_33
Prot_34
Prot_35
Prot_36
Prot_37
Prot_38
Prot_39
Prot_40
Prot_41
Prot_42
Prot_01_R2
Prot_01_R3
Prot_02_R2
Prot_02_R3
```

I guess one way to do this is to just loop over these accessions with our normal script. This is probably the easiest thing to do. 

```shell
#!/bin/bash
rawDataDirectory="/projects/ptx_results/OtherDataSets/dataset20201217_ccleProteomicsPmid31978347/"
############################################
for j in Prot_{01..42}  
do
  echo $j
  ##############################################################
  ##############################################################
  #users must edit the below variables
  #this is the text that will identify your files
  sampleTag=$j
  #this is the desired location for the output of the data processing process
  desiredBaseLocation=$rawDataDirectory
  #this is the base location where your raw data is stored
  dataStorageLocation=$rawDataDirectory
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
  searchGuiLocation="/projects/ptx_analysis/chughes/software/SearchGUI-4.0.4/SearchGUI-4.0.4.jar"
  peptideshakerLocation="/projects/ptx_analysis/chughes/software/PeptideShaker-2.0.5/PeptideShaker-2.0.5.jar"
  

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
  fixedModifications='"Carbamidomethylation of C, TMT 11-plex of peptide N-term, TMT 11-plex of K"'
  variableModifications='"Oxidation of M"'

  #####call the parameters file creation
  callingParameters='java -cp "$searchGuiLocation" eu.isas.searchgui.cmd.IdentificationParametersCLI -out ./databaseSearchParameters.par -prec_ppm 0 -prec_tol 1 -frag_ppm 0 -frag_tol 0.5 -fixed_mods "$fixedModifications" -variable_mods "$variableModifications" -msgf_instrument 0 -msgf_protocol 4 -msgf_fragmentation 1'

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
    if [ -f "${i::-4}searchgui_out.zip" ]; then
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

  ###now execute on each of these
  while ! [ -f "${sampleTag}_Default_PSM_Report.txt" ]; do
    printf "\nNo peptide shaker reports output exists, so it will be created.\n"
    #peptide shaker execution
    peptideshakerExecution="$javaLocation -Xmx120g -cp $peptideshakerLocation eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference $sampleTag -identification_files $PWD -spectrum_files $PWD -fasta_file ${PWD}/uniprotHumanCrap"$databaseExtension""$uniprotVersion".fasta -id_params ./databaseSearchParameters.par -out $PWD/$sampleTag.out.psdb"
    eval $peptideshakerExecution
    #reports output
    reportsExecution="$javaLocation -Xmx120g -cp $peptideshakerLocation eu.isas.peptideshaker.cmd.ReportCLI -in $sampleTag.out.psdb -out_reports $PWD -reports 0,3"
    eval $reportsExecution
  done
  printf "\nFinished all exporting reports for $sampleTag data set.\n\n"

  ####final output
  eval "rm ${PWD}/*.html"
  eval "scp ${startingDirectory}/dataProcessing.txt $PWD"
  eval "scp ${startingDirectory}/ccleDataProcessing_c57m16tmt11CidItHumanTrypsinCombined.sh $PWD"
  printf "\n\nCreating the processed data archive.\n\n"
  eval "tar -czvf ${sampleTag}_processedDataPackage.tar.gz $PWD/quantFiles *.fasta *.rds *.txt *.sh *.par"
  eval "mv ${sampleTag}_processedDataPackage.tar.gz $rawDataDirectory"
  printf "\n\nFinished processing a total of ${#fileList[@]} files!\n\n"

done
```

We can try this script and see how it works. I had to move the R2 and R3 sets and I will just process these after the fact on their own. 

### Update

I went through these data and I am getting far fewer peptide IDs for DLG2 specifically than what is originally reported in the manuscript. This seems odd. One error I found is that Comet is not working in the version of SearchGUI that I used, so this will need to be repeated. But, it makes me wonder if the search pipeline isn't the best. Below, I write a script for a new pipeline based on [Monocle](https://github.com/gygilab/Monocle), [Comet](http://comet-ms.sourceforge.net/release/release_202001/), and [Percolator](https://github.com/percolator).


```shell
#!/bin/bash
rawDataDirectory="/projects/ptx_results/OtherDataSets/dataset20201217_ccleProteomicsPmid31978347/searchTest/"

##software tools
monoclePath="/projects/ptx_analysis/chughes/software/monocle-0.3.37.19/Monocle.CLI/Monocle.CLI"

##processing files with Monocle
for i in *.raw
do
  echo $i
  eval $monoclePath -f $i -t "mzxml" -o ${rawDataDirectory}${i}.mzXML
done
```








############################################
for j in Prot_{01..42}  
do
  echo $j
  ##############################################################
  ##############################################################
  #users must edit the below variables
  #this is the text that will identify your files
  sampleTag=$j
  #this is the desired location for the output of the data processing process
  desiredBaseLocation=$rawDataDirectory
  #this is the base location where your raw data is stored
  dataStorageLocation=$rawDataDirectory
  #this is text that will be appended to your output folder that is created for this analysis
  folderAdapter="dataProcessing_"
  #users do not need to edit the statement below
  folderToCreate="$desiredBaseLocation$folderAdapter$sampleTag"
  #starting directory
  startingDirectory="${PWD}"




```
