## Reprocessing some proteomics data

This document describes the reprocessing of some proteomics data from a previous publication. Specifically:

"A deep proteome and transcriptome abundance atlas of 29 healthy human tissues"
Molecular Systems Biology, 2019, Pubmed ID: 30777892, ProteomeXchange: PXD010154

### Getting the raw data

They have already processed the data using MaxQuant, which I think is good enough. There is a lot of data and I do not want to re-process it all at this point. The files I am interested in on the server are:

```
P010747_Tonsil_fullproteome_RNAseq_txt.zip
P012502_Liver_fullproteome_RNAseq_txt.zip	
P013114_Spleenl_fullproteome_RNAseq_txt.zip	
P013128_Stomach_fullproteome_RNAseq_txt.zip
P013129_Brain_fullproteome_RNAseq_txt.zip	
P013163_Lung_fullproteome_RNAseq_txt.zip	
P013164_Testis_fullproteome_RNAseq_txt.zip	
P013187_Duodenum_fullproteome_RNAseq_txt.zip	
P013188_SmallIntestine_fullproteome_RNAseq_txt.zip	
P013196_UrinaryBladder_fullproteome_RNAseq_txt.zip
P013197_GallBladder_fullproteome_RNAseq_txt.zip	
P013198_Esophagus_fullproteome_RNAseq_txt.zip	
P013201_Heart_fullproteome_RNAseq_txt.zip	
P013385_Thyroid_fullproteome_RNAseq_txt.zip	
P013386_Endometrium_fullproteome_RNAseq_txt.zip	
P013387_Colon_fullproteome_RNAseq_txt.zip	
P013559_FallopianTube_fullproteome_RNAseq_txt.zip	
P013560_Kidney_fullproteome_RNAseq_txt.zip	
P013562_SmoothMuscle_fullproteome_RNAseq_txt.zip
P013675_Prostate_fullproteome_RNAseq_txt.zip
P013677_Appendix_fullproteome_RNAseq_txt.zip	
P013678_Pancreas_fullproteome_RNAseq_txt.zip
P013679_Ovary_fullproteome_RNAseq_txt.zip	
P013680_Placenta_fullproteome_RNAseq_txt.zip	
P013681_Rectum_fullproteome_RNAseq_txt.zip	
P015159_Fat_fullproteome_RNAseq_txt.zip	609 MB
P015160_LymophNode_fullproteome_RNAseq_txt.zip
P015161_SalivaryGland_fullproteome_RNAseq_txt.zip
P015424_AdrenalGland_fullproteome_RNAseq_txt.zip
```

I am just going to do this download and directory organization using a shell script. 



```shell
#!/bin/bash
rawDataOutputDirectory="/projects/ptx_results/OtherDataSets/dataset20190218_normalTissueProteomicsPmid30777892/maxquantResults/"
ftpLocation="ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/07/PXD010154/"
for i in P010747_Tonsil P012502_Liver P013114_Spleenl P013128_Stomach P013129_Brain P013163_Lung P013164_Testis P013187_Duodenum P013188_SmallIntestine P013196_UrinaryBladder P013197_GallBladder P013198_Esophagus P013201_Heart P013385_Thyroid P013386_Endometrium P013387_Colon P013559_FallopianTube P013560_Kidney P013562_SmoothMuscle P013675_Prostate P013677_Appendix P013678_Pancreas P013679_Ovary P013680_Placenta P013681_Rectum P015159_Fat P015160_LymophNode P015161_SalivaryGland P015424_AdrenalGland
do
  filenameBuild="${ftpLocation}${i}_fullproteome_RNAseq_txt.zip"
  echo $filenameBuild
  eval "cd ${rawDataOutputDirectory}"
  eval "wget $filenameBuild"
  eval "unzip ${i}_fullproteome_RNAseq_txt.zip"
  eval "mv ./txt ./${i}"
  eval "scp ${rawDataOutputDirectory}${i}/*_proteinGroups.txt ./"
done
```

Ok, this works fine. I will process the downloaded protein data in R in another script in this directory.
