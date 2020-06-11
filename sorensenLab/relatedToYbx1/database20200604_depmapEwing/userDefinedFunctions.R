#######################################################################
#######################################################################
##this function will get the correlation between a drug AUC values
##and the dependency scores on a gene-by-gene basis
#######################################################################
getDrugSensitivity = function(geneOfInterest, crisprData, annotationData, lineage, ...){
  ##
  crisprGene = crisprData[,c(1, which(grepl(geneOfInterest, 
                                            sub('(.*) \\([0-9]+\\)$', '\\1', colnames(crisprData)))))]
  colnames(crisprGene)[2] = geneOfInterest
  ##
  crisprLineage = crisprGene %>%
    left_join(annotationData) %>%
    dplyr::select(DepMap_ID, lineage2, area_under_curve, paste(geneOfInterest)) %>%
    filter(!is.na(area_under_curve) & grepl(lineage, lineage2))
  ##
  return(cor(crisprLineage[,c('area_under_curve', paste(geneOfInterest))], use = 'pairwise.complete.obs')[2])
}

#######################################################################
#######################################################################



