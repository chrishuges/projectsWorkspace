#######################################################################
#######################################################################
##this function creates survival curves for all input genes
##based on the input data specified in the file
##20200604File1_geoDataProcessing.Rmd
#######################################################################
#expressionData = esAnnotated
#phenotypeData = esPhenotype
#geneOfInterest = 'YBX1'

survivalAnalysis = function(geneOfInterest, expressionData, phenotypeData, ...){
  geneOfInterestExprs = expressionData %>%
    filter(grepl(geneOfInterest, symbol))
  
  for (i in 1:nrow(geneOfInterestExprs)){
    message(paste(nrow(geneOfInterestExprs), 'total probes for', geneOfInterest))
    message(paste('Working on probe',i))
    geneLong = geneOfInterestExprs[i,] %>%
      pivot_longer(cols = GSM439886:GSM439940, names_to = 'geo_accession', values_to = 'rnaExprs') %>%
      arrange(desc(rnaExprs))
    geneLong$geneLevel = 'medium'
    geneLong[1:round(nrow(geneLong) * 0.25), 'geneLevel'] = 'high'
    geneLong[round(nrow(geneLong) * 0.75):nrow(geneLong), 'geneLevel'] = 'low'
    geneLong$geneLevel = factor(geneLong$geneLevel, levels = c('low','medium','high'))
    geneSurvivalInput = phenotypeData %>%
      left_join(geneLong)
    survivalFit = survfit(Surv(ovs, status) ~ geneLevel, data = geneSurvivalInput)
    coxStats = coxph(Surv(ovs, status) ~ geneLevel, data = geneSurvivalInput)
    coxZScore = coef(coxStats)/sqrt(diag(vcov(coxStats)))
  }
}
#######################################################################
#######################################################################



message(coxStats$score)
ggsurvplot(survivalFit, pval = FALSE, conf.int = FALSE,
           linetype = "strata", # Change line type by groups
           ggtheme = theme_classic(), 
           palette = c("#E7B800", "#2E9FDF", 'red'))
ggsave(paste(baseWorkspace, '/microarray20200604_ewingSarcomaSavola/survival_', 
             geneSurvivalInput$symbol[1], '_', i, '.pdf', sep = ''), width = 4, height = 4, useDingbats = FALSE)
