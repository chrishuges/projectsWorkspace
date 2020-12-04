#######################################################################
#######################################################################
##This function performs survival analysis for input genes.
##For reference on the format of expression and phenotypeData,
##please see the example analysis in getGeoDataForSurvivalAnalysis.Rmd
##that can be found on this GitHub page.
#######################################################################
survivalAnalysis = function(geneOfInterest, expressionData, phenotypeData,
                            survivalPlot = FALSE, survivalData = FALSE, writeDirectory = '', printData = FALSE, ...){
  print(paste('Performing analysis for', geneOfInterest))
  ########
  #first we select our gene of interest, and sort based on median expression across all patients
  geneOfInterestExprs = expressionData %>%
    filter(symbol == geneOfInterest)
  ##check for data
  if (nrow(geneOfInterestExprs) < 1){
    print('Gene not found in first pass, searching secondary identifier annotation.')
    geneOfInterestExprs = expressionData %>%
      filter(arraySymbolOne == geneOfInterest)
  }
  if (nrow(geneOfInterestExprs) < 1){
    print('Gene not found in second pass, searching third identifier annotation.')
    geneOfInterestExprs = expressionData %>%
      filter(arraySymbolTwo == geneOfInterest)
  }
  if (nrow(geneOfInterestExprs) < 1){
    print('Gene not found. Moving on to next entry.')
    return(NA)
  }
  geneOfInterestExprs$medExprs = apply(geneOfInterestExprs[,which(grepl('GS',colnames(geneOfInterestExprs)))[1]:ncol(geneOfInterestExprs)], 
                                       1, function(x) mean(x, na.rm = TRUE))
  geneOfInterestSort = geneOfInterestExprs %>%
    arrange(desc(medExprs))
  ########
  #keep the top probe for the gene, based on the expression calculated above
  geneSurvivalInput = geneOfInterestSort[1,1:(ncol(geneOfInterestSort) - 1)] %>%
    pivot_longer(cols = colnames(geneOfInterestSort)[which(grepl('GS',colnames(geneOfInterestSort)))[1]]:colnames(geneOfInterestSort)[(ncol(geneOfInterestSort) - 1)], 
                 names_to = 'geo_accession', values_to = 'rnaExprs') %>%
    right_join(phenotypeData) %>%
    arrange(desc(rnaExprs)) 
  ########
  #sorting of the patients into low and high expression based on the top and bottom 25%
  geneSurvivalInput$geneLevel = ''
  geneSurvivalInput[1:round(nrow(geneSurvivalInput) * 0.25), 'geneLevel'] = 'high'
  geneSurvivalInput[round(nrow(geneSurvivalInput) * 0.75):nrow(geneSurvivalInput), 'geneLevel'] = 'low'
  geneSurvivalInput$geneLevel = factor(geneSurvivalInput$geneLevel, levels = c('low','high'))
  ########
  #calculation of the different survival metrics based on our data
  survivalFit = survfit(Surv(ovs, status) ~ geneLevel, data = geneSurvivalInput)
  survivalDiff = survdiff(Surv(ovs, status) ~ geneLevel, data = geneSurvivalInput)
  survivalPValue = 1 - pchisq(survivalDiff$chisq, length(survivalDiff$n) - 1)
  survivalSummary = surv_summary(survivalFit, data = geneSurvivalInput)
  coxStats = coxph(Surv(ovs, status) ~ geneLevel, data = geneSurvivalInput)
  coxZScore = coef(coxStats)/sqrt(diag(vcov(coxStats)))
  ########
  #this will output a survival plot
  if (survivalPlot == TRUE){
    ggsurvplot(survivalSummary, pval = survivalPValue, conf.int = FALSE,
               risk.table = FALSE, linetype = "strata", ggtheme = theme_classic(), palette = brewer.pal(8,'RdBu')[c(8,1)]
    )
    ggsave(paste(baseRepository, writeDirectory, '/survival_', 
                 geneOfInterest, '.pdf', sep = ''), width = 4, height = 4, useDingbats = FALSE)
  }
  ########
  #this will output a text file with the survival results
  if (survivalData == TRUE){
    survivalOutput = summary(survivalFit)$table
    write.table(survivalOutput, paste(baseRepository, writeDirectory, '/survival_', geneOfInterest, '.csv', sep = ''),
                col.names = TRUE, row.names = TRUE, quote = FALSE, sep = ',')
  }
  ########
  #print the data if desired
  if (printData == TRUE){
    print(survivalPValue)
    print(table(geneSurvivalInput$geneLevel))
  }
  ########
  #lastly we output the Z-score
  return(coxZScore)
  
}
#######################################################################
#######################################################################



#######################################################################
##This function is a variant of the original survivalAnalysis,
##but instead tries to determine the best expression cutoff between
##high and low cases
##more details here: http://r-addict.com/2016/11/21/Optimal-Cutpoint-maxstat.html
#######################################################################
postelExpression = readRDS(paste(baseRepository, '/microarray20200604_ewingSarcomaAllSets/dataset_postelVinayExpression.rds', sep = ''))
postelPhenotype = readRDS(paste(baseRepository, '/microarray20200604_ewingSarcomaAllSets/dataset_postelVinayPhenotype.rds', sep = ''))
expressionData = postelExpression
phenotypeData = postelPhenotype
geneOfInterest = 'YBX1'



survivalAnalysisV2 = function(geneOfInterest, expressionData, phenotypeData,
                              survivalPlot = FALSE, survivalData = FALSE, writeDirectory = '', printData = FALSE, ...){
  print(paste('Performing analysis for', geneOfInterest))
  ########
  #first we select our gene of interest, and sort based on median expression across all patients
  geneOfInterestExprs = expressionData %>%
    filter(symbol == geneOfInterest)
  ##check for data
  if (nrow(geneOfInterestExprs) < 1){
    print('Gene not found in first pass, searching secondary identifier annotation.')
    geneOfInterestExprs = expressionData %>%
      filter(arraySymbolOne == geneOfInterest)
  }
  if (nrow(geneOfInterestExprs) < 1){
    print('Gene not found in second pass, searching third identifier annotation.')
    geneOfInterestExprs = expressionData %>%
      filter(arraySymbolTwo == geneOfInterest)
  }
  if (nrow(geneOfInterestExprs) < 1){
    print('Gene not found. Moving on to next entry.')
    return(NA)
  }
  geneOfInterestExprs$medExprs = apply(geneOfInterestExprs[,which(grepl('GS',colnames(geneOfInterestExprs)))[1]:ncol(geneOfInterestExprs)], 
                                       1, function(x) mean(x, na.rm = TRUE))
  geneOfInterestSort = geneOfInterestExprs %>%
    arrange(desc(medExprs))
  ########
  #keep the top probe for the gene, based on the expression calculated above
  geneSurvivalInput = geneOfInterestSort[1,1:(ncol(geneOfInterestSort) - 1)] %>%
    pivot_longer(cols = colnames(geneOfInterestSort)[which(grepl('GS',colnames(geneOfInterestSort)))[1]]:colnames(geneOfInterestSort)[(ncol(geneOfInterestSort) - 1)], 
                 names_to = 'geo_accession', values_to = 'rnaExprs') %>%
    right_join(phenotypeData) %>%
    arrange(rnaExprs)
  ########
  geneCutPoint = surv_cutpoint(geneSurvivalInput,
                               time = 'ovs',
                               event = 'status',
                               variables = 'rnaExprs')
  geneSurvivalCat = surv_categorize(geneCutPoint)
  geneSurvivalInput$geneLevel = geneSurvivalCat$rnaExprs
  geneSurvivalInput$geneLevel = factor(geneSurvivalInput$geneLevel, levels = c('low','high'))
  #print(table(geneSurvivalInput$geneLevel))
  ########
  #calculation of the different survival metrics based on our data
  survivalFit = survfit(Surv(ovs, status) ~ geneLevel, data = geneSurvivalInput)
  survivalDiff = survdiff(Surv(ovs, status) ~ geneLevel, data = geneSurvivalInput)
  survivalPValue = 1 - pchisq(survivalDiff$chisq, length(survivalDiff$n) - 1)
  #print(survivalPValue)
  survivalSummary = surv_summary(survivalFit, data = geneSurvivalInput)
  coxStats = coxph(Surv(ovs, status) ~ geneLevel, data = geneSurvivalInput)
  coxZScore = coef(coxStats)/sqrt(diag(vcov(coxStats)))
  ########
  #this will output a survival plot
  if (survivalPlot == TRUE){
    ggsurvplot(survivalSummary, pval = survivalPValue, conf.int = FALSE,
               risk.table = FALSE, linetype = "strata", ggtheme = theme_classic(), palette = brewer.pal(8,'RdBu')[c(8,1)]
    )
    ggsave(paste(baseRepository, writeDirectory, '/survival_', 
                 geneOfInterest, '.pdf', sep = ''), width = 4, height = 4, useDingbats = FALSE)
  }
  ########
  #this will output a text file with the survival results
  if (survivalData == TRUE){
    survivalOutput = summary(survivalFit)$table
    write.table(survivalOutput, paste(baseRepository, writeDirectory, '/survival_', geneOfInterest, '.csv', sep = ''),
                col.names = TRUE, row.names = TRUE, quote = FALSE, sep = ',')
  }
  ########
  #print the data if desired
  if (printData == TRUE){
    print(survivalPValue)
    print(print(table(geneSurvivalInput$geneLevel)))
  }
  ########
  #lastly we output the Z-score
  return(coxZScore)
  
}

