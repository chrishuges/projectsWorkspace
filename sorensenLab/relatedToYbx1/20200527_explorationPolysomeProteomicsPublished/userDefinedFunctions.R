#######################################################################
#######################################################################
##this function takes a max quant data object specifically from the 
##polysome data and will parse it into a more usable data object
##this function is really designed for the Imami data
#######################################################################
parsingPolysomeMqData = function(mqData, fractionFilter = 1, addAnnotation = FALSE, annotationIndex, ...){
  #first we can rework the annotation columns
  message('Parsing the protein accessions.')
  mqData$accession = sapply(strsplit(as.character(mqData$`Protein IDs`),";"), '[', 1)
  mqData$accession = sapply(strsplit(as.character(mqData$accession),"-"), '[', 1)
  mqData$gene = sapply(strsplit(as.character(mqData$`Gene names`),";"), '[', 1)
  mqData = mqData[!is.na(mqData$gene), c('accession','gene', colnames(mqData)[3:74])]
  #separate out the two replicates
  message('Split the sample into replicates.')
  lightRep = mqData[,c(1,2,which(grepl('Ratio [M]\\/[L] HEK_fr', colnames(mqData))))]
  heavyRep = mqData[,c(1,2,which(grepl('Ratio [H]\\/[M] HEK_fr', colnames(mqData))))]
  #change the column names for each
  message('Change the column names to a better format.')
  colnames(lightRep) = c('accession','gene', paste('light_', sub('.*HEK_(.*)$', '\\1', colnames(lightRep)[3:38]), sep = ''))
  colnames(heavyRep) = c('accession','gene', paste('heavy_', sub('.*HEK_(.*)$', '\\1', colnames(heavyRep)[3:38]), sep = ''))
  #filter based on the fractionFilter
  lightRepFilter = subset(lightRep, rowSums(!is.na(lightRep[,3:38])) >= fractionFilter)
  heavyRepFilter = subset(heavyRep, rowSums(!is.na(heavyRep[,3:38])) >= fractionFilter)
  message(paste((nrow(lightRep) - nrow(lightRepFilter)), ' light proteins were removed based on the fractionFilter setting.', sep = ''))
  message(paste((nrow(heavyRep) - nrow(heavyRepFilter)), ' light proteins were removed based on the fractionFilter setting.', sep = ''))
  #now transform them to long format
  message('Transform the data into long format.')
  lightRepLong = lightRepFilter %>%
    pivot_longer(cols = light_fr05:light_fr40, names_to = 'sampleId', values_to = 'silacRatio') %>%
    mutate('fraction' = as.numeric(sub('.*fr(.*)$', '\\1', sampleId))) %>%
    mutate('replicate' = sub('(.*)_fr[0-9]+.*$', '\\1', sampleId))
  lightRepLong$fraction = factor(lightRepLong$fraction, levels = seq(5,40,1))
  heavyRepLong = heavyRepFilter %>%
    pivot_longer(cols = heavy_fr05:heavy_fr40, names_to = 'sampleId', values_to = 'silacRatio') %>%
    mutate('fraction' = as.numeric(sub('.*fr(.*)$', '\\1', sampleId))) %>%
    mutate('replicate' = sub('(.*)_fr[0-9]+.*$', '\\1', sampleId))
  heavyRepLong$fraction = factor(heavyRepLong$fraction, levels = seq(5,40,1))
  #the light replicate also needs to have its ratio flipped
  message('Flip the M/L ratio.')
  lightRepLong$silacRatio = 1 / lightRepLong$silacRatio
  #recombine the data
  polysomeData = rbind(lightRepLong, heavyRepLong)
  polysomeData$replicate = factor(polysomeData$replicate, levels = c('light','heavy'))
  #check if they want annotation
  if (addAnnotation == TRUE){
    message('Adding annotation.')
    polysomeDataAnnotated = polysomeData %>% left_join(annotationIndex)
    message('Output the final data.')
    return(polysomeDataAnnotated)
  }
  else {
    message('No annotation will be added.')
    message('Output the final data.')
    return(polysomeData)
  }
}
#######################################################################
#######################################################################



