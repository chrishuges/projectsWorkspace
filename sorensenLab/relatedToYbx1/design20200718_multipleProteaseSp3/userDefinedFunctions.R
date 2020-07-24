#######################################################################
#######################################################################
##this is a modified function from the OrgMassSpecR package to do
##peptide cleavage of other enzyme sites
#######################################################################
cleave = function(sequence, start, stop, misses, minLength, maxLength, proteinIndex, ...) {
  peptide = substring(sequence, start, stop)
  mc = rep(misses, times = length(peptide))
  pepLength = nchar(peptide)
  result = data.frame('seqnames' = paste('pro', proteinIndex, sep = ''), peptide, start, stop, mc, pepLength, stringsAsFactors = FALSE)
  resultSub = subset(result, result$pepLength >= minLength & result$pepLength <= maxLength)
  return(resultSub)
}
#######################################################################
#######################################################################



#######################################################################
#######################################################################
##this function is for calculating sequence coverage when
##using alternative sequence sites
#######################################################################
enzymeCoverage = function(sequence, cutSite, accessionCounter = 1, ...){
  cleavageData = data.frame()
  combinedCut = vector()
  seqVector = strsplit(sequence, split = '')[[1]]
  aaSeqEndPosition = length(seqVector)
  parentProtein = makeGRangesFromDataFrame(data.frame('seqnames' = paste('pro', accessionCounter, sep = ''), start = 1, stop = aaSeqEndPosition))
  cutAminoAcids = strsplit(cutSite, split = '')[[1]]
  ##
  for (seqCounter in 1:length(cutAminoAcids)){
    seqVector = strsplit(sequence, split = '')[[1]]
    aaSeqEndPosition = length(seqVector)
    aaSeqStop = grep(cutAminoAcids[seqCounter], seqVector)
    aaSeqStart = aaSeqStop + 1
    aaSeqStop = c(aaSeqStop, aaSeqEndPosition)
    aaSeqStart = c(1, aaSeqStart)
    seqCut = cleave(sequence, aaSeqStart, aaSeqStop, misses = 0, minLength = 6, maxLength = 30, proteinIndex = accessionCounter)
    if (nrow(seqCut) > 0){
      queryProtein = makeGRangesFromDataFrame(seqCut)
      proteinOverlap = as_tibble(GenomicRanges::setdiff(parentProtein, queryProtein))
      seqCut$coverage = round(100 - (sum(proteinOverlap$width, na.rm = TRUE) / aaSeqEndPosition) * 100, 1)
      seqCut$cutSite = cutAminoAcids[seqCounter]
      cleavageData = rbind(cleavageData, seqCut)
    }
    if (nrow(seqCut) == 0){
      print(paste('No peptides at cut site ', cutAminoAcids[seqCounter], '.', sep = ''))
    }
    combinedCut = paste(combinedCut, cutAminoAcids[seqCounter], '|', sep = '')
  }
  ##
  aaSeqStop = grep(substr(combinedCut, 1, nchar(combinedCut)-1), seqVector)
  aaSeqStart = aaSeqStop + 1
  aaSeqStop = c(aaSeqStop, aaSeqEndPosition)
  aaSeqStart = c(1, aaSeqStart)
  seqCut = cleave(sequence, aaSeqStart, aaSeqStop, misses = 0, minLength = 6, maxLength = 30, proteinIndex = accessionCounter)
  if (nrow(seqCut) > 0){
    queryProtein = makeGRangesFromDataFrame(seqCut)
    proteinOverlap = as_tibble(GenomicRanges::setdiff(parentProtein, queryProtein))
    seqCut$coverage = round(100 - (sum(proteinOverlap$width, na.rm = TRUE) / aaSeqEndPosition) * 100, 1)
    seqCut$cutSite = substr(combinedCut, 1, nchar(combinedCut)-1)
    cleavageData = rbind(cleavageData, seqCut)
  }
  if (nrow(seqCut) == 0){
    print(paste('No peptides at cut site ', substr(combinedCut, 1, nchar(combinedCut)-1), '.', sep = ''))
  }
  ##
  if (nrow(cleavageData) > 0){
    queryProtein = makeGRangesFromDataFrame(cleavageData)
    proteinOverlap = as_tibble(GenomicRanges::setdiff(parentProtein, queryProtein))
    cleavageData$combinedCoverage = round(100 - (sum(proteinOverlap$width, na.rm = TRUE) / aaSeqEndPosition) * 100, 1)
    return(unique(cleavageData[,c(1,7:9)]))
  }
  if (nrow(cleavageData == 0)){
    return(NA)
  }
  ##
  
}
#######################################################################
#######################################################################






















