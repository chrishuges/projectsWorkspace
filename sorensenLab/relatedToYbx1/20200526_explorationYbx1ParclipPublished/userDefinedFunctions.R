#######################################################################
#######################################################################
##this function takes a txdb as input and finds the longest
##transcript variant assigned to a single gene, and returns
##a dataframe of just the longest hits for each
#######################################################################
findLongTranscripts = function(txdbRaw, ...){
  message('Finding the longest transcript associated with each input gene in the txdb.')
  cdsRaw = transcripts(txdbRaw)
  cdsLengths = as_tibble(transcriptLengths(txdbRaw, with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE))
  cdsLengthsSet1 = cdsLengths %>%
    dplyr::select(tx_name, gene_id, tx_len) %>%
    arrange(gene_id, desc(tx_len))
  cdsDuplicated = !duplicated(cdsLengthsSet1[,2])
  cdsNoDuplicates = cdsLengthsSet1[cdsDuplicated,]
  message('Core set of long transcripts assembled. Data returned.')
  return(cdsNoDuplicates)
}
#######################################################################
#######################################################################


#######################################################################
#######################################################################
##this function performs annotation of an input txdb
##and will filter for longest transcripts based on a given set
#######################################################################
annotateTxdb = function(txdbInput, txdbRaw, longTranscripts, ...){
  message('Annotating the input txdb.')
  txdbKeys = as.character(unique(names(txdbInput)))
  txdbAnno1 = select(txdbRaw, keys = txdbKeys, columns = c('GENEID','TXTYPE'), keytype = 'TXNAME')
  txdbAnno1$SYMBOL = mapIds(org.Hs.eg.db, keys = txdbAnno1$GENEID, column = 'SYMBOL', keytype = 'ENSEMBL')
  colnames(txdbAnno1) = c('ensembltrans','ensembl','tx_type','symbol')
  txdbAnno2 = tibble('ensembltrans' = names(txdbInput)) %>% left_join(txdbAnno1)
  txdbAnno3 = txdbInput
  txdbAnno3$symbol = txdbAnno2$symbol
  txdbAnno3$ensg = txdbAnno2$ensembl
  txdbAnno3$tx_type = txdbAnno2$tx_type
  message('Filtering input txdb to keep only longest transcripts.')
  txdbSetSc = keepStandardChromosomes(txdbAnno3, pruning.mode = 'coarse')
  txdbSetPc = subset(txdbSetSc, names(txdbSetSc) %in% longTranscripts$tx_name)
  message('Returning the annotated txdb for alignment.')
  return(txdbSetPc)
}
#######################################################################
#######################################################################



#######################################################################
#######################################################################
##this function takes an output file from wavClusteR (wavClusters.rds)
##and will perform annotation of the peaks based on the Human ENSEMBL
##transcriptome database.
#######################################################################
wavclusterPeakAnnotation = function(queryDataset, userTxdb, removeMito = FALSE, ...){
  
  ####change the chromosome annotation and remove chrM if desired
  prepData = queryDataset
  prepData$seqnames = substr(prepData$seqnames,4,nchar(as.character(prepData$seqnames)))
  if (removeMito == TRUE){
    prepData = subset(prepData, !grepl('M', prepData$seqnames))
    message('Removing chromosome M hits.')
  } else {
    prepData$seqnames = ifelse(grepl('M', prepData$seqnames), 'MT', prepData$seqnames)
    message('Chromosome M hits are retained.')
  }
  ####now make the GRange for processing
  prepDataGrange = makeGRangesFromDataFrame(prepData)
  
  ####process 5UTR txdb index
  message('Annotating fiveUTR hits.')
  txdbSet = unlist(fiveUTRsByTranscript(userTxdb, use.names = TRUE))
  txdbLongTranscripts = findLongTranscripts(userTxdb)
  txdbSetPc = annotateTxdb(txdbSet, userTxdb, txdbLongTranscripts)
  message('Annotation index for fiveUTR created.')
  message('Annotating fiveUTR overlap hits.')
  overlapHits = findOverlaps(prepDataGrange, txdbSetPc, type = 'any')
  txdbPeakData = as.data.frame(txdbSetPc[subjectHits(overlapHits)], row.names = NULL)
  fivePeakData = prepData[queryHits(overlapHits),]
  fivePeakData$txdbStart = txdbPeakData$start
  fivePeakData$txdbEnd = txdbPeakData$end
  fivePeakData$txdbWidth = txdbPeakData$width
  fivePeakData$txdbExonId = txdbPeakData$exon_id
  fivePeakData$txdbSymbol = txdbPeakData$symbol
  fivePeakData$txdbEnsg = txdbPeakData$ensg
  fivePeakData$txdbType = txdbPeakData$tx_type
  fivePeakDataSorted = fivePeakData[order(fivePeakData$txdbSymbol, fivePeakData$txdbExonId),]
  fivePeakDataSorted$region = 'fiveUTR'
  message('Annotation for fiveUTR complete.')
  
  
  ####process CDS txdb index
  message('Annotating CDS hits.')
  txdbSet = unlist(cdsBy(userTxdb, use.names = TRUE))
  txdbSetPc = annotateTxdb(txdbSet, userTxdb, txdbLongTranscripts)
  message('Annotation index for CDS created.')
  message('Annotating CDS overlap hits.')
  overlapHits = findOverlaps(prepDataGrange, txdbSetPc, type = 'any')
  txdbPeakData = as.data.frame(txdbSetPc[subjectHits(overlapHits)], row.names = NULL)
  cdsPeakData = prepData[queryHits(overlapHits),]
  cdsPeakData$txdbStart = txdbPeakData$start
  cdsPeakData$txdbEnd = txdbPeakData$end
  cdsPeakData$txdbWidth = txdbPeakData$width
  cdsPeakData$txdbExonId = txdbPeakData$cds_id
  cdsPeakData$txdbSymbol = txdbPeakData$symbol
  cdsPeakData$txdbEnsg = txdbPeakData$ensg
  cdsPeakData$txdbType = txdbPeakData$tx_type
  cdsPeakDataSorted = cdsPeakData[order(cdsPeakData$txdbSymbol, cdsPeakData$txdbExonId),]
  cdsPeakDataSorted$region = 'CDS'
  message('Annotation for CDS complete.')
  
  
  
  ####process 3UTR txdb index
  message('Annotating threeUTR hits.')
  txdbSet = unlist(threeUTRsByTranscript(userTxdb, use.names = TRUE))
  txdbSetPc = annotateTxdb(txdbSet, userTxdb, txdbLongTranscripts)
  message('Annotation index for threeUTR created.')
  message('Annotating threeUTR overlap hits.')
  overlapHits = findOverlaps(prepDataGrange, txdbSetPc, type = 'any')
  txdbPeakData = as.data.frame(txdbSetPc[subjectHits(overlapHits)], row.names = NULL)
  threePeakData = prepData[queryHits(overlapHits),]
  threePeakData$txdbStart = txdbPeakData$start
  threePeakData$txdbEnd = txdbPeakData$end
  threePeakData$txdbWidth = txdbPeakData$width
  threePeakData$txdbExonId = txdbPeakData$exon_id
  threePeakData$txdbSymbol = txdbPeakData$symbol
  threePeakData$txdbEnsg = txdbPeakData$ensg
  threePeakData$txdbType = txdbPeakData$tx_type
  threePeakDataSorted = threePeakData[order(threePeakData$txdbSymbol, threePeakData$txdbExonId),]
  threePeakDataSorted$region = 'threeUTR'
  message('Annotation for threeUTR complete.')
  
  
  ####perform the location calculations
  message('Combining the data sets and performing the final calculations')
  gnSet = rbind(fivePeakDataSorted, cdsPeakDataSorted, threePeakDataSorted)
  gnSet$region = factor(gnSet$region, levels = c('fiveUTR','CDS','threeUTR'))
  gnSet$clusterMedian = apply(gnSet[,c('start','end')], 1, function(x) as.integer(median(x, na.rm = TRUE)))
  gnSet$regionLocation = (gnSet$clusterMedian - gnSet$txdbStart) / gnSet$txdbWidth
  gnSet$regionLocation = ifelse(grepl('-',gnSet$strand), 1 - gnSet$regionLocation, gnSet$regionLocation)
  gnSet$regionLocationAdjusted = ifelse(gnSet$region == 'CDS', gnSet$regionLocation + 1.1,
                                        ifelse(gnSet$region == 'threeUTR', gnSet$regionLocation + 2.2, gnSet$regionLocation))
  #gnSetSub1 = subset(gnSet, !(gnSet$clusterMedian <= gnSet$txdbStart) & !(gnSet$clusterMedian >= gnSet$txdbEnd))
  #revert the chromosomes to their original format
  gnSet$seqnames = paste('chr', gnSet$seqnames, sep = '')
  gnSetSub1 = subset(gnSet, !is.na(gnSet$txdbSymbol))
  message('Processing complete')
  return(gnSetSub1)
}
#######################################################################
#######################################################################

