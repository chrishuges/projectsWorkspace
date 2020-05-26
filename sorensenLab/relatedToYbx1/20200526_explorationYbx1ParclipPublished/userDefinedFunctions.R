#######################################################################
#######################################################################
##this function is a support function that performs mapping
##of parsed regions to the genes and ID's they belong to. It supports
## the wavclusterPeakAnnotation function
#######################################################################
prepTxDb = function(txdbInput, txdbRaw, ...){
  
  ########first we will prefilter our txdb for the longest mRNA
  message('Preparing a gold standard set of RNA CDS Ids')
  cdsRaw = transcripts(txdbRaw)
  cdsLengths = as_tibble(transcriptLengths(txdbRaw, with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE))
  cdsLengthsSet1 = cdsLengths %>%
    dplyr::select(tx_name, gene_id, tx_len) %>%
    arrange(gene_id, desc(tx_len))
  cdsDuplicated = !duplicated(cdsLengthsSet1[,2])
  cdsNoDuplicates = cdsLengthsSet1[cdsDuplicated,]
  message('Core CDS set assembled')
  
  ######now process the specific set of interest
  message('Now setting up your transcriptome database of interest')
  txdbKeys = as.character(unique(names(txdbInput)))
  txdbAnno1 = select(txdbRaw, keys = txdbKeys, columns = c('GENEID','TXTYPE'), keytype = 'TXNAME')
  txdbAnno1$SYMBOL = mapIds(org.Hs.eg.db, keys = txdbAnno1$GENEID, column = 'SYMBOL', keytype = 'ENSEMBL')
  colnames(txdbAnno1) = c('ensembltrans','ensembl','tx_type','symbol')
  txdbAnno2 = tibble('ensembltrans' = names(txdbInput)) %>% left_join(txdbAnno1)
  txdbAnno3 = txdbInput
  txdbAnno3$symbol = txdbAnno2$symbol
  txdbAnno3$ensg = txdbAnno2$ensembl
  txdbAnno3$tx_type = txdbAnno2$tx_type
  txdbSetSc = keepStandardChromosomes(txdbAnno3, pruning.mode = 'coarse')
  txdbSetPc = subset(txdbSetSc, txdbSetSc$tx_type == 'protein_coding' & names(txdbSetSc) %in% cdsNoDuplicates$tx_name & !is.na(txdbSetSc$symbol))
  
  #####return the data
  message('Preparing data output')
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
  message('Annotating fiveUTR hits')
  txdbSet = unlist(fiveUTRsByTranscript(userTxdb, use.names = TRUE))
  txdbSetPc = prepTxDb(txdbSet, userTxdb)
  message('Annotation index for fiveUTR created')
  message('Annotating fiveUTR overlap hits')
  overlapHits = findOverlaps(prepDataGrange, txdbSetPc, type = 'any')
  txdbPeakData = as.data.frame(txdbSetPc[subjectHits(overlapHits)], row.names = NULL)
  fivePeakData = queryDataset[queryHits(overlapHits),]
  fivePeakData$txdbStart = txdbPeakData$start
  fivePeakData$txdbEnd = txdbPeakData$end
  fivePeakData$txdbWidth = txdbPeakData$width
  fivePeakData$txdbExonId = txdbPeakData$exon_id
  fivePeakData$txdbSymbol = txdbPeakData$symbol
  fivePeakData$txdbEnsg = txdbPeakData$ensg
  fivePeakData$txdbType = txdbPeakData$tx_type
  fivePeakDataSorted = fivePeakData[order(fivePeakData$txdbSymbol, fivePeakData$txdbExonId),]
  fivePeakDataSorted$region = 'fiveUTR'
  message('Annotation for fiveUTR complete')
  
  
  ####process CDS txdb index
  message('Annotating CDS hits')
  txdbSet = unlist(cdsBy(userTxdb, use.names = TRUE))
  txdbSetPc = prepTxDb(txdbSet, userTxdb)
  message('Annotation index for CDS created')
  message('Annotating CDS overlap hits')
  overlapHits = findOverlaps(prepDataGrange, txdbSetPc, type = 'any')
  txdbPeakData = as.data.frame(txdbSetPc[subjectHits(overlapHits)], row.names = NULL)
  cdsPeakData = queryDataset[queryHits(overlapHits),]
  cdsPeakData$txdbStart = txdbPeakData$start
  cdsPeakData$txdbEnd = txdbPeakData$end
  cdsPeakData$txdbWidth = txdbPeakData$width
  cdsPeakData$txdbExonId = txdbPeakData$cds_id
  cdsPeakData$txdbSymbol = txdbPeakData$symbol
  cdsPeakData$txdbEnsg = txdbPeakData$ensg
  cdsPeakData$txdbType = txdbPeakData$tx_type
  cdsPeakDataSorted = cdsPeakData[order(cdsPeakData$txdbSymbol, cdsPeakData$txdbExonId),]
  cdsPeakDataSorted$region = 'CDS'
  message('Annotation for CDS complete')
  
  
  
  ####process 3UTR txdb index
  message('Annotating threeUTR hits')
  txdbSet = unlist(threeUTRsByTranscript(userTxdb, use.names = TRUE))
  txdbSetPc = prepTxDb(txdbSet, userTxdb)
  message('Annotation index for threeUTR created')
  message('Annotating threeUTR overlap hits')
  overlapHits = findOverlaps(prepDataGrange, txdbSetPc, type = 'any')
  txdbPeakData = as.data.frame(txdbSetPc[subjectHits(overlapHits)], row.names = NULL)
  threePeakData = queryDataset[queryHits(overlapHits),]
  threePeakData$txdbStart = txdbPeakData$start
  threePeakData$txdbEnd = txdbPeakData$end
  threePeakData$txdbWidth = txdbPeakData$width
  threePeakData$txdbExonId = txdbPeakData$exon_id
  threePeakData$txdbSymbol = txdbPeakData$symbol
  threePeakData$txdbEnsg = txdbPeakData$ensg
  threePeakData$txdbType = txdbPeakData$tx_type
  threePeakDataSorted = threePeakData[order(threePeakData$txdbSymbol, threePeakData$txdbExonId),]
  threePeakDataSorted$region = 'threeUTR'
  message('Annotation for threeUTR complete')
  
  
  ####perform the location calculations
  message('Combining the data sets and performing the final calculations')
  gnSet = rbind(fivePeakDataSorted, cdsPeakDataSorted, threePeakDataSorted)
  gnSet$region = factor(gnSet$region, levels = c('fiveUTR','CDS','threeUTR'))
  gnSet$clusterMedian = apply(gnSet[,c('start','end')], 1, function(x) as.integer(median(x, na.rm = TRUE)))
  gnSet$regionLocation = (gnSet$clusterMedian - gnSet$txdbStart) / gnSet$txdbWidth
  gnSet$regionLocation = ifelse(grepl('-',gnSet$strand), 1 - gnSet$regionLocation, gnSet$regionLocation)
  gnSet$regionLocationAdjusted = ifelse(gnSet$region == 'CDS', gnSet$regionLocation + 1.1,
                                        ifelse(gnSet$region == 'threeUTR', gnSet$regionLocation + 2.2, gnSet$regionLocation))
  gnSetSub1 = subset(gnSet, !(gnSet$clusterMedian <= gnSet$txdbStart) & !(gnSet$clusterMedian >= gnSet$txdbEnd))
  message('Processing complete')
  return(gnSetSub1)
}
#######################################################################
#######################################################################

