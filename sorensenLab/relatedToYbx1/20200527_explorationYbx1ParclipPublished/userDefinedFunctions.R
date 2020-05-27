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
##this function will calculate the length of all regions
##that are fed to it in the form of an annotated wavClusteR output
##it requires you to specify a txdb to use as a reference as well
##as a cutoff for readCounts, if desired
#######################################################################
geneLengthAnalysis = function(annotatedPeaks, readCountSum = 1, userTxdb, ...){
  
  txdbLongTranscripts = findLongTranscripts(userTxdb)
  ####process 5UTR txdb index
  message('Annotating fiveUTR hits.')
  gnSet = annotatedPeaks[annotatedPeaks$region == 'fiveUTR' & !is.na(annotatedPeaks$txdbSymbol) & 
                           annotatedPeaks$readCountSum >= readCountSum,]
  txdbSet = unlist(fiveUTRsByTranscript(userTxdb, use.names = TRUE))
  txdbSetPc = annotateTxdb(txdbSet, userTxdb, txdbLongTranscripts)
  txdbSetPc$status = ifelse(txdbSetPc$symbol %in% gnSet$txdbSymbol, 'bound', 'unbound')
  fiveUtrWidths = tibble('symbol' = txdbSetPc$symbol,
                         'status' = txdbSetPc$status,
                         'width' = width(txdbSetPc),
                         'region' = 'fiveUTR') %>%
    filter(!is.na(symbol)) %>%
    group_by(symbol, status, region) %>%
    summarize('medianWidth' = median(width, na.rm = TRUE))
  
  ####process CDS txdb index
  message('Annotating CDS hits.')
  gnSet = annotatedPeaks[annotatedPeaks$region == 'CDS' & !is.na(annotatedPeaks$txdbSymbol) & 
                           annotatedPeaks$readCountSum >= readCountSum,]
  txdbSet = unlist(cdsBy(userTxdb, use.names = TRUE))
  txdbSetPc = annotateTxdb(txdbSet, userTxdb, txdbLongTranscripts)
  txdbSetPc$status = ifelse(txdbSetPc$symbol %in% gnSet$txdbSymbol, 'bound', 'unbound')
  cdsWidths = tibble('symbol' = txdbSetPc$symbol,
                     'status' = txdbSetPc$status,
                     'width' = width(txdbSetPc),
                     'region' = 'CDS') %>%
    filter(!is.na(symbol)) %>%
    group_by(symbol, status, region) %>%
    summarize('medianWidth' = median(width, na.rm = TRUE))
  
  ####process 3UTR txdb index
  message('Annotating threeUTR hits.')
  gnSet = annotatedPeaks[annotatedPeaks$region == 'threeUTR' & !is.na(annotatedPeaks$txdbSymbol) & 
                           annotatedPeaks$readCountSum >= readCountSum,]
  txdbSet = unlist(threeUTRsByTranscript(userTxdb, use.names = TRUE))
  txdbSetPc = annotateTxdb(txdbSet, userTxdb, txdbLongTranscripts)
  txdbSetPc$status = ifelse(txdbSetPc$symbol %in% gnSet$txdbSymbol, 'bound', 'unbound')
  threeUtrWidths = tibble('symbol' = txdbSetPc$symbol,
                          'status' = txdbSetPc$status,
                          'width' = width(txdbSetPc),
                          'region' = 'threeUTR') %>%
    filter(!is.na(symbol)) %>%
    group_by(symbol, status, region) %>%
    summarize('medianWidth' = median(width, na.rm = TRUE))
  
  ####combine the data
  message('Performing final assembly.')
  allData = rbind(fiveUtrWidths, cdsWidths, threeUtrWidths)
  allData$region = factor(allData$region, levels = c('fiveUTR','CDS','threeUTR'))
  allData$status = factor(allData$status, levels = c('unbound','bound'))
  return(allData)
}
#######################################################################
#######################################################################



#######################################################################
#######################################################################
##this function will calculate the GC content of all regions
##that are fed to it in the form of an annotated wavClusteR output
##it requires you to specify a txdb to use as a reference as well
##as a cutoff for readCounts, if desired
#######################################################################
gcContentAnalysis = function(annotatedPeaks, readCountSum = 1, userTxdb, ...){
  
  txdbLongTranscripts = findLongTranscripts(userTxdb)
  ####process 5UTR txdb index
  message('Annotating GC content for fiveUTR hits')
  gnSet = annotatedPeaks[annotatedPeaks$region == 'fiveUTR' & !is.na(annotatedPeaks$txdbSymbol) & 
                           annotatedPeaks$readCountSum >= readCountSum,]
  txdbSet = unlist(fiveUTRsByTranscript(userTxdb, use.names = TRUE))
  txdbSetPc = annotateTxdb(txdbSet, userTxdb, txdbLongTranscripts)
  message('Annotation index created for fiveUTR')
  txdbEnsSeq = as_tibble(txdbSetPc)
  txdbEnsSeq$seqnames = paste('chr',txdbEnsSeq$seqnames, sep = '')
  txdbEnsSeq$seqnames = ifelse(grepl('chrMT',txdbEnsSeq$seqnames),'chrM',txdbEnsSeq$seqnames)
  txdbEnsSeqGRange = makeGRangesFromDataFrame(txdbEnsSeq) 
  txdbEnsSequences = getSeq(Hsapiens, txdbEnsSeqGRange)
  txdbSetPc$gcContent = round(letterFrequency(txdbEnsSequences, "GC", as.prob = TRUE),2)
  message('GC content for index calculated')
  fiveUtrGcContent = as_tibble(txdbSetPc) %>%
    dplyr::select(symbol, strand, G.C) %>%
    filter(!is.na(symbol)) %>%
    group_by(symbol) %>%
    summarize('gcContent' = median(G.C, na.rm = TRUE))
  fiveUtrGcContent$region = 'fiveUTR'
  fiveUtrGcContent$status = ifelse(fiveUtrGcContent$symbol %in% gnSet$txdbSymbol, 'bound', 'unbound')
  
  ####process CDS txdb index
  message('Annotating GC content for CDS hits')
  gnSet = annotatedPeaks[annotatedPeaks$region == 'CDS' & !is.na(annotatedPeaks$txdbSymbol) & 
                           annotatedPeaks$readCountSum >= readCountSum,]
  txdbSet = unlist(cdsBy(userTxdb, use.names = TRUE))
  txdbSetPc = annotateTxdb(txdbSet, userTxdb, txdbLongTranscripts)
  message('Annotation index created for CDS')
  txdbEnsSeq = as_tibble(txdbSetPc)
  txdbEnsSeq$seqnames = paste('chr',txdbEnsSeq$seqnames, sep = '')
  txdbEnsSeq$seqnames = ifelse(grepl('chrMT',txdbEnsSeq$seqnames),'chrM',txdbEnsSeq$seqnames)
  txdbEnsSeqGRange = makeGRangesFromDataFrame(txdbEnsSeq) 
  txdbEnsSequences = getSeq(Hsapiens, txdbEnsSeqGRange)
  txdbSetPc$gcContent = round(letterFrequency(txdbEnsSequences, "GC", as.prob = TRUE),2)
  message('GC content for index calculated')
  cdsGcContent = as_tibble(txdbSetPc) %>%
    dplyr::select(symbol, strand, G.C) %>%
    filter(!is.na(symbol)) %>%
    group_by(symbol) %>%
    summarize('gcContent' = median(G.C, na.rm = TRUE))
  cdsGcContent$region = 'CDS'
  cdsGcContent$status = ifelse(cdsGcContent$symbol %in% gnSet$txdbSymbol, 'bound', 'unbound')
  
  ####process 3UTR txdb index
  message('Annotating GC content for threeUTR hits')
  gnSet = annotatedPeaks[annotatedPeaks$region == 'threeUTR' & !is.na(annotatedPeaks$txdbSymbol) & 
                           annotatedPeaks$readCountSum >= readCountSum,]
  txdbSet = unlist(threeUTRsByTranscript(userTxdb, use.names = TRUE))
  txdbSetPc = annotateTxdb(txdbSet, userTxdb, txdbLongTranscripts)
  message('Annotation index created for threeUTR')
  txdbEnsSeq = as_tibble(txdbSetPc)
  txdbEnsSeq$seqnames = paste('chr',txdbEnsSeq$seqnames, sep = '')
  txdbEnsSeq$seqnames = ifelse(grepl('chrMT',txdbEnsSeq$seqnames),'chrM',txdbEnsSeq$seqnames)
  txdbEnsSeqGRange = makeGRangesFromDataFrame(txdbEnsSeq) 
  txdbEnsSequences = getSeq(Hsapiens, txdbEnsSeqGRange)
  txdbSetPc$gcContent = round(letterFrequency(txdbEnsSequences, "GC", as.prob = TRUE),2)
  message('GC content for index calculated')
  threeUtrGcContent = as_tibble(txdbSetPc) %>%
    dplyr::select(symbol, strand, G.C) %>%
    filter(!is.na(symbol)) %>%
    group_by(symbol) %>%
    summarize('gcContent' = median(G.C, na.rm = TRUE))
  threeUtrGcContent$region = 'threeUTR'
  threeUtrGcContent$status = ifelse(threeUtrGcContent$symbol %in% gnSet$txdbSymbol, 'bound', 'unbound')
  
  ####combine the data
  message('Performing final assembly')
  allData = rbind(fiveUtrGcContent, cdsGcContent, threeUtrGcContent)
  allData$region = factor(allData$region, levels = c('fiveUTR','CDS','threeUTR'))
  allData$status = factor(allData$status, levels = c('unbound','bound'))
  return(allData)
}
#######################################################################
#######################################################################

