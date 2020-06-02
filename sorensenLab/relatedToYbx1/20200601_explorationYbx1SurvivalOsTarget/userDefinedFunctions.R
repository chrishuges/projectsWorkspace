#######################################################################
#######################################################################
##this file does initial processing for the read count files that
##are read in order for them to be combined later on
#######################################################################
countsFileRead = function(inputDirPath, manifest, ...){
  fileIdentifier = sub('(.*)\\.htseq_counts\\.txt', '\\1', list.files(inputDirPath, pattern = '*\\.htseq_counts\\.txt'))
  inputFile = read_tsv(paste(inputDirPath, 
                             '\\/',
                             list.files(inputDirPath, pattern = '*\\.htseq_counts\\.txt'), 
                             sep = ''),
                       col_names = FALSE) %>%
    mutate(ensg = sub('(.*)\\.[0-9]+$', '\\1', X1)) %>%
    mutate(folderID = sub('.*Target\\/(.*)$', '\\1', inputDirPath)) %>%
    mutate(fileID = sub('(.*)\\.htseq_counts\\.txt', '\\1', list.files(inputDirPath, pattern = '*\\.htseq_counts\\.txt'))) %>%
    mutate(caseID = manifest[which(grepl(fileIdentifier, manifest$file_name)), 'cases'][[1]]$case_id) %>%
    filter(X2 >= 10) %>%
    dplyr::select(folderID, fileID, caseID, ensg, X2)
  colnames(inputFile)[5] = 'counts'
  inputFile$dePath = paste(inputFile$folderID, '/', inputFile$fileID, '.htseq_counts.txt', sep = '')
  return(inputFile)
}
#######################################################################
#######################################################################



