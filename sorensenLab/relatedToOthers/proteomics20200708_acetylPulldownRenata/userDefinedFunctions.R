#######################################################################
#######################################################################
##this function goes through the provided quant files and returns
##a parsed object for each. It returns signal to noise ratio for quant
#######################################################################
combineQuantFiles = function(filePath, ...){
  quantData = read_tsv(filePath) %>%
    dplyr::select(MS2ScanNumber, `126SignalToNoise`:`131CSignalToNoise`)
  colnames(quantData) = c('scan','tmt10plex_126','tmt10plex_127N','tmt10plex_127C','tmt10plex_128N',
                  'tmt10plex_128C','tmt10plex_129N','tmt10plex_129C','tmt10plex_130N','tmt10plex_130C','tmt10plex_131N','tmt10plex_131C')
  ##
  fraction = sub('.*ch_(.*)\\.raw_Matrix\\.txt$', '\\1', filePath)
  print(paste('Processing file ', fraction, '.', sep = ''))
  ##
  return(quantData)
}
#######################################################################
#######################################################################



