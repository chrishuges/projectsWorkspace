#!/usr/bin/env Rscript

##libraries
suppressMessages(require(tidyverse))
suppressMessages(require(superintronic))
suppressMessages(require(BiocParallel))
suppressMessages(require(ggplot2))
suppressMessages(require(GenomeInfoDb))
suppressMessages(require(GenomicRanges))

##get the input from the command line
commandLineInputs = commandArgs(trailingOnly=TRUE)

#do the superintronic processing
features = plyranges::read_gff('/projects/ptx_analysis/chughes/databases/refgenieIndexes/hg38/gencode_gtf/default/hg38.gtf') %>% 
  collect_parts() %>% 
  dplyr::filter(gene_name == "DLG2")
bamInputFile = commandLineInputs[1]
sampleId = sub('.*starResults\\/(.*)\\.sorted\\.bam$', '\\1', commandLineInputs[1])
design = tibble('sample' = sampleId) %>%
  dplyr::mutate(bam = bamInputFile)
cvg = compute_coverage_long(design, source = 'bam')
cvgStd = keepStandardChromosomes(cvg, pruning.mode = 'coarse')
cvg_over_features = cvgStd %>% 
  dplyr::select(-bam) %>% 
  join_parts(features)

#save the data output
saveRDS(cvg_over_features, paste(commandLineInputs[1],'.coverage.rds', sep = ''))

