#!/usr/bin/env Rscript

##libraries
suppressMessages(require(Biostrings))
suppressMessages(require(tidyverse))
suppressMessages(require(Peptides))
suppressMessages(require(OrgMassSpecR))

##get the input from the command line
commandLineInputs = commandArgs(trailingOnly=TRUE)

##process the fasta database
fastaDb = readAAStringSet(commandLineInputs[1])
#fastaDb = readAAStringSet('C:/Users/chris/OneDrive/Documents/bccrc/databases/uniprot/uniprotHuman202006.fasta')
fastaIndex = tibble('metadata' = names(fastaDb)) %>%
  mutate(accession = sub(".*[sptr]\\|(.*)\\|.*$", "\\1", metadata)) %>%
  mutate(length = width(fastaDb))
fastaIndex$gene = ifelse(grepl('GN=', fastaIndex$metadata), sub(".*GN=(.*) [PE].*$", "\\1", fastaIndex$metadata), NA)
fastaIndex$species = ifelse(grepl('sapiens', fastaIndex$metadata), 'human', 
                             ifelse(grepl('\\=9606', fastaIndex$metadata), 'human',
                                    ifelse(grepl('musculus', fastaIndex$metadata), 'mouse', 'other')))
fastaIndexFinal = dplyr::select(fastaIndex, accession, gene, species, length)

##get info about tryptic peptides
trypticPeptides = vector()
detectableLength = vector()
for (i in 1:length(fastaDb)){
  aaSeq = as.character(fastaDb[[i]])
  seqDigest = OrgMassSpecR::Digest(aaSeq, enzyme = 'trypsin', missed = 0, custom = list(code = c('X','U','Z','B'), mass = c(50, 60, 70, 80)))
  seqDigest$pepLength = (seqDigest$stop - seqDigest$start) + 1
  seqDigestSub = subset(seqDigest, (seqDigest$pepLength > 5) & (seqDigest$pepLength < 31))
  trypticPeptides = c(trypticPeptides, nrow(seqDigestSub))
  detectableLength = c(detectableLength, sum(seqDigestSub$pepLength, na.rm = TRUE))
  message(paste('Finished ', i, ' proteins.', sep = ''))
}

#now add to the database
fastaIndexFinal$detectablePeptides = trypticPeptides
fastaIndexFinal$detectableLength = detectableLength

##save the index
saveRDS(fastaIndexFinal, paste(commandLineInputs[1],'.annotated.rds',sep = ''))
