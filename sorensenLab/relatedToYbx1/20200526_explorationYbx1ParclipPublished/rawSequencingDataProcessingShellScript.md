# Reprocessing the raw sequencing data <!-- omit in toc -->

To reprocess the raw sequencing data, I wrote a shell script that can be found below with some explanations of the different code used in different sections. If you want to use this on your own system, you will need to modify the file paths.

## System and software information

These data were processed on a CentOS Linux release 7.4.1708 system with dual Intel(R) Xeon(R) CPU E5-2690 processors (16 total cores, 32 total threads) and 125GB of ram. For software, this script makes use of:

* BBTools [(version 38.61b)](https://sourceforge.net/projects/bbmap/files/)
* Bowtie [(version 1.2.3)](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.3/)
* Samtools [(version 1.9)](https://sourceforge.net/projects/samtools/files/samtools/)
* R language [(version 3.6.3)](https://www.r-project.org/)
* R package 'bamsignals' (I used version 1.18.0)
* R package 'BSgenome.Hsapiens.UCSC.hg38' (I used version 1.4.1)
* R package 'GenomicFeatures' (I used version 1.38.2)
* R package 'tidyverse' (I used version 1.3.0)
* R package 'wavClusteR' (I used version 2.20.0)

All R packages were installed using BiocManager::install('library').

## Building the index for alignment

To build a Bowtie index for hg38, I used sequence data from Ensembl. I downloaded all the individual chromosome files and extracted them:

```
$wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome*.gz
$gunzip *.gz
```

Then I built the Bowtie index:

```
$/projects/ptx_analysis/chughes/software/bowtie-1.2.3-linux-x86_64/bowtie-build --threads 24 Homo_sapiens.GRCh38.dna.chromosome.1.fa,Homo_sapiens.GRCh38.dna.chromosome.2.fa,Homo_sapiens.GRCh38.dna.chromosome.3.fa,Homo_sapiens.GRCh38.dna.chromosome.4.fa,Homo_sapiens.GRCh38.dna.chromosome.5.fa,Homo_sapiens.GRCh38.dna.chromosome.6.fa,Homo_sapiens.GRCh38.dna.chromosome.7.fa,Homo_sapiens.GRCh38.dna.chromosome.8.fa,Homo_sapiens.GRCh38.dna.chromosome.9.fa,Homo_sapiens.GRCh38.dna.chromosome.10.fa,Homo_sapiens.GRCh38.dna.chromosome.11.fa,Homo_sapiens.GRCh38.dna.chromosome.12.fa,Homo_sapiens.GRCh38.dna.chromosome.13.fa,Homo_sapiens.GRCh38.dna.chromosome.14.fa,Homo_sapiens.GRCh38.dna.chromosome.15.fa,Homo_sapiens.GRCh38.dna.chromosome.16.fa,Homo_sapiens.GRCh38.dna.chromosome.17.fa,Homo_sapiens.GRCh38.dna.chromosome.18.fa,Homo_sapiens.GRCh38.dna.chromosome.19.fa,Homo_sapiens.GRCh38.dna.chromosome.20.fa,Homo_sapiens.GRCh38.dna.chromosome.21.fa,Homo_sapiens.GRCh38.dna.chromosome.22.fa,Homo_sapiens.GRCh38.dna.chromosome.X.fa,Homo_sapiens.GRCh38.dna.chromosome.Y.fa,Homo_sapiens.GRCh38.dna.chromosome.MT.fa ./homoSapiensEnsemblGRCh38_rel100
```

Then I copied the created index files over to the Bowtie index directory to finish the install:

```
$scp homoSapiensEnsemblGRCh38_rel100.*.ebwt /projects/ptx_analysis/chughes/software/bowtie-1.2.3-linux-x86_64/indexes/
```

## Shell script

The shell script detailed below was run using the command:

```
$ ./parclipProcessing.sh |& tee data-processing.txt
```

Create a shell script and insert the below text.

```
#! /bin/bash

for i in SRR9623531 SRR9623532 SRR9623533 SRR9623534 SRR9623535 SRR9623536
do

echo $i
    
eval /projects/ptx_analysis/chughes/software/BBTools/bbmap/bbduk.sh in=${i}.fastq.gz out=${i}.trimmed.fastq.gz stats=${i}.trimming_stats.txt literal=AGATCGGAAGAGCACACGTCT,AAAAAAAAAAAA k=12 ktrim=r qtrim=rl ordered=t minlength=20 maxlength=40 tp=4

eval /projects/ptx_analysis/chughes/software/BBTools/bbmap/bbduk.sh in=${i}.trimmed.fastq.gz out=${i}.trimmed.rrna.fastq.gz stats=${i}.rRNA_stats.txt k=18 ordered=t ref=/projects/ptx_analysis/chughes/databases/genbank_rRNA/genbank_aug2019_rRNA.fasta forcetrimleft=3 forcetrimright=28

eval /projects/ptx_analysis/chughes/software/bowtie-1.2.3-linux-x86_64/bowtie homoSapiensEnsemblGRCh38_rel100 -S -q -v 2 -m 10 --no-unal --best --strata /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/${i}.trimmed.rrna.fastq.gz /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/${i}.trimmed.rrna.bowtie.sam

eval /projects/ptx_analysis/chughes/software/samtools-1.9/samtools view -S -b /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/${i}.trimmed.rrna.bowtie.sam > /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/${i}.trimmed.rrna.bowtie.bam

eval /projects/ptx_analysis/chughes/software/samtools-1.9/samtools view -b -F 4 /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/${i}.trimmed.rrna.bowtie.bam > /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/${i}.trimmed.rrna.bowtie.mapped.bam

eval /projects/ptx_analysis/chughes/software/samtools-1.9/samtools view -h /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/${i}.trimmed.rrna.bowtie.mapped.bam > /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/${i}.trimmed.rrna.bowtie.mapped.sam

eval /projects/ptx_analysis/chughes/software/samtools-1.9/samtools sort /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/${i}.trimmed.rrna.bowtie.mapped.bam > /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/${i}.trimmed.rrna.bowtie.mapped.sorted.bam

eval /projects/ptx_analysis/chughes/software/samtools-1.9/samtools index /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/${i}.trimmed.rrna.bowtie.mapped.sorted.bam > /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/${i}.trimmed.rrna.bowtie.mapped.sorted.bai

eval Rscript /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/wavClusterProcessing.R /projects/ptx_results/Sequencing/published-studies/201908_Chen_NatCellBiol_PMID31358969/parclip/ ${i}

done
```

There is also an R script that goes along with this shell script.

```
#!/usr/bin/env Rscript

#Load in the libraries
###################################################################################################################################################################

suppressMessages(require(tidyverse))
suppressMessages(require(wavClusteR))
suppressMessages(require(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(require(GenomicFeatures))
suppressMessages(require(bamsignals))

#######get the input directory
command_line_inputs = commandArgs(trailingOnly=TRUE)

###################################################################################################################################################################

base_dir = command_line_inputs[1]
filename = paste(base_dir, command_line_inputs[2], '.trimmed.rrna.bowtie.mapped.sorted.bam', sep = '')
message(base_dir)
message(filename)

##
message('reading bam file')
Bam <- readSortedBam(filename = filename)

message('preparing count table')
countTable <- getAllSub( Bam, minCov = 20 )
#head( countTable )

pdf(paste(base_dir, command_line_inputs[2], '.subs1.pdf', sep = ''))
plotSubstitutions( countTable, highlight = "TC" )
dev.off()

message('preparing model')
model <- fitMixtureModel(countTable, substitution = "TC")
data(model)
str(model)
(support <- getExpInterval( model, bayes = TRUE ) )

pdf(paste(base_dir, command_line_inputs[2], '.subs2.pdf', sep = ''))
plotSubstitutions( countTable, highlight = "TC", model )
dev.off()

message('getting high confidence substitutions')
highConfSub <- getHighConfSub( countTable, support = support, substitution = "TC" )
#head( highConfSub )

message('getting coverage')
coverage <- coverage( Bam )

message('getting clusters')
clusters <- getClusters( highConfSub = highConfSub, coverage = coverage, sortedBam = Bam, method = "mrn", cores = 12 )

message('filtering clusters')
wavclusters <- filterClusters( clusters = clusters,
                               highConfSub = highConfSub,
                               coverage = coverage,
                               model = model,
                               genome = Hsapiens,
                               refBase = "T",
                               minWidth = 12)

message('exporting high confidence substitutions')
exportHighConfSub( highConfSub = highConfSub,
                   filename = paste(base_dir, command_line_inputs[2], '.hcTC.bed', sep = ''),
                   trackname = "hcTC",                                                                                                                         description = "hcTC" )
saveRDS(highConfSub, paste(base_dir, command_line_inputs[2], '.hcTC.rds', sep = ''))


message('exporting clusters')
exportClusters( clusters = wavclusters,
                filename = paste(base_dir, command_line_inputs[2], '.wavClusters.bed', sep = ''),
                trackname = "wavClusters",
                description = "wavClusters" )
saveRDS(wavclusters, paste(base_dir, command_line_inputs[2], '.wavClusters.rds', sep = ''))


message('exporting coverage')
exportCoverage( coverage = coverage,
                filename = paste(base_dir, command_line_inputs[2], '.coverage.bigWig', sep = ''))
saveRDS(coverage, paste(base_dir, command_line_inputs[2], '.coverage.rds', sep = ''))

message('done processing')
```

### Some notes about the processing

* Although this sequencing was done in paired end format as stated in the manuscript, only the first read is provided in the data repository. According to the authors, only the first read was used in their analysis for the manuscript.
* I had to be really strict in the trimming of the reads. Even after adapter removal, there was still these long polyA runs that didn't look real and really had a negative impact on the overall alignment, so I trimmed them off for the most part. There were also many sequences where the adapter sequence was messed up, likely due to bad calls in those regions because of all the polyA, meaning I had to discard sequences where I couldn't detect the adapter. So, for a 150bp PE read set, these runs didn't yield a ton of sequencing data.
* Bowtie uses pretty standard parameters, allowing for 2 mismatches because of the expected T>C conversion due to the UV-crosslinking.
* The samtools processing is a bit convoluted, but is necessary to remove non-mapped, non-unique reads as wavClusteR will throw an error if you don't.
* In wavClusteR, I opted for a minimum coverage of 20 for T>C conversion sites. Because there was such a different depth of read coverage between the files, sometimes I could set this really high, and other times really low and be happy. However, I opted to set it somewhat low and to use replicates to filter out lower confidence sites (unfortunately they only had 2 replicates, but it is better than nothing).
