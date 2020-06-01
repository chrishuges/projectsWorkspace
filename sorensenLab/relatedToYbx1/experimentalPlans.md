# Experimental plans <!-- omit in toc -->

This document details experiments we are interested in performing related to the YBX1 project.

<hr style="height:6pt; visibility:hidden;" />

## Quick links <!-- omit in toc -->

- [General information](#general-information)
- [1. Phenotypes of hetKD lines](#1-phenotypes-of-hetkd-lines)
- [2. Global expression profiling](#2-global-expression-profiling)
- [3. CLIPseq profiling](#3-clipseq-profiling)
- [4. Polysome profiling atlas](#4-polysome-profiling-atlas)

<hr style="height:6pt; visibility:hidden;" />

<span id="general-information"></span>

## General information

First, just to note, these are the experiments I would like to start with because I think they will fill in gaps in our current understanding and data. We can branch out from here, but I think these are the most important first steps.

All of the below experiments are going to be performed in MG63.3 cells, which are a highly-metastatic variant of MG63 cells. Currently, we have multiple CRISPR clones where multiple copies of YBX1 has been knocked out (hetKD). We don't have a complete knockout because these cells do not survive. There is some literature suggesting that when YBX1 is knocked out that YBX3 will come in and compensate, but it appears this does not happen with MG63.3 cells. Our hetKD system is a bit nicer because the cells have not 'adapted' to the reduction of YBX1. It is not completely gone, so they are able to survive, but they are far more sensitive to certain stimuli than wild type cells and have a much slower growth rate.

To link the hetKD with function, we will use YBX1 mutants. From some studies of YBX1 in complex with RNA's it appears the cold shock domain is the most important region. Various residues are important in this region, but for our work, K64 and K81 will be the main targets due to potential links with acetylation. We have engineered constructs for both of these variants that can be introduced into the hetKD cells to measure their effects.

<span id="1-phenotypes-of-hetkd-lines"></span>

## 1. Phenotypes of hetKD lines

One thing we are missing right now is a nice clear cut analysis of the phenotypes of the knockdown lines (mostly because we only just generated them). My preliminary data suggests the hetKD lines are far more sensitive to oxidative stress induced by tBHP. But, we don't have anything specifically linking this to YBX1 function other than the fact that these are hetKD lines. The main experiment I would propose here initially would be Incucyte analysis for the two independent hetKD clones. This growth experiment will also be done with two tBHP concentrations (low and high). The logic here is that it will allow us to separate out potential SG mediated effects that will likely occur in high but not low concentrations.

The main questions here are:

* Do the hetKD lines behave differently in standard conditions?
* Do the hetKD lines respond differently to stress?
* Can the hetKD lines be rescued with different mutants of YBX1?

Samples will be:

1. wild type cells
2. hetKD
3. hetKD + wild type YB1
4. hetKD + K64A YB1
5. hetKD + K81A YB1
6. wild type cells + tBHP
7. hetKD + tBHP
8. hetKD + wild type YB1 + tBHP
9. hetKD + K64A YB1 + tBHP
10. hetKD + K81A YB1 + tBHP

Depending on the result of this, we move on to doing oxidative stress (flow cytometry), GSH levels, and perhaps GLO1 assays by metabolomics.

<span id="2-global-expression-profilin"></span>

## 2. Global expression profiling

Depending on the results of the phenotype analysis, profiling the global changes in these hetKD cells would be valuable. This objective is fairly self-explanatory, but it is valuable to acquire a 'base' dataset of RNA and protein expression for the hetKD lines that we can use as a reference for different observations we make in more targeted experiments. 

The main question to answer here would be:

* When YBX1 is reduced, what is changed in the cells at the RNA and protein levels?

Samples here would be:

1. wild-type x3 replicates
2. hetKD clone 1 x3 replicates
3. hetKD clone 2 x3 replicates

Matched sets will be sent for proteomics and RNAseq analysis. I am not convinced that we will see a lot of stuff here, but it is something that absolutely needs to be done to satisfy reviewers.

<span id="3-clipseq-profiling"></span>

## 3. CLIPseq profiling

Based on our preliminary data, I anticipate that we will see that the mutants demonstrate different abilities to recover hetKD cells (with or without stress). I would propose that this is due to YBX1 regulating different transcription events.

The main questions to answer here are:

* Which RNA's does YBX1 bind to (and where)?
* How is this modulated by YBX1 mutants?
* What role does stress play in the YBX1 RNA interactome?
  
So, we would do a CLIPseq analysis with the samples:

1. Base RNA interaction set
   1. wild type cells x3 replicates
   2. hetKD clone A x3 replicates
   3. hetKD clone B x3 replicates
2. hetKD RNA interaction set (pick one clone)
   1. hetKD x3 replicates
   2. hetKD + YB1 x3 replicates
   3. hetKD + K64A YB1 x3 replicates
   4. hetKD + K81A YB1 x3 replicates

Depending on the results of these first two runs, we would do additional experiments where we are using stress conditions. I don't want to do this right out of the gate because it ends up being a huge number of sequencing runs. If we can subset our sample set down based on our data, we can save a lot of time/money/effort. I envision the next set would be something like:

1. Follow-up set 1 (pick one clone)
   1. hetKD + wild type YB1 x3 replicates
   2. hetKD + wild type YB1 + tBHP x3 replicates
   3. hetKD + wild type YB1 + MS275 x3 replicates
   4. hetKD + wild type YB1 + tBHP + MS275 x3 replicates
2. Follow-up set 2 (pick one clone)
   1. hetKD + K81 YB1 x3 replicates
   2. hetKD + K81 YB1 + tBHP x3 replicates
   3. hetKD + K81 YB1 + MS275 x3 replicates
   4. hetKD + K81 YB1 + tBHP + MS275 x3 replicates

The assumption for the second follow-up set will be that K81A doesn't completely wipe out RNA binding, but rather just alters it (K64A most likely will wipe it out). From here, we will have a lot of information related to RNA binding and expression, that we can use to leverage our way into RNA modification studies, if we see the need. Or, even targeted modifications of validated RNA binding sites (such as on EIF4A1) to see what impact this has.

<span id="4-polysome-profiling-atlas"></span>

## 4. Polysome profiling atlas

One thing we have already done is to profile protein abundance across a fractionated sucrose gradient. Here we can see that YBX1 is present at all times and it does not disappear when polysome samples are EDTA treated (dissociates polysomes), indicating that we are likely seeing multiple different YBX1 containing complexes. Unfortunately, it is a bit hard to map this out with our current data because despite our sensitivity being much better than current studies, our resolution is a bit lower (resolution depends on how many fractions you collect). The ability to correlate elution profiles and thus establish potential interactions of proteins is really dependent on resolution. I originally did these experiments at a lower resolution because I did not anticipate we would see what we do for YBX1. Doing this with our hetKD lines would additionally allow us to answer some really specific questions about YBX1 and translation as we may be able to assign translation-specific interactions.

The main questions I want to answer here are:

* With a high-resolution profile of YBX1 protein abundance (that includes proper controls), can we link YBX1 abundance to other proteins (e.g. ribosomes)?
* When we remove polysomes (with EDTA treatment), which of these potential interactions drop away and which remain?
* Can we link those proteins that remain to a secondary (e.g. active translation independent) function for YBX1?
* Through integration with our protein interaction data, which of YBX1's interactions are translation-dependent?
* Which RNA's are present in the same regions as YBX1?
* Through integration with published CLIPseq data, can we establish links between YBX1 mRNA binding and active translation?

For this, I would propose we do:

1. Polysome proteome profiles (pick one hetKD clone)
   1. hetKD x3 replicates
   2. hetKD + wild type YB1 x3 replicates

From here, we would follow up to profile specific RNA's (targeted via our CLIPseq knowledge) by qPCR to see if when the distribution of YBX1 changes, their profile does as well.
