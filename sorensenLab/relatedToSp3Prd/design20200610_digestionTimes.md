# SP3 digestion duration <!-- omit in toc -->

This document describes the experiments related to the duration of digestion related to SP3 PRD.

**Goal:** Determine the length of time needed to achieve complete or near-complete digestion of proteins in a typical workflow for PRD. Ideally, we want to find the shortest time obtainable without additional hardware enhancers (e.g. pressure).

## Experiment design

Want to analyze a set of 8 samples:

* lysate, no digestion
* lysate, 5-minutes digestion
* lysate, 10-minutes digestion
* lysate, 20-minutes digestion
* lysate, 30-minutes digestion
* lysate, 60-minutes digestion
* lysate, 90-minutes digestion
* lysate, 120-minutes digestion

For each of these, the digestion volume will be 60uL and the trypsin ratio will be 1:10 (ug/ug of trypsin to protein). Because we are interested specifically in digestion times here, we will use a FastPrep disrupted lysate. This is to avoid any problems from chromatin. The lysis buffer will be the standard proteomic one that we use in the nuclease protocol. We will do a single SP3 on the samples before digestion, and equally aliquot the elution across the samples (60uL each). To measure the result, we can simply run an SDS-PAGE of a fraction of the digest and look for the disappearance of the intact proteins (take 20uL of sample and add 5uL of loading buffer, load 25uL on the gel). It is also a good idea to extract RNA here, so take 20uL and add 40uL of RNA binding buffer, and proceed with the RNA Clean and Concentrator kit. Perform qPCR on extracted RNA for GAPDH just to look at recovery. Lastly, we can also look at DNA, so take the remaining 20uL of volume, add 5uL of loading buffer and run on a 1% agarose gel.

Based on promising results, we can perform an MS-analysis of replicates of a specific time-point compared to an overnight digest.
