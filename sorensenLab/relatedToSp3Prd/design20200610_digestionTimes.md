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

For each of these, the digestion volume will be 80uL and the trypsin ratio will be 1:10 (ug/ug of trypsin to protein). Because we are interested specifically in digestion times here, we will use a FastPrep disrupted lysate. This is to avoid any problems from chromatin. The lysis buffer will be the standard proteomic one that we use in the nuclease protocol. We will do a single SP3 on the samples before digestion, and equally aliquot the elution across the samples (80uL each). To measure the result, we can simply run an SDS-PAGE of a fraction of the digest and look for the disappearance of the intact proteins (add 20uL of loading buffer, load 25uL on the gel). Based on promising results, we can perform an MS-analysis of replicates of a specific time-point compared to an overnight digest.
