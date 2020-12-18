## Reprocessing CCLE proteomics data

This document details reprocessing of the raw data for the CCLE proteomics performed by the group of Gygi at Harvard.

Quantitative Proteomics of the Cancer Cell Line Encyclopedia
Nusinow et al., Cell, PMID: 31978347, MASSIVE: MSV000085836

### Data download the preparation

The data are hosted on the MASSIVE ftp and I just downloaded them from there into the directory `/projects/ptx_results/OtherDataSets/dataset20201217_ccleProteomicsPmid31978347`. I had to make some changes to some of the raw file names in order to make it a cohesive set in terms of naming convention, but these changes were minor. 

We already have some good scripts for processing the actual data, but there are a total of 38 TMT batches here with 12 fractions in each, so I want it to loop through all of these batches individually. 

Batch names:

```
Prot_1
Prot_2
Prot_3
Prot_4
Prot_5
Prot_6
Prot_7
Prot_8
Prot_9
Prot_10
Prot_11
Prot_12
Prot_13
Prot_14
Prot_15
Prot_16
Prot_17
Prot_18
Prot_19
Prot_20
Prot_21
Prot_22
Prot_23
Prot_24
Prot_25
Prot_26
Prot_27
Prot_28
Prot_29
Prot_30
Prot_31
Prot_32
Prot_33
Prot_34
Prot_35
Prot_36
Prot_37
Prot_38
Prot_39
Prot_40
Prot_41
Prot_42
Prot_01_R2
Prot_01_R3
Prot_02_R2
Prot_02_R3
```
