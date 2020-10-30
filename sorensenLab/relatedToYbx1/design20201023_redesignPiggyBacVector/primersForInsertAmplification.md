# piggyBAC insert amplification

This document contains information about moving constructs from one piggyBAC vector to another. I found that the original vector with the EFS promoter didn't result in enough expression of the genes, so I am going to try moving to a different vector (specifically, CAT#PB510B-1, PB-CMV-MCS-EF1α-Puro, from SBI).

## Primers

In order to move the inserts we have already prepared and sequence-verified from gBlocks, we will just PCR them out and insert them in the MCS from the new vector. We will use NheI on the 5' end, and BstBI or BamHI on the 3' end. 


| Insert | 5' Enzyme | 3' Enzyme |
|------|------|------|
| 3xFLAG eGFP | NheI | BamHI |
| 3xFLAG YBX1 | NheI | BamHI |
| YBX1 3xFLAG | NheI | BamHI |
| APEX2 eGFP | NheI | BstBI |
| APEX2 YBX1 | NheI | BamHI |
| YBX1 APEX2 | NheI | BamHI |
| YBX1 | NheI | BstBI |
| YBX1 K64A | NheI | BstBI |
| YBX1 K81A | NheI | BstBI |


For the forward primer, we can use the same one for all sequences because they all use NheI.

moveInsertPiggyNheI

```
CAGAACACAGGCAAGTTTGTAC
```

For the reverse, we only need a primer with BamHI because the BstBI site is already there. We will just not cut at it when we don't want to.

moveInsertPiggyBamHI

```
CGGATCGTATGGATCCACCTGAGGATCACCACTTT
```

Anneal at +64C in PCR with Q5. Digest with the RE from the above table and anneal into vectors.

