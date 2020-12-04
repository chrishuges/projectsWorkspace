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

We can also try in another vector, pLVX-Puro. Here we need to use XhoI on the 5' end, so we need a different primer.

moveInsertPiggyXhoI

```
CAGCTACCTCGAGCAGAACACAGGCAAGTTTGTAC
```

But, we don't actually know where this vector was cut, so this is not so useful. Something we just discovered is that the new SBI vector doesn't include a poly-A sequence after the MCS, so we need to go back and add one without destroying the MCS.

I will target the BGH polyA sequence because it is strong and we have it in the original piggyBAC vector, so we can easily clone it out. 

BGH polyA

```
ctgtgccttctagttgccagccatctgttgtttgcccctcccccgtgccttccttgaccctggaaggtgccactcccactgtcctttcctaataaaatgaggaaattgcatcgcattgtctgagtaggtgtcattctattctggggggtggggtggggcaggacagcaagggggaggattgggaagacaatagcaggcatgctggggatgcggtgggctctatgg
```

bghBamHI_fwd

```
attcgacaGGATCCgcctcgactgtgccttctagttgc
```

bghNotI_rev

```
gatgcggtgggctctatggctcGCGGCCGCtctaacga
tcgttagaGCGGCCGCgagccatagagcccaccgcatc
```

combined insert

```
GGATCCgcctcgactgtgccttctagttgccagccatctgttgtttgcccctcccccgtgccttccttgaccctggaaggtgccactcccactgtcctttcctaataaaatgaggaaattgcatcgcattgtctgagtaggtgtcattctattctggggggtggggtggggcaggacagcaagggggaggattgggaagacaatagcaggcatgctggggatgcggtgggctctatggctcGCGGCCGC
```

NEB says to anneal at 72C. I think around 64C should be ok. 

Lastly, we will make an mCherry version for later use.

mCherryNheI_fwd

```
ACCATCgctagcaccATGgtgagcaagggcgagga
```

mCherryBstBi_rev

```
tagagggcccgtttaaacccgcttcgaaTTACCAGA
TCTGGTAAttcgaagcgggtttaaacgggccctcta
```




