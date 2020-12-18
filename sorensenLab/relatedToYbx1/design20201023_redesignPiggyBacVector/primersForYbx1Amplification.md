# ybx1 amplification and mutation

This document contains information related to moving wild-type YBX1 to the SBI vector with the BGH polyA sequence generated as described elsewhere in this folder. 

## Sequences and primers

wild type ybx1

```
ATGAGCAGCGAGGCCGAGACCCAGCAGCCGCCCGCCGCCCCCCCCGCCGCCCCCGCCCTCAGCGCCGCCG
ACACCAAGCCCGGCACTACGGGCAGCGGCGCAGGGAGCGGTGGCCCGGGCGGCCTCACATCGGCGGCGCC
TGCCGGCGGGGACAAGAAGGTCATCGCAACGAAGGTTTTGGGAACAGTAAAATGGTTCAATGTAAGGAAC
GGATATGGTTTCATCAACAGGAATGACACCAAGGAAGATGTATTTGTACACCAGACTGCCATAAAGAAGA
ATAACCCCAGGAAGTACCTTCGCAGTGTAGGAGATGGAGAGACTGTGGAGTTTGATGTTGTTGAAGGAGA
AAAGGGTGCGGAGGCAGCAAATGTTACAGGTCCTGGTGGTGTTCCAGTTCAAGGCAGTAAATATGCAGCA
GACCGTAACCATTATAGACGCTATCCACGTCGTAGGGGTCCTCCACGCAATTACCAGCAAAATTACCAGA
ATAGTGAGAGTGGGGAAAAGAACGAGGGATCGGAGAGTGCTCCCGAAGGCCAGGCCCAACAACGCCGGCC
CTACCGCAGGCGAAGGTTCCCACCTTACTACATGCGGAGACCCTATGGGCGTCGACCACAGTATTCCAAC
CCTCCTGTGCAGGGAGAAGTGATGGAGGGTGCTGACAACCAGGGTGCAGGAGAACAAGGTAGACCAGTGA
GGCAGAATATGTATCGGGGATATAGACCACGATTCCGCAGGGGCCCTCCTCGCCAAAGACAGCCTAGAGA
GGACGGCAATGAAGAAGATAAAGAAAATCAAGGAGATGAGACCCAAGGTCAGCAGCCACCTCAACGTCGG
TACCGCCGCAACTTCAATTACCGACGCAGACGCCCAGAAAACCCTAAACCACAAGATGGCAAAGAGACAA
AAGCAGCCGATCCACCAGCTGAGAATTCGTCCGCTCCCGAGGCTGAGCAGGGCGGGGCTGAGTAA
```

ybx1WtQC_fwd

```
ACCATCgctagcaccATGAGCAGCGAGGCCGAGAC
```

ybx1WtQC_rev (anneal at 64C)

```
GCAGGGCGGGGCTGAGTAAttcgaaTTACCAGA
TCTGGTAAttcgaaTTACTCAGCCCCGCCCTGC
```

combined sequence

```
ACCATCgctagcaccATGAGCAGCGAGGCCGAGACCCAGCAGCCGCCCGCCGCCCCCCCCGCCGCCCCCGCCCTCAGCGCCGCCGACACCAAGCCCGGCACTACGGGCAGCGGCGCAGGGAGCGGTGGCCCGGGCGGCCTCACATCGGCGGCGCCTGCCGGCGGGGACAAGAAGGTCATCGCAACGAAGGTTTTGGGAACAGTAAAATGGTTCAATGTAAGGAACGGATATGGTTTCATCAACAGGAATGACACCAAGGAAGATGTATTTGTACACCAGACTGCCATAAAGAAGAATAACCCCAGGAAGTACCTTCGCAGTGTAGGAGATGGAGAGACTGTGGAGTTTGATGTTGTTGAAGGAGAAAAGGGTGCGGAGGCAGCAAATGTTACAGGTCCTGGTGGTGTTCCAGTTCAAGGCAGTAAATATGCAGCAGACCGTAACCATTATAGACGCTATCCACGTCGTAGGGGTCCTCCACGCAATTACCAGCAAAATTACCAGAATAGTGAGAGTGGGGAAAAGAACGAGGGATCGGAGAGTGCTCCCGAAGGCCAGGCCCAACAACGCCGGCCCTACCGCAGGCGAAGGTTCCCACCTTACTACATGCGGAGACCCTATGGGCGTCGACCACAGTATTCCAACCCTCCTGTGCAGGGAGAAGTGATGGAGGGTGCTGACAACCAGGGTGCAGGAGAACAAGGTAGACCAGTGAGGCAGAATATGTATCGGGGATATAGACCACGATTCCGCAGGGGCCCTCCTCGCCAAAGACAGCCTAGAGAGGACGGCAATGAAGAAGATAAAGAAAATCAAGGAGATGAGACCCAAGGTCAGCAGCCACCTCAACGTCGGTACCGCCGCAACTTCAATTACCGACGCAGACGCCCAGAAAACCCTAAACCACAAGATGGCAAAGAGACAAAAGCAGCCGATCCACCAGCTGAGAATTCGTCCGCTCCCGAGGCTGAGCAGGGCGGGGCTGAGTAAttcgaaTTACCAGA
```

## Mutation work

The sequences below are related to generating mutated versions of YBX1 using the NEB Site Mutagenesis kit (CAT#E0554S). Primer sequenes were designed as directed with the kit using the NEB online tool.

ybx1 k64a

```
ATGAGCAGCGAGGCCGAGACCCAGCAGCCGCCCGCCGCCCCCCCCGCCGCCCCCGCCCTCAGCGCCGCCG
ACACCAAGCCCGGCACTACGGGCAGCGGCGCAGGGAGCGGTGGCCCGGGCGGCCTCACATCGGCGGCGCC
TGCCGGCGGGGACAAGAAGGTCATCGCAACGAAGGTTTTGGGAACAGTAGCCTGGTTCAATGTAAGGAAC
GGATATGGTTTCATCAACAGGAATGACACCAAGGAAGATGTATTTGTACACCAGACTGCCATAAAGAAGA
ATAACCCCAGGAAGTACCTTCGCAGTGTAGGAGATGGAGAGACTGTGGAGTTTGATGTTGTTGAAGGAGA
AAAGGGTGCGGAGGCAGCAAATGTTACAGGTCCTGGTGGTGTTCCAGTTCAAGGCAGTAAATATGCAGCA
GACCGTAACCATTATAGACGCTATCCACGTCGTAGGGGTCCTCCACGCAATTACCAGCAAAATTACCAGA
ATAGTGAGAGTGGGGAAAAGAACGAGGGATCGGAGAGTGCTCCCGAAGGCCAGGCCCAACAACGCCGGCC
CTACCGCAGGCGAAGGTTCCCACCTTACTACATGCGGAGACCCTATGGGCGTCGACCACAGTATTCCAAC
CCTCCTGTGCAGGGAGAAGTGATGGAGGGTGCTGACAACCAGGGTGCAGGAGAACAAGGTAGACCAGTGA
GGCAGAATATGTATCGGGGATATAGACCACGATTCCGCAGGGGCCCTCCTCGCCAAAGACAGCCTAGAGA
GGACGGCAATGAAGAAGATAAAGAAAATCAAGGAGATGAGACCCAAGGTCAGCAGCCACCTCAACGTCGG
TACCGCCGCAACTTCAATTACCGACGCAGACGCCCAGAAAACCCTAAACCACAAGATGGCAAAGAGACAA
AAGCAGCCGATCCACCAGCTGAGAATTCGTCCGCTCCCGAGGCTGAGCAGGGCGGGGCTGAGTAA
```

k64aQC_fwd

```
GGGAACAGTAgccTGGTTCAATGTAAGG
```

k64aQC_rev (anneal at 57C)

```
AAAACCTTCGTTGCGATG
```

ybx1 k81a

```
ATGAGCAGCGAGGCCGAGACCCAGCAGCCGCCCGCCGCCCCCCCCGCCGCCCCCGCCCTCAGCGCCGCCG
ACACCAAGCCCGGCACTACGGGCAGCGGCGCAGGGAGCGGTGGCCCGGGCGGCCTCACATCGGCGGCGCC
TGCCGGCGGGGACAAGAAGGTCATCGCAACGAAGGTTTTGGGAACAGTAAAATGGTTCAATGTAAGGAAC
GGATATGGTTTCATCAACAGGAATGACACCGCCGAAGATGTATTTGTACACCAGACTGCCATAAAGAAGA
ATAACCCCAGGAAGTACCTTCGCAGTGTAGGAGATGGAGAGACTGTGGAGTTTGATGTTGTTGAAGGAGA
AAAGGGTGCGGAGGCAGCAAATGTTACAGGTCCTGGTGGTGTTCCAGTTCAAGGCAGTAAATATGCAGCA
GACCGTAACCATTATAGACGCTATCCACGTCGTAGGGGTCCTCCACGCAATTACCAGCAAAATTACCAGA
ATAGTGAGAGTGGGGAAAAGAACGAGGGATCGGAGAGTGCTCCCGAAGGCCAGGCCCAACAACGCCGGCC
CTACCGCAGGCGAAGGTTCCCACCTTACTACATGCGGAGACCCTATGGGCGTCGACCACAGTATTCCAAC
CCTCCTGTGCAGGGAGAAGTGATGGAGGGTGCTGACAACCAGGGTGCAGGAGAACAAGGTAGACCAGTGA
GGCAGAATATGTATCGGGGATATAGACCACGATTCCGCAGGGGCCCTCCTCGCCAAAGACAGCCTAGAGA
GGACGGCAATGAAGAAGATAAAGAAAATCAAGGAGATGAGACCCAAGGTCAGCAGCCACCTCAACGTCGG
TACCGCCGCAACTTCAATTACCGACGCAGACGCCCAGAAAACCCTAAACCACAAGATGGCAAAGAGACAA
AAGCAGCCGATCCACCAGCTGAGAATTCGTCCGCTCCCGAGGCTGAGCAGGGCGGGGCTGAGTAA
```

k81aQC_fwd

```
GAATGACACCgccGAAGATGTATTTGTACAC
```

k81aQC_rev (anneal at 56C)

```
CTGTTGATGAAACCATATC
```

## Related to sanger sequencing

To sanger sequence and validate these constructs, we will use the primers below. The first sequence is the CMV foward primer, so can be used with any construct.

piggySbiYbx1Sanger
```
CGCAAATGGGCGGTAGGCGTG
GGGACAAGAAGGTCATCGCAAC
CAGCAGACCGTAACCATTATAGAC
AAGGTAGACCAGTGAGGCAG
CTCGACTGTGCCTTCTAGTTGC
```
