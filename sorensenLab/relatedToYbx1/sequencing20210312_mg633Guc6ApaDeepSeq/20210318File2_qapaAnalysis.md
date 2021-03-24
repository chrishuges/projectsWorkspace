## Processing with QAPA

The purpose of these sequencing runs was to look at alternative polyadenylation. To do this, we will use a tool called [QAPA](https://github.com/morrislab/qapa). We are pretty much going to follow their default workflow that they describe there. Getting it going on Python took a little bit of work, but we got there in the end.

The first you need is the annotation file. They give a good description of how to build one and also provide some pre-compiled files. For this first analysis, I think I will use the precompiled files.

First we need to process with our genome fasta.

```shell
qapa fasta -f /projects/ptx_analysis/chughes/databases/refgenieIndexes/alias/hg38/fasta/default/hg38.fa /projects/ptx_analysis/chughes/databases/qapa/qapa_3utrs.gencode_V31.hg38.bed /projects/ptx_analysis/chughes/databases/qapa/output_sequences.fa
```

Then generate the salmon index.

```
/projects/ptx_analysis/chughes/software/salmon-1.4/bin/salmon index -t ./output_sequences.fa -i utr_library
```

Then we need to run Salmon itself.

```shell
#!/bin/bash
salmonLocation="/projects/ptx_analysis/chughes/software/salmon-1.4/bin/salmon"
indexLocation="/projects/ptx_analysis/chughes/databases/qapa/utr_library/"
rawDataOutputDirectory="/projects/ptx_results/Sequencing/2021/sequencing20210312_mg633Ybx1Guc6ApaDeepSeq/"

###########################################
for i in ATCACG CGATGT
do
  echo $i
  ##
  salmonCall="$salmonLocation quant -i $indexLocation -l A -1 ${rawDataOutputDirectory}${i}_1.clean.fastq.gz -2 ${rawDataOutputDirectory}${i}_2.clean.fastq.gz -o ${rawDataOutputDirectory}${i}qapa_quant"
  eval $salmonCall
done
```

