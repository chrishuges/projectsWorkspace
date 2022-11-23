## Running MAJIQ

This document describes the process of running [MAJIQ](https://biociphers.bitbucket.io/majiq-docs-academic/getting-started-guide/quick-overview.html). I am following the instructions at that link. 

### Data pipeline

First we need to make our annotation file.

```shell
majiq build  -c settings_file.ini /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gff3 -o /home/chughes/databases/projectEwsDlg2/majiqFiles  -j 8

```

I actually ended up following [this document](https://biociphers.bitbucket.io/majiq-docs-academic/gallery/heterogen-vignette.html) more closely. I created a settings file as below.

```shell
# global parameters
[info]
# where do we find the input bam files?
bamdirs=/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/results
# this is used by voila to produce links to UCSC genome browser
genome=GRCh38
# if our experiments are stranded, we can get better results by specifying so
# (a little bit more work if we have mixed strandedness)
strandness=reverse

# we divide input experiments (bam or sj) into "build groups"
# we find these files in bamdirs or sjdirs
[experiments]
# EWS-FLI1 high
EWSFL1high=ATCACG_setA.sorted.bam,ATCACG_setB.sorted.bam,ATCACG_setC.sorted.bam

# EWS-FLI1 low
EWSFL1low=CGATGT_setA.sorted.bam,CGATGT_setB.sorted.bam,CGATGT_setC.sorted.bam
```

Then I ran MAJIQ with the command below.

```shell
 majiq build /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gff3 -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build -c /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/settingsNew.ini
```

This is the output.

```shell
majiq build /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gff3 -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build -c /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/settingsNew.ini
2022-11-22 03:06:07,972 (PID:1913756) - INFO - Majiq Build v2.4.dev3+g85d07819
2022-11-22 03:06:07,972 (PID:1913756) - INFO - Command: /home/chughes/virtualPython3810/bin/majiq build /home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gff3 -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build -c /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/settingsNew.ini
2022-11-22 03:06:07,973 (PID:1913756) - INFO - Parsing GFF3
2022-11-22 03:07:13,748 (PID:1913756) - INFO - Reading bamfiles
2022-11-22 03:07:13,760 (PID:1913756) - INFO - Group EWSFL1high, number of experiments: 3, minexperiments: 2
2022-11-22 03:07:13,776 (PID:1913756) - INFO - Reading bam file /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/results/ATCACG_setA.sorted.bam
2022-11-22 03:08:16,397 (PID:1913756) - INFO - Detect Intron retention ATCACG_setA.sorted
2022-11-22 03:09:20,963 (PID:1913756) - INFO - Done Reading file ATCACG_setA.sorted
2022-11-22 03:09:21,769 (PID:1913756) - INFO - Reading bam file /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/results/ATCACG_setB.sorted.bam
2022-11-22 03:10:12,465 (PID:1913756) - INFO - Detect Intron retention ATCACG_setB.sorted
2022-11-22 03:11:07,140 (PID:1913756) - INFO - Done Reading file ATCACG_setB.sorted
2022-11-22 03:11:07,901 (PID:1913756) - INFO - Reading bam file /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/results/ATCACG_setC.sorted.bam
2022-11-22 03:12:00,205 (PID:1913756) - INFO - Detect Intron retention ATCACG_setC.sorted
2022-11-22 03:12:54,795 (PID:1913756) - INFO - Done Reading file ATCACG_setC.sorted
2022-11-22 03:12:55,541 (PID:1913756) - INFO - Group EWSFL1low, number of experiments: 3, minexperiments: 2
2022-11-22 03:12:55,542 (PID:1913756) - INFO - Reading bam file /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/results/CGATGT_setA.sorted.bam
2022-11-22 03:13:52,683 (PID:1913756) - INFO - Detect Intron retention CGATGT_setA.sorted
2022-11-22 03:14:47,786 (PID:1913756) - INFO - Done Reading file CGATGT_setA.sorted
2022-11-22 03:14:48,548 (PID:1913756) - INFO - Reading bam file /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/results/CGATGT_setB.sorted.bam
2022-11-22 03:15:40,303 (PID:1913756) - INFO - Detect Intron retention CGATGT_setB.sorted
2022-11-22 03:16:32,933 (PID:1913756) - INFO - Done Reading file CGATGT_setB.sorted
2022-11-22 03:16:33,711 (PID:1913756) - INFO - Reading bam file /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/results/CGATGT_setC.sorted.bam
2022-11-22 03:17:36,756 (PID:1913756) - INFO - Detect Intron retention CGATGT_setC.sorted
2022-11-22 03:18:38,949 (PID:1913756) - INFO - Done Reading file CGATGT_setC.sorted
2022-11-22 03:18:39,988 (PID:1913756) - INFO - Detecting LSVs ngenes: 60649
2022-11-22 03:19:03,499 (PID:1913756) - INFO - 112728 LSV found
2022-11-22 03:19:03,559 (PID:1913756) - INFO - DUMP file b'ATCACG_setA.sorted'
2022-11-22 03:19:11,602 (PID:1913756) - INFO - Create majiq file
2022-11-22 03:19:14,093 (PID:1913756) - INFO - Dump majiq file
2022-11-22 03:19:14,915 (PID:1913756) - INFO - ATCACG_setA.sorted: 85767 LSVs
2022-11-22 03:19:15,018 (PID:1913756) - INFO - DUMP file b'ATCACG_setB.sorted'
2022-11-22 03:19:22,760 (PID:1913756) - INFO - Create majiq file
2022-11-22 03:19:25,375 (PID:1913756) - INFO - Dump majiq file
2022-11-22 03:19:26,149 (PID:1913756) - INFO - ATCACG_setB.sorted: 85767 LSVs
2022-11-22 03:19:26,251 (PID:1913756) - INFO - DUMP file b'ATCACG_setC.sorted'
2022-11-22 03:19:33,874 (PID:1913756) - INFO - Create majiq file
2022-11-22 03:19:36,518 (PID:1913756) - INFO - Dump majiq file
2022-11-22 03:19:37,378 (PID:1913756) - INFO - ATCACG_setC.sorted: 85767 LSVs
2022-11-22 03:19:37,479 (PID:1913756) - INFO - DUMP file b'CGATGT_setA.sorted'
2022-11-22 03:19:45,430 (PID:1913756) - INFO - Create majiq file
2022-11-22 03:19:47,995 (PID:1913756) - INFO - Dump majiq file
2022-11-22 03:19:48,811 (PID:1913756) - INFO - CGATGT_setA.sorted: 85767 LSVs
2022-11-22 03:19:48,913 (PID:1913756) - INFO - DUMP file b'CGATGT_setB.sorted'
2022-11-22 03:19:57,119 (PID:1913756) - INFO - Create majiq file
2022-11-22 03:19:59,658 (PID:1913756) - INFO - Dump majiq file
2022-11-22 03:20:00,463 (PID:1913756) - INFO - CGATGT_setB.sorted: 85767 LSVs
2022-11-22 03:20:00,567 (PID:1913756) - INFO - DUMP file b'CGATGT_setC.sorted'
2022-11-22 03:20:08,819 (PID:1913756) - INFO - Create majiq file
2022-11-22 03:20:11,458 (PID:1913756) - INFO - Dump majiq file
2022-11-22 03:20:12,270 (PID:1913756) - INFO - CGATGT_setC.sorted: 85767 LSVs
2022-11-22 03:20:13,302 (PID:1913756) - INFO - MAJIQ Builder is ended successfully!
```

Run the quantification.

```shell
majiq deltapsi -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/deltapsi -grp1 /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setA.sorted.majiq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setB.sorted.majiq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setC.sorted.majiq -grp2 /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setA.sorted.majiq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setB.sorted.majiq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setC.sorted.majiq  -n EWSFL1low EWSFL1high
```

This is the output.

```shell
majiq deltapsi -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/deltapsi -grp1 /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setA.sorted.majiq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setB.sorted.majiq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setC.sorted.majiq -grp2 /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setA.sorted.majiq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setB.sorted.majiq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setC.sorted.majiq  -n EWSFL1low EWSFL1high
2022-11-22 03:20:39,592 (PID:1914268) - INFO - Majiq deltapsi v2.4.dev3+g85d07819
2022-11-22 03:20:39,592 (PID:1914268) - INFO - Command: /home/chughes/virtualPython3810/bin/majiq deltapsi -o /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/deltapsi -grp1 /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setA.sorted.majiq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setB.sorted.majiq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setC.sorted.majiq -grp2 /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setA.sorted.majiq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setB.sorted.majiq /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setC.sorted.majiq -n EWSFL1low EWSFL1high
2022-11-22 03:20:39,593 (PID:1914268) - INFO - GROUP1: ['/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setA.sorted.majiq', '/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setB.sorted.majiq', '/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setC.sorted.majiq']
2022-11-22 03:20:39,593 (PID:1914268) - INFO - GROUP2: ['/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setA.sorted.majiq', '/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setB.sorted.majiq', '/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setC.sorted.majiq']
2022-11-22 03:20:39,606 (PID:1914268) - INFO - Parsing file: /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setA.sorted.majiq
2022-11-22 03:20:41,508 (PID:1914268) - INFO - Parsing file: /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setB.sorted.majiq
2022-11-22 03:20:43,708 (PID:1914268) - INFO - Parsing file: /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/CGATGT_setC.sorted.majiq
2022-11-22 03:20:46,708 (PID:1914268) - INFO - Group EWSFL1low: 74603 LSVs
2022-11-22 03:20:46,709 (PID:1914268) - INFO - Parsing file: /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setA.sorted.majiq
2022-11-22 03:20:48,753 (PID:1914268) - INFO - Parsing file: /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setB.sorted.majiq
2022-11-22 03:20:51,218 (PID:1914268) - INFO - Parsing file: /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/ATCACG_setC.sorted.majiq
2022-11-22 03:20:54,155 (PID:1914268) - INFO - Group EWSFL1high: 74603 LSVs
2022-11-22 03:20:54,176 (PID:1914268) - INFO - Number quantifiable LSVs: 69287
2022-11-22 03:20:54,176 (PID:1914268) - INFO - Calculating prior matrix...
BETA PARAMS: a: 25.9589 b: 25.9589 p: 51.9179 vari: 0.0047243 mean: 0.5
BETA PARAMS: a: 1281.49 b: 1281.49 p: 2562.98 vari: 9.75045e-05 mean: 0.5
BETA PARAMS: a: 25.4423 b: 25.4423 p: 50.8845 vari: 0.00481839 mean: 0.5
BETA PARAMS: a: 1589.04 b: 1589.04 p: 3178.08 vari: 7.8639e-05 mean: 0.5
2022-11-22 03:20:57,883 (PID:1914268) - INFO - Start Computing DeltaPSI
2022-11-22 03:21:30,653 (PID:1914268) - INFO - Computation done, saving results....
2022-11-22 03:24:08,252 (PID:1914268) - INFO - DeltaPSI calculation for EWSFL1low-EWSFL1high ended successfully! Result can be found at /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/deltapsi
```

Run the output command.

```shell
voila tsv -f /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/voila_output.tsv /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq
```

This is the output.

```shell
voila tsv -f /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/voila_output.tsv /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq
2022-11-22 03:24:42,564 (PID:1914388) - INFO - Command: /home/chughes/virtualPython3810/bin/voila tsv -f /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/voila_output.tsv /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq
2022-11-22 03:24:42,564 (PID:1914388) - INFO - Voila v2.4.dev3+g85d07819
2022-11-22 03:24:42,564 (PID:1914388) - INFO - config file: /tmp/tmpk61ob8ig
2022-11-22 03:24:42,599 (PID:1914388) - INFO - deltapsi TSV
2022-11-22 03:24:42,619 (PID:1914388) - INFO - Creating Tab-delimited output file
2022-11-22 03:24:59,871 (PID:1914388) - INFO - Delimited output file successfully created in: /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/voila_output.tsv
2022-11-22 03:24:59,881 (PID:1914388) - INFO - Duration: 0:00:17.265277
```

What if we output modulize.

```shell
voila modulize -d /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/modulize /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/build/splicegraph.sql /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput/majiq/deltapsi/EWSFL1low-EWSFL1high.deltapsi.voila -j8
```

Here is the output.

```shell

```