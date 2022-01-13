#header 

notes 

> output

```
code
```
add a backslash to start a new line


# Syntrophy Samples
Sample_A = &Delta;4H<sub>2</sub>ase pLD28.1 + &Delta;fdh 1&2 pLD28.2 \
Sample_B = &Delta;4H<sub>2</sub>ase pLD28.2 + &Delta;fdh 1&2 pLD28.1


# Barcodes
pLD28.1: GTATTGAA \
pLD28.2: ATACGTGG
```
perl -V
```
>Summary of my perl5 (revision 5 version 16 subversion 3)


Edited MultiCodes.pl file to look for barcodes with a length of 8 rather than 20 bp. Following code is looking at overall sequencing run barcode abundance, not yet dividing by sample using sequencing index. 

```
perl MultiCodes_8.pl -out A_barcodes -index something -preseq GAGGTCTCT -postseq CGT -nPreExpected 0:120 < Aamplicon_R12.fastq
```
>Reads 3615990 Multiplexed 3615990 Usable(8) 424706 (11.7%) unique codes 77 \
Wrong presequence position: 13 reads (0.000%)\
Wrote 77 unique barcodes to A_barcodes.codes \
Wrote number of barcodes seen a certain number of times to A_barcodes.counts; nOnce = 24 \
Wrote 145 off-by-1 pairs (1975 reads, fraction 0.005) to A_barcodes.close \
Wrote A_barcodes.good \
Aside from singletons and off-by-1s, see 3.0 barcodes (99.7% of reads) \

```
perl MultiCodes_8.pl -out B_barcodes -index something -preseq GAGGTCTCT -postseq CGT -nPreExpected 0:120 < Bamplicon_R12.fastq
```
>Reads 2985378 Multiplexed 2985378 Usable(8) 382443 (12.8%) unique codes 69 \
Wrong presequence position: 9 reads (0.000%) \
Wrote 69 unique barcodes to B_barcodes.codes \
Wrote number of barcodes seen a certain number of times to B_barcodes.counts; nOnce = 17 \
Wrote 130 off-by-1 pairs (2308 reads, fraction 0.006) to B_barcodes.close \
Wrote B_barcodes.good \
Aside from singletons and off-by-1s, see 2.0 barcodes (99.6% of reads) \

perl MultiCodes_8.pl -out A_idxbarcodes -primers index.txt -preseq gataactaataggtgaaatgcaGAGGTCTCT -postseq CGT -nPreExpected 0:120 < Aamplicon_R12.fastq
gataactaataggtgaaatgcaGAGGTCTCT


Looking at sequences in Geneious - doesn't seem like primer barcodes were sequenced much \
Use Bowtie2 to map sequences to pLD28 file 
Added in 8 Ns where amplicon barcode located and 8 Ns integrated barcode located

```
bowtie2-build pLD28.txt pLD28
bowtie2 -x pLD28 -q -p 4 --local -1 Aamplicon_S2_R1_001.fastq.gz -2 Aamplicon_S2_R2_001.fastq.gz -S A_amplicon.sam
```
>1807995 reads; of these:\
  1807995 (100.00%) were paired; of these:\
    13299 (0.74%) aligned concordantly 0 times\
    1794316 (99.24%) aligned concordantly exactly 1 time\
    380 (0.02%) aligned concordantly >1 times\
    ----\
    13299 pairs aligned concordantly 0 times; of these:\
      1127 (8.47%) aligned discordantly 1 time\
    ----\
    12172 pairs aligned 0 times concordantly or discordantly; of these:\
      24344 mates make up the pairs; of these:\
        21557 (88.55%) aligned 0 times\
        2766 (11.36%) aligned exactly 1 time\
        21 (0.09%) aligned >1 times\
99.40% overall alignment rate
```
bowtie2 -x pLD28 -q -p 4 --local -1 Bamplicon_S3_R1_001.fastq.gz -2 Bamplicon_S3_R2_001.fastq.gz -S B_amplicon.sam
```
>1492689 (100.00%) were paired; of these:\
    10379 (0.70%) aligned concordantly 0 times\
    1480899 (99.21%) aligned concordantly exactly 1 time\
    1411 (0.09%) aligned concordantly >1 times\
    ----\
    10379 pairs aligned concordantly 0 times; of these:/
      845 (8.14%) aligned discordantly 1 time\
    ----\
    9534 pairs aligned 0 times concordantly or discordantly; of these:\
      19068 mates make up the pairs; of these:\
        16867 (88.46%) aligned 0 times\
        2181 (11.44%) aligned exactly 1 time\
        20 (0.10%) aligned >1 times\
99.44% overall alignment rate\

```
module load fastqc 
mkdir fastqc
fastqc -o fastqc Aamplicon_S2_R1_001.fastq.gz
fastqc -o fastqc Aamplicon_S2_R2_001.fastq.gz
fastqc -o fastqc Bamplicon_S3_R1_001.fastq.gz
fastqc -o fastqc Bamplicon_S3_R2_001.fastq.gz
```
Shortened primer sequences to 4 bp
```
perl MultiCodes_8.pl -out A_idx5barcodes -index idx5 -preseq GCTAgataactaataggtgaaatgcaGAGGTCTCT -postseq CGT -nPreExpected 0:120 < Aamplicon_R12.fastq
```
perl MultiCodes_8.pl -out B_barcodes -primers index.txt -preseq GAGGTCTCT -postseq CGT -nPreExpected 21:23 < Bamplicon_R12.fastq

Not enough coverage of indexing primers - need to redesign indexed amplicon primers and send for sequencing again
