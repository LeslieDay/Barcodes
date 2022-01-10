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

