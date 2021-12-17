# Barcoded_transposons
### hello world 
### The following is script to access MSI, request node time, map amplicon sequences to the barcode plasmid, pull mapped sequences, and extract/count barcodes 
# Running MultiCodes.pl script on JJ and S2 library
### unzip .gz files and concatenate R1 and R2 into R12 file
```C++
gzip -c -d JJ_amplicon_S198_R1_001.fastq.gz > JJ_amplicon_S198_R1_001.fastq
gzip -c -d JJ_amplicon_S198_R2_001.fastq.gz > JJ_amplicon_S198_R2_001.fastq
gzip -c -d S2_amplicon_S199_R2_001.fastq.gz > S2_amplicon_S199_R1_001.fastq
gzip -c -d S2_amplicon_S199_R1_001.fastq.gz > S2_amplicon_S199_R2_001.fastq

cat S2_amplicon_S199_R1_001.fastq S2_amplicon_S199_R2_001.fastq > S2_amplicon_R12.fastq
cat JJ_amplicon_S198_R1_001.fastq JJ_amplicon_S198_R2_001.fastq > JJ_amplicon_R12.fastq
```
### Use MultiCodes.pl code to extract barcodes from R2 and R12 files
```perl
perl MultiCodes.pl -out JJ_barcodes -index something -preseq GAGGTCTCT -postseq CGT -nPreExpected 0:120 < JJ_amplicon_S198_R2_001.fastq
perl MultiCodes.pl -out S2_barcodes -index something -preseq GAGGTCTCT -postseq CGT -nPreExpected 0:120 < S2_amplicon_S199_R2_001.fastq
perl MultiCodes.pl -out S2_R12_barcodes -index something -preseq GAGGTCTCT -postseq CGT -nPreExpected 0:120 < S2_amplicon_R12.fastq
perl MultiCodes.pl -out JJ_R12_barcodes -index something -preseq GAGGTCTCT -postseq CGT -nPreExpected 0:120 < JJ_amplicon_R12.fastq
```

### Make a bowtie2 index of the barcode plasmid and use bowtie2 local alignment to map S2 reads 
#### updated barcode plasmid to add random primer overhang included to reach 400 bp length requirement for MIGs sequencing 
```
srun -N 1 --ntasks-per-node=4  --mem-per-cpu=6gb -t 4:00:00 -p interactive --pty bash
module load bowtie2
cd barcoded_tn
bowtie2-build barcode2_plasmid.fasta plasmid2_index
bowtie2 -x plasmid2_index -q -p 4 --local -1 S2_amplicon_S199_R1_001.fastq.gz -2 S2_amplicon_S199_R2_001.fastq.gz -S S2barcode_local.sam
```
### results log of bowtie2 alignment 
```
1974629 reads; of these:
  1974629 (100.00%) were paired; of these:
    219328 (11.11%) aligned concordantly 0 times
    1677764 (84.97%) aligned concordantly exactly 1 time
    77537 (3.93%) aligned concordantly >1 times
    ----
    219328 pairs aligned concordantly 0 times; of these:
      5020 (2.29%) aligned discordantly 1 time
    ----
    214308 pairs aligned 0 times concordantly or discordantly; of these:
      428616 mates make up the pairs; of these:
        421342 (98.30%) aligned 0 times
        6540 (1.53%) aligned exactly 1 time
        734 (0.17%) aligned >1 times
89.33% overall alignment rate

```
### convert .sam to .bam, sort, extract reads mapping to barcode region (262-281)

```
module load samtools 
# .sam to .bam
samtools view -S -b S2barcode_local.sam > S2barcode_local.bam

# sort bam file
samtools sort S2barcode_local.bam -o S2barcode_sorted.bam

# index the bam file
samtools index S2barcode_sorted.bam

#export reads mapping to the barcode
samtools view S2barcode_sorted.bam "barcode2_plasmid:262-281" > S2sorted_barcode_only.bam

```
### export fastq from sorted.bam using geneious 
### concatenate into single file to run barcode perl script on 
```
cat S2sorted_barcode_only_unpaired.fastq.gz S2sorted_barcode_onlyR2.fastq.gz S2sorted_barcode_onlyR1.fastq.gz > S2sorted_barcode_only.fastq.gz

gzip -c -d S2sorted_barcode_only.fastq.gz > S2sorted_barcode_only.fastq

perl MultiCodes.pl -out S2barcodes -index something -preseq GAGGTCTCT -postseq CGT -nPreExpected 0:120 < S2sorted_barcode_only.fastq
```
### S2 output 
```
Reads 561868 Multiplexed 561868 Usable(20) 499741 (88.9%) unique codes 52828 
Wrong presequence position: 58 reads (0.010%)
Wrote 52828 unique barcodes to S2barcodes.codes
Wrote number of barcodes seen a certain number of times to S2barcodes.counts; nOnce = 11687
Wrote 7732 off-by-1 pairs (8172 reads, fraction 0.016) to S2barcodes.close
Wrote S2barcodes.good
If 0.0% of reads are noise: diversity 71.8 K from total barcodes 52.8 K seen once 11.7 K seen twice 3.6 K
If 0.5% of reads are noise: diversity 62.0 K from total barcodes 50.3 K seen once 9.2 K seen twice 3.6 K
If 1.0% of reads are noise: diversity 54.0 K from total barcodes 47.8 K seen once 6.7 K seen twice 3.6 K
If 2.0% of reads are noise: diversity 43.2 K from total barcodes 42.8 K seen once 1.7 K seen twice 3.6 K
Aside from singletons and off-by-1s, see 40.8 K barcodes (97.5% of reads)
Barcodes with >= 45 reads each: 1.01% of codes (0.54 K), 7.98% of reads (39.9 K)
```


# repeat analysis with JJ sequences
```
bowtie2 -x plasmid2_index -q -p 4 --local -1 JJ_amplicon_S198_R1_001.fastq.gz -2 JJ_amplicon_S198_R2_001.fastq.gz -S JJbarcode_local.sam
```

# results log of JJ bowtie2 alignment 
```
2281225 reads; of these:
  2281225 (100.00%) were paired; of these:
    259471 (11.37%) aligned concordantly 0 times
    1921069 (84.21%) aligned concordantly exactly 1 time
    100685 (4.41%) aligned concordantly >1 times
    ----
    259471 pairs aligned concordantly 0 times; of these:
      3654 (1.41%) aligned discordantly 1 time
    ----
    255817 pairs aligned 0 times concordantly or discordantly; of these:
      511634 mates make up the pairs; of these:
        501091 (97.94%) aligned 0 times
        9515 (1.86%) aligned exactly 1 time
        1028 (0.20%) aligned >1 times
89.02% overall alignment rate
```
### convert .sam to .bam, sort, extract reads mapping to barcode region (262-281)

```
module load samtools 
# .sam to .bam
samtools view -S -b JJbarcode_local.sam > JJbarcode_local.bam

# sort bam file
samtools sort JJbarcode_local.bam -o JJbarcode_sorted.bam

# index the bam file
samtools index JJbarcode_sorted.bam

# export reads mapping to the barcode
samtools view JJbarcode_sorted.bam "barcode2_plasmid:262-281" > JJsorted_barcode_only.bam

```
# export fastq from sorted.bam using geneious 
# concatenate into single file to run barcode perl script on 
```
cat JJsorted_barcode_only_unpaired.fastq.gz JJsorted_barcode_onlyR2.fastq.gz JJsorted_barcode_onlyR1.fastq.gz > JJsorted_barcode_only.fastq.gz

gzip -c -d JJsorted_barcode_only.fastq.gz > JJsorted_barcode_only.fastq

perl MultiCodes.pl -out JJbarcodes -index something -preseq GAGGTCTCT -postseq CGT -nPreExpected 0:120 < JJsorted_barcode_only.fastq
```
### JJ output 
```
Reads 637992 Multiplexed 637992 Usable(20) 565387 (88.6%) unique codes 99513 
Wrong presequence position: 42 reads (0.007%)
Wrote 99513 unique barcodes to JJbarcodes.codes
Wrote number of barcodes seen a certain number of times to JJbarcodes.counts; nOnce = 19955
Wrote 5023 off-by-1 pairs (5423 reads, fraction 0.010) to JJbarcodes.close
Wrote JJbarcodes.good
If 0.0% of reads are noise: diversity 114.8 K from total barcodes 99.5 K seen once 20.0 K seen twice 13.0 K
If 0.5% of reads are noise: diversity 107.9 K from total barcodes 96.7 K seen once 17.1 K seen twice 13.0 K
If 1.0% of reads are noise: diversity 101.7 K from total barcodes 93.9 K seen once 14.3 K seen twice 13.0 K
If 2.0% of reads are noise: diversity 91.1 K from total barcodes 88.2 K seen once 8.6 K seen twice 13.0 K
Aside from singletons and off-by-1s, see 79.2 K barcodes (96.3% of reads)
Barcodes with >= 24 reads each: 1.11% of codes (1.11 K), 6.86% of reads (38.8 K)
```


