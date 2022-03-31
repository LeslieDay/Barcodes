
# Order to run scripts for analysis
### 1-4 Required for initial transposon library mapping
### 5-6 Required for BarSeq analysis 
1. MapTnSeq.pl
2. SetupOrg.pl
3. DesignRandomPool.pl 
4. PoolStats.R
6. MultiCodes.pl
7. combineBarSeq.pl

# Required Modules
module load blat
module load R
module load bioperl

# Required scripts 
## All scripts must be located in working directory & be sure file permissions are executable (--rxw)
## [RB-TnSeq Scripts](https://bitbucket.org/berkeleylab/feba/src/master/)
- [MapTnSeq.pl](https://bitbucket.org/berkeleylab/feba/src/master/bin/MapTnSeq.pl)
- [SetupOrg.pl](https://bitbucket.org/berkeleylab/feba/src/master/bin/SetupOrg.pl)
- [Edited -DesignRandomPool.pl](https://github.com/LeslieDay/Barcodes/blob/main/DesignRandomPool.pl) 
-- Removed line to automatically execute PoolStats.R 
-- Original version of script [DesignRandomPool.pl](https://bitbucket.org/berkeleylab/feba/src/master/bin/DesignRandomPool.pl)
- [PoolStats.R](https://bitbucket.org/berkeleylab/feba/src/master/lib/PoolStats.R)
-- To manually run if using edited DesignRandomPool.pl script provide DesignRandomPool.pl stat of $nMapped from output "Read $nMapped mapped reads for..."
- [MultiCodes.pl](https://bitbucket.org/berkeleylab/feba/src/master/bin/MultiCodes.pl)
- [combineBarSeq.pl](https://bitbucket.org/berkeleylab/feba/src/master/bin/combineBarSeq.pl)

## Scripts dependencies of above scripts i.e. must be present in directory and executable for main scripts to run
- [gbkToFaa.pl](https://bitbucket.org/berkeleylab/feba/src/master/bin/gbkToFaa.pl)
- [gbkToSeq.pl](https://bitbucket.org/berkeleylab/feba/src/master/bin/gbkToSeq.pl)
- [gbkToSeq2.pl](https://bitbucket.org/berkeleylab/feba/src/master/bin/gbkToSeq2.pl)

### From gfftools repository 
[genbank2gff.pl](https://github.com/ihh/gfftools/blob/master/genbank2gff.pl)
 

## Upload model file specific to our plasmid 
### line 1 = nnnnnnn (expect 6bp index) 20N (20 nucleotide barcode) followed by sequence to end of inverted repeat
### line 2 = plasmid backbone following inverted repeat, used to remove plasmid contamination sequences

- [model_file.txt](https://github.com/LeslieDay/Barcodes/files/8318817/model_file.txt)
nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGTGTCAGACCGGGGACTTATCAGCCAACCTGT
ATATCCATCACACTGGGCCGCTCGAGCATGC

[S2 Genome files](https://www.ncbi.nlm.nih.gov/datasets/genomes/?taxon=267377&utm_source=data-hub)

# ------------- S2 Tn Library mapping -------------
## 1. MapTnSeq.pl
```bash
# login to MSI and enter directory with sequencing data and scripts 
# my working directory = /home/barcodedTnLibrary
cd barcodedTnLibrary
# request interactive node time on MSI 
srun -N 1 --ntasks-per-node=4 --mem-per-cpu=12gb -t 4:00:00 -p interactive --pty bash
module load perl/5.28.1
module load blat 
module load bioperl/5.16.1
# Load required packages
# Run MapTnSeq.pl 
# Generates tab-delimited file with fields read,barcode,scaffold,pos,strand,uniq,qBeg,qEnd,score,identity where qBeg and qEnd are the positions in the read, after the transposon sequence is removed, that match the genome.

perl MapTnSeq.pl -genome S2_GCF_000011585.1_ASM1158v1_genomic.fna -model model_file -first S2_S1_R1_001.fastq.gz > S2_mapping

# Run SetupOrg.pl
# create directory /g for saving genome info 
mkdir g
perl SetupOrg.pl -gff S2_GCF000011585.1genomic.gff -fna S2_GCF_000011585.1_ASM1158v1_genomic.fna -name S2_GeneTable

module load perl/modules.centos7.5.26.1

perl DesignRandomPool.pl -pool S2_RandomPool -genes g/S2_GeneTable/genes.tab S2_mapping

#if script won't run because cant find DBI then run following line
perl -MCPAN -e 'install DBI'

# Run R Script analysis on S2 DesignRandomPool output
module load R
Rscript PoolStats.R S2_RandomPool g/S2_GeneTable/genes.tab 3041567

# Unzip sequencing files for MultiCodes.pl analysis
gunzip -c S2_S1_R1_001.fastq.gz >S2_S1_R1_001.fastq
#Analyze reads with MultiCodes.pl
perl MultiCodes.pl -out S2_barcodes -index S2_library -preseq GATGTCCACGAGGTCTCT -postseq CGT -nPreExpected 6 < S2_S1_R1_001.fastq

```
# ------------- S2 Tn Library mapping commands -------------
#    ------------- & corresponding output -------------
```
perl MapTnSeq.pl -genome S2_GCF_000011585.1_ASM1158v1_genomic.fna -model model_file -first S2_S1_R1_001.fastq.gz > S2_mapping
```
### MapTnSeq.pl S2 Output
>Parsed model model_file\
Barcodes of length 20, expected transposon region of 95\
Reads (gzipped) from S2_S1_R1_001.fastq.gz\
Read 5296936 reads\
Reads processed 5296936 Long-enough 5296936 Barcodes found 3320983\
Mapping attempted for 3207633 Mapped 2916841 Uniquely 2904674\
Hits past end of transposon: 135352 plus 20 weak/ambiguous; trumped hit to genome 16 times\
Proportions: Long-enough 1.000 Barcode 0.627 Attempted 0.606 Mapped 0.551 Past-end 0.026

```
perl SetupOrg.pl -gff S2_GCF000011585.1genomic.gff -fna S2_GCF_000011585.1_ASM1158v1_genomic.fna -name S2_GeneTable
```
### SetupOrg.pl S2 Output
>Replacement list is longer than search list at /soft/bioperl/5.16.1/lib/site_perl/5.16.1/Bio/Range.pm line 251.\
Replacement list is longer than search list at /soft/bioperl/5.16.1/lib/site_perl/5.16.1/Bio/Perl.pm line 627.\
Warning: stop codon within MMP_RS00795\
Warning: stop codon within MMP_RS02705\
Warning: stop codon within MMP_RS04695\
Warning: mstop codon within MMP_RS06685\
Warning: stop codon within MMP_RS07120\
Warning: stop codon within MMP_RS08715\
Warning: stop codon within MMP_RS09070\
Warning: stop codon within MMP_RS08735\
Warning: stop codon within MMP_RS08740

```
perl DesignRandomPool.pl -pool S2_RandomPool -genes g/S2_GeneTable/genes.tab S2_mapping
```
### DesignRandomPool.pl S2 output
>Reading mapping files:\
S2_mapping\
S2_mapping has 3,052,193 mappings, 99,909 barcodes, 49,736 non-unique barcodes\
Read 3041567 mapped reads for 94028 distinct barcodes\
(Skipped 10626 reads with qBeg > 3)\
40130 barcodes seen 10 or more times, map 38179 (minFrac 0.75 minRatio 8)\
FewGood	150	0.0037\
LoRatio	488	0.0122\
NoCons	1195	0.0298\
PastEnd	118	0.0029\
Usable	38179	0.9514\
Masked 6 off-by-1 barcodes (69 reads) leaving 38173 barcodes\
Reads for those barcodes: 2602847 of 3041567 (85.6%)\
Chao2 estimate of #barcodes present (may be inflated for sequencing error): 594490

```bash
Rscript PoolStats.R S2_RandomPool g/S2_GeneTable/genes.tab 3041567
```
### S2 PoolStats.R output 
>38173 insertions in genome are at 28673 different locations\
Found 21259 insertions (16573 distinct locations) in central 10-90% of genes\
Found central insertions for 1615 of 1741 protein-coding genes\
Hit rate in (crude) likely essentials: 0.42 other 0.92\
Wrote proteins of 300nt or more with no good insertions to S2_RandomPool.unhit\
Wrote 29 genes with surprising insertions in central 10-90% to S2_RandomPool.surprise\
Wrote read and strain counts for hit genes to S2_RandomPool.hit\
Strains per hit protein: median 8 mean 13.1\
Gene and transposon on same strand: 50.9%\
Reads per hit protein: median 451 mean 865.1 bias (ratio) 1.92\
Reads per million for hit proteins: median 148.28 mean 284.43

```bash
perl MultiCodes.pl -out S2_barcodes -index S2_library -preseq GATGTCCACGAGGTCTCT -postseq CGT -nPreExpected 6 < S2_S1_R1_001.fastq
```
### MultiCodes.pl S2 Output
>Reads 5296936 Multiplexed 5296936 Usable(20) 3316472 (62.6%) unique codes 131171 \
Wrong presequence position: 139 reads (0.003%)\
Wrote 131171 unique barcodes to S2_barcodes.codes\
Wrote number of barcodes seen a certain number of times to S2_barcodes.counts; nOnce = 77482\
Wrote 60938 off-by-1 pairs (90679 reads, fraction 0.027) to S2_barcodes.close\
Wrote S2_barcodes.good\
If 0.0% of reads are noise: diversity 1034.5 K from total barcodes 131.2 K seen once 77.5 K seen twice 3.3 K\
If 0.5% of reads are noise: diversity 672.6 K from total barcodes 114.6 K seen once 60.9 K seen twice 3.3 K\
If 1.0% of reads are noise: diversity 393.5 K from total barcodes 98.0 K seen once 44.3 K seen twice 3.3 K\
If 2.0% of reads are noise: diversity 83.6 K from total barcodes 64.8 K seen once 11.2 K seen twice 3.3 K\
Aside from singletons and off-by-1s, see 50.0 K barcodes (97.2% of reads)\
Barcodes with >= 209 reads each: 1.01% of codes (1.33 K), 12.23% of reads (405.4 K)

# ------------- Extra S2 analysis with outputs -------------
# ------------- Attempting to optimize mapping -------------
```
# DesignRandomPool.pl with minimum mapping number = 2 
# default is a minimum of ten sequences mapped in order to count barcode 
perl DesignRandomPool.pl -pool S2_RandomPool_minN2 -genes g/S2_GeneTable/genes.tab S2_mapping -minN 2
```
## DesignRandomPool S2 output
>Reading mapping files:\
S2_mapping\
S2_mapping has 3,052,193 mappings, 99,909 barcodes, 49,736 non-unique barcodes\
Read 3041567 mapped reads for 94028 distinct barcodes\
(Skipped 10626 reads with qBeg > 3)\
47332 barcodes seen 2 or more times, map 45207 (minFrac 0.75 minRatio 8)\
FewGood	175	0.0037\
LoRatio	522	0.0110\
NoCons	1298	0.0274\
PastEnd	130	0.0027\
Usable	45207	0.9551\
Masked 1805 off-by-1 barcodes (4619 reads) leaving 43409 barcodes\
Reads for those barcodes: 2635422 of 3041567 (86.6%)\
Chao2 estimate of #barcodes present (may be inflated for sequencing error): 594490


# ------------- JJ Tn Library mapping -------------
```bash
perl MapTnSeq.pl -genome JJ_fna_GCA_002945325.1_ASM294532v1_genomic.fna -model model_file -first JJ_S2_R1_001.fastq.gz > JJ_mapping

# Must have only bioperl/5.16.1 and perl/5.28.1 installed for this step to work so may have to uninstall and reinstall perl/modules.centos7.5.26.1
perl SetupOrg.pl -gff JJ_gff_GCA_002945325.1_ASM294532v1_genomic.gff -fna JJ_fna_GCA_002945325.1_ASM294532v1_genomic.fna -name JJ_GeneTable

module load perl/modules.centos7.5.26.1
perl DesignRandomPool.pl -pool JJ_RandomPool -genes g/JJ_GeneTable/genes.tab JJ_mapping

Rscript PoolStats.R JJ_RandomPool g/JJ_GeneTable/genes.tab 4277927

# fastq file must be unzipped for MultiCodes to run properly
gunzip -c JJ_S2_R1_001.fastq.gz >JJ_S2_R1_001.fastq

perl MultiCodes.pl -out JJ_barcodes -index JJ_library -preseq GATGTCCACGAGGTCTCT -postseq CGT -nPreExpected 6 < JJ_S2_R1_001.fastq
```
# ------------- JJ Tn Library mapping commands -------------
#    ------------- & corresponding output -------------

```bash
perl MapTnSeq.pl -genome JJ_fna_GCA_002945325.1_ASM294532v1_genomic.fna -model model_file -first JJ_S2_R1_001.fastq.gz > JJ_mapping
```
Parsed model model_file
Barcodes of length 20, expected transposon region of 95
Reads (gzipped) from JJ_S2_R1_001.fastq.gz
Read 5924086 reads
Reads processed 5924086 Long-enough 5924086 Barcodes found 4910193
Mapping attempted for 4674817 Mapped 4283950 Uniquely 4222641
Hits past end of transposon: 8426 plus 0 weak/ambiguous; trumped hit to genome 0 times
Proportions: Long-enough 1.000 Barcode 0.829 Attempted 0.789 Mapped 0.723 Past-end 0.001

```bash
perl SetupOrg.pl -gff JJ_gff_GCA_002945325.1_ASM294532v1_genomic.gff -fna JJ_fna_GCA_002945325.1_ASM294532v1_genomic.fna -name JJ_GeneTable
```
### SetupOrg.pl for JJ output

>Replacement list is longer than search list at /soft/bioperl/5.16.1/lib/site_perl/5.16.1/Bio/Range.pm line 251.\
Replacement list is longer than search list at /soft/bioperl/5.16.1/lib/site_perl/5.16.1/Bio/Perl.pm line 627.\
Warning: stop codon within MMJJ_00870\
Warning: stop codon within MMJJ_08070\
Warning: stop codon within MMJJ_09470\
Warning: stop codon within MMJJ_11320\
Warning: stop codon within MMJJ_11330\
Warning: stop codon within MMJJ_11360\
Warning: stop codon within MMJJ_11380\
Warning: stop codon within MMJJ_14570\
Warning: stop codon within MMJJ_15440

```bash
perl DesignRandomPool.pl -pool JJ_RandomPool -genes g/JJ_GeneTable/genes.tab JJ_mapping
```
### JJ DesignRandomPool.pl output
>Reading mapping files:\
JJ_mapping\
JJ_mapping has 4,292,376 mappings, 202,668 barcodes, 121,633 non-unique barcodes\
Read 4277927 mapped reads for 197929 distinct barcodes\
(Skipped 14449 reads with qBeg > 3)\
94143 barcodes seen 10 or more times, map 81382 (minFrac 0.75 minRatio 8)\
FewGood	1146	0.0122\
LoRatio	3324	0.0353\
NoCons	8283	0.0880\
PastEnd	8	0.0001\
Usable	81382	0.8645\
Masked 18 off-by-1 barcodes (221 reads) leaving 81364 barcodes\
Reads for those barcodes: 3218100 of 4277927 (75.2%)\
Chao2 estimate of #barcodes present (may be inflated for sequencing error): 779112


```bash
Rscript PoolStats.R JJ_RandomPool g/JJ_GeneTable/genes.tab 4277927
```

### JJ PoolStats.R output 
>81364 insertions in genome are at 50130 different locations\
Found 47387 insertions (30575 distinct locations) in central 10-90% of genes\
Found central insertions for 1522 of 1815 protein-coding genes\
Hit rate in (crude) likely essentials: 0.31 other 0.84\
Wrote proteins of 300nt or more with no good insertions to JJ_RandomPool.unhit\
Wrote 22 genes with surprising insertions in central 10-90% to JJ_RandomPool.surprise\
Wrote read and strain counts for hit genes to JJ_RandomPool.hit\
Strains per hit protein: median 21 mean 31.1\
Gene and transposon on same strand: 50.1%\
Reads per hit protein: median 758 mean 1220.1 bias (ratio) 1.61\
Reads per million for hit proteins: median 177.07 mean 285.21


```bash
perl MultiCodes.pl -out JJ_barcodes -index JJ_library -preseq GATGTCCACGAGGTCTCT -postseq CGT -nPreExpected 6 < JJ_S2_R1_001.fastq
```
### MultiCodes.pl JJ Output
Reads 5924086 Multiplexed 5924086 Usable(20) 4895594 (82.6%) unique codes 262975\
Wrong presequence position: 166 reads (0.003%)\
Wrote 262975 unique barcodes to JJ_barcodes.codes\
Wrote number of barcodes seen a certain number of times to JJ_barcodes.counts; nOnce = 128237\
Wrote 94791 off-by-1 pairs (115988 reads, fraction 0.024) to JJ_barcodes.close\
Wrote JJ_barcodes.good\
If 0.0% of reads are noise: diversity 1465.6 K from total barcodes 263.0 K seen once 128.2 K seen twice 6.8 K\
If 0.5% of reads are noise: diversity 1025.8 K from total barcodes 238.5 K seen once 103.8 K seen twice 6.8 K\
If 1.0% of reads are noise: diversity 673.7 K from total barcodes 214.0 K seen once 79.3 K seen twice 6.8 K\
If 2.0% of reads are noise: diversity 232.3 K from total barcodes 165.1 K seen once 30.3 K seen twice 6.8 K\
Aside from singletons and off-by-1s, see 128.3 K barcodes (96.9% of reads)\
Barcodes with >= 130 reads each: 1.03% of codes (2.71 K), 9.60% of reads (469.9 K)
