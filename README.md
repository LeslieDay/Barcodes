# Barcodes
This repository is for saving/organizing all code used for analyzing barcoded sequences 

### BarcodeCounts_pLD26.sh
Script for filtering and counting amplicon barcodes
- useful for determining # of random barcodes 
- based around MultiCodes.py from MagicPools https://journals.asm.org/doi/10.1128/mSystems.00143-17
  - manuscript code https://bitbucket.org/berkeleylab/feba/src/master/


Practicing using MapTnSeq.pl code from https://bitbucket.org/berkeleylab/feba/src/master/bin/MapTnSeq.pl
Description: MapTnSeq.pl identifies the random barcode, the junction between the transposon and the genome, and maps the remainder of the read to the genome. The result is a list of mapped reads with their barcodes and insertion locations.

```bash
module load blat
perl MapTnSeq.pl -genome GCF_000005845.2_ASM584v2_genomic.fna -model model_file -first Keio_ML9_index10_TAGCTT_L002_R1_001.fastq > Keio_test
 ```
 module_file is a single line file containing an example of what the tn/barcode insertion should look like
 Keio practice eg. nnnnnnCTAAGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACAGATGTGTATAAGAGACAG
 My samples actual- nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGTGTCAGACCGGGGACTTATCAGCCAACCTGTT
 Sequence following inverted repeat on plasmid- ATATCCATCACACTGGGCCGCTCGAGCATGC
 - The lower case n's indicate index that can be variable
 - Can add a second line for sequence that follows inverted repeat on plasmid to look for to filter out any sequences that had plasmid carry over past the IR 
 - outputs a file with barcode and location of insertion info 
  - Next use output file with DesignRandomPool.pl 
  - Description: uses the output of MapTnSeq.pl to identify barcodes that consistently map to a unique location in the genome. These are the useful barcodes. It also reports various metrics about the pool of mutants. (This step is done by invoking PoolStats.R.) Ideally, a mutant library has even insertions across the genome; has insertions in most of the protein-coding genes (except the essential ones); has a similar number of reads for insertions in most genes (i.e., no crazy positive selection for loss of a few genes); has insertions on both strands of genes (insertions on only one strand may indicate that the resistance marker's promoter is too weak); has tens of thousands or hundreds of thousands of useful barcodes; and the useful barcode s account for most of the reads.

In order for the DesignRandomPool.pl script to run properly had to move the module files FEBA_Utils.pm and Utils.pm into /home/perl5/lib/perl5

In addition to the MapTnSeq.pl output_file the DesignRandomPool.pl requries a gene_table which includes: scaffoldId, begin, end, strand, desc
and it should either include only protein-coding genes or it should
include the field 'type' with type=1 for protein-coding genes.
- use SetupOrg.pl to generate the gene_table
- For SetupOrg.pl to work had to import gffToGenes.pl, genesTabTranslation.pl, RegionGC.pl and update permissions to exicute (rw-- to rwx-)
- made a directory 'g' for saving genome information
/home/barcode_ratios/g
- creates a folder titled based on -name for outputting genome information
```bash
perl SetupOrg.pl -gff GCF_000005845.2_ASM584v2_genomic.gff -fna GCF_000005845.2_ASM584v2_genomic.fna -name Keio
```
For previous script had some issues with symbol errors think due to different versions of perl and bioperl communicating
- To see what versions of a module are loaded run
```bash
module info-loaded perl
module info-loaded bioperl
```
worked with perl/5.28.1 and bioperl/5.16.1
- if extra modules are loaded remove them using
```bash
module rm program/VERSIONINFO
```
Removed line for invoking PoolStats.R from DesignRandomPool.pl because communicating between perl MSI and me is dumb 

```bash
perl DesignRandomPool.pl -pool Keio_pool -genes g/Keio/genes.tab Keio_test
```

Pulled # of reads mapped from DesignRandomPool.pl output and ran R script independently
```bash
module load R
Rscript PoolStats.R Keio_pool g/Keio/genes.tab 3213323
```
Running these scripts on M Maripaludis sequencing data - S2 first
```bash
srun -N 1 --ntasks-per-node=4 --mem-per-cpu=12gb -t 4:00:00 -p interactive --pty bash
module load blat 
module load bioperl/5.16.1
perl MapTnSeq.pl -genome S2_GCF_000011585.1_ASM1158v1_genomic.fna -model model_file -first S2_S1_R1_001.fastq.gz > S2_mapping
```
### Output
>Parsed model model_file\
Barcodes of length 20, expected transposon region of 95\
Reads (gzipped) from S2_S1_R1_001.fastq.gz\
Read 5296936 reads\
Reads processed 5296936 Long-enough 5296936 Barcodes found 3320983\
Mapping attempted for 3207633 Mapped 2916841 Uniquely 2904674\
Hits past end of transposon: 135352 plus 20 weak/ambiguous; trumped hit to genome 16 times\
Proportions: Long-enough 1.000 Barcode 0.627 Attempted 0.606 Mapped 0.551 Past-end 0.026

[model_file.txt](https://github.com/LeslieDay/Barcodes/files/8318817/model_file.txt)

[MapTnSeq.pl script source](https://bitbucket.org/berkeleylab/feba/src/master/bin/MapTnSeq.pl)

[S2 Genome files](https://www.ncbi.nlm.nih.gov/datasets/genomes/?taxon=267377&utm_source=data-hub)
```
perl SetupOrg.pl -gff S2_GCF000011585.1genomic.gff -fna S2_GCF_000011585.1_ASM1158v1_genomic.fna -name S2_GeneTable
```
### Output
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

```bash
module load perl/modules.centos7.5.26.1

perl DesignRandomPool.pl -pool S2_RandomPool -genes g/S2_GeneTable/genes.tab S2_mapping

#if script won't run because cant find DBI then run following line
perl -MCPAN -e 'install DBI'
```
### DesignRandomPool.pl output
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

Analysis for JJ 
```bash
perl MapTnSeq.pl -genome JJ_fna_GCA_002945325.1_ASM294532v1_genomic.fna -model model_file -first JJ_S2_R1_001.fastq.gz > JJ_mapping
```
Must have only bioperl/5.16.1 and perl/5.28.1 installed for this step to work so may have to uninstall and reinstall perl/modules.centos7.5.26.1
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
Warning: stop codon within MMJJ_15440\

```bash
perl DesignRandomPool.pl -pool JJ_RandomPool -genes g/JJ_GeneTable/genes.tab JJ_mapping
```
### output
>Parsed model model_file\
Barcodes of length 20, expected transposon region of 95\
Reads (gzipped) from JJ_S2_R1_001.fastq.gz\
Read 5924086 reads\
Reads processed 5924086 Long-enough 5924086 Barcodes found 4910193\
Mapping attempted for 4674817 Mapped 4283950 Uniquely 4222641\
Hits past end of transposon: 8426 plus 0 weak/ambiguous; trumped hit to genome 0 times\
Proportions: Long-enough 1.000 Barcode 0.829 Attempted 0.789 Mapped 0.723 Past-end 0.001

```bash
module load perl/modules.centos7.5.26.1
perl DesignRandomPool.pl -pool JJ_RandomPool -genes g/JJ_GeneTable/genes.tab JJ_mapping
```
### output
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

Run R Script analysis on JJ and S2 DesignRandomPool output
```bash
module load R
Rscript PoolStats.R S2_RandomPool g/S2_GeneTable/genes.tab 3052193
```
### output 
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
Reads per million for hit proteins: median 147.76 mean 283.44

```bash
Rscript PoolStats.R JJ_RandomPool g/JJ_GeneTable/genes.tab 4292376
```

### output 
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
Reads per million for hit proteins: median 176.48 mean 284.25
```bash
perl MultiCodes.pl -out S2_barcodes -index S2_library -preseq GATGTCCACGAGGTCTCT -postseq CGT -nPreExpected 6 < S2_S1_R1_001.fastq
```
### output
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


```bash
# fastq file must be unzipped for MultiCodes to run properly
gunzip -c JJ_S2_R1_001.fastq.gz >JJ_S2_R1_001.fastq

perl MultiCodes.pl -out JJ_barcodes -index JJ_library -preseq GATGTCCACGAGGTCTCT -postseq CGT -nPreExpected 6 < JJ_S2_R1_001.fastq
```
### output
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
Barcodes with >= 130 reads each: 1.03% of codes (2.71 K), 9.60% of reads (469.9 K\
