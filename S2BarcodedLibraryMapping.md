
# ------------- S2 Tn Library mapping commands -------------
# -------------& corresponding output -------------
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

```bash
perl MultiCodes.pl -out S2_barcodes -index S2_library -preseq GATGTCCACGAGGTCTCT -postseq CGT -nPreExpected 6 < S2_S1_R1_001.fastq
```
### MultiCodes.pl S2 output
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
