>load fastq files and MultiCodes.pl script into ~/EH_barcodeRatio/tn_barcodes using filezilla
> 
>unzip files 

```
gzip -c -d EHH023_S236_R1_001.fastq.gz > EHH023_R1.fastq
gzip -c -d EHH021_S234_R1_001.fastq.gz > EHH021_R1.fastq
```
>extract barcodes 

```
perl MultiCodes.pl -out EHH023_barcodes -index EHH023 -preseq GAGGTCTCT -postseq CGT -nPreExpected 0:120 < EHH023_R1.fastq
```
Reads 5708042 Multiplexed 5708042 Usable(20) 263775 (4.6%) unique codes 63665 \
Wrong presequence position: 14469 reads (0.253%)\
Wrote 63665 unique barcodes to EHH023_barcodes.codes\
Wrote number of barcodes seen a certain number of times to EHH023_barcodes.counts; nOnce = 16509\
Wrote 2894 off-by-1 pairs (3262 reads, fraction 0.012) to EHH023_barcodes.close\
Wrote EHH023_barcodes.good\
If 0.0% of reads are noise: diversity 76.1 K from total barcodes 63.7 K seen once 16.5 K seen twice 11.0 K\
If 0.5% of reads are noise: diversity 72.9 K from total barcodes 62.3 K seen once 15.2 K seen twice 11.0 K\
If 1.0% of reads are noise: diversity 69.8 K from total barcodes 61.0 K seen once 13.9 K seen twice 11.0 K\
If 2.0% of reads are noise: diversity 64.1 K from total barcodes 58.4 K seen once 11.2 K seen twice 11.0 K\
Aside from singletons and off-by-1s, see 46.9 K barcodes (93.5% of reads)\
Barcodes with >= 18 reads each: 1.05% of codes (0.67 K), 5.68% of reads (15.0 K)

```
perl MultiCodes.pl -out EHH021_barcodes -index EHH021 -preseq GAGGTCTCT -postseq CGT -nPreExpected 0:120 < EHH021_R1.fastq
```

Reads 6280103 Multiplexed 6280103 Usable(20) 272724 (4.3%) unique codes 77811 \
Wrong presequence position: 15687 reads (0.250%) \
Wrote 77811 unique barcodes to EHH021_barcodes.codes\
Wrote number of barcodes seen a certain number of times to EHH021_barcodes.counts; nOnce = 23099\
Wrote 2458 off-by-1 pairs (2670 reads, fraction 0.010) to EHH021_barcodes.close\
Wrote EHH021_barcodes.good\
If 0.0% of reads are noise: diversity 95.3 K from total barcodes 77.8 K seen once 23.1 K seen twice 15.2 K\
If 0.5% of reads are noise: diversity 91.9 K from total barcodes 76.4 K seen once 21.7 K seen twice 15.2 K\
If 1.0% of reads are noise: diversity 88.7 K from total barcodes 75.1 K seen once 20.4 K seen twice 15.2 K\
If 2.0% of reads are noise: diversity 82.6 K from total barcodes 72.4 K seen once 17.6 K seen twice 15.2 K\
Aside from singletons and off-by-1s, see 54.6 K barcodes (91.4% of reads)\
Barcodes with >= 15 reads each: 1.02% of codes (0.79 K), 5.52% of reads (15.1 K)


This is the code to combine the amplicon barcode files with the random pool library (barcodes and their location from original library mapping). I cannot get this counting to work correctly. I think it is because of the order of the barcodes in the different files. The script works from top to bottom in each so if they are not ordered properly the code is assumed to not be present. Need to think more on this. 
```
perl combineBarSeq.pl JJ_combineOut JJ_RandomPool EHH023_barcodes.codes EHH021_barcodes.codes
```
