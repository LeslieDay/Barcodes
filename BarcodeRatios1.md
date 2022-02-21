Code used for analyzing samples sent for sequencing on February 7, 2022 
Has sequencing samples LD001-4

```bash
srun -N 1 --ntasks-per-node=4 --mem-per-cpu=12gb -t 4:00:00 -p interactive --pty bash
# unzip all sequencing files 
gunzip LD00*
# concatenate R1 and R2 sequencing files into single file
cat LD001* > LD001.fastq
cat LD002* > LD002.fastq
cat LD003* > LD003.fastq
cat LD004* > LD004.fastq

# python version 3.8.3 
python CountBarcodes.py LD001.fastq 30 index.txt barcodes.txt LD001
```
>
perl MultiCodes_8.pl -out LD003_barcodes -index something -preseq GAGGTCTCT -postseq CGT -nPreExpected 0:120 < LD003.fastq
