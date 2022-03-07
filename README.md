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
