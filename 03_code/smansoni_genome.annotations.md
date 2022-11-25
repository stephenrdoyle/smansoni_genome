# Schistosoma mansoni V10 genome: busco analysis

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk



## Comparison of V5 to V9 annotaitons
```bash
# working dir
/nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V5_v_V9_ANNOTATION

# get genome and annotation files
cp ../../papers/genome_improvement/smansoni/protasio_2012/schistosoma_mansoni.PRJEA36577.WBPS1.annotations.gff3 .
cp ../../papers/genome_improvement/smansoni/protasio_2012/schistosoma_mansoni.PRJEA36577.WBPS1.genomic.fa .
cp ../../papers/genome_improvement/smansoni/buddenborg_2021/SM_V9_ENA.fa .
cp ../SM_V9_16Mar.gff .

cp ../../papers/genome_improvement/smansoni/WBP_V7/schistosoma_mansoni.PRJEA36577.WBPS15.genomic.fa .
cp ../../papers/genome_improvement/smansoni/WBP_V7/schistosoma_mansoni.PRJEA36577.WBPS15.annotations.gff3 .

# extract CDS sequences
gffread -g schistosoma_mansoni.PRJEA36577.WBPS1.genomic.fa -x SM_V5_CDS.fa schistosoma_mansoni.PRJEA36577.WBPS1.annotations.gff3

gffread -g SM_V9_ENA.fa -x SM_V9_CDS.fa SM_V9_16Mar.gff

gffread -g schistosoma_mansoni.PRJEA36577.WBPS15.genomic.fa -x SM_V7_CDS.fa schistosoma_mansoni.PRJEA36577.WBPS15.annotations.gff3


# make a blast database of the V5 sequences
makeblastdb -in SM_V5_CDS.fa -parse_seqids -dbtype nucl

blastn -db SM_V5_CDS.fa -query SM_V9_CDS.fa -outfmt 6 -out V5vV9.result.txt

makeblastdb -in SM_V7_CDS.fa -parse_seqids -dbtype nucl

blastn -db SM_V7_CDS.fa -query SM_V9_CDS.fa -outfmt 6 -out V7vV9.result.txt

```
