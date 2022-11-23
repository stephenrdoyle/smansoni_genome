# S. mansoni: V5 to V10 comparisons

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk




## Gene set comparison
- Want to understand how the genesets have changed between V5 and V10.
- the rough workflow is as follows:

1. gene IDs
     - shared gene ids
     - unique to V5 - deleted in V10
     - unique to V10 - new

2. shared gene ids
     - transcripts - mRNA length
          - same - unchanged
          - updated - length and/or % id changed
     - proteins  - length / %ID
          - same - unchanged
          - updated - length and/or % id changed
               - small change < 20%
               - large change >= 20%


- note: focus on first transcript of each gene



## data
```bash
# working directory
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/V5_vs_V10

# V5
schistosoma_mansoni.PRJEA36577.WBPS1.genomic.fa
schistosoma_mansoni.PRJEA36577.WBPS1.annotations.gff3

#--- v5 proteins
gffread schistosoma_mansoni.PRJEA36577.WBPS1.annotations.gff3 -g schistosoma_mansoni.PRJEA36577.WBPS1.genomic.fa -y sm_v5.proteins.fa



#--- v5 spliced exons
gffread schistosoma_mansoni.PRJEA36577.WBPS1.annotations.gff3 -g schistosoma_mansoni.PRJEA36577.WBPS1.genomic.fa -w sm_v5.spliced-exons.fa




# v10
/nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/SM_V10.genome.preWBP18checked.fa
/nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/SM_V10.annotation.preWBP18checked.gff3

#--- v10 proteins
gffread /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/SM_V10.annotation.preWBP18checked.gff3 -g /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/SM_V10.genome.preWBP18checked.fa -y sm_v10.proteins.fa



#--- v10 spliced exons
gffread /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/SM_V10.annotation.preWBP18checked.gff3 -g /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/SM_V10.genome.preWBP18checked.fa -w sm_v10.spliced-exons.fa



# fix sequence files for the DNA or aa sequence is on one line, rather than split across lines (helps parsing file later)
fastaq to_fasta -l0 sm_v5.proteins.fa sm_v5.proteins.l0.fa
fastaq to_fasta -l0 sm_v5.spliced-exons.fa sm_v5.spliced-exons.l0.fa
fastaq to_fasta -l0 sm_v10.proteins.fa sm_v10.proteins.l0.fa
fastaq to_fasta -l0 sm_v10.spliced-exons.fa sm_v10.spliced-exons.l0.fa

```


## Comparing gene IDs between V5 and V10

```bash
# extract V10 ids
grep ">" sm_v10.proteins.l0.fa | cut -f2 -d " " | sort | uniq > sm_v10.geneids.txt

# extract V5 ids
grep ">" sm_v5.proteins.l0.fa | cut -f2 -d " " | sort | uniq > sm_v5.geneids.txt


# compare files using diff
#- if gene present in both files, they'll be reported
#- if present in V5 by missing in V10 (deleted) , will show "<", eg
#     gene=Smp_014270						      <
#- if new in V10, will show a ">", eg.
# 							      >	gene=Smp_013815

diff -y <(sort sm_v5.geneids.txt) <(sort sm_v10.geneids.txt) | grep "<" > sm_V5.unique_genes.txt && wc -l sm_V5.unique_genes.txt
#> 4234 genes in V5 missing from V10

diff -y <(sort sm_v5.geneids.txt) <(sort sm_v10.geneids.txt) | grep ">" > sm_V10.unique_genes.txt && wc -l sm_V10.unique_genes.txt
#> 3299 genes new in V10, not present in V5

diff -y <(sort sm_v5.geneids.txt) <(sort sm_v10.geneids.txt) | grep -E "gene.*gene" >
#> 6597 genes with the same ID between V5 and V9



```

## Checking for updates in transcripts
```bash
# shared gene ids
diff -y <(sort sm_v5.geneids.txt) <(sort sm_v10.geneids.txt) | grep -E "gene.*gene" | cut -f1 | sed 's/gene=//g' > sm_shared_genes.txt


# extract first transcript from V5 of shared genes
>sm_v5.spliced-exons.shared.fa

while read NAME; do
     grep -A1 -m1 "${NAME}" sm_v5.spliced-exons.l0.fa >> sm_v5.spliced-exons.shared.fa
done < sm_shared_genes.txt

# extract first transcript from V10 of shared genes
>sm_v10.spliced-exons.shared.fa

while read NAME; do
     grep -A1 -m1 "${NAME}" sm_v10.spliced-exons.l0.fa >> sm_v10.spliced-exons.shared.fa
done < sm_shared_genes.txt




makeblastdb -in sm_v5.spliced-exons.shared.fa -out  sm_v5.spliced-exons.shared -dbtype nucl

blastn -outfmt 6 -db sm_v5.spliced-exons.shared -query sm_v10.spliced-exons.shared.fa



>sm_v5.proteins.shared.fa

while read NAME; do
     grep -A1 -m1 "${NAME}" sm_v5.proteins.l0.fa >> sm_v5.proteins.shared.fa
done < sm_shared_genes.txt



>sm_v10.proteins.shared.fa

while read NAME; do
     grep -A1 -m1 "${NAME}" sm_v10.proteins.l0.fa >> sm_v10.proteins.shared.fa
done < sm_shared_genes.txt

# fix proteins files for diamond - remove stop codons "."
sed -i '/^>/!s/\.//g' sm_v5.proteins.shared.fa
sed -i '/^>/!s/\.//g' sm_v10.proteins.shared.fa


module load diamond/2.0.12

diamond makedb --in sm_v5.proteins.shared.fa --db sm_v5.proteins.shared

diamond blastp --db sm_v5.proteins.shared.dmnd --query sm_v10.proteins.shared.fa --outfmt "6 score" | sed 's/transcript://g' | awk '$1==$2 {print}' | more


>sm_v5.proteins.n1.fa

while read NAME; do
     grep -A1 -m1 "${NAME}" sm_v5.proteins.l0.fa >> sm_v5.proteins.n1.fa
done < sm_v5.proteins.ids.txt

>sm_v10.proteins.n1.fa

while read NAME; do
     grep -A1 -m1 "${NAME}" sm_v10.proteins.l0.fa >> sm_v10.proteins.n1.fa
done < sm_v10.proteins.ids.txt

sed -i '/^>/!s/\.//g' sm_v5.proteins.n1.fa
sed -i '/^>/!s/\.//g' sm_v10.proteins.n1.fa

diamond makedb --in sm_v5.proteins.n1.fa --db sm_v5.proteins.n1

diamond blastp --db sm_v5.proteins.n1.dmnd --query sm_v10.proteins.n1.fa --outfmt "6" | sed 's/transcript://g' | awk '$1==$2 {print}' | more
