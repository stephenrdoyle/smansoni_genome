# Normalising coverage to estimate copy number of gene family arrays: IPSE and OMEGA

Author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk

- IPSE and OMEGA are genes of high interest to the community
- in the new assembly, they are now represented as an array of multicopy genes, whereas perviously, they were single genes.
- however, there is still some degree of assembly collapse, evidence by a higher coverage in this region. This means that there are almost certainly more copies of these genes in reality, but are not present in the assembly
- want to try and normalise the coverage to estimate the true copy number of these genes.

- Matt B has carefully curated the gene models in artemis, and has provided a GFF of the annotatons
- will use this, together with read data from the single female mapped data, to estimate coverage.

### IPSE
```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/REPEATS

# ISPE coverage using samtools
cat mb_IPSE.gff | awk '{print $1,$4,$5,$9}' OFS="\t" > mb_IPSE.bed

samtools bedcov mb_IPSE.bed 6520_5_1_sorted.bam > mb_IPSE.coverage2

cat mb_IPSE.coverage2 | awk '{print $1,$2,$3,$4,$5/($3-$2)}' OFS="\t" | datamash median 5
#>53.7555

53.7555/44.6 = 1.21x higher coverage in region

16 * 1.21 = 20 likely copies of IPSE based on normalising coverage
```

### OMEGA
```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/REPEATS

# OMEGA coverage using samtools
cat mb_omega.gff | awk '{print $1,$4,$5,$9}' OFS="\t" > mb_omega.bed

samtools bedcov mb_omega.bed 6520_5_1_sorted.bam > mb_omega.coverage

cat mb_omega.coverage | awk '{print $1,$2,$3,$4,$5/($3-$2)}' OFS="\t" | datamash median 5

#> 78.35955/44.6 = 1.76x higher coverage in region

#> 5 copies * 1.76x = 8.8 likely copies of omega based on normalising coverage
#--- note - there are 7 copies of omega in this region, but only 5 are in a high coverage area. Two copies are not. See Supplementary Figure 1C
#--- therefore: 8.8 normalised copies + 2 addition copies =  10.8 copies, ~11 copies reported in the paper.

```
