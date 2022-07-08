# Updates to Schistosoma mansoni genome structure: v9 to v10


## remove the indel in the W copy U2AF
- the W copy of U2AF has a frameshift caused by an additional base present that shouldnt be there
- evidence that this is a technical error in the reference includes:
     - RNAseq read data that maps well to the gene, but all reads have an indel
     - genomic reads from the reference strain all contain the indel
     - SNP data from the global cohort show a fixed variant, suggesting it is not variable
- in the v9 paper, we descibed the presence of the indel, however, it caused some issue with the reviewer and will likely cause problems in the future, give its likely importance in sex determination .
- give some other large changes to the annotation, we thoughtv it best to fix it.


```bash
cd ~/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10

samtools faidx SM_V9_21Feb.fa SM_V9_WSR > SM_V9_WSR_v9.fa


# remove the indel by substituting the correct sequence
sed 's/CCCTACCTACGATCAC/CCCTACTACGATCAC/' SM_V9_WSR_v9.fa > SM_V9_WSR_v10.fa


# collate the other scaffolds, without the old WSR
cut -f1 SM_V9_21Feb.fa.fai | grep -v "SM_V9_WSR" | while read -r NAME; do
     samtools faidx SM_V9_21Feb.fa ${NAME} >> SM_v10.fa;
     done

# bring it together
cat SM_v10.fa SM_V9_WSR_v10.fa > tmp; mv tmp SM_v10.fa

```
