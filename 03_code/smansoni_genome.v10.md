# Schistosoma mansoni V10 genome: Updates to Schistosoma mansoni genome structure: v9 to v10

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk

- workflow describing changes per chromosome in the genome and annotation from V9 to V10

```bash
# working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10
```


## Chromosome 1
### merge original and updated annotations from apollo
```bash
GENOME=SM_V9_21Feb.fa
CHROMOSOME=SM_V9_1

mkdir /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}


samtools faidx ../${GENOME} ${CHROMOSOME} > ${CHROMOSOME}.fa


# apollo annotations saved as "chr1_v9.apollo.gff3.gz"
# ie. scp Annotations.gff3.gz sd21@farm5-head1:~/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/SM_V9_1/SM_V9_1_v9.apollo.gff3.gz
gunzip -f ${CHROMOSOME}_v9.apollo.gff3.gz

# get orignal annotations for chromosome 1
grep "${CHROMOSOME}" ../../SM_V9_16Mar.gff > ${CHROMOSOME}_v9.original.gff3



# step 1: fix gene ID and relationship in mRNA ID
awk -F'[\t]' '$3=="gene" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=..........|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_gene_IDs_NAMEs.txt


# step 2: fix mRNA IDs and all descendant IDs
awk -F'[\t]' '$3=="mRNA" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=............|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_mRNA_IDs_NAMEs.txt


# get product IDs
awk -F'[\t]'  '$3=="mRNA" {print $9}' ${CHROMOSOME}_v9.original.gff3 | sed 's/;/\t/g' | awk -F '[\t]' '{print $1,$3}' OFS="\t" | sed -e 's/ID=/Name=/g' | awk -F '[\t]' '{print $1,$1";"$2}' OFS="\t" > original_IDs_product.txt


cat apollo_gene_IDs_NAMEs.txt apollo_mRNA_IDs_NAMEs.txt > apollo_IDs_NAMEs.txt

# replace apollo unique ID with Smp ID in apollo GFF
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.renamed.gff3 apollo_IDs_NAMEs.txt ${CHROMOSOME}_v9.apollo.gff3


# split annotation into mRNA and other, so that product id's can be added back to mRNAs specifically. Initial testing showed that some CDSs were getting the product IDs added incorrectly, so to account for this, splitting them
awk '$3=="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp
awk '$3!="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp


# add product descriptions to mRNAs only
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 original_IDs_product.txt ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp


# bring the updated mRNAs and other annotations back together.
cat ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp | sort -k4,4n > ${CHROMOSOME}_v9.apollo.renamed2.gff3



# remove overlapping models which have been updated
bedtools intersect -s -v  -f 0.1 -a ${CHROMOSOME}_v9.original.gff3 -b ${CHROMOSOME}_v9.apollo.renamed2.gff3 > ${CHROMOSOME}_v9.original.gff3.no-overlaps


# remove updated gene IDs from original annotation. Might not be strictly necessary, but there were a couple of instances where there were old and new annotations kept.
while read ID GENE; do
     sed -i "/Liftoff.*${GENE}/d" ${CHROMOSOME}_v9.original.gff3.no-overlaps;
done < apollo_gene_IDs_NAMEs.txt



# bring apollo and unchanged original annotations back together.
cat ${CHROMOSOME}_v9.apollo.renamed2.gff3 ${CHROMOSOME}_v9.original.gff3.no-overlaps | grep -v "#" | sort -k4,4n > ${CHROMOSOME}_v9.updated.gff3



# manual check to see if there are duplicated gene IDs
cat ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="gene") print}' | sed -e 's/owner=irisadmin@local.host;//g' -e 's/owner=mb4@sanger.ac.uk//g' -e 's/owner=irisadmin@local.host,mb4@sanger.ac.uk//g' -e 's/,irisadmin@local.host//g' | sed 's/;/\t/g' | cut -f9 | sort | uniq -c | sort



# chromosome 1 specific changes - checked in apollo
2 ID=Smp_009760 - remove liftoff
2 ID=Smp_012750 - fixed
2 ID=Smp_032260 - remove liftoff
2 ID=Smp_037780 - remove liftoff
2 ID=Smp_200090 - remove liftoff
2 ID=Smp_200190 - remove liftoff
2 ID=Smp_205000 - fixed
2 ID=Smp_320500 - remove liftoff
2 ID=Smp_331740 - to delete
2 ID=Smp_336630 - fixed
2 ID=Smp_346560 - to delete

# removing old Liftoff versions
sed -i "/Liftoff.*Smp_009760/d" SM_V9_1_v9.updated.gff3
sed -i "/Liftoff.*Smp_032260/d" SM_V9_1_v9.updated.gff3
sed -i "/Liftoff.*Smp_037780/d" SM_V9_1_v9.updated.gff3
sed -i "/Liftoff.*Smp_200090/d" SM_V9_1_v9.updated.gff3
sed -i "/Liftoff.*Smp_200190/d" SM_V9_1_v9.updated.gff3
sed -i "/Liftoff.*Smp_320500/d" SM_V9_1_v9.updated.gff3

# deleting rubbish genes
sed -i "/Smp_331740/d" SM_V9_1_v9.updated.gff3
sed -i "/Smp_346560/d" SM_V9_1_v9.updated.gff3
sed -i "/Smp_003760/d" SM_V9_1_v9.updated.gff3
sed -i "/Smp_318520/d" SM_V9_1_v9.updated.gff3
sed -i "/Smp_301580/d" SM_V9_1_v9.updated.gff3
sed -i "/Smp_319280/d" SM_V9_1_v9.updated.gff3
sed -i "/Smp_334170/d" SM_V9_1_v9.updated.gff3
sed -i "/Smp_204580/d" SM_V9_1_v9.updated.gff3
sed -i "/Smp_200150/d" SM_V9_1_v9.updated.gff3




# checking for mini-introns placed in the gff in apollo to fix a gene model, but reflect an error in the genome that needs fixing
gt gff3 -tidy -addintrons -retainids SM_V9_1_v9.updated.gff3 | awk '{if($3=="intron") print $9,$5-$4}' | sort -k2,2nr | grep -v "warning"

# mini introns - gene ID, length
Parent=Smp_003340.1 1
Parent=Smp_007380.1 1
Parent=Smp_059790.1 1
Parent=Smp_059790.1 1
Parent=Smp_144140.1 1
Parent=Smp_166290.1 1
Parent=Smp_202170.1 1
Parent=Smp_243170.1 1
Parent=Smp_301670.1 1
Parent=Smp_317370.1 1
Parent=Smp_332790.1 1
Parent=Smp_334660.1 1
Parent=Smp_335410.1 1
Parent=Smp_337200.1 1
Parent=Smp_340320.1 1
Parent=Smp_342750.1 1
Parent=Smp_345080.1 1
Parent=Smp_345280.1 1
Parent=Smp_019380.1 0
Parent=Smp_059790.1 0
Parent=Smp_141160.1 0
Parent=Smp_337430.1 0

```

### repair indels
```bash
FASTA=SM_V9_1.fa
GFF=SM_V9_1_v9.updated.gff3


fastaq to_fasta -l0 ${FASTA} FASTA.repaired.tmp
cp ${GFF} GFF.repaired.tmp


# read in the list of genes. Note, list needs to be reverse sorted by position to make sure the corrections dont interfere with each other.
sed 1d repair_list.txt | sort -k3rn | while read -r GENE CHROMOSOME POSITION SIZE OLD_SEQUENCE NEW_SEQUENCE; do

# check if repair site is unique
if [ "$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)" -eq 1 ]; then
     # introduce the new sequence to the reference
     sed -i "s/${OLD_SEQUENCE}/${NEW_SEQUENCE}/" FASTA.repaired.tmp

     # update the coordinates in the GFF
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($4>POSITION) print $1,$2,$3,$4+SIZE,$5,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" GFF.repaired.tmp |
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($5>POSITION) print $1,$2,$3,$4,$5+SIZE,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" > GFF.tmp; mv GFF.tmp GFF.repaired.tmp

     echo "${GENE} target site ${POSITION} has been fixed, with a indel of ${SIZE} bp."

else
     COUNT=$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)
     echo "${GENE} target site ${POSITION} is not unique - it was found ${COUNT} times. Check and improve the target sequence so it is unique and try again."
fi;

done

mv FASTA.repaired.tmp ${FASTA%.fa.*$}.repaired.fa
mv GFF.repaired.tmp ${GFF%.gff.*$}.repaired.gff3
rm *.tmp*

```

### update from V9 to V10
```bash
# update chromosome name
sed 's/SM_V9_1/SM_V10_1/g' SM_V9_1.fa.repaired.fa > SM_V10_1.fa
samtools faidx SM_V10_1.fa

# fix chromosome name, and use gffread to close indels
sed 's/SM_V9_1/SM_V10_1/g' SM_V9_1_v9.updated.gff3.repaired.gff3 | gffread - -Z -F -O -o SM_V10_1.gff3

gffread SM_V10_1.gff3 -g SM_V10_1.fa -y SM_V10_1.proteins.fa

```




## Chromosome 2
### merge original and updated annotations from apollo
```bash
GENOME=SM_V9_21Feb.fa
CHROMOSOME=SM_V9_2

mkdir /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}


samtools faidx ../${GENOME} ${CHROMOSOME} > ${CHROMOSOME}.fa


# get orignal annotations for chromosome
grep "${CHROMOSOME}" ../../SM_V9_16Mar.gff > ${CHROMOSOME}_v9.original.gff3


# apollo annotations saved as "SM_V9_2_v9.apollo.gff3.gz"
# scp Annotations.gff3.gz sd21@farm5-head1:~/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/SM_V9_2/SM_V9_2_v9.apollo.gff3.gz
gunzip -f ${CHROMOSOME}_v9.apollo.gff3.gz


# step 1: fix gene ID and relationship in mRNA ID
awk -F'[\t]' '$3=="gene" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=..........|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_gene_IDs_NAMEs.txt


# step 2: fix mRNA IDs and all descendant IDs
awk -F'[\t]' '$3=="mRNA" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=............|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_mRNA_IDs_NAMEs.txt


# get product IDs
awk -F'[\t]'  '$3=="mRNA" {print $9}' ${CHROMOSOME}_v9.original.gff3 | sed 's/;/\t/g' | awk -F '[\t]' '{print $1,$3}' OFS="\t" | sed -e 's/ID=/Name=/g' | awk -F '[\t]' '{print $1,$1";"$2}' OFS="\t" > original_IDs_product.txt


cat apollo_gene_IDs_NAMEs.txt apollo_mRNA_IDs_NAMEs.txt > apollo_IDs_NAMEs.txt

# replace apollo unique ID with Smp ID in apollo GFF
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.renamed.gff3 apollo_IDs_NAMEs.txt ${CHROMOSOME}_v9.apollo.gff3


# split annotation into mRNA and other, so that product id's can be added back to mRNAs specifically. Initial testing showed that some CDSs were getting the product IDs added incorrectly, so to account for this, splitting them
awk '$3=="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp
awk '$3!="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp


# add product descriptions to mRNAs only
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 original_IDs_product.txt ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp


# bring the updated mRNAs and other annotations back together.
cat ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp | sort -k4,4n > ${CHROMOSOME}_v9.apollo.renamed2.gff3



# remove overlapping models which have been updated
bedtools intersect -s -v  -f 0.1 -a ${CHROMOSOME}_v9.original.gff3 -b ${CHROMOSOME}_v9.apollo.renamed2.gff3 > ${CHROMOSOME}_v9.original.gff3.no-overlaps


# remove updated gene IDs from original annotation. Might not be strictly necessary, but there were a couple of instances where there were old and new annotations kept.
while read ID GENE; do
     sed -i "/Liftoff.*${GENE}/d" ${CHROMOSOME}_v9.original.gff3.no-overlaps;
done < apollo_gene_IDs_NAMEs.txt



# bring apollo and unchanged original annotations back together.
cat ${CHROMOSOME}_v9.apollo.renamed2.gff3 ${CHROMOSOME}_v9.original.gff3.no-overlaps | grep -v "#" | sort -k4,4n > ${CHROMOSOME}_v9.updated.gff3



# manual check to see if there are duplicated gene IDs
cat ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="gene") print}' | sed -e 's/owner=irisadmin@local.host;//g' -e 's/owner=mb4@sanger.ac.uk//g' -e 's/owner=irisadmin@local.host,mb4@sanger.ac.uk//g' -e 's/,irisadmin@local.host//g' | sed 's/;/\t/g' | cut -f9 | sort | uniq -c | sort



# Chromosome 2 manual checks
# single duplicated ID
   2 ID=Smp_346010 - fixed, made new gene ID


# checking for mini-introns placed in the gff in apollo to fix a gene model, but reflect an error in the genome that needs fixing
gt gff3 -tidy -addintrons -retainids ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="intron") print $9,$5-$4}' | sort -k2,2nr | grep -v "warning"

# removing old Liftoff versions
sed -i "/Liftoff.*Smp_148390/d" SM_V9_2_v9.updated.gff3

# deleting rubbish genes
sed -i "/Smp_327345/d" SM_V9_2_v9.updated.gff3
sed -i "/Smp_340640/d" SM_V9_2_v9.updated.gff3
sed -i "/Smp_327390/d" SM_V9_2_v9.updated.gff3



# mini introns - gene ID, length
Parent=Smp_045410.1 1
Parent=Smp_045410.1 1
Parent=Smp_045410.2 1
Parent=Smp_045410.2 1
Parent=Smp_065290.1 1
Parent=Smp_130050.1 1
Parent=Smp_147330.1 1
Parent=Smp_159180.1 1
Parent=Smp_184350.3 1
Parent=Smp_242220.1 1
Parent=Smp_301160.1 1
Parent=Smp_311770.1 1
Parent=Smp_316380.1 1
Parent=Smp_343180.1 1
Parent=Smp_346650.1 1
Parent=Smp_346650.1 1
Parent=Smp_346650.1 1
Parent=Smp_045430.1 0
Parent=Smp_045430.2 0

```

### fix chromosome 2
```bash
FASTA=SM_V9_2.fa
GFF=SM_V9_2_v9.updated.gff3


fastaq to_fasta -l0 ${FASTA} FASTA.repaired.tmp
cp ${GFF} GFF.repaired.tmp


# read in the list of genes. Note, list needs to be reverse sorted by position to make sure the corrections dont interfere with each other.
sed 1d repair_list.txt | sort -k3rn | while read -r GENE CHROMOSOME POSITION SIZE OLD_SEQUENCE NEW_SEQUENCE; do

# check if repair site is unique
if [ "$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)" -eq 1 ]; then
     # introduce the new sequence to the reference
     sed -i "s/${OLD_SEQUENCE}/${NEW_SEQUENCE}/" FASTA.repaired.tmp

     # update the coordinates in the GFF
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($4>POSITION) print $1,$2,$3,$4+SIZE,$5,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" GFF.repaired.tmp |
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($5>POSITION) print $1,$2,$3,$4,$5+SIZE,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" > GFF.tmp; mv GFF.tmp GFF.repaired.tmp

     echo "${GENE} target site ${POSITION} has been fixed, with a indel of ${SIZE} bp."

else
     COUNT=$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)
     echo "${GENE} target site ${POSITION} is not unique - it was found ${COUNT} times. Check and improve the target sequence so it is unique and try again."
fi;

done

mv FASTA.repaired.tmp ${FASTA%.fa.*$}.repaired.fa
mv GFF.repaired.tmp ${GFF%.gff.*$}.repaired.gff3
rm *.tmp

```

### update from V9 to V10
```bash
# fix chromosome names in fasta
sed 's/SM_V9_2/SM_V10_2/g' SM_V9_2.fa.repaired.fa > SM_V10_2.fa
samtools faidx SM_V10_2.fa

# fix chromosome name in gff, and use gffread to close indels
sed 's/SM_V9_2/SM_V10_2/g' SM_V9_2_v9.updated.gff3.repaired.gff3 | gffread - -Z -F -O -o SM_V10_2.gff3

gffread SM_V10_2.gff3 -g SM_V10_2.fa -y SM_V10_2.proteins.fa

```












## Chromosome 3
### merge original and updated annotations from apollo

```bash  
GENOME=SM_V9_21Feb.fa
CHROMOSOME=SM_V9_3

mkdir /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}


samtools faidx ../${GENOME} ${CHROMOSOME} > ${CHROMOSOME}.fa


# get orignal annotations for chromosome
grep "${CHROMOSOME}" ../../SM_V9_16Mar.gff > ${CHROMOSOME}_v9.original.gff3


# apollo annotations saved as "SM_V9_2_v9.apollo.gff3.gz"
# scp Annotations.gff3.gz sd21@farm5-head1:~/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/SM_V9_3/SM_V9_3_v9.apollo.gff3.gz
gunzip -f ${CHROMOSOME}_v9.apollo.gff3.gz





# step 1: fix gene ID and relationship in mRNA ID
awk -F'[\t]' '$3=="gene" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=..........|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_gene_IDs_NAMEs.txt


# step 2: fix mRNA IDs and all descendant IDs
awk -F'[\t]' '$3=="mRNA" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=............|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_mRNA_IDs_NAMEs.txt


# get product IDs
awk -F'[\t]'  '$3=="mRNA" {print $9}' ${CHROMOSOME}_v9.original.gff3 | sed 's/;/\t/g' | awk -F '[\t]' '{print $1,$3}' OFS="\t" | sed -e 's/ID=/Name=/g' | awk -F '[\t]' '{print $1,$1";"$2}' OFS="\t" > original_IDs_product.txt


cat apollo_gene_IDs_NAMEs.txt apollo_mRNA_IDs_NAMEs.txt > apollo_IDs_NAMEs.txt

# replace apollo unique ID with Smp ID in apollo GFF
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.renamed.gff3 apollo_IDs_NAMEs.txt ${CHROMOSOME}_v9.apollo.gff3


# split annotation into mRNA and other, so that product id's can be added back to mRNAs specifically. Initial testing showed that some CDSs were getting the product IDs added incorrectly, so to account for this, splitting them
awk '$3=="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp
awk '$3!="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp


# add product descriptions to mRNAs only
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 original_IDs_product.txt ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp


# bring the updated mRNAs and other annotations back together.
cat ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp | sort -k4,4n > ${CHROMOSOME}_v9.apollo.renamed2.gff3



# remove overlapping models which have been updated
bedtools intersect -s -v  -f 0.1 -a ${CHROMOSOME}_v9.original.gff3 -b ${CHROMOSOME}_v9.apollo.renamed2.gff3 > ${CHROMOSOME}_v9.original.gff3.no-overlaps


# remove updated gene IDs from original annotation. Might not be strictly necessary, but there were a couple of instances where there were old and new annotations kept.
while read ID GENE; do
     sed -i "/Liftoff.*${GENE}/d" ${CHROMOSOME}_v9.original.gff3.no-overlaps;
done < apollo_gene_IDs_NAMEs.txt



# bring apollo and unchanged original annotations back together.
cat ${CHROMOSOME}_v9.apollo.renamed2.gff3 ${CHROMOSOME}_v9.original.gff3.no-overlaps | grep -v "#" | sort -k4,4n > ${CHROMOSOME}_v9.updated.gff3



# manual check to see if there are duplicated gene IDs
cat ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="gene") print}' | sed -e 's/owner=irisadmin@local.host;//g' -e 's/owner=mb4@sanger.ac.uk//g' -e 's/owner=irisadmin@local.host,mb4@sanger.ac.uk//g' -e 's/,irisadmin@local.host//g' | sed 's/;/\t/g' | cut -f9 | sort | uniq -c | sort


# Chromosome 2 manual checks
# single duplicated ID
2 ID=Smp_350180 - fixed - created new gene ID for one copy

# removing old Liftoff versions
sed -i "/Liftoff.*Smp_326780/d" SM_V9_3_v9.updated.gff3


# checking for mini-introns placed in the gff in apollo to fix a gene model, but reflect an error in the genome that needs fixing
gt gff3 -tidy -addintrons -retainids ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="intron") print $9,$5-$4}' | sort -k2,2nr | grep -v "warning"

# mini introns - gene ID, length
Parent=Smp_185180.1 3
Parent=Smp_074050.1 1
Parent=Smp_104500.1 1
Parent=Smp_104500.2 1
Parent=Smp_104500.3 1
Parent=Smp_104500.4 1
Parent=Smp_134910.1 1
Parent=Smp_335460.1 1
Parent=Smp_335480.1 1
Parent=Smp_347860.1 1
Parent=Smp_349810.1 1
Parent=Smp_317020.1 0
Parent=Smp_336790.1 0
```

### repair indels
```bash
FASTA=SM_V9_3.fa
GFF=SM_V9_3_v9.updated.gff3


fastaq to_fasta -l0 ${FASTA} FASTA.repaired.tmp
cp ${GFF} GFF.repaired.tmp


# read in the list of genes. Note, list needs to be reverse sorted by position to make sure the corrections dont interfere with each other.
sed 1d repair_list.txt | sort -k3rn | while read -r GENE CHROMOSOME POSITION SIZE OLD_SEQUENCE NEW_SEQUENCE; do

# check if repair site is unique
if [ "$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)" -eq 1 ]; then
     # introduce the new sequence to the reference
     sed -i "s/${OLD_SEQUENCE}/${NEW_SEQUENCE}/" FASTA.repaired.tmp

     # update the coordinates in the GFF
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($4>POSITION) print $1,$2,$3,$4+SIZE,$5,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" GFF.repaired.tmp |
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($5>POSITION) print $1,$2,$3,$4,$5+SIZE,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" > GFF.tmp; mv GFF.tmp GFF.repaired.tmp

     echo "${GENE} target site ${POSITION} has been fixed, with a indel of ${SIZE} bp."

else
     COUNT=$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)
     echo "${GENE} target site ${POSITION} is not unique - it was found ${COUNT} times. Check and improve the target sequence so it is unique and try again."
fi;

done

mv FASTA.repaired.tmp ${FASTA%.fa.*$}.repaired.fa
mv GFF.repaired.tmp ${GFF%.gff.*$}.repaired.gff3
rm *.tmp*

```

### update from V9 to V10
```bash
# update chromosome name
sed 's/SM_V9_3/SM_V10_3/g' SM_V9_3.fa.repaired.fa > SM_V10_3.fa
samtools faidx SM_V10_3.fa

# fix chromosome name, and use gffread to close indels
sed 's/SM_V9_3/SM_V10_3/g' SM_V9_3_v9.updated.gff3.repaired.gff3 | gffread - -Z -F -O -o SM_V10_3.gff3

gffread SM_V10_3.gff3 -g SM_V10_3.fa -y SM_V10_3.proteins.fa

```







## chromosome 4
### merge original and updated annotations from apollo
```bash
GENOME=SM_V9_21Feb.fa
CHROMOSOME=SM_V9_4

mkdir /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}


samtools faidx ../${GENOME} ${CHROMOSOME} > ${CHROMOSOME}.fa


# get orignal annotations for chromosome
grep "${CHROMOSOME}" ../../SM_V9_16Mar.gff > ${CHROMOSOME}_v9.original.gff3


# apollo annotations saved as "SM_V9_2_v9.apollo.gff3.gz"
# scp Annotations.gff3.gz sd21@farm5-head1:~/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/SM_V9_4/SM_V9_4_v9.apollo.gff3.gz
gunzip -f ${CHROMOSOME}_v9.apollo.gff3.gz





# step 1: fix gene ID and relationship in mRNA ID
awk -F'[\t]' '$3=="gene" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=..........|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_gene_IDs_NAMEs.txt


# step 2: fix mRNA IDs and all descendant IDs
awk -F'[\t]' '$3=="mRNA" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=............|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_mRNA_IDs_NAMEs.txt


# get product IDs
awk -F'[\t]'  '$3=="mRNA" {print $9}' ${CHROMOSOME}_v9.original.gff3 | sed 's/;/\t/g' | awk -F '[\t]' '{print $1,$3}' OFS="\t" | sed -e 's/ID=/Name=/g' | awk -F '[\t]' '{print $1,$1";"$2}' OFS="\t" > original_IDs_product.txt


cat apollo_gene_IDs_NAMEs.txt apollo_mRNA_IDs_NAMEs.txt > apollo_IDs_NAMEs.txt

# replace apollo unique ID with Smp ID in apollo GFF
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.renamed.gff3 apollo_IDs_NAMEs.txt ${CHROMOSOME}_v9.apollo.gff3


# split annotation into mRNA and other, so that product id's can be added back to mRNAs specifically. Initial testing showed that some CDSs were getting the product IDs added incorrectly, so to account for this, splitting them
awk '$3=="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp
awk '$3!="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp


# add product descriptions to mRNAs only
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 original_IDs_product.txt ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp


# bring the updated mRNAs and other annotations back together.
cat ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp | sort -k4,4n > ${CHROMOSOME}_v9.apollo.renamed2.gff3



# remove overlapping models which have been updated
bedtools intersect -s -v  -f 0.1 -a ${CHROMOSOME}_v9.original.gff3 -b ${CHROMOSOME}_v9.apollo.renamed2.gff3 > ${CHROMOSOME}_v9.original.gff3.no-overlaps


# remove updated gene IDs from original annotation. Might not be strictly necessary, but there were a couple of instances where there were old and new annotations kept.
while read ID GENE; do
     sed -i "/Liftoff.*${GENE}/d" ${CHROMOSOME}_v9.original.gff3.no-overlaps;
done < apollo_gene_IDs_NAMEs.txt



# bring apollo and unchanged original annotations back together.
cat ${CHROMOSOME}_v9.apollo.renamed2.gff3 ${CHROMOSOME}_v9.original.gff3.no-overlaps | grep -v "#" | sort -k4,4n > ${CHROMOSOME}_v9.updated.gff3



# manual check to see if there are duplicated gene IDs
cat ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="gene") print}' | sed -e 's/owner=irisadmin@local.host;//g' -e 's/owner=mb4@sanger.ac.uk//g' -e 's/owner=irisadmin@local.host,mb4@sanger.ac.uk//g' -e 's/,irisadmin@local.host//g' | sed 's/;/\t/g' | cut -f9 | sort | uniq -c | sort

#> no gene duplications

# removing old Liftoff versions
sed -i "/Liftoff.*Smp_312830/d" SM_V9_4_v9.updated.gff3
sed -i "/Liftoff.*Smp_344260/d" SM_V9_4_v9.updated.gff3
sed -i "/Liftoff.*Smp_328910/d" SM_V9_4_v9.updated.gff3
sed -i "/Liftoff.*Smp_192040/d" SM_V9_4_v9.updated.gff3




# checking for mini-introns placed in the gff in apollo to fix a gene model, but reflect an error in the genome that needs fixing
gt gff3 -tidy -addintrons -retainids ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="intron") print $9,$5-$4}' | sort -k2,2nr | grep -v "warning"

# mini introns - gene ID, length
Parent=Smp_120640.1 4
Parent=Smp_043300.1 3
Parent=Smp_062630.1 1
Parent=Smp_120640.1 1
Parent=Smp_149130.1 1
Parent=Smp_192040.1 1
Parent=Smp_312920.1 1
Parent=Smp_349730.1 1
Parent=Smp_349740.1 1
```

### repair indels
```bash
FASTA=SM_V9_4.fa
GFF=SM_V9_4_v9.updated.gff3


fastaq to_fasta -l0 ${FASTA} FASTA.repaired.tmp
cp ${GFF} GFF.repaired.tmp


# read in the list of genes. Note, list needs to be reverse sorted by position to make sure the corrections dont interfere with each other.
sed 1d repair_list.txt | sort -k3rn | while read -r GENE CHROMOSOME POSITION SIZE OLD_SEQUENCE NEW_SEQUENCE; do

# check if repair site is unique
if [ "$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)" -eq 1 ]; then
     # introduce the new sequence to the reference
     sed -i "s/${OLD_SEQUENCE}/${NEW_SEQUENCE}/" FASTA.repaired.tmp

     # update the coordinates in the GFF
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($4>POSITION) print $1,$2,$3,$4+SIZE,$5,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" GFF.repaired.tmp |
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($5>POSITION) print $1,$2,$3,$4,$5+SIZE,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" > GFF.tmp; mv GFF.tmp GFF.repaired.tmp

     echo "${GENE} target site ${POSITION} has been fixed, with a indel of ${SIZE} bp."

else
     COUNT=$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)
     echo "${GENE} target site ${POSITION} is not unique - it was found ${COUNT} times. Check and improve the target sequence so it is unique and try again."
fi;

done

mv FASTA.repaired.tmp ${FASTA%.fa.*$}.repaired.fa
mv GFF.repaired.tmp ${GFF%.gff.*$}.repaired.gff3
rm *.tmp*

```

### update from V9 to V10
```bash
# update chromosome name
sed 's/SM_V9_4/SM_V10_4/g' SM_V9_4.fa.repaired.fa > SM_V10_4.fa
samtools faidx SM_V10_4.fa

# fix chromosome name, and use gffread to close indels
sed 's/SM_V9_4/SM_V10_4/g' SM_V9_4_v9.updated.gff3.repaired.gff3 | gffread - -Z -F -O -o SM_V10_4.gff3

gffread SM_V10_4.gff3 -g SM_V10_4.fa -y SM_V10_4.proteins.fa

```






## Chromosome 5
### merge original and updated annotations from apollo
```bash
GENOME=SM_V9_21Feb.fa
CHROMOSOME=SM_V9_5

mkdir /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}


samtools faidx ../${GENOME} ${CHROMOSOME} > ${CHROMOSOME}.fa


# get orignal annotations for chromosome
grep "${CHROMOSOME}" ../../SM_V9_16Mar.gff > ${CHROMOSOME}_v9.original.gff3


# apollo annotations saved as "SM_V9_2_v9.apollo.gff3.gz"
# scp Annotations.gff3.gz sd21@farm5-head1:~/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/SM_V9_5/SM_V9_5_v9.apollo.gff3.gz
gunzip -f ${CHROMOSOME}_v9.apollo.gff3.gz





# step 1: fix gene ID and relationship in mRNA ID
awk -F'[\t]' '$3=="gene" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=..........|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_gene_IDs_NAMEs.txt


# step 2: fix mRNA IDs and all descendant IDs
awk -F'[\t]' '$3=="mRNA" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=............|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_mRNA_IDs_NAMEs.txt


# get product IDs
awk -F'[\t]'  '$3=="mRNA" {print $9}' ${CHROMOSOME}_v9.original.gff3 | sed 's/;/\t/g' | awk -F '[\t]' '{print $1,$3}' OFS="\t" | sed -e 's/ID=/Name=/g' | awk -F '[\t]' '{print $1,$1";"$2}' OFS="\t" > original_IDs_product.txt


cat apollo_gene_IDs_NAMEs.txt apollo_mRNA_IDs_NAMEs.txt > apollo_IDs_NAMEs.txt

# replace apollo unique ID with Smp ID in apollo GFF
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.renamed.gff3 apollo_IDs_NAMEs.txt ${CHROMOSOME}_v9.apollo.gff3


# split annotation into mRNA and other, so that product id's can be added back to mRNAs specifically. Initial testing showed that some CDSs were getting the product IDs added incorrectly, so to account for this, splitting them
awk '$3=="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp
awk '$3!="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp


# add product descriptions to mRNAs only
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 original_IDs_product.txt ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp


# bring the updated mRNAs and other annotations back together.
cat ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp | sort -k4,4n > ${CHROMOSOME}_v9.apollo.renamed2.gff3



# remove overlapping models which have been updated
bedtools intersect -s -v  -f 0.1 -a ${CHROMOSOME}_v9.original.gff3 -b ${CHROMOSOME}_v9.apollo.renamed2.gff3 > ${CHROMOSOME}_v9.original.gff3.no-overlaps


# remove updated gene IDs from original annotation. Might not be strictly necessary, but there were a couple of instances where there were old and new annotations kept.
while read ID GENE; do
     sed -i "/Liftoff.*${GENE}/d" ${CHROMOSOME}_v9.original.gff3.no-overlaps;
done < apollo_gene_IDs_NAMEs.txt



# bring apollo and unchanged original annotations back together.
cat ${CHROMOSOME}_v9.apollo.renamed2.gff3 ${CHROMOSOME}_v9.original.gff3.no-overlaps | grep -v "#" | sort -k4,4n > ${CHROMOSOME}_v9.updated.gff3



# manual check to see if there are duplicated gene IDs
cat ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="gene") print}' | sed -e 's/owner=irisadmin@local.host;//g' -e 's/owner=mb4@sanger.ac.uk//g' -e 's/owner=irisadmin@local.host,mb4@sanger.ac.uk//g' -e 's/,irisadmin@local.host//g' | sed 's/;/\t/g' | cut -f9 | sort | uniq -c | sort

#> no gene duplications

# checking for mini-introns placed in the gff in apollo to fix a gene model, but reflect an error in the genome that needs fixing
gt gff3 -tidy -addintrons -retainids ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="intron") print $9,$5-$4}' | sort -k2,2nr | grep -v "warning"


# removing old Liftoff versions
sed -i "/Liftoff.*Smp_329930/d" SM_V9_5_v9.updated.gff3
sed -i "/Liftoff.*Smp_332130/d" SM_V9_5_v9.updated.gff3

# mini introns - gene ID, length
Parent=Smp_330280.1 1
Parent=Smp_348380.1 1
Parent=Smp_349200.1 0
Parent=Smp_349290.1 0

```

### repair indels
```bash
FASTA=SM_V9_5.fa
GFF=SM_V9_5_v9.updated.gff3


fastaq to_fasta -l0 ${FASTA} FASTA.repaired.tmp
cp ${GFF} GFF.repaired.tmp


# read in the list of genes. Note, list needs to be reverse sorted by position to make sure the corrections dont interfere with each other.
sed 1d repair_list.txt | sort -k3rn | while read -r GENE CHROMOSOME POSITION SIZE OLD_SEQUENCE NEW_SEQUENCE; do

# check if repair site is unique
if [ "$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)" -eq 1 ]; then
     # introduce the new sequence to the reference
     sed -i "s/${OLD_SEQUENCE}/${NEW_SEQUENCE}/" FASTA.repaired.tmp

     # update the coordinates in the GFF
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($4>POSITION) print $1,$2,$3,$4+SIZE,$5,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" GFF.repaired.tmp |
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($5>POSITION) print $1,$2,$3,$4,$5+SIZE,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" > GFF.tmp; mv GFF.tmp GFF.repaired.tmp

     echo "${GENE} target site ${POSITION} has been fixed, with a indel of ${SIZE} bp."

else
     COUNT=$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)
     echo "${GENE} target site ${POSITION} is not unique - it was found ${COUNT} times. Check and improve the target sequence so it is unique and try again."
fi;

done

mv FASTA.repaired.tmp ${FASTA%.fa.*$}.repaired.fa
mv GFF.repaired.tmp ${GFF%.gff.*$}.repaired.gff3
rm *.tmp*

```

### update from V9 to V10
```bash
# update chromosome name
sed 's/SM_V9_5/SM_V10_5/g' SM_V9_5.fa.repaired.fa > SM_V10_5.fa
samtools faidx SM_V10_5.fa

# fix chromosome name, and use gffread to close indels
sed 's/SM_V9_5/SM_V10_5/g' SM_V9_5_v9.updated.gff3.repaired.gff3 | gffread - -Z -F -O -o SM_V10_5.gff3

gffread SM_V10_5.gff3 -g SM_V10_5.fa -y SM_V10_5.proteins.fa

```









## Chromosome 6
### merge original and updated annotations from apollo
```bash
GENOME=SM_V9_21Feb.fa
CHROMOSOME=SM_V9_6

mkdir /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}


samtools faidx ../${GENOME} ${CHROMOSOME} > ${CHROMOSOME}.fa


# get orignal annotations for chromosome
grep "${CHROMOSOME}" ../../SM_V9_16Mar.gff > ${CHROMOSOME}_v9.original.gff3


# apollo annotations saved as "SM_V9_2_v9.apollo.gff3.gz"
# scp Annotations.gff3.gz sd21@farm5-head1:~/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/SM_V9_6/SM_V9_6_v9.apollo.gff3.gz
gunzip -f ${CHROMOSOME}_v9.apollo.gff3.gz





# step 1: fix gene ID and relationship in mRNA ID
awk -F'[\t]' '$3=="gene" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=..........|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_gene_IDs_NAMEs.txt


# step 2: fix mRNA IDs and all descendant IDs
awk -F'[\t]' '$3=="mRNA" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=............|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_mRNA_IDs_NAMEs.txt


# get product IDs
awk -F'[\t]'  '$3=="mRNA" {print $9}' ${CHROMOSOME}_v9.original.gff3 | sed 's/;/\t/g' | awk -F '[\t]' '{print $1,$3}' OFS="\t" | sed -e 's/ID=/Name=/g' | awk -F '[\t]' '{print $1,$1";"$2}' OFS="\t" > original_IDs_product.txt


cat apollo_gene_IDs_NAMEs.txt apollo_mRNA_IDs_NAMEs.txt > apollo_IDs_NAMEs.txt

# replace apollo unique ID with Smp ID in apollo GFF
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.renamed.gff3 apollo_IDs_NAMEs.txt ${CHROMOSOME}_v9.apollo.gff3


# split annotation into mRNA and other, so that product id's can be added back to mRNAs specifically. Initial testing showed that some CDSs were getting the product IDs added incorrectly, so to account for this, splitting them
awk '$3=="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp
awk '$3!="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp


# add product descriptions to mRNAs only
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 original_IDs_product.txt ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp


# bring the updated mRNAs and other annotations back together.
cat ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp | sort -k4,4n > ${CHROMOSOME}_v9.apollo.renamed2.gff3



# remove overlapping models which have been updated
bedtools intersect -s -v  -f 0.1 -a ${CHROMOSOME}_v9.original.gff3 -b ${CHROMOSOME}_v9.apollo.renamed2.gff3 > ${CHROMOSOME}_v9.original.gff3.no-overlaps


# remove updated gene IDs from original annotation. Might not be strictly necessary, but there were a couple of instances where there were old and new annotations kept.
while read ID GENE; do
     sed -i "/Liftoff.*${GENE}/d" ${CHROMOSOME}_v9.original.gff3.no-overlaps;
done < apollo_gene_IDs_NAMEs.txt



# bring apollo and unchanged original annotations back together.
cat ${CHROMOSOME}_v9.apollo.renamed2.gff3 ${CHROMOSOME}_v9.original.gff3.no-overlaps | grep -v "#" | sort -k4,4n > ${CHROMOSOME}_v9.updated.gff3



# manual check to see if there are duplicated gene IDs
cat ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="gene") print}' | sed -e 's/owner=irisadmin@local.host;//g' -e 's/owner=mb4@sanger.ac.uk//g' -e 's/owner=irisadmin@local.host,mb4@sanger.ac.uk//g' -e 's/,irisadmin@local.host//g' | sed 's/;/\t/g' | cut -f9 | sort | uniq -c | sort

#> single gene duplication
   2 ID=Smp_131350 - fixed - incorrect labelling of another gene

# checking for mini-introns placed in the gff in apollo to fix a gene model, but reflect an error in the genome that needs fixing
gt gff3 -tidy -addintrons -retainids ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="intron") print $9,$5-$4}' | sort -k2,2nr | grep -v "warning"


# mini introns - gene ID, length
Parent=Smp_030860.1 1
Parent=Smp_173330.1 1
Parent=Smp_336880.1 1
Parent=Smp_123260.1 0


```

### repair indels
```bash
FASTA=SM_V9_6.fa
GFF=SM_V9_6_v9.updated.gff3


fastaq to_fasta -l0 ${FASTA} FASTA.repaired.tmp
cp ${GFF} GFF.repaired.tmp


# read in the list of genes. Note, list needs to be reverse sorted by position to make sure the corrections dont interfere with each other.
sed 1d repair_list.txt | sort -k3rn | while read -r GENE CHROMOSOME POSITION SIZE OLD_SEQUENCE NEW_SEQUENCE; do

# check if repair site is unique
if [ "$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)" -eq 1 ]; then
     # introduce the new sequence to the reference
     sed -i "s/${OLD_SEQUENCE}/${NEW_SEQUENCE}/" FASTA.repaired.tmp

     # update the coordinates in the GFF
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($4>POSITION) print $1,$2,$3,$4+SIZE,$5,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" GFF.repaired.tmp |
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($5>POSITION) print $1,$2,$3,$4,$5+SIZE,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" > GFF.tmp; mv GFF.tmp GFF.repaired.tmp

     echo "${GENE} target site ${POSITION} has been fixed, with a indel of ${SIZE} bp."

else
     COUNT=$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)
     echo "${GENE} target site ${POSITION} is not unique - it was found ${COUNT} times. Check and improve the target sequence so it is unique and try again."
fi;

done

mv FASTA.repaired.tmp ${FASTA%.fa.*$}.repaired.fa
mv GFF.repaired.tmp ${GFF%.gff.*$}.repaired.gff3
rm *.tmp*

```

### update from V9 to V10
```bash
# update chromosome name
sed 's/SM_V9_6/SM_V10_6/g' SM_V9_6.fa.repaired.fa > SM_V10_6.fa
samtools faidx SM_V10_6.fa

# fix chromosome name, and use gffread to close indels
sed 's/SM_V9_6/SM_V10_6/g' SM_V9_6_v9.updated.gff3.repaired.gff3 | gffread - -Z -F -O -o SM_V10_6.gff3

gffread SM_V10_6.gff3 -g SM_V10_6.fa -y SM_V10_6.proteins.fa

```





## chromosome 7
```bash
GENOME=SM_V9_21Feb.fa
CHROMOSOME=SM_V9_7

mkdir /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}


samtools faidx ../${GENOME} ${CHROMOSOME} > ${CHROMOSOME}.fa


# get orignal annotations for chromosome
grep "${CHROMOSOME}" ../../SM_V9_16Mar.gff > ${CHROMOSOME}_v9.original.gff3


# apollo annotations saved as "SM_V9_2_v9.apollo.gff3.gz"
# scp Annotations.gff3.gz sd21@farm5-head1:~/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/SM_V9_7/SM_V9_7_v9.apollo.gff3.gz
gunzip -f ${CHROMOSOME}_v9.apollo.gff3.gz


# step 1: fix gene ID and relationship in mRNA ID
awk -F'[\t]' '$3=="gene" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=..........|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_gene_IDs_NAMEs.txt


# step 2: fix mRNA IDs and all descendant IDs
awk -F'[\t]' '$3=="mRNA" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=............|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_mRNA_IDs_NAMEs.txt


# get product IDs
awk -F'[\t]'  '$3=="mRNA" {print $9}' ${CHROMOSOME}_v9.original.gff3 | sed 's/;/\t/g' | awk -F '[\t]' '{print $1,$3}' OFS="\t" | sed -e 's/ID=/Name=/g' | awk -F '[\t]' '{print $1,$1";"$2}' OFS="\t" > original_IDs_product.txt


cat apollo_gene_IDs_NAMEs.txt apollo_mRNA_IDs_NAMEs.txt > apollo_IDs_NAMEs.txt

# replace apollo unique ID with Smp ID in apollo GFF
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.renamed.gff3 apollo_IDs_NAMEs.txt ${CHROMOSOME}_v9.apollo.gff3


# split annotation into mRNA and other, so that product id's can be added back to mRNAs specifically. Initial testing showed that some CDSs were getting the product IDs added incorrectly, so to account for this, splitting them
awk '$3=="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp
awk '$3!="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp


# add product descriptions to mRNAs only
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 original_IDs_product.txt ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp


# bring the updated mRNAs and other annotations back together.
cat ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp | sort -k4,4n > ${CHROMOSOME}_v9.apollo.renamed2.gff3



# remove overlapping models which have been updated
bedtools intersect -s -v  -f 0.1 -a ${CHROMOSOME}_v9.original.gff3 -b ${CHROMOSOME}_v9.apollo.renamed2.gff3 > ${CHROMOSOME}_v9.original.gff3.no-overlaps


# remove updated gene IDs from original annotation. Might not be strictly necessary, but there were a couple of instances where there were old and new annotations kept.
while read ID GENE; do
     sed -i "/Liftoff.*${GENE}/d" ${CHROMOSOME}_v9.original.gff3.no-overlaps;
done < apollo_gene_IDs_NAMEs.txt



# bring apollo and unchanged original annotations back together.
cat ${CHROMOSOME}_v9.apollo.renamed2.gff3 ${CHROMOSOME}_v9.original.gff3.no-overlaps | grep -v "#" | sort -k4,4n > ${CHROMOSOME}_v9.updated.gff3



# manual check to see if there are duplicated gene IDs
cat ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="gene") print}' | sed -e 's/owner=irisadmin@local.host;//g' -e 's/owner=mb4@sanger.ac.uk//g' -e 's/owner=irisadmin@local.host,mb4@sanger.ac.uk//g' -e 's/,irisadmin@local.host//g' | sed 's/;/\t/g' | cut -f9 | sort | uniq -c | sort

#> no gene duplication


# checking for mini-introns placed in the gff in apollo to fix a gene model, but reflect an error in the genome that needs fixing
gt gff3 -tidy -addintrons -retainids ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="intron") print $9,$5-$4}' | sort -k2,2nr | grep -v "warning"

sed -i "/Smp_318440/d" SM_V9_7_v9.updated.gff3


# mini introns - gene ID, length
Parent=Smp_134510.1 9
Parent=Smp_181250.1 4
Parent=Smp_300490.1 1
Parent=Smp_181250.1 0
```

### repair indels
```bash
FASTA=SM_V9_7.fa
GFF=SM_V9_7_v9.updated.gff3


fastaq to_fasta -l0 ${FASTA} FASTA.repaired.tmp
cp ${GFF} GFF.repaired.tmp


# read in the list of genes. Note, list needs to be reverse sorted by position to make sure the corrections dont interfere with each other.
sed 1d repair_list.txt | sort -k3rn | while read -r GENE CHROMOSOME POSITION SIZE OLD_SEQUENCE NEW_SEQUENCE; do

# check if repair site is unique
if [ "$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)" -eq 1 ]; then
     # introduce the new sequence to the reference
     sed -i "s/${OLD_SEQUENCE}/${NEW_SEQUENCE}/" FASTA.repaired.tmp

     # update the coordinates in the GFF
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($4>POSITION) print $1,$2,$3,$4+SIZE,$5,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" GFF.repaired.tmp |
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($5>POSITION) print $1,$2,$3,$4,$5+SIZE,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" > GFF.tmp; mv GFF.tmp GFF.repaired.tmp

     echo "${GENE} target site ${POSITION} has been fixed, with a indel of ${SIZE} bp."

else
     COUNT=$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)
     echo "${GENE} target site ${POSITION} is not unique - it was found ${COUNT} times. Check and improve the target sequence so it is unique and try again."
fi;

done

mv FASTA.repaired.tmp ${FASTA%.fa.*$}.repaired.fa
mv GFF.repaired.tmp ${GFF%.gff.*$}.repaired.gff3
rm *.tmp*

```

### update from V9 to V10
```bash
# update chromosome name
sed 's/SM_V9_7/SM_V10_7/g' SM_V9_7.fa.repaired.fa > SM_V10_7.fa
samtools faidx SM_V10_7.fa

# fix chromosome name, and use gffread to close indels
sed 's/SM_V9_7/SM_V10_7/g' SM_V9_7_v9.updated.gff3.repaired.gff3 | gffread - -Z -F -O -o SM_V10_7.gff3

gffread SM_V10_7.gff3 -g SM_V10_7.fa -y SM_V10_7.proteins.fa

```






## PAR1
```bash
GENOME=SM_V9_21Feb.fa
CHROMOSOME=SM_V9_PAR1

mkdir /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}


samtools faidx ../${GENOME} ${CHROMOSOME} > ${CHROMOSOME}.fa


# get orignal annotations for chromosome
grep "${CHROMOSOME}" ../../SM_V9_16Mar.gff > ${CHROMOSOME}_v9.original.gff3


# apollo annotations saved as "SM_V9_2_v9.apollo.gff3.gz"
# scp Annotations.gff3.gz sd21@farm5-head1:~/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/SM_V9_PAR1/SM_V9_PAR1_v9.apollo.gff3.gz
gunzip -f ${CHROMOSOME}_v9.apollo.gff3.gz


# step 1: fix gene ID and relationship in mRNA ID
awk -F'[\t]' '$3=="gene" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=..........|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_gene_IDs_NAMEs.txt


# step 2: fix mRNA IDs and all descendant IDs
awk -F'[\t]' '$3=="mRNA" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=............|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_mRNA_IDs_NAMEs.txt


# get product IDs
awk -F'[\t]'  '$3=="mRNA" {print $9}' ${CHROMOSOME}_v9.original.gff3 | sed 's/;/\t/g' | awk -F '[\t]' '{print $1,$3}' OFS="\t" | sed -e 's/ID=/Name=/g' | awk -F '[\t]' '{print $1,$1";"$2}' OFS="\t" > original_IDs_product.txt


cat apollo_gene_IDs_NAMEs.txt apollo_mRNA_IDs_NAMEs.txt > apollo_IDs_NAMEs.txt

# replace apollo unique ID with Smp ID in apollo GFF
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.renamed.gff3 apollo_IDs_NAMEs.txt ${CHROMOSOME}_v9.apollo.gff3


# split annotation into mRNA and other, so that product id's can be added back to mRNAs specifically. Initial testing showed that some CDSs were getting the product IDs added incorrectly, so to account for this, splitting them
awk '$3=="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp
awk '$3!="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp


# add product descriptions to mRNAs only
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 original_IDs_product.txt ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp


# bring the updated mRNAs and other annotations back together.
cat ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp | sort -k4,4n > ${CHROMOSOME}_v9.apollo.renamed2.gff3



# remove overlapping models which have been updated
bedtools intersect -s -v  -f 0.1 -a ${CHROMOSOME}_v9.original.gff3 -b ${CHROMOSOME}_v9.apollo.renamed2.gff3 > ${CHROMOSOME}_v9.original.gff3.no-overlaps


# remove updated gene IDs from original annotation. Might not be strictly necessary, but there were a couple of instances where there were old and new annotations kept.
while read ID GENE; do
     sed -i "/Liftoff.*${GENE}/d" ${CHROMOSOME}_v9.original.gff3.no-overlaps;
done < apollo_gene_IDs_NAMEs.txt



# bring apollo and unchanged original annotations back together.
cat ${CHROMOSOME}_v9.apollo.renamed2.gff3 ${CHROMOSOME}_v9.original.gff3.no-overlaps | grep -v "#" | sort -k4,4n > ${CHROMOSOME}_v9.updated.gff3



# manual check to see if there are duplicated gene IDs
cat ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="gene") print}' | sed -e 's/owner=irisadmin@local.host;//g' -e 's/owner=mb4@sanger.ac.uk//g' -e 's/owner=irisadmin@local.host,mb4@sanger.ac.uk//g' -e 's/,irisadmin@local.host//g' | sed 's/;/\t/g' | cut -f9 | sort | uniq -c | sort

#> no gene duplication


# checking for mini-introns placed in the gff in apollo to fix a gene model, but reflect an error in the genome that needs fixing
gt gff3 -tidy -addintrons -retainids ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="intron") print $9,$5-$4}' | sort -k2,2nr | grep -v "warning"

sed -i "/Smp_349120/d" SM_V9_PAR1_v9.updated.gff3


```

### repair indels
```bash
FASTA=SM_V9_PAR1.fa
GFF=SM_V9_PAR1_v9.updated.gff3


fastaq to_fasta -l0 ${FASTA} FASTA.repaired.tmp
cp ${GFF} GFF.repaired.tmp


# read in the list of genes. Note, list needs to be reverse sorted by position to make sure the corrections dont interfere with each other.
sed 1d repair_list.txt | sort -k3rn | while read -r GENE CHROMOSOME POSITION SIZE OLD_SEQUENCE NEW_SEQUENCE; do

# check if repair site is unique
if [ "$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)" -eq 1 ]; then
     # introduce the new sequence to the reference
     sed -i "s/${OLD_SEQUENCE}/${NEW_SEQUENCE}/" FASTA.repaired.tmp

     # update the coordinates in the GFF
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($4>POSITION) print $1,$2,$3,$4+SIZE,$5,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" GFF.repaired.tmp |
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($5>POSITION) print $1,$2,$3,$4,$5+SIZE,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" > GFF.tmp; mv GFF.tmp GFF.repaired.tmp

     echo "${GENE} target site ${POSITION} has been fixed, with a indel of ${SIZE} bp."

else
     COUNT=$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)
     echo "${GENE} target site ${POSITION} is not unique - it was found ${COUNT} times. Check and improve the target sequence so it is unique and try again."
fi;

done

mv FASTA.repaired.tmp ${FASTA%.fa.*$}.repaired.fa
mv GFF.repaired.tmp ${GFF%.gff.*$}.repaired.gff3
rm *.tmp*

```

### update from V9 to V10
```bash
# update chromosome name
sed 's/SM_V9_PAR1/SM_V10_PAR1/g' SM_V9_PAR1.fa.repaired.fa > SM_V10_PAR1.fa
samtools faidx SM_V10_PAR1.fa

# fix chromosome name, and use gffread to close indels
sed 's/SM_V9_PAR1/SM_V10_PAR1/g' SM_V9_PAR1_v9.updated.gff3.repaired.gff3 | gffread - -Z -F -O -o SM_V10_PAR1.gff3

gffread SM_V10_PAR1.gff3 -g SM_V10_PAR1.fa -y SM_V10_PAR1.proteins.fa
```











## PAR2

```bash
GENOME=SM_V9_21Feb.fa
CHROMOSOME=SM_V9_PAR2

mkdir /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}


samtools faidx ../${GENOME} ${CHROMOSOME} > ${CHROMOSOME}.fa


# get orignal annotations for chromosome
grep "${CHROMOSOME}" ../../SM_V9_16Mar.gff > ${CHROMOSOME}_v9.original.gff3


# apollo annotations saved as "SM_V9_2_v9.apollo.gff3.gz"
# scp Annotations.gff3.gz sd21@farm5-head1:~/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/SM_V9_PAR2/SM_V9_PAR2_v9.apollo.gff3.gz
gunzip -f ${CHROMOSOME}_v9.apollo.gff3.gz


# step 1: fix gene ID and relationship in mRNA ID
awk -F'[\t]' '$3=="gene" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=..........|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_gene_IDs_NAMEs.txt


# step 2: fix mRNA IDs and all descendant IDs
awk -F'[\t]' '$3=="mRNA" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=............|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_mRNA_IDs_NAMEs.txt


# get product IDs
awk -F'[\t]'  '$3=="mRNA" {print $9}' ${CHROMOSOME}_v9.original.gff3 | sed 's/;/\t/g' | awk -F '[\t]' '{print $1,$3}' OFS="\t" | sed -e 's/ID=/Name=/g' | awk -F '[\t]' '{print $1,$1";"$2}' OFS="\t" > original_IDs_product.txt


cat apollo_gene_IDs_NAMEs.txt apollo_mRNA_IDs_NAMEs.txt > apollo_IDs_NAMEs.txt

# replace apollo unique ID with Smp ID in apollo GFF
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.renamed.gff3 apollo_IDs_NAMEs.txt ${CHROMOSOME}_v9.apollo.gff3


# split annotation into mRNA and other, so that product id's can be added back to mRNAs specifically. Initial testing showed that some CDSs were getting the product IDs added incorrectly, so to account for this, splitting them
awk '$3=="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp
awk '$3!="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp


# add product descriptions to mRNAs only
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 original_IDs_product.txt ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp


# bring the updated mRNAs and other annotations back together.
cat ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp | sort -k4,4n > ${CHROMOSOME}_v9.apollo.renamed2.gff3



# remove overlapping models which have been updated
bedtools intersect -s -v  -f 0.1 -a ${CHROMOSOME}_v9.original.gff3 -b ${CHROMOSOME}_v9.apollo.renamed2.gff3 > ${CHROMOSOME}_v9.original.gff3.no-overlaps


# remove updated gene IDs from original annotation. Might not be strictly necessary, but there were a couple of instances where there were old and new annotations kept.
while read ID GENE; do
     sed -i "/Liftoff.*${GENE}/d" ${CHROMOSOME}_v9.original.gff3.no-overlaps;
done < apollo_gene_IDs_NAMEs.txt



# bring apollo and unchanged original annotations back together.
cat ${CHROMOSOME}_v9.apollo.renamed2.gff3 ${CHROMOSOME}_v9.original.gff3.no-overlaps | grep -v "#" | sort -k4,4n > ${CHROMOSOME}_v9.updated.gff3



# manual check to see if there are duplicated gene IDs
cat ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="gene") print}' | sed -e 's/owner=irisadmin@local.host;//g' -e 's/owner=mb4@sanger.ac.uk//g' -e 's/owner=irisadmin@local.host,mb4@sanger.ac.uk//g' -e 's/,irisadmin@local.host//g' | sed 's/;/\t/g' | cut -f9 | sort | uniq -c | sort

#> no gene duplication


# checking for mini-introns placed in the gff in apollo to fix a gene model, but reflect an error in the genome that needs fixing
gt gff3 -tidy -addintrons -retainids ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="intron") print $9,$5-$4}' | sort -k2,2nr | grep -v "warning"



Parent=Smp_025370.1 1
Parent=Smp_146440.1 1
Parent=Smp_151285.1 1
Parent=Smp_243800.1 1
Parent=Smp_244080.1 1
Parent=Smp_244320.1 1
```

### repair indels
```bash
FASTA=SM_V9_PAR2.fa
GFF=SM_V9_PAR2_v9.updated.gff3


fastaq to_fasta -l0 ${FASTA} FASTA.repaired.tmp
cp ${GFF} GFF.repaired.tmp


# read in the list of genes. Note, list needs to be reverse sorted by position to make sure the corrections dont interfere with each other.
sed 1d repair_list.txt | sort -k3rn | while read -r GENE CHROMOSOME POSITION SIZE OLD_SEQUENCE NEW_SEQUENCE; do

# check if repair site is unique
if [ "$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)" -eq 1 ]; then
     # introduce the new sequence to the reference
     sed -i "s/${OLD_SEQUENCE}/${NEW_SEQUENCE}/" FASTA.repaired.tmp

     # update the coordinates in the GFF
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($4>POSITION) print $1,$2,$3,$4+SIZE,$5,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" GFF.repaired.tmp |
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($5>POSITION) print $1,$2,$3,$4,$5+SIZE,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" > GFF.tmp; mv GFF.tmp GFF.repaired.tmp

     echo "${GENE} target site ${POSITION} has been fixed, with a indel of ${SIZE} bp."

else
     COUNT=$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)
     echo "${GENE} target site ${POSITION} is not unique - it was found ${COUNT} times. Check and improve the target sequence so it is unique and try again."
fi;

done

mv FASTA.repaired.tmp ${FASTA%.fa.*$}.repaired.fa
mv GFF.repaired.tmp ${GFF%.gff.*$}.repaired.gff3
rm *.tmp*

```

### update from V9 to V10
```bash
# update chromosome name
sed 's/SM_V9_PAR2/SM_V10_PAR2/g' SM_V9_PAR2.fa.repaired.fa > SM_V10_PAR2.fa
samtools faidx SM_V10_PAR2.fa

# fix chromosome name, and use gffread to close indels
sed 's/SM_V9_PAR2/SM_V10_PAR2/g' SM_V9_PAR2_v9.updated.gff3.repaired.gff3 | gffread - -Z -F -O -o SM_V10_PAR2.gff3

gffread SM_V10_PAR2.gff3 -g SM_V10_PAR2.fa -y SM_V10_PAR2.proteins.fa

```






## ZSR
```bash
GENOME=SM_V9_21Feb.fa
CHROMOSOME=SM_V9_ZSR

mkdir /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}


samtools faidx ../${GENOME} ${CHROMOSOME} > ${CHROMOSOME}.fa


# get orignal annotations for chromosome
grep "${CHROMOSOME}" ../../SM_V9_16Mar.gff > ${CHROMOSOME}_v9.original.gff3


# apollo annotations saved as "SM_V9_2_v9.apollo.gff3.gz"
# scp Annotations.gff3.gz sd21@farm5-head1:~/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/SM_V9_ZSR/SM_V9_ZSR_v9.apollo.gff3.gz
gunzip -f ${CHROMOSOME}_v9.apollo.gff3.gz


# step 1: fix gene ID and relationship in mRNA ID
awk -F'[\t]' '$3=="gene" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=..........|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_gene_IDs_NAMEs.txt


# step 2: fix mRNA IDs and all descendant IDs
awk -F'[\t]' '$3=="mRNA" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=............|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_mRNA_IDs_NAMEs.txt


# get product IDs
awk -F'[\t]'  '$3=="mRNA" {print $9}' ${CHROMOSOME}_v9.original.gff3 | sed 's/;/\t/g' | awk -F '[\t]' '{print $1,$3}' OFS="\t" | sed -e 's/ID=/Name=/g' | awk -F '[\t]' '{print $1,$1";"$2}' OFS="\t" > original_IDs_product.txt


cat apollo_gene_IDs_NAMEs.txt apollo_mRNA_IDs_NAMEs.txt > apollo_IDs_NAMEs.txt

# replace apollo unique ID with Smp ID in apollo GFF
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.renamed.gff3 apollo_IDs_NAMEs.txt ${CHROMOSOME}_v9.apollo.gff3


# split annotation into mRNA and other, so that product id's can be added back to mRNAs specifically. Initial testing showed that some CDSs were getting the product IDs added incorrectly, so to account for this, splitting them
awk '$3=="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp
awk '$3!="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp


# add product descriptions to mRNAs only
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 original_IDs_product.txt ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp


# bring the updated mRNAs and other annotations back together.
cat ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp | sort -k4,4n > ${CHROMOSOME}_v9.apollo.renamed2.gff3



# remove overlapping models which have been updated
bedtools intersect -s -v  -f 0.1 -a ${CHROMOSOME}_v9.original.gff3 -b ${CHROMOSOME}_v9.apollo.renamed2.gff3 > ${CHROMOSOME}_v9.original.gff3.no-overlaps


# remove updated gene IDs from original annotation. Might not be strictly necessary, but there were a couple of instances where there were old and new annotations kept.
while read ID GENE; do
     sed -i "/Liftoff.*${GENE}/d" ${CHROMOSOME}_v9.original.gff3.no-overlaps;
done < apollo_gene_IDs_NAMEs.txt



# bring apollo and unchanged original annotations back together.
cat ${CHROMOSOME}_v9.apollo.renamed2.gff3 ${CHROMOSOME}_v9.original.gff3.no-overlaps | grep -v "#" | sort -k4,4n > ${CHROMOSOME}_v9.updated.gff3



# manual check to see if there are duplicated gene IDs
cat ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="gene") print}' | sed -e 's/owner=irisadmin@local.host;//g' -e 's/owner=mb4@sanger.ac.uk//g' -e 's/owner=irisadmin@local.host,mb4@sanger.ac.uk//g' -e 's/,irisadmin@local.host//g' | sed 's/;/\t/g' | cut -f9 | sort | uniq -c | sort

#> two gene duplication
2 ID=Smp_044470 - fixed
2 ID=Smp_304500 - fixed

sed -i "/Liftoff.*Smp_093110/d" SM_V9_ZSR_v9.updated.gff3


sed -i "/Smp_349430/d" SM_V9_ZSR_v9.updated.gff3



# checking for mini-introns placed in the gff in apollo to fix a gene model, but reflect an error in the genome that needs fixing
gt gff3 -tidy -addintrons -retainids ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="intron") print $9,$5-$4}' | sort -k2,2nr | grep -v "warning"


Parent=Smp_343920.1 3
Parent=Smp_120320.1 1
Parent=Smp_242530.1 1
Parent=Smp_342890.1 1
Parent=Smp_134160.1 0
Parent=Smp_142170.1 0
Parent=Smp_193440.1 0
Parent=Smp_203410.1 0
Parent=Smp_346520.1 0

```

### repair indels
```bash
FASTA=SM_V9_ZSR.fa
GFF=SM_V9_ZSR_v9.updated.gff3


fastaq to_fasta -l0 ${FASTA} FASTA.repaired.tmp
cp ${GFF} GFF.repaired.tmp


# read in the list of genes. Note, list needs to be reverse sorted by position to make sure the corrections dont interfere with each other.
sed 1d repair_list.txt | sort -k3rn | while read -r GENE CHROMOSOME POSITION SIZE OLD_SEQUENCE NEW_SEQUENCE; do

# check if repair site is unique
if [ "$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)" -eq 1 ]; then
     # introduce the new sequence to the reference
     sed -i "s/${OLD_SEQUENCE}/${NEW_SEQUENCE}/" FASTA.repaired.tmp

     # update the coordinates in the GFF
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($4>POSITION) print $1,$2,$3,$4+SIZE,$5,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" GFF.repaired.tmp |
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($5>POSITION) print $1,$2,$3,$4,$5+SIZE,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" > GFF.tmp; mv GFF.tmp GFF.repaired.tmp

     echo "${GENE} target site ${POSITION} has been fixed, with a indel of ${SIZE} bp."

else
     COUNT=$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)
     echo "${GENE} target site ${POSITION} is not unique - it was found ${COUNT} times. Check and improve the target sequence so it is unique and try again."
fi;

done

mv FASTA.repaired.tmp ${FASTA%.fa.*$}.repaired.fa
mv GFF.repaired.tmp ${GFF%.gff.*$}.repaired.gff3
rm *.tmp*

```

### update from V9 to V10
```bash
# update chromosome name
sed 's/SM_V9_ZSR/SM_V10_ZSR/g' SM_V9_ZSR.fa.repaired.fa > SM_V10_ZSR.fa
samtools faidx SM_V10_ZSR.fa

# fix chromosome name, and use gffread to close indels
sed 's/SM_V9_ZSR/SM_V10_ZSR/g' SM_V9_ZSR_v9.updated.gff3.repaired.gff3 | gffread - -Z -F -O -o SM_V10_ZSR.gff3

gffread SM_V10_ZSR.gff3 -g SM_V10_ZSR.fa -y SM_V10_ZSR.proteins.fa

```







## WSR

```bash
GENOME=SM_V9_21Feb.fa
CHROMOSOME=SM_V9_WSR

mkdir /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/${CHROMOSOME}


samtools faidx ../${GENOME} ${CHROMOSOME} > ${CHROMOSOME}.fa


# get orignal annotations for chromosome
grep "${CHROMOSOME}" ../../SM_V9_16Mar.gff > ${CHROMOSOME}_v9.original.gff3


# apollo annotations saved as "SM_V9_2_v9.apollo.gff3.gz"
# scp Annotations.gff3.gz sd21@farm5-head1:~/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10/SM_V9_WSR/SM_V9_WSR_v9.apollo.gff3.gz
gunzip -f ${CHROMOSOME}_v9.apollo.gff3.gz


# step 1: fix gene ID and relationship in mRNA ID
awk -F'[\t]' '$3=="gene" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=..........|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_gene_IDs_NAMEs.txt


# step 2: fix mRNA IDs and all descendant IDs
awk -F'[\t]' '$3=="mRNA" {print $9}' OFS="\t" ${CHROMOSOME}_v9.apollo.gff3 | grep -Po "Name=............|ID=...................................." | paste - - | sed -e 's/ID=//g' -e 's/Name=/\t/g' | awk '{print $1,$2}' OFS="\t" > apollo_mRNA_IDs_NAMEs.txt


# get product IDs
awk -F'[\t]'  '$3=="mRNA" {print $9}' ${CHROMOSOME}_v9.original.gff3 | sed 's/;/\t/g' | awk -F '[\t]' '{print $1,$3}' OFS="\t" | sed -e 's/ID=/Name=/g' | awk -F '[\t]' '{print $1,$1";"$2}' OFS="\t" > original_IDs_product.txt


cat apollo_gene_IDs_NAMEs.txt apollo_mRNA_IDs_NAMEs.txt > apollo_IDs_NAMEs.txt

# replace apollo unique ID with Smp ID in apollo GFF
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.renamed.gff3 apollo_IDs_NAMEs.txt ${CHROMOSOME}_v9.apollo.gff3


# split annotation into mRNA and other, so that product id's can be added back to mRNAs specifically. Initial testing showed that some CDSs were getting the product IDs added incorrectly, so to account for this, splitting them
awk '$3=="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp
awk '$3!="mRNA" {print}' ${CHROMOSOME}_v9.apollo.renamed.gff3 > ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp


# add product descriptions to mRNAs only
fsed  --pattern-format=tsv --output ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 original_IDs_product.txt ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp


# bring the updated mRNAs and other annotations back together.
cat ${CHROMOSOME}_v9.apollo.gff3.mRNA.tmp.2 ${CHROMOSOME}_v9.apollo.gff3.no-mRNA.tmp | sort -k4,4n > ${CHROMOSOME}_v9.apollo.renamed2.gff3



# remove overlapping models which have been updated
bedtools intersect -s -v  -f 0.1 -a ${CHROMOSOME}_v9.original.gff3 -b ${CHROMOSOME}_v9.apollo.renamed2.gff3 > ${CHROMOSOME}_v9.original.gff3.no-overlaps


# remove updated gene IDs from original annotation. Might not be strictly necessary, but there were a couple of instances where there were old and new annotations kept.
while read ID GENE; do
     sed -i "/Liftoff.*${GENE}/d" ${CHROMOSOME}_v9.original.gff3.no-overlaps;
done < apollo_gene_IDs_NAMEs.txt



# bring apollo and unchanged original annotations back together.
cat ${CHROMOSOME}_v9.apollo.renamed2.gff3 ${CHROMOSOME}_v9.original.gff3.no-overlaps | grep -v "#" | sort -k4,4n > ${CHROMOSOME}_v9.updated.gff3



# manual check to see if there are duplicated gene IDs
cat ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="gene") print}' | sed -e 's/owner=irisadmin@local.host;//g' -e 's/owner=mb4@sanger.ac.uk//g' -e 's/owner=irisadmin@local.host,mb4@sanger.ac.uk//g' -e 's/,irisadmin@local.host//g' | sed 's/;/\t/g' | cut -f9 | sort | uniq -c | sort

#>  gene duplication


# checking for mini-introns placed in the gff in apollo to fix a gene model, but reflect an error in the genome that needs fixing
gt gff3 -tidy -addintrons -retainids ${CHROMOSOME}_v9.updated.gff3 | awk '{if($3=="intron") print $9,$5-$4}' | sort -k2,2nr | grep -v "warning"

Parent=Smp_318650.1 2
Parent=Smp_348790.1 0
Parent=Smp_348790.2 0
```


### repair indels
```bash
FASTA=SM_V9_WSR.fa
GFF=SM_V9_WSR_v9.updated.gff3


fastaq to_fasta -l0 ${FASTA} FASTA.repaired.tmp
cp ${GFF} GFF.repaired.tmp


# read in the list of genes. Note, list needs to be reverse sorted by position to make sure the corrections dont interfere with each other.
sed 1d repair_list.txt | sort -k3rn | while read -r GENE CHROMOSOME POSITION SIZE OLD_SEQUENCE NEW_SEQUENCE; do

# check if repair site is unique
if [ "$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)" -eq 1 ]; then
     # introduce the new sequence to the reference
     sed -i "s/${OLD_SEQUENCE}/${NEW_SEQUENCE}/" FASTA.repaired.tmp

     # update the coordinates in the GFF
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($4>POSITION) print $1,$2,$3,$4+SIZE,$5,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" GFF.repaired.tmp |
     awk -F '[\t]' -v POSITION=$POSITION -v SIZE=$SIZE '{if($5>POSITION) print $1,$2,$3,$4,$5+SIZE,$6,$7,$8,$9; else print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" > GFF.tmp; mv GFF.tmp GFF.repaired.tmp

     echo "${GENE} target site ${POSITION} has been fixed, with a indel of ${SIZE} bp."

else
     COUNT=$(grep -o "${OLD_SEQUENCE}" FASTA.repaired.tmp | wc -l)
     echo "${GENE} target site ${POSITION} is not unique - it was found ${COUNT} times. Check and improve the target sequence so it is unique and try again."
fi;

done

mv FASTA.repaired.tmp ${FASTA%.fa.*$}.repaired.fa
mv GFF.repaired.tmp ${GFF%.gff.*$}.repaired.gff3
rm *.tmp*

```

### update from V9 to V10
```bash
# update chromosome name
sed 's/SM_V9_WSR/SM_V10_WSR/g' SM_V9_WSR.fa.repaired.fa > SM_V10_WSR.fa
samtools faidx SM_V10_WSR.fa

# fix chromosome name, and use gffread to close indels
sed 's/SM_V9_WSR/SM_V10_WSR/g' SM_V9_WSR_v9.updated.gff3.repaired.gff3 | gffread - -Z -F -O -o SM_V10_WSR.gff3

gffread SM_V10_WSR.gff3 -g SM_V10_WSR.fa -y SM_V10_WSR.proteins.fa

```




## extract V9 gff
```bash  

cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/GENOME_FIX_V10

# extract GFF from individual chromosomes, and parse it via gffread into a collated GFF
cat SM_V9_*/*_v9.updated.gff3 | gffread - -F -O -o SM_V9_updated_221114.gff3

gffread SM_V9_updated_221114.gff3 -g SM_V9_21Feb.fa -y SM_V9_updated_221114.proteins.fa
```

## Fixing Z chromosome
## Merge PAR1, ZSR and PAR2 into a contiguous Z chromosome
- PAR1, ZSR and PAR2 were kept separate in V9 due to the fact that PAR1 and PAR2 are shared with both ZSR and WSR
- however, it is slightly confusing, and is not very useful for downstream uses of the resource for mapping etc
- Matt, Sarah and I concluded that making a contiguous Z would fix this, which we would then describe the PARs as being shared, and keep WSR as a separate sequence

```bash
echo ">SM_V10_Z" > SM_V10_Z.fa.tmp

# extract the sequences in the order in which they will be placed in the chromosome. Note that the orientaiton is already correct
cat SM_V9_PAR1/SM_V10_PAR1.fa > tmp
cat SM_V9_ZSR/SM_V10_ZSR.fa >> tmp
cat SM_V9_PAR2/SM_V10_PAR2.fa >> tmp

# remove the header and whitespaces
grep -v ">" tmp | sed 's/^[[:space:]]*//g' >> SM_V10_Z.fa.tmp

# clean up
fastaq to_fasta SM_V10_Z.fa.tmp SM_V10_Z.fa

rm *tmp

```

##  Merge PAR1, ZSR and PAR2 annotations into into a contiguous Z chromosome annotation
- lengths of the individual chromosomes after correction
     - PAR1 = 10680111
     - ZSR = 33063211
     - PAR2 = 42949106
     - TOTAL = 86692428
- need to join the annotations, using the length of the chromosomes to adjust the annotation positions

```bash
# print the gff columns, adding the chromosome lengths to the ZSR and PAR2 to account for the positional offset
awk -F '[\t]' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS="\t" SM_V9_PAR1/SM_V10_PAR1.gff3 | grep -v "#" > tmp
awk -F '[\t]' '{print $1,$2,$3,$4+10680111,$5+10680111,$6,$7,$8,$9}' OFS="\t" SM_V9_ZSR/SM_V10_ZSR.gff3 | grep -v "#" >> tmp
awk -F '[\t]' '{print $1,$2,$3,$4+43743322,$5+43743322,$6,$7,$8,$9}' OFS="\t" SM_V9_PAR2/SM_V10_PAR2.gff3 | grep -v "#" >> tmp

sort -k1,1 -k 4,4n tmp | gffread - -F -O -o SM_V10_Z.gff3

# change the individual chromosome names to Z
sed -i -e 's/SM_V10_PAR1/SM_V10_Z/g' -e 's/SM_V10_ZSR/SM_V10_Z/g' -e 's/SM_V10_PAR2/SM_V10_Z/g' SM_V10_Z.gff3

gffread SM_V10_Z.gff3 -g SM_V10_Z.fa -y SM_V10_Z.proteins.fa

```

## Bring genome and annotation together
```bash
# genome
cat SM_V9_1/SM_V10_1.fa \
SM_V9_2/SM_V10_2.fa \
SM_V9_3/SM_V10_3.fa \
SM_V9_4/SM_V10_4.fa \
SM_V9_5/SM_V10_5.fa \
SM_V9_6/SM_V10_6.fa \
SM_V9_7/SM_V10_7.fa \
SM_V10_Z.fa \
SM_V9_WSR/SM_V10_WSR.fa > SM_V10.fa

samtools faidx SM_V10.fa

# annotation
echo "##gff-version 3" > SM_V10.gff3
cat SM_V9_1/SM_V10_1.gff3 \
SM_V9_2/SM_V10_2.gff3 \
SM_V9_3/SM_V10_3.gff3 \
SM_V9_4/SM_V10_4.gff3 \
SM_V9_5/SM_V10_5.gff3 \
SM_V9_6/SM_V10_6.gff3 \
SM_V9_7/SM_V10_7.gff3 \
SM_V10_Z.gff3 \
SM_V9_WSR/SM_V10_WSR.gff3 | grep -v "#" >> SM_V10.gff3


|  gffread - -F -O -o SM_V10.gff3


gffread SM_V10.gff3 -g SM_V10.fa -y SM_V10.proteins.fa


```






PROTIENS=*.proteins.fa

fastaq to_fasta -l0 ${PROTIENS} protein.fa.tmp

grep ">" protein.fa.tmp | sed "s/>//g" > protein.ids

>gene_check.txt
while read ID GENE; do \
     grep -A1 "${ID}" protein.fa.tmp > protein.tmp; \
     SEQUENCE=$( tail -n 1 protein.tmp ); \
     NUMBER_STOPS=$( echo ${SEQUENCE} | grep -o "\." | wc -l ); \
     M_START=$(if [[ $SEQUENCE =~ ^M.* ]]; then echo yes; else echo no; fi); \
     echo -e $ID"\t"$M_START"\t"$NUMBER_STOPS >> gene_check.txt ; \
     done < protein.ids


cat SM_V9_*/gene_check.txt | awk '{if($2=="no") print $0}' OFS="\t" > Sm_genes_noM.txt
cat SM_V9_*/gene_check.txt | awk '{if($3=="0") print $0}' OFS="\t" > Sm_genes_noSTOP.txt



cat SM_V9_1/SM_V10_1.fa \
SM_V9_2/SM_V10_2.fa \
SM_V9_3/SM_V10_3.fa \
SM_V9_4/SM_V10_4.fa \
SM_V9_5/SM_V10_5.fa \
SM_V9_6/SM_V10_6.fa \
SM_V9_7/SM_V10_7.fa \
SM_V10_Z.fa \
SM_V9_WSR/SM_V10_WSR.fa > SM_V10.tmp.fa


cat SM_V9_1/SM_V10_1.fa \
SM_V9_2/SM_V10_2.fa \
SM_V9_3/SM_V10_3.fa \
SM_V9_4/SM_V10_4.fa \
SM_V9_5/SM_V10_5.fa \
SM_V9_6/SM_V10_6.fa \
SM_V9_7/SM_V10_7.fa \
SM_V10_Z.fa \
SM_V9_WSR/SM_V10_WSR.fa > SM_V10.fa

samtools faidx SM_V10.fa
#





# annotation
echo "##gff-version 3" > SM_V10.tmp.gff3

cat SM_V9_1/SM_V10_1.gff3 \
SM_V9_2/SM_V10_2.gff3 \
SM_V9_3/SM_V10_3.gff3 \
SM_V9_4/SM_V10_4.gff3 \
SM_V9_5/SM_V10_5.gff3 \
SM_V9_6/SM_V10_6.gff3 \
SM_V9_7/SM_V10_7.gff3 \
SM_V10_Z.gff3 \
SM_V9_WSR/SM_V10_WSR.gff3 | grep -v "#" >> SM_V10.tmp.gff3




echo "##gff-version 3" > SM_V10.gff3

cat SM_V9_1/SM_V10_1.gff3 \
SM_V9_2/SM_V10_2.gff3 \
SM_V9_3/SM_V10_3.gff3 \
SM_V9_4/SM_V10_4.gff3 \
SM_V9_5/SM_V10_5.gff3 \
SM_V9_6/SM_V10_6.gff3 \
SM_V9_7/SM_V10_7.gff3 \
SM_V10_Z.gff3 \
SM_V9_WSR/SM_V10_WSR.gff3 | grep -v "#" >> SM_V10.gff3
#
#
#
#

#



gffread SM_V10.tmp.gff3 -g SM_V10.tmp.fa -y SM_V10.proteins.tmp.fa
gffread SM_V10.gff3 -g SM_V10.fa -y SM_V10.proteins.fa
