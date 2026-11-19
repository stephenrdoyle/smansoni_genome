#!/usr/bin/env bash
# S. mansoni v10 reference files shared by every species.
#
# usage: prepare_ref.sh WORKDIR
#
# Writes WORKDIR/ref/:
#   sman.annotations.gff3.gz
#   sman.genome.fa, SM_V10_Z.fa   v10 genome and its Z (target for vs_SmZ.paf)
#   sman_longest_tx.txt      longest protein-coding isoform per gene
#   sman_proteins.fa(.fai)   proteins of those isoforms (miniprot queries)
#   sman_genes.tsv           gene, chromosome, start, end in v10 coordinates
set -euo pipefail

WORKDIR=$1
WBPS=https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/schistosoma_mansoni/PRJEA36577
mkdir -p ${WORKDIR}/ref
cd ${WORKDIR}/ref

wget -q -O sman.annotations.gff3.gz ${WBPS}/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3.gz
wget -q -O sman.protein.fa.gz ${WBPS}/schistosoma_mansoni.PRJEA36577.WBPS19.protein.fa.gz
wget -q -O - ${WBPS}/schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa.gz | gunzip -c > sman.genome.fa
samtools faidx sman.genome.fa
samtools faidx sman.genome.fa SM_V10_Z > SM_V10_Z.fa

seqkit fx2tab -n -i -l sman.protein.fa.gz |\
    awk -F'\t' '{g=$1; sub(/\.[0-9]+$/,"",g); if($2>L[g]){L[g]=$2; T[g]=$1}} END{for(g in T) print T[g]}' |\
    sort > sman_longest_tx.txt
seqkit grep -f sman_longest_tx.txt sman.protein.fa.gz | seqkit seq -i | sed '/^[^>]/s/\*$//' > sman_proteins.fa
samtools faidx sman_proteins.fa

zcat sman.annotations.gff3.gz |\
    awk -F'\t' '$3=="gene" && $9~/biotype=protein_coding/{match($9,/ID=gene:[^;]+/); print substr($9,RSTART+8,RLENGTH-8),$1,$4,$5}' OFS='\t' |\
    sort -k2,2 -k3,3n |\
    awk 'BEGIN{print "gene\tchrom\tstart\tend"}{print}' OFS='\t' > sman_genes.tsv
