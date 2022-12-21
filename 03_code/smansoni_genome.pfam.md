# Schistosoma mansoni V10 genome: Pfam domain and clustering analysis

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk

- the original V9 paper described an analysis of Pfam domains to determine if there was evidence of clustering of genes 
- here is the update on V10
- note it relies on a script to do the clustering step written by Zhigang Liu


## Get some data
```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/PFAM

ln -s /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/SM_V10.annotation.preWBP18checked.gff3
/nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/SM_V10.annotation.preWBP18checked.fa

gffread -y PROTEINS.fa -g SM_V10.annotation.preWBP18checked.fa SM_V10.annotation.preWBP18checked.gff3

```


## Get Pfam data and annotate Sm proteins using HMMER
```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases//Pfam35.0/Pfam-A.hmm.gz

module load hmmer/3.2.1=hfc679d8_0-c1

gunzip Pfam-A.hmm.gz

hmmpress Pfam-A.hmm

bsub.py 2 hmm2 "hmmsearch --tblout hmm_out_tbl.txt -E 1e-5 --cpu 2 Pfam-A.hmm PROTEINS.fa"
```

## Make input files for zhigang's script
```bash
# 1. gene-chr-start.txt (must list all genes)
cat SM_V10.annotation.preWBP18checked.gff3 | awk -F '[\t;]' '$3=="gene" {print $9, $1, $4}' OFS="\t" | sed 's/ID=//g' > gene-chr-start.txt
#grep "Smp_" gene-chr-start.txt | sort -k2,2 -k3,3n > tmp; mv tmp gene-chr-start.txt
grep "Smp_" gene-chr-start.txt | sort -k1,1 > tmp; mv tmp gene-chr-start.txt

# 2. gene-func.txt (gene and domain ids separated by ,)
>gene-func.txt
while read GENE CHR POS; do
    PFAM_LIST=$(grep "${GENE}" hmm_out_tbl.txt | awk '{print $4}' | cut -c-7 | sed -z 's/\n/,/g;s/,$/\n/')
    echo -e ${GENE}"\t"${PFAM_LIST} >> gene-func.txt;
    done < gene-chr-start.txt

sort -k1,1 gene-func.txt > tmp; mv tmp gene-func.txt

# 3. func-names.txt (domain id and name)
cat hmm_out_tbl.txt | grep -v "#" | awk '{print $4,$3}' OFS="\t" | sort | uniq | awk -F '[.\t]' '{print $1,$3}' OFS="\t" > func-names.txt

# 4. chr-length.txt (chromosome lengths for plotting, seprated by " ")
cut -f1,2 SM_V10.genome.preWBP18checked.fa.fai | sort > chr-length.txt

```
## Run Zhigang's tool: "_functionalClusters.sh" 
```bash
# get Zhigangs code 
git clone https://github.com/zglu/FunctionalClusters_adjacentGenes.git

cd FunctionalClusters_adjacentGenes

ln -f -s ../gene-chr-start.txt
ln -f -s ../gene-func.txt
ln -f -s ../func-names.txt
ln -f -s ../chr-length.txt

./_functionalClusters.sh

```


## Repeat for V5 analysis
```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/PFAM/V5

ln -s /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/V5_vs_V10/schistosoma_mansoni.PRJEA36577.WBPS1.genomic.fa
ln -s /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/V5_vs_V10/schistosoma_mansoni.PRJEA36577.WBPS1.annotations.gff3
ln -s /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/V5_vs_V10/schistosoma_mansoni.PRJEA36577.WBPS1.protein.fa
ln -s ../Pfam-A.hmm

module load hmmer/3.2.1=hfc679d8_0-c1

hmmpress Pfam-A.hmm
bsub.py 2 --threads 10 hmm "hmmsearch --tblout hmm_out_tbl.txt -E 1e-5 --cpu 10 Pfam-A.hmm schistosoma_mansoni.PRJEA36577.WBPS1.protein.fa"

# 1. gene-chr-start.txt (must list all genes)
cat schistosoma_mansoni.PRJEA36577.WBPS1.annotations.gff3 | awk -F '[\t;]' '$3=="gene" {print $9, $1, $4}' OFS="\t" | sed -e 's/gene://g' -e 's/ID=//g' > gene-chr-start.txt
#grep "Smp_" gene-chr-start.txt | sort -k2,2 -k3,3n > tmp; mv tmp gene-chr-start.txt
grep "Smp_" gene-chr-start.txt | sort -k1,1 > tmp; mv tmp gene-chr-start.txt


# 2. gene-func.txt (gene and domain ids separated by ,)
>gene-func.txt
while read GENE CHR POS; do
    PFAM_LIST=$(grep "${GENE}" hmm_out_tbl.txt | awk '{print $4}' | cut -c-7 | sed -z 's/\n/,/g;s/,$/\n/')
    echo -e ${GENE}"\t"${PFAM_LIST} >> gene-func.txt;
    done < gene-chr-start.txt

sort -k1,1 gene-func.txt > tmp; mv tmp gene-func.txt

# 3. func-names.txt (domain id and name)
cat hmm_out_tbl.txt | grep -v "#" | awk '{print $4,$3}' OFS="\t" | sort | uniq | awk -F '[.\t]' '{print $1,$3}' OFS="\t" > func-names.txt



# 4. chr-length.txt (chromosome lengths for plotting, seprated by " ")
samtools faidx schistosoma_mansoni.PRJEA36577.WBPS1.genomic.fa
cut -f1,2 schistosoma_mansoni.PRJEA36577.WBPS1.genomic.fa.fai | sort > chr-length.txt


# get Zhigangs code 
git clone https://github.com/zglu/FunctionalClusters_adjacentGenes.git

cd FunctionalClusters_adjacentGenes

ln -f -s ../gene-chr-start.txt
ln -f -s ../gene-func.txt
ln -f -s ../func-names.txt
ln -f -s ../chr-length.txt

./_functionalClusters.sh
```