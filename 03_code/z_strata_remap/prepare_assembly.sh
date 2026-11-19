#!/usr/bin/env bash
# Download one assembly and place S. mansoni v10 genes on it.
#
# usage: prepare_assembly.sh WORKDIR SPECIES_DIR ASSEMBLY [THREADS]
#
# ASSEMBLY is a WBPS19 bioproject (e.g. TD1_PRJEB44434, fetched from WormBase
# ParaSite as SPECIES_DIR.ASSEMBLY) or an NCBI accession (GCA_..., fetched
# from the NCBI FTP site).
#
# Writes WORKDIR/assemblies/<species_dir>.<assembly>/:
#   genome.fa(.fai)    the assembly
#   loci.bed           CDS intervals of each S. mansoni gene's best placement (name = Smp gene)
#   liftover.tsv       one row per placed gene: locus, identity, coverage, ambiguity flag
#   vs_SmZ.paf         whole-genome alignment to S. mansoni v10 Z (minimap2 asm20), for
#                      finding the species' Z and checking rearrangements at boundaries
#
# For S. mansoni itself the v10 annotation is used directly. Other species get
# the v10 proteins aligned with miniprot; a gene is kept when its best hit
# covers >=50% of the protein at >=40% identity, and flagged ambiguous when a
# second hit is nearly as good (likely paralog or collapsed duplicate).
set -euo pipefail

WORKDIR=$1; SPECIES=$2; ASSEMBLY=$3; THREADS=${4:-16}
WBPS=https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all
REF=${WORKDIR}/ref
OUT=${WORKDIR}/assemblies/${SPECIES}.${ASSEMBLY}
mkdir -p ${OUT}
cd ${OUT}

if [[ ! -s genome.fa.fai ]]; then
    if [[ ${ASSEMBLY} == GCA_* ]]; then
        # GCA_025215515.1 -> GCA/025/215/515/GCA_025215515.1_<assembly name>/
        DIGITS=${ASSEMBLY#GCA_}; DIGITS=${DIGITS%.*}
        DIR=${NCBI}/GCA/${DIGITS:0:3}/${DIGITS:3:3}/${DIGITS:6:3}
        NAME=$(wget -q -O - ${DIR}/ | grep -o "${ASSEMBLY}_[^/\"]*" | head -1)
        URL=${DIR}/${NAME}/${NAME}_genomic.fna.gz
    else
        URL=${WBPS}/${SPECIES}/${ASSEMBLY}/${SPECIES}.${ASSEMBLY}.WBPS19.genomic.fa.gz
    fi
    wget -q -O - ${URL} | gunzip -c | sed '/^>/s/ .*//' > genome.fa
    samtools faidx genome.fa
fi

if [[ ! -s vs_SmZ.paf ]]; then
    minimap2 -x asm20 -t ${THREADS} ${REF}/SM_V10_Z.fa genome.fa > vs_SmZ.paf
fi

if [[ ${SPECIES} == "schistosoma_mansoni" ]]; then
    # CDS of the longest isoform per gene, straight from the annotation
    zcat ${REF}/sman.annotations.gff3.gz |\
        awk -F'\t' -v keep=${REF}/sman_longest_tx.txt '
            BEGIN{while((getline t < keep) > 0) K["transcript:"t]=1}
            $3=="CDS"{match($9,/Parent=[^;]+/); p=substr($9,RSTART+7,RLENGTH-7);
                if(p in K){g=p; sub(/^transcript:/,"",g); sub(/\.[0-9]+$/,"",g); print $1,$4-1,$5,g}}' OFS='\t' |\
        sort -k1,1 -k2,2n > loci.bed
    awk -F'\t' '{print $4, $1, 1, 1, 0}' OFS='\t' loci.bed | sort -u |\
        awk 'BEGIN{print "gene\tcontig\tidentity\tcoverage\tambiguous"}{print}' OFS='\t' > liftover.tsv
    exit 0
fi

if [[ ! -s miniprot.gff ]]; then
    miniprot -t ${THREADS} --gff -N 5 --outs=0.5 genome.fa ${REF}/sman_proteins.fa > miniprot.gff.tmp
    mv miniprot.gff.tmp miniprot.gff
fi

# mRNA lines carry Rank, Identity, Positive and Target=<protein> <start> <end>
awk -F'\t' -v fai=${REF}/sman_proteins.fa.fai '
    BEGIN{OFS="\t"; while((getline l < fai) > 0){split(l,f,"\t"); L[f[1]]=f[2]}}
    $3=="mRNA"{
        id=""; rank=0; ident=0; pos=0; tgt=""; ts=0; te=0
        n=split($9,a,";")
        for(i=1;i<=n;i++){
            split(a[i],kv,"=")
            if(kv[1]=="ID") id=kv[2]
            else if(kv[1]=="Rank") rank=kv[2]
            else if(kv[1]=="Identity") ident=kv[2]
            else if(kv[1]=="Positive") pos=kv[2]
            else if(kv[1]=="Target"){split(kv[2],t," "); tgt=t[1]; ts=t[2]; te=t[3]}
        }
        g=tgt; sub(/\.[0-9]+$/,"",g)
        cov=(te-ts+1)/L[tgt]
        if(rank==1){best[g]=id; ID1[g]=ident; P1[g]=pos; C1[g]=cov; CTG[g]=$1}
        else if(rank==2){P2[g]=pos; C2[g]=cov}
    }
    END{
        print "gene\tcontig\tidentity\tcoverage\tambiguous" > "liftover.tsv"
        for(g in best){
            if(ID1[g] < 0.4 || C1[g] < 0.5) continue
            amb=(g in P2 && P2[g] >= P1[g]-0.02 && C2[g] >= 0.5*C1[g]) ? 1 : 0
            print g, CTG[g], ID1[g], C1[g], amb >> "liftover.tsv"
            print best[g], g > "keep_ids.tmp"
        }
    }' miniprot.gff

awk -F'\t' 'NR==FNR{K[$1]=$2; next}
    $3=="CDS"{match($9,/Parent=[^;]+/); p=substr($9,RSTART+7,RLENGTH-7); if(p in K) print $1,$4-1,$5,K[p]}' OFS='\t' \
    keep_ids.tmp miniprot.gff | sort -k1,1 -k2,2n > loci.bed
rm keep_ids.tmp
