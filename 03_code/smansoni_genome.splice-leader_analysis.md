# splice-leader analysis

# author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk

```
#!/usr/bin/env bash
########################################################################

# SPLICE LEADER (SL) FINDER

########################################################################

# Usage: ~sd21/bash_scripts/run_spliceleader_finder.sh <PREFIX> <REFERENCE> <GFF> <splice leader sequence> <SL_MIN_LENGTH> <R1.fastq> <R2.fastq>

# Requirements in path
#--- samtools-1.6
#--- bedtools v2.17.0
#--- hisat2


# Notes
#--- use bsub to run
#--- relies on a Augustus-like GFF to detect overlap of SL with transcription start sites, and may fail to to the last step if not compatible. Can hack if needed.


# Modifications
#- 180712: added a minimum length of SL parameter for cutadapt to work with - originally set to 10 bp, but probably should be set to ~14 based on Alans script to find SLrelated seq by chance in genome
#- 180712: changed length of upstream window from start codon from 50 to 100 bp when looking for overlap between SL sequence and transcripts


########################################################################

PREFIX=$1
REFERENCE=$2
GFF=$3
SL_SEQ=$4
SL_MIN_LENGTH=$5
R1=$6
R2=$7


if [ "$#" -eq 0 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
	echo ""
    echo "Usage: ~sd21/bash_scripts/run_spliceleader_finder.sh <PREFIX> <REFERENCE> <GFF> <splice leader sequence> <SL_MIN_LENGTH> <R1.fastq> <R2.fastq>"
    echo ""
    exit 0
fi

if [ -d "${PREFIX}_SL_ANALYSIS_out" ]; then
	echo -e "\nThere is already a directory started with this sample name. Rename and start again\n"
    exit 0
fi


mkdir ${PREFIX}_SL_ANALYSIS_out
cd ${PREFIX}_SL_ANALYSIS_out
ln -s ${R1}
ln -s ${R2}
ln -s ${REFERENCE}
ln -s ${GFF}


# cutadapt - http://cutadapt.readthedocs.io/en/stable/

# R1 reads with mates

#
#cutadapt -j 7 \
#	--prefix=${PREFIX}_SL_trimmed_ \
#	--overlap=${SL_MIN_LENGTH} \
#	--error-rate=0.1 \
#	-g ${SL_SEQ} -G ${SL_SEQ} \
#	-o ${PREFIX}.trimmed.1.tmp.fastq \
#	-p ${PREFIX}.trimmed.2.tmp.fastq \
#	${R1} ${R2} >> ${PREFIX}_cutadapt_trim.summary

# R1 reads with mates
cutadapt -j 7 --trimmed-only --prefix=${PREFIX}_SL_trimmed_ --overlap=${SL_MIN_LENGTH} --error-rate=0.1 -g ${SL_SEQ} -o ${PREFIX}.trimmed.1_1.tmp.fastq -p ${PREFIX}.trimmed.1_2.tmp.fastq ${R1} ${R2} > ${PREFIX}_cutadapt_trim.summary


# get the reverse reads by switching R1 and R2 around
cutadapt -j 7 --trimmed-only --prefix=${PREFIX}_SL_trimmed_ --overlap=${SL_MIN_LENGTH} --error-rate=0.1 -g ${SL_SEQ} -o ${PREFIX}.trimmed.2_2.tmp.fastq -p ${PREFIX}.trimmed.2_1.tmp.fastq ${R2} ${R1} >> ${PREFIX}_cutadapt_trim.summary


# merge the two datasets together
cat ${PREFIX}.trimmed.1_1.tmp.fastq ${PREFIX}.trimmed.2_1.tmp.fastq > ${PREFIX}.trimmed.merged_R1.tmp.fastq
cat ${PREFIX}.trimmed.1_2.tmp.fastq ${PREFIX}.trimmed.2_2.tmp.fastq > ${PREFIX}.trimmed.merged_R2.tmp.fastq



# map reads using histat
hisat2-build -p 1 ${REFERENCE} ${PREFIX}.REFidx.tmp
hisat2 ${PREFIX}.REFidx.tmp -1 ${PREFIX}.trimmed.merged_R1.tmp.fastq -2 ${PREFIX}.trimmed.merged_R2.tmp.fastq -q -S ${PREFIX}.tmp.sam


# extract trimmed and mapped SL reads
samtools view ${PREFIX}.tmp.sam -H > ${PREFIX}.tmp.sam.header
grep "${PREFIX}_SL_trimmed_" ${PREFIX}.tmp.sam > ${PREFIX}.tmp.sam.body
cat ${PREFIX}.tmp.sam.header ${PREFIX}.tmp.sam.body > ${PREFIX}.SL-only.tmp.sam
samtools view -bS ${PREFIX}.SL-only.tmp.sam > ${PREFIX}.SL-only.tmp.bam
samtools sort ${PREFIX}.SL-only.tmp.bam -o ${PREFIX}.SL-only.sorted.bam
samtools index -b ${PREFIX}.SL-only.sorted.bam

# count SL reads using feature counts
featureCounts -a ${GFF} -o ${PREFIX}.featurecounts.out -g "ID" -M -O -t exon ${PREFIX}.SL-only.sorted.bam



# for schisto
#- extract first and second exons of genes (if UTR is present, sometimes is in 2nd exon), count if 1 or more SL reads are mapped
grep  -e $'\-1\t' -e $'\-2\t'  sm.featurecounts.out | awk '{if($7>=1) print}' | cut -f1 | cut -c-10 | sort | uniq | wc -l
#> 4271

# check for additonal SL sequences in other exons
grep  -v -e $'\-1\t' -e $'\-2\t'  sm.featurecounts.out | awk '{if($7>=1) print}' | cut -f1 | cut -c-10 | sort | uniq | wc -l
#> 5280

# 5280 - 4271 = 1009 additional genes with an internal SL

cat  sm.featurecounts.out | awk '{if($7>=1) print}' | cut -f1 | cut -c-10 | sort | uniq | wc -l
#> 6909 genes with an SL

# 6909 - 4271 = 2638 additional genes with an internal SL

# summarise SL reads per transcript, summing over exon 1 and 2 whee needed.
grep  -e $'\-1\t' -e $'\-2\t'  sm.featurecounts.out | awk '{if($7>=1) print}' | awk  '{ a = substr($1,1,10);b = substr($1,1,12) ;  print b,a,$2,$3,$4,$5,$7 }' OFS="\t" | datamash --full groupby 1 sum 7 | wc -l
```
