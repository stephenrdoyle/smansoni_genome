# Schistosoma mansoni V10 genome: Splice Leader analysis

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk


```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/SPLICE_LEADER

bsub.py 1 test "bash ./run_spliceleader_finder.sh \
    SM_V10 /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/SM_V10.genome.preWBP18checked.fa \
    /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/SM_V10.annotation.preWBP18checked.gff3 \ AACCGTCACGGTTTTACTCTTGTGATTTGTTGCATG \
    10 \
    $PWD/merged_1.fastq.gz \
    $PWD/merged_2.fastq.gz"
```

where "run_spliceleader_finder.sh" is
```bash
#!/usr/bin/env bash
########################################################################

# SPLICE LEADER (SL) FINDER

########################################################################

# Usage: ~sd21/bash_scripts/run_spliceleader_finder.sh <PREFIX> <REFERENCE> <GFF> <splice leader sequence> <SL_MIN_LENGTH> <R1.fastq> <R2.fastq>

# Requirements in path
#--- samtools-1.6
#--- bedtools v2.17.0
#--- hisat2

eval "$(conda shell.bash hook)"
conda activate cutadaptenv
module load hisat2/2.1.0--py36pl5.22.0_0
module load subread/2.0.1--hed695b0_0
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
featureCounts -a ${GFF} -o ${PREFIX}.featurecounts.out -g "Parent" -M -O -t exon ${PREFIX}.SL-only.sorted.bam

```



## Analysis
- need to identify transcripts in that contain SL sequence at their 5' and genes that have internal SL sequences
```bash
>SM_V10.transcripts_w_5prime_SLs.txt
>SM_V10.transcripts_w_SLs_anywhere.txt
cut -f1 SM_V10.featurecounts.out | cut -f2 -d ":" | sort | uniq | grep "Smp" | while read TRANSCRIPT; do

    grep "${TRANSCRIPT}" SM_V10.featurecounts.out > ${TRANSCRIPT}.data.tmp
    
    cat ${TRANSCRIPT}.data.tmp | 
        STRAND=$( head -n 1 ${TRANSCRIPT}.data.tmp | cut -f5 )
        if [ ${STRAND}=='+' ]; then 
            COUNT=$(head -n2 ${TRANSCRIPT}.data.tmp | datamash sum 7)
            GENE=$(echo ${TRANSCRIPT} | cut -c-10)
            if [ ${COUNT} -ge 1 ]; then
            echo -e ${TRANSCRIPT}"\t"${GENE}"\t"${COUNT} >> SM_V10.transcripts_w_5prime_SLs.txt;
            fi
        else 
            COUNT=$(tail -n2 ${TRANSCRIPT}.data.tmp | datamash sum 7)
            if [ ${COUNT} -ge 1 ]; then
            echo -e ${TRANSCRIPT}"\t"${GENE}"\t"${COUNT} >> SM_V10.transcripts_w_5prime_SLs.txt;
            fi
        fi;

    cat ${TRANSCRIPT}.data.tmp | 
        
        if [ ${STRAND}=='+' ]; then 
            COUNT=$(cat ${TRANSCRIPT}.data.tmp | datamash sum 7)
            GENE=$(echo ${TRANSCRIPT} | cut -c-10 )
            if [ ${COUNT} -ge "1" ]; then
            echo -e ${TRANSCRIPT}"\t"${GENE}"\t"${COUNT} >> SM_V10.transcripts_w_SLs_anywhere.txt;
            fi
        else 
            COUNT=$(cat ${TRANSCRIPT}.data.tmp | datamash sum 7)
            if [ ${COUNT} -ge "1" ]; then
            echo -e ${TRANSCRIPT}"\t"${GENE}"\t"${COUNT} >> SM_V10.transcripts_w_SLs_anywhere.txt;
            fi
        fi;
        rm ${TRANSCRIPT}.data.tmp
    done


# count number of transcripts with 5' SLs
wc -l SM_V10.transcripts_w_5prime_SLs.txt
#> 4793 SM_V10.transcripts_w_5prime_SLs.txt

cat SM_V10.transcripts_w_5prime_SLs.txt | cut -c-10 | sort | uniq | wc -l
#> 4326 genes

# count number of transcripts with SLs anywhere in the gene
wc -l SM_V10.transcripts_w_5prime_SLs.txt
#> 7711 SM_V10.transcripts_w_SLs_anywhere.txt

cat SM_V10.transcripts_w_SLs_anywhere.txt | cut -c-10 | sort | uniq | wc -l
#> 6897 genes

#> 6897 - 4326 = 2571 with cryptic SLs internally

# percentages
cut -f1 SM_V10.featurecounts.out | cut -f2 -d ":" | sort | uniq | grep "Smp" | cut -c-10 | sort | uniq | wc -l
#> 9918 total genes

cut -f1 SM_V10.featurecounts.out | cut -f2 -d ":" | sort | uniq | grep "Smp" | sort | uniq | wc -l
#> 10958 transcripts

#> 4326 5' SL genes /  9918 total genes
#> = 0.436 or 43.6%

#> 2571 extra / 9918 total genes
#> 0.259 or 25.9%

```