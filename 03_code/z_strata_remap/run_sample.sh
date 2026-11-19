#!/usr/bin/env bash
# Download a sample's runs, downsample, map to its own species' assembly and measure depth.
#
# usage: run_sample.sh WORKDIR SAMPLE_ID [THREADS] [TARGET_COVERAGE]
#
# Reads the row for SAMPLE_ID from samples.tsv (next to this script). Needs
# prepare_assembly.sh to have been run for that sample's assembly. A sample can
# have several runs (e.g. two lanes of one individual): their reads are pooled.
# Platform "bam" takes an existing BAM against the same assembly, given as the
# fastq field. Writes WORKDIR/samples/<SAMPLE_ID>/:
#   <SAMPLE_ID>.bam(.bai)                primary alignments, MAPQ >= 20
#   <SAMPLE_ID>.loci.regions.bed.gz      mean depth per CDS interval of each placed S. mansoni gene
#   <SAMPLE_ID>.w10k.regions.bed.gz      mean depth in 10 kb windows over the whole assembly
#
# MAPQ 20 keeps reads that could sit on either the PARs/ZSR or W out of the depth.
set -euo pipefail

WORKDIR=$1; SAMPLE=$2; THREADS=${3:-16}; TARGET=${4:-30}
MAPQ=20
SHEET=$(cd "$(dirname "$0")" && pwd)/samples.tsv

read -r SPECIES ASSEMBLY PLATFORM RUNS BASES FASTQ < <(
    awk -F'\t' -v s=${SAMPLE} '$1==s{print $3, $4, $8, $9, $10, $11}' ${SHEET})
[[ -n ${RUNS:-} ]] || { echo "no sample ${SAMPLE} in ${SHEET}" >&2; exit 1; }

ASM=${WORKDIR}/assemblies/${SPECIES}.${ASSEMBLY}
OUT=${WORKDIR}/samples/${SAMPLE}
mkdir -p ${OUT}
cd ${OUT}

if [[ ${PLATFORM} == "bam" ]]; then
    [[ -s ${SAMPLE}.bam ]] || { ln -s ${FASTQ} ${SAMPLE}.bam; samtools index ${SAMPLE}.bam; }
fi

if [[ ! -s ${SAMPLE}.bam.bai ]]; then
    # keep about TARGET x of the assembly; exome reads are all kept.
    # seqtk with the same seed keeps mates paired
    GSIZE=$(awk '{s+=$2} END{print s}' ${ASM}/genome.fa.fai)
    FRAC=$(awk -v t=${TARGET} -v g=${GSIZE} -v b=${BASES} -v p=${PLATFORM} \
        'BEGIN{f=(p=="illumina_exome") ? 1 : t*g/b; print (f>1 ? 1 : f)}')
    echo "${SAMPLE}: ${RUNS} ${PLATFORM}, ${BASES} bp over ${GSIZE} bp assembly, keeping fraction ${FRAC}"

    # seqtk reads a fraction of exactly 1 as "keep 1 read", so skip it then
    subsample() {
        if awk -v f=${FRAC} 'BEGIN{exit !(f >= 1)}'; then zcat -f; else seqtk sample -s 11 - ${FRAC}; fi
    }

    # pool all runs: *_1 into reads_1, *_2 into reads_2, unpaired (PacBio) into reads_1.
    # ENA lists some SRA-mirrored runs it won't serve (403); fetch those from NCBI instead
    rm -f reads_1.fq.gz reads_2.fq.gz part.fq.gz
    IFS=';' read -r -a URLS <<< "${FASTQ}"
    for url in "${URLS[@]}"; do
        file=$(basename ${url})
        case ${file} in
            *_2.fastq.gz) MATE=2 ;;
            *) MATE=1 ;;
        esac
        if ! wget -q -O part.fq.gz "https://${url}"; then
            run=${file%%[_.]*}
            echo "${SAMPLE}: ENA download failed for ${file}, using fasterq-dump ${run}" >&2
            [[ -s sra/${file%.gz} ]] || fasterq-dump --split-files -e ${THREADS} -O sra ${run}
            gzip -1 -c sra/${file%.gz} > part.fq.gz
        fi
        subsample < part.fq.gz | gzip -1 >> reads_${MATE}.fq.gz
        rm -f part.fq.gz
    done
    rm -rf sra

    case ${PLATFORM} in
        illumina|illumina_exome) PRESET=sr ;;
        pacbio_clr)  PRESET=map-pb ;;
        pacbio_hifi) PRESET=map-hifi ;;
        ont)         PRESET=map-ont ;;
        *) echo "unknown platform ${PLATFORM}" >&2; exit 1 ;;
    esac

    minimap2 -t ${THREADS} -ax ${PRESET} -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" ${ASM}/genome.fa reads_*.fq.gz |\
        samtools view -u -q ${MAPQ} -F 0x904 - |\
        samtools sort -@ 4 -m 2G -T ${SAMPLE}.sort -o ${SAMPLE}.bam -
    rm -f reads_*.fq.gz

    # stop before indexing if almost nothing mapped, so a rerun doesn't skip this sample
    MAPPED=$(samtools view -c ${SAMPLE}.bam)
    echo "${SAMPLE}: ${MAPPED} reads mapped at MAPQ >= ${MAPQ}"
    (( MAPPED >= 100000 )) || { echo "${SAMPLE}: too few reads mapped" >&2; exit 1; }
    samtools index ${SAMPLE}.bam
fi

mosdepth -t 4 -n --mapq ${MAPQ} --by ${ASM}/loci.bed ${SAMPLE}.loci ${SAMPLE}.bam
mosdepth -t 4 -n --mapq ${MAPQ} --by 10000 ${SAMPLE}.w10k ${SAMPLE}.bam
