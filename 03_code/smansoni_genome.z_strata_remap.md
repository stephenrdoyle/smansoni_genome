# Schistosoma mansoni V10 genome: Z chromosome strata, remapped with defined boundaries

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk

Status: run on tier 1 (4 sexed species, 18 samples) and three tier-2 outgroups, 2026-11-19 to 21. Results in
section 4.
Replaces the coverage approach in `smansoni_genome.z_strata.md` (Fig. 3a). Merges the earlier
`smansoni_genome.z_strata_v2_plan.md` (sexed data) with the one-run-per-species remap.


## Why
- Fig. 3a strata come from read depth over S. mansoni Z orthologs, one or two libraries per species, normalised per
  species. Boundaries were drawn by eye (African 10,680,112-23,250,000; ancestral 23,250,001-43,743,322).
- Only S. japonicum (SRR6841388 female, SRR6841389 male) had a recorded sex. S. haematobium (ERR036251, ERR037800;
  sample ERS023537) and S. rodhaini (ERR310938; "ERS09443" in the Methods is ERS094843 / ERR114786) have no sex in ENA.
- Z-W gametologue dS on V10 (protein-coding W genes, 2026-11-19) does not fit the current two strata:
    - dS ~0.14 (median, n=21) for Z genes at 11.3-18.4 Mb vs ~0.61 (n=4) at 19.1-22.8 Mb; best single breakpoint
      between 18.42 and 19.14 Mb (MWU p=0.0002)
    - U2AF65 (Smp_019690 / Smp_348790, 29.5 Mb) saturated (ML dS ~2.6-2.8) even in the conserved C-terminus
    - the DnaJ pair at 23.46 Mb (Smp_336770 / Smp_020920) has dS ~0.27, i.e. not ancestral-like, despite falling
      inside the "ancestral" stratum
- aims:
    1. define ZSR / strata boundaries from female:male depth, with confidence intervals, per species
    2. test whether the 18.4-19.1 Mb dS step is a stratum boundary shared with other species
    3. date the strata relative to the African-Asian split using S. japonicum's assembled W


## Design
- each sample mapped to **its own species' assembly** (WBPS19, or NCBI for S. japonicum GCA_025215515.1, which has a W)
- S. mansoni v10 proteins placed on each assembly with miniprot, so every species is read gene by gene on v10 Z
  coordinates; minimap2 asm20 of each assembly against v10 Z is kept to find the species' Z and check rearrangements
- MAPQ >= 20 so reads that could sit on PAR/ZSR or W don't count
- profiles segmented:
    - **tier 1, sexed (S. mansoni, S. japonicum, S. rodhaini, S. haematobium):** female:male gene-depth ratio per species, replicates averaged; S. mansoni also in 10 kb windows
      (no liftover needed). Exome capture (S. rodhaini PRJNA773498) is used only as F/M, where capture bias cancels
    - **tier 2, breadth:** one run per species, depth / median over genes from S. mansoni autosomes; hemizygous ≈ 0.5.
      Sex inferred from the data (most runs have none recorded); a mixed pool gives a shallower dip, still usable
- boundaries by circular binary segmentation (DNAcopy), 95% intervals from a block bootstrap; each boundary flagged if
  its interval overlaps the 18.42-19.14 Mb dS step, and every profile summarised over the old strata with the African
  stratum split at the dS step


## Data (`z_strata_remap/samples.tsv`)

| species | assembly | tier 1: female | tier 1: male | tier 2 (sex inferred) |
|---|---|---|---|---|
| S. mansoni | V10 (WBPS19) | ERR3467675 (PRJEB2320); 6520_5 PCR-free female BAM optional | ERR3467676 (PRJEB2320) | |
| S. japonicum | GCA_025215515.1 (Xu et al. 2023): Z = CM046069.1, W = CM046068.1 | SRR17422882, -84, -86 (PRJNA793934) + SRR6841388 | SRR17422883, -85, -87 + SRR6841389 | |
| S. rodhaini | TD1_PRJEB44434 (GCA_944470415.2), no W | ERR7978144 (WGS); exome SRR16526441-44 (2 individuals × 2 lanes) | exome SRR16526437-40 | ERR310938, ERR114786 |
| S. haematobium | TD1_PRJEB44434 | ERR037800 + ERR036251 (SAMEA980563, ShaemFemaleNoPCR)* | SRR15400745 (Senegal), SRR15400771 (Mali) (PRJNA78265)*; backups SRR15400770, SRR15400772 | |
| haematobium group | PRJEB44434 / TD1 | | | bovis ERR119622, curassoni ERR119623, guineensis ERR119612, intercalatum ERR119613, mattheei ERR103051, margrebowiei ERR310940 |
| S. spindale, S. turkestanicum | PRJEB44434 | | | PacBio ERR10166842, ERR10166850 (turkestanicum's only Illumina run is male) |
| Heterobilharzia americana | TD1_PRJEB44434 | | | PacBio ERR10166852 |
| T. regenti / T. szidati | PRJEB44434 | | | single-miracidium cercarial clones ERR266615 / ERR266616 |
| E. caproni | PRJEB1207 | | | ERR065032 (hermaphrodite control) |

- tier-2 species have a backup run each, for when the primary turns out male
- \* S. haematobium sexes come from sample names only (ENA sex field empty): the female from the alias
  "ShaemFemaleNoPCR-sc-2011-02-16-1080706" (the female used in the submitted Fig. 3a), the males from the titles
  "Schistosoma haematobium DNA single male <country>". `sex_calls.tsv` checks both from the coverage
- sex in older Sanger samples is in the "ArrayExpress-Sex" attribute, not the ENA sex field
- none of the PRJEB44434 TD assemblies has a recorded sex

Scripts in `03_code/z_strata_remap/`:

| file | what it does |
|---|---|
| `samples.tsv` | samples, tier, sex, platform, runs (lanes pooled), assembly, ENA fastq URLs |
| `env.yaml` | conda environment (minimap2, miniprot, samtools, mosdepth, seqtk, seqkit, R + DNAcopy) |
| `prepare_ref.sh` | v10 genome, Z, annotation, longest-isoform proteins, gene coordinates |
| `prepare_assembly.sh` | download an assembly (WBPS or NCBI); miniprot v10 proteins onto it; asm20 to v10 Z |
| `run_sample.sh` | stream + downsample reads to ~30x (exome kept whole), map (minimap2 sr / map-pb), mosdepth over loci and 10 kb windows |
| `segment_strata.R` | normalise, F/M ratios, infer sex, segment, bootstrap boundaries, interval summary, plot |
| `checks/fm_check.R` | S. mansoni F:M in 10 kb windows per region (validation) |
| `checks/wsr_check.R` | S. mansoni W depth, windows vs W gene loci |
| `checks/sjap_z_order.R` | S. japonicum F:M read along its own Z |
| `checks/inversion_check.R` | orientation and junctions of the 20.9-23.25 Mb segment in each assembly (gene order only) |
| `checks/outgroup_z.R` | outgroup depth per S. mansoni chromosome and along the outgroup's own Z |


## Setup
```bash
WORKDIR=/nfs/users/nfs_s/sd21/lustre_link/schistosoma_mansoni/project-genome-v10/04_analysis/z_strata_remap
CODE=/nfs/users/nfs_s/sd21/lustre_link/schistosoma_mansoni/project-genome-v10/smansoni_genome/03_code/z_strata_remap

mkdir -p ${WORKDIR} && cd ${WORKDIR}

mamba env create -f ${CODE}/env.yaml
conda activate z_strata_remap

# every new login needs the env and both paths; jobs inherit them from the submitting shell
cat > ${WORKDIR}/env.sh <<EOF
conda activate z_strata_remap
WORKDIR=${WORKDIR}
CODE=${CODE}
EOF

bsub.py 8 ref_sman "bash ${CODE}/prepare_ref.sh ${WORKDIR}"
```

## 1. Place S. mansoni genes on each assembly
```bash
# one job per assembly (16); miniprot on a ~1 Gb Trichobilharzia genome takes a few hours
awk -F'\t' 'NR>1{print $3, $4}' ${CODE}/samples.tsv | sort -u | while read SPECIES ASSEMBLY; do
     bsub.py --threads 16 40 asm_${SPECIES} "bash ${CODE}/prepare_assembly.sh ${WORKDIR} ${SPECIES} ${ASSEMBLY} 16"
done

# how many v10 genes were placed, and how many are usable (unambiguous)
for i in assemblies/*/liftover.tsv; do
     echo -e "$(dirname ${i} | xargs basename)\t$(tail -n +2 ${i} | wc -l)\t$(awk '$5==0' ${i} | wc -l)";
done

# which contigs carry S. mansoni Z in each assembly
for i in assemblies/*/vs_SmZ.paf; do
     echo ${i}; awk '$12>=20{a[$1]+=$10} END{for(c in a) print c, a[c]}' ${i} | sort -k2,2nr | head -3;
done
```

## 2. Map
```bash
# primary samples (tier 1 and 2); run a backup only if segment_strata.R calls the primary male
awk -F'\t' 'NR>1 && $7=="primary"{print $1}' ${CODE}/samples.tsv | while read SAMPLE; do
     bsub.py --threads 16 40 map_${SAMPLE} "bash ${CODE}/run_sample.sh ${WORKDIR} ${SAMPLE} 16 30"
done

# optional: the paper's 6520_5 female. Set its fastq field in samples.tsv to the BAM path first
# (BAM must be against V10 with the WBPS sequence names)
bsub.py --threads 4 8 map_sman_6520_5 "bash ${CODE}/run_sample.sh ${WORKDIR} sman_6520_5 4"
```

## 3. Segment, define boundaries, plot
```bash
# third argument: Z-W gametologue dS table, drawn as a track above the coverage
# 04_analysis/zw_gametologue_ds.tsv: 35 protein-coding pairs + 5 W pseudogenes (smansoni_genome.gametologue_dS.md)
bsub.py --threads 2 16 strata "Rscript ${CODE}/segment_strata.R ${WORKDIR} 1000 ${CODE}/../../04_analysis/zw_gametologue_ds.tsv"
```

Outputs in `${WORKDIR}/strata/`:

- `sex_calls.tsv`: inferred sex per tier-2 (and tier-1) sample, next to the ENA metadata
- `segments.tsv`: every segment per profile, its extent on v10 Z and median value
- `breakpoints.tsv`: every boundary: position (midpoint of the gap between flanking genes/windows), the gap itself
  (resolution limit), 95% bootstrap interval, support, whether it is a full-to-half step (`hemizygous_step`), and
  whether its interval overlaps the 18.42-19.14 Mb dS step (`in_dS_gap`)
- `interval_summary.tsv`: median per profile over PAR1, African_a (10.68-18.42), African_b (19.14-23.25), ancestral,
  Asian, PAR2, with a Wilcoxon test of African_a vs African_b. With hundreds of autocorrelated points the p-value is
  descriptive only; read it with the segment medians.
- `z_strata_remap.pdf`: dS track, then one panel per profile on v10 Z (F/M profiles first within each species)

What to look for (aim 2): a species whose hemizygous region starts at ~19.1 Mb rather than ~10.7 Mb, i.e. a boundary
with `in_dS_gap` TRUE. In S. mansoni itself both halves of the ZSR are hemizygous, so the step can only show as a
boundary in another species.


## 4. Results (2026-11-19 to 21)

```bash
# all from ${WORKDIR}; checks/*.R take WORKDIR as their only argument
Rscript ${CODE}/checks/fm_check.R ${WORKDIR}
Rscript ${CODE}/checks/wsr_check.R ${WORKDIR}
Rscript ${CODE}/checks/sjap_z_order.R ${WORKDIR}     # needs strata/gene_depth.tsv
Rscript ${CODE}/checks/inversion_check.R ${WORKDIR}  # gene order only, no reads
Rscript ${CODE}/checks/outgroup_z.R ${WORKDIR}       # needs strata/gene_depth.tsv
```

### Mapping
- every tier-1 sample and the three primary outgroups mapped: S. mansoni 59-62 M reads at MAPQ >= 20, S. rodhaini
  WGS 53 M, exome 30-35 M per sample, treg_1 134 M, tszi_1 112 M, hame_1 2.1 M PacBio reads (~30x)
- download failures on the first submission (all at the same minute, likely ENA refusing many connections at once)
  and an exome subsampling bug (seqtk reads a fraction of 1 as "keep 1 read") are fixed in `run_sample.sh`

### S. mansoni validation (`fm_check.R`, `wsr_check.R`)
| region | windows | F | M | F/M |
|---|---|---|---|---|
| autosomes | 29,877 | 1 | 1 | 1 |
| Z PARs | 5,363 | 1.00 | 1.00 | 1.00 |
| ZSR African 10.7-18.4 | 774 | 0.52 | 1.01 | 0.51 |
| ZSR African 19.1-23.25 | 411 | 0.52 | 1.02 | 0.51 |
| ZSR ancestral | 2,049 | 0.50 | 0.99 | 0.51 |
| WSR | 597 | 0.03 | 0 | - |

- no W reads leak onto the Z at MAPQ 20 (ZSR F/M 0.51)
- WSR: window median 0.03 in the female because most of the W is repeat; the 160 W gene loci are at 0.50 in the
  female and 0 in the male

### Sex calls
- all 14 tier-1 WGS samples called as their metadata say
- outgroups: hame_1 (H. americana) female, tszi_1 (T. szidati) female, treg_1 (T. regenti) male

### Strata in four species (`strata/interval_summary.tsv`, `breakpoints.tsv`)
F/M (median per interval, S. mansoni v10 Z coordinates):

| interval | S. mansoni | S. rodhaini WGS / exome | S. haematobium | S. japonicum (4F 4M) |
|---|---|---|---|---|
| PAR1 (< 10.68) | 1.01 | 1.00 / 0.98 | 0.95 | 1.00 |
| African_a 10.68-18.42 | 0.52 | 0.52 / 0.53 | 0.55 | 0.99 |
| African_b 19.14-23.25 | 0.51 | 0.50 / 0.52 | 0.52 | 0.97 |
| ancestral 23.25-43.74 | 0.51 | 0.51 / 0.52 | 0.51 | 0.51 |
| Asian 43.74-51.3 | 1.00 | 0.98 / 0.97 | 1.03 | 0.53 |
| PAR2 (> 51.3) | 1.00 | 1.00 / 1.00 | 1.00 | 1.00 |

Full-to-half boundaries (bootstrap support 0.94-1):
- S. mansoni 10.6 (CI 10.4-10.7) and 43.7 Mb; S. rodhaini 10.7 and 43.7 Mb (WGS and exome agree)
- S. haematobium 10.1 (CI 9.6-10.2) and 43.9 Mb: left edge ~0.5 Mb beyond the others; PAR1 F/M slightly low (0.95)
- S. japonicum, all four females: 20.9-21.8, 23.2-44.5 and 46.3-51.2 Mb. Ancestral boundary between Smp_347690
  (23.2 Mb, F/M 1.02) and Smp_340330 (23.3 Mb, 0.46), i.e. the paper's 23.25 Mb, now from sexed data

What this shows:
- ancestral stratum hemizygous in all four species (before the African-Asian split); African stratum only in the
  African species; Asian stratum only in S. japonicum - the paper's model, now from F/M in each species
- **depth does not split the African stratum**: no boundary in any species falls in the 18.42-19.14 Mb dS step
  (`in_dS_gap` FALSE throughout), and African_a = African_b in every African species. The two-step history rests on
  dS alone; no sampled species diverged between the two steps
- U2AF65 (dS ~2.6) is inside the stratum shared by all four species, consistent with it being oldest
- the DnaJ pair Smp_336770 (23.46 Mb, dS 0.27) is hemizygous in S. japonicum too (F/M 0.63): a low-dS gametologue
  inside the ancestral stratum, still unexplained (gene conversion, or a later W copy)
- per-sample Wilcoxon p-values are not meaningful: shaem_F's p = 2e-5 (African_a 0.51 vs African_b 0.54) is matched
  by the same shift in its males, i.e. coverage bias, which the F/M ratio cancels (p = 0.08)

### S. japonicum: the 20.9-21.8 Mb block is an inversion, not a stratum (`sjap_z_order.R`)
- every gene in S. mansoni 19-24 Mb sits on CM046069.1 (the S. japonicum Z), hemizygous or not
- read in S. japonicum order, the hemizygous genes form one block, S. japonicum Z ~37.5-62.0 Mb (~25 Mb). S. mansoni
  positions along it are shuffled: S. japonicum 37-43 Mb holds S. mansoni 42-50 Mb (end of ancestral + Asian),
  44-60 Mb holds the ancestral stratum in blocks
- at the end of the block:

| S. mansoni genes | S. japonicum Z | F/M |
|---|---|---|
| 24 -> 23.3 (start of ancestral) | 59-61 Mb | 0.5 |
| 20.9 -> 21.8 (Smp_084270 ... Smp_166610) | 61.2-62.0 Mb | 0.5 |
| 21.8 -> 23.1 (Smp_104720 ...) | 62.3-63 Mb | 1.0 |
| 20.9 -> 20.5 (Smp_167870, Smp_092880) | 63.9-64.1 Mb | 1.0 |

- i.e. S. mansoni 20.9-23.2 Mb is inverted in S. japonicum. Its 20.9 end now abuts the ancestral stratum, and the
  S. japonicum hemizygous region runs ~0.9 Mb into it, ending between Smp_166610 and Smp_104720 (neighbours in both
  species). The "island" is an artefact of drawing S. japonicum on S. mansoni coordinates
- the high-dS genes are not this block: Smp_007630 (19.14), Smp_340400 (20.27), Smp_092880 (20.51) and
  Smp_158740 (22.85) are diploid in S. japonicum; only Smp_166600 (21.71) falls in the hemizygous edge
- the 44.5-46.3 Mb "diploid gap" is probably the same kind of artefact (its genes appear inside and just outside
  the block)

### Polarising the inversion (`inversion_check.R`)
Distance between segment ends on each assembly's Z (median of 3 end genes; A = 20.0 Mb to Smp_084270,
B = Smp_084270 to 23.25 Mb, C = 23.25-24.5 Mb):

| assembly | A_end-B_start | B_end-C_start | A_end-B_end | B_start-C_start | arrangement |
|---|---|---|---|---|---|
| S. haematobium (control) | 48 kb | 195 kb | 2.4 Mb | 2.6 Mb | S. mansoni |
| S. japonicum (control) | 2.7 Mb | 2.7 Mb | 141 kb | 159 kb | S. japonicum |
| H. americana | 131 kb | 505 kb | 4.3 Mb | 3.7 Mb | S. mansoni |
| T. regenti | 85 Mb | 39 Mb | 1.3 Mb | 45 Mb | shuffled |
| T. szidati | 92 Mb | 17 Mb | 1.0 Mb | 77 Mb | shuffled |

- columns 2-3 are S. mansoni-order junctions, 4-5 S. japonicum-order junctions
- orientation (Spearman rho per segment) agrees: S. haematobium all +1; S. japonicum flanks -1 / -0.86, B +1
- H. americana keeps both S. mansoni junctions, which points to S. mansoni order being ancestral and the inversion
  arising in the Asian lineage. If so, the ancestral boundary at 23.25 Mb is an inversion breakpoint only in
  S. japonicum, and the inversion did not stop recombination there in S. mansoni
- in both Trichobilharzia the region is scattered across their Z (rho -0.3 to 0.6); the one ~1 Mb S. japonicum-type
  junction is likely chance or breakpoint reuse
- H. americana and Trichobilharzia are equally distant from Schistosoma (sister clades, see below), so this rests on
  H. americana alone, the one outgroup whose junctions survive. In the paper: "H. americana shares the S. mansoni
  arrangement at both breakpoints, suggesting the inversion arose in the Asian lineage"
- the minimap2 asm20 alignments (`vs_SmZ.paf`) are useless at this distance (10-40 kb aligned); gene placement is
  what counts

### Outgroups: the whole Z is hemizygous (`outgroup_z.R`)
- the Z is the same chromosome across the family: 90-100% of genes on each outgroup's CHR_Z are on the S. mansoni Z
  (~150-180 Mb vs 86.7 Mb), and ~99/101 (T. regenti), 83/84 (T. szidati), 104/108 (H. americana) genes from
  S. mansoni 20-24.5 Mb are on it
- normalisation is fine: autosomes 0.93-1.03 in both females; treg_1 (male) ~1 everywhere

| sample | S. mansoni autosomes | SM_V10_Z | SM_V10_WSR | own CHR_Z, 5 Mb bins |
|---|---|---|---|---|
| hame_1 (female) | 0.95-1.03 | 0.53 | 0.66 | 0.42-0.63 along all 150 Mb |
| tszi_1 (female) | 0.93-1.02 | 0.53 | 0.51 | 0.48-0.55 from 20 Mb on; 0.78-0.80 in some bins at 0-20 Mb (6-11 genes each) |
| treg_1 (male) | 0.99-1.01 | 0.99 | 0.95 | 0.95-1.02 |

- in H. americana and T. szidati females every S. mansoni interval, PAR1 and PAR2 included, is at 0.51-0.56: no
  PAR at 5 Mb resolution
- T. szidati CHR_Z 0-20 Mb, gene by gene (2026-11-21): 34 genes at 2.3-19.6 Mb, median 0.75, range 0.31-2.55;
  back to ~0.48 from 21.5 Mb. Not a PAR (that would be a tight run at ~1). Either the youngest part of its ZW
  (W copies not yet diverged, so W reads map onto Z) or an assembly artefact at the chromosome end (this assembly
  has the most ambiguous placements). 10 kb windows can't help: repeats push CHR_Z windows to 0.5-2.0 even where
  genes are at 0.5. Unresolved without a male; backup tszi_2 (ERR119617) has no recorded sex
- phylogeny (UCE phylogenomics, Pathogens 2022, 11:769): Heterobilharzia + Schistosomatium (HS) is sister to the
  derived avian schistosomes (DAS, incl. Trichobilharzia), strongly supported; Schistosoma is sister to
  Bivitellobilharzia (SB), and (SB + Macrobilharzia) is sister to (HS + DAS); Austrobilharzia + Ornithobilharzia
  is sister to all. So both outgroups here belong to one clade, equally distant from Schistosoma
- the genes of the Schistosoma PARs were Z-linked in the common ancestor and are hemizygous in both outgroups. A W
  can't regain lost genes, so they were still PAR in the ancestor of Schistosoma and the outgroups, and
  recombination stopped along the whole Z later on the HS + DAS side. One event on the HS + DAS stem is enough;
  separate events in each lineage can't be excluded. Schistosoma is the lineage that kept large PARs; its strata
  boundaries are its own
- depth can't date the ancestral stratum against the outgroups (everything is hemizygous there). That would need
  dS on any outgroup Z-W gametologues that survive, or a female from the SB side (Bivitellobilharzia) or from
  Austrobilharzia / Ornithobilharzia

### Is the outgroup W assembled? (2026-11-21)
- window depth misleads for W: W windows are mostly repeat and sit near 0 even in a female (S. mansoni female WSR
  window median 0.03), and duplicated contigs drop to 0 through the MAPQ filter. Gene-level depth and
  self-alignment are needed
- H. americana: no W assembled (CHR_Z at 0.5, almost no unplaced sequence)
- T. szidati (female): T. szidati CHR_7 looked half-depth in windows but its genes sit at 0.93 (homologue of
  S. mansoni chr 7; repeat artefact, not a neo-Z). Its zero-depth unplaced contigs are mostly >= 99% copies of
  other sequence (MAPQ artefact). No W identifiable from this sample
- T. regenti (treg_1, male): 67 unplaced contigs >= 20 kb at zero depth, 64 with no close copy elsewhere (minimap2
  asm20). Genes on them are enriched for second copies of Z genes. Five contigs (UNPLACED_00115, 00117, 00373,
  00423, 00478; 21 genes) paired with CHR_Z genes and given dS (`03_code/zw_gametologue_dS/treg_zw.py`,
  `04_analysis/treg_zw_gametologue_ds.tsv`, T. regenti WBPS19 annotation):
    - most copies near-identical to their Z partner (dS 0-0.06); one gene in two pieces diverged (TREG1_137610/620
      vs Z TREG1_1250, dS ~0.13, aa 92%); five hits not Z copies (excluded)
    - depth of each partner in the male decides whether the W copy's zero is real:

| W-contig gene | Z partner | Z depth (male) | W depth (male) | reading |
|---|---|---|---|---|
| TREG1_137610/620 (00117), dS 0.13 | TREG1_1250 | 1.04 | 0.00 | W copy truly absent from the male: best W evidence |
| TREG1_137570/580 (00115), dS 0.01-0.03 | TREG1_25140 / 25130 | 0.50 / 0.80 | 0.00 / 0.02 | mostly real absence; young copies |
| TREG1_142340 (00423) | TREG1_28520 | 0.30 | 0.01 | partly masked |
| TREG1_143810 (00478) | TREG1_1910 | 0.24 | 0.01 | mostly masked |
| TREG1_141010 (00373) | TREG1_13090 | 0.07 | 0.00 | identical copies; uninformative |

- so UNPLACED_00117 (and probably 00115) is sequence missing from the male that carries a diverged copy of a Z gene:
  the T. regenti W, or sequence specific to the reference individual. The four-point test on protein distances
  puts each W copy of the Smp_120320 homologue with its own species' Z copy (S. mansoni pair Smp_120320 /
  Smp_310950, dS 0.17): consistent with independent suppression, one gene only
- correction: UNPLACED_00373 is not evidence of W. Its gene is 99.8% identical to its Z partner and both have
  zero depth, so the S. mansoni W protein (Smp_320940) landed there only because the two copies are identical
- the near-identical copies can't date the outgroup's whole-Z suppression. A T. regenti female (treg_2, sex
  unknown) would test the diverged parts of 00117 / 00115: ~0.5 if W

### Against published work (checked 2026-11-21)
- Picard et al. 2018, eLife 7:e35684 (PMC6089595): Z-specific genes from F/M coverage in S. mansoni and
  S. japonicum (pooled sexed worms) and S. haematobium (females only). Defined S0 (shared, 367 genes), S1mans
  (African, 180) and S1jap (Asian, 66). Also Z dosage: partial compensation by Z upregulation in both sexes
- Elkrewi et al. 2021, MBE 38:5345 (PMC8662593): W-linked genes in S. mansoni and S. japonicum. W copies of u2af2
  (= U2AF65) from both species cluster, so S0 stopped recombining before the African-Asian split. Median Z-W dS:
  S0 0.88 (Sm) / 0.61 (Sj), S1man 0.18, S1jap 0.47. A very young S. japonicum stratum, S2jap (28.8 Mb, 720 genes,
  dS 0.085), at the distal end of its Z, found from female-male FST, not coverage. u2af2 proposed as the
  sex-determination candidate
- Xu et al. 2023 (online 2022), Mol Ecol Resour 23(1) (doi 10.1111/1755-0998.13689): the S. japonicum V3 assembly
  used here (chrZ = CM046069.1). Strata from F/M depth, female SNP density, W-Z identity and orthology to
  S. mansoni V7 (fused chrZW). On their chrZ (Fig. 4, approx.): PAR 0-11 Mb, S4 ~11-33, S2 ~33, S3 ~34-35,
  S1 ~35-40, S0 ~40-62, PAR ~62-87. Only S1 and S0 (and S2, small) show 2x M/F depth; S3 and S4 are young,
  equal depth, female-biased SNPs. S0 = the part of their SDR shared with the S. mansoni SDR (20.5% of it).
  Gene trees for 7 genes (S4, S2) support strata forming independently after the split. No W fragments
  assembled for S0; u2af2 found only on chrZ. They note "Z-linked rearrangements" between the species and a
  Z-linked inversion that "probably reversed the order of SjapS2 and SjapS3", but don't describe our inversion
- agreement with Xu: our S. japonicum hemizygous block (its Z 37.5-62.0 Mb) is their S1 + S0, with the same right
  edge (62.0-62.3 Mb). Our "Asian stratum" (S. mansoni 43.7-51.3) is their SjapS1. Their S2-S4 are invisible to
  depth, as expected
- difference: by their orthology rule, the inverted 0.9 Mb piece (S. mansoni 20.9-21.8, on their chrZ 61.2-62.0)
  falls in S0, because it is hemizygous in both species. We read it as an Asian-lineage extension: in S. mansoni it
  sits in the African stratum, diploid in S. japonicum on both sides of it, its Z-W dS (Smp_166600, 0.53) is below
  S0 (0.88, Elkrewi), and the inversion put it next to S0. Worth stating in the paper
- Grossman, Short & Cain 1981, Chromosoma 84:413-430 (doi 10.1007/BF00286030): C-banded karyotypes. African
  Schistosoma: Z and W similar in size, W partly heterochromatic. Asian (S. japonicum, S. mekongi): W smaller than
  Z (deletion). Schistosomatium douthitti and H. americana: grossly heteromorphic ZW, Z the largest chromosome,
  W about a third of its length (S. douthitti relative lengths 33.2 vs 12.3); H. americana W also heterochromatic.
  They infer translocation of W material onto the Z. Trichobilharzia not studied
- no genome-level work found on Z coverage in Heterobilharzia or Trichobilharzia

| our result | status |
|---|---|
| three strata (ancestral, African, Asian) | known (Picard 2018); we add sexed F/M per species, S. haematobium males, S. rodhaini, boundaries with CIs |
| ancestral stratum predates the African-Asian split | known (Elkrewi 2021, u2af2 gene trees); our shared hemizygosity agrees |
| U2AF65 the most diverged pair | known (u2af2 in S0); fits their sex-determination candidate |
| African stratum in two steps (dS 0.15 vs 0.63 at 18.4-19.1 Mb) | not reported; our lower block matches their S1man dS (0.18). Rests on 5 pairs |
| S. japonicum hemizygous block and its edges | agrees with Xu 2023 (their S1 + S0, right edge ~62 Mb) |
| 20.9-23.25 Mb inversion, Asian lineage | new: Xu's Fig. 4a shows the Z rearranged but neither describes nor polarises it, and counts the inverted 0.9 Mb as S0 |
| whole Z hemizygous in H. americana and T. szidati | new at genome level; fits the heteromorphic ZW of H. americana in Grossman et al. 1981 |

- coverage only sees strata where the W has diverged enough to lose reads: Xu's S3 and S4 (and Elkrewi's S2jap)
  are invisible here. Say "strata visible in coverage" for S. japonicum
- in the paper: cite Picard 2018 and Elkrewi 2021 for the three-strata model and S0's age; the new claims are the
  sexed per-species boundaries, the African two-step dS, the inversion and the outgroups

### Current model
- ancestral stratum (23.25-43.74 Mb): hemizygous in all four Schistosoma; contains U2AF65 (dS ~2.6)
- African stratum (10.68-23.25 Mb): after the African-Asian split; dS suggests two steps (~18.4-19.1 Mb), which
  depth can't test
- Asian stratum: S. japonicum only, one block on its own Z, including ~0.9 Mb past the ancestral boundary brought
  next to it by the 20.9-23.25 Mb inversion (Asian lineage)
- outgroups (H. americana, T. szidati; one clade, HS + DAS): the whole Z hemizygous, Schistosoma PAR genes
  included - likely one event on that clade's stem


## 5. Date the strata with S. japonicum's W (to do; lower priority)
- Elkrewi et al. 2021 already did this for u2af2 (S0 predates the split). What's left: other S0 genes, and the
  DnaJ pair
- find Z-W gametologue pairs in S. japonicum (W CM046068.1) and compute dS as for S. mansoni (zw_ds.py)
- gene trees for gametologues present in both species: if Sm-W and Sj-W cluster together, recombination stopped
  before the African-Asian split (ancestral); if Sm-W clusters with Sm-Z, it stopped after (lineage-specific)
- U2AF65 and the 23.46 Mb DnaJ pair are the priority genes


## Notes on the method
- Gene depth is the length-weighted mean depth over the CDS of each placed gene. Genes with < 300 bp CDS, normalised
  depth < 0.05 or > 3, or an ambiguous placement (a second hit almost as good) are dropped. F/M ratios also need raw
  depth >= 5 in every replicate and male normalised depth > 0.2 (drops uncaptured genes in the exome).
- Segments shorter than 10 points, or differing from a neighbour by < 0.2 log2, are merged after CBS. A female's
  full-to-half step is 1 log2, so this only removes noise.
- Bootstrap: fitted segment medians plus residuals resampled in blocks (10 genes, or 20 windows), re-segmented;
  each original boundary is matched to the nearest bootstrap boundary within 50 genes (100 windows).
- Other species are projected onto v10 Z gene order, so an inversion relative to S. mansoni shows as extra boundaries
  rather than a moved one. Check `vs_SmZ.paf` / `liftover.tsv` at any boundary before reading it as a stratum
  (this is also the inversion test Hirai et al. proposed and v5 couldn't support).
- Tested on simulated coverage: females, males, a mixed pool, replicated F/M, exome F/M with capture bias, and a
  species whose hemizygous region starts at 19.1 Mb. Every boundary was recovered within its interval, males gave
  no boundaries, and the 19.1 Mb boundary was flagged `in_dS_gap`.


## Output for the paper
- new Fig. 3a: dS track plus per-species female:male (or normalised) depth along S. mansoni Z with segmented
  boundaries (± CI)
- table of boundaries per species in S. mansoni coordinates
- Methods: sample sexes and accessions (fix "ERS09443" -> ERS094843)

## Open
- T. szidati CHR_Z 0-20 Mb: intermediate depth (section 4); map tszi_2 and, if male, compare F/M
- script the model figure (strata bar, dS track, hemizygous blocks per species, inversion) from `strata/` outputs
- results tables in section 4 were transcribed from the farm output; copy `strata/*.tsv` into `04_analysis/` if
  they go in the paper
- S. mansoni clonal male PCR-free library ID (Table S1), to pair with 6520_5 instead of ERR3467675/76
- S. bovis mated pair (PRJNA451066, mixed) vs male ERR3061640: weaker ~0.75 signal; not in samples.tsv yet
- request female data for the haematobium group / new single-sex sequencing for Heterobilharzia, S. spindale,
  S. turkestanicum?
