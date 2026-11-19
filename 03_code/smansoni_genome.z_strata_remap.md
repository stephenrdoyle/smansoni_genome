# Schistosoma mansoni V10 genome: Z chromosome strata, remapped with defined boundaries

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk

Status: scripts written and tested on simulated coverage; nothing run on real data yet.
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


## Setup
```bash
WORKDIR=/nfs/users/nfs_s/sd21/lustre_link/schistosoma_mansoni/V10/STRATA_REMAP
CODE=/nfs/users/nfs_s/sd21/lustre_link/schistosoma_mansoni/V10/smansoni_genome/03_code/z_strata_remap

mkdir -p ${WORKDIR} && cd ${WORKDIR}

mamba env create -f ${CODE}/env.yaml
conda activate z_strata_remap

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


## 4. Date the strata with S. japonicum's W (to do)
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
- S. mansoni clonal male PCR-free library ID (Table S1), to pair with 6520_5 instead of ERR3467675/76
- S. bovis mated pair (PRJNA451066, mixed) vs male ERR3061640: weaker ~0.75 signal; not in samples.tsv yet
- request female data for the haematobium group / new single-sex sequencing for Heterobilharzia, S. spindale,
  S. turkestanicum?
