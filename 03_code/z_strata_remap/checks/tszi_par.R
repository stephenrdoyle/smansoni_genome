# Where the T. szidati PAR genes sit on the S. mansoni Z.
#
# usage: Rscript tszi_par.R WORKDIR
#
# Needs strata/gene_depth.tsv with tszi_1 (female) and tszi_2 (mostly male pool).
# For each gene on T. szidati CHR_Z: depth in both samples, their ratio (~1 in a
# PAR, ~0.6 where hemizygous for this pool), and its S. mansoni position and
# interval. Lists the genes in the first 25 Mb, then tabulates T. szidati region
# (PAR, < 20 Mb, vs the rest) against S. mansoni interval.

suppressMessages(library(tidyverse))
setwd(commandArgs(trailingOnly = TRUE)[1])

sm_interval <- function(chrom, pos) case_when(
    chrom != "SM_V10_Z" ~ chrom,
    pos < 10680112 ~ "PAR1",
    pos < 23250000 ~ "African",
    pos < 43743322 ~ "ancestral",
    pos < 51300000 ~ "Asian (Sj only)",
    TRUE ~ "PAR2")

gd <- read_tsv("strata/gene_depth.tsv", show_col_types = FALSE)
asm <- Sys.glob("assemblies/trichobilharzia_szidati.*")
pos <- read_tsv(file.path(asm, "loci.bed"), col_names = c("contig", "s", "e", "gene"), col_types = "cddc") |>
    group_by(gene, contig) |> summarise(og = min(s) / 1e6, .groups = "drop")

z <- gd |>
    filter(sample_id %in% c("tszi_1", "tszi_2")) |>
    select(gene, chrom, start, sample_id, norm) |>
    pivot_wider(names_from = sample_id, values_from = norm) |>
    inner_join(pos, by = "gene") |>
    filter(contig == "CHR_Z", tszi_1 > 0.05, tszi_2 > 0.05) |>
    mutate(ratio = tszi_1 / tszi_2, sm_Mb = round(start / 1e6, 2), sm_interval = sm_interval(chrom, start),
           region = if_else(og < 20, "Tsz PAR (< 20 Mb)", "Tsz rest")) |>
    arrange(og)

cat("T. szidati CHR_Z genes in the first 25 Mb\n")
z |> filter(og < 25) |>
    transmute(gene, og_Mb = round(og, 2), tszi_1 = round(tszi_1, 2), tszi_2 = round(tszi_2, 2),
              ratio = round(ratio, 2), chrom, sm_Mb, sm_interval) |>
    print(n = Inf)

cat("\nT. szidati region x S. mansoni interval: genes and median ratio\n")
z |> group_by(region, sm_interval) |> summarise(n = n(), ratio = round(median(ratio), 2), .groups = "drop") |>
    print(n = Inf)

cat("\nS. mansoni Z genes by interval: where they sit on T. szidati CHR_Z\n")
z |> filter(chrom == "SM_V10_Z") |> group_by(sm_interval) |>
    summarise(n = n(), in_Tsz_PAR = sum(og < 20), og_min = round(min(og), 1), og_median = round(median(og), 1),
              ratio = round(median(ratio), 2)) |>
    print()
