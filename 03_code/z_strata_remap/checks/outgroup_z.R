# Outgroup depth per S. mansoni chromosome, and along each outgroup's own Z.
#
# usage: Rscript outgroup_z.R WORKDIR
#
# Needs strata/gene_depth.tsv from segment_strata.R. For each tier-2 outgroup
# sample: median normalised depth per S. mansoni chromosome (autosomes should be
# ~1), then the outgroup CHR_Z in 5 Mb bins with depth and the share of genes
# that are on the S. mansoni Z (sm_Z). A PAR would show as bins at ~1 in a female.
# Where a species has two samples, also their gene-by-gene depth ratio along its
# CHR_Z (first sample / second): ~0.5 where one is female and the other male.

suppressMessages(library(tidyverse))
setwd(commandArgs(trailingOnly = TRUE)[1])

gd <- read_tsv("strata/gene_depth.tsv", show_col_types = FALSE)
dirs <- c(hame_1 = "heterobilharzia_americana", tszi_1 = "trichobilharzia_szidati", tszi_2 = "trichobilharzia_szidati",
          treg_1 = "trichobilharzia_regenti", treg_2 = "trichobilharzia_regenti")
dirs <- dirs[names(dirs) %in% gd$sample_id]
own_z <- list()
for (s in names(dirs)) {
    asm <- Sys.glob(paste0("assemblies/", dirs[[s]], ".*"))
    pos <- read_tsv(file.path(asm, "loci.bed"), col_names = c("contig", "s", "e", "gene"), col_types = "cddc") |>
        group_by(gene, contig) |> summarise(og = min(s) / 1e6, .groups = "drop")
    d <- gd |> filter(sample_id == s) |> inner_join(pos, by = "gene")
    cat("\n", s, ": median depth per S. mansoni chromosome\n", sep = "")
    d |> group_by(chrom) |> summarise(n = n(), norm = round(median(norm), 2)) |> print()
    cat(s, ": along its own CHR_Z in 5 Mb bins (sm_Z = share of genes on the S. mansoni Z)\n", sep = "")
    d |> filter(contig == "CHR_Z") |> group_by(bin_Mb = 5 * floor(og / 5)) |>
        summarise(n = n(), norm = round(median(norm), 2), sm_Z = round(mean(chrom == "SM_V10_Z"), 2)) |> print(n = Inf)
    own_z[[s]] <- d |> select(gene, contig, og, chrom, norm) |> mutate(sample_id = s, species = dirs[[s]])
}

pairs <- bind_rows(own_z) |> group_by(species) |> filter(n_distinct(sample_id) == 2) |> ungroup()
for (sp in unique(pairs$species)) {
    ids <- sort(unique(pairs$sample_id[pairs$species == sp]))
    r <- pairs |> filter(species == sp) |> select(gene, contig, og, chrom, sample_id, norm) |>
        pivot_wider(names_from = sample_id, values_from = norm) |>
        filter(.data[[ids[1]]] > 0.05, .data[[ids[2]]] > 0.05) |>
        mutate(ratio = .data[[ids[1]]] / .data[[ids[2]]])
    cat("\n", sp, ": ", ids[1], "/", ids[2], " per S. mansoni chromosome\n", sep = "")
    r |> group_by(chrom) |> summarise(n = n(), ratio = round(median(ratio), 2)) |> print()
    cat(sp, ": ", ids[1], "/", ids[2], " along its own CHR_Z in 5 Mb bins\n", sep = "")
    r |> filter(contig == "CHR_Z") |> group_by(bin_Mb = 5 * floor(og / 5)) |>
        summarise(n = n(), ratio = round(median(ratio), 2), sm_Z = round(mean(chrom == "SM_V10_Z"), 2)) |>
        print(n = Inf)
}
