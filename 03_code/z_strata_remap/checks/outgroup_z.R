# Outgroup depth per S. mansoni chromosome, and along each outgroup's own Z.
#
# usage: Rscript outgroup_z.R WORKDIR
#
# Needs strata/gene_depth.tsv from segment_strata.R. For each tier-2 outgroup
# sample: median normalised depth per S. mansoni chromosome (autosomes should be
# ~1), then the outgroup CHR_Z in 5 Mb bins with depth and the share of genes
# that are on the S. mansoni Z (sm_Z). A PAR would show as bins at ~1 in a female.

suppressMessages(library(tidyverse))
setwd(commandArgs(trailingOnly = TRUE)[1])

gd <- read_tsv("strata/gene_depth.tsv", show_col_types = FALSE)
dirs <- c(hame_1 = "heterobilharzia_americana", tszi_1 = "trichobilharzia_szidati", treg_1 = "trichobilharzia_regenti")
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
}
