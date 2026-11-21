# S. japonicum female:male depth read along its own Z (CM046069.1), not S. mansoni's.
#
# usage: Rscript sjap_z_order.R WORKDIR
#
# Needs strata/gene_depth.tsv from segment_strata.R. Genes are placed on the
# S. japonicum Z from loci.bed; F/M is the mean of the four females over the
# mean of the four males; hemizygous = F/M < 0.7. Prints:
#   1. runs of hemizygous / diploid genes in S. japonicum order (runs of >= 3)
#   2. the S. mansoni 20.5-22.2 Mb genes with their S. japonicum position
#   3. S. japonicum Z 34-64 Mb in 1 Mb bins, with the S. mansoni Z position of
#      the genes in each bin

suppressMessages(library(tidyverse))
setwd(commandArgs(trailingOnly = TRUE)[1])

asm <- Sys.glob("assemblies/schistosoma_japonicum.*")
pos <- read_tsv(file.path(asm, "loci.bed"), col_names = c("contig", "s", "e", "gene"), col_types = "cddc") |>
    filter(contig == "CM046069.1") |> group_by(gene) |> summarise(sj_Mb = min(s) / 1e6)
fm <- read_tsv("strata/gene_depth.tsv", show_col_types = FALSE) |>
    filter(str_detect(sample_id, "^sjap_")) |>
    mutate(sex = str_sub(sample_id, 6, 6)) |>
    group_by(gene, chrom, start, sex) |> summarise(norm = mean(norm), .groups = "drop") |>
    pivot_wider(names_from = sex, values_from = norm) |> mutate(FM = F / M)
d <- inner_join(pos, fm, by = "gene") |> filter(M > 0.2) |> arrange(sj_Mb) |>
    mutate(state = if_else(FM < 0.7, "hemi", "dip"), run = consecutive_id(state))

cat("\nRuns along the S. japonicum Z (sm_ columns mix chromosomes; use the bins below for S. mansoni positions)\n")
d |> group_by(run, state) |>
    summarise(sj_from = min(sj_Mb), sj_to = max(sj_Mb), n = n(), FM = median(FM),
              sm_chrom = paste(unique(str_remove(chrom, "SM_V10_")), collapse = ","),
              sm_from = min(start) / 1e6, sm_to = max(start) / 1e6, .groups = "drop") |>
    filter(n >= 3) |> mutate(across(where(is.double), \(x) round(x, 2))) |> print(n = Inf, width = Inf)

cat("\nS. mansoni 20.5-22.2 Mb genes on the S. japonicum Z\n")
d |> filter(chrom == "SM_V10_Z", start > 20.5e6, start < 22.2e6) |>
    transmute(gene, sm_Mb = round(start / 1e6, 2), sj_Mb = round(sj_Mb, 2), FM = round(FM, 2)) |>
    arrange(sm_Mb) |> print(n = Inf)

cat("\nS. japonicum Z 34-64 Mb in 1 Mb bins, S. mansoni Z genes only\n")
d |> filter(chrom == "SM_V10_Z", sj_Mb >= 34, sj_Mb < 64) |>
    group_by(sj_bin = floor(sj_Mb)) |>
    summarise(n = n(), FM = round(median(FM), 2),
              sm_q10 = round(quantile(start, .1) / 1e6, 1), sm_median = round(median(start) / 1e6, 1),
              sm_q90 = round(quantile(start, .9) / 1e6, 1)) |>
    print(n = Inf)
