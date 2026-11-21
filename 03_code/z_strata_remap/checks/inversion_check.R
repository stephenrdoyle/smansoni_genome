# Polarise the S. mansoni / S. japonicum inversion of S. mansoni Z 20.9-23.25 Mb.
#
# usage: Rscript inversion_check.R WORKDIR
#
# Gene order only (loci.bed, unambiguous placements), no reads. Segments by
# S. mansoni position: A = left flank (20.0 Mb to Smp_084270), B = the inverted
# segment (Smp_084270 to 23.25 Mb), C = start of the ancestral stratum
# (23.25-24.5 Mb). For each assembly, on the contig holding most of these genes:
#   1. each segment's position and orientation (Spearman rho against S. mansoni)
#   2. the distance between segment ends (median of the 3 end genes):
#      S. mansoni order puts A_end next to B_start and B_end next to C_start;
#      S. japonicum order puts A_end next to B_end and B_start next to C_start
# S. haematobium and S. japonicum are the controls.

suppressMessages(library(tidyverse))
setwd(commandArgs(trailingOnly = TRUE)[1])

genes <- read_tsv("ref/sman_genes.tsv", show_col_types = FALSE) |>
    filter(chrom == "SM_V10_Z", start >= 20e6, start < 24.5e6) |> arrange(start)
cut <- genes$start[genes$gene == "Smp_084270"]
genes <- genes |> mutate(seg = case_when(start < cut ~ "A", start < 23.25e6 ~ "B", TRUE ~ "C"))

asms <- unlist(map(c("schistosoma_haematobium", "schistosoma_japonicum", "trichobilharzia", "heterobilharzia"),
                   \(x) Sys.glob(paste0("assemblies/", x, "*"))))
for (asm in asms) {
    ok <- read_tsv(file.path(asm, "liftover.tsv"), show_col_types = FALSE) |> filter(ambiguous == 0) |> pull(gene)
    pos <- read_tsv(file.path(asm, "loci.bed"), col_names = c("contig", "s", "e", "gene"), col_types = "cddc") |>
        filter(gene %in% ok) |> group_by(gene, contig) |> summarise(og = min(s), .groups = "drop")
    d <- inner_join(genes, pos, by = "gene")
    main <- names(which.max(table(d$contig)))
    cat("\n", basename(asm), ": ", nrow(d), " genes placed, ", sum(d$contig == main), " on ", main, "\n", sep = "")
    d <- filter(d, contig == main)
    d |> group_by(seg) |>
        summarise(n = n(), og_from = round(min(og) / 1e6, 2), og_median = round(median(og) / 1e6, 2),
                  og_to = round(max(og) / 1e6, 2), rho = round(cor(start, og, method = "spearman"), 2)) |> print()

    A <- filter(d, seg == "A"); B <- filter(d, seg == "B"); C <- filter(d, seg == "C")
    end3 <- \(x, first) median((if (first) head(x, 3) else tail(x, 3))$og)
    j <- c(A_end = end3(A, FALSE), B_start = end3(B, TRUE), B_end = end3(B, FALSE), C_start = end3(C, TRUE))
    kb <- \(a, b) round(abs(j[[a]] - j[[b]]) / 1e3)
    cat(sprintf("mansoni-type: A_end-B_start %d kb, B_end-C_start %d kb | japonicum-type: A_end-B_end %d kb, B_start-C_start %d kb\n",
        kb("A_end", "B_start"), kb("B_end", "C_start"), kb("A_end", "B_end"), kb("B_start", "C_start")))
    cat("end genes:", paste(tail(A$gene, 3), collapse = ","), "|", paste(head(B$gene, 3), collapse = ","), "...",
        paste(tail(B$gene, 3), collapse = ","), "|", paste(head(C$gene, 3), collapse = ","), "\n")
}
