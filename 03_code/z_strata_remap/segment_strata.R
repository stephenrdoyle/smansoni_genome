# Segment Z-chromosome coverage into strata, with bootstrap intervals on each boundary.
#
# usage: Rscript segment_strata.R WORKDIR [N_BOOTSTRAP] [ZW_DS_TSV]
#
# Every sample is placed on S. mansoni v10 Z coordinates through the
# S. mansoni genes lifted onto its own assembly (prepare_assembly.sh).
#
# Profiles segmented:
#   - each whole-genome sample: gene depth / median over genes from S. mansoni
#     autosomes, so a hemizygous region in a female sits at ~0.5
#   - each species with sexed samples (tier 1): female:male gene depth ratio,
#     replicates averaged, WGS and exome kept apart (capture bias cancels in
#     the ratio, so exome is only used this way)
#   - S. mansoni female:male in 10 kb windows, which needs no liftover
#
# Segmentation is circular binary segmentation (DNAcopy) on log2 values.
# Boundary intervals come from a block bootstrap: fitted segment medians plus
# residuals resampled in blocks (keeping local autocorrelation), re-segmented,
# and each original boundary matched to the nearest bootstrap boundary.
#
# ZW_DS_TSV, if given, is zw_gametologue_ds.tsv from zw_ds.py: Z-W gametologue
# dS is drawn as a track above the coverage panels.
#
# Writes WORKDIR/strata/: gene_depth.tsv, segments.tsv, breakpoints.tsv,
# sex_calls.tsv, interval_summary.tsv and z_strata_remap.pdf/png.

suppressPackageStartupMessages({
    library(tidyverse)
    library(DNAcopy)
    library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
workdir <- args[1]
n_boot <- if (length(args) > 1) as.integer(args[2]) else 500
ds_file <- if (length(args) > 2) args[3] else NA
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
outdir <- file.path(workdir, "strata")
dir.create(outdir, showWarnings = FALSE)
set.seed(11)

AUTOSOMES <- paste0("SM_V10_", 1:7)
Z <- "SM_V10_Z"
# boundaries used in the manuscript (Fig. 3a), drawn for comparison only
OLD_BOUNDARIES <- c(10680112, 23250000, 43743322, 51300000)
# intervals summarised per profile; 18.42-19.14 Mb is the gap holding the
# step in Z-W gametologue dS (~0.14 below, ~0.61 above)
INTERVALS <- tribble(
    ~interval,             ~from,      ~to,
    "PAR1",                0,          10680112,
    "African_a (dS low)",  10680112,   18420000,
    "African_b (dS high)", 19140000,   23250000,
    "ancestral",           23250000,   43743322,
    "Asian (Sj only)",     43743322,   51300000,
    "PAR2",                51300000,   Inf
)

samples <- read_tsv(file.path(script_dir, "samples.tsv"), show_col_types = FALSE,
                    col_types = cols(base_count = col_character()))
genes <- read_tsv(file.path(workdir, "ref", "sman_genes.tsv"), show_col_types = FALSE)


# --- segmentation helpers ------------------------------------------------------

# CBS, then merge segments too short to trust (< min_n points) or too similar
# to their neighbour to matter (< min_step in log2; a female's full-to-half
# step is 1, a mixed pool's is smaller but still > 0.2)
cbs <- function(y, pos, alpha = 0.001, min_n = 10, min_step = 0.2) {
    cna <- CNA(y, chrom = rep("Z", length(y)), maploc = pos, data.type = "logratio",
               sampleid = "s", presorted = TRUE)
    seg <- segment(smooth.CNA(cna), alpha = alpha, min.width = 5, verbose = 0)
    n <- seg$output$num.mark
    repeat {
        ends <- cumsum(n)
        means <- map2_dbl(c(1, head(ends, -1) + 1), ends, \(a, b) median(y[a:b]))
        if (length(n) == 1) break
        steps <- abs(diff(means))
        short <- which(n < min_n)
        if (length(short)) {
            k <- short[which.min(n[short])]
            # join the short segment to its more similar neighbour
            left <- if (k > 1) abs(means[k] - means[k - 1]) else Inf
            right <- if (k < length(n)) abs(means[k] - means[k + 1]) else Inf
            j <- if (left <= right) k - 1 else k + 1
        } else if (min(steps) < min_step) {
            k <- which.min(steps); j <- k + 1
        } else break
        n[min(k, j)] <- n[k] + n[j]
        n <- n[-max(k, j)]
    }
    ends <- cumsum(n)
    data.frame(num.mark = n, end_idx = ends,
               seg.mean = map2_dbl(c(1, head(ends, -1) + 1), ends, \(a, b) median(y[a:b])))
}

# boundary between segments k and k+1 sits between points end_idx[k] and end_idx[k]+1
boundary_index <- function(seg) head(seg$end_idx, -1)

block_resample <- function(r, block) {
    n <- length(r)
    starts <- sample.int(max(n - block + 1, 1), ceiling(n / block), replace = TRUE)
    idx <- as.vector(outer(0:(block - 1), starts, "+"))
    r[pmin(idx, n)][seq_len(n)]
}

segment_with_ci <- function(y, pos, block, n_boot, max_shift) {
    fit <- cbs(y, pos)
    fitted <- rep(fit$seg.mean, fit$num.mark)
    resid <- y - fitted
    b0 <- boundary_index(fit)
    hits <- matrix(NA_real_, nrow = n_boot, ncol = length(b0))
    if (length(b0) > 0) {
        for (b in seq_len(n_boot)) {
            yb <- fitted + block_resample(resid, block)
            bb <- boundary_index(cbs(yb, pos))
            for (k in seq_along(b0)) {
                if (length(bb) == 0) next
                d <- bb - b0[k]
                j <- which.min(abs(d))
                if (abs(d[j]) <= max_shift) hits[b, k] <- bb[j]
            }
        }
    }
    # a boundary at index i lies between pos[i] and pos[i + 1]
    mid <- function(i) (pos[i] + pos[i + 1]) / 2
    segs <- tibble(seg = seq_len(nrow(fit)), start = pos[c(1, b0 + 1)], end = pos[fit$end_idx],
                   n = fit$num.mark, log2_mean = fit$seg.mean, ratio = 2^fit$seg.mean)
    bps <- tibble(
        boundary = seq_along(b0),
        position = mid(b0),
        gap_start = pos[b0], gap_end = pos[b0 + 1],
        ci_low = map_dbl(seq_along(b0), \(k) { h <- hits[, k]; h <- h[!is.na(h)]; if (length(h)) mid(quantile(h, 0.025, type = 1)) else NA }),
        ci_high = map_dbl(seq_along(b0), \(k) { h <- hits[, k]; h <- h[!is.na(h)]; if (length(h)) mid(quantile(h, 0.975, type = 1)) else NA }),
        support = colMeans(!is.na(hits)),
        left_ratio = 2^fit$seg.mean[seq_along(b0)],
        right_ratio = 2^fit$seg.mean[seq_along(b0) + 1]
    )
    list(segments = segs, breakpoints = bps)
}

# female: a long segment near half depth; pooled: shallower dip; male: neither
call_sex <- function(segs) {
    long <- filter(segs, n >= 50)
    if (any(between(long$ratio, 0.35, 0.65))) "female"
    else if (any(between(long$ratio, 0.65, 0.85))) "pooled/mixed"
    else "male (no hemizygous region)"
}



# --- per-sample gene depth on S. mansoni coordinates ------------------------------

read_loci <- function(sample_id, species_dir, assembly) {
    asm <- file.path(workdir, "assemblies", paste(species_dir, assembly, sep = "."))
    f <- file.path(workdir, "samples", sample_id, paste0(sample_id, ".loci.regions.bed.gz"))
    if (!file.exists(f)) return(NULL)
    lift <- read_tsv(file.path(asm, "liftover.tsv"), show_col_types = FALSE)
    read_tsv(f, col_names = c("contig", "start", "end", "gene", "depth"), col_types = "cddcd") |>
        mutate(len = end - start) |>
        group_by(gene) |>
        summarise(depth = weighted.mean(depth, len), cds_len = sum(len), .groups = "drop") |>
        inner_join(filter(lift, ambiguous == 0) |> select(gene), by = "gene") |>
        mutate(sample_id = sample_id)
}

depth <- samples |>
    pmap(\(sample_id, species_dir, assembly, ...) read_loci(sample_id, species_dir, assembly)) |>
    bind_rows() |>
    inner_join(genes, by = "gene") |>
    filter(cds_len >= 300) |>
    left_join(samples |> select(sample_id, species, tier, platform, sex_metadata), by = "sample_id") |>
    group_by(sample_id) |>
    mutate(norm = depth / median(depth[chrom %in% AUTOSOMES & depth > 0])) |>
    ungroup()

write_tsv(depth, file.path(outdir, "gene_depth.tsv"))

# single-sample profiles (not exome: off-target genes make its depth uninterpretable alone)
single <- depth |>
    filter(chrom == Z, platform != "illumina_exome", norm > 0.05, norm < 3) |>
    transmute(profile = sample_id, species, kind = "sample", pos = start, value = norm)

# female:male profiles, per species and data type, replicates averaged
sexed <- depth |>
    filter(tier == 1, str_detect(sex_metadata, "^(fe)?male$")) |>
    mutate(sex = if_else(str_detect(sex_metadata, "^female"), "F", "M"),
           data = if_else(platform == "illumina_exome", "exome", "WGS"))

fm_genes <- sexed |>
    group_by(species, data, gene, chrom, start, sex) |>
    summarise(norm = mean(norm), raw = min(depth), n = n(), .groups = "drop") |>
    pivot_wider(names_from = sex, values_from = c(norm, raw, n)) |>
    filter(!is.na(norm_F), !is.na(norm_M), raw_F >= 5, raw_M >= 5, norm_M > 0.2) |>
    group_by(species, data) |>
    mutate(ratio = (norm_F / norm_M) / median((norm_F / norm_M)[chrom %in% AUTOSOMES])) |>
    ungroup()

fm <- fm_genes |>
    filter(chrom == Z, ratio > 0.05, ratio < 3) |>
    transmute(profile = paste0(species, " F/M (", data, ", ", n_F, "F ", n_M, "M)"),
              species, kind = "F/M genes", pos = start, value = ratio)

# S. mansoni female:male in 10 kb windows, each sample normalised by its autosomal median
win <- function(id) {
    f <- file.path(workdir, "samples", id, paste0(id, ".w10k.regions.bed.gz"))
    if (!file.exists(f)) return(NULL)
    read_tsv(f, col_names = c("chrom", "start", "end", "depth"), col_types = "cddd") |>
        mutate(norm = depth / median(depth[chrom %in% AUTOSOMES & depth > 0]), sample_id = id)
}
sman_sexed <- sexed |> filter(species == "Schistosoma mansoni", data == "WGS") |> distinct(sample_id, sex)
fm_win <- NULL
if (all(c("F", "M") %in% sman_sexed$sex)) {
    w <- map(sman_sexed$sample_id, win) |> bind_rows() |> inner_join(sman_sexed, by = "sample_id")
    if (nrow(w)) {
        fm_win <- w |>
            group_by(chrom, start, end, sex) |>
            summarise(norm = mean(norm), .groups = "drop") |>
            pivot_wider(names_from = sex, values_from = norm) |>
            filter(chrom == Z, F > 0.05, M > 0.2, F < 3, M < 3) |>
            transmute(profile = sprintf("Schistosoma mansoni F/M (10 kb windows, %dF %dM)",
                                        sum(sman_sexed$sex == "F"), sum(sman_sexed$sex == "M")),
                      species = "Schistosoma mansoni", kind = "F/M windows",
                      pos = (start + end) / 2, value = F / M)
    }
}

profiles <- bind_rows(fm_win, fm, single) |> arrange(profile, pos)


# --- segment each profile --------------------------------------------------------

results <- profiles |>
    group_split(profile) |>
    map(\(d) {
        windows <- d$kind[1] == "F/M windows"
        r <- segment_with_ci(log2(d$value), d$pos,
                             block = if (windows) 20 else 10, n_boot = n_boot,
                             max_shift = if (windows) 100 else 50)
        key <- distinct(d, profile, species, kind)
        list(segments = bind_cols(key, r$segments), breakpoints = bind_cols(key, r$breakpoints))
    })

segments <- map(results, "segments") |> bind_rows()
breakpoints <- map(results, "breakpoints") |> bind_rows() |>
    mutate(hemizygous_step = (pmin(left_ratio, right_ratio) < 0.75) & (pmax(left_ratio, right_ratio) >= 0.75),
           in_dS_gap = ci_high >= 18420000 & ci_low <= 19140000)

sex_calls <- segments |>
    filter(kind == "sample") |>
    group_by(sample_id = profile) |>
    group_modify(\(s, k) tibble(inferred_sex = call_sex(s))) |>
    ungroup() |>
    left_join(samples |> select(sample_id, species, tier, runs, sex_metadata), by = "sample_id")

# median value per interval, and whether the two halves of the old African stratum differ
interval_summary <- profiles |>
    mutate(interval = cut(pos, c(INTERVALS$from, Inf), labels = INTERVALS$interval, right = FALSE)) |>
    filter(!is.na(interval)) |>
    group_by(profile, species, kind, interval) |>
    summarise(n = n(), median = median(value), .groups = "drop_last") |>
    left_join(
        profiles |>
            group_by(profile) |>
            summarise(p_African_a_vs_b = {
                a <- value[pos >= 10680112 & pos < 18420000]
                b <- value[pos >= 19140000 & pos < 23250000]
                if (length(a) >= 5 && length(b) >= 5) wilcox.test(a, b)$p.value else NA_real_
            }),
        by = "profile") |>
    ungroup()

write_tsv(segments, file.path(outdir, "segments.tsv"))
write_tsv(breakpoints, file.path(outdir, "breakpoints.tsv"))
write_tsv(sex_calls, file.path(outdir, "sex_calls.tsv"))
write_tsv(interval_summary, file.path(outdir, "interval_summary.tsv"))
print(sex_calls, n = Inf)
print(breakpoints |> filter(hemizygous_step) |>
      select(profile, position, ci_low, ci_high, support, left_ratio, right_ratio, in_dS_gap), n = Inf)


# --- figure --------------------------------------------------------------------

species_order <- unique(samples$species)
panel_order <- profiles |>
    distinct(profile, species, kind) |>
    mutate(s = match(species, species_order), k = match(kind, c("F/M windows", "F/M genes", "sample"))) |>
    arrange(s, k, profile) |>
    pull(profile)
facet <- \(d) mutate(d, profile = factor(profile, levels = panel_order))

p <- ggplot(facet(profiles)) +
    annotate("rect", xmin = 18.42, xmax = 19.14, ymin = 0, ymax = 2, fill = "#4E79C2", alpha = 0.15) +
    geom_rect(data = facet(filter(breakpoints, hemizygous_step)),
              aes(xmin = ci_low / 1e6, xmax = ci_high / 1e6, ymin = 0, ymax = 2), fill = "#C24E9C", alpha = 0.25) +
    geom_vline(xintercept = OLD_BOUNDARIES / 1e6, linetype = "dotted", colour = "grey40") +
    geom_hline(yintercept = c(0.5, 1), linetype = "dashed", colour = "grey70") +
    geom_point(aes(pos / 1e6, value), size = 0.3, alpha = 0.5) +
    geom_segment(data = facet(segments),
                 aes(x = start / 1e6, xend = end / 1e6, y = ratio, yend = ratio), colour = "#C24E9C", linewidth = 0.8) +
    facet_grid(profile ~ ., switch = "y", labeller = label_wrap_gen(30)) +
    coord_cartesian(ylim = c(0, 2)) +
    theme_bw() +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 0, size = 6), panel.spacing = unit(0.1, "lines")) +
    labs(x = "S. mansoni v10 Z chromosome position (Mb)",
         y = "Normalised coverage, or female:male ratio",
         caption = paste("Lines: CBS segments. Pink: 95% bootstrap interval of each full-to-half boundary.",
                         "Blue: 18.42-19.14 Mb gametologue dS step. Dotted: boundaries in the submitted Fig. 3a.", sep = "\n"))

if (!is.na(ds_file) && file.exists(ds_file)) {
    ds <- read_tsv(ds_file, show_col_types = FALSE) |> filter(Z_chrom == Z)
    p_ds <- ggplot(ds, aes(Z_start / 1e6, pmin(dS_YN00, 3), shape = W_type)) +
        annotate("rect", xmin = 18.42, xmax = 19.14, ymin = 0, ymax = 3, fill = "#4E79C2", alpha = 0.15) +
        geom_vline(xintercept = OLD_BOUNDARIES / 1e6, linetype = "dotted", colour = "grey40") +
        geom_point() +
        coord_cartesian(xlim = range(profiles$pos) / 1e6) +
        theme_bw() + theme(legend.position = "top", axis.title.x = element_blank()) +
        labs(y = "Z-W dS (YN00,\ncapped at 3)", shape = NULL)
    n_panels <- n_distinct(profiles$profile)
    p <- p_ds / p + plot_layout(heights = c(1.5, n_panels))
}

n_panels <- n_distinct(profiles$profile)
ggsave(file.path(outdir, "z_strata_remap.pdf"), p, height = 0.8 * n_panels + 2.5, width = 7.5, limitsize = FALSE)
ggsave(file.path(outdir, "z_strata_remap.png"), p, height = 0.8 * n_panels + 2.5, width = 7.5, dpi = 200, limitsize = FALSE)
