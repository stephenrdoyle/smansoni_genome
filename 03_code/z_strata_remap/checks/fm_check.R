# S. mansoni female:male depth in 10 kb windows, summarised per region of v10.
#
# usage: Rscript fm_check.R WORKDIR
#
# Each sample is normalised to its own autosomal median. Expected: autosomes and
# Z PARs ~1, every ZSR block F ~0.5 / M ~1, WSR M ~0. The ZSR is split at the
# 18.42-19.14 Mb dS step ("other" is that gap).

suppressMessages(library(tidyverse))
setwd(commandArgs(trailingOnly = TRUE)[1])

rd <- function(s) read_tsv(sprintf("samples/%s/%s.w10k.regions.bed.gz", s, s),
    col_names = c("chr", "start", "end", "depth"), col_types = "cddd")
d <- inner_join(rd("sman_F"), rd("sman_M"), by = c("chr", "start", "end"), suffix = c("_F", "_M"))
aut <- str_detect(d$chr, "^SM_V10_[1-7]$")
d <- d %>% mutate(F = depth_F / median(depth_F[aut]), M = depth_M / median(depth_M[aut]),
    region = case_when(
        aut ~ "autosomes",
        chr == "SM_V10_WSR" ~ "WSR",
        chr == "SM_V10_Z" & (end <= 10680112 | start > 43743322) ~ "Z PARs",
        chr == "SM_V10_Z" & end <= 18420000 ~ "ZSR African 10.7-18.4",
        chr == "SM_V10_Z" & start >= 19140000 & end <= 23250000 ~ "ZSR African 19.1-23.25",
        chr == "SM_V10_Z" & start > 23250000 ~ "ZSR ancestral",
        TRUE ~ "other"))
d %>% group_by(region) %>%
    summarise(windows = n(), F = median(F), M = median(M), F_over_M = median(F / M, na.rm = TRUE)) %>%
    print(width = Inf)
