# S. mansoni W (SM_V10_WSR) depth in the female and male, windows vs gene loci.
#
# usage: Rscript wsr_check.R WORKDIR
#
# Most of the WSR is repeat, so few reads pass MAPQ 20 and the window median is
# near 0 even in the female. W genes are single copy and should sit at ~0.5 in
# the female and ~0 in the male.

suppressMessages(library(tidyverse))
setwd(commandArgs(trailingOnly = TRUE)[1])

rd <- function(s, w) read_tsv(sprintf("samples/%s/%s.%s.regions.bed.gz", s, s, w),
    col_names = FALSE, col_types = cols(.default = "c")) %>%
    mutate(depth = as.numeric(.[[ncol(.)]]))
for (s in c("sman_F", "sman_M")) {
    w <- rd(s, "w10k")
    a <- median(w$depth[str_detect(w$X1, "^SM_V10_[1-7]$")])
    cat("\n", s, ": WSR 10 kb windows, depth / autosome median, quantiles 10-90%\n")
    print(round(quantile(w$depth[w$X1 == "SM_V10_WSR"] / a, c(.1, .25, .5, .75, .9)), 2))
    l <- rd(s, "loci") %>% filter(X1 == "SM_V10_WSR")
    cat(s, ": WSR gene loci (n =", nrow(l), "), median depth / autosome median =", round(median(l$depth / a), 2), "\n")
}
