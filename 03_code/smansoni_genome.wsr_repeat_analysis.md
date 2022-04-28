# W specific repeat analysis

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk


# TO-DO - add nucmer analyses for repeats



```bash
# bedtools to calculate coverage in repeats
bedtools multicov -bams 6520_5_1_sorted.bam -bed combined_repeats_SD220223.bed > combined_repeats_SD220223.coverage

# multiply the number of reads by read length, and divide by length of array
cat combined_repeats_SD220223.coverage | awk '{print $5*100/($3-$2)}
```



# plot of WSR repeats
```R
library(tidyverse)
library(RColorBrewer)

data <- read.table("combined_repeats_SD22022_to_v9.WSR.coords", header=F)

# set colours
n <- 39
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(data, aes(V1,fill=V13)) +
     geom_histogram(binwidth=10000) +
     theme_bw() +
     scale_fill_manual(values=col_vector) +
     labs(fill="Repeat", x="Chromosome position (bp)", y="Repeat count")

ggsave("WSR_repeats.pdf", width=7, height=5)
```



## Normalising coverage to estimate copy number of gene family arrays: IPSE and OMEGA

### IPSE
```bash

# ISPE coverage using samtools
cat mb_IPSE.gff | awk '{print $1,$4,$5,$9}' OFS="\t" > mb_IPSE.bed

samtools bedcov mb_IPSE.bed 6520_5_1_sorted.bam > mb_IPSE.coverage2

cat mb_IPSE.coverage2 | awk '{print $1,$2,$3,$4,$5/($3-$2)}' OFS="\t" | datamash median 5
#>53.7555

53.7555/44.6 = 1.21x higher coverage in region

16 * 1.21 = 20 likely copies of IPSE based on normalising coverage
```

### OMEGA
```
# OMEGA coverage using samtools
cat mb_omega.gff | awk '{print $1,$4,$5,$9}' OFS="\t" > mb_omega.bed

samtools bedcov mb_omega.bed 6520_5_1_sorted.bam > mb_omega.coverage

cat mb_omega.coverage | awk '{print $1,$2,$3,$4,$5/($3-$2)}' OFS="\t" | datamash median 5

#> 78.35955/44.6 = 1.76x higher coverage in region

#> 5 copies * 1.76x = 8.8 likely copies of omega based on normalising coverage
#--- note - there are 7 copies of omega in this region, but only 5 are in a high coverage area. Two copies are not. See Supplementary Figure 1C
#--- therefore: 8.8 normalised copies + 2 addition copies =  10.8 copies, ~11 copies reported in the paper.

```




```bash

library(tidyverse)
library(patchwork)

data <- read.table("V7W_vs_V9_mum.filtered.coords", header=F)
rep_data <- read.table("W_repeat_data.txt", header=T, sep="\t")
genes_data <- read.table("../WSR_genes.bed", header=F)


# dotplot
ggplot(data) + geom_segment(aes(x=V1, y=V1+V3, xend=V2, yend=V1+V4, col=V11), size=1) + theme_bw() + labs(x="SM_V9_WSR", y="W scaffolds from V7", col="W scaffolds from V7")






data2 <- data %>% group_by(V11) %>% mutate(group_id = median(V1))
rep_data2 <- rep_data %>% group_by(Repeat.name) %>% mutate(group_id = median(V1))

plot_v7w <- ggplot(data2) +
     geom_segment(aes(x=V1, xend=V2, y=as.factor(group_id), yend=as.factor(group_id), col=V11), size=3) +
     theme_bw() + theme(legend.position="none", axis.text.y = element_blank()) +
     labs(title="V7 W scaffolds", y="", x="Genomic position on V9 W scaffold (bp)")

plot_rep <- ggplot(rep_data2) +
     geom_segment(aes(x=START, xend=END, y=reorder(group_id,START), yend=reorder(group_id,START), col=Repeat.name), size=3) +
     theme_bw() + theme(legend.position="none", axis.text.y = element_blank()) +
     labs(title="Repeats", y="", x="Genomic position on V9 W scaffold (bp)")

plot_genes <- ggplot(genes_data) +
     geom_segment(aes(x=V2, xend=V3, y=V4, yend=V4, col=V5), size=3) +
     theme_bw() + theme(legend.position="none") +
     labs(title="V9 W genes", y="Strand", x="Genomic position on V9 W scaffold (bp)")

plot_v7w + plot_rep + plot_genes + plot_layout(ncol=1)



```
