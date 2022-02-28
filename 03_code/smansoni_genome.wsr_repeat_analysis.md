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
library(RColorBrewer)
n <- 39
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(a, aes(V1,fill=V13)) + geom_histogram(bins=1000) + theme_bw() + scale_fill_manual(values=col_vector) + labs(fill="Repeat", x="Chromosome position (bp)", y="Repeat count")
```





```bash
# IPSE coverage
#bedtools multicov  -bams 6520_5_1_sorted.bam -bed mb_IPSE.gff > mb_IPSE.coverage

#cat mb_IPSE.coverage | awk '{print $9,$10*100/($5-$4)}' OFS="\t" | awk '{print $1,$2,$2/44.6}' OFS="\t" | grep -v "colour=13" | datamash median 2
#79.476


# OMEGA coverage
#bedtools multicov  -bams 6520_5_1_sorted.bam -bed mb_omega.gff > mb_omega.coverage

#cat mb_omega.coverage | awk '{print $9,$10*100/($5-$4)}' OFS="\t" | awk '{print $1,$2,$2/44.6}' OFS="\t" | datamash median 2
#151.5255




# ISPE coverage using samtools
cat mb_IPSE.gff | awk '{print $1,$4,$5,$9}' OFS="\t" > mb_IPSE.bed

samtools bedcov mb_IPSE.bed 6520_5_1_sorted.bam > mb_IPSE.coverage2

cat mb_IPSE.coverage2 | awk '{print $1,$2,$3,$4,$5/($3-$2)}' OFS="\t" | datamash median 5
#>53.7555

53.7555/44.6 = 1.21x higher coverage in region

16 * 1.21 = 20 likely copies of IPSE based on normalising coverage

# OMEGA coverage using samtools
cat mb_omega.gff | awk '{print $1,$4,$5,$9}' OFS="\t" > mb_omega.bed

samtools bedcov mb_omega.bed 6520_5_1_sorted.bam > mb_omega.coverage

cat mb_omega.coverage | awk '{print $1,$2,$3,$4,$5/($3-$2)}' OFS="\t" | datamash median 5

#> 78.35955/44.6 = 1.76x higher coverage in region

#> 5 * 1.76 = 8.8 likely copies of omega based on normalising coverage
```
