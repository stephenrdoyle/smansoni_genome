# W specific repeat analysis

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk


# TO-DO - add nucmer analyses for repeats



```bash
# bedtools to calculate coverage in repeats
bedtools multicov -bams 6520_5_1_sorted.bam -bed combined_repeats_SD220223.bed > combined_repeats_SD220223.coverage

# multiply the number of reads by read length, and divide by length of array
cat combined_repeats_SD220223.coverage | awk '{print $5*100/($3-$2)}'

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




```R

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


```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/ALT_CONTIGS

minimap2 -Y --secondary=no -x asm20 SM_V9_21Feb.fa SM_V7_Wcontigs.fa -o V7w_to_V9_minimap_v2.paf    

grep "SM_V9_WSR" V7w_to_V9_minimap_v2.paf | sort -k8n | cut -f1-11 > V7w_to_V9_minimap_v2.txt
```



```R
library(tidyverse)

data <- read.table("V7w_to_V9_minimap_v2.txt", header=F)
data2 <- data %>% group_by(V1) %>% mutate(group_id = min(V8))

plot_V7W <- ggplot(data2) +
     geom_segment(aes(x=V8/1e6, xend=V9/1e6, y=reorder(V1,(group_id)), yend=reorder(V1,(group_id)), col=as.factor(V5)), size=3) +
     labs(x="Genomic position on WSR in V9 assembly (Mb)", y="W scaffolds from V7 assembly") +
     theme_bw()
```


```bash
junction_subreads.bed

wspecific_pilon	947666	947667
wspecific_pilon	950283	950284
wspecific_pilon	1970116	1970117
wspecific_pilon	4754002	4754003
wspecific_pilon	6634601	6634602
wspecific_pilon	6641096	6641097
wspecific_pilon	8062142	8062143
wspecific_pilon	13127900	13127901
wspecific_pilon	13754066	13754067
wspecific_pilon	14901326	14901327
wspecific_pilon	15344882	15344883
wspecific_pilon	16081118	16081119
wspecific_pilon	16164607	16164608
wspecific_pilon	17231003	17231004
wspecific_pilon	17338728	17338729


samtools view -b -L junction_subreads.bed out.sorted.markdup.bam | samtools fasta - > junction_subreads.fasta

minimap2 -Y --secondary=no -x map-pb SM_V9_21Feb.fa junction_subreads.fasta > junction_subreads.paf

cat junction_subreads.paf | grep "SM_V9_WSR" > junction_subreads.WSR.paf
```


```R
library(tidyverse)

reads <- read.table("junction_subreads.WSR.paf", header=F, sep="\t")
reads2 <- reads %>% group_by(V1) %>% mutate(group_id = min(V8))

plot_reads <- ggplot(reads2) +
     geom_segment(aes(x=V8/1e6, xend=V9/1e6, y=reorder(V1,(group_id)), yend=reorder(V1,(group_id)), col=as.factor(V5)), size=3) +
     labs(x="Genomic position on WSR in V9 assembly (Mb)", y="W scaffolds from V7 assembly") +
     theme_bw()



ggplot() +     
     geom_segment(aes(x=reads2$V8/1e6, xend=reads2$V9/1e6, y=as.factor(1), yend=as.factor(1), col=as.factor(reads2$V5)), size=3) + geom_segment(aes(x=data2$V8/1e6, xend=data2$V9/1e6, y=reorder(data2$V1,(data2$group_id)), yend=reorder(data2$V1,(data2$group_id)), col=as.factor(data2$V5)), size=3) +
     geom_vline(aes(xintercept=reads2$V8/1e6))

vline=reads2$V8
ggplot(data2) +
     geom_vline(xintercept=vline/1e6) +
     geom_segment(aes(x=V8/1e6, xend=V9/1e6, y=V1, yend=V1, col=as.factor(V5)), size=3)
```


```bash
# new list from Alan
while read NAME; do
     samtools faidx pb_reads.fasta $NAME;
done < junction_subreads.list > junction_subreads.v2.fasta

minimap2 -Y --secondary=no -x map-pb SM_V9_21Feb.fa junction_subreads.v2.fasta > junction_subreads.v2.paf

cat junction_subreads.v2.paf | grep "SM_V9_WSR" > junction_subreads.v2.WSR.paf

```

```R
library(tidyverse)

reads <- read.table("junction_subreads.v2.WSR.paf", header=F, sep="\t")
reads2 <- reads %>% group_by(V1) %>% mutate(group_id = min(V8))

plot_reads <- ggplot(reads2) +
     geom_segment(aes(x=V8/1e6, xend=V9/1e6, y=reorder(V1,(group_id)), yend=reorder(V1,(group_id)), col=as.factor(V5)), size=3) +
     labs(x="Genomic position on WSR in V9 assembly (Mb)", y="W scaffolds from V7 assembly") +
     theme_bw()
plot_reads

```
