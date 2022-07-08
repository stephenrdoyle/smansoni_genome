### Coverage plots of repetitive regions: NOR, IPSE and Omega

```bash
# working directory
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/COV

# get the reference and index it, using the index file to make a genome file for bedtools
cp ../SM_V9_21Feb.fa .

samtools faidx SM_V9_21Feb.fa

cut -f1,2 SM_V9_21Feb.fa.fai > SM_V9.genome

# extract bed coordinates for the genome in 500 bp windows
bedtools makewindows -g SM_V9.genome -w 500 > SM_V9.500b.bed



# nucleotide frequency
bedtools nuc -fi SM_V9_21Feb.fa -bed SM_V9.500b.bed > SM_V9.500b.nucfreq

#--- extract nucleotide frequency data for chromosome 1
grep 'SM_V9_1' SM_V9.500b.nucfreq > SM_V9_1.500b.nucfreq

#--- extract nucleotide frequency data for  chromosome 2
grep 'SM_V9_2' SM_V9.500b.nucfreq > SM_V9_2.500b.nucfreq



# get mapped data from the single female
ln -s ../BAMS/6520_5_1_sorted.bam
ln -s ../BAMS/6520_5_1_sorted.bam.bai

# calculate genome coverage in 500 bp windows
samtools bedcov SM_V9.500b.bed 6520_5_1_sorted.bam > SM_V9.500b.coverage &

#--- extract coverage data for chromosome 1
grep 'SM_V9_1' SM_V9.500b.coverage > SM_V9_1.500b.coverage

#--- extract  coverage data for chromosome 2
grep 'SM_V9_2' SM_V9.500b.coverage > SM_V9_2.500b.coverage

#--- extract  coverage data for chromosome 2
grep 'SM_V9_Z' SM_V9.500b.coverage > SM_V9_Z.500b.coverage



# extract coordinate data for CDS for IPSE and omega
cd ~/lustre118_link/schistosoma_mansoni/REPEATS
cat mb_IPSE.gff | awk -F '[\t;]' '{print $1,$4,$5,$7,$9,$11}' OFS="\t" | sed 's/ID=//g' > mb_IPSE.bed

cat mb_omega.gff | awk -F '[\t;]' '{print $1,$4,$5,$7,$9,$11}' OFS="\t" | sed 's/ID=//g' > mb_omega.bed

```




### NOR plot
```R
# load libraries
library(tidyverse)
library(patchwork)

# load nucleotide frequency data, and calculate median GC
nuc <- read.table("../COV/SM_V9_2.500b.nucfreq", header=F)
nuc_median <- median(nuc$V5)

# load coverage data, and calculate median coverage
cov <- read.table("../COV/SM_V9_2.500b.coverage", header=F)
cov_median <- median(cov$V4/500)


# make some plots
plot_nuc <- ggplot(nuc) +
     geom_line(aes(V2, V5), size=0.75) +
     xlim(39.4e6, 39.75e6) +
     theme_bw() + guides(x = "none") +
     labs(x="", y="GC Content", title="NOR") +
     geom_hline(yintercept=nuc_median, linetype="dashed")

plot_cov <- ggplot(cov) +
     geom_line(aes(V2, V4/500), size=0.75) +
     xlim(39.4e6, 39.75e6) +
     theme_bw() + guides(x = "none") +
     labs(x="", y="Coverage")  +
     geom_hline(yintercept=cov_median, linetype="dashed")

plot_feature <- ggplot() +
     geom_rect(aes(xmin=39.445e6, xmax=39.47e6, ymin=0.5, ymax=1.5), fill="red") +
     geom_rect(aes(xmin=39.52e6, xmax=39.675e6, ymin=0.5, ymax=1.5), fill="cornflowerblue") +
     geom_rect(aes(xmin=39.545e6, xmax=39.555e6, ymin=0.5, ymax=1.5), fill="green") +
     theme_classic() +
     guides(y = "none") + theme(legend.position="none") +
     xlim(39.4e6, 39.75e6) +
     labs(x="Genome position (bp)", y="")


# combine into multipanel
plot_nuc + plot_cov + plot_feature + plot_layout(ncol=1, heights=c(3,3,1))
```

### IPSE plot
```R
# load libraries
library(tidyverse)
library(patchwork)

# load gene data, and determine min and max coords per gene
gene <- read.table("mb_IPSE.bed", header=F)
gene_summary <- data %>% group_by(V5) %>% summarise(min=min(V2), max=max(V3), strand=V4)

# load nucleotide frequency data, and calculate median GC
nuc <- read.table("../COV/SM_V9_1.500b.nucfreq", header=F)
nuc_median <- median(nuc$V5)

# load coverage data, and calculate median coverage
cov <- read.table("../COV/SM_V9_1.500b.coverage", header=F)
cov_median <- median(cov$V4/500)

# make some plots
plot_nuc <- ggplot(nuc) +
     geom_line(aes(V2, V5), size=0.75) +
     xlim(7090000, 7280000) +
     theme_bw() + guides(x = "none") + labs(x="", y="GC Content", title="IPSE") +
     geom_hline(yintercept=nuc_median, linetype="dashed")

plot_cov <- ggplot(cov) +
     geom_line(aes(V2, V4/500), size=0.75) +
     xlim(7090000, 7280000) + ylim(0,200) +
     theme_bw() + guides(x = "none") +
     labs(x="", y="Coverage")  +
     geom_hline(yintercept=cov_median, linetype="dashed")

plot_gene <- ggplot() +
     geom_segment(data=gene_summary, aes(x=min, xend=max, y=1, yend=1), size=1, col="grey") +
     geom_rect(data=gene, aes(xmin=V2, xmax=V3, ymin=0.5, ymax=1.5, col=as.factor(V6))) +
     theme_classic() +
     guides(y = "none") + theme(legend.position="none") +
     xlim(7090000, 7280000) +
     labs(x="Genome position (bp)", y="")

# combine into multipanel
plot_nuc + plot_cov + plot_gene + plot_layout(ncol=1, heights=c(3,3,1))
```


# OMEGA plot
```R
# load libraries
library(tidyverse)
library(patchwork)


# load gene data, and determine min and max coords per gene
gene <- read.table("mb_omega.bed", header=F)
gene_summary <- data %>% group_by(V5) %>% summarise(min=min(V2), max=max(V3), strand=V4)


# load nucleotide frequency data, and calculate median GC
nuc <- read.table("../COV/SM_V9_1.500b.nucfreq", header=F)
nuc_median <- median(nuc$V5)

# load coverage data, and calculate median coverage
cov <- read.table("../COV/SM_V9_1.500b.coverage", header=F)
cov_median <- median(cov$V4/500)

# make some plots
plot_nuc <- ggplot(nuc) +
     geom_line(aes(V2, V5), size=0.75) +
     xlim(3620000, 3880000) +
     theme_bw() + guides(x = "none") + labs(x="", y="GC Content", title="omega") +
     geom_hline(yintercept=nuc_median, linetype="dashed")

plot_cov <- ggplot(cov) +
     geom_line(aes(V2, V4/500), size=0.75) +
     xlim(3620000, 3880000) + ylim(0,200) +
     theme_bw() + guides(x = "none") + labs(x="", y="Coverage")  +
     geom_hline(yintercept=cov_median, linetype="dashed")

plot_gene <- ggplot() +
     geom_segment(data=gene_summary, aes(x=min, xend=max, y=1, yend=1), size=1, col="grey") +
     geom_rect(data=gene, aes(xmin=V2, xmax=V3, ymin=0.5, ymax=1.5, col=as.factor(V6))) +
     theme_classic() +
     guides(y = "none") + theme(legend.position="none") +
     xlim(3620000, 3880000) +
     labs(x="Genome position (bp)", y="")

# combine into multipanel
plot_nuc + plot_cov + plot_gene + plot_layout(ncol=1, heights=c(3,3,1))
```




cp /lustre/scratch118/infgen/team133/alt/SCHISTO/V5/schistosoma_mansoni.PRJEA36577.WBPS7.genomic.fa SM_V5.fa
samtools faidx SM_V5.fa

cut -f1,2 SM_V5.fa.fai > SM_V5.genome

# extract bed coordinates for the genome in 500 bp windows
bedtools makewindows -g SM_V5.genome -w 500 > SM_V5.500b.bed
bedtools makewindows -g SM_V5.genome -w 10000 > SM_V5.10000b.bed

bedtools makewindows -g SM_V9.genome -w 10000 > SM_V9.10000b.bed
bedtools makewindows -g SM_V9.genome -w 50000 > SM_V9.50000b.bed

# get mapped data from the single female
ln -s /lustre/scratch118/infgen/team133/alt/SCHISTO/V5/6520_5_coverage/out.sorted.markdup.bam V5.6520_5_1_sorted.bam
ln -s /lustre/scratch118/infgen/team133/alt/SCHISTO/V5/6520_5_coverage/out.sorted.markdup.bam.bai V5.6520_5_1_sorted.bam.bai

# calculate genome coverage in 500 bp windows
samtools bedcov SM_V5.500b.bed V5.6520_5_1_sorted.bam > SM_V5.500b.coverage &
samtools bedcov SM_V5.10000b.bed V5.6520_5_1_sorted.bam > SM_V5.10000b.coverage &

samtools bedcov SM_V9.10000b.bed 6520_5_1_sorted.bam > SM_V9.10000b.coverage &

#--- extract coverage data for Z chromosome
grep 'Smp.Chr_ZW' SM_V5.500b.coverage > SM_V5_ZW.500b.coverage
grep 'Smp.Chr_ZW' SM_V5.10000b.coverage > SM_V5_ZW.10000b.coverage


grep 'SM_V9_PAR1' SM_V9.10000b.coverage > SM_V9_Z.10000b.coverage
grep 'SM_V9_ZSR' SM_V9.10000b.coverage >> SM_V9_Z.10000b.coverage
grep 'SM_V9_PAR2' SM_V9.10000b.coverage >> SM_V9_Z.10000b.coverage

grep 'SM_V9_PAR1' SM_V9.50000b.coverage > SM_V9_Z.50000b.coverage
grep 'SM_V9_ZSR' SM_V9.50000b.coverage >> SM_V9_Z.50000b.coverage
grep 'SM_V9_PAR2' SM_V9.50000b.coverage >> SM_V9_Z.50000b.coverage
```

```R
library(tidyverse)
library(patchwork)


cov_V5 <- read.table("SM_V5_ZW.10000b.coverage", header=F)
cov_V5$cov_median <- median(cov_V5$V4/10000)
cov_V5$genome <- "V5"
cov_V5 <- mutate(cov_V5, id = row_number())


cov_V9 <- read.table("SM_V9_Z.10000b.coverage", header=F)
cov_V9$cov_median <- median(cov_V9$V4/10000)
cov_V9$genome <- "V9"
cov_V9 <- mutate(cov_V9, id = row_number())

cov_V9 <- read.table("SM_V9_Z.50000b.coverage", header=F)
cov_V9$cov_median <- median(cov_V9$V4/50000)
cov_V9$genome <- "V9"
cov_V9 <- mutate(cov_V9, id = row_number())



cov_V5_median <- median(cov_V5$V4/10000)
plot_cov_v5 <- ggplot(cov_V5) +
     geom_point(aes(1:nrow(cov_V5)*10000, V4/10000), size=0.1) +
     theme_bw() +
     labs(x="Genomic position (bp)", y="Coverage")  +
     ylim(0,100) +
     geom_hline(yintercept=cov_V5_median, linetype="dashed")

cov_V9_median <- median(cov_V9$V4/10000)
plot_cov_v9 <- ggplot(cov_V9) +
     geom_point(aes(1:nrow(cov_V9)*10000, V4/10000), size=0.1) +
     theme_bw() +
     labs(x="Genomic position (bp)", y="Coverage")  +
     ylim(0,100) +
     geom_hline(yintercept=cov_V9_median, linetype="dashed") +
     facet_grid(scales="free", space = "free")


plot_cov_v5 + plot_cov_v9 + plot_layout(ncol=1)
```

# comparison track between v5 and v9
```bash
samtools faidx SM_V9_21Feb.fa SM_V9_PAR1 > SM_V9_Z.fa
samtools faidx SM_V9_21Feb.fa SM_V9_ZSR >> SM_V9_Z.fa
samtools faidx SM_V9_21Feb.fa SM_V9_PAR2 >> SM_V9_Z.fa
cat SM_V9_Z.fa | sed -e 's/>SM_V9_ZSR//g' -e 's/>SM_V9_PAR2//g' -e 's/_PAR1/_Z/g' > SM_V9_Z.merge.fa


samtools faidx SM_V5.fa Smp.Chr_ZW > SM_V5_Z.fa


bsub.py 10 promer_ZV5_v_ZV9 "promer -p ZV5_v_ZV9 SM_V5_Z.fa SM_V9_Z.fa"
bsub.py 10 promer_ZV5_v_ZV9.2 "promer -p ZV5_v_ZV9 SM_V5_Z.fa SM_V9_Z.merge.fa"

show-coords -k -lTH -I 98 -L 5000 ZV5_v_ZV9.delta | awk '{print $1,$2,$3,$4,$5,$6,$7,$10,$11,$14,$15}' OFS="\t" > ZV5_v_ZV9.filtered.coords
```


<!--
# minimap2 to dotter
```
minimap2 -cx asm20 --cs SM_V9_Z.fa SM_V5_Z.fa > ZV5_v_ZV9.paf
minimap2 SM_V9_Z.fa SM_V5_Z.fa > ZV5_v_ZV9.paf

cat ZV5_v_ZV9.paf | /nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_ASSEMBLY/minimap/utils/bin/layout > layout.txt

/software/R-4.0.3/bin/Rscript /nfs/users/nfs_s/sd21/lustre118_link/software/GENOME_ASSEMBLY/minimap/utils/plot/plotLayout.R -f layout.txt -p out.pdf
``` -->

# load libraries
library(tidyverse)
library(patchwork)

# V5 coverage
cov_V5 <- read.table("SM_V5_ZW.10000b.coverage", header=F)
cov_V5$cov_median <- median(cov_V5$V4/10000)
cov_V5$genome <- "V5"
cov_V5 <- mutate(cov_V5, id = row_number())

cov_V5_median <- median(cov_V5$V4/10000)

plot_cov_v5 <- ggplot(cov_V5) +
     geom_point(aes(1:nrow(cov_V5)*10000/1e6, V4/10000), size=0.05, alpha=0.25) +
     theme_bw() +
     labs(x="Genomic position (bp)", y="Coverage")  +
     ylim(0,100) +
     geom_hline(yintercept=cov_V5_median, linetype="dashed") +
     facet_grid(scales="free", space = "free")

# V9 coverage     
cov_V9 <- read.table("SM_V9_Z.10000b.coverage", header=F)
cov_V9$cov_median <- median(cov_V9$V4/10000)
cov_V9$genome <- "V9"
cov_V9 <- mutate(cov_V9, id = row_number())

cov_V9_median <- median(cov_V9$V4/10000)
plot_cov_v9 <- ggplot(cov_V9) +
     geom_point(aes(1:nrow(cov_V9)*10000/1e6, V4/10000), size=0.05, alpha=0.25) +
     theme_bw() +
     labs(x="Genomic position (bp)", y="Coverage")  +
     ylim(0,100) +
     xlim(0,90) +
     geom_hline(yintercept=cov_V9_median, linetype="dashed")

# didnt use this dotter in the end
<!-- # dotter
data <- read.table("ZV5_v_ZV9.filtered.coords", header=F)
data <- mutate(data, col = if_else(V4>V3 & V2>V1,1,if_else(V4<V3 & V2<V1,1,0)))

plot_dotter <- ggplot(data) +
     geom_segment(aes(x=V3/1e6,xend=V4/1e6,y=V1/1e6,yend=V2/1e6, col=as.factor(col)), lineend = "round", linejoin = "round", size=1, show.legend = FALSE) +
     theme(legend.position="none") + theme_bw() +
     xlim(0,90) +
     guides(x="none", y="none") + labs(x="", y="")

# make multipanel
plot_cov_v5 + plot_dotter + plot_spacer() + plot_cov_v9 + plot_layout(ncol=2, widths = c(1, 4), heights = c(4, 1)) -->

plot_cov_v5 + plot_spacer() + plot_cov_v9 + plot_layout(ncol=1)

ggsave("updated_figure2_coverageplots.pdf", height=5, width=7)
```










# recreating HIC maps of chromosome 1 and WSR
```bash
bsub.py --queue yesterday --threads 7 10 minimap_sr "minimap2 -t 7 -ax sr ../V9/smansoni_buddenborg2021_genome.fa  32442_4#1_1.fastq.gz 32442_4#1_2.fastq.gz -o minimap.sam"

samtools view -H minimap.sam > tmp

cat minimap.sam | grep "WSR" > WSR.sam
cat minimap.sam | grep "SM_V9_1" > SM_V9_1.sam

samtools view -H minimap.sam > tmp
cat tmp WSR.sam > tmp2; mv tmp2 WSR.sam
cat tmp SM_V9_1.sam > tmp2; mv tmp2 SM_V9_1.sam

samtools sort -n -o WSR.sorted.bam WSR.bam
samtools sort -n -o SM_V9_1.sorted.sam SM_V9_1.sam

samtools view -f 64 WSR.sorted.bam | awk '{if($3=="SM_V9_WSR" && $7=="=") print $1,$4,$5,$9}' OFS="\t" > WSR.sorted.txt
samtools view -f 64 SM_V9_1.sorted.sam | awk '{if($3=="SM_V9_1" && $7=="=") print $1,$4,$5,$9}' OFS="\t" > SM_V9_1.sorted.txt
```

```R
library(tidyverse)
library(viridis)
library(patchwork)

#data <- read.table("WSR.sorted.txt", header=F, sep="\t")
data <- read.table("SM_V9_1.sorted.txt", header=F, sep="\t")
data <- filter(data, V3>0)
data2 <- mutate(data, midpoint = if_else(V4>0, V2+(V4/2), V2-(-V4/2)))
data2 <- mutate(data2, distance = abs(V4))

#ggplot(data2) + geom_hex(aes(V2+(V4/2),V4),bins = 75) + scale_fill_viridis(direction = -1)
plot_1 <- ggplot(data2) + geom_hex(aes(midpoint,distance),bins = 75) + scale_fill_viridis(direction = -1)
plot_2 <- ggplot(data2) + geom_hex(aes(midpoint,log10(distance)),bins = 75) + scale_fill_viridis(direction = -1)

plot_1 + plot_2 + plot_layout(ncol=1)




data <- read.table("WSR.sorted.txt", header=F, sep="\t")
#data <- read.table("SM_V9_1.sorted.txt", header=F, sep="\t")
data <- filter(data, V3>0)
data2 <- mutate(data, midpoint = if_else(V4>0, V2+(V4/2), V2-(-V4/2)))
data2 <- mutate(data2, distance = abs(V4))

#ggplot(data2) + geom_hex(aes(V2+(V4/2),V4),bins = 75) + scale_fill_viridis(direction = -1)
plot_3 <- ggplot(data2) + geom_hex(aes(midpoint,distance),bins = 75) + scale_fill_viridis(direction = -1)
plot_4 <- ggplot(data2) + geom_hex(aes(midpoint,log10(distance)),bins = 75) + scale_fill_viridis(direction = -1)

plot_1 + plot_2 + plot_3 + plot_4 + plot_layout(ncol=1)
```
