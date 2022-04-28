


cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/COV

cp ../SM_V9_21Feb.fa .

samtools faidx SM_V9_21Feb.fa

cut -f1,2 SM_V9_21Feb.fa.fai > SM_V9.genome

bedtools makewindows -g SM_V9.genome -w 100 > SM_V9.100b.bed



# nucleotide frequency
bedtools nuc -fi SM_V9_21Feb.fa -bed SM_V9.500b.bed > SM_V9.500b.nucfreq

#--- extract chromosome 1
grep 'SM_V9_1' SM_V9.500b.nucfreq > SM_V9_1.500b.nucfreq

#--- extract chromosome 2
grep 'SM_V9_2' SM_V9.500b.nucfreq > SM_V9_2.500b.nucfreq



# bam coverage
ln -s ../BAMS/6520_5_1_sorted.bam
ln -s ../BAMS/6520_5_1_sorted.bam.bai

samtools bedcov SM_V9.100b.bed 6520_5_1_sorted.bam > SM_V9.100b.coverage &

samtools bedcov SM_V9.500b.bed 6520_5_1_sorted.bam > SM_V9.500b.coverage &

#--- extract chromosome 1
grep 'SM_V9_1' SM_V9.500b.coverage > SM_V9_1.500b.coverage

#--- extract chromosome 2
grep 'SM_V9_2' SM_V9.500b.coverage > SM_V9_2.500b.coverage



# genes
cd ~/lustre118_link/schistosoma_mansoni/REPEATS
cat mb_IPSE.gff | awk -F '[\t;]' '{print $1,$4,$5,$7,$9,$11}' OFS="\t" | sed 's/ID=//g' > mb_IPSE.bed

cat mb_omega.gff | awk -F '[\t;]' '{print $1,$4,$5,$7,$9,$11}' OFS="\t" | sed 's/ID=//g' > mb_omega.bed



library(tidyverse)
data <- read.table("mb_omega.bed", header=F)
data_summary <- data %>% group_by(V5) %>% summarise(min=min(V2), max=max(V3), strand=V4)


ggplot() +
     geom_rect(data=data, aes(xmin=V2, xmax=V3, ymin=0, ymax=1, col=as.factor(V6))) +
     geom_segment(data=data_summary, aes(x=min, xend=max, y=0.5, yend=0.5)) +
     theme_classic()



library(tidyverse)
data <- read.table("mb_IPSE.bed", header=F)
data_summary <- data %>% group_by(V5) %>% summarise(min=min(V2), max=max(V3), strand=V4)


ggplot() +
     geom_rect(data=data, aes(xmin=V2, xmax=V3, ymin=0, ymax=1, col=as.factor(V6))) +
     geom_segment(data=data_summary, aes(x=min, xend=max, y=0.5, yend=0.5)) +
     theme_classic() +
     guides(y = "none")




# NOR plot
nuc <- read.table("../COV/SM_V9_2.500b.nucfreq", header=F)
nuc_median <- median(nuc$V5)

cov <- read.table("../COV/SM_V9_2.500b.coverage", header=F)
cov_median <- median(cov$V4/500)


plot_nuc <- ggplot(nuc) + geom_line(aes(V2, V5), size=0.75) + xlim(39.4e6, 39.75e6) + theme_bw() + guides(x = "none") + labs(x="", y="GC Content", title="NOR") + geom_hline(yintercept=nuc_median, linetype="dashed")

plot_cov <- ggplot(cov) + geom_line(aes(V2, V4/500), size=0.75) + xlim(39.4e6, 39.75e6) + theme_bw() + guides(x = "none") + labs(x="", y="Coverage")  + geom_hline(yintercept=cov_median, linetype="dashed")

plot_feature <- ggplot() +
     geom_rect(aes(xmin=39.445e6, xmax=39.47e6, ymin=0.5, ymax=1.5), fill="red") +
     geom_rect(aes(xmin=39.52e6, xmax=39.675e6, ymin=0.5, ymax=1.5), fill="cornflowerblue") +
     geom_rect(aes(xmin=39.545e6, xmax=39.555e6, ymin=0.5, ymax=1.5), fill="green") +
     theme_classic() +
     guides(y = "none") + theme(legend.position="none") +
     xlim(39.4e6, 39.75e6) + labs(x="Genome position (bp)", y="")

plot_nuc + plot_cov + plot_feature + plot_layout(ncol=1, heights=c(3,3,1))


# IPSE plot

library(tidyverse)
library(patchwork)


gene <- read.table("mb_IPSE.bed", header=F)
gene_summary <- data %>% group_by(V5) %>% summarise(min=min(V2), max=max(V3), strand=V4)

nuc <- read.table("../COV/SM_V9_1.500b.nucfreq", header=F)
nuc_median <- median(nuc$V5)
cov <- read.table("../COV/SM_V9_1.500b.coverage", header=F)
cov_median <- median(cov$V4/500)


plot_nuc <- ggplot(nuc) + geom_line(aes(V2, V5), size=0.75) + xlim(7090000, 7280000) + theme_bw() + guides(x = "none") + labs(x="", y="GC Content", title="IPSE") + geom_hline(yintercept=nuc_median, linetype="dashed")

plot_cov <- ggplot(cov) + geom_line(aes(V2, V4/500), size=0.75) + xlim(7090000, 7280000) + theme_bw() + ylim(0,200) + guides(x = "none") + labs(x="", y="Coverage")  + geom_hline(yintercept=cov_median, linetype="dashed")

plot_gene <- ggplot() +
     geom_segment(data=gene_summary, aes(x=min, xend=max, y=1, yend=1), size=1, col="grey") +
     geom_rect(data=gene, aes(xmin=V2, xmax=V3, ymin=0.5, ymax=1.5, col=as.factor(V6))) +
     theme_classic() +
     guides(y = "none") + theme(legend.position="none") +
     xlim(7090000, 7280000) + labs(x="Genome position (bp)", y="")

plot_nuc + plot_cov + plot_gene + plot_layout(ncol=1, heights=c(3,3,1))



# OMEGA plot

library(tidyverse)
library(patchwork)


gene <- read.table("mb_omega.bed", header=F)
gene_summary <- data %>% group_by(V5) %>% summarise(min=min(V2), max=max(V3), strand=V4)

nuc <- read.table("../COV/SM_V9_1.500b.nucfreq", header=F)
nuc_median <- median(nuc$V5)
cov <- read.table("../COV/SM_V9_1.500b.coverage", header=F)
cov_median <- median(cov$V4/500)


plot_nuc <- ggplot(nuc) + geom_line(aes(V2, V5), size=0.75) + xlim(3620000, 3880000) + theme_bw() + guides(x = "none") + labs(x="", y="GC Content", title="omega") + geom_hline(yintercept=nuc_median, linetype="dashed")

plot_cov <- ggplot(cov) + geom_line(aes(V2, V4/500), size=0.75) + xlim(3620000, 3880000) + theme_bw() + ylim(0,200) + guides(x = "none") + labs(x="", y="Coverage")  + geom_hline(yintercept=cov_median, linetype="dashed")

plot_gene <- ggplot() +
     geom_segment(data=gene_summary, aes(x=min, xend=max, y=1, yend=1), size=1, col="grey") +
     geom_rect(data=gene, aes(xmin=V2, xmax=V3, ymin=0.5, ymax=1.5, col=as.factor(V6))) + 
     theme_classic() +
     guides(y = "none") + theme(legend.position="none") +
     xlim(3620000, 3880000) + labs(x="Genome position (bp)", y="")

plot_nuc + plot_cov + plot_gene + plot_layout(ncol=1, heights=c(3,3,1))
