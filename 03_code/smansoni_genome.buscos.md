# Schistosoma mansoni V10 genome: busco analysis

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk


## SM_V10
```bash
# genome
bsub.py --threads 10 20 busco_smv10_genome "busco --in /lustre/scratch118/infgen/team333/sd21/schistosoma_mansoni/V10/REF/SM_V10.genome.preWBP18checked.fa --out SM_V10_genome --mode genome --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:89.4%[S:88.6%,D:0.8%],F:4.3%,M:6.3%,n:255
	228	Complete BUSCOs (C)
	226	Complete and single-copy BUSCOs (S)
	2	Complete and duplicated BUSCOs (D)
	11	Fragmented BUSCOs (F)
	16	Missing BUSCOs (M)
	255	Total BUSCO groups searched

# protein
bsub.py --threads 10 20 busco_smv10_protein "busco --in /lustre/scratch118/infgen/team333/sd21/schistosoma_mansoni/V10/REF/SM_V10.genome.preWBP18checked.proteins.fa --out SM_V10_protein --mode protein --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:96.1%[S:91.0%,D:5.1%],F:0.8%,M:3.1%,n:255
	245	Complete BUSCOs (C)
	232	Complete and single-copy BUSCOs (S)
	13	Complete and duplicated BUSCOs (D)
	2	Fragmented BUSCOs (F)
	8	Missing BUSCOs (M)
	255	Total BUSCO groups searched
```

## SM_V5
```bash
# genome
bsub.py --threads 10 20 busco_smv5_genome "busco --in /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/V5_vs_V10/schistosoma_mansoni.PRJEA36577.WBPS1.genomic.fa --out SM_V5_genome --mode genome --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:88.6%[S:88.2%,D:0.4%],F:4.3%,M:7.1%,n:255
	226	Complete BUSCOs (C)
	225	Complete and single-copy BUSCOs (S)
	1	Complete and duplicated BUSCOs (D)
	11	Fragmented BUSCOs (F)
	18	Missing BUSCOs (M)
	255	Total BUSCO groups searched

# protein
bsub.py --threads 10 20 busco_smv5_protein "busco --in /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/V5_vs_V10/schistosoma_mansoni.PRJEA36577.WBPS1.protein.fa --out SM_V5_protein --mode protein --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:87.5%[S:80.8%,D:6.7%],F:6.7%,M:5.8%,n:255
	223	Complete BUSCOs (C)
	206	Complete and single-copy BUSCOs (S)
	17	Complete and duplicated BUSCOs (D)
	17	Fragmented BUSCOs (F)
	15	Missing BUSCOs (M)
	255	Total BUSCO groups searched
```

## SBOS_v1.1_pilon
```bash
# genome
fastaq enumerate_names --suffix _scaffold /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/SBOS_v1.1_pilon.fasta SBOS_v1.1_pilon.renamed.fasta

/nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/SBOS_v1.1_pilon.fasta
bsub.py --threads 10 20 busco_SBOS_v1.1_pilon_genome "busco --in SBOS_v1.1_pilon.renamed.fasta --out SBOS_v1.1_pilon_genome --mode genome --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:83.9%[S:80.4%,D:3.5%],F:4.7%,M:11.4%,n:255
	214	Complete BUSCOs (C)
	205	Complete and single-copy BUSCOs (S)
	9	Complete and duplicated BUSCOs (D)
	12	Fragmented BUSCOs (F)
	29	Missing BUSCOs (M)
	255	Total BUSCO groups searched

```

## sbovis_GCA_003958945
```bash
fastaq enumerate_names --suffix _scaffold /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/sbovis_GCA_003958945.fa sbovis_GCA_003958945.renamed.fa

bsub.py --threads 10 20 busco_sbovis_GCA_003958945_genome "busco --in sbovis_GCA_003958945.renamed.fa --out sbovis_GCA_003958945_genome --mode genome --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:84.3%[S:82.7%,D:1.6%],F:7.8%,M:7.9%,n:255
	215	Complete BUSCOs (C)
	211	Complete and single-copy BUSCOs (S)
	4	Complete and duplicated BUSCOs (D)
	20	Fragmented BUSCOs (F)
	20	Missing BUSCOs (M)
	255	Total BUSCO groups searched
```

## scurasoni_GCA_900618015
```bash
fastaq enumerate_names --suffix _scaffold /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/scurasoni_GCA_900618015.1.fa scurasoni_GCA_900618015.1.renamed.fa

bsub.py --threads 10 20 busco_scurasoni_GCA_900618015_genome "busco --in scurasoni_GCA_900618015.1.renamed.fa --out scurasoni_GCA_900618015_genome --mode genome --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:56.5%[S:56.5%,D:0.0%],F:27.1%,M:16.4%,n:255
	144	Complete BUSCOs (C)
	144	Complete and single-copy BUSCOs (S)
	0	Complete and duplicated BUSCOs (D)
	69	Fragmented BUSCOs (F)
	42	Missing BUSCOs (M)
	255	Total BUSCO groups searched

```

## shaematobium_GCA_000699445.3
```bash
fastaq enumerate_names --suffix _scaffold /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/shaematobium_GCA_000699445.3.fasta shaematobium_GCA_000699445.3.renamed.fasta

bsub.py --threads 10 20 busco_shaematobium_GCA_000699445.3_genome "busco --in shaematobium_GCA_000699445.3.renamed.fasta --out shaematobium_GCA_000699445.3_genome --mode genome --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:88.3%[S:86.7%,D:1.6%],F:5.5%,M:6.2%,n:255
	225	Complete BUSCOs (C)
	221	Complete and single-copy BUSCOs (S)
	4	Complete and duplicated BUSCOs (D)
	14	Fragmented BUSCOs (F)
	16	Missing BUSCOs (M)
	255	Total BUSCO groups searched

```


## sjaponicum_GCA_006368765.1
```bash
fastaq enumerate_names --suffix _scaffold /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/sjaponicum_GCA_006368765.1.fasta sjaponicum_GCA_006368765.1.renamed.fasta

bsub.py --threads 10 20 busco_sjaponicum_GCA_006368765.1_genome "busco --in sjaponicum_GCA_006368765.1.renamed.fasta --out sjaponicum_GCA_006368765.1_genome --mode genome --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:83.5%[S:82.7%,D:0.8%],F:7.5%,M:9.0%,n:255
	213	Complete BUSCOs (C)
	211	Complete and single-copy BUSCOs (S)
	2	Complete and duplicated BUSCOs (D)
	19	Fragmented BUSCOs (F)
	23	Missing BUSCOs (M)
	255	Total BUSCO groups searched

```

## smargrebowiei_GCA_900618395.1
```bash
fastaq enumerate_names --suffix _scaffold /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/smargrebowiei_GCA_900618395.1.fasta smargrebowiei_GCA_900618395.1.renamed.fasta

bsub.py --threads 10 20 busco_smargrebowiei_GCA_900618395.1_genome "busco --in smargrebowiei_GCA_900618395.1.renamed.fasta --out smargrebowiei_GCA_900618395.1_genome --mode genome --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:69.8%[S:69.8%,D:0.0%],F:20.0%,M:10.2%,n:255
	178	Complete BUSCOs (C)
	178	Complete and single-copy BUSCOs (S)
	0	Complete and duplicated BUSCOs (D)
	51	Fragmented BUSCOs (F)
	26	Missing BUSCOs (M)
	255	Total BUSCO groups searched

```

## smattheei_GCA_900617995.1
```bash
fastaq enumerate_names --suffix _scaffold /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/smattheei_GCA_900617995.1.fasta smattheei_GCA_900617995.1.renamed.fasta

bsub.py --threads 10 20 busco_smattheei_GCA_900617995.1_genome "busco --in smattheei_GCA_900617995.1.renamed.fasta --out smattheei_GCA_900617995.1_genome --mode genome --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:54.5%[S:54.5%,D:0.0%],F:31.8%,M:13.7%,n:255
	139	Complete BUSCOs (C)
	139	Complete and single-copy BUSCOs (S)
	0	Complete and duplicated BUSCOs (D)
	81	Fragmented BUSCOs (F)
	35	Missing BUSCOs (M)
	255	Total BUSCO groups searched

```

## schistosoma_rodhaini
```bash
fastaq enumerate_names --suffix _scaffold /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/schistosoma_rodhaini.PRJEB526.WBPS17.genomic.fa schistosoma_rodhaini.PRJEB526.WBPS17.genomic.renamed.fa

bsub.py --threads 10 20 busco_srodhaini.PRJEB526.WBPS17_genome "busco --in schistosoma_rodhaini.PRJEB526.WBPS17.genomic.renamed.fa --out srodhaini.PRJEB526.WBPS17_genome --mode genome --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:56.5%[S:56.5%,D:0.0%],F:29.8%,M:13.7%,n:255
	144	Complete BUSCOs (C)
	144	Complete and single-copy BUSCOs (S)
	0	Complete and duplicated BUSCOs (D)
	76	Fragmented BUSCOs (F)
	35	Missing BUSCOs (M)
	255	Total BUSCO groups searched

```

## schistosoma_japonicum
```bash 
# paper: https://www.sciencedirect.com/science/article/pii/S2211124722003904

fastaq enumerate_names --suffix _scaffold /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/sjaponicum_GCA_021461655.1_ASM2146165v1_genomic.fasta sjaponicum_GCA_021461655.1_ASM2146165v1_genomic.renamed.fasta

gffread -y sjaponicum_GCA_021461655.1_ASM2146165v1.protein.fa -g /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/sjaponicum_GCA_021461655.1_ASM2146165v1_genomic.fasta /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/sjaponicum_GCA_021461655.1_ASM2146165v1_genomic.gff


# genome
bsub.py --threads 10 20 busco_sjaponicum_GCA_021461655.1_ASM2146165v1_genome "busco --in /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/BUSCO/sjaponicum_GCA_021461655.1_ASM2146165v1_genomic.renamed.fasta --out sjaponicum_GCA_021461655.1_ASM2146165v1_genome --mode genome --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:87.9%[S:85.9%,D:2.0%],F:7.1%,M:5.0%,n:255
	224	Complete BUSCOs (C)
	219	Complete and single-copy BUSCOs (S)
	5	Complete and duplicated BUSCOs (D)
	18	Fragmented BUSCOs (F)
	13	Missing BUSCOs (M)
	255	Total BUSCO groups searched

# protein
bsub.py --threads 10 20 busco_sjaponicum_GCA_021461655.1_ASM2146165v1_protein "busco --in sjaponicum_GCA_021461655.1_ASM2146165v1.protein.fa --out sjaponicum_GCA_021461655.1_ASM2146165v1_protein --mode protein --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:89.8%[S:30.2%,D:59.6%],F:4.7%,M:5.5%,n:255
	229	Complete BUSCOs (C)
	77	Complete and single-copy BUSCOs (S)
	152	Complete and duplicated BUSCOs (D)
	12	Fragmented BUSCOs (F)
	14	Missing BUSCOs (M)
	255	Total BUSCO groups searched
```



## schistosoma mansoni V9 - WBPSv17
```bash

fastaq enumerate_names --suffix _scaffold /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/schistosoma_mansoni.PRJEA36577.WBPS17.genomic.fa schistosoma_mansoni.PRJEA36577.WBPS17.renamed.fasta

gffread -y schistosoma_mansoni.PRJEA36577.WBPS17.protein.fa -g /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/schistosoma_mansoni.PRJEA36577.WBPS17.genomic.fa /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/schistosoma_mansoni.PRJEA36577.WBPS17.annotations.gff3

# genome
bsub.py --threads 10 20 busco_schistosoma_mansoni.PRJEA36577.WBPS17_genome "busco --in /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/BUSCO/schistosoma_mansoni.PRJEA36577.WBPS17.renamed.fasta --out schistosoma_mansoni.PRJEA36577.WBPS17_genome --mode genome --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:89.4%[S:88.6%,D:0.8%],F:4.3%,M:6.3%,n:255
	228	Complete BUSCOs (C)
	226	Complete and single-copy BUSCOs (S)
	2	Complete and duplicated BUSCOs (D)
	11	Fragmented BUSCOs (F)
	16	Missing BUSCOs (M)
	255	Total BUSCO groups searched

# protein
bsub.py --threads 10 20 busco_schistosoma_mansoni.PRJEA36577.WBPS17_protein "busco --in schistosoma_mansoni.PRJEA36577.WBPS17.protein.fa --out schistosoma_mansoni.PRJEA36577.WBPS17_protein --mode protein --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:94.9%[S:77.3%,D:17.6%],F:2.0%,M:3.1%,n:255
	242	Complete BUSCOs (C)
	197	Complete and single-copy BUSCOs (S)
	45	Complete and duplicated BUSCOs (D)
	5	Fragmented BUSCOs (F)
	8	Missing BUSCOs (M)
	255	Total BUSCO groups searched
```


## schistosoma mansoni V7 - WBPSv15
```bash

fastaq enumerate_names --suffix _scaffold /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/schistosoma_mansoni.PRJEA36577.WBPS15.genomic.fa schistosoma_mansoni.PRJEA36577.WBPS15.renamed.fasta

gffread -y schistosoma_mansoni.PRJEA36577.WBPS15.protein.fa -g /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/schistosoma_mansoni.PRJEA36577.WBPS15.genomic.fa /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/REF/OTHER_SCHISTO_GENOMES/schistosoma_mansoni.PRJEA36577.WBPS15.annotations.gff3

# genome
bsub.py --threads 10 20 busco_schistosoma_mansoni.PRJEA36577.WBPS15_genome "busco --in /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/BUSCO/schistosoma_mansoni.PRJEA36577.WBPS15.renamed.fasta --out schistosoma_mansoni.PRJEA36577.WBPS15_genome --mode genome --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"


# protein
bsub.py --threads 10 20 busco_schistosoma_mansoni.PRJEA36577.WBPS15_protein "busco --in schistosoma_mansoni.PRJEA36577.WBPS15.protein.fa --out schistosoma_mansoni.PRJEA36577.WBPS15_protein --mode protein --lineage_dataset ~sd21/lustre118_link/databases/busco/eukaryota_odb10 --cpu 10 -f"

	C:95.3%[S:76.5%,D:18.8%],F:2.0%,M:2.7%,n:255
	243	Complete BUSCOs (C)
	195	Complete and single-copy BUSCOs (S)
	48	Complete and duplicated BUSCOs (D)
	5	Fragmented BUSCOs (F)
	7	Missing BUSCOs (M)
	255	Total BUSCO groups searched
```
