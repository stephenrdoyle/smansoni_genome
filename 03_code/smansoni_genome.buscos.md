# Schistosoma mansoni V10 genome: busco analysis

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk


## Comparison of V9 genome and proteome
```bash
# genome
/nfs/users/nfs_s/sd21/bash_scripts/run_busco_metazoa.sh busco_buddenborg_2021 SM_V9_ENA.fa

#INFO	C:73.8%[S:72.0%,D:1.8%],F:5.8%,M:20.4%,n:978
#INFO	722 Complete BUSCOs (C)
#INFO	704 Complete and single-copy BUSCOs (S)
#INFO	18 Complete and duplicated BUSCOs (D)
#INFO	57 Fragmented BUSCOs (F)
#INFO	199 Missing BUSCOs (M)
#INFO	978 Total BUSCO groups searched

# proteins
/nfs/users/nfs_s/sd21/lustre118_link/software/ASSEMBLY_QC/busco_v3/scripts/run_BUSCO.py --in proteins.fa --mode proteins --out protein --lineage_path /nfs/users/nfs_s/sd21/lustre118_link/databases/busco/metazoa_odb9/ --force

#INFO	C:85.4%[S:68.5%,D:16.9%],F:3.8%,M:10.8%,n:978
#INFO	835 Complete BUSCOs (C)
#INFO	670 Complete and single-copy BUSCOs (S)
#INFO	165 Complete and duplicated BUSCOs (D)
#INFO	37 Fragmented BUSCOs (F)
#INFO	106 Missing BUSCOs (M)
#INFO	978 Total BUSCO groups searched
```
- where "run_busco_metasoa.sh" is
```bash
#!/usr/local/bin/bash

# run nematode busco

export AUGUSTUS_CONFIG_PATH=/nfs/users/nfs_s/sd21/software/augustus-3.2.1/config
export PATH=$PATH:/nfs/users/nfs_s/sd21/software/augustus-3.2.1/bin/

PREFIX=$1
REF=$2


/nfs/users/nfs_s/sd21/lustre118_link/software/ASSEMBLY_QC/busco_v3/scripts/run_BUSCO.py \
--in ${REF} \
--out ${PREFIX}_busco3.02.metazoa \
--mode genome \
--lineage_path  ~sd21/lustre118_link/databases/busco/metazoa_odb9/ \
--species caenorhabditis \
--cpu 7 --tarzip --force --restart --long --blast_single_core \
--tmp_path ${REF}.tmp
```



## 
```bash
# genome
/nfs/users/nfs_s/sd21/bash_scripts/run_busco_metazoa.sh busco_protasio_2012 schistosoma_mansoni.PRJEA36577.WBPS1.genomic.fa

#INFO	C:72.5%[S:71.1%,D:1.4%],F:6.1%,M:21.4%,n:978
#INFO	709 Complete BUSCOs (C)
#INFO	695 Complete and single-copy BUSCOs (S)
#INFO	14 Complete and duplicated BUSCOs (D)
#INFO	60 Fragmented BUSCOs (F)
#INFO	209 Missing BUSCOs (M)
#INFO	978 Total BUSCO groups searched

# proteins
/nfs/users/nfs_s/sd21/lustre118_link/software/ASSEMBLY_QC/busco_v3/scripts/run_BUSCO.py --in proteins.fa --mode proteins --out protein --lineage_path /nfs/users/nfs_s/sd21/lustre118_link/databases/busco/metazoa_odb9/ --force


#INFO	C:80.4%[S:71.7%,D:8.7%],F:6.6%,M:13.0%,n:978
#INFO	786 Complete BUSCOs (C)
#INFO	701 Complete and single-copy BUSCOs (S)
#INFO	85 Complete and duplicated BUSCOs (D)
#INFO	65 Fragmented BUSCOs (F)
#INFO	127 Missing BUSCOs (M)
#INFO	978 Total BUSCO groups searched
```
