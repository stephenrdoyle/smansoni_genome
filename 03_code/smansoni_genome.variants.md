# Schistosoma mansoni V10 genome: variant analysis

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk


```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/V10/VARIANTS

ln -s ../REF/SM_V10.genome.preWBP18checked.fa
ln -s ../../SM_V9_21Feb.fa

ln -s /lustre/scratch118/infgen/team133/db22/software/snpEff/ALL.normed.MANSONI.snpeff.vcf


cat SM_V9_21Feb.fa > V9ref.mod.fa
sed -i 's/SM_V9_//g' V9ref.mod.fa
samtools faidx V9ref.mod.fa
```


## Converting variants on V9 to V10
- inital analyses performed on the V9 genome
- want to convert variants on V9 to V10, to avoid having to remap everything
- following protocol here: https://github.com/informationsea/transanno

```bash

module load minimap2/2.16=h84994c4_1-c1

bsub.py --threads 20 20 minimap "minimap2 -t 20 -cx asm5 --cs SM_V10.genome.preWBP18checked.fa V9ref.mod.fa  -o PAF_FILE.paf"


# Run transanno to create chain file
/nfs/users/nfs_s/sd21/lustre118_link/software/SNP_CALLING/transanno/target/release/transanno minimap2chain PAF_FILE.paf --output CHAINFILE.chain



bsub.py --queue long 10 liftvcf "/nfs/users/nfs_s/sd21/lustre118_link/software/SNP_CALLING/transanno/target/release/transanno liftvcf \
    --chain CHAINFILE.chain \
    --fail FAILED.vcf.gz \
    --output SUCCEEDED.vcf.gz \
    --new-assembly SM_V10.genome.preWBP18checked.fa \
    --original-assembly V9ref.mod.fa \
    --vcf ALL.normed.MANSONI.snpeff.vcf"



/nfs/users/nfs_s/sd21/lustre118_link/schistosoma_mansoni/NUCDIV/samples.male.list
```