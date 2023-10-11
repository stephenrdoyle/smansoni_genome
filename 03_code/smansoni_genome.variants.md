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


### Calculated sequencing coverage of PAR and ZSR to determine sex
- PAR1 to ZSR ratio will be sufficient to infer the sex of the samples
- sex: 
    - male = ZZ = ZSR/PAR1 = ~1
    - female = ZW = ZSR/PAR1 = ~0.5


```bash
cd COVERAGE

for i in $( cd ../ ; ls -1d 2*  ) ; do 
    mosdepth --use-median --fast-mode --by z-windows.bed --no-per-base ${i} ../${i}/${i}.bam; 
done


echo -e sample"\t"PAR1_coverage"\t"ZSR_coverage"\t"Inferred_sex > coverage.stats

for i in $( cd ../ ; ls -1d 2*  ) ; do 
    # Extract the first and last values from the file
    PAR1=$(zcat "${i}.regions.bed.gz" | cut -f4 | head -n1)
    ZSR=$(zcat "${i}.regions.bed.gz" | cut -f4 | tail -n1)

    # Perform the division and store the result in SEX
    SEX=$(echo "scale=2; $ZSR / $PAR1" | bc)

    # Check if SEX is greater than 0.8
    if (( $(echo "$SEX > 0.8" | bc -l) )); then
        SEX="male"
    else
        SEX="female"
    fi

    # Append the data to coverage.stats
    echo -e "${i}\t${PAR1}\t${ZSR}\t${SEX}" >> coverage.stats
done


cat coverage.stats | awk '{if($2>=10 && $4=="female") print}' | wc -l
#> 54 females

cat coverage.stats | awk '{if($2>=10 && $4=="male") print}' | wc -l
#> 58 males


# select the females with sufficient coverage for SNP calling
cat coverage.stats | awk '{if($2>=10 && $4=="female") print $1}'  > female_samples_w_coverage.list

# select samples with sufficient coverage for SNP calling
cat coverage.stats | awk '{if($2>=10) print $1}' > samples_w_coverage.list
```

### Call variants
```bash

cd /nfs/users/nfs_s/sd21/lustre_link/schistosoma_mansoni/V10/VARIANTS

cp COVERAGE/female_samples_w_coverage.list .

# while read SAMPLE; do 
#     ls -1 $PWD/$SAMPLE/$SAMPLE.bam; 
#     done < female_samples_w_coverage.list > bam.list


while read SAMPLE; do 
    ls -1 $PWD/$SAMPLE/$SAMPLE.bam; 
    done < samples_w_coverage.list > bam.list


 ./run_gatk_hc.sh

 ```

- where "run_gatk_hc.sh" is:
 ```bash
#!/bin/bash

#-------------------------------------------------------------------------------
# run_gatk_hc_allsites.sh
#-------------------------------------------------------------------------------

# stephen doyle
# Jan 2023

# Export environment variables
export PREFIX=smansoni_female_allsites  # prefix for output files
export REFERENCE=/nfs/users/nfs_s/sd21/lustre_link/schistosoma_mansoni/V10/VARIANTS/SM_V10.genome.preWBP18checked.fa  # path to reference genome
export BAM_LIST=/nfs/users/nfs_s/sd21/lustre_link/schistosoma_mansoni/V10/VARIANTS/bam.list  # path to list of BAM files

# Load GATK module
module load gatk/4.1.4.1

# Define file locations
export LOG_FILES="${PWD}/gatk_hc_${PREFIX}/LOG_FILES"  # directory for log files
export REFERENCE_FILES="${PWD}/gatk_hc_${PREFIX}/REFERENCE_FILES"  # directory for reference files
export GATK_HC_GVCFs="${PWD}/gatk_hc_${PREFIX}/GATK_HC_GVCFs"  # directory for GATK HC GVCF files
export GATK_HC_MERGED="${PWD}/gatk_hc_${PREFIX}/GATK_HC_MERGED"  # directory for merged haplotype caller files

# Create directories if they don't exist
[ -d ${LOG_FILES} ] || mkdir -p ${LOG_FILES}
[ -d ${REFERENCE_FILES} ] || mkdir -p ${REFERENCE_FILES}
[ -d ${GATK_HC_GVCFs} ] || mkdir -p ${GATK_HC_GVCFs}
[ -d ${GATK_HC_MERGED} ] || mkdir -p ${GATK_HC_MERGED}



# Save current script in run folder to reproduce the exact output
cp ${PWD}/run_gatk_hc.sh ${PWD}/gatk_hc_${PREFIX}/commands.$(date -Iminutes).txt


#-------------------------------------------------------------------------------
### 01. Prepare reference files
#-------------------------------------------------------------------------------

func_build_reference() {
    # Check if the reference genome file already exists
    if [ -f "${REFERENCE_FILES}/REF.fa" ]; then
        echo "Reference is already setup. Moving on."
        exit 0
    else
        # Copy the reference genome file to the REFERENCE_FILES directory
        cp "${REFERENCE}" "${REFERENCE_FILES}/REF.fa"
        # Create an index file for the reference genome
        samtools faidx "${REFERENCE_FILES}/REF.fa"
        # Create a dictionary file for the reference genome
        samtools dict "${REFERENCE_FILES}/REF.fa" > "${REFERENCE_FILES}/REF.dict"

        # Append the BAM_LIST file to the bam.list file in the REFERENCE_FILES directory
        cat "${BAM_LIST}" >> "${REFERENCE_FILES}/bam.list"

        # Split the reference genome file into chunks of approximately 10 Mb in size
        fastaq split_by_base_count "${REFERENCE_FILES}/REF.fa" "${REFERENCE_FILES}/REFsplit" 10000000

        # For each chunked genome section, create a list of contig/scaffold names
        for i in "${REFERENCE_FILES}/REFsplit"* ; do
            # Extract the name of the chunked genome section
            NAME=$( echo "${i}" | awk -F '/' '{print $NF}' )
            # Extract the contig/scaffold names from the chunked genome section
            grep ">" "${i}" | sed 's/>//g' > "${REFERENCE_FILES}/${NAME}.list"
        done
    fi
}

export -f func_build_reference


#-------------------------------------------------------------------------------
### 02. Make GVCF per sample
#-------------------------------------------------------------------------------
func_make_gvcf() {

# make jobs
COUNT=0
while read BAM; do
	n=1

    SAMPLE=$( echo ${BAM} | awk -F '/' '{print $NF}' | sed -e 's/.bam//g' )

    # check if the sample directory exists already - if yes, stop and move on
    if [ -d "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_complete" ] ; then
        echo -e "\nThere is already a run started/completed with this sample name. Rename and start again, or move on to the enxt sample\n"
        exit 0
    fi

    # make sample directories
    mkdir ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started
	mkdir ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES

     echo "gatk GatherVcfsCloud \\" > ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_gather_${SAMPLE}_gvcf

    for SEQUENCE in ${REFERENCE_FILES}/REFsplit*list; do
        SEQUENCE=$( echo ${SEQUENCE} | awk -F '/' '{print $NF}' )
	    echo -e "gatk HaplotypeCaller --input ${BAM} --output ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/${n}.${SAMPLE}.tmp.g.vcf.gz --reference ${REFERENCE_FILES}/REF.fa --intervals ${REFERENCE_FILES}/${SEQUENCE} --emit-ref-confidence GVCF " > ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_hc_${SAMPLE}.tmp.job_${n};
	    echo -e "--input ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/${n}.${SAMPLE}.tmp.g.vcf.gz \\" >> ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_gather_${SAMPLE}_gvcf;
	    let "n+=1";
        done;

	echo -e "--output ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/${SAMPLE}.g.vcf.gz; tabix -p vcf ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/${SAMPLE}.g.vcf.gz" >> ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_gather_${SAMPLE}_gvcf;

	echo -e "rm ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/*.tmp.* && mv ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_complete" > ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_clean_${SAMPLE};

	chmod a+x ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_*

	# setup job conditions
	JOBS=$( ls -1 ${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_hc_* | wc -l )
	ID="U$(date +%s)"

	#submit job array to call variants put scaffold / contig
	bsub -q long -R'span[hosts=1] select[mem>15000] rusage[mem=15000]' -n 6 -M15000 -J "gatk_make_gvcf_${ID}_[1-$JOBS]%100" -e "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES/gatk_make_gvcf_${ID}_[1-$JOBS].e" -o "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES/gatk_make_gvcf_${ID}_[1-$JOBS].o" "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_hc_${SAMPLE}.tmp.job_\$LSB_JOBINDEX"

	#submit job to gather gvcfs into a single, per sample gvcf
	bsub -q normal -w "done(gatk_make_gvcf_${ID}_[1-$JOBS])" -R'span[hosts=1] select[mem>500] rusage[mem=500]' -n 1 -M500 -J "gatk_gather_gvcf_${ID}" -e "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES/gatk_gather_gvcf_${ID}.e" -o "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES/gatk_gather_gvcf_${ID}s.o" "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_gather_${SAMPLE}_gvcf"

	# clean up
	bsub -q normal -w "done(gatk_gather_gvcf_${ID})" -R'span[hosts=1] select[mem>500] rusage[mem=500]' -n 1 -M500 -J "gatk_clean_gvcf_${ID}" -e "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES/gatk_clean_gvcf_${ID}.e" -o "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/LOGFILES/gatk_clean_gvcf_${ID}.o" "${GATK_HC_GVCFs}/${SAMPLE}_GATK_HC_GVCF_started/run_clean_${SAMPLE}"

	sleep 1
done < ${BAM_LIST}


# check that GVCF directories are complete before finishing this step
while true; do
  found=0
  count=$(cat $BAM_LIST | wc -l)
while read -r NAME; do
    if [ -d "${GATK_HC_GVCFs}/${NAME}_GATK_HC_GVCF_complete" ]; then
      found=$((found + 1))
    fi
    done  < <(cat "$BAM_LIST" | awk -F '/' '{print $NF}' | sed 's/.bam//g' )
  if [ ${found} -eq ${count} ]; then
    echo "Directories are complete. Moving on."
    break
  fi
  echo "Directories not found, waiting..."
  sleep 20
done

}

export -f func_make_gvcf



#-------------------------------------------------------------------------------
### 03. Merge GVCFs
#-------------------------------------------------------------------------------
func_merge_gvcf() {

ls -1 ${GATK_HC_GVCFs}/*complete/*gz > ${GATK_HC_MERGED}/gvcf.list

[ -d ${GATK_HC_MERGED}/LOGFILES ] || mkdir -p ${GATK_HC_MERGED}/LOGFILES


n=1
for SEQUENCE in ${REFERENCE_FILES}/REFsplit*list; do
    SEQUENCE=$( echo ${SEQUENCE} | awk -F '/' '{print $NF}' )
    echo -e "gatk CombineGVCFs -R ${REFERENCE_FILES}/REF.fa --intervals ${REFERENCE_FILES}/${SEQUENCE} \\" > ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_${n}
    while read SAMPLE; do
        echo -e "--variant ${SAMPLE} \\" >> ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_${n};
   done < ${GATK_HC_MERGED}/gvcf.list
   echo -e "--output ${GATK_HC_MERGED}/${n}.cohort.tmp.g.vcf.gz" >> ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_${n};
   let "n+=1";
done

chmod a+x ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_*

# setup job conditions
JOBS=$( ls -1 ${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_* | wc -l )
ID="U$(date +%s)"

#submit job array to call variants put scaffold / contig
bsub -q long -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -n 10 -M10000 -J "gatk_merge_gvcf_[1-$JOBS]%100" -e "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_[1-$JOBS].e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_[1-$JOBS].o" "${GATK_HC_MERGED}/run_merge_gvcfs.tmp.job_\$LSB_JOBINDEX"

rm ${GATK_HC_MERGED}/MERGE_ARRAY_FINISHED
bsub -w "done(gatk_merge_gvcf_)" -q normal -R'span[hosts=1] select[mem>100] rusage[mem=100]' -n 1 -M100 -J "gatk_merge_gvcf_finish" -e "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_finish.e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_merge_gvcf_finish.o" "touch ${GATK_HC_MERGED}/MERGE_ARRAY_FINISHED"

until [ -f "${GATK_HC_MERGED}/MERGE_ARRAY_FINISHED" ]
do
     sleep 10
done

}

export -f func_merge_gvcf





#-------------------------------------------------------------------------------
### 04. Genotype GVCFs
#-------------------------------------------------------------------------------

func_genotype_gvcfs() {

# split each chromosome up into separate jobs, and run genotyping on each individually.
n=1
for SEQUENCE in ${REFERENCE_FILES}/REFsplit*list; do
    SEQUENCE=$( echo ${SEQUENCE} | awk -F '/' '{print $NF}' )
    echo -e "gatk GenotypeGVCFs \
    -R ${REFERENCE_FILES}/REF.fa \
    -V ${GATK_HC_MERGED}/${n}.cohort.tmp.g.vcf.gz \
    --intervals ${REFERENCE_FILES}/${SEQUENCE} \
    -all-sites \
    -O ${GATK_HC_MERGED}/${n}.cohort.tmp.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation" > ${GATK_HC_MERGED}/run_hc_genotype.tmp.job_${n};
    let "n+=1";
done

chmod a+x ${GATK_HC_MERGED}/run_hc_genotype*

# setup job conditions
JOBS=$( ls -1 ${GATK_HC_MERGED}/run_hc_genotype* | wc -l )
ID="U$(date +%s)"

bsub -q long -R'span[hosts=1] select[mem>20000] rusage[mem=20000]' -n 6 -M20000 -J "gatk_genotype_cohort_gvcf_[1-$JOBS]" -e "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_[1-$JOBS].e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_[1-$JOBS].o" "${GATK_HC_MERGED}/run_hc_genotype.tmp.job_*\$LSB_JOBINDEX"

rm ${GATK_HC_MERGED}/GENOTYPE_ARRAY_FINISHED
bsub -w "done(gatk_genotype_cohort_gvcf_)" -q normal -R'span[hosts=1] select[mem>100] rusage[mem=100]' -n 1 -M100 -J "gatk_genotype_cohort_gvcf_finish" -e "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_finish.e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_finish.o" "touch ${GATK_HC_MERGED}/GENOTYPE_ARRAY_FINISHED"

until [ -f "${GATK_HC_MERGED}/GENOTYPE_ARRAY_FINISHED" ]
do
     sleep 10
done

}

export -f func_genotype_gvcfs



#-------------------------------------------------------------------------------
### 05. Finish making VCF and cleanup
#-------------------------------------------------------------------------------


func_finish_vcf() {

    #
    ls ${GATK_HC_MERGED}/*.cohort.tmp.vcf.gz > ${GATK_HC_MERGED}/cohort.vcf.list

    # concatenate the vcf files in the list
    vcf-concat --files ${GATK_HC_MERGED}/cohort.vcf.list > ${GATK_HC_MERGED}/${PREFIX}.cohort.$(date -I).vcf

    # Compress the combined VCF file with bgzip
    bgzip -f ${GATK_HC_MERGED}/${PREFIX}.cohort.$(date -I).vcf

    # Create a tabix index for the compressed combined VCF file
    tabix -f ${GATK_HC_MERGED}/${PREFIX}.cohort.$(date -I).vcf.gz

    # Remove all files in the directory specified by GATK_HC_MERGED that match the pattern *tmp*
    rm ${GATK_HC_MERGED}/*tmp*

}

export -f func_finish_vcf





#-------------------------------------------------------------------------------
# running the pipeline
#-------------------------------------------------------------------------------

# func_build_reference
bsub -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>1000] rusage[mem=1000]" -M1000 -o ${LOG_FILES}/gatk_01_build_reference.o -e ${LOG_FILES}/gatk_01_build_reference.e -J gatk_01_build_reference_${PREFIX} func_build_reference

# func_make_gvcf
bsub -w "done(gatk_01_build_reference_${PREFIX})" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n20 -o ${LOG_FILES}/gatk_02_make_gvcf.o -e ${LOG_FILES}/gatk_02_make_gvcf.e -J gatk_02_make_gvcf_${PREFIX} func_make_gvcf

# func_merge_gvcf
bsub -w "done(gatk_02_make_gvcf_${PREFIX})" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n20 -o ${LOG_FILES}/gatk_03_merge_gvcf.o -e ${LOG_FILES}/gatk_03_merge_gvcf.e -J gatk_03_merge_gvcf_${PREFIX} func_merge_gvcf

# func_genotype_gvcfs
bsub -w "done(gatk_03_merge_gvcf_${PREFIX})" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n20 -o ${LOG_FILES}/gatk_04_genotype_gvcfs.o -e ${LOG_FILES}/gatk_04_genotype_gvcfs.e -J gatk_04_genotype_gvcfs_${PREFIX} func_genotype_gvcfs

# func_finish_vcf
bsub -w "done(gatk_04_genotype_gvcfs_${PREFIX})" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>1000] rusage[mem=1000]" -q long -M1000 -n1 -o ${LOG_FILES}/gatk_05_finish_vcf.o -e ${LOG_FILES}/gatk_05_finish_vcf.e -J gatk_05_finish_vcf_${PREFIX} func_finish_vcf

 ```