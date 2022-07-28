#!/bin/bash
#SBATCH --account=def-sanrehan 
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=48:00:00        # time
#SBATCH --job-name=gatkHpCaller
#SBATCH --output=slurm-%x_%j.out

export OMP_NUM_THREADS=$SLURM_NTASKS

# 1. GATK Variant Calling 

echo "Task started for Set ${SLURM_ARRAY_TASK_ID}"
#Parallelization setup similar to the one in mapping script 
# We will be using GNU parallel, and we need to record the environment
parallel --record-env # this has to go just after the sbatch directives, it copies all user-defined environmental variables to the remote system

STARTTIME=$(date +%s)
date_stamp=$(date +"%F %H:%M:%S") && echo "Job started on $date_stamp"
echo "START-------------------------"

## Packages to load

module load gatk
module load samtools
module load vcftools
module load bcftools

## One way to increase efficiency of your runs is to "stage" or copy the IO
## heavy files into localscratch (within a jobscript can be accessed with $SLURM_TMPDIR)
REFPATH=/home/sfarida/projects/def-sanrehan/sfarida/reference_fasta
parallel --will-cite -j ${SLURM_NTASKS} cp {} ${SLURM_TMPDIR} ::: ${REFPATH}/jasmine-uni1041-mb-hirise-teril_10-15-2020__final_assembly.fasta
REF=${SLURM_TMPDIR}/jasmine-uni1041-mb-hirise-teril_10-15-2020__final_assembly.fasta #using new reference genome file 

#bam files have been mapped to the reference genome, reads have been assigned, and duplicates were marker prior to variant calling 

## Staging the inputs might be good as well
STAGED_PATH=${SLURM_TMPDIR}/inputs
TO_STAGE=/home/sfarida/scratch/Mapping_2/Mapped_bam_files/MappedRGMD_bam_files_for_analysis/Subset${SLURM_ARRAY_TASK_ID}
mkdir -p ${STAGED_PATH}
parallel --will-cite -j ${SLURM_NTASKS} cp {} ${STAGED_PATH} ::: ${TO_STAGE}/*.ba*
INPUTS=$(ls ${STAGED_PATH}/*.bam)

# Setup Outputs
OUTPUTDIR=${SLURM_TMPDIR}/output/GATK_output${SLURM_ARRAY_TASK_ID}
mkdir -p ${OUTPUTDIR}
cd ${OUTPUTDIR}

## INDEX REFERENCE FILES and create Dictionary (This only has to be done once and then this step isn't used again for all samples)
if [ -s ${REF}.fai ]
then
        echo "1. $(basename ${REF}) already indexed, skipping"
else
        echo "1. Creating index and dictionary for ${REF}"
        samtools faidx ${REF}
        gatk CreateSequenceDictionary -R ${REF}
fi

for bamfile in ${INPUTS}

do
        PREF="${bamfile/_mapped_sorted*}"
        PREF=$(basename ${PREF})
        if [ -s "${PREF}.g.vcf.gz" ]
        then
                echo "GVCF already exists for $(basename ${PREF}.g.vcf.gz), skipping."
        else
                echo "2. Processing ${PREF} to ${OUTPUTDIR}"
                ( time bash -c \
                "gatk HaplotypeCaller \
                -R $REF \
                -I "$bamfile" \
                 -ERC GVCF \
                --native-pair-hmm-threads ${SLURM_NTASKS} \
                -O "${PREF}".g.vcf.gz")
        fi

        # COPY RESULTS OF THIS STEP TO PROJECT (IN CASE IT TIMES OUT)
        cp ${PREF}.g.vcf.gz* ${TO_STAGE}
done

echo "FINISH--------------------------"
ENDTIME=$(date +%s)
date_stamp=$(date +"%F %H:%M:%S") && echo "Job finished on $date_stamp"
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"


rsync -raz --progress * ${TO_STAGE}
rsync -raz --progress ${STAGED_PATH}/* ${TO_STAGE}
echo "All transfers done"

# --- once the above script finishes running, veryfy that all files were transferred successfully, and proceed to submit the next chunk as a separate SLURM job ---

### Creating GATK database

        gatk GenomicsDBImport \
        --genomicsdb-workspace-path ~/scratch/VariantCallingPipeline/GATK_pipeline/DB/DB_workspace_mergedINT \
        --sample-name-map GVCF_map.txt \ #this is a file listing all sample IDs and their corresponding GVCF files 
        --batch-size 50 \
        --L scaffolds.list \ #this is a file listing all the scaffolds 
        --merge-contigs-into-num-partitions 25

# --- once the above script finishes running, veryfy that all files were transferred successfully, and proceed to submit the next chunk as a separate SLURM job ---

### Genotyping individual GVCFs

REF=/home/sfarida/projects/def-sanrehan/sfarida/reference_fasta/jasmine-uni1041-mb-hirise-teril_10-15-2020__final_assembly.fasta

if [ -s ${REF}.fai ]
then
        echo "1. $(basename ${REF}) already indexed, skipping"
else
        echo "1. Creating index and dictionary for ${REF}"
        samtools faidx ${REF}
fi


if [-s ${REF}.dict]
then
        echo "(basename ${REF}) already has associated dictionary, skiping"
else
        echo "Creating dicttionary for ${REF}"
        gatk CreateSequenceDictionary -R ${REF}
fi

        gatk GenotypeGVCFs \
        -R $REF \
        -V gendb:///home/sfarida/scratch/VariantCallingPipeline/GATK_pipeline/DB/DB_workspace_mergedINT \
        -O gatk_output.vcf.gz

### Variant Calling complete 
### gatk_output.vcf.gz is used for downstream analysis

# --- once the above script finishes running, veryfy that all files were transferred successfully, and proceed to submit the next chunk as a separate SLURM job ---


bcftools filter --SnpGap 10 ../gatk_output.vcf.gz -O z -o gatk_output_gap10.vcf.gz

bcftools index gatk_output_gap10.vcf.gz -t

bcftools view gatk_output_gap10.vcf.gz --types snps -O z -o gatk_SNPs.vcf.gz

bcftools index gatk_SNPs.vcf.gz -t

gatk VariantFiltration \
                -V gatk_SNPs.vcf.gz \
                -filter "QD < 2.0" --filter-name "QD2" \
                -filter "QUAL < 30.0" --filter-name "QUAL30" \
                -filter "SOR > 3.0" --filter-name "SOR3" \
                -filter "FS > 60.0" --filter-name "FS60" \
                -filter "MQ < 40.0" --filter-name "MQ40" \
                -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
                -O gatk_filtered_snps.vcf.gz


bcftools view -f 'PASS,.' gatk_filtered_snps.vcf.gz -O z -o gatk_snps_filterpass.vcf.gz
bcftools stats gatk_snps_filterpass.vcf.gz > ./stats/snps_after_hard_filtering.txt

vcftools --gzvcf gatk_snps_filterpass.vcf.gz --min-meanDP 5  --max-meanDP 100 --mac 8 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out  gatk_snps_fpass_mf1

#Remove individuals with less than 3X coverage

vcftools --vcf gatk_snps_fpass_mf1.recode.vcf --remove samples_to_remove.txt --minDP 3 --min-meanDP 8 --max-missing 0.9 --recode --recode-INFO-all --out gatk_snps_fpass_mf1_178indv_minDP3_minmeanDP8_maxMissing90

bgzip gatk_snps_fpass_mf1_178indv_minDP3_minmeanDP8_maxMissing90.recode.vcf

# --- once the above script finishes running, veryfy that all files were transferred successfully, and proceed to submit the next chunk as a separate SLURM job ---

# 2. BCFTools mpileup Variant Calling 

bcftools mpileup -Ou -f $REF -b list.txt -d 8000 -q 40 -Q 20 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
        | bcftools call --threads 4 -f GQ -m -v -O v -o bcf_mp_call_180_jasMap.vcf

mkdir -p stats


bcfools filter --SnpGap 10 ../bcf_mp_call_180_jasMap.vcf -O z -o mp_variants_gap10.vcf.gz


#echo "Selecting SNPS -----------------------------"
        bcftools view  mp_variants_gap10.vcf.gz --types snps -O z -o mp_SNPs.vcf.gz
        gatk VariantFiltration \
                -V mp_SNPs.vcf.gz \
                -filter "QUAL < 30.0" --filter-name "QUAL30" \
                -filter "MQ < 40.0" --filter-name "MQ40" \
                -O mp_filtered_snps.vcf.gz


bcftools view -f 'PASS,.' mp_filtered_snps.vcf.gz -O z -o mp_snps_filterpass.vcf.gz

bcftools stats mp_snps_filterpass.vcf.gz > ./stats/SNPs_after_hardfiltering.txt

mkdir -p  manual_filtering

cd manual_filtering

mkdir -p stats

vcftools --gzvcf ../mp_snps_filterpass.vcf.gz --min-meanDP 5  --max-meanDP 100 --mac 8 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out  mp_snps_filterpass_manualfilter1
 
bgzip mp_snps_filterpass_manualfilter1.recode.vcf

bcftools index mp_snps_filterpass_manualfilter1.recode.vcf.gz -t

bcftools stats mp_snps_filterpass_manualfilter1.recode.vcf.gz  > ./stats/snps_after_manual_filter1.txt

#Removing low coverage individuals <3X 

vcftools --gzvcf mp_snps_filterpass_manualfilter1.recode.vcf.gz --remove samples_to_remove.txt  --max-missing 0.8 --recode --recode-INFO-all --out mp_snps_fpass_mf1_178indv_filter2

bgzip mp_snps_fpass_mf1_178indv_filter2.recode.vcf


# --- once the above script finishes running, veryfy that all files were transferred successfully, and proceed to submit the next chunk as a separate SLURM job ---


# 3. Intersecting GATK and BCFTools VCF files

VCFgatk=./gatk_snps_fpass_mf1_178indv_minDP3_minmeanDP8_maxMissing90.recode.vcf.gz
VCFmp=./mp_snps_fpass_mf1_178indv_filter2.recode.vcf.gz


module load bcftools
module load vcftools

bcftools isec -p SNPs_in_common $VCFgatk $VCFmp

#Intersected GATK file was used for downstream analysis 








