#!/bin/bash
#SBATCH --account=def-sanrehan #Your PI's Account
#SBATCH --cpus-per-task=48     # Ask for full nodes is better (32 graham and Beluga / 48 in Cedar)
#SBATCH --nodes=1              # Make sure is only one node (nothing is MPI)
#SBATCH --mem=0                # Reserve all memory in that node
#SBATCH --time=20:00:00        # time
#SBATCH --array=1-12           # Generate the array. This is because I personally split my datafiles into groups of 30 (15 sets of paired end fastq files) for easier management.

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK # set number of threads which would be 48

# We will be using GNU parallel, and we need to record the environment
parallel --record-env # this has to go just after the sbatch directives, it copies all user-defined environmental variables to the remote system

## Packages to load for bamtools/2.5.1
module load StdEnv/2020

## Main modules used for analysis
module load samtools/1.10
module load bwa/0.7.17
module load bamtools/2.5.1
module load bedtools/2.29.2
module load sambamba/0.8.0

## One way to increase efficiency of your runs is to "stage" or copy the IO 
## heavy files into localscratch (within a jobscript can be accessed with $SLURM_TMPDIR)
REFPATH=/scratch/chauk/project1_2021/annotation_12.30.2020
parallel --will-cite -j ${SLURM_CPUS_PER_TASK} cp {} ${SLURM_TMPDIR} ::: ${REFPATH}/jasmine-uni1041-mb-hirise-teril_10-15-2020__final_assembly.fasta
REF=${SLURM_TMPDIR}/jasmine-uni1041-mb-hirise-teril_10-15-2020__final_assembly.fasta

## Staging the inputs might be good as well
STAGED_PATH=${SLURM_TMPDIR}/inputs
TO_STAGE=/home/chauk/scratch/project1_2021/Ceratina${SLURM_ARRAY_TASK_ID}/ceratina${SLURM_ARRAY_TASK_ID}_paired
mkdir -p ${STAGED_PATH}
parallel --will-cite -j ${SLURM_CPUS_PER_TASK} cp {} ${STAGED_PATH} ::: ${TO_STAGE}/*.fastq.gz
INPUTS=$(ls ${STAGED_PATH}/*_R1_trimmed_paired.fastq.gz)

## Setup Outputs
OUTPUTDIR=${SLURM_TMPDIR}/output/Ceratina${SLURM_ARRAY_TASK_ID}
mkdir -p ${OUTPUTDIR}
cd ${OUTPUTDIR}

## INDEX REFERENCE FILES (This only has to be done once and then this step isn't used again for all samples)
if [ -s ${REF}.bwt ]
then
	echo "1. $(basename ${REF}) already indexed, skipping"
else
	echo "1. Creating index for ${REF}"
	bwa index ${REF}
fi

## BWA MEM IN TRIMMED PAIRED FOLDER TO SORTED FILES
## SINCE THIS RUNS ON FOLDER WITH 15 PAIRS IS MORE EFFICIENT TO USE THREADS OF THE PROGRAMS
for R1 in ${INPUTS}
do
	PREF="${R1/_R1*}"
	R2="${PREF}_R2_trimmed_paired.fastq.gz"
	PREF=$(basename ${PREF})
	# MAP, CONVERT, SORT, AND INDEX
	if [ -s "${STAGED_PATH}/${PREF}_aln_sorted.bam" ]
	then
		echo "2. SAM, BAM or SORTED BAM file already exists for $(basename ${R1/_R1*}_aln.sam), skipping."
	else
		echo "2. Processing ${PREF} to ${OUTPUTDIR}"
		( time bash -c \
		"bwa mem -t ${SLURM_CPUS_PER_TASK} "${REF}" "${R1}" "${R2}" -M -B 2 | \
		sambamba view -t ${SLURM_CPUS_PER_TASK} -S -f bam /dev/stdin | \
		sambamba sort -t ${SLURM_CPUS_PER_TASK} -o "${PREF}_aln_sorted.bam" /dev/stdin | \
		sambamba index -t ${SLURM_CPUS_PER_TASK} "${PREF}_aln_sorted.bam"" ) 2> ${PREF}_STEPS2to4.log
		bamfile="${PREF}_aln_sorted.bam"
		# COPY RESULTS OF THIS STEP TO PROJECT (IN CASE IT TIMES OUT)
		cp ${PREF}_aln_sorted* ${TO_STAGE}
	fi
	
  
	# GET STATS
	if [ -s "${STAGED_PATH}/${PREF}.[sb]amtools.stats.txt" ]
	then
		echo "5. Mapping stats for ${PREF} already done, skipping."  # If STATS files already exist, skip.
	else
		echo "5. Gathering mapping stats for $(basename ${PREF})"
		( time bamtools stats -in ${bamfile} -insert > ${PREF}.bamtools.stats.txt ) 2> ${PREF}_stats.log # THIS ONE RUNS SEQUENTIALLY DO YOU REALLY NEED IT?
		( time bash -c \
		"sambamba flagstat -t ${SLURM_CPUS_PER_TASK} "${bamfile}" > "${PREF}.samtools.stats.txt"" ) 2>> ${PREF}_stats.log
		# COPY RESULTS OF THIS STEP TO PROJECT (IN CASE IT TIMES OUT)
		cp ${PREF}*.stats.txt ${TO_STAGE}
	fi
	
  
	# GET UNMAPPED and MAPPED READS
	if [ -s "${STAGED_PATH}/${PREF}_*mapped_sorted.bam" ] 
	then
		echo "${PREF} SORTED MAPPED and SORTED UNMAPPED BAM files exist, skipping"
	else
		echo "Extracting UNMAPPED and MAPPED reads from $(basename ${bamfile})"
		# YOU CAN ADD MORE FILTERS... CHECK https://github.com/biod/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax
		( time bash -c \
		"sambamba view -f bam -F 'unmapped or mate_is_unmapped and not secondary_alignment' -t ${SLURM_CPUS_PER_TASK} ${bamfile} | \
		sambamba sort -t ${SLURM_CPUS_PER_TASK} -o ${PREF}_unmapped_sorted.bam /dev/stdin" ) 2> ${PREF}_unmapped.log
		( time bash -c \
		"sambamba view -f bam -F 'not (unmapped or mate_is_unmapped) and not secondary_alignment' -t ${SLURM_CPUS_PER_TASK} ${bamfile} | \
		sambamba sort -t ${SLURM_CPUS_PER_TASK} -o ${PREF}_mapped_sorted.bam /dev/stdin" ) 2> ${PREF}_mapped.log
		# COPY RESULTS OF THIS STEP TO PROJECT (IN CASE IT TIMES OUT)
		cp ${PREF}*mapped_sorted.bam ${TO_STAGE}
		bamfile_mapped=${PREF}_mapped_sorted.bam
		bamfile_unmapped=${PREF}_unmapped_sorted.bam
	fi
	
	
	# CONVERT UNMAPPED BAM TO FASTQ
	if [ -s ${STAGED_PATH}/${PREF}.fastq ]
	then
		echo "FASTQ for ${bamfile_unmapped} reads already exists, skipping."
	else
		echo "Extracting FASTQ for the reads in $(basename ${bamfile})."
		bamToFastq -i ${bamfile_unmapped} -fq ${PREFIX}_unmapped.R1.fastq -fq2 ${PREFIX}_unmapped.R2.fastq  
		cp ${PREF}_unmapped.R[12].fastq ${TO_STAGE}
	fi
	
	# CONVERT MAPPED BAM TO FASTQ/FASTA
	if [ -s ${STAGED_PATH}/${PREF}.fastq ]
	then
		echo "FASTQ for ${bamfile_mapped} reads already exists, skipping."
	else
		echo "Extracting FASTQ for the reads in $(basename ${bamfile})."
		bamToFastq -i ${bamfile_mapped} -fq ${PREFIX}_mapped.R1.fastq -fq2 ${PREFIX}_mapped.R2.fastq  
		cp ${PREF}_mapped.R[12].fastq ${TO_STAGE}
	fi 
done

echo "Pipeline Done"
# MAKE SURE EVERYTHING HAS BEEN MOVED, BUT DO NOT COPY THINGS ALREADY THERE
rsync -raz --progress * ${TO_STAGE}
rsync -raz --progress ${STAGED_PATH}/* ${TO_STAGE}
echo "All transfers done"
