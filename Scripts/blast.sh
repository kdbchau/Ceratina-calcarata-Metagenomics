#!/bin/bash
#SBATCH --account=def-sanrehan
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --array=0-180  # I had 180 samples to go through

module load StdEnv/2020
module load gcc/9.3.0
module load blast+/2.11.0

# Make access to the blast database faster through slurm tmpdir
parallel --will-cite --j "${SLURM_CPUS_PER_TASK}" rsync -lrtv {} "${SLURM_TMPDIR}" ::: /datashare/BLASTDB/nt*  # Load the nt database onto slurm tmpdir

# make access to input files faster by passing them through slurm tmpdir
STAGED_PATH=${SLURM_TMPDIR}/inputs
# declare lets you change the variable attribute without having to reassign it each time you loop through them.
# this is esp useful when the variable (i.e. filename) is identical in each instance.
# For example out of my 180 metagenome_results folders, the contig file is always named "contigs.fasta", so the name is the same but the actual file is different
declare -a ALL_INPUTS=( /scratch/chauk/4-Proper/3.10-metaspades/*_metagenome_results/contigs.fasta )
TO_STAGE=${ALL_INPUTS[${SLURM_ARRAY_TASK_ID}]}
mkdir -p "${STAGED_PATH}"

cp "${TO_STAGE}" "${STAGED_PATH}"  # copy over the contigs.fasta files to staged path

FILE="${STAGED_PATH}"/contigs.fasta
PREFIX=$(basename "$(dirname "${TO_STAGE}")")
PREFIX="${PREFIX%%_metagenome_results}"
OUTPUT="${SLURM_TMPDIR}/${PREFIX}_blastn_singlehit_output.txt"

blastn -db "${SLURM_TMPDIR}"/nt -num_threads 32 -max_target_seqs 1 -max_hsps 1 \  # -max_target_seqs 1 and -max_hsps 1 to get just the top hit for each read
        -outfmt "6 qseqid sseqid qlen slen pident evalue score staxids stitle" \
        -query "${FILE}" > "${OUTPUT}" && wait && echo "DONE" >> "${OUTPUT}"
wait

cp "${OUTPUT}" "${PWD}"
