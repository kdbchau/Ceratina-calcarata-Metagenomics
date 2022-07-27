#!/bin/bash
#SBATCH --account=def-sanrehan
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=0G
#SBATCH --time=24:00:00

module load StdEnv/2020
module load spades/3.10.1

# Run metaspades.py with different kmers for the metagenomic Ceratina dataset
INPUTS=$(ls /scratch/chauk/4-Proper/0-unmapped_fastq/*R1.fastq)

for files in ${INPUTS}
do
        PREFIX="${files%%_unmapped*}"
        R2="${PREFIX}_unmapped_R2.fastq"
        PREFIX="$(basename ${PREFIX})"
        metaspades.py -k 21,33,55,77 -t 32 -1 $files -2 $R2 -o ${PREFIX}_metagenome_results/
done

## the contig files are found in the metagenome_results folders
