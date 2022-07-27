#!/bin/bash
#SBATCH --account=def-sanrehan
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=16G
#SBATCH --time=03:00:00

module load java
module load StdEnv/2020
module load trimmomatic

# Save the forward and reverse reads as separate arrays
for_array=(*R1.fastq.gz) #all forward R1 files in current directory
rev_array=(*R2.fastq.gz) #all reverse R2 files in current directory

# Call each forward and reverse (R1 and R2) read for the same sample. Run trimmomatic
for ((i=0; i<${#for_array[@]}; i++));
do
        echo "RUNNING R1 ${for_array[$i]} and R2 ${rev_array[$i]}" >> trim_0.39_output.log
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 16 -phred33\
        ${for_array[$i]} ${rev_array[$i]}\
        ${for_array[$i]%%.fastq.gz}_trimmed_paired.fastq.gz\
        ${for_array[$i]%%.fastq.gz}_trimmed_unpaired.fastq.gz\
        ${rev_array[$i]%%.fastq.gz}_trimmed_paired.fastq.gz\
        ${rev_array[$i]%%.fastq.gz}_trimmed_unpaired.fastq.gz\
        ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq2-PE.fa:2:30:10 \
        SLIDINGWINDOW:4:20 \
        HEADCROP:5 \
        MINLEN:36 2>> trim_0.39_output.log
done
