#!/bin/bash
#SBATCH --account=def-sanrehan
#SBATCH --nodes=1
#SBATCH --mem=16G 
#SBATCH --time=10:00:00

module load fastqc

# I saved the names of the raw fastq files in a textfile as a list
# so in case some files did not complete I could just edit the textfile 
# and work on those files instead of accidentally rerunning fastqc on all files in the directory

# NOTE: the output directory fastqc_raw was created prior to executing this script! Must be done or fastqc won't work

while read line;
do
  echo "Starting fastqc for $line"
  fastqc --kmers 7 --outdir fastqc_raw $line
done < fastq_files.txt
