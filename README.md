# Integrative metagenomics of wild bees reveals urbanization tied to reduced genetic diversity and increased pathogen loads

This repository hosts the scripts and brief methods from our manuscript (title above). This project focused on using mapped and unmapped reads from whole genome sequencing of the small carpenter bee (_Ceratina calcarata_) to obtain population genetic and metagenomic data, respectively.

Authors:
* [Katherine D. Chau (me)](https://www.linkedin.com/in/balasink/)
* [Farida Samad-zada](https://www.linkedin.com/in/faridasamadzada/?originalSubdomain=ca)
* [Evan P. Kelemen](https://scholar.google.com/citations?user=eg7ziI8AAAAJ&hl=en)
* [Sandra M. Rehan](http://www.rehanlab.com/people.html)

Acknowledgements:
* [Rehan Lab](http://www.rehanlab.com/) for constructive feedback on earlier versions of our manuscript
* Alessia Schembri, Phuong Nguyen, and Yousaf Malik for assistance with field collection of bee specimens.
* [SHARCNET](https://www.sharcnet.ca/my/front/) community for technical support on high-performance computing (HPC) clusters ([Graham](https://docs.alliancecan.ca/wiki/Graham), [Beluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga), and [Cedar](https://docs.alliancecan.ca/wiki/Cedar))
    * Jose Sergio Hleap for assistance with the HPC clusters and overall pipeline performance.

Funding:
* [Mitacs Elevate Postdoctoral Fellowship](https://www.mitacs.ca/en/programs/elevate) to KDC
* [NSF Biological Collections Postdoctoral Fellowship](https://www.nsf.gov/) to EPK
* [Foundation for Food and Agricultural Research Pollinator Health](https://foundationfar.org/programs/pollinator-health-fund/) to SMR
* [Weston Family Foundation Microbiome Initiative](https://westonfoundation.ca/weston-family-microbiome-initiative/) to SMR
* [NSERC Discovery Grant](https://www.nserc-crsng.gc.ca/professors-professeurs/grants-subs/dgigp-psigp_eng.asp) to SMR
* [Supplement and E.W.R Steacie Memorial Fellowship](https://www.nserc-crsng.gc.ca/prizes-prix/steacie-steacie/about-apropos_eng.asp) to SMR.

# Table of Contents
1. [Install software](#1-install-software)
2. [Quality Check](#2-quality-check)
3. [Separating Mapped and Unmapped Reads](#3-separating-mapped-and-unmapped-reads)
4. [Population Genetics](#4-population-genetics)
    *
5. [Metagenomics](#5-metagenomics)  
    * [metaSPADES](#51-metaspades)
    * [Diversity Statistics](#52-diversity-statistics)
    * [Random Forests](#53-random-forests)
    * [DESeq2 and WGCNA](#54-deseq2-and-wgcna)
    * [Functional Analyses](#55-functional-analyses)
        * [GhostKOALA](#551-ghostkoala)
        * [eggNOG-mapper](#552-eggnog-mapper)
6. [CoNet - Integration of Pop Gen and Metagenomics](#6-conet-integration-of-pop-gen-and-metagenomics)

# 1. Install software
We will be using several software in order to take our unmapped whole genome sequence (WGS) reads and extract taxonomic information from them.

1. [ __FastQC__](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (version 0.11.9)
2. [__Trimmomatic__](http://www.usadellab.org/cms/?page=trimmomatic) (version 0.39)
3. [__BWA__](https://github.com/lh3/bwa) (version 0.7.17)
4. [__SAMtools__](http://www.htslib.org/) (version 1.10)
5. [__Sambamba__](https://lomereiter.github.io/sambamba/) (version 0.8.0)
6. [__BAMtools__](https://github.com/pezmaster31/bamtools) (version 2.5.1)
7. [__Bedtools__](https://bedtools.readthedocs.io/en/latest/) (version 2.29.2)
8. [__Metaspades__](https://cab.spbu.ru/software/meta-spades/) (version 3.10.1)
9. [__Kraken2__](https://ccb.jhu.edu/software/kraken2/) (version 2.1.1)
10. [__Diamond__](https://github.com/bbuchfink/diamond) (version 2.0.13)
11. [__Blast+__](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (version 2.12.0)
12. [__Entrez Direct EUtilities (Efetch)__](https://www.ncbi.nlm.nih.gov/books/NBK179288/) (version 15.3)
13. [__seqtk__](https://github.com/lh3/seqtk) (version 1.3-r106)
14. [__FragGeneScan__](https://sourceforge.net/projects/fraggenescan/) (version 1.31)

Note: The above software are already in the clusters and just need to be loaded using "module load". For example:

```
   module load fastqc # will load the latest version in the cluster
   module load samtools/1.10 # specify the version to use a specific one
   module load bwa/0.7.17
```

Loading these modules as such will appear in some of our scripts.


# 2. Quality Check
Run FASTQC on the raw reads to check their quality, followed by trimming the adaptors and poor quality bases. After the trim, run FASTQC again to ensure adaptors and base quality is acceptable.
[__FastQC script__](https://github.com/kdbchau/Ceratina-calcarata-Metagenomics/blob/main/Scripts/fastqc.sh)

Adaptors are there and there are some poor quality bases. Will run trimmomatic to improve quality.
[__Trimmomatic script__](https://github.com/kdbchau/Ceratina-calcarata-Metagenomics/blob/main/Scripts/trimscript.sh)

After this trim, run FASTQC again to check all is fine. Then move on to the read mapping! Can even use [__MultiQC__](https://multiqc.info/) here to get a nice comprehensive report of the trimmed results.

# 3. Separating Mapped and Unmapped Reads

The following script uses BWA to map reads to a reference genome, and incorporates sambamba, samtools, and bedtools for extraction of mapped and unmapped reads.
[__Mapping script__](https://github.com/kdbchau/Ceratina-calcarata-Metagenomics/blob/main/Scripts/mapping.sh).

We used the **_Ceratina calcarata_ (BioProject: PRJNA791561)** genome to map the reads. The data will be held until December 31, 2023.

This will separate mapped and unmapped reads into fastq files which can then be used for population genetics [section 4](#4-population-genetics) or metagenomics [section 5](#5-metagenomics), respectively.

# 4. Population Genetics

## 4.1.

# 5. Metagenomics

## 5.1. metaSPADES

[__Metaspades__](https://cab.spbu.ru/software/meta-spades/) (version 3.10.1) is another option to use for gathering information on the sequences. In this case, contigs are produced which can then be BLASTed. [__Metaspades script__](https://github.com/kdbchau/Ceratina-calcarata-Metagenomics/blob/main/Scripts/metaspades.sh).

Once the ```contigs.fasta``` files are obtained, we can BLAST the contigs to obtain taxonomy information. To do this, we use BLAST+ and the [__BLAST+ script__](https://github.com/kdbchau/Ceratina-calcarata-Metagenomics/blob/main/Scripts/blast.sh).

## 5.2. Diversity Statistics

Here are a list of mini R scripts or R code that I used to run diversity stats for my data.

PERMANOVA via _adonis2_ using Bray-Curtis method, followed by ANOVA for beta dispersion via _betadisper_. This is then followed by TUKEY HSD via _TukeyHSD_ if PERMANOVA was significant and beta dispersion was not significant. [Diversity R Script]().
