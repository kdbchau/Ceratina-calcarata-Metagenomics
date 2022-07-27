# Integrative metagenomics of wild bees reveals urbanization tied to reduced genetic diversity and increased pathogen loads

This repository hosts the scripts and brief methods from our [__manuscript__](). This project focused on using mapped and unmapped reads from whole genome sequencing of the small carpenter bee (_Ceratina calcarata_) to obtain population genetic and metagenomic data, respectively. Our main goal was to understand how the genetic and metagenomic components of the bee are influenced by an urbanization gradient.

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
15. [__Past__](https://past.en.lo4d.com/windows) (version 4.06)
16. [__GhostKOALA__](https://www.kegg.jp/ghostkoala/)
17. [__eggNOG-mapper__](https://github.com/eggnogdb/eggnog-mapper) (version 2.1.6)
18. [__R__](https://www.r-project.org/) (version 4.2.0)

Note: Majority of the above software are already in the Compute Canada clusters and just need to be loaded using "module load". For example:

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

PERMANOVA via _adonis2_ using Bray-Curtis method, followed by ANOVA for beta dispersion via _betadisper_. This is then followed by TUKEY HSD via _TukeyHSD_ if PERMANOVA was significant and beta dispersion was not significant. [__Diversity R Script__](https://github.com/kdbchau/Ceratina-calcarata-Metagenomics/blob/main/Scripts/diversity.R).

SIMPER analyses was done using [__Past__](https://past.en.lo4d.com/windows) and simply uploading the dataframe which includes taxa relative contig abundances and environmental feature bin information. The dataframe should be set up where each sample is the row and each variable (landscape feature and family/genus) are the columns. Then simply selecting SIMPER analysis while highlighting all the taxa and all the landscape features, selecting pool data, and using Bray-Curtis dissimilarities.

## 5.3. Random Forests

To test if random forests would be worthwhile on this data, I used the following [__random forest R script__](https://github.com/kdbchau/Ceratina-calcarata-Metagenomics/blob/main/Scripts/randomforest.R) to simply assess whether classification or regression random forest analyses would produce low out-of-box errors or high variance explained %.

## 5.4. DESeq2 and WGCNA

The taxa dataframe of relative contig abundances is first normalized using the [__DESeq2 R script__](https://github.com/kdbchau/Ceratina-calcarata-Metagenomics/blob/main/Scripts/deseq2.R), and then the obtained normalized dataframe is used for WGCNA analysis using the [__WGCNA R script__](https://github.com/kdbchau/Ceratina-calcarata-Metagenomics/blob/main/Scripts/wgcna.R).

DESeq2 will identify in which category are taxa overrepresented (positive, significant Log2foldC). If the Log2FoldC is negative, then that taxa is more abundant in the second, contrasted sample.

WGCNA aims to find significant modules or clusters of taxa with respect to each landscape variable, and from those clusters we can identify hub genes by highlighting those that have an absolute gene significance (GS) value ≥ 0.2 and an absolute module membership (MM) value ≥ 0.8 (this can be done in Excel).

## 5.5. Functional Analyses

For the contigs that correspond to our domains of interest (these sequences were isolated from the metaspades ```contigs.fasta``` output files), we run them through **FragGeneScan** which finds fragmented genes in short reads. 


```
# Run FragGeneScan on the cluster like this, or in an sbatch script
for file in *extracted_contigs.fasta; do name=${file%%_*}; /scratch/chauk/06-EXTRACTED_CONTIGS/FragGeneScan1.31/FragGeneScan -s /scratch/chauk/06-EXTRACTED_CONTIGS/${name}_extracted_contigs.fasta -o /scratch/chauk/06-EXTRACTED_CONTIGS/${name}_output -p 1 -w 0 -t illumina_5; echo "done for $file"; done
# -w 0 means fasta file has short sequence reads (-w 1 would mean complete genomic sequence)
# -p 1 means use 1 thread (default is already 1, can omit this)
# -t illumina_5 means train the data using Illumina sequencing reads with about 0.5% error rate.
```

This will produce three types of output files but we will be using the ```output.faa``` files for further functional analysis.


### 5.5.1. GhostKOALA

The ```output.faa``` files were concantenated by bin for each of the variables and run through [__GhostKoala__](https://www.kegg.jp/ghostkoala/). So for example, for agricultural percent I had concantenated all the ```output.faa``` samples that belonged to bin 1, then to bin 2, then bin 3, etc. I ended up with 5 ```.faa``` files. I ran these through GhostKoala to get the KEGG reconstructed pathways list.


### 5.5.2. eggNOG-mapper

eggNOG-mapper will use Diamond (the nr database) to blast the fragmented sequences to obtain protein information. The ```output.faa``` files are used. eggNOG-mapper needs to be manually installed into the cluster. After downloading the ```tar.gz``` file and it is uploaded into the cluster, run the following:


```
module load python
tar xvzf 2.1.6.tar.gz
rm 2.1.6.tar.gz
cd eggnog-mapper-2.1.6
pip install -r requirements.txt
python download_eggnog_data.py  ## note have python module loaded before, but make sure python version 3.7.7 is actually used for eggnog mapper analysis
```

Once complete, then we can run eggNOG-mapper as an sbatch script or in the terminal if the samples are small enough.

```
module load StdEnv/2020 gcc/9.3.0
module load hmmer diamond prodigal
module load python/3.7.7
module load mmseqs2/13-45111
for file in *_output.faa; do name=${file%%_*}; ../../08-EGGNOGGMAPPER/eggnog-mapper-2.1.6/emapper.py -i $file -o ../../08-EGGNOGGMAPPER/${name}; done
# Running it simply like this with emapper.py and no other commands means it will use diamond database as default
```

This will produce ```.emapper.hits```, ```.emapper.annotations```, and ```.emapper.seed_orthologs``` files. KEGG ID information is found within ```.emapper.annotations```. From this file, extract the KEGG ids which is in the KEGG_KO column (they look something like "ko:K01100" or with multiple ids). Then once we have our last of KEGG ids from eggNOG-mapper, we can input them into the [__KEGG Database__](https://www.genome.jp/kegg/ko.html) to see what the functional names are associated with these KEGG terms.

# 6. CoNet - Integration of Pop Gen and Metagenomics

Finally, to integrate the population genetic and metagenomic component to infer potential risks to wild bee health, we analyzed the correlation between population and landscape-level data with potenital bee and plant pathogen diversity. 

We created a dataframe that had samples as the rows, and the columns included: pathogenicity (Bray-Curtis dissimilarity for just identified bee and plant pathogens), Yang's relatedness, environmental distance, and resistance distance. We ran this with our [__CoNet script__](https://github.com/kdbchau/Ceratina-calcarata-Metagenomics/blob/main/Scripts/conet.R).
