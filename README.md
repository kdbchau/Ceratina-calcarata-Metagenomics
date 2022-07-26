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
* [SHARCNET](https://www.sharcnet.ca/my/front/) for technical support on high-performance computing (HPC) clusters ([Graham](https://docs.alliancecan.ca/wiki/Graham), [Beluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga), and [Cedar](https://docs.alliancecan.ca/wiki/Cedar))
    * [Jose Sergio Hleap]() for assistance with the HPC clusters and overall pipeline performance.

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
4. [Population Genetics]
    *
5. [Metagenomics]    
    * [metaSPADES](#51-metaspades)
    * [Diversity Statistics](#52-diversity-statistics)
    * [Random Forests](#53-random-forests)
    * [WGCNA](#54-wgcna)
    * [Functional Analyses](#55-functional-analyses)
        * [GhostKOALA](#551-ghostkoala)
        * [eggNOG-mapper](#552-eggnog-mapper)
6. [CoNet - Integration of Pop Gen and Metagenomics](#6-conet-integration-of-pop-gen-and-metagenomics)


