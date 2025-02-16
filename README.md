# *Trans-Specific Polymorphisms Between Cryptic Daphnia Species Affect Fitness and Behavior*

![License: MIT](https://img.shields.io/badge/License-MIT-blue)
![Repo Size](https://img.shields.io/github/repo-size/connor122721/SharedPolymorphismsDaphnia)
![# Languages](https://img.shields.io/github/languages/count/connor122721/SharedPolymorphismsDaphnia)
![Top Language](https://img.shields.io/github/languages/top/connor122721/SharedPolymorphismsDaphnia)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2Fmec.17632-blue)](https://doi.org/10.1111/mec.17632)

## This repository hosts scripts used to generate and analyze the shared polymorphism dataset of *Daphnia* species.

### Author
Connor S. Murray  
Email: [csm6hg@virginia.edu](mailto:csm6hg@virginia.edu)

### Repository Overview

#### Figures
Contains the primary figures along with data and scripts used for plotting.

#### Supplemental
Contains supplemental figures, tables, data, and associated scripts.

#### VCF Generation
This section includes the pipeline for generating Variant Call Format (VCF) files:

- **VCF Creation**:  
  - Downloading short-read Illumina FASTQ files.  
  - Mapping to the European *Daphnia pulex* genome (D84A; GenBank assembly: GCA_023526725.1).  
  - BAM generation and merging.  
  - SNP calling using HaplotypeCaller and GATK.  
  - VCF and genomic data structure (GDS) generation.  

- **VCF Filtration**:  
  - Filtering to benchmark universal single-copy ortholog (BUSCO) genes.  
  - Multi-locus genotype (MLG) classification.  
  - SNP filtering based on recommendations for species lacking reference SNP panels.  

- **Quality Control**:  
  - Assessing missingness, quality scores, and coverage.  
  - MultiQC on BAMs.

- **Read-Based Phasing**:  
  - WhatsHap phasing of samples and phased VCF creation.

- **LiftOver**:  
  - Lifting over the European and North American *Daphnia pulex* genomes (KAP4; GenBank assembly: GCF_021134715.1).  
  - Repeat masking, single-copy orthologous genes, and genotype concordance.

## Notes
- I would really love to make the feature of identifying trans-specific polymorphisms a package in R or as a standalone sfotware. I am currently editing this feature so it is more user-friendly, but I would appreciate any help or tips!
  
## License
This project is licensed under the MIT License.
