# *Trans-Specific Polymorphisms Between Cryptic Daphnia Species Affect Fitness and Behavior*
## [Molecular Ecology - DOI: 10.1101/2024.04.16.589693](https://doi.org/10.1101/2024.04.16.589693)

This repository hosts scripts used to generate and analyze the shared polymorphism dataset of *Daphnia* species.

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

