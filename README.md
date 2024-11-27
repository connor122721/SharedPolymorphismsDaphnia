# *Balancing selection and the functional effects of shared polymorphism in cryptic Daphnia species*
## BioRxiv DOI: https://doi.org/10.1101/2024.04.16.589693

This repository hosts the scripts utilized for generating and analyzing the shared polymorphism dataset of *Daphnia* species.

### Author:
Connor S. Murray  
Email: csm6hg@virginia.edu

### Script Repository:

#### Figures:
Contains the main figures along with data and scripts used for plotting.

#### Supplemental:
Contains supplemental figures, tables, data, and associated scripts.

#### VCF Generation:
This section encompasses the pipeline for generating Variant Call Format (VCF) files:

- **VCF Creation**: 
  - Downloading short-read Illumina fastqs.
  - Mapping to European *Daphnia pulex* (D84A; GenBank assembly: GCA_023526725.1) genome.
  - BAM generation and merging.
  - HaplotypeCaller, GATK processing, merging, and SNP calling.
  - VCF and genomic data structure (GDS) generation.

- **VCF Filtration**:
  - Filter to benchmark universal single-copy ortholog (BUSCO) genes.
  - Multi-locus genotype (MLG) classification.
  - Filtering SNPs based on recommendations for species lacking reference SNP panels.
    
- **Quality Control**:
  - Missingness, quality scores, coverage, and MultiQC on BAMs.

- **Read Based Phasing**
  - WhatsHap phasing of samples and phased VCF creation.
 
- **Lift Over**
  - LiftOver of European and North American *Daphnia pulex* (KAP4; GenBank assembly: GCF_021134715.1) genomes.
  - Repeat masking, single-copy orthologous genes, and genotype concordance.
