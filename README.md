# *Trans-Specific Polymorphisms Between Cryptic Daphnia Species Affect Fitness and Behavior*

![License: MIT](https://img.shields.io/badge/License-MIT-blue)
![Repo Size](https://img.shields.io/github/repo-size/connor122721/SharedPolymorphismsDaphnia)
![# Languages](https://img.shields.io/github/languages/count/connor122721/SharedPolymorphismsDaphnia)
![Top Language](https://img.shields.io/github/languages/top/connor122721/SharedPolymorphismsDaphnia)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2Fmec.17632-blue)](https://doi.org/10.1111/mec.17632)

## Overview

This repository hosts scripts used to generate and analyze the shared polymorphism dataset of *Daphnia* species.

### Author
Connor S. Murray  
Email: [csm6hg@virginia.edu](mailto:csm6hg@virginia.edu)

### Repository Structure

#### Figures
This directory contains the primary figures along with data and scripts used for plotting.

#### Supplemental
This directory contains supplemental figures, tables, data, and associated scripts.

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

### Container Setup and Usage

To improve reproducibility and ease of deployment, we have added an Apptainer container that includes the necessary dependencies for running R and Python scripts in this repository.

#### Installing Apptainer

If you do not have Apptainer installed, follow these steps to install it:

1. Install dependencies:
   ```sh
   sudo apt-get update
   sudo apt-get install -y build-essential libseccomp-dev pkg-config squashfs-tools cryptsetup
   ```

2. Download and install Go:
   ```sh
   wget https://dl.google.com/go/go1.16.5.linux-amd64.tar.gz
   sudo tar -C /usr/local -xzf go1.16.5.linux-amd64.tar.gz
   export PATH=$PATH:/usr/local/go/bin
   ```

3. Download the Apptainer source code:
   ```sh
   git clone https://github.com/apptainer/apptainer.git
   cd apptainer
   ```

4. Compile and install Apptainer:
   ```sh
   ./mconfig
   make -C builddir
   sudo make -C builddir install
   ```

#### Building the Apptainer Container

1. Clone this repository to your local machine:
   ```sh
   git clone https://github.com/connor122721/SharedPolymorphismsDaphnia.git
   cd SharedPolymorphismsDaphnia
   ```

2. Build the Apptainer image:
   ```sh
   apptainer build daphnia-container.sif Apptainer.def
   ```

#### Running the Apptainer Container

1. Start the container:
   ```sh
   apptainer exec daphnia-container.sif R
   ```

## Notes
- I would really love to make the feature of identifying trans-specific polymorphisms a package in R or as a standalone software. I am currently editing this feature so it is more user-friendly, but I would appreciate any collaboration or tips!
  
## License
This project is licensed under the MIT License.
