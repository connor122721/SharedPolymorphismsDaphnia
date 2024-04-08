#!/bin/bash
# Script used to prepare VCF and run snapp for phylogenetic splits
# 11.1.2022
# Connor Murray

# Load modules
module load bcftools
module load ruby

# Working directory
wd="/project/berglandlab/connor/snapp5"

# VCF
vcf="/project/berglandlab/connor/new_vcf2/daphnia.filtered.chr.busco.vcf.gz"

# Beast2 directory
beast2="/home/csm6hg/beast/bin/beast"

# 5 indviduals per species - highest mean coverage
samps=${wd}/highdep.2inds.species.continent.list

# Starting tree: base_tree.nwk
#(((((Daphnia.pulex.NorthAmerica,Daphnia.pulexcaria.NorthAmerica,Daphnia.pulicaria.NorthAmerica),Daphnia.pulicaria.Europe),Daphnia.pulex.Europe),Daphnia.obtusa.NorthAmerica),Daphnia.obtusa.Europe)

# Extract 1 individual per species
srun --ntasks=1 --cpus-per-task=15 --time=02:00:00 --account=berglandlab -p standard \
bcftools view \
--threads 15 \
-S ${samps} \
-Ov \
-o ${wd}/daphnia.genome.2inds.vcf \
${vcf}

# Remove biallelic SNPs and no alternative alleles
srun --ntasks=1 --cpus-per-task=15 --time=02:00:00 --account=berglandlab -p standard \
bcftools view \
--threads 15 \
-e 'AC==0 || AC==AN' \
-m2 \
-M2 \
-Ov \
-o ${wd}/daphnia.genome.2inds.biallelic.vcf \
${wd}/daphnia.genome.2inds.vcf

# Run ruby input script
# From: https://github.com/mmatschiner/snapp_prep
ruby ${wd}/snapp_prep.rb \
-v ${wd}/daphnia.genome.2inds.biallelic.vcf \
-t ${wd}/individuals.2inds.continent.withObtusa.txt \
-c ${wd}/constraints.withobtusa.txt \
-s ${wd}/base_tree.nwk \
-x snapp.0.1.xml \
-o snapp.0.1 \
-m 1000 \
-l 1000000

# Open resutls with Tracer
/home/csm6hg/tracer/bin/tracer

# Open density tree with beast2
/home/csm6hg/beast/bin/densitree
