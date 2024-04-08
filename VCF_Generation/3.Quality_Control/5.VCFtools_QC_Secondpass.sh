#! /bin/bash
#
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab_standard
#SBATCH -o /project/berglandlab/connor/err/FinQCVCF.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/FinQCVCF.err # Standard error

# This script will check various QC statistics from a raw VCF file

# Load Modules
module load vcftools

# Working folder is core folder where this pipeline is being run.
wd=/project/berglandlab/connor/new_vcf2

# Input file
IN_GZVCF=${wd}/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.vcf.gz

# Move to working directory
cd ${wd}

# Survey quality of final VCF
analyses=("--weir-fst-pop species.pulex.nam.pop.fst.txt --weir-fst-pop species.pulex.euro.pop.fst.txt" "--site-mean-depth --keep species.pop.fst.txt")

# For loop to run through 6 QC steps
for i in {2}

do
i=2
echo "Now processing" ${analyses[${i}]}

# Run VCFTools
vcftools \
--gzvcf $IN_GZVCF \
`echo ${analyses[${i}]}` \
--out ${wd}/vcftools

# Finish i
done

echo "VCF QC completed" $(date)
