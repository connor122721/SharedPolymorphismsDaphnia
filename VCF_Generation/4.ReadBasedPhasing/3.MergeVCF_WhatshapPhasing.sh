#!/usr/bin/env bash
##SBATCH -J phasing_mergeVCF # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -t 0-02:00:00 # Running time of 2 hours
#SBATCH --mem 50G # Memory request of 10GB
#SBATCH -o /project/berglandlab/connor/err/popPhasing_mergeVCF.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/popPhasing_mergeVCF.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load Modules
module load gcc/9.2.0
module load htslib/1.10.2
module load tabix/0.2.6
module load bcftools

# Parameter file
intervals="/project/berglandlab/connor/metadata/goodChrom.txt"

# Working directory for main VCF
wd="/project/berglandlab/connor/new_phase"

# Chromosomes to analyze (1:12)
chrom=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

# CPUs per job
threads=10

# Progress message
echo ${chrom}

# Whatshap list of VCFs
ls -d ${wd}/${chrom}/*vcf.gz > \
${wd}/${chrom}.whatshap.list

# Merge using bcftools into chromosomal vcf
bcftools \
merge \
-l ${wd}/${chrom}.whatshap.list \
-o ${wd}/${chrom}.filtered.whatshap.vcf \
-O v \
--threads ${threads}

# Index bcf
bcftools \
index \
-f \
--threads ${threads} \
${wd}/${chrom}.filtered.whatshap.vcf

# Finish
echo "Finish"

done
