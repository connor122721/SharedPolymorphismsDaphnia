#!/usr/bin/env bash
##SBATCH -J phasing_mergeVCF # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -t 0-02:00:00 # Running time of 2 hours
#SBATCH --mem 50G # Memory request of 10GB
#SBATCH -o /project/berglandlab/connor/err/popPhasing_concat.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/popPhasing_concat.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load Modules
module load gcc/9.2.0
module load htslib/1.10.2
module load tabix/0.2.6
module load bcftools

# Working directory for main VCF
wd="/project/berglandlab/connor/new_phase"

# CPUs per job
threads=10

# Whatshap list of VCFs
ls -d ${wd}/*filtered.whatshap.vcf > \
${wd}/whatshap.list

# Concatinate all whatshap.vcf files into common vcf
bcftools \
concat \
--threads ${threads} \
-f ${wd}/whatshap.list \
-Ov \
-o ${wd}/daphnia.whatshap.vcf

# Bgzip and tabix
bgzip ${wd}/daphnia.whatshap.vcf
tabix -p vcf ${wd}/daphnia.whatshap.vcf.gz

# Combined VCF name
vcf="daphnia.whatshap.vcf.gz"

# Output annotated VCF name
out_vcf="daphnia.whatshap.ann.vcf"

# java parameters
JAVAMEM=50G

# Run snpEFF on raw vcf
cd ${wd}
java -Xmx${JAVAMEM} -jar \
/home/csm6hg/SNPEFF/snpEff.jar ann \
dpgenome \
${wd}/${vcf} \
-o vcf \
-t \
-htmlStats ${wd}/snpeff_summary_whatshap.html > \
${out_vcf}

# BGzip and tab index annotated vcf
bgzip ${wd}/${out_vcf}
tabix -p vcf ${wd}/${out_vcf}.gz

# Finish
echo "Finish"
