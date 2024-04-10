#!/usr/bin/env bash
#SBATCH -J smcpp-barnacle # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -t 1-00:00 # Running time of 1 hour
#SBATCH --mem 50G # Memory request of 100GB
#SBATCH -o /project/berglandlab/connor/msmc/err/smcpp.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/msmc/err/smcpp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# Load SMC++ & Modules
module load htslib
module load gnuplot/5.2.2

# Working directory
wd="/project/berglandlab/connor"

# Samples
paramFile=${wd}/msmc/mlg.1.species.country.list

# Extract constants from parameter file
sample=$( sed -n ${SLURM_ARRAY_TASK_ID}p $paramFile )

# Chromosome file
intervals=${wd}/metadata/goodChrom.txt

# Reference genome
ref=${wd}/totalHiCwithallbestgapclosed.fa

# Cores
threads=10

# Create combined chromosome samples
while read -r i stop; do
  # i="Scaffold_9201_HRSCAF_10758"

  # Progress message
  echo ${i}

  # Index VCFs
  gunzip ${wd}/msmc/Scaffold*/${sample}.${i}.phase.vcf.gz
  bgzip ${wd}/msmc/Scaffold*/${sample}.${i}.phase.vcf
  tabix -p vcf -f ${wd}/msmc/Scaffold*/${sample}.${i}.phase.vcf.gz

  # Sample 1D sfs - create input file
  /project/berglandlab/connor/smc++ vcf2smc \
  ${wd}/msmc/Scaffold*/${sample}.${i}.phase.vcf.gz \
  --cores ${threads} \
  --mask ${wd}/data/miss10.daphnia.pulex.merged.RMoutHiCGM.final.bed.gz \
  ${wd}/msmc/smcpp/${sample}.${i}.phase.smc.gz \
  ${i} \
  ${sample}:${sample}

done < ${intervals}

# Estimate demographic history
/project/berglandlab/connor/smc++ estimate \
-o ${wd}/msmc/smcpp \
--base ${sample}.phase \
--cores ${threads} \
--em-iterations 30 \
--timepoints 1000 5e6 \
5.69e-09 \
${wd}/msmc/smcpp/${sample}.*.phase.smc.gz

# Plot output
/project/berglandlab/connor/smc++ plot \
${wd}/msmc/smcpp/${sample}.filt.phase.smcpp.pdf \
${wd}/msmc/smcpp/${sample}.phase.final.json \
-g 1 \
-c

# Finish
echo "Finish"
