#!/usr/bin/env bash
#SBATCH -J extractFASTA # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-05:00 # 4 hours
#SBATCH --mem 50G
#SBATCH -o /project/berglandlab/connor/err/extractFASTA.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/extractFASTA.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Make metadata file
cd /scratch/csm6hg/data
module load bcftools
bcftools query -l  combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.vcf.gz \
| grep "SRR" -v | grep "ERR" -v > \
euro.samples.filt

# Make bed file
echo "Scaffold_1931_HRSCAF_2197:6350219-6351796" > blop.bed

SLURM_ARRAY_TASK_ID=1

# Extract chromosome from each individual
getRegionInd () {

### load modules
module load samtools
module load bcftools
# iteration=30

# Working directory
wd="/scratch/csm6hg"

# Parameter file
# 1st column = sample
parameterFile="${wd}/data/euro.samples.filt"

# Bed file
chromFile="${wd}/data/blop.bed"

# Job id information
samp=$( sed "${iteration}q;d" ${parameterFile} | cut -f1 )

# Reference genome
ref="/scratch/csm6hg/ref/totalHiCwithallbestgapclosed.fa"

# Message for sample
echo "Sample:" ${samp} "id:" ${iteration} "Pop:" ${pop}
echo "${iteration}"

# Runs function for each haplotype
getHaplos () {
whichHaplo=${1}
samp=${2}
chrom=${3}
#whichHaplo=2;chrom=Scaffold_6786_HRSCAF_7541:13195666-13199143

# Message
echo ${samp} ${whichHaplo} ${chrom} ${pop}

# Parameters
chr=$( echo ${chrom} | tr ":" '\t' | cut -f1 )
start=$( echo ${chrom} | tr ":" '\t' | cut -f2 | tr "-" '\t' | cut -f1)
stop=$( echo ${chrom} | tr ":" '\t' | cut -f2 | tr "-" '\t' | cut -f2 )

# Message for chromosome
echo "Chromosome:" ${chr} ${start} ${stop}

# Creates a consensus fasta file for each haplotype
bcftools \
consensus \
-f <(samtools faidx ${ref} ${chrom}) \
-H ${whichHaplo} \
-s ${samp} \
-o ${wd}/blop/fasta/${samp}.${whichHaplo}.${chr}.${start}.${stop}.fa \
${wd}/data/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.vcf.gz

# Add sample name and haplotype to header
f="${samp}.${whichHaplo}.${chr}.${start}.${stop}"

# Extracts chromosome via bed file
sed -i "s/^>*${chrom}/>${f}/g" \
${wd}/blop/fasta/${samp}.${whichHaplo}.${chr}.${start}.${stop}.fa

# Index new fasta
samtools faidx ${wd}/blop/fasta/${samp}.${whichHaplo}.${chr}.${start}.${stop}.fa

# Show header of new FASTA
head -n 1 ${wd}/blop/fasta/${samp}.${whichHaplo}.${chr}.${start}.${stop}.fa

}

# Exports function
export -f getHaplos

# Parallelize tasks
parChr() {

# Runs "Extract haplotype" function
getHaplos 1 ${samp} ${i}
getHaplos 2 ${samp} ${i}

chr=$( echo ${chrom} | tr ":" '\t' | cut -f1 )
start=$( echo ${i} | tr ":" '\t' | cut -f2 | tr "-" '\t' | cut -f1 )
stop=$( echo ${i} | tr ":" '\t' | cut -f2 | tr "-" '\t' | cut -f2 )

}

# For loop through each window
i=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${chromFile} )
parChr "${i}"

}

# Export individual function
export getRegionInd

# For loop through each individual
for iteration in {1..611}; do getRegionInd "${iteration}"; done

# Finish script
echo "Finish"
