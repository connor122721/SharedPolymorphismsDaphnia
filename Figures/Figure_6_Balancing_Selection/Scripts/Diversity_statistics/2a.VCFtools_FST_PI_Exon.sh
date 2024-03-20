#!/usr/bin/env bash
#SBATCH -J exon_fst_pi_vcftools
#SBATCH --ntasks-per-node=1
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00 # hours
#SBATCH --mem 1G
#SBATCH -o /project/berglandlab/connor/err/exon_vcf.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/exon_vcf.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# start
echo "Start" $(date)

# Exon Functions
exonFunc () {

# Load Modules
module load vcftools
module load tabix

# Working folder is core folder where this pipeline is being run.
wd=/project/berglandlab/connor/new_vcf2

# Input file
IN_GZVCF=${wd}/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.vcf.gz

# Master Bed file
bed=/project/berglandlab/connor/data/miss10.daphnia.pulex.merged.RMoutHiCGM.NsandDepthandChrEnd.final.merge50.bed

# Meta Exon
exon=($( ls /project/berglandlab/connor/data/exon/* ))

# Move to working directory
cd ${wd}

# Go through each exon
while read -r i; do
#i="Scaffold_7757_HRSCAF_8726       3472398 3472582 exon    Daphnia04572-RA"

# Start
echo ${i}
date

# Extract metadata
chrom=$( echo ${i} | cut -f1 -d " " )
start=$( echo ${i} | cut -f2 -d " " )
stop=$( echo ${i} | cut -f3 -d " " )
gene=$( echo ${i} | cut -f5 -d " " )

# Calculate length
len=$(echo $i | awk '{$6 = $3-$2 } 1' | cut -f6 -d " ")

# VCF functions
analy="--weir-fst-pop species.pulex.euro.150samps_dbunk_d8.pop.fst.txt \
--weir-fst-pop species.pulex.euro.85samps_gbr.pop.fst.txt --fst-window-size ${len}"

# Output name
out_namey="euro_fst_nomiss"

# Subset VCF and Run VCFTools
tabix -h ${IN_GZVCF} ${chrom}:${start}-${stop} |
vcftools \
--vcf - \
--max-missing 1 \
--chr ${chrom} \
--from-bp ${start} \
--to-bp ${stop} \
--exclude-positions ${bed} \
`echo ${analy}` \
-c | \
awk 'NR>1 {print $0 "\t" "'"$start"'" "\t" "'"$stop"'" "\t" "'"$gene"'"}' \
>> ${wd}/vcftools/${out_namey}_${SLURM_ARRAY_TASK_ID}

# Finish exon
done < ${exon[${SLURM_ARRAY_TASK_ID}]}
}

# Exports function
export -f exonFunc

# Go through each function
exonFunc ${SLURM_ARRAY_TASK_ID}

echo "VCF completed" $(date)
