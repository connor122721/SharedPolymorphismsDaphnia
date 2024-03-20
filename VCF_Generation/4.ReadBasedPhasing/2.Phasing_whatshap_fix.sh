#!/usr/bin/env bash
#SBATCH -J whatshap-phase-chrom # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 day
#SBATCH --mem 10G
#SBATCH -o /project/berglandlab/connor/err/phylo.whatshap.chrom.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/phylo.whatshap.chrom.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/2020.11-py3.8 samtools htslib bcftools/1.9 gparallel/20170822

# Install whatshap (run once)
# pip3 install --user whatshap
export PATH=$HOME/.local/bin:$PATH

# Working directory
wd="/project/berglandlab/connor"

# Reference genome
ref="${wd}/totalHiCwithallbestgapclosed.fa"

# Parameter file
# 1st column = sample, 2nd column = chromosome, seperated = ","
parameterFile="${wd}/new_phase/phasing_paramList_phylo"

# Working directory for bams
bam="/scratch/csm6hg/from_old_scratch/all_bam"

# Extract sample name
samp=$( cat ${parameterFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )

# Chromosomes to analyze
chrom=$( cat ${parameterFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )

# CPUs per job
threads=1

# Message for sample and chromosome
echo "Read-back phasing -" "Sample:" ${samp} $SLURM_ARRAY_TASK_ID
echo "Chromosome:" ${chrom}

# Create temporary folders
if [[ -d "${wd}/new_phase" ]]
then
	echo "Working tmp folder exist"
	echo "lets move on"
	date
else
	echo "Folder doesnt exist. Let us fix that."
	mkdir ${wd}/new_phase
	date
fi

# Create temporary folders per chromosome
if [[ -d "${wd}/new_phase/${chrom}" ]]
then
	echo "Working tmp chromosome folder exist"
	echo "lets move on"
	date
else
	echo "Folder doesnt exist. Let us fix that."
	mkdir ${wd}/new_phase/${chrom}
	date
fi

echo "Extracting chromosome:" ${chrom} "from:" ${samp} "#" ${SLURM_ARRAY_TASK_ID}

# bcftools to extract individual and chromosome into an individual VCF
bcftools view \
    -s ${samp} \
    -r ${chrom} \
    -O v \
    ${wd}/new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.vcf.gz | \
    awk '{
      a=0
      if(substr($0, 0, 1)=="#") {
        print $0
      } else {
        printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t.\tGT\t"
        split($10, sp, ":")
        printf sp[1]"\n"
      }
    }' > \
    ${wd}/new_phase/${chrom}/${samp}_${chrom}.vcf

# Extract read groups from bam file per chromosome per sample
echo "Getting bam"

# Dorthe dataset
if [ -f ${bam}/bam_dorthe_fin/${samp}*RG.bam ]; then
samtools view \
    -b \
    -q 20 \
    -M \
    -O BAM \
    -@ ${threads} \
    ${bam}/bam_dorthe_fin/${samp}_finalmap_RG.bam \
    ${chrom} > \
    ${wd}/new_phase/${chrom}/${samp}_${chrom}.bam

# Index bam
samtools index ${wd}/new_phase/${chrom}/${samp}_${chrom}.bam

# Add read groups - hyphens are illegal!
Group_library="SingleInds"
Library_Platform="illumina"
Group_platform="place"
module load picard/2.23.4

# Replace - with _
samp1=$( echo "$samp" | tr '-' '_' )

echo $samp $samp1 > \
${wd}/new_phase/${chrom}/$samp.rename

# Rename sample in VCF without hyphen
bcftools reheader \
${wd}/new_phase/${chrom}/${samp}_${chrom}.vcf \
-s ${wd}/new_phase/${chrom}/$samp.rename \
-o ${wd}/new_phase/${chrom}/${samp1}_${chrom}.vcf

# Add corrected read groups
java -jar $PICARD AddOrReplaceReadGroups \
I=${wd}/new_phase/${chrom}/${samp}_${chrom}.bam \
O=${wd}/new_phase/${chrom}/${samp}_${chrom}.sm.bam \
RGID=4 \
RGLB=$Group_library \
RGPL=$Library_Platform \
RGPU=$Group_platform \
RGSM=$samp1

mv ${wd}/new_phase/${chrom}/${samp}_${chrom}.sm.bam \
${wd}/new_phase/${chrom}/${samp1}_${chrom}.bam

samp=$samp1

fi

# SRA dataset
if [ -f ${bam}/final_bam/${samp}*_RG.bam ]; then
samtools view \
    -b \
    -q 20 \
    -M \
    -O BAM \
    -@ ${threads} \
    ${bam}/final_bam/${samp}_finalmap_RG.bam \
    ${chrom} > \
    ${wd}/new_phase/${chrom}/${samp}_${chrom}.bam
fi

# European dataset
if [ -f ${bam}/Euro_bams/${samp}_finalmap_mdup.bam ]; then
samtools view \
    -b \
    -q 20 \
    -M \
    -O BAM \
    -@ ${threads} \
    ${bam}/Euro_bams/${samp}_finalmap_mdup.bam \
    ${chrom} > \
    ${wd}/new_phase/${chrom}/${samp}_${chrom}.bam
fi

# Index bam
samtools index ${wd}/new_phase/${chrom}/${samp}_${chrom}.bam

echo "Running whatshap"

# Run whatshap
whatshap \
  phase \
  -r ${ref} \
  -o ${wd}/new_phase/${chrom}/${samp}.${chrom}.phase.vcf.gz \
  --chromosome ${chrom} \
  --sample ${samp} \
  ${wd}/new_phase/${chrom}/${samp}_${chrom}.vcf \
  ${wd}/new_phase/${chrom}/${samp}_${chrom}.bam

# Tabix vcf
module load tabix

# tabix
tabix -p vcf ${wd}/new_phase/${chrom}/${samp}.${chrom}.phase.vcf.gz

# Check if output files are below 1Mb threshold
#find ${wd}/new_phase/${chrom}/${samp}.${chrom}.phase.vcf.gz -type f -size -1M | sed -r 's/.{39}$//' > ../failed.samples

# Clean up temporary files
rm ${wd}/new_phase/${chrom}/${samp}_${chrom}.vcf
rm ${wd}/new_phase/${chrom}/${samp}_${chrom}*.bam
rm ${wd}/new_phase/${chrom}/${samp}_${chrom}*.bam.bai

# Finish
echo "Finish:" ${samp} ${chrom} $SLURM_ARRAY_TASK_ID

# Lists number of files lower than 10Mb per subfolder
# for i in */ .*/ ; do echo -n $i": " ; (find "$i" -type f -size -1M | wc -l) ;  done

# Lists extensions for array jobs that failed
# for i in */ ; do chrom=$( echo $i | sed 's/.$//' ); find $i*.phase.vcf -type f -size -5M >> failed.samples.new1;  done

# Remove files for array jobs that failed
# for i in */ ; do chrom=$( echo $i | sed 's/.$//' ); find $i*.bam -type f -exec rm {} +; done
