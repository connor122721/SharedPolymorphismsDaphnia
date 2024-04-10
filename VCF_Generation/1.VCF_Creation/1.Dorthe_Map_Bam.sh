#!/usr/bin/env bash
#
#SBATCH -J Dorthe_bam # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 3-00:00 # 3 days
#SBATCH --mem 100G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/err/map_dorthe.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/err/map.dorthe.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load gcc/9.2.0 htslib/1.10.2 
module load trimmomatic/0.36
module load bwa
module load samtools
module load picard

# Working directory
wd="/scratch/csm6hg/daphnia_phylo/bam_dorthe"
parameterFile="/scratch/csm6hg/daphnia_phylo/dorthe.samples"
outfq="/scratch/csm6hg/daphnia_phylo/bam_dorthe"
outbam="/scratch/csm6hg/daphnia_phylo/bam_dorthe_fin"
forbam="/scratch/csm6hg/daphnia_phylo/dorthe.forward.fastq"
revbam="/scratch/csm6hg/daphnia_phylo/dorthe.reverse.fastq"

# Parameters
threads=10

# Extract sample
samp=`sed -n ${SLURM_ARRAY_TASK_ID}p $parameterFile`
forw=`sed -n ${SLURM_ARRAY_TASK_ID}p $forbam`
rev=`sed -n ${SLURM_ARRAY_TASK_ID}p $revbam`

# Trim illumina adaptors
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads ${threads} \
${outfq}/${forw}.fastq.gz \
${outfq}/${rev}.fastq.gz \
${outfq}/${samp}_1.P.trimm.fastq.gz \
${outfq}/${samp}_1.U.trimm.fastq.gz \
${outfq}/${samp}_2.P.trimm.fastq.gz \
${outfq}/${samp}_2.U.trimm.fastq.gz \
ILLUMINACLIP:/scratch/csm6hg/daphnia_phylo/data/CombinedPE-PE.fa:2:30:10:8:true 
 
# Merge overlapping reads
/scratch/csm6hg/daphnia_phylo/pear \
-f ${outfq}/${samp}_1.P.trimm.fastq.gz \
-r ${outfq}/${samp}_2.P.trimm.fastq.gz \
-o ${outfq}/${samp} \
-j ${threads}
   
# Map to reference genome
bwa mem -t ${threads} -K 100000000 -Y \
${wd}/totalHiCwithallbestgapclosed.fa \
${outfq}/${samp}.assembled.fastq | \
samtools view -Suh -q 20 -F 0x100 | \
samtools sort --threads ${threads} -o ${outfq}/${samp}.sort.bam
samtools index ${outfq}/${samp}.sort.bam

# Unassembled reads
bwa mem -t ${threads} -K 100000000 -Y \
${wd}/totalHiCwithallbestgapclosed.fa \
${outfq}/${samp}.unassembled.forward.fastq \
${outfq}/${samp}.unassembled.reverse.fastq | \
samtools view -Suh -q 20 -F 0x100 | \
samtools sort --threads ${threads} -o ${outfq}/${samp}.filt.unassembled.sort.bam
samtools index ${outfq}/${samp}.filt.unassembled.sort.bam
	
# Merge assembled and unassembled bam files and mark duplicates
samtools merge ${outfq}/${samp}.filt.merged.bam \
${outfq}/${samp}.sort.bam \
${outfq}/${samp}.filt.unassembled.sort.bam

# Index merged bam
samtools index ${outfq}/${samp}.filt.merged.bam

# Mark duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=10000 \
REMOVE_DUPLICATES=true \
INPUT=${outfq}/${samp}.filt.merged.bam \
OUTPUT=${outfq}/${samp}_finalmap.bam \
METRICS_FILE=${outfq}/${samp}_finalmap_mdups.metrics \
CREATE_INDEX=true           

# Move output
mv ${outfq}/${samp}_finalmap* ${outbam}/

# Remove previous steps
rm ${outfq}/${samp}*

# Finish
echo "Finish"