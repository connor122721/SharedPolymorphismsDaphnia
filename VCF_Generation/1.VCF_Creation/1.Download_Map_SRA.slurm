#!/usr/bin/env bash
#
#SBATCH -J DownloadMap # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/err/down.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/err/down.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load sratoolkit/2.10.5
module load gcc/9.2.0 htslib/1.10.2
module load trimmomatic/0.36
module load bwa
module load samtools
module load picard

# Working directory
wd="/scratch/csm6hg/daphnia_phylo"
parameterFile="/scratch/csm6hg/daphnia_phylo/Daphnia.SRR.csv"
outfq="/scratch/csm6hg/daphnia_phylo/fqs"
outbam="/scratch/csm6hg/daphnia_phylo/bam"
threads=10

# Extract sample
samp=$( cat ${parameterFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )

# Download samples
fasterq-dump ${samp} --outdir ${outfq} --threads 10
bgzip ${outfq}/${samp}_1.fastq --threads 10
bgzip ${outfq}/${samp}_2.fastq --threads 10

# Fastq Quality control
module load singularity
singularity run /home/yey2sn/software/multiqc.sif ${outfq} -o ${wd}/mapping_stats

# Trim illumina adaptors
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads ${threads} \
${outfq}/${samp}_1.fastq.gz \
${outfq}/${samp}_2.fastq.gz \
${outfq}/${samp}_1.P.trimm.fastq.gz \
${outfq}/${samp}_1.U.trimm.fastq.gz \
${outfq}/${samp}_2.P.trimm.fastq.gz \
${outfq}/${samp}_2.U.trimm.fastq.gz \
ILLUMINACLIP:/project/berglandlab/connor/TrimmomaticAdaptors/CombinedPE-PE.fa:2:30:10:8:true

# Merge overlapping reads
${wd}/pear \
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
REMOVE_DUPLICATES=true \
INPUT=${outfq}/${samp}.filt.merged.bam \
OUTPUT=${outfq}/${samp}_finalmap.bam \
METRICS_FILE=${outfq}/${samp}_finalmap_mdups.metrics \
CREATE_INDEX=true

# Move output
mv ${outfq}/${samp}_finalmap* ${outbam}/

# Remove previous steps
rm ${outfq}/${samp}

# Finish
echo "Finish"
