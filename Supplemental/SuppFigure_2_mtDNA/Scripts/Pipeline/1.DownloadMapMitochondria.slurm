#!/usr/bin/env bash
#SBATCH -J MapMito # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 days
#SBATCH --mem 100G
#SBATCH -o /scratch/csm6hg/mito/err/map.mito.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/mito/err/map.mito.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load sratoolkit/2.10.5
module load gcc/9.2.0 htslib/1.10.2
module load trimmomatic/0.36
module load bwa/0.7.17
module load samtools
module load picard
module load fastqc
module load bcftools

# Working directory
wd="/scratch/csm6hg/mito"

# Sample names
paramFile="${wd}/samples_mito"

# Final output - mitochondrial bam and fasta
outbam="${wd}/bam"
outfasta="${wd}/mito.fasta"
outvcf="${wd}/vcf"

# Number of CPUs
threads=10

# Extract constants from parameter file
slurmID=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f1 )
samp=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
species=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )
continent=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )
bam=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f5 )
fastaname=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f6 )
fastawild=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f7 )
Origin=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f8 )
fastfor=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f9 )
fastrev=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f10 )

# Conduct VCF calling
if [ ! -f "${outvcf}/${samp}*filt.mito.vcf.gz" ]; then
    echo "${samp} does not exist."

# Reference genome
if [[ (${species} == "Daphnia pulex" && ${continent} == "NorthAmerica") || \
  (${species} == "Daphnia pulicaria" && ${continent} == "NorthAmerica") ]]; then
  echo ${samp} "mapping to NAm. Daphnia pulex"
  ref="/scratch/csm6hg/Daphnia_mtdna/nam_pulex/NC_000844.1_Daphnia_pulex_mitochondrion.fasta"
fi

# Reference genome
if [[ (${species} == "Daphnia pulex" && ${continent} == "Europe") || \
  (${species} == "Daphnia pulicaria" && ${continent} == "Europe") ]]; then
  echo ${samp} "mapping to Euro. Daphnia pulex"
  ref="/scratch/csm6hg/Daphnia_mtdna/euro_pulex/D84Amtdna.fa"
fi

# Reference genome
if [[ ${species} == "Daphnia obtusa" ]]; then
  echo ${samp} "mapping to Daphnia obtusa"
  ref="/scratch/csm6hg/Daphnia_mtdna/obtusa/nam_obtusa_mtdna.fasta"
fi

# Sample in SRA
if [[ ${Origin} == "SRA" ]]; then

  # Progress
  echo "Sample is from SRA"

  # Output directory
  outfq="/scratch/csm6hg/all_fqs/SRA_fqs"
  cd ${outfq}

  # Conditional download
  if [ ! -f "${outfq}/${samp}_1.fastq.gz" ]; then
      echo "${samp} does not exist. Downloading now."

    # Get Fastq
    fasterq-dump ${samp} \
    --outdir ${outfq} \
    --threads ${threads}

    # FASTQC
    fastqc ${outfq}/${samp}_1.fastq \
    -t ${threads}

    fastqc ${outfq}/${samp}_2.fastq \
    -t ${threads}

    mv ${outfq}/${samp}*fastqc* \
    ${outfq}/qc

    # BGzip fastq
    bgzip ${outfq}/${samp}_1.fastq \
    --threads ${threads}

    bgzip ${outfq}/${samp}_2.fastq \
    --threads ${threads}

  fi

  # Message
  echo "Finished fastq download and QC"

  # Trim illumina adaptors
  java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE \
  -threads ${threads} \
  ${outfq}/${samp}_1.fastq.gz \
  ${outfq}/${samp}_2.fastq.gz \
  ${outfq}/${samp}_1.P.trimm.fastq.gz \
  ${outfq}/${samp}_1.U.trimm.fastq.gz \
  ${outfq}/${samp}_2.P.trimm.fastq.gz \
  ${outfq}/${samp}_2.U.trimm.fastq.gz \
  ILLUMINACLIP:${wd}/data/CombinedPE-PE.fa:2:30:10:8:true

  # Merge overlapping reads
  /scratch/csm6hg/pear \
  -f ${outfq}/${samp}_1.P.trimm.fastq.gz \
  -r ${outfq}/${samp}_2.P.trimm.fastq.gz \
  -o ${outfq}/${samp} \
  -j ${threads}

  # Message
  echo "Mapping now"

  # Map to mitochondrial genome
  bwa mem -t ${threads} -K 100000000 -Y \
  ${ref} \
  ${outfq}/${samp}.assembled.fastq | \
  samtools view -Suh -q 20 -F 0x100 | \
  samtools sort --threads ${threads} -o ${outfq}/${samp}.sort.bam
  samtools index ${outfq}/${samp}.sort.bam

  # Unassembled reads
  bwa mem -t ${threads} -K 100000000 -Y \
  ${ref} \
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

  # Call snps
  bcftools mpileup -q 30 -Q 20 -Ou -f ${ref} \
  ${outfq}/${samp}_finalmap.bam | \
  bcftools call -mv -V indels | bcftools filter \
  -i '%QUAL>20 & DP>=100' -Oz -o ${outfq}/${samp}.filt.mito.vcf.gz

  # Index VCF
  tabix -p vcf ${outfq}/${samp}.filt.mito.vcf.gz

  # apply variants to create consensus sequence
  cat ${ref} | \
  bcftools consensus ${outfq}/${samp}.filt.mito.vcf.gz \
  --sample ${outfq}/${samp}_finalmap.bam > \
  ${outfq}/${samp}.filt.consensus.mito.fa

  # Chnage name of header
  name=$( echo ${samp}.filt.consensus.mito.fa | tr "." '\t' | cut -f1 | tr "/" '\t' | cut -f1 )
  echo ${name}

  # Extracts chromosome via bed file
  sed -i "s/^>*"mtdna"/>"mtdna."${name}/g" \
  ${outfq}/${samp}.filt.consensus.mito.fa

  # Index fasta
  samtools faidx ${outfq}/${samp}.filt.consensus.mito.fa

  # Move output
  mv ${outfq}/${samp}_finalmap* ${outbam}/
  mv ${outfq}/${samp}*mito.fa* ${outfasta}/
  mv ${outfq}/${samp}*mito.vcf* ${outvcf}/

  # Remove previous steps
  rm ${outfq}/${samp}*trimm*
  rm ${outfq}/${samp}*unassembled*
  rm ${outfq}/${samp}*assembled*
  rm ${outfq}/${samp}*discarded*
  rm ${outfq}/${samp}*merged*
  rm ${outfq}/${samp}*sort.bam*

else
  echo "Sample not in SRA dataset"
fi

# Sample in Dorthe dataset
if [[ ${Origin} == "Dorthe" ]]; then

  # Progress
  echo "Sample is from Dorthe"
  echo ${samp}
  outfq="/scratch/csm6hg/all_fqs/dorthe"

  # Trim illumina adaptors
  java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE \
  -threads ${threads} \
  ${outfq}/${fastfor}.fastq.gz \
  ${outfq}/${fastrev}.fastq.gz \
  ${outfq}/${samp}_1.P.trimm.fastq.gz \
  ${outfq}/${samp}_1.U.trimm.fastq.gz \
  ${outfq}/${samp}_2.P.trimm.fastq.gz \
  ${outfq}/${samp}_2.U.trimm.fastq.gz \
  ILLUMINACLIP:${wd}/data/CombinedPE-PE.fa:2:30:10:8:true

  # Merge overlapping reads
  /scratch/csm6hg/pear \
  -f ${outfq}/${samp}_1.P.trimm.fastq.gz \
  -r ${outfq}/${samp}_2.P.trimm.fastq.gz \
  -o ${outfq}/${samp} \
  -j ${threads}

  # Message
  echo "Mapping now"

  # Map to mitochondrial genome
  bwa mem -t ${threads} -K 100000000 -Y \
  ${ref} \
  ${outfq}/${samp}.assembled.fastq | \
  samtools view -Suh -q 20 -F 0x100 | \
  samtools sort --threads ${threads} -o ${outfq}/${samp}.sort.bam
  samtools index ${outfq}/${samp}.sort.bam

  # Unassembled reads
  bwa mem -t ${threads} -K 100000000 -Y \
  ${ref} \
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

  # Call snps
  bcftools mpileup -q 30 -Q 20 -Ou -f ${ref} \
  ${outfq}/${samp}_finalmap.bam | \
  bcftools call -mv -V indels | bcftools filter \
  -i '%QUAL>20 & DP>=100' -Oz -o ${outfq}/${samp}.filt.mito.vcf.gz

  # Index VCF
  tabix -p vcf ${outfq}/${samp}.filt.mito.vcf.gz

  # apply variants to create consensus sequence
  cat ${ref} | \
  bcftools consensus ${outfq}/${samp}.filt.mito.vcf.gz \
  --sample ${outfq}/${samp}_finalmap.bam > \
  ${outfq}/${samp}.filt.consensus.mito.fa

	# Chnage name of header
	name=$( echo ${samp}.filt.consensus.mito.fa | tr "." '\t' | cut -f1 | tr "/" '\t' | cut -f1 )
	echo ${name}

	# Rename fasta chromosome
	sed -i "s/^>*"mtdna"/>"mtdna."${name}/g" \
  ${outfq}/${samp}.filt.consensus.mito.fa

	# Index fasta
	samtools faidx ${outfq}/${samp}.filt.consensus.mito.fa

  # Move output
  mv ${outfq}/${samp}_finalmap* ${outbam}/
  mv ${outfq}/${samp}*mito.fa* ${outfasta}/
  mv ${outfq}/${samp}*mito.vcf* ${outvcf}/

  # Remove previous steps
  rm ${outfq}/${samp}*trimm*
  rm ${outfq}/${samp}*unassembled*
  rm ${outfq}/${samp}*assembled*
  rm ${outfq}/${samp}*discarded*
  rm ${outfq}/${samp}*merged*
  rm ${outfq}/${samp}*sort.bam*

else
  echo "Sample not in Dorthe dataset"
fi

# Sample in Euro dataset
if [[ ${Origin} == "Europe" ]]; then

  # Progress
  echo "Sample is from Europe"
  echo ${samp}
  outfq="/scratch/csm6hg/all_fqs/Euro_fqs"

  # Move to Karen's directory
  cd /project/berglandlab/Karen

  # Map reads to Mitochondria genome
  lanes=($( find . -type f -name *${fastawild}* | awk 'FS="_" {print $2}' | sort | uniq ))

  # Forloop through different lanes
  for lane in ${lanes[@]}; do

    # Progress lane
    echo "${lane}"

    # Identify forward and reverse
    fast1=($( find . -type f -name *${lane}*_1_*${fastawild}* ))
    fast2=($( find . -type f -name *${lane}*_2_*${fastawild}* ))

    # Trim illumina adaptors
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE \
    -threads ${threads} \
    ${fast1} \
    ${fast2} \
    ${outfq}/${samp}_${lane}_1.P.trimm.fastq.gz \
    ${outfq}/${samp}_${lane}_1.U.trimm.fastq.gz \
    ${outfq}/${samp}_${lane}_2.P.trimm.fastq.gz \
    ${outfq}/${samp}_${lane}_2.U.trimm.fastq.gz \
    ILLUMINACLIP:${wd}/data/CombinedPE-PE.fa:2:30:10:8:true

    # Merge overlapping reads
    /scratch/csm6hg/pear \
    -f ${outfq}/${samp}_${lane}_1.P.trimm.fastq.gz \
    -r ${outfq}/${samp}_${lane}_2.P.trimm.fastq.gz \
    -o ${outfq}/${samp}_${lane} \
    -j ${threads}

    # Message
    echo "Mapping now"

    # Map to mitochondrial genome
    bwa mem -t ${threads} -K 100000000 -Y \
    ${ref} \
    ${outfq}/${samp}_${lane}.assembled.fastq | \
    samtools view -Suh -q 20 -F 0x100 | \
    samtools sort --threads ${threads} -o ${outfq}/${samp}_${lane}.sort.bam
    samtools index ${outfq}/${samp}_${lane}.sort.bam

    # Unassembled reads
    bwa mem -t ${threads} -K 100000000 -Y \
    ${ref} \
    ${outfq}/${samp}_${lane}.unassembled.forward.fastq \
    ${outfq}/${samp}_${lane}.unassembled.reverse.fastq | \
    samtools view -Suh -q 20 -F 0x100 | \
    samtools sort --threads ${threads} -o ${outfq}/${samp}_${lane}.filt.unassembled.sort.bam
    samtools index ${outfq}/${samp}_${lane}.filt.unassembled.sort.bam

    # Merge assembled and unassembled bam files and mark duplicates
    samtools merge ${outfq}/${samp}_${lane}.filt.merged.bam \
    ${outfq}/${samp}_${lane}.sort.bam \
    ${outfq}/${samp}_${lane}.filt.unassembled.sort.bam

    # Index merged bam
    samtools index ${outfq}/${samp}_${lane}.filt.merged.bam

    echo "Finish" ${lane}

  # Finish lanes
  done

  ### next, merge bam files to single bam file
  samtools merge ${outfq}/${samp}_finalmap.merged.bam \
  ${outfq}/${samp}_*.filt.merged.bam
  samtools index ${outfq}/${samp}_finalmap.merged.bam

  # Mark duplicates
  java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  REMOVE_DUPLICATES=true \
  INPUT=${outfq}/${samp}_finalmap.merged.bam \
  OUTPUT=${outfq}/${samp}_finalmap.bam \
  METRICS_FILE=${outfq}/${samp}_finalmap_mdups.metrics \
  CREATE_INDEX=true

  # Call snps
  bcftools mpileup -q 30 -Q 20 -Ou -f ${ref} \
  ${outfq}/${samp}_finalmap.bam | \
  bcftools call -mv -V indels | bcftools filter \
  -i '%QUAL>20 & DP>=100' -Oz -o ${outfq}/${samp}.filt.mito.vcf.gz

  # Index VCF
  tabix -p vcf ${outfq}/${samp}.filt.mito.vcf.gz

  # apply variants to create consensus sequence
  cat ${ref} | \
  bcftools consensus ${outfq}/${samp}.filt.mito.vcf.gz \
  --sample ${outfq}/${samp}_finalmap.bam > \
  ${outfq}/${samp}.filt.consensus.mito.fa

  # Chnage name of header
  name=$( echo ${samp}.filt.consensus.mito.fa | tr "." '\t' | cut -f1 | tr "/" '\t' | cut -f1 )
  echo ${name}

  # Extracts chromosome via bed file
  sed -i "s/^>*"mtdna"/>"mtdna."${name}/g" \
  ${outfq}/${samp}.filt.consensus.mito.fa

  # Index fasta
  samtools faidx ${outfq}/${samp}.filt.consensus.mito.fa

  # Move output
  mv ${outfq}/${samp}_finalmap* ${outbam}/
  mv ${outfq}/${samp}*mito.fa* ${outfasta}/
  mv ${outfq}/${samp}*mito.vcf* ${outvcf}/

  # Remove previous steps
  rm ${outfq}/${samp}*trimm*
  rm ${outfq}/${samp}*unassembled*
  rm ${outfq}/${samp}*assembled*
  rm ${outfq}/${samp}*discarded*
  rm ${outfq}/${samp}*merged*
  rm ${outfq}/${samp}*sort.bam*

else
  echo "Sample not in Euro dataset"
fi

else
    echo "${samp} exists."
fi

# Finish
echo "Finish"
