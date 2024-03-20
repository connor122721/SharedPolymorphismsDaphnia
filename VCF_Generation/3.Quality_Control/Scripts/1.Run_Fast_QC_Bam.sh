#!/usr/bin/env bash
#SBATCH -J fastqc # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -t 1-00:00 # Running time of 1 day
#SBATCH --mem 20G # Memory request of 20 GB
#SBATCH -o /scratch/csm6hg/daphnia_phylo/fastqc.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/fastqc.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load fastqc
module load singularity

# Parameters
threads=20

# Bam locations
wd="/scratch/csm6hg/daphnia_phylo/all_bam/"

# Output QC statistics
out="/scratch/csm6hg/daphnia_phylo/mapping_stats/fastqc_output"

# Create temporary folders
if [[ -d "${out}" ]]
then
	echo "Working tmp folder exist"
	echo "lets move on"
	date
else
	echo "Folder doesnt exist. Let us fix that."
	mkdir ${out}
	date
fi

# For loop for every bam - run fastqc
for f in $wd**/*bam; do
#f="/scratch/csm6hg/daphnia_phylo/all_bam/Euro_bams/April_2017_DBunk_147_finalmap_mdup.bam"

echo "Processing:" ${f}
file=$( echo ${f} | sed 's/....$//' )

# Run samtools
fastqc ${f} -t ${threads}

# Now consolidate all files into a single folder
find ${file}*fastqc* -exec mv {} ${out} \;

# Finish loop
done

# Message for transition into multiqc
echo "Finished running fastqc"
date
echo "Now consolidating fastqc output using multiqc"

# Move to fastqc output folder
cd ${out}

# Now run multi QC
singularity run /home/csm6hg/multiqc.sif ./

echo "Completed running multiqc"
date
