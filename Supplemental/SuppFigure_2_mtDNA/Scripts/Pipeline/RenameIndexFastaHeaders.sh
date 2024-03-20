#!/usr/bin/env bash
#SBATCH -J geneFasta
#SBATCH --ntasks-per-node=1 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 day
#SBATCH --mem 10G
#SBATCH -o /scratch/csm6hg/mito/err/renameFasta.out # Standard output
#SBATCH -e /scratch/csm6hg/mito/err/renameFasta.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load Modules
module load samtools

# Working directory
wd="/scratch/csm6hg/mito"

# Get samples
cd ${wd}/mito.fasta
ls *filt.consensus.mito.fa > ${wd}/samps.list

# Go through each gene
while read i; do
	# i=SRR7592710.filt.consensus.mito.fa

	echo ${i}
	cd $wd/mito.fasta

	samp=$( echo ${i} | cut -f1 -d'.')
	echo ${samp}
	# done < ${wd}/samps.list

	sed -i "s/^>NC_000844.1 Daphnia pulex mitochondrion, complete genome/>"mtdna."${samp}/g" \
	${wd}/mito.fasta/${samp}.filt.consensus.mito.fa

	sed -i "s/^>NC_026914.1*/>"mtdna."${samp}/g" \
	${wd}/mito.fasta/${samp}.filt.consensus.mito.fa

	sed -i "s/^>NC_026914.1 Daphnia magna mitochondrion, complete genome/>"mtdna."${samp}/g" \
	${wd}/mito.fasta/${samp}.filt.consensus.mito.fa

	sed -i "s/^>NC_000844.1/>"mtdna."${samp}/g" \
	${wd}/mito.fasta/${samp}.filt.consensus.mito.fa

	sed -i "s/^>CM028013.1 Daphnia obtusa isolate FS6 mitochondrion, complete sequence, whole genome shotgun sequence/>"mtdna."${samp}/g" \
	${wd}/mito.fasta/${samp}.filt.consensus.mito.fa

	sed -i "s/^>CM028013.1/>"mtdna."${samp}/g" \
	${wd}/mito.fasta/${samp}.filt.consensus.mito.fa

	samtools faidx ${wd}/mito.fasta/${samp}.filt.consensus.mito.fa

done < ${wd}/samps.list
