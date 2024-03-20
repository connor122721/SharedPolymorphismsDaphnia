#!/usr/bin/env bash
#SBATCH -J ext.Exon # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-05:00 # 4 hours
#SBATCH --mem 5G
#SBATCH -o /project/berglandlab/connor/err/ext.Exon.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/ext.Exon.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Extract chromosome from each individual
getRegionInd () {

### load modules
module load samtools
module load bcftools

# Working directory
wd="/project/berglandlab/connor"

# Parameter file
# 1st column = sample
parameterFile="${wd}/candgene/samples_phase_var1"

# Bed file
chromFile="${wd}/candgene/exons.genome.list.tspset6"

# Job id information
samp=$( sed "${iteration}q;d" ${parameterFile} | cut -f1 )
pop=$( sed "${iteration}q;d" ${parameterFile} | cut -f2 )

# Reference genome
ref="${wd}/totalHiCwithallbestgapclosed.fa"

# Message for sample
echo "Sample:" ${samp} "id:" ${iteration} "Pop:" ${pop}
echo "${iteration}"

# Runs function for each haplotype
getHaplos () {
whichHaplo=${1}
samp=${2}
chrom=${3}
#whichHaplo=2;chrom=Scaffold_7627_HRSCAF_8561:1-943

# Message
echo ${samp} ${whichHaplo} ${chrom} ${pop}

# Parameters
chr=$( echo ${chrom} | tr ":" '\t' | cut -f1 -d"," | cut -f1 )
start=$( echo ${chrom} | tr ":" '\t' | cut -f1 -d"," | tr "-" '\t' | cut -f2 )
stop=$( echo ${chrom} | tr ":" '\t' | cut -f1 -d"," | tr "-" '\t' | cut -f3 )

# Message for chromosome
echo "Chromosome:" ${chr} ${start} ${stop}

# Creates a consensus fasta file for each haplotype
bcftools \
consensus \
-f <(samtools faidx ${ref} $( echo $chrom | cut -f1 -d"," )) \
-H ${whichHaplo} \
-s ${samp} \
-o ${wd}/candgene/fasta_exon/${samp}.${whichHaplo}.${chr}.${start}.${stop}.fa \
${wd}/new_phase/daphnia.whatshap.ann.vcf.gz

# Add sample name and haplotype to header
f="${pop}_${samp}.${whichHaplo}.${chr}.${start}.${stop}"

# Extracts chromosome via bed file
sed -i "s/${chrom}/${f}/g" \
${wd}/candgene/fasta_exon/${samp}.${whichHaplo}.${chr}.${start}.${stop}.fa

# Show header of new FASTA
head -n 1 ${wd}/candgene/fasta_exon/${samp}.${whichHaplo}.${chr}.${start}.${stop}.fa

}

# Exports function
export -f getHaplos

# Parallelize tasks
parChr() {

# Runs "Extract haplotype" function
getHaplos 1 ${samp} ${i}
getHaplos 2 ${samp} ${i}

chr=$( echo ${chrom} | tr ":" '\t' | cut -f1 )
start=$( echo ${i} | tr ":" '\t' | cut -f1 -d"," | tr "-" '\t' | cut -f1 )
stop=$( echo ${i} | tr ":" '\t' | cut -f1 -d"," | tr "-" '\t' | cut -f2 )

}

# For loop through each window
i=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${chromFile} | cut -f1 -d"," )
parChr "${i}"

}

# Export individual function
export getRegionInd

# For loop through each individual
for iteration in {1..62}; do getRegionInd "${iteration}"; done

# Finish script
echo "Finish fasta extraction"

# Load modules
module purge
module load gcc/9.2.0 openmpi/3.1.6 mafft/7.475

# Working directory
wd="/project/berglandlab/connor"

# Bed file
chromFile="/project/berglandlab/connor/candgene/exons.genome.list.tspset6"

# Go through each gene through each window
i=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${chromFile} | cut -f1 -d"," )

# Parameters
chr=$( echo ${i} | tr ":" '\t' | cut -f1 )
start=$( echo ${i} | tr ":" '\t' | cut -f2 | tr "-" '\t' | cut -f1)
stop=$( echo ${i} | tr ":" '\t' | cut -f2 | tr "-" '\t' | cut -f2 )
gene=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${chromFile} | cut -f2 -d"," )
class=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${chromFile} | cut -f3 -d"," )

# Progress message
echo "Gene:" $gene
echo "Classification:" $class
echo "Coordinates:" $chr $start $stop

# Add sample name and haplotype to header
f="${chr}.${start}.${stop}.${gene}.${class}"

# European Daphnia pulex reference
ref_europul="/project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa"

# Euro D. pulex gff
gff="/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gff"

# Finding fastas
f=${chr}.${start}.${stop}

# Concat every fasta
cat ${wd}/candgene/fasta_exon/*${f}.*fa > \
${wd}/candgene/tree_exon_tsp_250/${f}.fa

# Align fasta
mafft --auto \
${wd}/candgene/tree_exon_tsp_250/${f}.fa > \
${wd}/candgene/tree_exon_tsp_250/${f}.aln.fa

# Run iqtree2
~/iqtree2 \
-redo \
-bb 1000 \
-s ${wd}/candgene/tree_exon_tsp_250/${f}.aln.fa \
-T 1 \
--prefix ${wd}/candgene/tree_exon_tsp_250/${f}

# Remove temp files
rm ${wd}/candgene/tree_exon_tsp_250/${f}.splits.nex
rm ${wd}/candgene/tree_exon_tsp_250/${f}.model.gz
rm ${wd}/candgene/tree_exon_tsp_250/${f}.varsites.phy
rm ${wd}/candgene/tree_exon_tsp_250/${f}.iqtree
rm ${wd}/candgene/tree_exon_tsp_250/${f}.bionj
rm ${wd}/candgene/tree_exon_tsp_250/${f}.contree
rm ${wd}/candgene/tree_exon_tsp_250/${f}.mldist
rm ${wd}/candgene/tree_exon_tsp_250/${f}.ufboot
rm ${wd}/candgene/tree_exon_tsp_250/${f}.ckp.gz
rm ${wd}/candgene/tree_exon_tsp_250/${f}.log
rm ${wd}/candgene/tree_exon_tsp_250/${f}.fa
rm ${wd}/candgene/fasta_exon_250/*${f}.*fa*

# Finish tree
echo "Finished tree"
