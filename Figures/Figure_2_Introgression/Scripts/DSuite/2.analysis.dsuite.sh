#!/usr/bin/env bash
#
#SBATCH -J run_dsuite # A single job name for the array
#SBATCH --ntasks-per-node=1 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 0-05:00 # 5 hours
#SBATCH --mem 2G
#SBATCH -o /project/berglandlab/connor/backup_project/chapter1/dsuite4/err/dsuite.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/backup_project/chapter1/dsuite4/err/dsuite.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules: SLURM_ARRAY_TASK_ID=2
module load gcc

### PARAMETERS FOR INPUT AND OUTPUT ###

# Working directory
wd="/project/berglandlab/connor/backup_project/chapter1/dsuite4/"

# Cluster file
clust="pop.dsuite.country.ameuropul.clust"
clust2="pop.dsuite.country.europul.clust"
clust3="pop.ameuropul.clust"

# Parameter file
intervals="interval_dsuite_paramList"

# Cluster population file
echo "Array job:" ${SLURM_ARRAY_TASK_ID}

# Start of chromosome window (position)
start=$( cat ${wd}/${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )

# Stop of chromosome window (position)
stop=$( cat ${wd}/${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )

# Length of chromsome window
len=$( cat ${wd}/${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )

# Progress message
echo "Start:" ${start} "Stop:" ${stop} "Width:" ${len}

# Tree 1 - Euro pulic
#(((Daphnia.pulex.NorthAmerica,Daphnia.pulex.Europe),Daphnia.pulicaria.Europe),Outgroup)

# Tree 2 - Euro pulic
#(((Daphnia.pulex.NorthAmerica,Daphnia.pulicaria.Europe),Daphnia.pulex.Europe),Outgroup)

# Dsuite - with continent and species as populations Euro Pulicaria
~/dsuite Dtrios \
-t ${wd}/tree1_europulic.nwk \
-r ${start},${len} \
-o ${wd}/data5/pulicaria_euro_P3.${start}.${stop} \
/project/berglandlab/connor/backup_project/new_vcf2/daphnia.filt.mlg.genome.11.18.22.vcf.gz \
${wd}/${clust3}

# Dsuite - with continent and species as populations Euro Pulicaria
~/dsuite Dtrios \
-t ${wd}/tree2_europulic.nwk \
-r ${start},${len} \
-o ${wd}/data5/pulex_euro_P3.${start}.${stop} \
/project/berglandlab/connor/backup_project/new_vcf2/daphnia.filt.mlg.genome.11.18.22.vcf.gz \
${wd}/${clust3}

# Finish
echo "Finish" ${start} "-" ${stop}