#!/usr/bin/env bash
#
#SBATCH -J run_dsuite # A single job name for the array
#SBATCH --ntasks-per-node=1 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 0-05:00 # 5 hours
#SBATCH --mem 2G
#SBATCH -o /project/berglandlab/connor/dsuite4/err/dsuite.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/dsuite4/err/dsuite.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# Load modules
module load gcc

### PARAMETERS FOR INPUT AND OUTPUT ###

# Working directory
wd="/project/berglandlab/connor/dsuite4"

# Cluster file
clust="pop.dsuite.country.ameuropul.clust"
clust2="pop.dsuite.country.europul.clust"

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

# Tree 1
#(((Daphnia.pulex.NorthAmerica,Daphnia.pulexcaria.NorthAmerica),Daphnia.pulicaria.NorthAmerica),Outgroup)

# Tree 2
#(((Daphnia.pulex.NorthAmerica,Daphnia.pulicaria.Europe),Daphnia.pulicaria.NorthAmerica),Outgroup)

# Tree 3
#(((Daphnia.pulex.NorthAmerica,Daphnia.pulicaria.NorthAmerica),Daphnia.pulicaria.Europe),Outgroup)

# Tree 4
#(((Daphnia.pulex.NorthAmerica,Daphnia.pulicaria.NorthAmerica),Daphnia.pulex.Europe),Outgroup)

# Dsuite - with continent and species as populations Euro Pulicaria
~/dsuite Dtrios \
-t ${wd}/tree3.nwk \
-r ${start},${len} \
-o ${wd}/data4/pulicaria_euro_P3.${start}.${stop} \
/project/berglandlab/connor/new_vcf2/daphnia.filt.qual.miss.rep.ann.vcf.gz \
${wd}/${clust}

# Dsuite - with continent and species as populations Euro Pulicaria
~/dsuite Dtrios \
-t ${wd}/tree4.nwk \
-r ${start},${len} \
-o ${wd}/data4/pulex_euro_P3.${start}.${stop} \
/project/berglandlab/connor/new_vcf2/daphnia.filt.qual.miss.rep.ann.vcf.gz \
${wd}/${clust2}

# Finish
echo "Finish" ${start} "-" ${stop}
