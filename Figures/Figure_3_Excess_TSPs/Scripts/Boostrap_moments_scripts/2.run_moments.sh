#!/usr/bin/env bash
#
#SBATCH -J moments_boot # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-00:30 # 30 mins
#SBATCH --mem 2G
#SBATCH -o /scratch/csm6hg/err/moments.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/err/moments.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# This Slurm script runs an array job of moments runs with one job per line of instructions_run_moments.txt,
# which should contain one SFS projection size per line, e.g. 20 corresponds to 20x20 SFS projections in moments

# Working directory: SLURM_ARRAY_TASK_ID=5
wd="/scratch/csm6hg/phd/moments/"
cd ${wd}

# Modules
module load anaconda/2023.07-py3.11
conda activate msprime_env

# Options file should contain output of "ls [input SFS directory]"
options_file="instructions.metadata.txt"
script_file="run_moments.py"

# Feed parameters
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file} | cut -f1)
modei=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file} | cut -f3)

# Progress + quality of life
echo "SLURM Job:" ${SLURM_ARRAY_TASK_ID}
echo "SFS Projection size for run:" ${OPTS}
echo "Model:" ${modei}

# Run moments script
python \
${script_file} \
${OPTS} \
${SLURM_ARRAY_TASK_ID} \
${modei}
