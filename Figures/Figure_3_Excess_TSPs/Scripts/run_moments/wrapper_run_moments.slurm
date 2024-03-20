#!/usr/bin/env bash
#
#SBATCH -J moments_boot # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-05:00 # 5 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/csm6hg/err/moments.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/err/moments.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# This Slurm script runs an array job of moments runs with one job per line of instructions_run_moments.txt,
# which should contain one SFS projection size per line, e.g. 20 corresponds to 20x20 SFS projections in moments

# Options file should contain output of "ls [input SFS directory]": SLURM_ARRAY_TASK_ID=283
options_file="/scratch/csm6hg/moments/instructions.metadata.txt"
script_file="/scratch/csm6hg/moments/2.run_moments.py"

# Feed parameters
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file} | cut -f1 )

# Progress + quality of life
echo "SLURM Job:" ${SLURM_ARRAY_TASK_ID}
echo "SFS Projection size for run:" ${OPTS}

# Run moments script
python ${script_file} \
${OPTS} \
${SLURM_ARRAY_TASK_ID}
