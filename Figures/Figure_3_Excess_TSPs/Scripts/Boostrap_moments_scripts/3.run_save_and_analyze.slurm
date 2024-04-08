#!/usr/bin/env bash
#
#SBATCH -J moments_boot # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00 # 1 hour
#SBATCH --mem 10G
#SBATCH -o /scratch/csm6hg/err/moments.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/err/moments.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# This Slurm script runs an array job of moments runs with one job per line of instructions_run_moments.txt,
# which should contain one SFS projection size per line, e.g. 20 corresponds to 20x20 SFS projections in moments

# Working directory: SLURM_ARRAY_TASK_ID=1748
cd /scratch/csm6hg/moments/

# Modules
module purge
module load anaconda
DIR=/home/$USER/.conda/envs/msprime_env
conda activate msprime_env
export PATH=$DIR/bin:$PATH
export LD_LIBRARY_PATH=$DIR/lib:$PATH
export PYTHONPATH=$DIR/lib/python3.8/site-packages:$PATH

# Options file should contain output of "ls [input SFS directory]"
options_file="/scratch/csm6hg/moments/optimizedRuns_moments.txt"
script_file="/scratch/csm6hg/moments/3.save_and_analyze_sfss.py"

# Feed parameters
modei=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file} | cut -f2)

# Conditional make parameter
if [ ${modei} == "split_mig" ]
then
    OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file} | \
    cut -f3-6 | sed 's/\t/,/g' )
elif [ ${modei} == "split_no_mig" ]
then
    OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file} | \
    cut -f3-5 | sed 's/\t/,/g' ),0
fi

# Progress + quality of life
echo "SLURM Job:" ${SLURM_ARRAY_TASK_ID}
echo "SFS Parameters:" ${OPTS}
echo "Model:" ${modei}

# Run moments script
python \
${script_file} \
${modei} \
${OPTS} \
${SLURM_ARRAY_TASK_ID}
