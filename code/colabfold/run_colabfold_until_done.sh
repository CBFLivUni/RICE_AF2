#!/bin/bash -l

#SBATCH --open-mode=append  # if job rerunning, append to existing output file
#SBATCH --export=ALL
#SBATCH --job-name colab_GPU

#SBATCH -o colabfold_gpu.out
# request first GPU parition available
#SBATCH -p gpu,gpulowbig,gpulowmed,gpulowsmall
#SBATCH -N 1
#SBATCH --gres=gpu:1
#SBATCH -n 6
#SBATCH --array=1-8  # 8 is the max GPU's that can be requested at once on barkla
#SBATCH -t 1-00:00:00

#SBATCH --mail-user=aroths@liverpool.ac.uk
#SBATCH --mail-type=ALL

# MAKE SURE STARTS WITH CONDA DEACTIVATED
module purge
for i in $(seq ${CONDA_SHLVL}); do
  conda deactivate
done

# run for 23 hours, then timeout, so not killed by SLURM
timeout 23h ./run_colabfold_barkla.sh

# if terminates due to timeout, requeue
if [[ $? == 124 ]]; then
  scontrol requeue $SLURM_JOB_ID
fi
