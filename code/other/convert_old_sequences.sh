#!/bin/bash -l

# Use the current environment for this job
#SBATCH --export=ALL
#SBATCH --job-name gen_convert_old_sequences

#SBATCH -o gen_convert_old_sequences.out
#SBATCH -p nodes
#SBATCH -N 1

#SBATCH -n 8
#SBATCH -t 2-00:00:00

#SBATCH --mail-user=aroths@liverpool.ac.uk
#SBATCH --mail-type=ALL

module purge
module load apps/anaconda3/2023.03-poetry
source /mnt/data1/users/software/anaconda/anaconda3-2023.03/etc/profile.d/conda.sh

conda activate base
conda activate AF_sequence_data_env

python convert_old_sequences.py