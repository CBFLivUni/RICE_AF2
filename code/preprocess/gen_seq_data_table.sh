#!/bin/bash -l

# Use the current environment for this job
#SBATCH --export=ALL
#SBATCH --job-name gen_seq_data

#SBATCH -o gen_seq_data.out
#SBATCH -p nodes
#SBATCH -N 1

#SBATCH -n 6
#SBATCH -t 2-00:00:00

#SBATCH --mail-user=aroths@liverpool.ac.uk
#SBATCH --mail-type=ALL

source af2venv/bin/activate

python gen_seq_data_table.py