#!/bin/bash -l
# Define job name
#SBATCH -J get_results
# Define a standard output file. When the job is running, %u will be replaced by user name,
#SBATCH -o get_results.out
# Request the partition
#SBATCH -p low
# Request the number of nodes
#SBATCH -N 1
# Request the number of cores
#SBATCH -n 16
# Specify time limit in format a-bb:cc:dd, where a is days, b is hours, c is minutes, and d is seconds.
#SBATCH -t 3-00:00:00
# Insert your own username to get e-mail notifications (note: keep just one "#" before SBATCH)
#SBATCH --mail-user=aroths@liverpool.ac.uk
# Notify user by email when certain event types occur
#SBATCH --mail-type=ALL
# Define job array
##SBATCH --array=1-8

source af2venv/bin/activate

python get_results.py