#!/bin/bash

fasta_path="../../data/fastas_for_colabfold/priority_1/"
fasta_names=($(ls $fasta_path))

# count all gpu's available
gpu_count=$(nvidia-smi -L | wc -l)

# write commands to be run to txt file for scheduler
for f in "${fasta_names[@]}";
do
	echo "colabfold_batch --num-models 3 --use-gpu-relax --amber ${fasta_path}${f} ../data/colabfold_output/" >> colabfold_commands.txt
done

# ensure colabfold in path
PATH="../../localcolabfold/colabfold-conda/bin:$PATH"

# run using gpu scheduler using all gpu's
simple_gpu_scheduler --gpus $(seq -s ',' 0 $(($gpu_count-1))) < colabfold_commands.txt
