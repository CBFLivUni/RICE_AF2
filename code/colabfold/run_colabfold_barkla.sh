module purge

# load tensorflow and relevant modules
module load apps/anaconda3/5.2.0

source activate base  # activate base first as sometimes conda activate doesn't work
conda activate my_gpu

# SET LD_LIBRARY_PATH TO THE CORRECT PATH FOR FINDING THE CORRECT GCC INSTALL
# ADDED TO THE START OF THE LIBRARY PATH OTHERWISE USES DIFFERENT DEFAULT PATH AND WON'T WORK CONSISTENTLY
# SHOULDN'T BE NECESSARY, BUT THIS IS THE BEST WAY I'VE FOUND SO FAR TO MAKE SURE IT USES THE RIGHT GCC VERSION AND DOESN'T FAIL
export LD_LIBRARY_PATH=/mnt/data1/users/aroths/.conda/envs/my_gpu/lib:$LD_LIBRARY_PATH

# check GPU visible
echo "CUDA_VISIBLE_DEVICES : $CUDA_VISIBLE_DEVICES"
echo "GPU_DEVICE_ORDINAL   : $GPU_DEVICE_ORDINAL"

# array of file names from chromosome 7
#fasta_names=($(ls main/fasta_output_unique/*7_*))

fasta_path="main/fastas_for_colabfold/priority_1/"
fasta_names=($(ls $fasta_path))

# add an idx offset to make sure each job looks at different fastas, last job should mop up all remaining sequences
# if sequence already written is automatically skipped as .done.txt file is created

for f in "${!fasta_names[@]}";
do
	echo "job_number: {$SLURM_ARRAY_TASK_ID} on ${fasta_names[$f]}"
	echo $(($f + (($SLURM_ARRAY_TASK_ID * 2))))
	colabfold_batch --num-models 3 --use-gpu-relax --amber ${fasta_path}${fasta_names[(($f + (($SLURM_ARRAY_TASK_ID * 2))))]} colabfold_output/
done

source deactivate my_gpu
