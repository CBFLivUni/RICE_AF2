#!/bin/bash -l
#SBATCH -J foldseek
#SBATCH -o foldseek.out
#SBATCH -p low
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 3-00:00:00
#SBATCH --mail-user=aroths@liverpool.ac.uk
# Notify user by email when certain event types occur
#SBATCH --mail-type=ALL

# use relaxed models as generally more accurate, and best ranked
model_names=($(ls ../../data/colabfold_output/*_relaxed_rank_001_*.pdb))

for m in "${!model_names[@]}";
do
  # get id of model before "_"
  id=$(basename "${model_names[$m]}" | cut -d'_' -f1)

  echo "Querying ${id}"
  echo ${model_names[$m]}
  foldseek easy-search ${model_names[$m]} `# input structure` \
  ../../foldseek_databases/pdb `# db to query` \
  ../../data/foldseek_output/"$id".m8 tmp `# output file` \
  --format-output "query,target,alntmscore,lddt,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage" \
  --exhaustive-search \
  --threads $SLURM_NTASKS

done