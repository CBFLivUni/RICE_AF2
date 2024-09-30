import pandas as pd
from tqdm import tqdm
import os
import glob
import json
import statistics
from time import time
import os

source_df = pd.read_csv(os.path.join("..", "..", "data", "source_unique_df_preprocessed_w_prior1.csv"))

results_chr7_df = source_df[source_df['maj_vot_chr'] == 7][['RAP-DB - official', 'MSU - official', 'OsNip']]
results_chr7_df['idx'] = results_chr7_df.index

idx_to_get = results_chr7_df.index.tolist()

# colabfold results
colabfold_results_list = []

colabfold_results_path = 'colabfold_output'

# get path for colabfold results, faster to get first
idx_to_colabfold_path = {}

print('Getting colabfold paths')

best_model_paths = glob.glob(os.path.join("..", "..", "data", colabfold_results_path, '*_scores_rank_001_*.json'))

for path in tqdm(best_model_paths):
    base_path = os.path.basename(path)
    idx_to_colabfold_path[base_path.split("_")[0]] = base_path

# get colabfold results
print('Getting colabfold results')

# colabfold results
colabfold_results_list = []

for idx in tqdm(idx_to_get):
    
	results_dict = {}
    
    # get the path of the best ranking model
	try:
		f = idx_to_colabfold_path[str(idx)]
	except KeyError:
		print(f"Index {idx} not found in colabfold results")
		continue
 
	with open(os.path.join(colabfold_results_path, f)) as json_file:
		json_dict = json.load(json_file)

	results_dict['ptm'] = json_dict['ptm'] 
	results_dict['mean_plddt'] = statistics.mean(json_dict['plddt'])  # faster than np
 
	# count plddt greater than x
	results_dict['plddt_50'] = sum(i > 50 for i in json_dict['plddt'])
	results_dict['plddt_60'] = sum(i > 60 for i in json_dict['plddt'])
	results_dict['plddt_70'] = sum(i > 70 for i in json_dict['plddt'])
	results_dict['plddt_80'] = sum(i > 80 for i in json_dict['plddt'])
	results_dict['plddt_90'] = sum(i > 90 for i in json_dict['plddt'])

	results_dict['idx'] = idx

	colabfold_results_list.append(results_dict)

colabfold_results_df = pd.DataFrame(colabfold_results_list)

# get foldseek results
foldseek_results_list = []

foldseek_results_path = 'foldseek_output'

print('Getting foldseek results')
for idx in tqdm(idx_to_get):
    
    # if file empty, or not file, then return 'NA' and continue
    
	df_cols = ["query","target","alntmscore","lddt","fident","alnlen",
			"mismatch","gapopen","qstart","qend","tstart","tend","evalue",
			"bits","taxid","taxname","taxlineage"]
 
	results_dict = {}
	results_dict['idx'] = idx

	try:
		if not os.path.isfile(os.path.join("..", "..", "data", foldseek_results_path, f"{idx}.m8")):
			raise pd.errors.EmptyDataError
		df = pd.read_csv(os.path.join('foldseek_output', f"{idx}.m8"), sep='\t', header=None)
	except pd.errors.EmptyDataError:

		for col in df_cols:
			results_dict[col] = 'NA'
		foldseek_results_list.append(results_dict)

		continue

	# get columns from foldseek input
	df.columns = df_cols

	best_hit_df = df[df['evalue'] == df['evalue'].min()]

	# loop over cols and add vals to dict
	for col in df_cols:
		results_dict[col] = best_hit_df[col].values[0]

	foldseek_results_list.append(results_dict)

foldseek_results_df = pd.DataFrame(foldseek_results_list)

# then match on idx
results_df = results_chr7_df.merge(colabfold_results_df, on='idx').merge(foldseek_results_df, on='idx')

results_df.to_csv('results_chr7.csv', index=False)
