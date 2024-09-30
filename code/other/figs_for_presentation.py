import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import os

results_df = pd.read_csv('results_chr7.csv')

# note if present in RAP-DB, MSU, OsNip
present_in = []

for idx, row in tqdm(results_df.iterrows(), total=results_df.shape[0]):
	
	present_for_row = []

	# check if exists
	if pd.notna(row['RAP-DB - official']):
		present_for_row.append('RAP-DB')
	
	if pd.notna(row['MSU - official']):
		present_for_row.append('MSU')
	
	if pd.notna(row['OsNip']):
		present_for_row.append('OsNip')
  
	if len(present_for_row) == 0:
		present_for_row = ['None']
 
	present_in.append(' + '.join(present_for_row))
 
results_df['present_in'] = present_in

# violin plots for each column based on present_in

# set order of columns for violin plots
results_df['present_in'] = pd.Categorical(results_df['present_in'],
                            categories=['RAP-DB', 'MSU', 'OsNip', 'RAP-DB + MSU', 'RAP-DB + OsNip', 'MSU + OsNip', 'RAP-DB + MSU + OsNip', 'None'],
                            ordered=True)

for col in results_df.columns[4:len(results_df.columns)-1]:
	if col in ['query', 'target', 'taxid', 'taxname', 'taxlineage']:
		continue

	plt.figure(figsize=(12, 8))
	sns.violinplot(x='present_in', y=col, data=results_df, cut=0)

	plt.title(col) 
	plt.gca().spines['top'].set_visible(False)
	plt.gca().spines['right'].set_visible(False)
	plt.xticks(rotation=45)
	plt.ylim(0, None)
	plt.tight_layout()
	plt.xlabel('Present in')
 
	plt.savefig(os.path.join('plots', f'{col}.png'))
	plt.close()

# ptm score group by present_in group by agg
results_df['ptm'].groupby(results_df['present_in']).agg(['mean', 'std', 'min', 'max', 'median', 'count'])
