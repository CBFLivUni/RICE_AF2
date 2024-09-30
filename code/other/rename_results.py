import os
import pandas as pd
from tqdm import tqdm

files = os.listdir('colabfold_output')
source_df = pd.read_csv("source_unique_df_preprocessed_w_prior1.csv")

priority_df = source_df[source_df['priority_1'] == 'T']

# replace the count iterating with the index, which is what actually want.
old_to_new = {}
for count, (index, row) in enumerate(priority_df.iterrows()):
    old_to_new[str(count)] = str(index)


# if not "_" in name, then split on "." and take the first part

for name in tqdm(files):
    
    # ignore "cite.bibtex", "config.json" and "log.txt"
	if name == "cite.bibtex" or name == "config.json" or name == "log.txt":
		continue

	else:
    
		if not "_" in name:
			new_name = name.replace(name.split(".")[0], old_to_new[name.split(".")[0]])
			os.rename('colabfold_output/' + name, 'colabfold_output/' + new_name)
	
		else:
			new_name = name.replace(name.split("_")[0], old_to_new[name.split("_")[0]])
			os.rename('colabfold_output/' + name, 'colabfold_output/' + new_name)
	
		print(name + " -> " + new_name)
