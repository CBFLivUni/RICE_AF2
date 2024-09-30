import pandas as pd
import shutil
import os
from tqdm import tqdm
import json

# import old sequences
old_source_df = pd.read_csv("source_df_preprocessed.csv")
new_source_df = pd.read_csv("source_unique_df_preprocessed.csv")

old_source_df = old_source_df.drop(columns=['Unnamed: 0'])
new_source_df = new_source_df.drop(columns=['Unnamed: 0'])

# map old indexes to new indexes
# preallocate for speed
old_idx = list(range(len(old_source_df)))
new_idx = [None] * len(old_source_df)
old_seq_list = old_source_df['seq'].tolist()
new_seq_list = new_source_df['seq'].tolist()

for i in tqdm(range(len(old_seq_list))):
	old_seq = old_seq_list[i]
	for j in range(len(new_seq_list)):
		new_seq = new_seq_list[j]
		if old_seq == new_seq:
			new_idx[i] = j
			break
   
if None in new_idx:
    raise ValueError("Some sequences in old_source_df are not in new_source_df")

old_to_new = dict(zip(old_idx, new_idx))

# write to json
with open("old_to_new.json", "w") as outfile: 
    json.dump(old_to_new, outfile)


with open("old_to_new.json") as json_file:
    old_to_new = json.load(json_file)

# loop over keys. i.e. old indices.
old_files = os.listdir(os.path.join("..", "colabfold_output"))

for old_idx in tqdm(old_to_new.keys()):
    
    # needs to be checked each time, in case folder has been created
	new_files = os.listdir(os.path.join("..", "colabfold_output_new"))
 
	# new idxs with a file
	new_files_idx = [file.split("_")[0] for file in new_files]
 
	new_idx = old_to_new[old_idx]
 
	if new_idx is None:
		# if doesn't exist, then skip, but this should never happen
		continue
 
	# if new idx already has a file, then skip
	if str(new_idx) in new_files_idx:
		print(f"new_idx: {new_idx} already exists")
		continue

	# if new idx doesn't have a file, check if old_idx file has been created and needs moving
	else:
		print(f"new_idx: {new_idx} doesn't exist")
		# find all the files with that name
		# ie. by separating on "." and "_"
  
		old_files_with_old_idx = [file for file in old_files if file.split(".")[0] == str(old_idx) or file.split("_")[0] == str(old_idx)]
  
		if len(old_files_with_old_idx) == 0:
			# files not processed yet, so skip
			print(f"no files with old_idx: {old_idx}")
			continue

		else:
			# have been processed, so move and rename from old_idx to new_idx
			for old_file in old_files_with_old_idx:
				# if is a directory, then copytree
				if os.path.isdir(os.path.join("..", "colabfold_output", old_file)):
					print("copying directory:")
					print(os.path.join("..", "colabfold_output", old_file))
					print(os.path.join("..", "colabfold_output_new", old_file.replace(old_idx, str(new_idx))))
					shutil.copytree(os.path.join("..", "colabfold_output", old_file), os.path.join("..", "colabfold_output_new", old_file.replace(old_idx, str(new_idx))))
				# else is file, so copy that
				else:
					print("copying file:")
					print(os.path.join("..", "colabfold_output", old_file))
					print(os.path.join("..", "colabfold_output_new", old_file.replace(old_idx, str(new_idx))))
					shutil.copy(os.path.join("..", "colabfold_output", old_file), os.path.join("..", "colabfold_output_new", old_file.replace(old_idx, str(new_idx))))
