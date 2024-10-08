import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

# list of priority 2 and 3 genes (+ all proteins assigned to those genes)
# with check that they've not been in the priority 1 list

# priority 1 file doesn't match to the original priority 1 file i used to generate the fastas, so I've used that instead
#priority_1_df = pd.read_csv(os.path.join('..', '..', 'data', 'raw', 'Nipponbare_priority_1.tab'), sep='\t')
priority_1_df = pd.read_csv(os.path.join('..', '..', 'data', 'raw', 'Nipponbare_new_filter_Priority_1.tab'), sep='\t')
priority_2_df = pd.read_csv(os.path.join('..', '..', 'data', 'raw', 'Nipponbare_priority_2.tab'), sep='\t')
priority_3_df = pd.read_csv(os.path.join('..', '..', 'data', 'raw', 'Nipponbare_priority_3.tab'), sep='\t')
priority_4_df = pd.read_csv(os.path.join('..', '..', 'data', 'raw', 'Nipponbare_priority_4.tab'), sep='\t')

priority_1_trans = priority_1_df['child_transcripts'].tolist()
priority_2_trans = priority_2_df['transcript'].tolist()
priority_3_trans = priority_3_df['transcript'].tolist()
priority_4_trans = priority_4_df['transcript'].tolist()

# get all sequences for those transcripts
source_df = pd.read_csv(os.path.join('..', '..', 'data', 'source_unique_df_preprocessed_w_prior1.csv'))

# drop predicted
source_df = source_df.drop(columns=['RAP-DB - predicted', 'Ensembl - predicted'])

def get_sequences(source_df, priority_list):
	
	# k: transcript, v: seq
	priority_seqs = {}
 
	# check every column for a match.
	id_cols = ['RAP-DB - official',
		'Ensembl - official', 'MSU - official',
		'MSU - TE', 'OsNip', 'UniProt']

	for i, row in tqdm(source_df.iterrows(), total=source_df.shape[0]):
     
		# turn row into string
		source_string = '|'.join(row[id_cols].values.astype(str))
    
		# search row for sequence, if found, append seq to list
		for id in priority_list:
			if id in source_string:
				priority_seqs[id] = row['seq']
	
	return priority_seqs

# dicts of transcript: seq
priority_1_seqs = get_sequences(source_df, priority_1_trans)
priority_2_seqs = get_sequences(source_df, priority_2_trans)
priority_3_seqs = get_sequences(source_df, priority_3_trans)
priority_4_seqs = get_sequences(source_df, priority_4_trans)

# then check sequences against each other

priority_2_not_inc_1 = {}
priority_3_not_inc_1_2 = {}

# check not in 2 for 1
priority_1_set_seqs = list(set(priority_1_seqs.values()))

for k, v in tqdm(priority_2_seqs.items()):
	if v not in priority_1_set_seqs:
		priority_2_not_inc_1[k] = v

print(len(priority_2_seqs))
print(len(priority_2_not_inc_1))

# check not in 3 for 1 and 2
priority_2_set_seqs = list(set(priority_2_seqs.values()))

for k, v in tqdm(priority_3_seqs.items()):
	if v not in priority_1_set_seqs and v not in priority_2_set_seqs:
		priority_3_not_inc_1_2[k] = v

print(len(priority_3_seqs))
print(len(priority_3_not_inc_1_2))

# then filter out eshan's lists based on transcripts
priority_2_filtered_df = priority_2_df[priority_2_df['transcript'].isin(priority_2_not_inc_1.keys())]
priority_3_filtered_df = priority_3_df[priority_3_df['transcript'].isin(priority_3_not_inc_1_2.keys())]

# save to file
priority_2_filtered_df.to_csv(os.path.join('..', '..', 'data', 'raw', 'Nipponbare_priority_2_filtered_out_prio1.tab'), sep='\t', index=False)
priority_3_filtered_df.to_csv(os.path.join('..', '..', 'data', 'raw', 'Nipponbare_priority_3_filtered_out_prio1_2.tab'), sep='\t', index=False)
