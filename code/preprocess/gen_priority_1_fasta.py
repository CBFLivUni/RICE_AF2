import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

priority_df = pd.read_csv(os.path.join('..', '..', 'data', 'raw', 'Nipponbare_new_filter_Priority_1.tab'), sep='\t')

# child transcript is the only column I can match on
priority_list = priority_df['child_transcripts'].tolist()

# master source df
source_df = pd.read_csv(os.path.join('..', '..', 'data', 'source_unique_df_preprocessed.csv'))

# add marker if in priority list
# check every column for a match.
id_cols = ['RAP-DB - official', 'RAP-DB - predicted',
       'Ensembl - official', 'Ensembl - predicted', 'MSU - official',
       'MSU - TE', 'OsNip', 'UniProt'] 

source_string = source_df[id_cols].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)

priority_1_col = []

# join the id's and search for priority sequences in the source df
for row in tqdm(source_string):
    
    # set default values to false
	is_priority = 'F'
    
	for id in priority_list:
		if id in row:
			is_priority = 'T'
			break
	
	priority_1_col.append(is_priority)
  
source_df['priority_1'] = priority_1_col
source_df.to_csv(os.path.join('..', '..', 'data', 'source_unique_df_preprocessed_w_prior1.csv'), index=False)

priority_1_df = source_df[source_df['priority_1'] == 'T']
priority_1_df['idx'] = priority_1_df.index

JOB_NAME = "priority_1"

def save_fasta(current_chr, fasta_chunk, records_for_fasta):
    
    # chr changed to reorder, reset to original
	if current_chr == 999:
		current_chr = 'UNKNOWN'
	elif current_chr == 998:
		current_chr = 'ERROR'
    
	fasta_name = os.path.join('..', '..', 'data', 'fastas_for_colabfold', JOB_NAME, f'{current_chr}_{fasta_chunk}.fasta')
	with open(fasta_name, "w") as output_handle:
		SeqIO.write(records_for_fasta, output_handle, "fasta")

# loop over rows in table

# initial values
fasta_chunk = 1  # n fasta has been chunked per chr
current_chr = 1 # current chromosome being used
count_per_fasta = 0  # records in fasta
records_for_fasta = []  # records to be written to fasta

for i in tqdm(range(len(priority_1_df)), total=len(priority_1_df)):
    
    # itterows breaks because idx has been changed
	row = priority_1_df.iloc[i]

	if i == 0:
		print(i)
		print(row['seq'])
    
	new_chr = row['maj_vot_chr']
    
	if new_chr == current_chr:

		one_sequence = Seq(data=row['seq'])
		one_record = SeqRecord(one_sequence, id=str(row['idx']), description="")
		records_for_fasta.append(one_record)
  
		count_per_fasta += 1

		# max 20 seqs per fasta file
		if count_per_fasta == 20:
			
			save_fasta(current_chr, fasta_chunk, records_for_fasta)

			fasta_chunk += 1

			# reset the records going into the fasta
			current_chr = new_chr
			count_per_fasta = 0
			records_for_fasta = []

	else:
		# must be new chromosome
		# check there are records to save, then save final fasta for that chr
  
		if len(records_for_fasta) != 0:
			save_fasta(current_chr, fasta_chunk, records_for_fasta)

		# reset chunk for the new chromosome
		fasta_chunk = 1
  
		# add the first record for the new chromosome
		records_for_fasta = []
		one_sequence = Seq(data=row['seq'])
		one_record = SeqRecord(one_sequence, id=str(row['idx']), description="")
		records_for_fasta.append(one_record)

		# reset the records going into the fasta
		current_chr = new_chr
		count_per_fasta = 1

# save final record
save_fasta(current_chr, fasta_chunk, records_for_fasta)
