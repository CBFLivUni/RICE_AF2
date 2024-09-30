import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

# table of sequences and chromosomes
source_df = pd.read_csv("source_unique_df.csv")

source_df = source_df.drop(columns=['Unnamed: 0'])

# sorted on chr
# make UNKNOWN and ERROR chromosomes that don't exist and are last
source_df['maj_vot_chr'] = source_df['maj_vot_chr'].replace({'UNKNOWN': 999})
source_df['maj_vot_chr'] = source_df['maj_vot_chr'].astype(int)

source_df = source_df.sort_values(by=['maj_vot_chr'])

# reset idx
source_df = source_df.reset_index(drop=True)

# save preprocessed source df
source_to_save = source_df.copy()
source_to_save['maj_vot_chr'] = source_to_save['maj_vot_chr'].replace({999: 'UNKNOWN'})
source_to_save.to_csv("source_unique_df_preprocessed.csv")

def save_fasta(current_chr, fasta_chunk, records_for_fasta):
    
    # chr changed to reorder, reset to original
	if current_chr == 999:
		current_chr = 'UNKNOWN'
	elif current_chr == 998:
		current_chr = 'ERROR'
    
	fasta_name = f"fasta_output_unique/{current_chr}_{fasta_chunk}.fasta"
	with open(fasta_name, "w") as output_handle:
		SeqIO.write(records_for_fasta, output_handle, "fasta")

# loop over rows in table

# initial values
fasta_chunk = 1  # n fasta has been chunked per chr
current_chr = 1 # current chromosome being used
count_per_fasta = 0  # records in fasta
records_for_fasta = []  # records to be written to fasta

for i in tqdm(range(len(source_df)), total=len(source_df)):
    
    # itterows breaks because idx has been changed
	row = source_df.iloc[i]

	if i == 0:
		print(i)
		print(row['seq'])
    
	new_chr = row['maj_vot_chr']
    
	if new_chr == current_chr:

		one_sequence = Seq(data=row['seq'])
		one_record = SeqRecord(one_sequence, id=str(i), description="")
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
		one_record = SeqRecord(one_sequence, id=str(i), description="")
		records_for_fasta.append(one_record)

		# reset the records going into the fasta
		current_chr = new_chr
		count_per_fasta = 1

# save final record
save_fasta(current_chr, fasta_chunk, records_for_fasta)