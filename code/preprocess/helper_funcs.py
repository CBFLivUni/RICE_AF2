import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import seaborn as sns
import matplotlib.pyplot as plt
import upsetplot as up
import json
from tqdm import tqdm
import os

def read_fasta_to_dict(file_path):
	seq_dict = {}
	records = list(SeqIO.parse(file_path, "fasta"))
	for record in records:
		seq_dict[str(record.id)] = record.seq

	return seq_dict


def count_plot(source_pro):

	len_pro = {}
	[len_pro.update({p: len(source_pro[p])}) for p in source_pro.keys()]

	# count
	count_df = pd.DataFrame(data={'Source': list(len_pro.keys()),
					'Count': list(len_pro.values())})

	count_df['Count (k)'] = count_df['Count'] / 1000

	fig, ax = plt.subplots(figsize=(9, 7))

	sns.barplot(data=count_df,
				x="Count (k)",
				y="Source",
				ax=ax)

	plt.tight_layout()

	# Hide the right and top spines
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	# Only show ticks on the left and bottom spines
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')

	fig.savefig(fname="plots/barplot_count.png")


def remove_stop(source_pro):
	# remove stop codon if present for consistency across sources
	for source in source_pro.keys():
		
		for fasta_header, seq in source_pro[source].items():
			if seq.endswith('*'):
				source_pro[source][fasta_header] = seq[:-1]

	return source_pro


def gen_upset_plot(rm_stop_source_pro):

	# Call set on the list of sequences, and convert to string
	# both of those are required to work with upsetplot
	# set explains discrepancy between the number of sequences in seq record and string records
	set_dict_as_string = {}

	for source, records in rm_stop_source_pro.items():

		list_seqs= list(records.values())
		unique_seqs = list(set(list_seqs))
		unique_seqs_str = [str(seq) for seq in unique_seqs]

		set_dict_as_string[source] = unique_seqs_str

	# upset plot
	fig = plt.figure()
	pro_set = up.from_contents(set_dict_as_string)
	upset_plt = up.plot(pro_set,
						fig=fig,
						element_size=50,
						sort_by="cardinality",
						facecolor="#434549",
						min_subset_size=20,
						show_counts='%,d')
	plt.savefig("plots/upset_plot_all.png") 

	# RAP-DB and Emsemble only
	rapdb_ensemble = {key: set_dict_as_string[key] for key in
					["RAP-DB - official", "RAP-DB - predicted", "Ensembl - official", "Ensembl - predicted"]}

	fig = plt.figure()
	pro_set = up.from_contents(rapdb_ensemble)
	upset_plt = up.plot(pro_set,
						fig=fig,
						element_size=70,
						sort_by="cardinality",
						facecolor="#434549",
						min_subset_size=20,
						show_counts='%,d')
	plt.savefig("plots/upset_plot_rapdb_ensemble.png")


class ChromosomeGetter:

	"""
	Use class to avoid passing large uniprot json around
	"""
	def __init__(self):
		# uniprot taxonomy to get chromosome
		with open(os.path.join("..", "..", "data", "uniprotkb_taxonomy_id_4530_2024_07_17.json")) as f:
			self.uniprot_json = json.load(f)

	def getchromosome(self, record_id):

		try:

			if "OsNip" in record_id:
				# then osnip
				chromosome = str(int(record_id.split("OsNip_")[1][:2]))

			elif "Os" in record_id:
				# then msu, ens, rap
				chromosome = str(int(record_id.split("Os")[1][:2]))

			elif "|" in record_id:
				# then uniprot
				# find record id in uniprot json

				chromosome = "UNKNOWN" # set default and only change if can get record

				accession = record_id.split("|")[1]
				
				uniprot_rec_idx = None  # set default
				for idx, uniprot_records in enumerate(self.uniprot_json['results']):
					if uniprot_records['primaryAccession'] == accession:
						uniprot_rec_idx = idx

				if uniprot_rec_idx is not None:
					unipro_db_records = self.uniprot_json['results'][uniprot_rec_idx]['uniProtKBCrossReferences']
					# get first proteomes db record and use that as chromosome
					for unipro_db_record in unipro_db_records:
						if unipro_db_record['database'] == "Proteomes":
							# check record is for chromosome
							db_record_value = unipro_db_record['properties'][0]['value']
							if "Chromosome" in db_record_value:
								chromosome = db_record_value.split("Chromosome ")[1]

			else:
				# predicted and therefore unknown
				chromosome = "UNKNOWN"

			return chromosome
		
		except:
			print("ERROR")
			print(record_id)
			return "ERROR"
	
def most_common(lst):
	# remove "UNKNOWN" from list
	flatten_list = []
	for l in lst:
		if ',' in l:
			flatten_list.extend(l.split(","))
		else:
			flatten_list.append(l)

	sanitise_list = [x for x in flatten_list if x != 'UNKNOWN' and x != 'ERROR']
	# if list is empty return "UNKNOWN"
	if len(sanitise_list) == 0:
		return "UNKNOWN"

	# then return majority vote
	return max(set(sanitise_list), key=sanitise_list.count)


def gen_seq_source_chr_table(rm_stop_source_pro, all_seqs, filter=None):
	# all_seqs is all sequences from all sources
	# rm_stop_source_pro is the dictionary of sources and their sequences
	# filter is the number of sequences to filter to

	# encode which sequence is from which source in a df
	# each data source is a column, each row is a sequence
	# each value is a list of headers that are in that source

	if filter is not None:
		all_seqs = all_seqs[:filter]

	getChr = ChromosomeGetter()

	source_df = pd.DataFrame(columns=list(rm_stop_source_pro.keys()) +
							['seq'] +
							[source + " - chr" for source in list(rm_stop_source_pro.keys())] +
							['maj_vot_chr'])
	for seq in tqdm(all_seqs):
		row = {}
		for source in rm_stop_source_pro.keys():
			# get idx of seq in source
			sequences_from_source = list(rm_stop_source_pro[source].values())
			seq_idx = [i for i in range(len(sequences_from_source)) if sequences_from_source[i] == seq]

			# if seq in source
			if len(seq_idx) != 0:
				# append headers from those sequences
				headers_from_source = list(rm_stop_source_pro[source].keys())
				headers = [headers_from_source[idx] for idx in seq_idx]
				row[source] = ",".join(headers)

				# get chromosome from header
				chromosomes = [getChr.getchromosome(header) for header in headers]
				row[source + " - chr"] = ",".join(chromosomes)

			row['seq'] = str(seq)
			maj_vot_chr = most_common([v for k, v in row.items() if "chr" in k])  # take majority vote from all chromosome rows
			row['maj_vot_chr'] = maj_vot_chr

		source_df = pd.concat([source_df, pd.DataFrame(row, index=[0])])
	
	return source_df


def make_priority_list(source_df, priority_list):
	"""
	Takes list of transcripts and checks if they are in the source DataFrame, returning a list of 'T' or 'F' for each row.
	Args:
		source_df (pd.DataFrame): The source DataFrame containing the data to be checked.
		priority_list (list): A list of transcript names that are considered priority.
	Returns:
		list: A list of 'T' or 'F' indicating whether each row in the source DataFrame contains any of the priority IDs.
	"""
    
	priority_T_F = []
    
	print('Making priority list...')
	for i, row in tqdm(source_df.iterrows(), total=source_df.shape[0]):
        
        # check every column for a match.
		id_cols = ['RAP-DB - official',
			'Ensembl - official', 'MSU - official',
			'MSU - TE', 'OsNip', 'UniProt']

		# turn row into string
		source_string = '|'.join(row[id_cols].values.astype(str))

		# set default to False
		is_priority = 'F'
  
		for id in priority_list:
			if id in source_string:
				is_priority = 'T'
				break

		priority_T_F.append(is_priority)
  
	return priority_T_F


def save_fasta(current_chr, fasta_chunk, records_for_fasta, JOB_NAME):
	"""
	Saves the records to a fasta file, with file name as "chr_chunk.fasta".
	"""
    
    # chr changed to reorder, reset to original
	if current_chr == 999:
		current_chr = 'UNKNOWN'
	elif current_chr == 998:
		current_chr = 'ERROR'
    
	fasta_name = os.path.join('fastas_for_colabfold', JOB_NAME, f'{current_chr}_{fasta_chunk}.fasta')
	with open(fasta_name, "w") as output_handle:
		SeqIO.write(records_for_fasta, output_handle, "fasta")


def to_priority_df(source_df, priority):
	"""
	Returns a DataFrame containing only the rows where the priority column is 'T'.
	"""
    
	priority_df = source_df[source_df[priority] == 'T']
	priority_df['idx'] = priority_df.index
	return priority_df


def make_fastas(priority_df, JOB_NAME):
	"""
	Splits the priority DataFrame into chunks of 20 and saves them to a fasta file.
	"""
    
	# initial values
	fasta_chunk = 1  # n fasta has been chunked per chr
	current_chr = 1 # current chromosome being used
	count_per_fasta = 0  # records in fasta
	records_for_fasta = []  # records to be written to fasta

	print('Making fastas...')
	for i in tqdm(range(len(priority_df)), total=len(priority_df)):
		
		# itterows breaks because idx has been changed
		row = priority_df.iloc[i]

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
				
				save_fasta(current_chr, fasta_chunk, records_for_fasta, JOB_NAME)

				fasta_chunk += 1

				# reset the records going into the fasta
				current_chr = new_chr
				count_per_fasta = 0
				records_for_fasta = []

		else:
			# must be new chromosome
			# check there are records to save, then save final fasta for that chr
	
			if len(records_for_fasta) != 0:
				save_fasta(current_chr, fasta_chunk, records_for_fasta, JOB_NAME)

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
	save_fasta(current_chr, fasta_chunk, records_for_fasta, JOB_NAME)


def get_sequences(source_df, priority_list):
	"""
	Given a list of 'T' and 'F' values, returns a list of sequences from the source DataFrame where the value is 'T'.
	"""
	
	priority_seqs = []
 
	# check every column for a match.
	id_cols = ['RAP-DB - official',
		'Ensembl - official', 'MSU - official',
		'MSU - TE', 'OsNip', 'UniProt']

	print('Getting sequences...')
 
	all_seqs = source_df['seq'].tolist()
 
	for i in tqdm(range(len(priority_list))):
		if priority_list[i] == 'T':
			priority_seqs.append(all_seqs[i])
	
	return priority_seqs


def check_for_duplicates(old_source_df, priority_1_all, priority_2_all, priority_3_all, priority_4_all):
	"""
	Checks for duplicates between the priority lists and returns a list of sequences that are unique to each priority list.
	If a sequence is in an earlier priority list, it will not be included in the later priority lists.
	"""
    
	priority_1_seqs = get_sequences(old_source_df, priority_1_all)
	priority_2_seqs = get_sequences(old_source_df, priority_2_all)
	priority_3_seqs = get_sequences(old_source_df, priority_3_all)
	priority_4_seqs = get_sequences(old_source_df, priority_4_all)
    
    # check sequences against each other
	priority_2_not_inc_1 = []
	priority_3_not_inc_1_2 = []
	priority_4_not_inc_1_2_3 = []
 
	print('Checking for duplicates...')
 
	# check not in 2 for 1
	priority_1_set_seqs = list(set(priority_1_seqs))
 
	for seq in priority_2_seqs:
		if seq not in priority_1_set_seqs:
			priority_2_not_inc_1.append(seq)
   
	# check not in 3 for 1, 2
	priority_2_set_seqs = list(set(priority_2_seqs))
 
	for seq in priority_3_seqs:
		if seq not in priority_1_set_seqs and seq not in priority_2_set_seqs:
			priority_3_not_inc_1_2.append(seq)
   
	# check not in 4 for 1, 2, 3
	priority_3_set_seqs = list(set(priority_3_seqs))
 
	for seq in priority_4_seqs:
		if seq not in priority_1_set_seqs and seq not in priority_2_set_seqs and seq not in priority_3_set_seqs:
			priority_4_not_inc_1_2_3.append(seq)
   
	priority_4_set_seqs = list(set(priority_4_seqs))
 
	return priority_1_set_seqs, priority_2_set_seqs, priority_3_set_seqs, priority_4_set_seqs


def add_cols_to_df(source_df, priority_1_set_seqs, priority_2_set_seqs, priority_3_set_seqs, priority_4_seqs):
	"""
	Adds columns to the source DataFrame indicating whether a sequence is in a priority list.
	"""

	priority_1_col = []
	priority_2_col = []
	priority_3_col = []
	priority_4_col = []

	seqs = source_df['seq'].tolist()
 
	print('Adding columns to df...')
 
	for seq in tqdm(seqs):
		if seq in priority_1_set_seqs:
			priority_1_col.append('T')
		else:
			priority_1_col.append('F')
 
		if seq in priority_2_set_seqs:
			priority_2_col.append('T')
		else:
			priority_2_col.append('F')
 
		if seq in priority_3_set_seqs:
			priority_3_col.append('T')
		else:
			priority_3_col.append('F')
 
		if seq in priority_4_seqs:
			priority_4_col.append('T')
		else:
			priority_4_col.append('F')
   
	source_df['priority_1'] = priority_1_col
	source_df['priority_2'] = priority_2_col
	source_df['priority_3'] = priority_3_col
	source_df['priority_4'] = priority_4_col
 
	return source_df
