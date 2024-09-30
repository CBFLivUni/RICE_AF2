import pandas as pd
import helper_funcs as funcs
import os

# import
rap_db_off = funcs.read_fasta_to_dict(os.path.join("..", "..", "data", "raw", "IRGSP-1.0_protein_2024-01-11.fasta"))
rap_db_pred = funcs.read_fasta_to_dict(os.path.join("..", "..", "data", "raw", "IRGSP-1.0_predicted-protein_2024-01-11.fasta"))
ens_off = funcs.read_fasta_to_dict(os.path.join("..", "..", "data", "raw", "Oryza_sativa.IRGSP-1.0.pep.all.fa"))
ens_pred = funcs.read_fasta_to_dict(os.path.join("..", "..", "data", "raw", "Oryza_sativa.IRGSP-1.0.pep.abinitio.fa"))
msu_off = funcs.read_fasta_to_dict(os.path.join("..", "..", "data", "raw", "all.pep"))
osnip = funcs.read_fasta_to_dict(os.path.join("..", "..", "data", "raw", "oryza_sativa_Nip.protein.fasta"))
gram_pred = funcs.read_fasta_to_dict(os.path.join("..", "..", "data", "raw", "gram_Oryza_sativa.IRGSP-1.0.pep.abinitio.fa"))
uniprot = funcs.read_fasta_to_dict(os.path.join("..", "..", "data", "raw", "uniprotkb_proteome_UP000059680_2024_07_04.fasta"))

# msu transposable elements
msu_df = pd.read_csv(os.path.join("..", "..", "data", "raw", "all.locus_brief_info.7.0"), sep='\t')
te_msu_list = list(set(msu_df[msu_df['is_TE'] == "Y"]['model'].to_list()))
te_msu_off = {key: msu_off[key] for key in te_msu_list}

source_pro = {"RAP-DB - official": rap_db_off,
        	"RAP-DB - predicted": rap_db_pred,
        	"Ensembl - official": ens_off,
        	"Ensembl - predicted": ens_pred,
        	"MSU - official": msu_off,
            "MSU - TE": te_msu_off,
        	"OsNip": osnip,
      		"UniProt" : uniprot}

funcs.count_plot(source_pro)

# remove stop codon if present for consistency across sources
rm_stop_source_pro = funcs.remove_stop(source_pro)

# check got all sequences from all sources
assert sum(len(v) for v in source_pro.values()) == sum(len(v) for v in rm_stop_source_pro.values())

funcs.gen_upset_plot(rm_stop_source_pro)

# get all proteins
all_seqs = []
for records in rm_stop_source_pro.values():
	for seq in records.values():
		all_seqs.append(seq)

# check got all sequences from all sources
assert sum(len(v) for v in source_pro.values()) == len(all_seqs)

# only use unique sequences
unique_seqs = list(set(all_seqs))

source_df = funcs.gen_seq_source_chr_table(rm_stop_source_pro, unique_seqs, filter=None)

source_df.to_csv(os.path.join("..", "..", "data", "source_unique_df.csv"))
