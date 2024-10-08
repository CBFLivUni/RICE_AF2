import pandas as pd
import helper_funcs as funcs
import os

# loads the priority lists
priority_1_df = pd.read_csv(os.path.join('Nipponbare_new_filter_Priority_1.tab'), sep='\t')
priority_2_df = pd.read_csv(os.path.join('Nipponbare_priority_2.tab'), sep='\t')
priority_3_df = pd.read_csv(os.path.join('Nipponbare_priority_3.tab'), sep='\t')
priority_4_df = pd.read_csv(os.path.join('Nipponbare_priority_4.tab'), sep='\t')

# get transcript name
priority_1_list = priority_1_df['child_transcripts'].tolist()
priority_2_list = priority_2_df['transcript'].tolist()
priority_3_list = priority_3_df['transcript'].tolist()
priority_4_list = priority_4_df['transcript'].tolist()

# master source df
source_df = pd.read_csv(os.path.join('source_unique_df.csv'))

# drop predicted
# note this will lead into slight differences with previous versions, as previously they were included
source_df = source_df.drop(columns=['RAP-DB - predicted', 'Ensembl - predicted'])

# generate lists for each priority inc duplicates
priority_1_all = funcs.make_priority_list(source_df, priority_1_list)
priority_2_all = funcs.make_priority_list(source_df, priority_2_list)
priority_3_all = funcs.make_priority_list(source_df, priority_3_list)
priority_4_all = funcs.make_priority_list(source_df, priority_4_list)

# check for duplicates across priorities then add as columns
priority_1_set_seqs, priority_2_set_seqs, priority_3_set_seqs, priority_4_set_seqs = funcs.check_for_duplicates(
    source_df, priority_1_all, priority_2_all, priority_3_all, priority_4_all
    )

# add to source df
source_df = funcs.add_cols_to_df(source_df, priority_1_set_seqs, priority_2_set_seqs, priority_3_set_seqs, priority_4_set_seqs)	

source_df.to_csv(os.path.join('source_unique_df_w_prior.csv'), index=False)

# generate fastas for each priority
priority_1_df = funcs.to_priority_df(source_df, 'priority_1')
priority_2_df = funcs.to_priority_df(source_df, 'priority_2')
priority_3_df = funcs.to_priority_df(source_df, 'priority_3')
priority_4_df = funcs.to_priority_df(source_df, 'priority_4')

funcs.make_fastas(priority_1_df, "priority_1")
funcs.make_fastas(priority_2_df, "priority_2")
funcs.make_fastas(priority_3_df, "priority_3")
funcs.make_fastas(priority_4_df, "priority_4")
