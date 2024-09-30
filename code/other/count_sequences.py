import pandas as pd
import os

source_df = pd.read_csv('main/source_unique_df_preprocessed_w_prior1.csv')

priority_1_df = source_df[source_df['priority_1'] == 'T']

# rename Unnamed column
priority_1_df.rename(columns={'Unnamed: 0': 'Name'}, inplace=True)


# list currently processed
list_current = os.listdir('colabfold_output')

list_env = [f for f in list_current if '_env' in f]

len(list_env)

# check if the file is already processed
count = 0
for i in range(len(priority_1_df)):
	if str(priority_1_df.iloc[i]['Name']) + '_env' in list_env:
		count += 1

print('Left to process:', len(priority_1_df) - count)