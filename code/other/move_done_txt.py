import os

# get done.txt from colabfold_output
list_done = os.listdir('colabfold_output')
list_done_txt = [f for f in list_done if 'done.txt' in f]

# copy to new folder
os.system('mkdir done_txt')

for f in list_done_txt:
	os.system(f'cp colabfold_output/{f} done_txt/')