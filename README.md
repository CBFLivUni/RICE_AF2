# Rice AF2

## Prerequisites

### Download the data
- Depending on how far though the pipeline you want to run the code from, download the relevant data from [Dropbox](https://www.dropbox.com/scl/fo/ow3quj34qnc9ssx1t3s48/ABPHfAW83K1I51vurVwVkqo?rlkey=yh2468xp4xlhqez6cnyyh77z4&st=qu9mn7py&dl=0).
- `data/colabfold_output` contains colabfold output files
- `data/fastas_for_colabfold` contains the input fasta files to colabfold
- `data/foldseek_output` contains the foldseek output files
- `data/raw` contains the protein sequences from each source

### Install the software
- Install [colabfold](https://github.com/YoshitakaMo/localcolabfold?tab=readme-ov-file#for-linux)
- Install [foldseek](https://github.com/steineggerlab/foldseek?tab=readme-ov-file#installation)
- Install [foldseek databases](https://github.com/steineggerlab/foldseek?tab=readme-ov-file#databases) to `data/foldseek_databases`. In this example only PDB was used.
- If using multiple local GPUs, such as via a compute engine on GCP then `pip install simple-gpu-scheduler`

### Other setup
- `other/GCP_doc.pdf` includes setup of compute engines on GCP.
- If some sequences have already been run on other machines and you want to avoid duplicating effort, rather than editing fastas, it's easier to copy the `.done.txt` from the `data/colabfold_output` folder. These files tell colabfold if predicted have already been generated for a sequence and colabfold will skip them.
  - `other/move_done_txt.py` can assist with this.
  - The scripts assume there's a conda environment on Barkla called `my_gpu` that is a copy of `gpu`, and a venv called `af2venv`.
 
## Running the pipeline

### Generate table of unique sequences from sources
- `preprocess/gen_seq_data_table.sh`
  - Merge sources and get unique sequences
  - Saves as `source_unique_df.csv`

### Generate fastas of these sequences
- `preprocess/gen_priority_1_fastas.sh`
  - For priority one sequences, get chromosome of sequences and put into fasta's of 20 sequences, where the header is the corresponding id in `source_unique_df.csv`
  - See `other/gen_fasta.sh` and `gen_fasta.py` for an example of how to generate fastas for all sequences, not just those that were priority for this project.
 
### Run colabfold
-
