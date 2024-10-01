# Rice AF2

## Prerequisites

### Download the data
- Depending on how far though the pipeline you want to run the code from, download the relevant data from [Dropbox](https://www.dropbox.com/scl/fo/ow3quj34qnc9ssx1t3s48/ABPHfAW83K1I51vurVwVkqo?rlkey=yh2468xp4xlhqez6cnyyh77z4&st=qu9mn7py&dl=0).
- `data/colabfold_output` contains colabfold output files
- `data/fastas_for_colabfold` contains the input fasta files to colabfold
- `data/foldseek_output` contains the foldseek output files
- `data/raw` contains the protein sequences from each source

#### Note on rclone
- Data is backed up on Dropbox [here](https://www.dropbox.com/scl/fo/ow3quj34qnc9ssx1t3s48/ABPHfAW83K1I51vurVwVkqo?rlkey=yh2468xp4xlhqez6cnyyh77z4&st=qu9mn7py&dl=0).
- [rclone](https://rclone.org/install/) can make backing up on Dropbox easier.
- `sudo -v ; curl https://rclone.org/install.sh | sudo bash` to install
- `rsync config` to configure
- `rclone copy colabfold_output "remote:/Shared Data/AF2/colabfold_output" -P --ignore-existing` as an example, to copy the colabfold_output folder to Dropbox, observing progressing and ignoring duplicates.


### Install the software
- Install [colabfold](https://github.com/YoshitakaMo/localcolabfold?tab=readme-ov-file#for-linux)
- Install [foldseek](https://github.com/steineggerlab/foldseek?tab=readme-ov-file#installation)
- Install [foldseek databases](https://github.com/steineggerlab/foldseek?tab=readme-ov-file#databases) to `data/foldseek_databases`. In this example only PDB was used.
- If using multiple local GPUs, such as via a compute engine on GCP then `pip install simple-gpu-scheduler`

### Other setup
- `other/GCP_doc.pdf` includes setup of compute engines on GCP.
- If some sequences have already been run on other machines and you want to avoid duplicating effort, rather than editing fastas, it's easier to copy the `.done.txt` from the `data/colabfold_output` folder. These files tell colabfold if predicted have already been generated for a sequence and colabfold will skip them.
  - `coda/other/move_done_txt.py` can assist with this.
  - The scripts assume there's a conda environment on Barkla called `my_gpu` that is a copy of `gpu`, and a venv called `af2venv`.
 
## Running the pipeline

### Generate table of unique sequences from sources
- `code/preprocess/gen_seq_data_table.sh`
  - Merge sources and get unique sequences
  - Saves as `source_unique_df.csv`

### Generate fastas of these sequences
- `code/preprocess/gen_priority_1_fastas.sh`
  - For priority one sequences, get chromosome of sequences and put into fasta's of 20 sequences, where the header is the corresponding id in `source_unique_df.csv`
  - See `code/other/gen_fasta.sh` and `gen_fasta.py` for an example of how to generate fastas for all sequences, not just those that were priority for this project.
 
### Run colabfold
- `code/colabfold/run_colabfold_until_done.sh`
  - Will run colabfold on Barkla.
  - Sets up an array of 8 jobs (the max that can be requested at one time on Barkla) on the first GPU partittion that becomes available.
  - Edits the vars `fasta_path` and/ or `fasta_names` in `colabfold/run_colabfold_barkla.sh` if you want to use different fastas.
  - Save to `data/colabfold_output/`
- `code/colabfolabfold/run_colabfold_gcp.sh`
  - Will run colabfold on a compute engine on GCP.
  - Edits the vars `fasta_path` and/ or `fasta_names` in `colabfold/run_colabfold_barkla.sh` if you want to use different fastas.
  - Save to `data/colabfold_output/`
 
### Run foldseek
- `code/foldseek/foldseek_query.sh`
  - Run foldseek against the highest ranked relaxed models for each sequence in `data/colabfold_output/`
  - Outputs to `data/foldseek_output/`
  - This script ran on the PGB HPC.
 
### Process results
- `data/process_output/get_results.sh`
  - calculate some results from colabfold and foldseek outputs and save as csv.
 
## Other scripts

- Don't expect these scripts to work as part of the the pipeline, but here are some other scripts I used to solve problems I ran into.
- `code/other/move_done_txt.py`
  - copy `*done.txt` files from `colabfold_output` and move them to a new directory called `done_txt`.
  - `done.txt` files indicate to colabfold that a sequence has already been predicted and that these sequences can be skipped.
  - These can then be moved to `colabfold_output` on a different system to ensure that those fastas as skipped and results aren't duplicated when processing.
  - This allows all fastas to be moved, rather than selecting only those that need moving.
- `code/other/convert_old_sequences.sh`
  - Used to change the name of files that had already been generated, but incorrectly named
  - This opens `.csv`'s with the incorrect names (`old_source_df`) and correct names (`new_source_df`) and maps the old name to the new one.
  - Then renames the files in the incorrectly named fold (`colabfold_output`), and copies them to a new folder ('new_colabfold_output').
- `code/other/figs_for_presentation.py`
  - Opens the csv generated by `code/process_output/get_results.py`
  - In a folder called `plots`, violin plots are generated of the distribution of each feature by which dataset the sequence was obtained from.
- `code/other/count_sequences.py`
  - Get names of sequences that need running, check if results have been generated for them in `colabfold_output`.
  - Print count of sequences left to process.
- `code/other/gen_fasta.sh`
  - Similar to `code/preprocess/gen_priority_1_fasta.sh`, though doesn't check against being a priority 1 sequence.
- `code/other/rename_results.py`
  - Some fasta records were incorrectly named based on count, rather than a sequence id. This script renames those files.
