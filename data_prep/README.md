## Dataset 
Data were prepared based on [tutorial.ipynb](data_prep/tutorial.ipynb)

- data are already splitted to train and test (see cryptobench/cryptobench-dataset/folds folder)
- 20 train examples.json is here, in data_prep folder 
- to get the data:
1) activate any venv with BioPython 
2) run get_data.sh

./get_data.sh does: 
1) downloads and unzips cryptobench 
2) extracts sequences from APO structures, for which pockets are defined only on 1 chain. No multichain. 


### Details 
- 20_train_sequences.json  and 20_train_structures.json are 20 APO proteins I chose for initial experimenting. There are no APO structures with multichain pockets.
- these 2 files correspond to sequences.json and dataset.json in cryptobench folder. However, there are no multichain structures in sequences.json (while they are in dataset.json). 
- sequences.json contains sequences, author residue indices and label residue indices for all 946 APO structures without multichain pockets. 

- data_prep.ipynb -> my testing of some utils for cryptobench to get chains, sequences and so on. 
- tutorial.ipynb -> just copied from cryptobench repo. Useful, but you can skip it, if you don't need to retrive some data from dataset that don't have function in utils/

## P2Rank
First, download p2rank via
```
get_p2rank.sh
```
- this downloads it from the repo and extracts it into the folder `../p2rank`

Then, make sure you have a Java (OpenJDK) 17 (or higher) installed
- on MetaCentrum, run `module add openjdk/17`

You can run the p2rank predictions with
```
run_p2rank.sh
```
