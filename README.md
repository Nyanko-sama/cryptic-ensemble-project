# Introduction

# Repo structure
```
|---data_prep   
|---data
|   |---cryptobench             # cryptobench dataset, contains dataset.json, sequences.json and cif-files/
|   |---bioemu
|---bioemu                      # scripts for running bioemu, installation guide later 
|---utils
|   |---working_with_cryptobenchpy          # some functions for cryptobench. Stolen from tutorial :) 
```

# Dataset
 CryptoBench benchmark, containing PDB IDs of apo-holo structure pairs along with the corresponding binding site residue annotations (ground truth labels), together with the corresponding CIF files. The full dataset is also available from the [CryptoBench OSF](https://osf.io/pz4a9/wiki?wiki=e7tdg) project site.

To obtain dataset: 
- activate any venv with BioPython
- run get_data.sh

# BioEmu
I had prepared env on metacentrum. So no intallation guide. I will add it later. 
