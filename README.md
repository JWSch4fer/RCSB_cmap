# Create a contact map from RCSB

This python-CLI enables users to create contact maps of proteins directly from the [RCSB](https://www.rcsb.org/) Protein Data Bank (PDB). By specifying a protein's RCSB ID, one can download the structure, select specific chains, and generate contact maps.

Key Features:
- Download by RCSB ID: Retrieve protein structures from the RCSB website using Python's requests module.
- Local structures (pdb/cif) can also be used to generate contact maps
- Chain Specification: Specify individual protein chains for focused analysis.
- Dual-Fold Contact Maps: allow the creation of dual-fold contact maps (one protein is displayed on the upper triangle and another on the lower triangle, making it easier to visualize and compare structural differences)

## Installation
```
pip install -r requirements.txt

bash test.sh
```

## Usage:
```
 -pdb ####      |  RCSB pdb id for the protein of interest (example: 1fha)
                |  if ####.pdb(cif) is in cwd the local file will be used (example: 1fha.pdb)
 -chain #       |  Chain id of interest (example: A)
 -pdb2 ####     |  used to create a dualfold comparison contact map with -pdb
 -chain2 #      |  specify the chain of -pdb_2 that will be used for the comparison
 -oligomer      |  If -pdb is an oligomer this flag will create a superposition of of all contacts
 -chains_like # |  retain only chains that are similar to the selected chain. Useful for hetero-oligomer (example: C)
 -leven #       |
```
**NOTE**: chains_like calculates the levenshtein distance between chains and retains chains that are within 30 deletions/insertions/mutations

### Examples
Full contact associated with RCSB ID 1A5M.
```
python3 RCSB_cmap.py -pdb 1a5m
```
![temporary text](/img/1a5m.png)


Use -chains_like flag to find all similar amino acid chains in a protein complex

| All chains with sequences like chain C  | Collapsed version                 |
| ------------------------------- | ----------------------------------------- |
|![](/img/1a5m_chains_like_C.png) | ![](/img/1a5m_chains_like_C_collapse.png) |
|```python3 RCSB_cmap.py -pdb 1a5m -chains_like C``` |```python3 RCSB_cmap.py -pdb 1a5m -chains_like C -oligomer ```|

Specify both -pdb and -pdb2 to create a dual-fold contact map to compare similar structures
```
python3 RCSB_cmap.py -pdb 1rep -oligomer -pdb2 2z9o
```
![](/img/comp_1rep_2z9o_collapse.png)


### 🤝 Contributing

Contributions are welcome! Please fork the repository and create a pull request with your changes. For major changes, please open an issue to discuss what you would like to change.
