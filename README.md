# RCSB_cmap
Create contact maps of RCSB structures
- Convient tool for quickly creating contact maps of proteins of interest from [rcsb.org](https://www.rcsb.org/)
- This can also generate contact maps of local .pdb files

### Install
```
pip install -r requirements.txt

bash test.sh
```

## Available flags:
- -pdb ####      |  RCSB pdb id for the protein of interest (example: 1fha)
                 |  if ####.pdb(cif) is in cwd the local file will be used (example: 1fha.pdb)
- -chain #       |  Chain id of interest (example: A)
- -pdb2 ####     |  used to create a dualfold comparison contact map with -pdb
- -chain2 #      |  specify the chain of -pdb_2 that will be used for the comparison
- -oligomer      |  If -pdb is an oligomer this flag will create a superposition of of all contacts
- -chains_like # |  retain only chains that are similar to the selected chain. Useful for hetero-oligomer (example: C)
- -leven #   |
**NOTE**: chains_like calculates the levenshtein distance between chains and retains chains that are within 30
**NOTE**: of the adjust -leven if this is to restrictive/permissive

![temporary text](/png/1a5m.png)
