# RCSB_cmap
Create contact maps of RCSB structures
- Convient tool for quickly creating contact maps of proteins of interest from [rcsb.org](https://www.rcsb.org/)
- This can also generate contact maps of local .pdb files

### Build a MiniConda environment
bash conda_setup.sh
bash test.sh

## Usage
python RCSB_cmap.py -pdb ####(.pdb) -chain # -pdb2 ####(.pdb) -chain2 # -oligomer
- pdb "flag to specify protein pdbid (to run a local file pass the extension otherwise this will download the structure)"
- chain "specify a chain of interest"
- pdb2 "flag to specify protein pdbid (to run a local file pass the extension otherwise this will download the structure)"
- chain2 "specify a chain of interest"
- oligomer "If a protein is a homo-oligomer this flag will collapse the redundant chains down to one with inter/intra chains denoted"


