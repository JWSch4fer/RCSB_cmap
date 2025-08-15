# Create a contact map from RCSB

This python-CLI enables users to create contact maps of proteins directly from the [RCSB](https://www.rcsb.org/) Protein Data Bank (PDB). By specifying a protein's RCSB ID, one can download the structure, select specific chains, and generate contact maps.

Key Features:
- Download by RCSB ID: Retrieve protein structures from the RCSB website using Python's requests module.
- Local structures (pdb/cif) can also be used to generate contact maps
- Chain Specification: Specify individual protein chains for focused analysis.
- *Chains like* Specification: Specify a chain and the CLI will search the protein complex for chains like the one specified(similarity is calculated using Levenshtein distance).
- Dual-Fold Contact Maps: allow the creation of dual-fold contact maps (one protein is displayed on the upper triangle and another on the lower triangle, making it easier to visualize and compare structural differences)

## Installation
```
pip install -r requirements.txt
pip install .

```

## Usage:
```
usage: rcsb_cmap [-h] -p PDB [-c CHAIN] [-p2 PDB2] [-c2 CHAIN2] [-o] [--cutoff CUTOFF] [--chains_like CHAINS_LIKE] [--levenshtein LEVENSHTEIN]

Generate residue–residue contact maps from PDB/mmCIF structures.

options:
  -h, --help            show this help message and exit
  -p PDB, --pdb PDB     PDB ID or local file path (.pdb or .cif)
  -c CHAIN, --chain CHAIN
                        Chain ID (e.g. A)
  -p2 PDB2, --pdb2 PDB2
                        PDB ID or local file path (.pdb or .cif)
  -c2 CHAIN2, --chain2 CHAIN2
                        Chain ID (e.g. A)
  -o, --oligomer        Collapse homo‑oligomer contacts
  --cutoff CUTOFF       Distance cutoff in Å for contacts
  --chains_like CHAINS_LIKE
                        retain only chains that are similar to the selected chain. Useful for hetero-oligomer (example: C)
  --levenshtein LEVENSHTEIN
                        chains_like calculates the levenshtein distance between chains and retains chains that are within 30 of the adjust if this is to restrictive/permissive```
**NOTE**: chains_like calculates the levenshtein distance between chains and retains chains that are within 30 deletions/insertions/mutations

### Examples
Full contact associated with RCSB ID 1A5M.
```
rcsb_cmap -p 1a5m
```
![temporary text](/img/1a5m.png)


Use -chains_like flag to find all similar amino acid chains in a protein complex

| All chains with sequences like chain C  | Collapsed version                 |
| ------------------------------- | ----------------------------------------- |
|![](/img/1a5m_chains_like_C.png) | ![](/img/1a5m_chains_like_C_collapse.png) |
|```rcsb_cmap -p 1a5m --chains_like C``` |```rcsb_cmap -p 1a5m -chains_like C -o ```|

Specify both --pdb and --pdb2 to create a dual-fold contact map to compare similar structures
```
rcsb_cmap -p 1rep -o -p2 2z9o
```
![](/img/comp_1rep_2z9o_collapse.png)


### 🤝 Contributing

Contributions are welcome! Please fork the repository and create a pull request with your changes. For major changes, please open an issue to discuss what you would like to change.
