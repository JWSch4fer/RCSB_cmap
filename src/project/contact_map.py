### src/project/contact_map.py

import numpy as np
import itertools
from typing import Tuple, Optional
from Bio.PDB import Selection
from scipy.spatial import distance_matrix
from .utils import levenshtein_distance

TRANSLATE_AA = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}


class ContactMap:
    """Compute residue–residue contact maps from a Biopython Structure."""

    def __init__(self, structure, cutoff: float = 8.0) -> None:
        """
        Args:
            structure: Biopython Structure object.
            cutoff: Å distance threshold.
        """
        self.structure = structure
        self.cutoff = cutoff
        self.coordinates = self._extract_ca_coordinates()

    def _extract_ca_coordinates(self) -> np.ndarray:
        # some RCSB enteries are poorly written, we'll deal with empty chains here
        detach_list = []

        for chain in self.structure[0].get_chains():
            residues = list(Selection.unfold_entities(chain, "R"))
            if not residues:
                detach_list.append(chain)
                continue

            # Iterate over each consecutive pair of residues
            for prev_res, curr_res in itertools.pairwise(residues):
                prev_idx = prev_res.get_id()[1]
                curr_idx = curr_res.get_id()[1]

                # Generate missing indices lazily
                for missing_idx in range(prev_idx + 1, curr_idx):
                    new_res = curr_res.copy()
                    # Detach parent and set placeholder coords
                    for atom in new_res:
                        atom.detach_parent()
                        atom.set_coord(np.array([9999, 9999, 9999], dtype="float32"))
                    new_res.id = (" ", missing_idx, chain.get_id())
                    chain.add(new_res)

        # Remove any chains that had zero residues
        for ch in detach_list:
            self.structure[0].detach_child(ch.id)

        return np.array(coords, dtype=float)

    def compute_contact_map(self) -> np.ndarray:
        dm = distance_matrix(self.coordinates, self.coordinates)
        cmap = (dm < self.cutoff).astype(int)
        np.fill_diagonal(cmap, 0)
        return cmap

    def collapse_homo(self, cmap: np.ndarray, n_chains: int) -> np.ndarray:
        L = cmap.shape[0] // n_chains
        intra = sum(
            cmap[i * L : (i + 1) * L, i * L : (i + 1) * L] for i in range(n_chains)
        )
        inter = sum(
            cmap[i * L : (i + 1) * L, j * L : (j + 1) * L]
            + cmap[j * L : (j + 1) * L, i * L : (i + 1) * L]
            for i, j in itertools.permutations(range(n_chains), 2)
        )
        intra = (intra > 0).astype(int)
        inter = (inter > 0).astype(int) * 2
        return intra + inter
