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
        coords = []
        for chain in self.structure.get_chains():
            residues = sorted(
                Selection.unfold_entities(chain, "R"), key=lambda r: r.get_id()[1]
            )
            for res in residues:
                if "CA" in res:
                    coords.append(res["CA"].get_coord())
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
