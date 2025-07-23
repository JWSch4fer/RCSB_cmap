### src/project/contact_map.py
import numpy as np
from sklearn.neighbors import radius_neighbors_graph
from scipy.sparse import csr_matrix
import itertools
from typing import Tuple, Optional
from Bio.PDB import Selection
from scipy.spatial import distance_matrix
from .utils import levenshtein_distance, plot_contact_map
import sys

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
        self._process_protein()
        self.all_coords, self.all_coords_ids = self._extract_coordinates()

    def _process_protein(self) -> None:
        atom_number = 1
        aa_return = []
        aa = list(TRANSLATE_AA.keys())
        for chain in self.structure.get_chains():
            detach_residue_list = set()
            residues = Selection.unfold_entities(chain, "R")
            for residue in residues:
                detach_id_list = []
                # remove residues that don't have a CA
                if "CA" not in residue.child_dict.keys():
                    detach_residue_list.add(residue._id)

                # remove anything that is not protein
                if residue.full_id[-1][0] != " ":
                    detach_residue_list.add(residue._id)
                if residue.resname not in aa:
                    detach_residue_list.add(residue._id)

                if residue.resname in aa:
                    aa_return.append(residue.resname)
                    for atom in residue:
                        # remove hydrogens
                        if atom.element.strip() == "H":
                            detach_id_list.append(atom.id.strip())
                        if atom.element != "H":
                            atom.serial_number = atom_number
                            atom_number += 1
                # remove all the atoms we're not interested in
                for atom_id in detach_id_list:
                    residue.detach_child(atom_id)

            # remove everything other than the protein
            for residue in detach_residue_list:
                chain.detach_child(residue)

    def _extract_coordinates(self) -> Tuple[np.ndarray, np.ndarray]:
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
                if prev_idx + 1 != curr_idx:
                    self.fill_gap(
                        res_start=prev_idx,
                        res_stop=curr_idx,
                        chain=chain,
                        curr_res=curr_res,
                    )

        # Remove any chains that had zero residues
        for ch in detach_list:
            self.structure[0].detach_child(ch.id)

        all_coords = [
            (atom.coord, chain._id + str(residue.get_id()[1]))
            for chain in self.structure[0].get_chains()
            for residue in Selection.unfold_entities(chain, "R")
            for atom in residue.get_list()
        ]

        # return coordinates and identifiers
        _ = np.array(all_coords, dtype=object)
        all_coords = np.stack(_[:, 0], axis=0)
        all_coords_ids = np.stack(_[:, 1], axis=0)
        return all_coords, all_coords_ids

    def fill_gap(self, res_start: int, res_stop: int, chain, curr_res) -> None:
        # Generate missing indices lazily
        for missing_idx in range(res_start + 1, res_stop):
            new_res = curr_res.copy()
            # Detach parent and set placeholder coords
            for atom in new_res:
                atom.detach_parent()
                atom.set_coord(np.array([9999, 9999, 9999], dtype="float32"))
            new_res.id = (" ", missing_idx, chain.get_id())
            chain.add(new_res)

    def compute_residue_contact_map(self) -> np.ndarray:
        """
        Compute a residue‐level contact map by aggregating heavy‐atom distances.

        Args:
            structure: Biopython Structure object containing chains and atoms.
            cutoff: distance threshold (in Å) defining a contact.

        Returns:
            R×R binary contact map, where R is the total number of residues.
        """
        coords = self.all_coords
        resid_idx = self.all_coords_ids

        # Compute H×H pairwise distances and threshold all atoms
        A = radius_neighbors_graph(
            X=coords,
            radius=self.cutoff,
            mode="connectivity",
            metric="euclidean",
            include_self=False,
            n_jobs=-1,
        )

        valid = ~np.all(coords == 9999.0, axis=1)  # deal with placeholder atoms
        mask2d = valid[:, None] & valid[None, :]
        A = A.multiply(mask2d)

        # Build one‐hot encoding (R × H)
        unique_res = np.unique(resid_idx)
        R = unique_res.size
        # Map residue sequence numbers to 0…R‑1 indices
        res_to_idx = {res: i for i, res in enumerate(unique_res)}
        mapped = np.vectorize(res_to_idx.get)(resid_idx)  # length H
        one_hot = np.zeros((R, coords.shape[0]), dtype=int)
        one_hot[mapped, np.arange(coords.shape[0])] = 1  # O(H) assignment
        one_hot_sparse = csr_matrix(one_hot)

        # Aggregate via matrix multiplication
        # this stip must be sparse matrix operation
        contact_counts = one_hot_sparse @ A @ one_hot_sparse.T

        contact_map = (contact_counts > 0).astype(int)  # Binarize >0

        # Remove diagonal and near‐diagonal hits
        idxs = np.arange(R)
        diag_mask = np.abs(idxs[:, None] - idxs[None, :]) <= 3
        contact_map[diag_mask] = 0

        plot_contact_map(contact_map.toarray())
        sys.exit()
        return contact_map

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
