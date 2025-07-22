### src/project/contact_map.py
import numpy as np
from sklearn.neighbors import radius_neighbors_graph
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
        self.all_coords, self.all_coords_ids = self._extract_ca_coordinates()

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

            # remove everythin other than the protein
            for residue in detach_residue_list:
                chain.detach_child(residue)

    def _extract_ca_coordinates(self) -> Tuple[np.ndarray, np.ndarray]:
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

    # def compute_contact_map(self) -> np.ndarray:
    # dm = distance_matrix(self.coordinates, self.coordinates)
    # cmap = (dm < self.cutoff).astype(int)
    # np.fill_diagonal(cmap, 0)
    # return cmap

    def compute_residue_contact_map(self) -> np.ndarray:
        """
        Compute a residue‐level contact map by aggregating heavy‐atom distances.

        Args:
            structure: Biopython Structure object containing chains and atoms.
            cutoff: distance threshold (in Å) defining a contact.

        Returns:
            R×R binary contact map, where R is the total number of residues.
        """
        # 1. Extract heavy‐atom coords and residue indices
        # coords = []
        # resid_idx = []
        # for chain in self.structure[0]:
        #     for residue in sorted(
        #         Selection.unfold_entities(chain, "R"), key=lambda r: r.get_id()[1]
        #     ):
        #         idx = residue.get_id()[1]  # residue sequence number
        #         for atom in residue:
        #             if atom.element != "H":  # heavy atom check
        #                 coords.append(atom.coord)
        #                 resid_idx.append(idx)
        # coords = np.array(coords, dtype=float)  # shape: (H, 3)
        # resid_idx = np.array(resid_idx, dtype=int)  # length H
        coords = self.all_coords
        resid_idx = self.all_coords_ids

        # 2. Compute H×H pairwise distances and threshold
        print(coords.shape)
        print(resid_idx.shape)
        dm = distance_matrix(coords, coords, p=2)  # SciPy fast C‑loop
        mask = (dm < self.cutoff) & (dm > 0.0)  # Boolean contact mask
        print("1")
        # 3. Build one‐hot encoding (R × H)
        unique_res = np.unique(resid_idx)
        R = unique_res.size
        # Map residue sequence numbers to 0…R‑1 indices
        res_to_idx = {res: i for i, res in enumerate(unique_res)}

        print("2")
        # 3. Build one‐hot encoding (R × H)
        mapped = np.vectorize(res_to_idx.get)(resid_idx)  # length H
        one_hot = np.zeros((R, coords.shape[0]), dtype=int)
        one_hot[mapped, np.arange(coords.shape[0])] = 1  # O(H) assignment

        print("3")
        # 3. Build one‐hot encoding (R × H)
        # 4. Aggregate via matrix multiplication
        # contact_counts = one_hot @ mask.astype(int) @ one_hot.T

        # X.shape == (H, 3)
        A = radius_neighbors_graph(
            X=coords,
            radius=self.cutoff,
            mode="connectivity",
            metric="euclidean",
            include_self=False,
            n_jobs=-1,
        )
        # A is a CSR matrix of shape (H, H); next map to residues same as above
        print(A)
        plot_contact_map(A.toarray())
        sys.exit()
        contact_map = (contact_counts > 0).astype(int)  # Binarize >0

        # 5. Remove diagonal and near‐diagonal hits
        idxs = np.arange(R)
        diag_mask = np.abs(idxs[:, None] - idxs[None, :]) <= 3
        contact_map[diag_mask] = 0

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
