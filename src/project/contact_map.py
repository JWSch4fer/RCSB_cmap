### src/project/contact_map.py
import numpy as np
from sklearn.neighbors import radius_neighbors_graph
from sklearn.preprocessing import LabelBinarizer
from scipy.sparse import csr_matrix
import itertools, operator
from typing import Tuple, Optional, List

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
            all_coords: every heavy atom x,y,z
            all_coords_res_ids: the chain id:res id
            all_coords_int_ids: integer index from 0 to # of residues
        """
        self.structure = structure
        self.cutoff = cutoff
        self._process_protein()
        self._check_length()
        for chain in self.structure.get_chains():
            residues = sorted(
                Selection.unfold_entities(chain, "R"), key=lambda r: r.get_id()[1]
            )
            print(len(residues))

        self.all_coords, self.all_coords_res_ids, self.all_coords_int_ids = (
            self._extract_coordinates()
        )

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

    def _extract_coordinates(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        # some RCSB enteries are poorly written, we'll deal with empty chains here
        detach_list = []

        for chain in self.structure[0].get_chains():
            residues = list(Selection.unfold_entities(chain, "R"))
            if not residues:
                detach_list.append(chain)
                continue

            # Iterate over each consecutive pair of residues
            gap_list = []
            copy_res = residues[0]
            for prev_res, curr_res in itertools.pairwise(residues):
                prev_idx = prev_res.get_id()[1]
                curr_idx = curr_res.get_id()[1]
                if prev_idx + 1 != curr_idx:
                    gap_list += list(range(prev_idx, curr_idx))
            self.fill_gap(
                res_list=gap_list,
                chain=chain,
                curr_res=copy_res,
            )

        # Remove any chains that had zero residues
        for ch in detach_list:
            self.structure[0].detach_child(ch.id)

        _temp_id = 0
        all_coords = []
        for chain in self.structure.get_chains():
            for residue in sorted(
                Selection.unfold_entities(chain, "R"), key=lambda r: r.get_id()[1]
            ):
                for atom in residue.get_list():
                    all_coords.append(
                        (
                            atom.coord,
                            chain._id + ":" + str(residue.get_id()[1]),
                            _temp_id,
                        )
                    )
                _temp_id += 1

        # return coordinates and identifiers
        _ = np.array(all_coords, dtype=object)
        all_coords = np.stack(_[:, 0], axis=0)
        all_coords_res_ids = np.stack(_[:, 1], axis=0)
        all_coords_int_ids = np.stack(_[:, 2], axis=0)
        return all_coords, all_coords_res_ids, all_coords_int_ids

    def fill_gap(self, res_list: list[int], chain, curr_res) -> None:
        # Generate missing indices lazily
        for missing_idx in res_list:
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
        resid_idx = self.all_coords_int_ids

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
        lb = LabelBinarizer()
        lb.fit(resid_idx)
        one_hot = lb.transform(resid_idx)
        one_hot_sparse = csr_matrix(one_hot)

        # Aggregate via matrix multiplication
        # this step must be sparse matrix operation

        contact_counts = one_hot_sparse.T @ A @ one_hot_sparse

        # contact_counts = one_hot_sparse @ A @ one_hot_sparse.T

        contact_map = (contact_counts > 0).astype(int)  # Binarize >0
        contact_map = contact_map.toarray()

        # Remove diagonal and near‐diagonal hits
        idxs = np.arange(contact_map.shape[0])
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

    def _check_length(self):
        """
        Standardize homo‑oligomer chain lengths by padding missing residues
        at the N‑terminus (NTD) and C‑terminus (CTD) based on Levenshtein alignments.
        """
        # get list of chains
        chains = list(self.structure[0].get_chains())
        if len(chains) < 2:
            return

        # build the AA string for each chain once
        aa_strings = []
        for chain in chains:
            residues = sorted(
                Selection.unfold_entities(chain, "R"), key=lambda r: r.get_id()[1]
            )
            aa_strings.append("".join(TRANSLATE_AA[res.resname] for res in residues))
        for a in aa_strings:
            print(len(a))

        # # helper: compute matrices for a given domain
        def compute_matrices() -> Tuple[set, Tuple[int, List[str]]]:
            matrices = []
            groups = set()
            for ref in aa_strings:
                row = [levenshtein_distance(ref, aa) for aa in aa_strings]
                g = tuple(idx_r for idx_r, r in enumerate(row) if r[0] <= 30)
                groups.add(g)
                matrices.append(row)
            return groups, matrices

        #     if domain == "NTD":
        #         # for NTD, stop at first insertion within first 10
        #         matrices = []
        #         for ref in aa_strings:
        #             row = [levenshtein_distance(ref, aa) for aa in aa_strings]
        #             matrices.append(row)
        #             # if any insertion ('I') in the first 10 ops of this alignment, break
        #             if any("I" in ops[-1][:10] for _, ops in row):
        #                 break
        #         return matrices

        #     elif domain == "CTD":
        #         # for CTD, align everything to the longest sequence
        #         ref = max(aa_strings, key=len)
        #         return [[levenshtein_distance(ref, aa) for aa in aa_strings]]

        #     else:
        #         return []

        # helper: pad chains based on matrices and domain
        def pad_chains(groups_matrices: Tuple[set, Tuple[int, List[str]]]):
            groups, matrices = groups_matrices
            for group in groups:
                getter = operator.itemgetter(*group)
                for mtx_row, chain in zip(getter(matrices), getter(chains)):
                    # boundary residue (first for NTD, last for CTD)
                    residues = sorted(
                        Selection.unfold_entities(chain, "R"),
                        key=lambda r: r.get_id()[1],
                    )
                    if not residues:
                        continue


"""
groups are a good idea to reduce useless computation
still need to use longest in a group as reference
This should be the only member of the group that has to be iterated fully
"""
                    gap_res = residues[0]
                    missing = 

                    print(mtx_row)

                sys.exit(1)
            for mtx_row, chain in zip(matrices, chains):
                # boundary residue (first for NTD, last for CTD)
                residues = sorted(
                    Selection.unfold_entities(chain, "R"), key=lambda r: r.get_id()[1]
                )
                if not residues:
                    continue

                if domain == "NTD":
                    for leven_score, mem in mtx_row:
                        if leven_score <= 30:
                            boundary = residues[0]
                            end_idx = boundary.get_id()[1]
                            ops = np.array(list(mtx_row[-1][-1]), dtype="<U1")
                            start_idx = end_idx - np.argmax(ops == "M")

                            # print("ntd", mem)
                            # print([idx for idx, L in enumerate(mem) if L == "D"])

                    missing = list(range(start_idx, end_idx))

                else:  # CTD
                    for leven_score, mem in mtx_row:
                        if leven_score <= 30:

                            boundary = residues[-1]
                            start_idx = boundary.get_id()[1] + 1
                            ops = np.array(list(mtx_row[-1][-1]), dtype="<U1")
                            end_idx = start_idx + np.count_nonzero(ops == "D")
                            # print("ctd", mem)
                            # print([idx for idx, L in enumerate(mem) if L == "I"])
                    missing = list(range(start_idx, end_idx))

                self.fill_gap(res_list=missing, chain=chain, curr_res=boundary)

        # pad NTD then CTD
        pad_chains(compute_matrices())
        # pad_chains(compute_matrices("CTD"), "CTD")
