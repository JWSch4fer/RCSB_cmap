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

    def __init__(
        self,
        structure,
        chains_like: bool = False,
        levenshtein_cutoff: int = 30,
        cutoff: float = 8.0,
    ) -> None:
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
        self._remove_empty_chains()
        if chains_like:
            self._only_keep_chains_like(chains_like, levenshtein_cutoff)
        self._check_length()

        # for chain in self.structure.get_chains():
        #     residues = sorted(
        #         Selection.unfold_entities(chain, "R"), key=lambda r: r.get_id()[1]
        #     )
        #     print(len(residues))

        self.all_coords, self.all_coords_res_ids, self.all_coords_int_ids = (
            self._extract_coordinates()
        )

        print(self.all_coords)
        print(self.all_coords_int_ids)
        print(self.all_coords_res_ids)

    def _remove_empty_chains(self):
        """
        some RCSB entries have empty chains
        """
        detach_list = []
        for chain in self.structure[0].get_chains():
            residues = Selection.unfold_entities(chain, "R")
            if len(residues) == 0:
                detach_list.append(chain)
        if detach_list:
            for chain in detach_list:
                self.structure[0].detach_child(chain.id)

    def _only_keep_chains_like(self, selected_chain: str, levenshtein_cutoff: int):
        # remove chains that are not related enough to the target sequence
        chains = list(self.structure[0].get_chains())
        ref = chains.index(
            [chain for chain in chains if chain.full_id[-1] == selected_chain][0]
        )

        # build the AA string for each chain once
        aa_strings = []
        for chain in chains:
            residues = sorted(
                Selection.unfold_entities(chain, "R"), key=lambda r: r.get_id()[1]
            )
            aa_strings.append("".join(TRANSLATE_AA[res.resname] for res in residues))

        # calculate the levenshtein_distance to decide which chains aren't related
        leven_score = [levenshtein_distance(aa_strings[ref], aa) for aa in aa_strings]
        detach_list = [
            idx
            for idx, score in enumerate(leven_score)
            if score[0] > levenshtein_cutoff
        ]

        # detach chains
        keys = list(self.structure[0].child_dict.keys())
        for chain in detach_list:
            self.structure[0].detach_child(keys[chain])

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
            residues = sorted(
                Selection.unfold_entities(chain, "R"),
                key=lambda r: r.get_id()[1],
            )

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
                    gap_list += list(range(prev_idx + 1, curr_idx))
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

    def get_list_of_indicies(self):
        # assume `chains` is an ordered dict or list of your chains’ CA‐atom lists
        # coor_idx = [list(range(len(atom_list))) for atom_list in chains.values()]
        coor_idx = []
        chains = list(self.structure[0].get_chains())
        for chain in chains:
            residues = Selection.unfold_entities(chain, "R")
            coor_idx.append([idx for idx, _ in enumerate(residues)])

        tot_idx = []
        offset = 0
        for local_idxs in coor_idx:
            global_idxs = [i + offset for i in local_idxs]
            tot_idx.append(global_idxs)
            offset += len(local_idxs)
        return tot_idx

    def get_list_of_chain_ids(self):
        return [chain.full_id[-1] for chain in self.structure[0].get_chains()]

    def collapse_homo(self, cmap: np.ndarray, n_chains: int) -> np.ndarray:
        # n_chains = len(list(self.structure[0].get_chains()))
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
        if n_chains > 1:
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

        # helper: pad chains based on matrices and domain
        def pad_chains(groups_matrices: Tuple[set, Tuple[int, List[str]]]):
            groups, matrices = groups_matrices
            for group in groups:
                if len(group) < 2:
                    continue

                getter = operator.itemgetter(*group)
                getter_aa_strings = getter(aa_strings)
                getter_matrices = getter(matrices)
                getter_chains = getter(chains)

                ref = aa_strings.index(max(getter_aa_strings, key=len))

                for mtx_row, chain, aa_s, g in zip(
                    getter_matrices, getter_chains, getter_aa_strings, group
                ):
                    # define reference as the longest chain in the group
                    if g == ref:
                        continue  # don't edit the ref

                    # residues = sorted(
                    #     Selection.unfold_entities(chain, "R"),
                    #     key=lambda r: r.get_id()[1],
                    # )
                    # if not residues:
                    #     continue  # skip empty chains

                    # grad indices that need to be inserted

                    # deletions from reference change indexing
                    idx_deletion_offset = sum([1 for L in mtx_row[ref][1] if L == "D"])

                    idx_for_update = [
                        idx - idx_deletion_offset
                        for idx, L in enumerate(mtx_row[ref][1])
                        if L == "I"
                    ]

                    if idx_for_update:
                        ref_residues = sorted(
                            Selection.unfold_entities(chains[ref], "R"),
                            key=lambda r: r.get_id()[1],
                        )

                        # have to make sure to match indicies of chains
                        if len(idx_for_update) == 1:
                            get_resids = [
                                ref_residues[idx_for_update[0]].full_id[-1][1]
                            ]
                        else:
                            getter_resids = operator.itemgetter(*idx_for_update)
                            get_resids = getter_resids(ref_residues)
                            get_resids = [r.full_id[-1][1] for r in get_resids]

                        self.fill_gap(
                            res_list=get_resids, chain=chain, curr_res=ref_residues[0]
                        )

        pad_chains(compute_matrices())
