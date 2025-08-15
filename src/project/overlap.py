### src/project/overlap.py
import numpy as np
from dataclasses import dataclass


def overlap_definition(A, B, mtx_return=False):
    # pad edges for sliding window consistency
    padded_A = np.pad(A, ((1, 1), (1, 1)), mode="constant", constant_values=0)
    # Extract all possible windows
    windows = np.lib.stride_tricks.sliding_window_view(padded_A, (3, 3))
    # positions with contacts
    mask_1 = np.where(B == 1)
    mask_2 = np.where(B == 2)
    mask_3 = np.where(B == 3)
    # find number of contacts in B that match a (3,3) window in A that has a nonzero elements
    matches_1 = np.any(windows[mask_1] != 0, axis=(-1, -2))
    matches_2 = np.any(windows[mask_2] != 0, axis=(-1, -2))
    matches_3 = np.any(windows[mask_3] != 0, axis=(-1, -2))

    # return B array with +s where a +-1 match in A exists and -s where they don't exist
    # 0s where no contact exists, this function retains the type of contact
    if mtx_return == True:
        B_true = np.copy(B)
        B_true[mask_1] = matches_1 * 1
        B_true[mask_2] = matches_2 * 2
        B_true[mask_3] = matches_3 * 3

        B_false = np.copy(B)
        B_false[mask_1] = ~matches_1 * 1
        B_false[mask_2] = ~matches_2 * 2
        B_false[mask_3] = ~matches_3 * 3
        B_false = B_false * -1
        return B_true + B_false
    return sum([*matches_1, *matches_2, *matches_3])


@dataclass
class comparison_info:
    cmap_comp: np.ndarray
    x_name: str
    padding_offset: int
    y_name: str
    align_offset: int
    max_overlap: int


def compare_contact_maps(cmap1, cmap2, pdb1_name, pdb2_name) -> comparison_info:
    """
    Compare two contact maps and return a contact map alignment based on number of symmetric contacts
    """
    A = cmap1
    A_name = pdb1_name
    B = cmap2
    B_name = pdb2_name

    if B.shape[0] < A.shape[0]:
        B = np.pad(
            B, ((0, A.shape[0] - B.shape[0]), (0, A.shape[0] - B.shape[0])), "constant"
        )

    # pad first matrix to allow for sliding window
    n = max(1, int(np.ceil(B.shape[0] * 0.1)))
    B = np.pad(B, ((n, n), (n, n)), mode="constant", constant_values=0)

    mask = np.triu(np.ones(B.shape))
    B = B * mask

    # initialize the maximum number of overlapping 1s
    max_overlap = 0

    # initialize the starting idices of the optimal alignment
    start_i, start_j = 0, 0

    beg = range(B.shape[0] - A.shape[0] + 1)
    end = range(B.shape[0] - A.shape[0] + 1)

    for i, j in zip(beg, end):
        subset = B.copy()
        subset = subset[i : i + A.shape[0], j : j + A.shape[1]]
        # only consider exact mathches
        # overlap = len(np.argwhere(np.multiply(A, subset) == 1))

        # consider matches within +- 1 | good luck
        overlap = overlap_definition(A, subset)
        # update the maximum number of overlapping 1s and the best starting indices if needed
        if overlap > max_overlap:
            max_overlap = overlap
            start_i, start_j = i, j

    # print the maximum number of overlapping 1s and the starting indicies of the optimal alignment
    print(f"Maximum number of overlapping 1s: {max_overlap}")
    print(f"Starting indices of the optimal alignment: ({start_i}, {start_j})")

    lower_triangle_idx = np.tril_indices(A.shape[0], 0)

    # offset indices to create the aligned dualfold cmap
    lower_triangle_idx_align = np.transpose(lower_triangle_idx)
    lower_triangle_idx_align = lower_triangle_idx_align + start_i
    lower_triangle_idx_align = np.transpose(lower_triangle_idx_align)
    lower_triangle_idx_align = tuple(np.array(i) for i in lower_triangle_idx_align)

    # create the duafold cmap
    B[lower_triangle_idx_align] = A[lower_triangle_idx]

    # return B, B_name, n, A_name, start_i, max_overlap
    return comparison_info(
        cmap_comp=B,
        x_name=B_name,
        padding_offset=n,
        y_name=A_name,
        align_offset=start_i,
        max_overlap=max_overlap,
    )
