
import numpy as np
import matplotlib.pyplot as plt

from project.utils import (
    levenshtein_distance,
    cell_scatter_size,
    create_oligomer_mask,
)


def test_levenshtein_distance_basic():
    # classic example; we only assert the distance (ops can vary)
    dist, ops = levenshtein_distance("kitten", "sitting")
    assert dist == 3
    assert isinstance(ops, list) and len(ops) >= dist


def test_cell_scatter_size_positive_and_scales():
    fig, ax = plt.subplots(figsize=(4, 4), dpi=100)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    s1 = cell_scatter_size(ax, units_per_cell=1.0)
    assert s1 > 0

    # Bigger axes (more pixels per unit) should increase marker area
    fig2, ax2 = plt.subplots(figsize=(8, 8), dpi=100)
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 10)
    s2 = cell_scatter_size(ax2, units_per_cell=1.0)
    assert s2 > s1


def test_create_oligomer_mask_and_map_recoloring():
    # 2 chains of length 2 each → indices [0,1] and [2,3]
    chains_idx = [[0, 1], [2, 3]]
    contact_map = np.zeros((4, 4), dtype=int)

    # Put inter-chain contacts in the off-diagonal blocks
    contact_map[0, 2] = 1
    contact_map[2, 0] = 1
    contact_map[1, 3] = 1
    contact_map[3, 1] = 1

    mask = create_oligomer_mask(contact_map, chains_idx, value_step=0.2)

    # Mask should shade off-diagonal blocks and be symmetric
    assert mask.shape == (4, 4)
    assert np.allclose(mask[:2, :2], 0)  # intra-block unshaded
    assert np.all(mask[:2, 2:] > 0) and np.all(mask[2:, :2] > 0)

    # Function recolors *existing* inter-chain contacts to 2 in the map
    assert contact_map[0, 2] == 2 and contact_map[2, 0] == 2
    assert contact_map[1, 3] == 2 and contact_map[3, 1] == 2
