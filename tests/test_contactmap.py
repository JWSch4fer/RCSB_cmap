### tests/test_contact_map.py
import numpy as np
from project.contact_map import ContactMap


def test_compute_contact_map_simple():
    # Create dummy ContactMap with preset coordinates
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ]
    )
    cmap_obj = ContactMap.__new__(ContactMap)
    cmap_obj.all_coords = coords
    cmap_obj.all_coords_int_ids = np.array([1, 2, 3, 4, 5])
    cmap_obj.cutoff = 1.5

    cmap = cmap_obj.compute_residue_contact_map()
    # Expect contact between 0-1 only
    expected = np.array(
        [
            [0, 0, 0, 0, 1],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
        ]
    )
    assert np.array_equal(cmap, expected)


def test_collapse_homo():
    # Simulate a 2-chain homo-oligomer map of length 2 each
    cmap = np.array(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0],
        ]
    )
    cmap_obj = ContactMap.__new__(ContactMap)
    collapsed = cmap_obj.collapse_homo(cmap, n_chains=2)
    # Intrachain contacts become 1, interchain become 2
    assert collapsed.shape == (2, 2)
    # All intrachain summed: chain length=2, so indices 0-1 and 2-3
    assert np.array_equal(collapsed, np.array([[0, 1], [1, 0]]))
