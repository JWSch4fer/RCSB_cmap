import numpy as np
from project.overlap import compare_contact_maps, comparison_info, overlap_definition


def test_overlap_definition_window_match_matrix():
    # overlap_definition is called by COMPARE
    # COMPARE garuntees that A is larger
    A = np.zeros((7, 7), dtype=int)
    A[1, :] = 1
    # A[:,3] = 1
    A[0, 3] = 1
    A[3, 0] = 1

    B = np.zeros((7, 7), dtype=int)
    B[:, 2] = 1
    out = overlap_definition(A, B, mtx_return=True)
    # Where there was a '1' in B:
    assert np.array_equal(out[0:3, 2], np.array([1, 1, 1]))  # matched
    assert np.array_equal(out[3:, 2], np.array([-1, -1, -1, -1]))  # matched


    comp_info = compare_contact_maps(A, B, "A", "B")
    result = np.array([[0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 1., 0., 0., 0., 0., 0.],
       [0., 1., 1., 1., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 1., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.]])
    assert np.array_equal(comp_info.cmap_comp, result)
