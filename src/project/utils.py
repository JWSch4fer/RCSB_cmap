### src/project/utils.py
import numpy as np
from typing import Tuple, List
import matplotlib.pyplot as plt


def levenshtein_distance(s1: str, s2: str) -> Tuple[int, List[str]]:
    """
    Compute Levenshtein distance and edit operations from s2 → s1.

    Returns:
        (distance, list of edit ops).
    """
    len1, len2 = len(s1), len(s2)
    dp = [[(0, []) for _ in range(len2 + 1)] for _ in range(len1 + 1)]

    for i in range(1, len1 + 1):
        dp[i][0] = (i, ["D"] * i)
    for j in range(1, len2 + 1):
        dp[0][j] = (j, ["I"] * j)

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            cost, op = (0, "M") if s1[i - 1] == s2[j - 1] else (1, "S")
            d_cost, d_ops = dp[i - 1][j]
            i_cost, i_ops = dp[i][j - 1]
            s_cost, s_ops = dp[i - 1][j - 1]
            choices = [
                (d_cost + 1, d_ops + ["D"]),
                (i_cost + 1, i_ops + ["I"]),
                (s_cost + cost, s_ops + [op]),
            ]
            dp[i][j] = min(choices, key=lambda x: x[0])

    dist, ops = dp[len1][len2]
    return dist, ops


def pad_with(vector: np.ndarray, pad_width: Tuple[int, int], iaxis: int, kwargs):
    """
    Custom padder: fill pad edges with a constant.
    Used as callback for np.pad(mode='constant').
    """
    pad_value = kwargs.get("padder", 0)
    vector[: pad_width[0]] = pad_value
    vector[-pad_width[1] :] = pad_value


def plot_contact_map(contact_map: np.ndarray) -> None:
    """
    Display a residue–residue contact map.

    Args:
        contact_map: 2D binary NumPy array (shape: R×R).
    """
    plt.figure(figsize=(8, 8))
    # Render the contact map as an image
    img = plt.imshow(
        contact_map,
        cmap="viridis",  # perceptually uniform colormap
        interpolation="none",  # no smoothing between cells
        origin="lower",  # put [0,0] in the bottom-left
    )
    plt.colorbar(
        img, label="Contact presence"
    )  # show scale :contentReference[oaicite:2]{index=2}
    plt.xlabel("Residue index")
    plt.ylabel("Residue index")
    plt.title("Protein Contact Map")
    plt.tight_layout()
    plt.savefig("temp.png")
    plt.clf()
