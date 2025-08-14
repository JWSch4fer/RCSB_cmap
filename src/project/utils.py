### src/project/utils.py
import numpy as np
import pandas as pd
from typing import Tuple, List, Sequence, Optional
import itertools
import matplotlib.pyplot as plt

from .overlap import comparison_info, overlap_definition


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


def plot_contact_map(contact_map: np.ndarray, mask: Optional[np.ndarray] = []) -> None:
    """
    Display a residue–residue contact map.

    Args:
        contact_map: 2D binary NumPy array (shape: R×R).
    """

    f, ax = plt.subplots(1, 1, figsize=(9, 9))
    axes = []

    # Render the contact map as an image
    plt.imshow(
        contact_map,
        cmap="viridis",  # perceptually uniform colormap
        interpolation="none",  # no smoothing between cells
        origin="lower",  # put [0,0] in the bottom-left
    )
    if len(mask) > 0:
        axes.append(ax.imshow(mask, cmap="Greys", interpolation="none", alpha=0.11))
    plt.xlabel("Residue index")
    plt.ylabel("Residue index")
    plt.title("Protein Contact Map")
    plt.tight_layout()
    plt.savefig("temp.png")
    plt.clf()


def create_contact_map_plot(
    contact_map: np.ndarray,
    chain_labels: list[str],
    chain_midpoints: list[int],
    mask: Optional[np.ndarray] = np.ndarray([]),
    name: Optional[str] = "prot",
    comp: Optional[comparison_info] = None,
) -> None:
    """
    Create pandas data frame of contact map for visualizing cmap
    """

    low = np.tril(contact_map)
    low[np.triu_indices(low.shape[0])] = low.T[np.triu_indices(low.shape[0])]

    up = np.triu(contact_map)
    up[np.tril_indices(up.shape[0])] = up.T[np.tril_indices(up.shape[0])]

    low = overlap_definition(up, low, mtx_return=True)
    up = overlap_definition(low, up, mtx_return=True)
    colors = {
        "INTRA_COMMON": "#383838",
        "INTRA_UNIQUE": "#377c2b",
        "INTER_COMMON": "#7b7b7b",
        "INTER_UNIQUE": "#b9dcb3",
    }  # green version

    # check for all contacts to find matches to ACE contacts
    contact_types = {}
    contact_types["low_true_intra"] = np.argwhere(
        (np.tril(low) == 1) | (np.tril(low) == 3)
    ).tolist()
    contact_types["low_false_intra"] = np.argwhere(
        (np.tril(low) == -1) | (np.tril(low) == -3)
    ).tolist()
    contact_types["low_true_inter"] = np.argwhere((np.tril(low) == 2)).tolist()
    contact_types["low_false_inter"] = np.argwhere((np.tril(low) == -2)).tolist()
    contact_types["up_true_intra"] = np.argwhere(
        (np.triu(up) == 1) | (np.triu(up) == 3)
    ).tolist()
    contact_types["up_false_intra"] = np.argwhere(
        (np.triu(up) == -1) | (np.triu(up) == -3)
    ).tolist()
    contact_types["up_true_inter"] = np.argwhere((np.triu(up) == 2)).tolist()
    contact_types["up_false_inter"] = np.argwhere((np.triu(up) == -2)).tolist()

    df = {"i": [], "j": [], "type": []}
    for key in contact_types.keys():
        for i, j in contact_types[key]:
            if key == "low_true_intra":
                df["i"].append(i)
                df["j"].append(j)
                df["type"].append("INTRA_COMMON")
            if key == "low_false_intra":
                df["i"].append(i)
                df["j"].append(j)
                df["type"].append("INTRA_UNIQUE")
            if key == "low_true_inter":
                df["i"].append(i)
                df["j"].append(j)
                df["type"].append("INTER_COMMON")
            if key == "low_false_inter":
                df["i"].append(i)
                df["j"].append(j)
                df["type"].append("INTER_UNIQUE")
            if key == "up_true_intra":
                df["i"].append(i)
                df["j"].append(j)
                df["type"].append("INTRA_COMMON")
            if key == "up_false_intra":
                df["i"].append(i)
                df["j"].append(j)
                df["type"].append("INTRA_UNIQUE")
            if key == "up_true_inter":
                df["i"].append(i)
                df["j"].append(j)
                df["type"].append("INTER_COMMON")
            if key == "up_false_inter":
                df["i"].append(i)
                df["j"].append(j)
                df["type"].append("INTER_UNIQUE")

    df = pd.DataFrame.from_dict(df)

    f, ax = plt.subplots(1, 1, figsize=(9, 9), dpi=300)
    axes = []

    # first, if you have a mask, draw it as a background in axes coords
    if len(mask) > 0 and "collapse" not in name:
        ax.imshow(
            mask,
            cmap="Greys",
            interpolation="none",
            alpha=0.11,
            origin="lower",
            zorder=0,  # behind scatter
            clip_on=False,  # ensure nothing clips it at the data border
        )

    # predefine limits for marker size calculations
    ax.set_xlim(-0.5, contact_map.shape[0] + 0.5)
    ax.set_ylim(-0.5, contact_map.shape[0] + 0.5)
    # standardize size of markers no matter the size of the protein
    r_0 = cell_scatter_size(ax)  # radius of markers

    # each 'type' needs its own handle for matplotlib to give unique legend elements
    for t in df["type"].unique():
        axes.append(
            ax.scatter(
                x=df.loc[df["type"].eq(t), "i"],
                y=df.loc[df["type"].eq(t), "j"],
                c=df.loc[df["type"].eq(t), "type"].map(colors),
                s=r_0,
                linewidth=0,
                linestyle="None",
            )
        )

    # ax.set_ylabel(y_name)
    ax.set_xlabel(name.replace(".pdb", ""))
    if comp:
        ax.set_ylabel(comp.y_name.replace(".pdb", ""))
        ax.set_xlabel(comp.x_name.replace(".pdb", ""))

    if "collapse" not in name and comp == None:
        # chain names at midpoints
        ax.set_yticks(chain_midpoints)
        ax.set_yticklabels(chain_labels)

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height * 0.85])

    # Put a legend to the right of the current axis
    lgnd = ax.legend(df["type"].unique(), loc="center left", bbox_to_anchor=(1, 0.5))
    for handle in lgnd.legend_handles:
        handle.set_sizes([15.0])

    # name = name.replace(re.search(r"(.cif|.pdb)", name).group(0), '') if bool(re.search(r"(.cif|.pdb)", name)) else name
    # print("Writing csv...", f"{name}_df.csv")
    # df.to_csv(f'{name}_df.csv')
    print("Writing image...", "{:}.png".format(name.replace(".pdb", "")))
    plt.savefig("{:}.png".format(name.replace(".pdb", "")))
    plt.clf()
    plt.close()


def cell_scatter_size(ax, frac=1.0, units_per_cell=1.0):
    """
    Compute matplotlib.scatter(s=...) so a marker's DIAMETER is `frac` of one grid cell.
    This standardizes the size of a point no matter the area of the image

    NOTE:
        1 point = 1/72 inch (the PostScript point that matplotlib uses).

        inches = pixels / dpi
        points = inches × 72

        Combined: points = pixels * 72 / dpi

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes the scatter will be drawn on (must have correct xlim/ylim set).
    frac : float in (0, 1]
        Desired marker diameter as a fraction of one cell (e.g., 0.8 = 80% of a cell).
    units_per_cell : float
        Size of a cell in data units (1.0 for typical contact maps on integer grids).

    Returns
    -------
    s : float
        Area in points^2 to pass to scatter(..., s=s).
    """
    # axes pixel size
    fig = ax.figure
    dpi = fig.dpi
    fig_w_in, fig_h_in = fig.get_size_inches()
    pos = ax.get_position()  # axes bbox in figure-relative coords (0..1)
    ax_w_px = fig_w_in * pos.width * dpi
    ax_h_px = fig_h_in * pos.height * dpi

    # Data span (works even if axes are inverted)
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    span_x = abs(x1 - x0)
    span_y = abs(y1 - y0)

    # Pixels per data unit along each axis
    px_per_unit_x = ax_w_px / span_x
    px_per_unit_y = ax_h_px / span_y

    # One cell (units_per_cell wide/tall) in pixels; should always be a square for proteins...
    cell_px = min(px_per_unit_x, px_per_unit_y) * units_per_cell

    # Desired marker diameter in pixels NOTE: 1 = entire cell
    d_px = frac * cell_px

    # Convert diameter in pixels -> points, then area (pt^2) for scatter's 's'
    d_pt = d_px * 72.0 / dpi
    s = np.pi * (d_pt / 2.0) ** 2
    return max(s, 1e-6)  # guard against zero if something degenerate happens


def create_oligomer_mask(
    contact_map: np.ndarray,
    chains_idx: Sequence[Sequence[int]],
    value_step: float = 0.2,
) -> np.ndarray:
    """
    Build an inter‑chain mask for a homo‑oligomer contact map.

    Args:
        contact_map:        square (N×N) contact map array (modified in place).
        chains_idx:         sequence of index lists, one per chain,
                            where each list contains the flat indices for that chain.
        value_step:         how much darker each successive inter‑chain block is.

    Returns:
        mask:               N×N float array with inter‑chain shading levels.
                            Also mutates `contact_map`, marking any existing
                            contacts in these regions with a value of 2.
    """
    n = contact_map.shape[0]
    mask = np.zeros((n, n), dtype=float)

    # for each unordered pair of chains (i < j)
    for i, j in itertools.combinations(range(len(chains_idx)), 2):
        rows = chains_idx[i]
        cols = chains_idx[j]
        level = 1.0 - (j - i) * value_step

        # Shade the off‑diagonal block [rows, cols] and its transpose
        mask[np.ix_(rows, cols)] = level
        mask[np.ix_(cols, rows)] = level

        # Wherever there was a contact (==1) in the original map, recolor to 2
        block = contact_map[np.ix_(rows, cols)]
        hits = block == 1
        if np.any(hits):
            # get the positions of those hits within the block
            hit_rows, hit_cols = np.nonzero(hits)
            # map back to global indices
            global_rows = np.array(rows)[hit_rows]
            global_cols = np.array(cols)[hit_cols]
            contact_map[global_rows, global_cols] = 2
            contact_map[global_cols, global_rows] = 2

    return mask
