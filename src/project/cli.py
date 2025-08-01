### src/project/cli.py
import argparse
import logging
import sys
from pathlib import Path
import numpy as np

from .rcsb import RCSBClient
from .contact_map import ContactMap
from .utils import plot_contact_map, create_oligomer_mask, create_contact_map_plot

# Configure logging
logging.basicConfig(
    filename="cmap.log",
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="rcsb_cmap",
        description="Generate residue–residue contact maps from PDB/mmCIF structures.",
    )
    parser.add_argument(
        "-p", "--pdb", required=True, help="PDB ID or local file path (.pdb or .cif)"
    )
    parser.add_argument("-c", "--chain", help="Chain ID (e.g. A)")
    parser.add_argument(
        "-p2", "--pdb2", help="PDB ID or local file path (.pdb or .cif)"
    )
    parser.add_argument("-c2", "--chain2", help="Chain ID (e.g. A)")
    parser.add_argument(
        "-o",
        "--oligomer",
        action="store_true",
        help="Collapse homo‑oligomer contacts",
    )
    parser.add_argument(
        "--cutoff", type=float, default=8.0, help="Distance cutoff in Å for contacts"
    )
    parser.add_argument(
        "--chains_like",
        help="retain only chains that are similar to the selected chain. Useful for hetero-oligomer (example: C)",
    )
    parser.add_argument(
        "--levenshtein",
        default=30,
        type=int,
        help="chains_like calculates the levenshtein distance between chains and retains chains that are within 30 of the adjust if this is to restrictive/permissive",
    )

    return parser.parse_args()


def main() -> None:

    # grab cli information
    args = parse_args()

    # read information about protein from file or RCSB.org
    client = RCSBClient()
    pdb_input = Path(args.pdb)
    if pdb_input.is_file():
        raw, fmt = pdb_input.read_bytes(), pdb_input.suffix.lstrip(".")
        target_name = pdb_input.stem
    else:
        raw, fmt = client.download_pdb(args.pdb, chain_id=args.chain)
        target_name = args.pdb

    # take the raw protein information and clean it for visualization
    structure = client.parse_structure(raw, fmt, target_name)
    cmap_obj = ContactMap(
        structure,
        cutoff=args.cutoff,
        chains_like=args.chains_like,
        levenshtein_cutoff=args.levenshtein,
    )
    cmap = cmap_obj.compute_residue_contact_map()
    mask = create_oligomer_mask(
        contact_map=cmap, chains_idx=cmap_obj.get_list_of_indicies()
    )
    plot_contact_map(cmap, mask)
    create_contact_map_plot(
        contact_map=cmap,
        chain_labels=cmap_obj.get_list_of_chain_ids(),
        chain_midpoints=[np.mean(idxs) for idxs in cmap_obj.get_list_of_indicies()],
        mask=mask,
        name=args.pdb,
    )

    if args.oligomer:
        cmap_collapsed = cmap_obj.collapse_homo(cmap)
        plot_contact_map(cmap_collapsed)

    # TODO: integrate visualization module for saving CSV & figures

    # if args.compare:
    #     # Repeat for comparison mode...
    #     pass


if __name__ == "__main__":
    try:
        main()
    except Exception:
        logging.exception("Failed to generate contact map:")
        sys.exit(1)
