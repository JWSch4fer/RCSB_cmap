### src/project/cli.py
import argparse
import logging
import sys
from pathlib import Path

from .rcsb import RCSBClient
from .contact_map import ContactMap

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
    args = parse_args()
    client = RCSBClient()
    print("}}}")
    pdb_input = Path(args.pdb)
    if pdb_input.is_file():
        raw, fmt = pdb_input.read_bytes(), pdb_input.suffix.lstrip(".")
        target_name = pdb_input.stem
    else:
        raw, fmt = client.download_pdb(args.pdb, chain_id=args.chain)
        target_name = args.pdb

    # print(raw, fmt, target_name)
    structure = client.parse_structure(raw, fmt, target_name)
    print("//", structure)
    cmap_obj = ContactMap(structure, cutoff=args.cutoff)
    cmap = cmap_obj.compute_contact_map()

    if args.oligomer:
        n_chains = len(structure[0].child_list)
        cmap = cmap_obj.collapse_homo(cmap, n_chains)

    # TODO: integrate visualization module for saving CSV & figures

    if args.compare:
        # Repeat for comparison mode...
        pass


if __name__ == "__main__":
    try:
        main()
    except Exception:
        logging.exception("Failed to generate contact map:")
        sys.exit(1)
