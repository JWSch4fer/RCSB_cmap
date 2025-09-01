### src/project/cli.py
import argparse
from dataclasses import dataclass
import logging
import sys
from pathlib import Path
import numpy as np

from .rcsb import RCSBClient
from .contact_map import ContactMap
from .utils import create_oligomer_mask, create_contact_map_plot
from .overlap import compare_contact_maps

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
        "-p2", "--pdb2", default="", help="PDB ID or local file path (.pdb or .cif)"
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
    create_contact_map_plot(
        contact_map=cmap,
        chain_labels=cmap_obj.get_list_of_chain_ids(),
        chain_midpoints=[np.mean(idxs) for idxs in cmap_obj.get_list_of_indicies()],
        mask=mask,
        name=args.pdb,
    )

    if args.oligomer:
        cmap_collapsed = cmap_obj.collapse_homo(
            cmap, len(cmap_obj.get_list_of_chain_ids())
        )
        create_contact_map_plot(
            contact_map=cmap_collapsed,
            chain_labels=cmap_obj.get_list_of_chain_ids(),
            chain_midpoints=[np.mean(idxs) for idxs in cmap_obj.get_list_of_indicies()],
            mask=mask,
            name=args.pdb + "-collapsed",
        )

    if args.pdb2:
        # Repeat for comparison mode...
        # read information about protein from file or RCSB.org
        client_2 = RCSBClient()
        pdb_input_2 = Path(args.pdb2)
        if pdb_input_2.is_file():
            raw_2, fmt_2 = pdb_input_2.read_bytes(), pdb_input_2.suffix.lstrip(".")
            target_name_2 = pdb_input_2.stem
        else:
            raw_2, fmt_2 = client_2.download_pdb(args.pdb2, chain_id=args.chain2)
            target_name_2 = args.pdb2

        # take the raw protein information and clean it for visualization
        structure_2 = client_2.parse_structure(raw_2, fmt_2, target_name_2)
        cmap_obj_2 = ContactMap(
            structure_2,
            cutoff=args.cutoff,
            chains_like=args.chains_like,
            levenshtein_cutoff=args.levenshtein,
        )
        cmap_2 = cmap_obj_2.compute_residue_contact_map()
        mask_2 = create_oligomer_mask(
            contact_map=cmap_2, chains_idx=cmap_obj_2.get_list_of_indicies()
        )

        create_contact_map_plot(
            contact_map=cmap_2,
            chain_labels=cmap_obj_2.get_list_of_chain_ids(),
            chain_midpoints=[
                np.mean(idxs) for idxs in cmap_obj_2.get_list_of_indicies()
            ],
            mask=mask_2,
            name=args.pdb2,
        )

        comp_info = compare_contact_maps(cmap, cmap_2, args.pdb, args.pdb2)
        create_contact_map_plot(
            contact_map=comp_info.cmap_comp,
            chain_labels=cmap_obj_2.get_list_of_chain_ids(),
            chain_midpoints=[
                np.mean(idxs) for idxs in cmap_obj_2.get_list_of_indicies()
            ],
            mask=np.pad(
                mask_2,
                ((comp_info.padding_offset, 0), (comp_info.padding_offset, 0)),
                mode="edge",
            ),
            name=args.pdb + "-" + args.pdb2 + "-comparison",
        )

        if args.oligomer:
            cmap_collapsed_2 = cmap_obj_2.collapse_homo(
                cmap_2, len(cmap_obj_2.get_list_of_chain_ids())
            )
            create_contact_map_plot(
                contact_map=cmap_collapsed_2,
                chain_labels=cmap_obj_2.get_list_of_chain_ids(),
                chain_midpoints=[
                    np.mean(idxs) for idxs in cmap_obj_2.get_list_of_indicies()
                ],
                mask=mask_2,
                name=args.pdb2 + "-collapsed",
            )

            comp_info_collapsed = compare_contact_maps(
                cmap_collapsed, cmap_collapsed_2, args.pdb, args.pdb2
            )

            create_contact_map_plot(
                contact_map=comp_info_collapsed.cmap_comp,
                chain_labels=cmap_obj.get_list_of_chain_ids(),
                chain_midpoints=[
                    np.mean(idxs) for idxs in cmap_obj.get_list_of_indicies()
                ],
                mask=mask_2,
                name=args.pdb + "-" + args.pdb2 + "-comparison" + "-collapsed",
            )


if __name__ == "__main__":
    try:
        main()
    except Exception:
        logging.exception("Failed to generate contact map:")
        sys.exit(1)
