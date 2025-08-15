import types
import numpy as np
import builtins
import matplotlib.pyplot as plt
from pathlib import Path
import sys

import project.cli as cli_mod


class DummyContactMap:
    def __init__(
        self,
        structure,
        chain_id=None,
        cutoff=8.0,
        chains_like=None,
        levenshtein_cutoff=None,
        **kwargs,
    ):
        self.structure = structure
        self.chain_id = chain_id
        self.cutoff = cutoff
        self.chains_like = chains_like
        self.levenshtein_cutoff = levenshtein_cutoff

    def compute_residue_contact_map(self):
        # Simple 4x4 with two contacts
        cm = np.zeros((4, 4), dtype=int)
        cm[0, 1] = cm[1, 0] = 1  # intra
        cm[2, 3] = cm[3, 2] = 1  # intra
        cm[0, 2] = cm[2, 0] = 1  # inter (will be recolored by mask in utils)
        return cm

    def get_list_of_indicies(self):
        # two chains of length 2
        return [[0, 1], [2, 3]]

    def get_list_of_chain_ids(self):
        return ["A", "B"]

    def collapse_homo(self, cmap):
        # Collapse to 2x2 (intra=1, inter=2 when present)
        collapsed = np.array([[1, 1], [1, 1]], dtype=int)
        return collapsed


class DummyClient:
    def __init__(self):
        pass

    def download_pdb(self, pdb_id, chain_id=None):
        return b"", "pdb"

    def parse_structure(self, raw, fmt, target):
        return object()  # just a sentinel


def test_cli_main_smoke(monkeypatch, tmp_path):
    # Patch in our dummies
    monkeypatch.setattr(cli_mod, "RCSBClient", DummyClient)
    monkeypatch.setattr(cli_mod, "ContactMap", DummyContactMap)

    # Avoid actually rendering big figures
    saves = []

    def fake_savefig(path, *args, **kwargs):
        saves.append(Path(path))

    monkeypatch.setattr(plt, "savefig", fake_savefig, raising=False)

    argv = ["prog", "--pdb", "1abc", "--cutoff", "8.0"]
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    monkeypatch.setattr(sys, "argv", argv)

    # Should run end-to-end without exceptions and produce at least one figure save
    cli_mod.main()
    assert len(saves) >= 1
