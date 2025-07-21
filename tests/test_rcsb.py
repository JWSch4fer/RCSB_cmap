### tests/test_rcsb.py
import io
import pytest
from project.rcsb import RCSBClient

PDB_SAMPLE = b"""
ATOM      1  N   ALA A   1      11.104  13.207   9.923  1.00 20.00           N  
ATOM      2  CA  ALA A   1      12.560  13.207   9.923  1.00 20.00           C  
ATOM      3  C   ALA A   1      12.904  14.640   9.523  1.00 20.00           C  
ATOM      4  O   ALA A   1      12.104  15.540   9.223  1.00 20.00           O  
TER
END
"""


def test_parse_structure_pdb():
    client = RCSBClient()
    # Parse in-memory PDB bytes
    structure = client.parse_structure(PDB_SAMPLE, fmt="pdb", target="test")
    models = list(structure.get_models())
    assert len(models) == 1
    chains = list(models[0].get_chains())
    assert len(chains) == 1
    residues = list(chains[0].get_residues())
    assert residues[0].get_resname() == "ALA"


@pytest.mark.parametrize(
    "chain_id,expected_prefix",
    [
        (None, b"ATOM"),
        ("A", b"ATOM"),
    ],
)
def test_download_pdb_fallback(monkeypatch, chain_id, expected_prefix):
    # Monkeypatch session.get to simulate a valid PDB download
    client = RCSBClient()

    class DummyResp:
        def __init__(self, content):
            self.content = content
            self.ok = True

        def raise_for_status(self):
            pass

    monkeypatch.setattr(
        client.session, "get", lambda url, params=None: DummyResp(PDB_SAMPLE)
    )
    data, fmt = client.download_pdb("xxxx", 1, chain_id)
    assert data.startswith(expected_prefix)
    assert fmt == "pdb"
