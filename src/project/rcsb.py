### src/project/rcsb.py
from pathlib import Path
import io, json
import logging
from typing import Optional, Tuple, Union
import tempfile, gzip

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from Bio.PDB import PDBParser, MMCIFParser

# Suppress Biopython warnings to WARNING level
logging.getLogger("Bio.PDB").setLevel(logging.WARNING)


class RCSBClient:
    """Client to download and parse PDB/mmCIF structures from RCSB."""

    BASE_PDB_URL = "https://www.rcsb.org/pdb/download/downloadFile.do"
    BASE_CIF_URL = "https://files.rcsb.org/pub/pdb/data/assemblies/mmCIF/divided"

    def __init__(self, retries: int = 5, backoff_factor: float = 0.5) -> None:
        """
        Initialize a requests.Session with retry logic.

        Args:
            retries: number of retry attempts for transient errors.
            backoff_factor: sleep factor between retries.
        """
        self.session = requests.Session()
        retry_strategy = Retry(
            total=retries,
            backoff_factor=backoff_factor,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET"],
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

    def download_pdb(
        self, pdb_id: str, assembly_id: int = 1, chain_id: Optional[str] = None
    ) -> Tuple[bytes, str]:
        """
        Download a PDB or mmCIF file for the given PDB ID.

        Args:
            pdb_id: 4-character PDB identifier.
            assembly_id: assembly number (usually 1).
            chain_id: optional chain filter for PDB format.

        Returns:
            Tuple of raw content bytes and format string ('pdb' or 'cif').

        Raises:
            HTTPError: if neither download succeeds.
        """
        pdb_id = pdb_id.lower()

        # Fallback to mmCIF (.cif.gz)
        cif_path = f"{pdb_id[1:3]}/{pdb_id}-assembly{assembly_id}.cif.gz"
        url = f"{self.BASE_CIF_URL}/{cif_path}"
        response = self.session.get(url)
        response.raise_for_status()
        logging.info(f"mmCIF download successful for {pdb_id}")

        if response.ok:
            return response.content, "cif"

        # Try PDB format (allows chain filtering)
        params = {
            "fileFormat": "pdb",
            "compression": "NO",
            "structureId": pdb_id.upper(),
            "assemblyId": assembly_id,
        }
        if chain_id:
            params["chainID"] = chain_id
        response = self.session.get(self.BASE_PDB_URL, params=params)
        if response.ok:
            logging.info(f"PDB download successful for {pdb_id}")
            return response.content, "pdb"

    def parse_structure(self, raw: bytes, fmt: Optional[str], target: str):
        """
        Parse PDB/mmCIF from bytes. Robust to gzip regardless of `fmt`.
        If `fmt` is wrong, we sniff and fallback.
        """
        if not isinstance(raw, (bytes, bytearray)):
            raise TypeError("parse_structure expects raw bytes")

        # --- 1) Decompress if gzipped (regardless of fmt hint) ---
        if len(raw) >= 2 and raw[0] == 0x1F and raw[1] == 0x8B:  # gzip magic
            raw = gzip.decompress(raw)

        # normalize fmt hint
        fmt_l = (fmt or "").lower().rstrip()
        if fmt_l.endswith(".gz"):
            fmt_l = fmt_l[:-3]
        if fmt_l in ("cif",):
            fmt_l = "mmcif"

        # --- 2) Get a text handle ---
        # mmCIF/PDB are ASCII; use utf-8 with strict to avoid injecting BOM-like chars.
        text = raw.decode("utf-8", errors="strict")
        # allow leading whitespace/newlines (rare)
        first_nonempty = next((ln for ln in text.splitlines() if ln.strip()), "")
        looks_cif = first_nonempty.startswith("data_")
        looks_pdb = first_nonempty.startswith(("HEADER"))

        # --- 3) Choose parser: prefer fmt hint, but fallback on sniffed content ---
        def _parse_cif(txt: str):
            handle = io.StringIO(txt)
            parser = MMCIFParser(QUIET=True)
            return parser.get_structure(target, handle)

        def _parse_pdb(txt: str):
            handle = io.StringIO(txt)
            parser = PDBParser(QUIET=True)
            return parser.get_structure(target, handle)

        # Try in order: hinted fmt → sniffed fallback
        if fmt_l == "mmcif":
            try:
                return _parse_cif(text)
            except ValueError as e:
                # common case: actually PDB text or truncated cif; try PDB
                if looks_pdb:
                    return _parse_pdb(text)
                raise
        elif fmt_l == "pdb":
            try:
                return _parse_pdb(text)
            except Exception:
                if looks_cif:
                    return _parse_cif(text)
                raise
        else:
            # No hint: use sniffed format
            if looks_cif:
                return _parse_cif(text)
            # Default to PDB if we can't confidently identify mmCIF
            return _parse_pdb(text)
