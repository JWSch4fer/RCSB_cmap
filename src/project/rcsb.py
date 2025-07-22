### src/project/rcsb.py
from pathlib import Path
import io
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
        if response.ok and response.content.strip().startswith(b"ATOM"):
            logging.info(f"PDB download successful for {pdb_id}")
            return response.content, "pdb"

    def parse_structure(self, data: bytes, fmt: str, target: str):
        """
        Parse raw PDB/mmCIF bytes into a Biopython Structure.

        Args:
            data: raw bytes or file-like object.
            fmt: 'pdb' or 'cif'.
            target: structure identifier.

        Returns:
            Biopython Structure object.
        """
        with tempfile.NamedTemporaryFile() as temp_file:
            if fmt == "pdb":
                parser = PDBParser(QUIET=True)
                temp_file.write(data)
                structure = parser.get_structure(target, temp_file.name)
                return structure
            if fmt == "cif":
                parser = MMCIFParser(QUIET=True)
                # create an in-memory stream to read the gzipped data
                with gzip.open(io.BytesIO(data), "rb") as gz:
                    gz_bytes = gz.read()  # .decode('utf-8')
                temp_file.write(gz_bytes)
                structure = parser.get_structure(target, temp_file.name)
                return structure
