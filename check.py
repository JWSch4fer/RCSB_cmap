from rcsbapi.search import TextQuery
import requests

def list_assemblies(pdb_id: str) -> list[str]:
    """
    Return all biological assembly identifiers for a PDB entry,
    e.g. ["1HFA-1", "1HFA-2", ...].
    """
    # Build and execute a full‑text query on the PDB ID
    query = TextQuery(value=pdb_id)            # TextQuery does full-text search
    assemblies = list(query(return_type="assembly"))  
    # return_type="assembly" yields IDs like "1HFA-1", "1HFA-2"
    return assemblies

def download_assembly_cifs(pdb_id: str):
    """
    For each biological assembly of pdb_id, download and save the mmCIF file.
    """
    assemblies = list_assemblies(pdb_id)
    for asm in assemblies:
        # Split "1HFA-1" into entry and assembly number
        entry, asm_id = asm.split('-', maxsplit=1)
        # Construct the CIF download URL
        url = (
            f"https://files.rcsb.org/download/"
            f"{entry.lower()}-assembly{asm_id}.cif.gz"
        )  # Biological Assembly PDBx/mmCIF download URL
        
        # Fetch the raw mmCIF file
        resp = requests.get(url)
        resp.raise_for_status()  # ensure we got a 200 OK
        
        # Save to disk
        filename = f"{entry}-assembly{asm_id}.cif.gz"
        with open(filename, "wb") as f:
            f.write(resp.content)
        print(f"Saved {filename}")

if __name__ == "__main__":
    pdb_id = "1HFA"
    download_assembly_cifs(pdb_id)

