# src/spestimator/database.py

import logging
import shutil
import subprocess
import requests
import gzip
from pathlib import Path

logger = logging.getLogger(__name__)

# URL for RefSeq 16S rRNA FASTA
REFSEQ_16S_FASTA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"
)

def get_bundled_db_prefix():
    """
    Returns the path to the bundled 16S BLAST database prefix.
    """
    package_dir = Path(__file__).parent
    data_dir = package_dir / "data"
    db_prefix = data_dir / "bacteria.16SrRNA"
    return db_prefix

def get_db_info(db_path):
    """
    Runs 'blastdbcmd -info' to get metadata about the database.
    """
    cmd = ["blastdbcmd", "-db", str(db_path), "-info"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        return "Could not retrieve database info (blastdbcmd failed or DB missing)."
    except FileNotFoundError:
        return "blastdbcmd not found in PATH."

def download_database(target_dir):
    """
    Downloads the RefSeq 16S rRNA FASTA, decompresses it, and builds a BLAST database.
    """
    target_path = Path(target_dir)
    target_path.mkdir(parents=True, exist_ok=True)

    logger.info(f"Target directory: {target_path}")

    gz_fasta = target_path / "bacteria.16SrRNA.fna.gz"
    fasta = target_path / "bacteria.16SrRNA.fna"
    db_prefix = target_path / "bacteria.16SrRNA"

    try:
        # --- 1. Download FASTA ---
        logger.info(f"Downloading RefSeq 16S FASTA from {REFSEQ_16S_FASTA_URL}...")

        with requests.get(REFSEQ_16S_FASTA_URL, stream=True) as r:
            r.raise_for_status()
            with open(gz_fasta, 'wb') as f:
                shutil.copyfileobj(r.raw, f)

        logger.info("Download complete. Decompressing...")

        # --- 2. Decompress .gz file ---
        with gzip.open(gz_fasta, "rb") as f_in, open(fasta, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        gz_fasta.unlink()  # cleanup

        # --- 3. Build BLAST DB ---
        logger.info("Building BLAST database using makeblastdb...")

        cmd = [
            "makeblastdb",
            "-in", str(fasta),
            "-dbtype", "nucl",
            "-parse_seqids",
            "-out", str(db_prefix)
        ]

        subprocess.run(cmd, check=True)

        logger.info(f"RefSeq 16S BLAST database successfully created at {db_prefix}")

    except Exception as e:
        logger.error(f"Failed to download or build database: {e}")
        if gz_fasta.exists():
            gz_fasta.unlink()
        if fasta.exists():
            fasta.unlink()
