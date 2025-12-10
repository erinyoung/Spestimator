# src/spestimator/database.py

import logging
import shutil
import tarfile
import subprocess
import requests
from pathlib import Path

logger = logging.getLogger(__name__)

# URL for the official NCBI 16S rRNA BLAST database
NCBI_16S_URL = "https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz"

def get_bundled_db_prefix():
    """
    Returns the path to the bundled 16S BLAST database prefix.
    Assumes file structure: src/spestimator/data/16S_ribosomal_RNA.*
    """
    # Locate the 'data' directory relative to this file
    package_dir = Path(__file__).parent
    data_dir = package_dir / "data"
    db_prefix = data_dir / "16S_ribosomal_RNA"
    return db_prefix

def get_db_info(db_path):
    """
    Runs 'blastdbcmd -info' to get metadata about the database (date, sequences, etc).
    """
    cmd = ["blastdbcmd", "-db", str(db_path), "-info"]
    try:
        # Runs command and captures output
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        # This returns if the BLAST DB is incomplete or corrupt
        return "Could not retrieve database info (blastdbcmd failed or DB missing)."
    except FileNotFoundError:
        # This returns if the 'blastdbcmd' binary is missing from PATH
        return "blastdbcmd not found in PATH."

def download_database(target_dir):
    """
    Downloads and extracts the latest 16S_ribosomal_RNA database from NCBI.
    
    This function deliberately avoids downloading the large taxdb.tar.gz file
    because metadata (TaxID, GCF) is built separately via an API query.
    """
    target_path = Path(target_dir)
    target_path.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Target directory: {target_path}")

    # --- 1. Download 16S Database ---
    tar_16s = target_path / "16S_ribosomal_RNA.tar.gz"
    
    try:
        logger.info(f"Downloading 16S database from {NCBI_16S_URL}...")
        
        # Download with streaming
        with requests.get(NCBI_16S_URL, stream=True) as r:
            r.raise_for_status()
            with open(tar_16s, 'wb') as f:
                shutil.copyfileobj(r.raw, f)
        
        logger.info("Download complete. Extracting...")
        
        # 2. Extract
        with tarfile.open(tar_16s, "r:gz") as tar:
            # This extracts the .nsq, .nhr, .nin files etc.
            tar.extractall(path=target_path)
        
        logger.info(f"16S BLAST database successfully updated in {target_path}")
        
        # 3. Cleanup
        tar_16s.unlink() 
        
    except Exception as e:
        logger.error(f"Failed to download/extract database: {e}")
        # Clean up partial file if it failed
        if tar_16s.exists():
            tar_16s.unlink()