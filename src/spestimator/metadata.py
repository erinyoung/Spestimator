# src/spestimator/metadata.py

import pandas as pd
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

def get_metadata_path():
    """
    Returns the path to the bundled metadata.csv.gz file.
    Assumes it lives in src/spestimator/data/
    """
    # __file__ is .../spestimator/metadata.py
    # We want .../spestimator/data/metadata.csv.gz
    return Path(__file__).parent / "data" / "metadata.csv.gz"

def load_metadata():
    """
    Loads the bundled metadata CSV into a Pandas DataFrame.
    
    Expected Columns in CSV:
      - blast_sacc 
      - taxid
      - organism
      - refseq_accession
      
    Returns:
        pd.DataFrame: Returns empty DataFrame if file is missing or corrupt.
    """
    path = get_metadata_path()
    
    if not path.exists():
        logger.warning(f"Metadata file not found at {path}. Results will use raw BLAST names.")
        return pd.DataFrame()
    
    try:
        # Load compressed CSV
        df = pd.read_csv(path, compression='gzip', dtype=str)
        
        # Strip whitespace from column names just in case
        df.columns = [c.strip() for c in df.columns]
        
        # Verify Key Column Exists
        if 'blast_sacc' not in df.columns:
            logger.error(f"Metadata file at {path} is invalid (missing 'blast_sacc' column).")
            return pd.DataFrame()
        
        # We merge on the 'blast_sacc' column directly in cli.py, so we return the full DataFrame.
        return df
        
    except Exception as e:
        logger.error(f"Failed to load metadata file: {e}")
        return pd.DataFrame()
