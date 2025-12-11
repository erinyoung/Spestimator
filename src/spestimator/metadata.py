import pandas as pd
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def get_bundled_metadata_path():
    """
    Returns the path to the bundled metadata.csv.gz file.
    Assumes it lives in src/spestimator/data/
    """
    return Path(__file__).parent / "data" / "metadata.csv.gz"


def load_metadata(fpath=None):
    """
    Loads the metadata CSV into a Pandas DataFrame.

    Args:
        fpath (Path, optional): Path to the metadata.csv.gz file.
                                If None, defaults to the bundled package file.

    Expected Columns in CSV:
      - blast_sacc
      - taxid
      - organism
      - refseq_accession

    Returns:
        pd.DataFrame: Returns empty DataFrame if file is missing or corrupt.
    """
    # Determine which path to use
    if fpath:
        path = Path(fpath)
    else:
        path = get_bundled_metadata_path()

    if not path.exists():
        logger.warning(
            f"Metadata file not found at {path}. Results will use raw BLAST names."
        )
        return pd.DataFrame()

    try:
        # Load compressed CSV
        df = pd.read_csv(path, compression="gzip", dtype=str)

        # Strip whitespace from column names just in case
        df.columns = [c.strip() for c in df.columns]

        # Verify Key Column Exists
        # Note: We check for 'blast_sacc' OR 'accession' to be safe, though our builder makes 'blast_sacc'
        if "blast_sacc" not in df.columns and "accession" not in df.columns:
            logger.error(
                f"Metadata file at {path} is invalid (missing 'blast_sacc' column)."
            )
            return pd.DataFrame()

        return df

    except Exception as e:
        logger.error(f"Failed to load metadata file: {e}")
        return pd.DataFrame()
