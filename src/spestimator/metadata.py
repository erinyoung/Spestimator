import pandas as pd
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def load_metadata(fpath=None):
    """
    Loads the metadata CSV into a Pandas DataFrame.

    Args:
        fpath (Path, optional): Path to the metadata.csv.gz file.
                                If None, defaults to the bundled package file.

    Returns:
        pd.DataFrame: Returns empty DataFrame if file is missing or corrupt.
    """
    # Determine which path to use
    if fpath:
        path = Path(fpath)

    if not path.exists():
        logger.fatal(f"Metadata file not found at {path}.")
        return pd.DataFrame()

    try:
        df = pd.read_csv(path, compression="gzip", dtype=str)
        df["sacc"] = df["blast_sacc"].str.split(".").str[0]

        return df

    except Exception as e:
        logger.fatal(f"Failed to load metadata file: {e}")
        return pd.DataFrame()
