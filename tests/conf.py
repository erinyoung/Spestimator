import pytest
import pandas as pd
from pathlib import Path

@pytest.fixture
def sample_blast_df():
    """
    Creates a dummy BLAST DataFrame mimicking the output of estimation.run_blast.
    Columns: qseqid, sacc, stitle, pident, length, qlen, evalue, bitscore
    """
    data = [
        # Good Hit (NR_001, count 2)
        ["read1", "NR_001", "Bacteria A strain 1", 99.0, 1500, 1500, 0.0, 2000],
        # Short Alignment (Should be filtered by min_len=1000)
        ["read2", "NR_002", "Bacteria B", 100.0, 50, 1500, 1e-5, 100],
        # Low Identity (Should be filtered by min_ident=90.0)
        ["read3", "NR_003", "Bacteria C", 80.0, 1500, 1500, 0.0, 1800],
        # Duplicate organism for Top-K testing (count 2 for NR_001)
        ["read4", "NR_001", "Bacteria A strain 2", 98.0, 1490, 1500, 0.0, 1950],
        # Medium Hit (NR_004)
        ["read5", "NR_004", "Bacteria D", 95.0, 1500, 1500, 0.0, 1900],
    ]
    cols = ["qseqid", "sacc", "stitle", "pident", "length", "qlen", "evalue", "bitscore"]
    return pd.DataFrame(data, columns=cols)

@pytest.fixture
def test_data_dir():
    """Returns the path to the tests/data directory (You may need to recreate this directory)."""
    return Path(__file__).parent / "data"
