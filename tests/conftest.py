import pytest
import pandas as pd
from pathlib import Path

# --- Fixtures for Sample Files ---


@pytest.fixture
def fasta_positive(tmp_path):
    """Creates a FASTA file with a valid E. coli sequence."""
    p = tmp_path / "sample_positive.fasta"
    seq = (
        "GCTTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAAT"
        "GTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAG"
    )
    p.write_text(f">Test_Positive|E_coli\n{seq}\n")
    return p


@pytest.fixture
def fasta_nomatch(tmp_path):
    """Creates a FASTA file that yields no BLAST hits."""
    p = tmp_path / "sample_nomatch.fasta"
    seq = "TGCAT" * 50
    p.write_text(f">Test_NoMatch|Alien\n{seq}\n")
    return p


@pytest.fixture
def fasta_error(tmp_path):
    """Creates a FASTA file that triggers a BLAST error (invalid chars)."""
    p = tmp_path / "sample_error.fasta"
    p.write_text(">Test_Error|Invalid\nATGCATGCAT!!@@##$$%%^^&&**\n")
    return p


# --- Fixtures for Data Structures ---


@pytest.fixture
def sample_blast_df():
    """
    Creates a dummy BLAST DataFrame mimicking output for estimation.py tests.
    """
    data = [
        # Hit 1: Good Match (NR_001) - Count 2 (read1, read4)
        ["read1", "NR_001", "Bacteria A strain 1", 99.0, 1500, 1500, 0.0, 2000],
        ["read2", "NR_002", "Bacteria B", 100.0, 50, 1500, 1e-5, 100],
        ["read3", "NR_003", "Bacteria C", 80.0, 1500, 1500, 0.0, 1800],
        ["read4", "NR_001", "Bacteria A strain 2", 98.0, 1490, 1500, 0.0, 1950],
        ["read5", "NR_004", "Bacteria D", 95.0, 1500, 1500, 0.0, 1700],
    ]
    cols = [
        "qseqid",
        "sacc",
        "stitle",
        "pident",
        "length",
        "qlen",
        "evalue",
        "bitscore",
    ]
    return pd.DataFrame(data, columns=cols)


@pytest.fixture
def test_data_dir():
    """Returns the path to the tests/data directory."""
    # Ensure this directory actually exists or create it
    d = Path(__file__).parent / "data"
    d.mkdir(exist_ok=True)
    return d


@pytest.fixture
def mock_db_path(tmp_path):
    """Creates dummy BLAST database files so existence checks pass."""
    db_dir = tmp_path / "mock_db"
    db_dir.mkdir()
    prefix = db_dir / "bacteria.16SrRNA"

    # Touch required extensions
    for ext in [".nsq", ".nin", ".nhr"]:
        (db_dir / f"bacteria.16SrRNA{ext}").touch()

    return prefix
