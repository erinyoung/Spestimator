import pytest
import subprocess
from unittest.mock import patch, MagicMock
from pathlib import Path
from spestimator.database import get_db_info, download_database, clean_16s_title

# --- Test clean_16s_title ---


def test_clean_16s_title():
    """Test the organism name cleaning logic."""
    assert (
        clean_16s_title("Alitibacter langaaensis strain ATCC 43328")
        == "Alitibacter langaaensis"
    )
    assert (
        clean_16s_title("Roseovarius maritimus 16S ribosomal RNA")
        == "Roseovarius maritimus"
    )
    assert clean_16s_title("E. coli strain K12 16S ribosomal RNA") == "E. coli"
    assert clean_16s_title("Unknown bacteria partial sequence") == "Unknown bacteria"
    assert clean_16s_title("Bacillus subtilis") == "Bacillus subtilis"


# --- Test get_db_info ---


@patch("spestimator.database.subprocess.run")
def test_get_db_info_success(mock_run, tmp_path):
    mock_run.return_value.stdout = (
        "Database: 16S rRNA\n10,000 sequences\nDate: Jan 2024"
    )
    mock_run.return_value.returncode = 0

    db_path = tmp_path / "bacteria.16SrRNA"
    info = get_db_info(db_path)

    assert "10,000 sequences" in info
    mock_run.assert_called_once()


@patch("spestimator.database.subprocess.run")
def test_get_db_info_missing_blast(mock_run, tmp_path):
    mock_run.side_effect = FileNotFoundError
    db_path = tmp_path / "bacteria.16SrRNA"
    info = get_db_info(db_path)
    assert "blastdbcmd not found" in info


# --- Test download_database ---


@patch("spestimator.database.download_file_with_progress")
@patch("spestimator.database.gzip.open")
@patch("spestimator.database.subprocess.run")
def test_download_database_success(mock_run, mock_gzip, mock_download, tmp_path):
    """
    Test the full flow: Download -> Decompress -> Build DB.
    """
    target_dir = tmp_path / "spestimator_db"

    # 1. Setup Mock for Download Side Effect
    # This creates the dummy .gz file so that .unlink() doesn't crash later.
    def create_dummy_file(url, dest_path, desc=None):
        with open(dest_path, "wb") as f:
            f.write(b"fake gzip content")
        return True

    mock_download.side_effect = create_dummy_file

    # 2. Mock Gzip
    # Configure the file handle to return empty bytes to simulate EOF for shutil.copyfileobj
    mock_file_handle = MagicMock()
    mock_file_handle.read.return_value = b""
    mock_gzip.return_value.__enter__.return_value = mock_file_handle

    # 3. Mock Subprocess
    mock_run.return_value.returncode = 0

    # Run
    download_database(target_dir, force=True)

    # Assertions
    mock_download.assert_called_once()
    assert "bacteria.16SrRNA.fna.gz" in str(mock_download.call_args[0][1])

    mock_gzip.assert_called_once()
    mock_run.assert_called_once()


@patch("spestimator.database.download_file_with_progress")
@patch("spestimator.database.gzip.open")
@patch("spestimator.database.subprocess.run")
def test_download_database_skips_if_exists(
    mock_run, mock_gzip, mock_download, tmp_path
):
    """
    Test that we skip download if file exists and force=False.
    """
    target_dir = tmp_path / "spestimator_db"
    target_dir.mkdir(parents=True)

    # Create the FASTA file so it 'exists'
    fasta = target_dir / "bacteria.16SrRNA.fna"
    fasta.touch()

    download_database(target_dir, force=False)

    # Assertions
    mock_download.assert_not_called()
    mock_gzip.assert_not_called()
    mock_run.assert_called_once()
