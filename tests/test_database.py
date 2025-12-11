import pytest
import subprocess
from unittest.mock import patch, MagicMock
from pathlib import Path
from spestimator.database import get_db_info, download_database

# --- Test get_db_info ---

@patch('spestimator.database.subprocess.run')
def test_get_db_info_success(mock_run, tmp_path):
    """Test parsing of blastdbcmd output when it works."""
    # Mock successful output from blastdbcmd
    mock_run.return_value.stdout = "Database: 16S rRNA\n10,000 sequences\nDate: Jan 2024"
    mock_run.return_value.returncode = 0
    
    db_path = tmp_path / "bacteria.16SrRNA"
    info = get_db_info(db_path)
    
    assert "10,000 sequences" in info
    mock_run.assert_called_once()
    # Ensure correct args were passed
    cmd_args = mock_run.call_args[0][0]
    assert cmd_args[0] == "blastdbcmd"
    assert cmd_args[2] == str(db_path)

@patch('spestimator.database.subprocess.run')
def test_get_db_info_missing_blast(mock_run, tmp_path):
    """Test handling of FileNotFoundError if blast isn't installed."""
    mock_run.side_effect = FileNotFoundError
    
    db_path = tmp_path / "bacteria.16SrRNA"
    info = get_db_info(db_path)
    
    assert "blastdbcmd not found" in info

# --- Test download_database ---

@patch('spestimator.database.requests.get')
@patch('spestimator.database.gzip.open')
@patch('spestimator.database.subprocess.run')
def test_download_database_success(mock_run, mock_gzip, mock_requests, tmp_path):
    """
    Test the full flow: Download -> Decompress -> Build DB.
    """
    # 1. Mock Requests (Streaming response context manager)
    mock_response = MagicMock()
    mock_response.raise_for_status = MagicMock()
    mock_response.raw = MagicMock()
    mock_requests.return_value.__enter__.return_value = mock_response
    
    # 2. Mock Gzip (Read/Write context manager)
    mock_gzip.return_value.__enter__.return_value = MagicMock() # f_in
    
    # 3. Mock Subprocess (makeblastdb)
    mock_run.return_value.returncode = 0
    
    target_dir = tmp_path / "spestimator_db"
    
    download_database(target_dir)
    
    # Assertions:
    
    # A. Check Download
    mock_requests.assert_called_once()
    assert "ftp.ncbi.nlm.nih.gov" in mock_requests.call_args[0][0]
    
    # B. Check Decompression (gzip open called)
    mock_gzip.assert_called_once()
    
    # C. Check DB Build (makeblastdb called)
    mock_run.assert_called_once()
    cmd_args = mock_run.call_args[0][0]
    assert cmd_args[0] == "makeblastdb"
    assert cmd_args[4] == "nucl" # dbtype
    assert str(target_dir / "bacteria.16SrRNA.fna") in cmd_args[2] # input file