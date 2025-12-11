import subprocess
import sys
import pytest
from pathlib import Path


# Helper to run spestimator via subprocess
def run_cli(args):
    # Ensure the script runs with the correct package context
    cmd = [sys.executable, "-m", "spestimator.cli"] + args
    return subprocess.run(cmd, capture_output=True, text=True)


def test_cli_help():
    """Ensure --help works."""
    result = run_cli(["--help"])
    assert result.returncode == 0
    assert "Spestimator" in result.stdout


def test_cli_version():
    """Ensure --version works."""
    result = run_cli(["--version"])
    assert result.returncode == 0
    assert "spestimator" in result.stdout


def test_missing_input_file():
    """Ensure it warns gracefully if an input file doesn't exist."""
    # Since the processing loop continues if one file is missing, we check for a warning
    result = run_cli(["-i", "non_existent_ghost_file.fasta"])
    # If the database isn't installed, it might exit 1, but we look for the warning message
    assert "Skipping" in result.stderr or "Database not found" in result.stderr


def test_download_genome_flag_default(test_data_dir):
    """Test that --download-genomes without a path defaults correctly."""
    input_file = test_data_dir / "sample_input.fasta"
    # Create a dummy input file for the test
    if not input_file.exists():
        input_file.write_text(">test_read\nAGCT")

    # This tests the nargs='?'
    result = run_cli(["-i", str(input_file), "--download-genomes"])

    # If the DB is missing, it will fail, but the key is that argparse passed without error
    assert "expected one argument" not in result.stderr
