import pandas as pd
import pytest
from spestimator.estimation import process_results


def test_process_results_basic(sample_blast_df):
    """Test basic aggregation without filters."""
    filters = {}
    results = process_results(sample_blast_df, "test.fasta", filters)

    # We expect 4 unique accessions: NR_001, NR_002, NR_003, NR_004
    assert len(results) == 4

    # NR_001 should be top because it has a count of 2.
    assert results[0]["sacc"] == "NR_001"
    assert results[0]["count"] == 2

    # Verify the sum of bitscores for NR_001 (2000 + 1950)
    assert results[0]["total_bitscore"] == 3950

    # Verify organism name cleaning (Bacteria A strain 1 -> Bacteria A)
    # Based on our cleaning logic, "Bacteria A strain 1" splits at " strain "
    assert results[0]["organism"] == "Bacteria A"


def test_filter_min_len(sample_blast_df):
    """Test that short alignments are removed (--min-alignment-len)."""
    filters = {"min_alignment_len": 1000}
    results = process_results(sample_blast_df, "test.fasta", filters)

    # Read2 (NR_002, length 50) should be gone.
    accessions = [r["sacc"] for r in results]
    assert "NR_002" not in accessions

    # NR_001, NR_003, NR_004 should remain
    assert len(results) == 3


def test_filter_min_identity(sample_blast_df):
    """Test that low identity hits are removed (--min-identity)."""
    filters = {"min_identity": 90.0}
    results = process_results(sample_blast_df, "test.fasta", filters)

    # Read3 (NR_003, 80%) should be gone.
    accessions = [r["sacc"] for r in results]
    assert "NR_003" not in accessions

    # NR_001, NR_002, NR_004 should remain
    assert len(results) == 3


def test_top_k_filter(sample_blast_df):
    """Test that only the top K results are returned (--top-k-taxa)."""
    # The expected order based on count/score in sample_blast_df is:
    # 1. NR_001 (Count 2, Score 3950)
    # 2. NR_003 (Count 1, Score 1800)
    # 3. NR_004 (Count 1, Score 1700)
    # 4. NR_002 (Count 1, Score 100)

    filters = {"top_k_taxa": 2}
    results = process_results(sample_blast_df, "test.fasta", filters)

    assert len(results) == 2
    assert results[0]["sacc"] == "NR_001"
    assert results[1]["sacc"] == "NR_003"
