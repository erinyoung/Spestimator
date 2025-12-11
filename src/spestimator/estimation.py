import subprocess
import pandas as pd
import logging
import sys
from pathlib import Path
from io import StringIO

logger = logging.getLogger(__name__)


def run_blast(query_path, db_path, threads, max_target_seqs=10):
    """
    Runs BLASTN locally.
    Output Format 6 columns:
    1. qseqid (Query ID)
    2. sacc (Subject Accession - used for metadata merge)
    3. stitle (Subject Title - fallback organism name)
    4. pident (Percent Identity)
    5. length (Alignment Length)
    6. qlen (Query Length - used for coverage)
    7. evalue
    8. bitscore
    """
    # We ask for stitle to have a fallback name if metadata lookup fails
    outfmt = "6 qseqid sacc stitle pident length qlen evalue bitscore"

    cmd = [
        "blastn",
        "-query",
        str(query_path),
        "-db",
        str(db_path),
        "-outfmt",
        outfmt,
        "-num_threads",
        str(threads),
        "-max_target_seqs",
        str(max_target_seqs),
    ]

    try:
        # Run BLAST
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        if not result.stdout:
            return pd.DataFrame()

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
        df = pd.read_csv(StringIO(result.stdout), sep="\t", names=cols)
        return df

    except subprocess.CalledProcessError as e:
        logger.error(f"BLAST failed: {e.stderr}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error parsing BLAST output: {e}")
        return pd.DataFrame()


def clean_organism_name(text):
    """
    Cleans up the BLAST 'stitle' to get a concise organism name.
    Matches the logic used in database.py for consistency.
    """
    if pd.isna(text):
        return "Unknown"

    # Sometimes stitle starts with the Accession (e.g., "NR_123.1 Organism...").
    # If the first word looks like an accession (contains underscore), drop it.
    first_word = text.split()[0]
    if "_" in first_word and (
        first_word.startswith("NR_") or first_word.startswith("NZ_")
    ):
        parts = text.split(maxsplit=1)
        if len(parts) > 1:
            text = parts[1]

    # Stop at specific delimiters that denote strain info
    if " strain " in text:
        text = text.split(" strain ")[0]

    # Remove common suffixes
    # Note: Order matters. Remove " 16S..." first.
    replacements = [
        " 16S ribosomal RNA",
        " 16S rRNA",
        " partial sequence",
        " complete sequence",
        " gene for 16S rRNA",
    ]

    for r in replacements:
        text = text.replace(r, "")

    # Remove trailing commas or spaces
    return text.strip().strip(",")


def process_results(df, fasta_file, filters):
    """
    Filters raw BLAST hits and aggregates them by Organism/Accession.
    """
    if df.empty:
        return []

    # 1. Calculate Coverage (length / qlen * 100)
    # Ensure no division by zero (though qlen should be > 0 from BLAST)
    df["qcov"] = (df["length"] / df["qlen"]) * 100

    # 2. Apply Filters (Row-wise)
    if filters.get("min_identity"):
        df = df[df["pident"] >= filters["min_identity"]]

    if filters.get("min_coverage"):
        df = df[df["qcov"] >= filters["min_coverage"]]

    if filters.get("min_alignment_len"):
        df = df[df["length"] >= filters["min_alignment_len"]]

    if df.empty:
        return []

    df = df.copy()

    # 3. Clean Organism Name from 'stitle'
    # This provides a clean name for the report even if metadata lookup fails later
    df["organism_clean"] = df["stitle"].apply(clean_organism_name)

    # 4. Aggregate by Accession (sacc)
    # Group by 'sacc' to get unique organisms identified by 16S ID.
    stats = (
        df.groupby("sacc")
        .agg(
            organism=("organism_clean", "first"),
            count=("qseqid", "nunique"),  # How many distinct input reads hit this?
            total_bitscore=("bitscore", "sum"),  # Sum score = abundance * quality
            avg_bitscore=("bitscore", "mean"),
            avg_pident=("pident", "mean"),
            max_pident=("pident", "max"),
            avg_qcov=("qcov", "mean"),
            best_evalue=("evalue", "min"),
        )
        .reset_index()
    )

    # 5. Filter by Min Hits (Read Count)
    if filters.get("min_hits"):
        stats = stats[stats["count"] >= filters["min_hits"]]

    # 6. Sort by Count (Descending) AND Total Bitscore
    # This bubbles the most abundant/confident matches to the top
    stats = stats.sort_values(by=["total_bitscore", "count"], ascending=[False, False])

    # 7. Apply Top-K Taxa Filter
    # e.g., only show the top 10 organisms found in the sample
    if filters.get("top_k_taxa") and filters["top_k_taxa"] > 0:
        stats = stats.head(filters["top_k_taxa"])

    # 8. Format for Output
    stats["input file"] = Path(fasta_file).name

    # Return list of dicts for easy DataFrame conversion later
    return stats.to_dict("records")


def run_estimation(fasta_file, db_path, threads, blast_args, filter_args):
    """
    Main entry point used by CLI.
    """
    df_raw = run_blast(
        fasta_file, db_path, threads, blast_args.get("max_target_seqs", 10)
    )
    results = process_results(df_raw, fasta_file, filter_args)

    return results
