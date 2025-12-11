import argparse
import sys
import logging
import pandas as pd
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError

# --- Internal Imports ---
from spestimator.estimation import run_estimation
from spestimator.database import (
    download_database,
    get_bundled_db_prefix,
    get_db_info,
    create_metadata_table,
)
from spestimator.metadata import load_metadata
from spestimator.genome import download_genomes_bulk


def main():
    try:
        pkg_version = version("spestimator")
    except PackageNotFoundError:
        pkg_version = "unknown"

    parser = argparse.ArgumentParser(
        description="Spestimator: Predict bacterial TaxIDs from 16S and download genomes."
    )

    parser.add_argument(
        "-v", "--version", action="version", version=f"spestimator {pkg_version}"
    )

    # --- Input/Output ---
    parser.add_argument("-i", "--input", nargs="+", help="Input FASTA files")
    parser.add_argument("-o", "--output", default="results.csv", help="Output CSV file")
    parser.add_argument(
        "-d",
        "--download-genomes",
        type=Path,
        nargs="?",
        const=Path("genomes"),
        default=None,
        metavar="DIR",
        help="Download found genomes. Defaults to 'genomes/' if flag is used without a path.",
    )

    # --- Database Args ---
    parser.add_argument(
        "--db-dir", type=Path, help="Override path to BLAST database directory"
    )
    parser.add_argument(
        "--db-name",
        type=str,
        help="Custom name for the database to appear in results (Default: DB filename)",
    )

    # --- Update & API Args ---
    parser.add_argument(
        "-u",
        "--update-db",
        action="store_true",
        help="Download database and generate metadata",
    )
    parser.add_argument(
        "--api-key", type=str, help="NCBI API Key (Speeds up metadata generation)"
    )

    parser.add_argument("-t", "--threads", type=int, default=4, help="BLAST threads")

    # --- Filtering Options ---
    filter_group = parser.add_argument_group("Filtering Options")
    filter_group.add_argument(
        "--max-target-seqs",
        type=int,
        default=10,
        help="BLAST: Hits to keep per read (Default: 10)",
    )
    filter_group.add_argument(
        "--min-identity",
        type=float,
        default=90.0,
        help="Filter: Minimum Percent Identity (0-100). Default: 90.0",
    )
    filter_group.add_argument(
        "--min-coverage",
        type=float,
        default=0.0,
        help="Filter: Minimum Query Coverage (0-100). Default: 0.0",
    )
    filter_group.add_argument(
        "--min-hits",
        type=int,
        default=1,
        help="Filter: Minimum reads required to report an organism",
    )
    filter_group.add_argument(
        "--min-alignment-len",
        type=int,
        default=0,
        help="Filter: Minimum Alignment Length in bp (Default: 0/No Filter)",
    )
    filter_group.add_argument(
        "--top-k-taxa",
        type=int,
        default=10,
        help="Report: Only keep the top K unique organisms per file (Default: 10)",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(message)s", datefmt="%H:%M:%S"
    )

    # --- Mode 1: Update Database & Metadata ---
    if args.update_db:
        target_dir = args.db_dir if args.db_dir else Path("spestimator_db")

        # 1. Download & Build BLAST DB
        # force=True ensures we re-download if the file exists but might be old/corrupt
        download_database(target_dir, force=True)

        # 2. Generate Metadata
        db_prefix = target_dir / "bacteria.16SrRNA"
        metadata_output = target_dir / "metadata.csv.gz"

        # Safe check for DB existence (avoiding pathlib suffix replacement issues)
        nsq_path = Path(str(db_prefix) + ".nsq")

        if not nsq_path.exists():
            logging.error(
                f"BLAST DB creation failed (Checked for {nsq_path}). Skipping metadata generation."
            )
            sys.exit(1)

        logging.info(
            "Generating metadata table (this uses NCBI API and may take a minute)..."
        )
        create_metadata_table(db_prefix, metadata_output, api_key=args.api_key)

        logging.info("Update complete.")
        sys.exit(0)

    if not args.input:
        parser.print_help()
        sys.exit(0)

    # --- Mode 2: Locate Database ---
    if args.db_dir:
        db_prefix = args.db_dir / "bacteria.16SrRNA"
    else:
        db_prefix = get_bundled_db_prefix()

    # Check for core BLAST files
    if not Path(str(db_prefix) + ".nsq").exists():
        logging.error(f"Database not found at {db_prefix}")
        logging.error("Please run with --update-db or check your installation.")
        sys.exit(1)

    report_db_name = args.db_name if args.db_name else db_prefix.parent.name
    logging.info(f"Using database: {db_prefix} (ID: {report_db_name})")

    # --- Check Database Info ---
    db_info = get_db_info(db_prefix)
    logging.info(f"Database Metadata:\n{db_info}")

    # --- Mode 3: Estimation (Phase 1: BLAST & Filter) ---
    blast_args = {"max_target_seqs": args.max_target_seqs}
    filter_args = {
        "min_hits": args.min_hits,
        "min_identity": args.min_identity,
        "min_coverage": args.min_coverage,
        "min_alignment_len": args.min_alignment_len,
        "top_k_taxa": args.top_k_taxa,
    }

    all_raw_results = []
    logging.info(f"Processing {len(args.input)} input files...")

    for fasta_file in args.input:
        fpath = Path(fasta_file)
        if not fpath.exists():
            logging.warning(f"Skipping {fpath.name} (File not found)")
            continue

        file_results = run_estimation(
            fasta_file, db_prefix, args.threads, blast_args, filter_args
        )

        if file_results:
            top_hit = file_results[0]
            organism = top_hit.get("organism", "Unknown")
            count = top_hit.get("count", 0)
            logging.info(f"File: {fpath.name} -> Top Hit: {organism} ({count} reads)")
            all_raw_results.extend(file_results)
        else:
            logging.warning(f"File: {fpath.name} -> No matches found.")
            all_raw_results.append(
                {
                    "input file": fpath.name,
                    "organism": "No Match",
                    "count": 0,
                    "sacc": "",
                    "total_bitscore": 0,
                    "avg_bitscore": 0,
                    "avg_pident": 0,
                    "max_pident": 0,
                    "avg_qcov": 0,
                    "best_evalue": "",
                }
            )

    if not all_raw_results:
        logging.error("No results generated.")
        sys.exit(1)

    df = pd.DataFrame(all_raw_results)

    # --- Mode 3: Estimation (Phase 2: Local Metadata Merge) ---
    # Determine metadata path based on current DB location
    metadata_path = db_prefix.parent / "metadata.csv.gz"

    logging.info(f"Loading metadata from {metadata_path}...")
    meta_df = load_metadata(metadata_path)

    if not meta_df.empty:
        # 1. Normalize BLAST Dataframe Accession Column
        for col in ["sacc", "accession"]:
            if col in df.columns:
                if "blast_sacc" not in df.columns:
                    df.rename(columns={col: "blast_sacc"}, inplace=True)
                else:
                    df["blast_sacc"] = df["blast_sacc"].fillna(df[col])
                    df.drop(columns=[col], inplace=True)

        if "blast_sacc" in df.columns:
            # Create merge keys (strip version numbers if necessary, though BLAST usually returns them)
            df["merge_key"] = df["blast_sacc"].astype(str).str.split(".").str[0]

            # Normalize Metadata Accession Column
            col_to_clean = (
                "blast_sacc" if "blast_sacc" in meta_df.columns else meta_df.columns[0]
            )
            meta_df["merge_key"] = (
                meta_df[col_to_clean].astype(str).str.split(".").str[0]
            )

            logging.info(f"Merging results with {len(meta_df)} metadata records...")

            # 2. Merge
            # We include 'organism' from metadata to get the clean RefSeq name
            df = df.merge(
                meta_df[["merge_key", "taxid", "refseq_assembly", "organism"]],
                on="merge_key",
                how="left",
                suffixes=("", "_clean"),
            )

            # 3. Prioritize Clean Names
            # If metadata gave us a clean organism name, use it. Otherwise keep BLAST result.
            if "organism_clean" in df.columns:
                df["organism"] = df["organism_clean"].fillna(df["organism"])
                df.drop(columns=["organism_clean"], inplace=True)

            # Rename refseq assembly -> refseq accession for output clarity
            if "refseq_assembly" in df.columns:
                df.rename(columns={"refseq_assembly": "refseq_accession"}, inplace=True)

            logging.info("Metadata merged successfully.")
        else:
            logging.warning(
                "BLAST results missing 'blast_sacc' column. Skipping metadata merge."
            )

    # --- Mode 3: Estimation (Phase 3: Formatting & Save) ---
    cols = [
        "input file",
        "organism",
        "taxid",
        "refseq_accession",
        "blast_sacc",
        "count",
        "total_bitscore",
        "avg_bitscore",
        "avg_pident",
        "max_pident",
        "avg_qcov",
        "best_evalue",
    ]

    final_cols = [c for c in cols if c in df.columns]
    df = df[final_cols].fillna("")

    df.to_csv(args.output, index=False)
    logging.info(f"Results saved to {args.output}")

    # --- Mode 4: Genome Download (Bulk) ---
    if args.download_genomes:
        genome_dir = args.download_genomes
        logging.info("--- Starting Genome Downloads ---")

        if "refseq_accession" in df.columns:
            valid_gcfs = df[df["refseq_accession"].astype(str).str.startswith("GCF")]
            unique_gcfs = valid_gcfs["refseq_accession"].unique().tolist()

            if unique_gcfs:
                download_genomes_bulk(unique_gcfs, genome_dir)
            else:
                logging.warning("No valid RefSeq GCFs found in results to download.")
        else:
            logging.warning("RefSeq accession column missing (Metadata merge failed).")


if __name__ == "__main__":
    main()
