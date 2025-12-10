# src/spestimator/cli.py

import argparse
import sys
import logging
import pandas as pd
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError

# --- Internal Imports ---
from spestimator.estimation import run_estimation
from spestimator.database import download_database, get_bundled_db_prefix, get_db_info
from spestimator.metadata import load_metadata
from spestimator.genome import download_genomes_bulk

def main():
    try:
        pkg_version = version("spestimator")
    except PackageNotFoundError:
        pkg_version = "unknown"

    parser = argparse.ArgumentParser(description="Spestimator: Predict bacterial TaxIDs from 16S and download genomes.")
    
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {pkg_version}")
    
    # --- Input/Output ---
    parser.add_argument("-i", "--input", nargs="+", help="Input FASTA files")
    parser.add_argument("--output", default="results.csv", help="Output CSV file")
    # FIX: Use nargs='?' and metavar="DIR" for flexible path handling
    parser.add_argument("--download-genomes", type=Path, nargs='?', const=Path("genomes"), default=None, metavar="DIR",
                        help="Download found genomes. Defaults to 'genomes/' if flag is used without a path.")

    # --- Database Args ---
    parser.add_argument("--db-dir", type=Path, help="Override path to BLAST database directory")
    
    # Optional: Allow user to name the run/db in the output
    parser.add_argument("--db-name", type=str, 
                        help="Custom name for the database to appear in results (Default: DB filename)")
    
    parser.add_argument("--update-db", action="store_true", help="Download database")
    parser.add_argument("-t", "--threads", type=int, default=4, help="BLAST threads")

    # --- Filtering Options ---
    filter_group = parser.add_argument_group("Filtering Options")
    filter_group.add_argument("--max-target-seqs", type=int, default=10, 
                              help="BLAST: Hits to keep per read (Default: 10)")
    filter_group.add_argument("--min-identity", type=float, default=90.0, 
                              help="Filter: Minimum Percent Identity (0-100). Default: 90.0")
    filter_group.add_argument("--min-coverage", type=float, default=0.0, 
                              help="Filter: Minimum Query Coverage (0-100). Default: 0.0")
    filter_group.add_argument("--min-hits", type=int, default=1, 
                              help="Filter: Minimum reads required to report an organism")
    
    filter_group.add_argument("--min-alignment-len", type=int, default=0,
                              help="Filter: Minimum Alignment Length in bp (Default: 0/No Filter)")
    filter_group.add_argument("--top-k-taxa", type=int, default=10,
                              help="Report: Only keep the top K unique organisms per file (Default: 10)")
    
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s", datefmt="%H:%M:%S")

    # --- Mode 1: Update Database ---
    if args.update_db:
        target_dir = args.db_dir if args.db_dir else Path("spestimator_db")
        download_database(target_dir)
        sys.exit(0)

    if not args.input:
        parser.print_help()
        sys.exit(0)

    # --- Mode 2: Locate Database ---
    if args.db_dir:
        db_prefix = args.db_dir / "16S_ribosomal_RNA"
    else:
        db_prefix = get_bundled_db_prefix()

    if not Path(str(db_prefix) + ".nsq").exists():
        logging.error(f"Database not found at {db_prefix}")
        logging.error("Please run with --update-db or check your installation.")
        sys.exit(1)
        
    # Set DB Name for Report
    report_db_name = args.db_name if args.db_name else db_prefix.parent.name
    logging.info(f"Using database: {db_prefix} (ID: {report_db_name})")
    
    # --- Check Database Info ---
    db_info = get_db_info(db_prefix)
    logging.info(f"Database Metadata:\n{db_info}")

    # --- Mode 3: Estimation (Phase 1: BLAST & Filter) ---
    blast_args = {'max_target_seqs': args.max_target_seqs}
    filter_args = {
        'min_hits': args.min_hits, 
        'min_identity': args.min_identity,
        'min_coverage': args.min_coverage,
        'min_alignment_len': args.min_alignment_len,
        'top_k_taxa': args.top_k_taxa
    }

    all_raw_results = []
    logging.info(f"Processing {len(args.input)} input files...")
    
    for fasta_file in args.input:
        fpath = Path(fasta_file)
        if not fpath.exists():
            logging.warning(f"Skipping {fpath.name} (File not found)")
            continue
        
        # Run estimation (Returns list of dicts)
        file_results = run_estimation(fasta_file, db_prefix, args.threads, blast_args, filter_args)
        
        if file_results:
            top_hit = file_results[0]
            organism = top_hit.get('organism', 'Unknown')
            count = top_hit.get('count', 0)
            logging.info(f"File: {fpath.name} -> Top Hit: {organism} ({count} reads)")
            all_raw_results.extend(file_results)
        else:
            # FIX: Handle Empty Results (Zero Hits) - Use 'sacc' key to match estimation.py output
            logging.warning(f"File: {fpath.name} -> No matches found.")
            all_raw_results.append({
                "input file": fpath.name,
                "organism": "No Match",
                "count": 0,
                "sacc": "",           
                "total_bitscore": 0,
                "avg_bitscore": 0,
                "avg_pident": 0,
                "max_pident": 0,
                "avg_qcov": 0,
                "best_evalue": ""
            })

    if not all_raw_results:
        logging.error("No results generated.")
        sys.exit(1)

    # Convert results to DataFrame
    df = pd.DataFrame(all_raw_results)

    # --- Mode 3: Estimation (Phase 2: Local Metadata Merge) ---
    logging.info("Loading local metadata...")
    meta_df = load_metadata()

    if not meta_df.empty:
        # FIX: Robust Column Standardization
        for col in ['sacc', 'accession']:
            if col in df.columns:
                if 'blast_sacc' in df.columns:
                    df['blast_sacc'] = df['blast_sacc'].replace("", pd.NA).fillna(df[col])
                    df.drop(columns=[col], inplace=True)
                else:
                    df.rename(columns={col: 'blast_sacc'}, inplace=True)

        if 'blast_sacc' in df.columns:
            # Clean Version Numbers for Robust Merging (e.g. NR_123.1 -> NR_123)
            df['merge_key'] = df['blast_sacc'].astype(str).str.split('.').str[0]
            
            meta_clean = meta_df.reset_index()
            col_to_clean = 'blast_sacc' if 'blast_sacc' in meta_clean.columns else meta_clean.columns[0]
            meta_clean['merge_key'] = meta_clean[col_to_clean].astype(str).str.split('.').str[0]
            
            logging.info(f"Merging results with {len(meta_clean)} metadata records...")
            
            df = df.merge(
                meta_clean[['merge_key', 'taxid', 'organism', 'refseq_accession']], 
                on='merge_key', 
                how='left', 
                suffixes=('_blast', '')
            )
            
            if 'organism' in df.columns and 'organism_blast' in df.columns:
                df['organism'] = df['organism'].fillna(df['organism_blast'])
                df.drop(columns=['organism_blast'], inplace=True)
            elif 'organism_blast' in df.columns:
                 df.rename(columns={'organism_blast': 'organism'}, inplace=True)
                
            logging.info("Metadata merged successfully.")
        else:
            logging.warning("BLAST results missing 'blast_sacc' column. Skipping metadata merge.")
    else:
        logging.warning("Skipping metadata merge (File empty or not found).")

    # --- Mode 3: Estimation (Phase 3: Formatting & Save) ---
    
    # Define Output Column Order (Database column removed)
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
        "best_evalue"
    ]
    
    final_cols = [c for c in cols if c in df.columns]
    df = df[final_cols]
    
    # Fill NAs for cleanliness
    df = df.fillna("")
    
    df.to_csv(args.output, index=False)
    logging.info(f"Results saved to {args.output}")

    # --- Mode 4: Genome Download (Bulk) ---
    if args.download_genomes:
        genome_dir = args.download_genomes
        logging.info("--- Starting Genome Downloads ---")
        
        if 'refseq_accession' in df.columns:
            valid_gcfs = df[df['refseq_accession'].astype(str).str.startswith("GCF")]
            unique_gcfs = valid_gcfs['refseq_accession'].unique().tolist()
            
            if unique_gcfs:
                download_genomes_bulk(unique_gcfs, genome_dir)
            else:
                logging.warning("No valid RefSeq GCFs found in results to download.")
        else:
             logging.warning("RefSeq accession column missing (Metadata merge failed).")

if __name__ == "__main__":
    main()
