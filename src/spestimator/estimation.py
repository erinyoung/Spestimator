# src/spestimator/estimation.py

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
    outfmt = "6 qseqid sacc stitle pident length qlen evalue bitscore"
    
    cmd = [
        "blastn",
        "-query", str(query_path),
        "-db", str(db_path),
        "-outfmt", outfmt,
        "-num_threads", str(threads),
        "-max_target_seqs", str(max_target_seqs)
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if not result.stdout:
            return pd.DataFrame()

        cols = ["qseqid", "sacc", "stitle", "pident", "length", "qlen", "evalue", "bitscore"]
        df = pd.read_csv(StringIO(result.stdout), sep="\t", names=cols)
        return df

    except subprocess.CalledProcessError as e:
        logger.error(f"BLAST failed: {e.stderr}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error parsing BLAST output: {e}")
        return pd.DataFrame()

def process_results(df, fasta_file, filters):
    """
    Filters raw BLAST hits and aggregates them by Organism/Accession.
    """
    if df.empty:
        return []

    # 1. Calculate Coverage (length / qlen * 100)
    df['qcov'] = (df['length'] / df['qlen']) * 100

    # 2. Apply Filters (Row-wise)
    if filters.get('min_identity'):
        df = df[df['pident'] >= filters['min_identity']]
    
    if filters.get('min_coverage'):
        df = df[df['qcov'] >= filters['min_coverage']]

    if filters.get('min_alignment_len'):
        df = df[df['length'] >= filters['min_alignment_len']]

    if df.empty:
        return []

    # FIX: Avoid SettingWithCopyWarning
    df = df.copy()

    # 3. Clean Organism Name from 'stitle'
    def clean_organism_name(text):
        if pd.isna(text): return "Unknown"
        parts = text.split(maxsplit=1)
        if len(parts) > 1:
            text = parts[1] # Drop "NR_xxxx.1"
        
        delimiters = [" strain", " 16S", " partial", " complete", ","]
        for d in delimiters:
            if d in text:
                text = text.split(d)[0]
        return text.strip()

    df['organism_clean'] = df['stitle'].apply(clean_organism_name)

    # 4. Aggregate by Accession (sacc)
    stats = df.groupby('sacc').agg(
        organism=('organism_clean', 'first'),
        count=('qseqid', 'nunique'),
        total_bitscore=('bitscore', 'sum'),
        avg_bitscore=('bitscore', 'mean'),
        avg_pident=('pident', 'mean'),
        max_pident=('pident', 'max'),
        avg_qcov=('qcov', 'mean'),
        best_evalue=('evalue', 'min')
    ).reset_index()

    # 5. Filter by Min Hits (Read Count)
    if filters.get('min_hits'):
        stats = stats[stats['count'] >= filters['min_hits']]

    # 6. Sort by Count (Descending) AND Total Bitscore
    stats = stats.sort_values(by=["count", "total_bitscore"], ascending=[False, False])

    # 7. Apply Top-K Taxa Filter
    if filters.get('top_k_taxa') and filters['top_k_taxa'] > 0:
        stats = stats.head(filters['top_k_taxa'])

    # 8. Format for Output
    stats['input file'] = Path(fasta_file).name
    
    return stats.to_dict('records')

def run_estimation(fasta_file, db_path, threads, blast_args, filter_args):
    """
    Main entry point used by CLI.
    """
    df_raw = run_blast(fasta_file, db_path, threads, blast_args.get('max_target_seqs', 10))
    results = process_results(df_raw, fasta_file, filter_args)
    
    return results
