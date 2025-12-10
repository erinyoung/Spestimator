# scripts/build_metadata.py

import argparse
import subprocess
import pandas as pd
import logging
import time
import requests
from pathlib import Path

# Set up logging for console feedback
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s", datefmt="%H:%M:%S")

# --- NCBI API Constants ---
# Use the E-utilities API for quick assembly summary lookups
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EUTILS_DB = "assembly" # The target database for finding RefSeq assemblies

def get_blast_accessions(db_path):
    """
    Runs 'blastdbcmd' to extract all accessions from the 16S database.
    This provides the list of NR_... IDs we need to query the API with.
    """
    logging.info("Running blastdbcmd to extract accessions...")
    # Get all sequence IDs (accessions) in the database
    cmd = ["blastdbcmd", "-db", str(db_path), "-entry", "all", "-outfmt", "%a"]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Split by newline and filter out any empty strings
        accessions = [acc.strip() for acc in result.stdout.split('\n') if acc.strip()]
        logging.info(f"Extracted {len(accessions)} accessions.")
        return accessions
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running blastdbcmd: {e.stderr}")
        return []
    except FileNotFoundError:
        logging.error("BLAST+ command 'blastdbcmd' not found. Is BLAST+ installed and in your PATH?")
        return []

def fetch_assembly_metadata(accessions, api_key=None):
    """
    Fetches taxonomic and assembly metadata for a batch of 16S accessions (NR_...).
    This is the slow, API-dependent step.
    """
    # Create a DataFrame for results
    results = []
    
    # Process in batches of 200 (E-Utils limit)
    batch_size = 200
    
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]
        acc_list = ','.join(batch)
        
        logging.info(f"Fetching metadata for batch {int(i/batch_size) + 1} (Accessions {i} to {i + len(batch) - 1})...")

        # --- Step 1: Query NCBI to find assembly IDs (UIDs) based on 16S accessions ---
        # We query the 'nucleotide' database first to get the TaxID and then the 'assembly' database.
        
        # Querying nucleotide to get the TaxID (primary goal)
        esearch_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=nucleotide&term={acc_list}&usehistory=y&retmax=200&retmode=json"
        )
        if api_key:
            esearch_url += f"&api_key={api_key}"

        try:
            r = requests.get(esearch_url, timeout=30)
            r.raise_for_status()
            data = r.json()
            
            # Extract the WebEnv and QueryKey for history-based fetching
            webenv = data['esearchresult']['webenv']
            query_key = data['esearchresult']['querykey']
            
            # --- Step 2: Fetch the summary (esummary) using the history context ---
            esummary_url = (
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                f"?db=nucleotide&query_key={query_key}&WebEnv={webenv}&retmode=json"
            )
            if api_key:
                esummary_url += f"&api_key={api_key}"
                
            r = requests.get(esummary_url, timeout=30)
            r.raise_for_status()
            summary_data = r.json()
            
            # The result is a dictionary mapping UIDs to summary data
            uids_to_summaries = summary_data['result']['uids']
            
            for uid in uids_to_summaries:
                summary = summary_data['result'][uid]
                
                # Extract core information
                accession = summary.get('accessionversion', '').split('.')[0] # NR_XXXXX
                taxid = str(summary.get('taxid', ''))
                organism = summary.get('title', '').split(maxsplit=1)[-1].strip()
                
                # --- Step 3: Find the corresponding RefSeq Assembly ID (GCF_) ---
                # This requires an *additional* query to the Assembly database using the TaxID
                gcf_id = ""
                if taxid:
                    assembly_search_url = (
                        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                        f"?db=assembly&term=txid{taxid}[Organism] AND latest[filter] AND refseq[filter] AND 'complete genome'[Title]&retmax=1&retmode=json"
                    )
                    if api_key:
                        assembly_search_url += f"&api_key={api_key}"
                        
                    r_asm = requests.get(assembly_search_url, timeout=30)
                    r_asm.raise_for_status()
                    asm_data = r_asm.json()
                    
                    assembly_uids = asm_data['esearchresult']['idlist']
                    
                    if assembly_uids:
                        # Fetch the actual accession for the assembly UID
                        assembly_esummary_url = (
                            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                            f"?db=assembly&id={assembly_uids[0]}&retmode=json"
                        )
                        if api_key:
                            assembly_esummary_url += f"&api_key={api_key}"
                            
                        r_asm_sum = requests.get(assembly_esummary_url, timeout=30)
                        r_asm_sum.raise_for_status()
                        asm_sum_data = r_asm_sum.json()
                        
                        # Extract the RefSeq Assembly ID (GCF_...)
                        if assembly_uids[0] in asm_sum_data['result']:
                            asm_summary = asm_sum_data['result'][assembly_uids[0]]['model']
                            gcf_id = asm_summary.get('refseq', '')
                
                results.append({
                    "blast_sacc": accession,
                    "taxid": taxid,
                    "organism": organism,
                    "refseq_accession": gcf_id
                })
                
            time.sleep(0.35) # Respect the NCBI limit (3 queries/second)

        except requests.exceptions.RequestException as e:
            logging.error(f"API request failed for batch starting at {i}: {e}. Skipping batch.")
            time.sleep(5) # Wait longer after failure
        except Exception as e:
            logging.error(f"General error processing batch starting at {i}: {e}. Skipping batch.")
            time.sleep(5)

    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description="Builds the metadata CSV file for Spestimator from NCBI APIs.")
    parser.add_argument("--db-path", required=True, type=Path, help="Path to the 16S BLAST database prefix (e.g., 'data/16S_ribosomal_RNA')")
    parser.add_argument("--output", required=True, type=Path, help="Output path for the compressed CSV (e.g., 'data/metadata.csv.gz')")
    parser.add_argument("--api-key", type=str, help="NCBI API key for faster, higher-volume queries.")
    
    args = parser.parse_args()
    
    if not args.db_path.exists():
        logging.error(f"Database prefix not found at {args.db_path}. Run 'spestimator --update-db' first.")
        return

    # 1. Get all accessions from the local BLAST database
    all_accessions = get_blast_accessions(args.db_path)
    
    if not all_accessions:
        logging.error("No accessions found in the database. Aborting metadata build.")
        return

    # 2. Query NCBI APIs for metadata
    df_metadata = fetch_assembly_metadata(all_accessions, args.api_key)
    
    if df_metadata.empty:
        logging.error("Failed to retrieve any metadata from NCBI. Output file will be empty.")
    else:
        # 3. Save to compressed CSV
        df_metadata.to_csv(args.output, index=False, compression='gzip')
        logging.info(f"Successfully saved {len(df_metadata)} metadata records to {args.output}")

if __name__ == "__main__":
    main()
