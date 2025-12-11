import argparse
import subprocess
import pandas as pd
import logging
import time
import requests
from pathlib import Path

# --- Logging ---
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(message)s",
                    datefmt="%H:%M:%S")

# --- NCBI API constants ---
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EUTILS_DB = "assembly"

# --- Functions ---


def get_blast_accessions(db_path: Path):
    """
    Runs 'blastdbcmd' to extract all accessions from the 16S database
    using '%a' format (archive accession).
    """
    logging.info("Running blastdbcmd to extract accessions from BLAST DB...")

    cmd = ["blastdbcmd", "-db", str(db_path), "-entry", "all", "-outfmt", "%a"]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        accessions = [acc.strip() for acc in result.stdout.splitlines() if acc.strip()]
        logging.info(f"Extracted {len(accessions)} accessions.")
        return accessions

    except subprocess.CalledProcessError as e:
        logging.error(f"Error running blastdbcmd: {e.stderr}")
        return []
    except FileNotFoundError:
        logging.error("BLAST+ command 'blastdbcmd' not found. Install BLAST+ and ensure it's in PATH.")
        return []


def fetch_assembly_metadata(accessions, api_key=None):
    """
    Fetch taxonomic and assembly metadata for 16S accessions from NCBI.
    """
    results = []
    batch_size = 200  # E-Utils limit

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]
        acc_list = ','.join(batch)
        logging.info(f"Fetching metadata for batch {int(i / batch_size) + 1} ({i}-{i + len(batch) - 1})")

        try:
            # Step 1: nucleotide esearch
            esearch_url = (
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                f"?db=nucleotide&term={acc_list}&usehistory=y&retmax=200&retmode=json"
            )
            if api_key:
                esearch_url += f"&api_key={api_key}"

            r = requests.get(esearch_url, timeout=30)
            r.raise_for_status()
            data = r.json()
            webenv = data['esearchresult']['webenv']
            query_key = data['esearchresult']['querykey']

            # Step 2: nucleotide esummary
            esummary_url = (
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                f"?db=nucleotide&query_key={query_key}&WebEnv={webenv}&retmode=json"
            )
            if api_key:
                esummary_url += f"&api_key={api_key}"

            r = requests.get(esummary_url, timeout=30)
            r.raise_for_status()
            summary_data = r.json()
            uids_to_summaries = summary_data['result']['uids']

            for uid in uids_to_summaries:
                summary = summary_data['result'][uid]

                accession = summary.get('accessionversion', '').split('.')[0]

                # Clean organism name
                title = summary.get("title", "")
                organism = title.replace("16S ribosomal RNA", "") \
                                .replace("partial sequence", "") \
                                .replace("complete sequence", "").strip().strip(",")

                taxid = str(summary.get("taxid", ""))

                # Step 3: find RefSeq assembly (GCF_)
                gcf_id = ""
                if taxid:
                    assembly_search_url = (
                        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                        f"?db=assembly&term=txid{taxid}[Organism] AND latest[filter] AND refseq[filter]&retmax=1&retmode=json"
                    )
                    if api_key:
                        assembly_search_url += f"&api_key={api_key}"

                    r_asm = requests.get(assembly_search_url, timeout=30)
                    r_asm.raise_for_status()
                    asm_data = r_asm.json()
                    assembly_uids = asm_data['esearchresult']['idlist']

                    if assembly_uids:
                        assembly_esummary_url = (
                            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                            f"?db=assembly&id={assembly_uids[0]}&retmode=json"
                        )
                        if api_key:
                            assembly_esummary_url += f"&api_key={api_key}"

                        r_asm_sum = requests.get(assembly_esummary_url, timeout=30)
                        r_asm_sum.raise_for_status()
                        asm_sum_data = r_asm_sum.json()

                        if assembly_uids[0] in asm_sum_data['result']:
                            asm_summary = asm_sum_data['result'][assembly_uids[0]]['model']
                            gcf_id = asm_summary.get('refseq', '')

                results.append({
                    "blast_sacc": accession,
                    "taxid": taxid,
                    "organism": organism,
                    "refseq_accession": gcf_id
                })

            time.sleep(0.35)

        except requests.exceptions.RequestException as e:
            logging.error(f"API request failed for batch starting at {i}: {e}. Skipping batch.")
            time.sleep(5)
        except Exception as e:
            logging.error(f"General error processing batch starting at {i}: {e}. Skipping batch.")
            time.sleep(5)

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(
        description="Build metadata CSV for Spestimator from NCBI APIs."
    )
    parser.add_argument("--db-path", required=True, type=Path,
                        help="Path to 16S BLAST DB prefix (e.g., 'data/bacteria.16SrRNA')")
    parser.add_argument("--output", required=True, type=Path,
                        help="Output path for compressed CSV (e.g., 'data/metadata.csv.gz')")
    parser.add_argument("--api-key", type=str, help="NCBI API key")

    args = parser.parse_args()

    # --- BLAST DB prefix check ---
    db_prefix = args.db_path
    required_exts = [".nhr", ".nin", ".nsq"]
    found_files = [db_prefix.parent / (db_prefix.name + ext) for ext in required_exts if (db_prefix.parent / (db_prefix.name + ext)).exists()]

    if not found_files:
        logging.error(
            f"BLAST database not found using prefix '{db_prefix.resolve()}'. Checked files: "
            f"{', '.join([db_prefix.name + ext for ext in required_exts])}. Run 'spestimator --update-db' first."
        )
        return
    else:
        logging.info(f"Found BLAST DB files: {', '.join([f.name for f in found_files])}")



    # --- Extract accessions ---
    all_accessions = get_blast_accessions(db_prefix)
    if not all_accessions:
        logging.error("No accessions found in the BLAST DB. Aborting.")
        return

    # --- Fetch metadata from NCBI ---
    df_metadata = fetch_assembly_metadata(all_accessions, args.api_key)

    if df_metadata.empty:
        logging.error("Failed to retrieve any metadata. Output will be empty.")
    else:
        df_metadata.to_csv(args.output, index=False, compression='gzip')
        logging.info(f"Saved {len(df_metadata)} metadata records to {args.output}")


if __name__ == "__main__":
    main()
