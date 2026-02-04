import logging
import shutil
import subprocess
import requests
import gzip
import pandas as pd
import time
import tempfile
import tarfile
import os
from pathlib import Path
from tqdm import tqdm

logger = logging.getLogger(__name__)

# --- URLs ---
REFSEQ_16S_FASTA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"
)
REFSEQ_ASSEMBLY_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
)

TAXDUMP_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
)

# --- NCBI API Constants ---
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def get_bundled_db_prefix():
    """Returns the path to the bundled 16S BLAST database prefix."""
    package_dir = Path(__file__).parent
    data_dir = package_dir / "data"
    db_prefix = data_dir / "bacteria.16SrRNA"
    return db_prefix


def get_db_info(db_path):
    """Runs 'blastdbcmd -info' to get metadata about the database."""
    cmd = ["blastdbcmd", "-db", str(db_path), "-info"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        return "Could not retrieve database info (blastdbcmd failed or DB missing)."
    except FileNotFoundError:
        return "blastdbcmd not found in PATH."


def download_file_with_progress(url, dest_path, desc="Downloading"):
    """Helper to download a file with a tqdm progress bar and timeout."""
    try:
        with requests.get(url, stream=True, timeout=30) as r:
            r.raise_for_status()
            total_size = int(r.headers.get("content-length", 0))

            with (
                open(dest_path, "wb") as f,
                tqdm(
                    desc=desc,
                    total=total_size,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    leave=True,
                ) as bar,
            ):
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        size = f.write(chunk)
                        bar.update(size)
        return True
    except Exception as e:
        logger.error(f"Failed to download {url}: {e}")
        return False


def download_database(target_dir, force=False):
    """Downloads the RefSeq 16S rRNA FASTA, decompresses it, and builds a BLAST database."""
    target_path = Path(target_dir)
    target_path.mkdir(parents=True, exist_ok=True)

    logger.info(f"Target directory: {target_path}")

    gz_fasta = target_path / "bacteria.16SrRNA.fna.gz"
    fasta = target_path / "bacteria.16SrRNA.fna"
    db_prefix = target_path / "bacteria.16SrRNA"

    try:
        # --- 1. Download FASTA ---
        if force or not fasta.exists():
            logger.info("Starting download of 16S FASTA...")
            if fasta.exists():
                fasta.unlink()

            success = download_file_with_progress(
                REFSEQ_16S_FASTA_URL, gz_fasta, desc="Downloading 16S FASTA"
            )
            if not success:
                return

            logger.info("Decompressing FASTA...")
            with gzip.open(gz_fasta, "rb") as f_in, open(fasta, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            gz_fasta.unlink()
        else:
            logger.info("FASTA already exists. Skipping download.")

        # --- 2. Build BLAST DB ---
        logger.info("Building BLAST database...")
        cmd = [
            "makeblastdb",
            "-in",
            str(fasta),
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-out",
            str(db_prefix),
        ]
        subprocess.run(cmd, check=True)
        logger.info(f"BLAST database created at {db_prefix}")

    except Exception as e:
        logger.error(f"Failed to build database: {e}")
        if gz_fasta.exists():
            gz_fasta.unlink()


def get_blast_accessions(db_path):
    """Runs 'blastdbcmd' to extract all accessions."""
    logger.info("Extracting accessions from local BLAST database...")
    # Note: This blast database is missing many fields, even if TAXDB is present
    cmd = ["blastdbcmd", "-db", str(db_path), "-entry", "all", "-outfmt", "%a"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        accessions = [acc.strip() for acc in result.stdout.splitlines() if acc.strip()]
        logger.info(f"Extracted {len(accessions)} accessions.")
        return accessions
    except Exception as e:
        logger.error(f"Error extracting accessions: {e}")
        return []


def fetch_taxids_for_accessions(accessions, taxdump_dict, api_key=None):
    """Batches accessions and fetches their TaxIDs using NCBI eSummary."""
    base_url = f"{EUTILS_BASE}esummary.fcgi"
    results = []
    batch_size = 300

    pbar = tqdm(total=len(accessions), desc="Fetching TaxIDs (API)")

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        ids = ",".join(batch)

        params = {"db": "nucleotide", "id": ids, "retmode": "json", "version": "2.0"}
        if api_key:
            params["api_key"] = api_key

        try:
            r = requests.post(base_url, data=params, timeout=30)
            r.raise_for_status()
            data = r.json()

            if "result" in data and "uids" in data["result"]:
                uids = data["result"]["uids"]
                for uid in uids:
                    item = data["result"][uid]
                    taxid = str(item.get("taxid",""))
                    if taxid:
                        species_taxid = get_species_id(taxid, taxdump_dict)
                    else:
                        species_taxid = ""
                    results.append(
                        {
                            "blast_sacc": item.get("accessionversion", ""),
                            "blast_taxid": taxid,
                            "blast_species_taxid": species_taxid,
                            "blast_organism": item.get("title", "").split(",")[0],
                        }
                    )

            pbar.update(len(batch))
            time.sleep(0.35 if api_key else 1.0)

        except Exception as e:
            logger.error(f"Error fetching batch {i}: {e}")
            pbar.update(len(batch))
            time.sleep(5)

    pbar.close()
    
    df = pd.DataFrame(results)
    df["blast_organism"] = df["blast_organism"].apply(clean_16s_title)

    return df


def get_species_id(query_taxid, node_dict):
    query_taxid = str(query_taxid)
    
    # Check if ID exists
    if query_taxid not in node_dict:
        return f"Error: TaxID {query_taxid} not found in database."

    current_id = query_taxid

    while True:
        node = node_dict[current_id]
        rank = node['rank']
        parent = node['parent']
        
        if rank == 'species':
            return current_id
            
        if current_id == parent or current_id == '1':
            # We reached the top without hitting 'species'. 
            # This happens if the input is higher than species (e.g. a Genus)
            return None
        
        # 3. Move up the tree
        current_id = parent

def fetch_refseq_assembly_summary(tax_list):
    """
    Downloads and parses the RefSeq Assembly Summary.
    Captures: GCF Accession, TaxID, and Organism Name.
    """
    tmp_dir = Path(tempfile.gettempdir())
    tmp_path = tmp_dir / "assembly_summary.txt"

    logger.info("Retrieving Genome Assembly Summary...")

    try:
        if not tmp_path.exists():
            success = download_file_with_progress(
                REFSEQ_ASSEMBLY_SUMMARY_URL,
                tmp_path,
                desc="Downloading Assembly Summary",
            )
            if not success:
                return pd.DataFrame()
            
        else:
            logger.info("Using cached assembly summary from temp...")

        logger.info("Parsing Assembly Summary...")

        df = pd.read_csv(
            tmp_path,
            sep="\t",
            dtype=str,
            header=1,
            usecols=["#assembly_accession", "refseq_category", "taxid", "species_taxid", "organism_name"],
        )

        df.rename(columns={
            '#assembly_accession': 'assembly_accession',
            'organism_name': 'refseq_organism_name'
        }, inplace=True)
        df = df[df['refseq_category'] != 'na']

        # note: some of these will match on taxid and some will match on species_taxid
        df['any_taxid'] = df['taxid']

        df_species = df.copy()
        df_species['any_taxid'] = df_species['species_taxid']

        df = pd.concat([df, df_species], ignore_index=True)
        # for when species_taxid == strain_taxid
        df = df.drop_duplicates()

        df = df[df['any_taxid'].isin(tax_list)]

        found_taxids = set(df['any_taxid'])
        missing_taxids = list(set(tax_list) - found_taxids)

        logging.info(f"Count of missing IDs: {len(missing_taxids)}")
        logging.info(f"Missing IDs: {missing_taxids}")

        logger.info(f"Identified {len(df)} reference genomes for {len(tax_list)} taxids.")
        return df

    except Exception as e:
        logger.error(f"Failed to process assembly summary: {e}")
        return pd.DataFrame()
    finally:
        pass


def fetch_taxdump_files():
    """
    Downloads and parses the Taxdump files for taxonomy walking.
    """
    tmp_dir = Path(tempfile.gettempdir())
    tmp_path = tmp_dir / "taxdump.tar.gz"

    nodes = {}

    logger.info("Retrieving Taxdump files...")

    try:
        if not tmp_path.exists():
            success = download_file_with_progress(
                TAXDUMP_SUMMARY_URL,
                tmp_path,
                desc="Downloading TaxDump files",
            )
            if not success:
                return pd.DataFrame()
            
        else:
            logger.info("Using taxdump files from temp...")

        logger.info("Parsing Taxdump files for taxonomy walking...")

        with tarfile.open(tmp_path, "r:gz") as tar:
            # Extract nodes.dmp as a file-like object
            with tar.extractfile("nodes.dmp") as f:
                # Read line by line
                for line in f:
                    # Decode bytes to string and strip whitespace
                    line = line.decode('utf-8').strip()
                    parts = line.split('|')

                    taxid = parts[0].strip()
                    parent_taxid = parts[1].strip()
                    rank = parts[2].strip()
                    
                    # Store in a dict: tax_id -> (parent, rank)
                    nodes[taxid] = {'parent': parent_taxid, 'rank': rank}
        return nodes

    except Exception as e:
        logger.error(f"Failed to download and process taxdump files: {e}")
        return {}
    finally:
        pass

def clean_16s_title(text):
    """
    Cleans up a raw 16S sequence title.
    1. Cuts off at ' strain ' to remove specific strain info.
    2. Removes common suffixes like ' 16S ribosomal RNA'.
    """
    if not isinstance(text, str):
        return text

    # Logic 1: Stop at " strain "
    if " strain " in text:
        text = text.split(" strain ")[0]

    # Logic 2: Remove " 16S..." if it remains (or if "strain" wasn't there)
    text = text.replace(" 16S ribosomal RNA", "")
    text = text.replace(" partial sequence", "")
    text = text.replace(" complete sequence", "")
    text = text.replace(" gene for 16S rRNA", "")

    return text.strip()

import random
def create_metadata_table(db_path, output_path, api_key=None):
    """Main function to generate the metadata table."""

    # 1. Get Accessions from blast_db
    accessions = get_blast_accessions(db_path)

    # keep for testing
    # accessions = random.sample(accessions, 50)

    # 2. Get Taxdump files
    taxdump_dict = fetch_taxdump_files()

    # 3. Get TaxIDs for each Accession
    df = fetch_taxids_for_accessions(accessions, taxdump_dict, api_key)

    # 4. Get RefSeq Reference Accesions for each TaxID
    taxid_list = pd.concat([df['blast_taxid'], df['blast_species_taxid']]).unique().tolist()
    df_assemblies = fetch_refseq_assembly_summary(taxid_list)

    # 5. Combine Dataframes
    df_taxid = df.merge(df_assemblies, left_on='blast_species_taxid', right_on='any_taxid', how = 'left')
    df_species = df.merge(df_assemblies, left_on='blast_taxid', right_on='any_taxid', how = 'left')

    final_df = pd.concat([df_taxid, df_species], ignore_index=True)
    final_df = final_df.sort_values(['blast_organism', 'blast_sacc'])
    final_df = final_df.drop_duplicates()

    # 6. Save
    output_path = Path(output_path)
    final_df.to_csv(output_path, index=False, compression="gzip")
    logger.info(f"DB Metadata saved to {output_path}")

