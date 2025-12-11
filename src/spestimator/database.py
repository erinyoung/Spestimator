import logging
import shutil
import subprocess
import requests
import gzip
import pandas as pd
import time
import tempfile
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


def fetch_taxids_for_accessions(accessions, api_key=None):
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
                    results.append(
                        {
                            "blast_sacc": item.get("accessionversion", ""),
                            "taxid": str(item.get("taxid", "")),
                            "organism": item.get("title", "").split(",")[0],
                        }
                    )

            pbar.update(len(batch))
            time.sleep(0.35 if api_key else 1.0)

        except Exception as e:
            logger.error(f"Error fetching batch {i}: {e}")
            pbar.update(len(batch))
            time.sleep(5)

    pbar.close()
    return pd.DataFrame(results)


def fetch_refseq_assembly_summary():
    """
    Downloads and parses the RefSeq Assembly Summary.
    Captures: GCF Accession, TaxID, and Organism Name.
    """
    tmp_dir = Path(tempfile.gettempdir())
    tmp_path = tmp_dir / "spestimator_assembly_summary.txt"

    logger.info("Retrieving Genome Assembly Summary (~100MB+)...")

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

        # Col 0: assembly_accession, Col 4: refseq_category, Col 6: species_taxid, Col 7: organism_name, Col 11: assembly_level
        df = pd.read_csv(
            tmp_path,
            sep="\t",
            dtype=str,
            comment="#",
            header=None,
            usecols=[0, 4, 6, 7, 11],
            names=["refseq_assembly", "category", "taxid", "organism_clean", "level"],
        )

        cat_map = {"reference genome": 1, "representative genome": 2}
        df["prio_cat"] = df["category"].map(cat_map).fillna(3)
        lvl_map = {"Complete Genome": 1, "Chromosome": 2, "Scaffold": 3, "Contig": 4}
        df["prio_lvl"] = df["level"].map(lvl_map).fillna(5)

        df_sorted = df.sort_values(by=["taxid", "prio_cat", "prio_lvl"])
        df_unique = df_sorted.drop_duplicates(subset=["taxid"], keep="first")

        logger.info(f"Loaded {len(df_unique)} unique species genomes.")
        return df_unique[["taxid", "refseq_assembly", "organism_clean"]]

    except Exception as e:
        logger.error(f"Failed to process assembly summary: {e}")
        return pd.DataFrame()
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


def create_metadata_table(db_path, output_path, api_key=None):
    """Main function to generate the metadata table."""

    # 1. Get Accessions
    accessions = get_blast_accessions(db_path)
    if not accessions:
        return

    # 2. Get TaxIDs
    df_tax = fetch_taxids_for_accessions(accessions, api_key)
    if df_tax.empty:
        logger.error("Could not fetch TaxIDs. Aborting.")
        return

    # --- CLEANUP STEP ---
    # Apply cleaning to the 16S 'organism' name immediately.
    # This acts as our fallback if RefSeq doesn't have a better name.
    logger.info("Cleaning organism names...")
    df_tax["organism"] = df_tax["organism"].apply(clean_16s_title)

    # 3. Get Assembly Map
    df_assemblies = fetch_refseq_assembly_summary()

    # 4. Merge
    logger.info("Merging 16S data with Reference Genomes...")

    df_tax["taxid"] = df_tax["taxid"].astype(str)

    if not df_assemblies.empty:
        df_assemblies["taxid"] = df_assemblies["taxid"].astype(str)
        df_final = pd.merge(df_tax, df_assemblies, on="taxid", how="left")

        # PREFERENCE LOGIC:
        # 1. Use 'organism_clean' (from RefSeq Assembly) if available.
        # 2. If NA, fallback to 'organism' (which we just cleaned above).
        df_final["organism"] = df_final["organism_clean"].fillna(df_final["organism"])

        df_final.drop(columns=["organism_clean"], inplace=True)
    else:
        df_final = df_tax
        df_final["refseq_assembly"] = pd.NA

    df_final["refseq_assembly"] = df_final["refseq_assembly"].fillna("NA")

    # 5. Save
    output_path = Path(output_path)
    df_final.to_csv(output_path, index=False, compression="gzip")
    logger.info(f"Metadata saved to {output_path}")
