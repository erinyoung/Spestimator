import logging
import time
import zipfile
import io
import shutil
from pathlib import Path

# --- Client Library Imports ---
import ncbi.datasets.openapi
from ncbi.datasets.openapi.api import GenomeApi
from ncbi.datasets.openapi.rest import ApiException

logger = logging.getLogger(__name__)

# Configure the client globally
CONFIGURATION = ncbi.datasets.openapi.Configuration(
    host="https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
)


def download_genomes_bulk(accession_list, output_dir):
    """
    Uses the ncbi-datasets client to download a batch of genomes as a ZIP.
    Extracts the .fna (FASTA) files to the output directory.
    """
    if not accession_list:
        return

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Filter: Only download what is missing locally
    to_download = [
        acc
        for acc in set(accession_list)
        if acc and not (output_dir / f"{acc}.fasta").exists()
    ]

    if not to_download:
        logger.info("All requested genomes are already present.")
        return

    logger.info(f"Downloading {len(to_download)} genomes using NCBI Client...")

    # Use the client as a context manager
    with ncbi.datasets.openapi.ApiClient(CONFIGURATION) as api_client:
        api_instance = GenomeApi(api_client)

        try:
            # Respect rate limits slightly
            time.sleep(0.35)

            # FIX: Removed '_preload_content=False'.
            # Recent versions of the library return the bytes directly by default for this endpoint.
            response = api_instance.download_assembly_package(
                accessions=to_download, include_annotation_type=["GENOME_FASTA"]
            )

            # Handle potential response wrappers (just in case of version differences)
            if hasattr(response, "data"):
                zip_bytes = response.data
            else:
                zip_bytes = response

            # Wrap the raw bytes in a buffer
            zip_buffer = io.BytesIO(zip_bytes)

            # Extract Files
            with zipfile.ZipFile(zip_buffer) as z:
                extracted_count = 0
                for filename in z.namelist():
                    # Only extract the FASTA files
                    if filename.endswith(".fna"):
                        # Extract the GCF accession from the path
                        parts = filename.split("/")
                        gcf_id = next((p for p in parts if p.startswith("GCF_")), None)

                        if gcf_id:
                            target_path = output_dir / f"{gcf_id}.fasta"

                            with (
                                z.open(filename) as source,
                                open(target_path, "wb") as target,
                            ):
                                shutil.copyfileobj(source, target)

                            extracted_count += 1

            logger.info(
                f"Successfully downloaded and extracted {extracted_count} genomes."
            )

        except ApiException as e:
            logger.error(f"NCBI Client Error (Download): {e.status} - {e.reason}")
            logger.error("Double check that these accessions exist in RefSeq.")
        except Exception as e:
            logger.error(f"Unexpected error during download: {e}")
