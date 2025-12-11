import pytest
from unittest.mock import MagicMock, patch
from spestimator.genome import download_genomes_bulk
import zipfile
import io


def test_download_genomes_mocked(tmp_path):
    """
    Test download logic without hitting NCBI API (using patch/mock).
    """
    # 1. Create a Fake ZIP file in memory containing a fake FASTA
    fake_zip_buffer = io.BytesIO()
    with zipfile.ZipFile(fake_zip_buffer, "w") as z:
        # Simulate NCBI folder structure
        z.writestr("ncbi_dataset/data/GCF_001/GCF_001_Assembly.fna", ">seq1\nATGC")
        z.writestr("ncbi_dataset/data/GCF_002/GCF_002_Assembly.fna", ">seq2\nTGCA")
        z.writestr("ncbi_dataset/data_report.jsonl", "{}")

    fake_zip_bytes = fake_zip_buffer.getvalue()

    # 2. Mock the ApiClient and GenomeApi
    with patch("spestimator.genome.ncbi.datasets.openapi.ApiClient") as MockClient:
        with patch("spestimator.genome.GenomeApi") as MockApi:
            mock_instance = MockApi.return_value

            # FIX: Create a mock response object that has a .data attribute
            # because your code calls: api_instance.download_assembly_package(...).data
            mock_response = MagicMock()
            mock_response.data = fake_zip_bytes

            mock_instance.download_assembly_package.return_value = mock_response

            # 3. Run the function
            output_dir = tmp_path / "genomes"
            download_genomes_bulk(["GCF_001", "GCF_002"], output_dir)

            # 4. Verify results
            expected_file_1 = output_dir / "GCF_001.fasta"
            expected_file_2 = output_dir / "GCF_002.fasta"

            assert expected_file_1.exists()
            assert expected_file_2.exists()
            assert expected_file_1.read_text() == ">seq1\nATGC"

            # Verify API was called correctly
            mock_instance.download_assembly_package.assert_called_once()
