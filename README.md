<div align="center">
  <img src="assets/logo.png" width="200" alt="Spestimator Logo">
  <h1>Spestimator</h1>
</div>

![PyPI - Version](https://img.shields.io/pypi/v/Spestimator)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/Spestimator)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/erinyoung/Spestimator/test.yml?label=tests)
![GitHub License](https://img.shields.io/github/license/erinyoung/Spestimator)

**Species Estimator & Genome Retriever**

Species + Estimator = Spestimator

Spestimator is a lightweight Python command-line tool designed to quickly identify bacterial species from input FASTA sequences (contigs or reads) and automatically retrieve their corresponding reference genomes from NCBI RefSeq.


## Why Spestimator?

I needed a tool that could:
1. Quickly estimate a range of potential reference genomes for an input sample.
2. Use a database maintained by a trustworthy, authoritative group (NCBI RefSeq) rather than maintaining a custom, static database that gets outdated.

---

## Installation

### From Pip/From Source
Spestimator requires [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata) to be installed and available in your system path.

```bash
sudo apt-get install ncbi-blast+
```

```bash
# installation from pypi
pip install spestimator
```

```bash
# from source
git clone https://github.com/erinyoung/Spestimator.git
cd Spestimator
pip install .
```

### From Conda (includes blast)
```bash
conda install -c bioconda spestimator
```




---

## Quick Start

Identify organisms in your FASTA file(s)

```bash
spestimator -i *.fasta -o results.csv
```

Identify organisms and download the matched RefSeq genomes to a folder:

```bash
spestimator -i *.fasta -o results.csv -d genomes_dir/
```

### Updating the database

RefSeq updates quarterly. There is an attempt to keep this package on a similar update schedule, but this may not be feasible. To get a current version of the database.

```bash
# downloads https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz for blast database
# uses eutils to get taxids for blast database (NCBI api-key will speed this step up)
# downloads https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt for refseq accesions
spestimator --update-db --db-dir database --api-key NCBI_API_KEY
```
This database can then be used by

```bash
spestimator -i *.fasta -o results.csv -d genomes_dir/ --db-dir database
```

---

## Usage & Options

```
usage: spestimator [-h] [-v] [-i INPUT [INPUT ...]] [-o OUTPUT] [-d [DIR]] [--db-dir DB_DIR] [--db-name DB_NAME] [-u] [--api-key API_KEY] [-t THREADS]
                   [--max-target-seqs MAX_TARGET_SEQS] [--min-identity MIN_IDENTITY] [--min-coverage MIN_COVERAGE] [--min-hits MIN_HITS]
                   [--min-alignment-len MIN_ALIGNMENT_LEN] [--top-k-taxa TOP_K_TAXA]

Spestimator: Predict bacterial TaxIDs from 16S and download genomes.

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i, --input INPUT [INPUT ...]
                        Input FASTA files
  -o, --output OUTPUT   Output CSV file
  -d, --download-genomes [DIR]
                        Download found genomes. Defaults to 'genomes/' if flag is used without a path.
  --db-dir DB_DIR       Override path to BLAST database directory
  --db-name DB_NAME     Custom name for the database to appear in results (Default: DB filename)
  -u, --update-db       Download database and generate metadata
  --api-key API_KEY     NCBI API Key (Speeds up metadata generation)
  -t, --threads THREADS
                        BLAST threads

Filtering Options:
  --max-target-seqs MAX_TARGET_SEQS
                        BLAST: Hits to keep per read (Default: 10)
  --min-identity MIN_IDENTITY
                        Filter: Minimum Percent Identity (0-100). Default: 90.0
  --min-coverage MIN_COVERAGE
                        Filter: Minimum Query Coverage (0-100). Default: 0.0
  --min-hits MIN_HITS   Filter: Minimum reads required to report an organism
  --min-alignment-len MIN_ALIGNMENT_LEN
                        Filter: Minimum Alignment Length in bp (Default: 0/No Filter)
  --top-k-taxa TOP_K_TAXA
                        Report: Only keep the top K unique organisms per file (Default: 10)
```

## How It Works
- Runs [blastn](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata) of your input against the 16S database in the repo.
- Filters hits based on identity and coverage.
- Aggregates hits to identify the most likely species present.
- (Optional): Uses [ncbi-datasets-pyclient](https://github.com/misialq/ncbi-datasets-pyclient) to download the reference genome assembly (`GCF_xxxx`) associated with the identified species.

---

## Output Format

The results CSV contains the following columns:

- **organism:** The clean species name (e.g., *E. coli*).
- **taxid:** NCBI Taxonomy ID.
- **refseq_accession:** The Assembly ID (e.g., `GCF_000005845.2`) used for downloading.
- **blast_sacc:** The specific 16S sequence accession hit.
- **count:** Number of input sequences matching this organism.
- **avg_pident:** Average percent identity of the matches.
- **avg_bitscore:** A metric combining match quality and length.

## Example

```bash
spestimator -i tests/sample_positive.fasta
```

The results should look like this

```
input file,organism,taxid,refseq_accession,blast_sacc,count,total_bitscore,avg_bitscore,avg_pident,max_pident,avg_qcov,best_evalue
sample_positive.fasta,Streptococcus pyogenes,1314,GCF_900475035.1,NR_112088,1,2772,2772.0,100.0,100.0,33.25947263461113,0.0
sample_positive.fasta,Streptococcus canis,1329,GCF_900636575.1,NR_115729,1,2545,2545.0,98.085,98.085,32.395302459561265,0.0
sample_positive.fasta,Streptococcus canis,1329,GCF_900636575.1,NR_024633,1,2536,2536.0,97.202,97.202,33.25947263461113,0.0
sample_positive.fasta,Streptococcus urinalis,149016,GCF_900636885.1,NR_037101,1,2534,2534.0,97.135,97.135,33.25947263461113,0.0
sample_positive.fasta,Streptococcus dysgalactiae subsp. equisimilis,119602,,NR_043661,1,2534,2534.0,97.443,97.443,32.92709949036118,0.0
sample_positive.fasta,Streptococcus iniae,1346,GCF_000831485.1,NR_025148,1,2532,2532.0,96.422,96.422,34.05716818081099,0.0
sample_positive.fasta,Streptococcus ictaluri 707-05,764299,,NR_115802,1,2523,2523.0,97.627,97.627,32.68335918457789,0.0
sample_positive.fasta,Streptococcus penaeicida,1765960,GCF_965139625.1,NR_178901,1,2514,2514.0,96.217,96.217,33.96853534234434,0.0
sample_positive.fasta,Streptococcus tangpeifui,2709400,GCF_011038795.1,NR_199923,1,2494,2494.0,95.858,95.858,34.23443385774429,0.0
sample_positive.fasta,Streptococcus agalactiae ATCC 13813,888745,,NR_040821,1,2484,2484.0,96.536,96.536,33.25947263461113,0.0
```

## Disclaimer: NCBI Datasets

Spestimator relies on [ncbi-datasets-pyclient](https://github.com/misialq/ncbi-datasets-pyclient) to retrieve genomes. This library updates frequently, as does the underlying [datasets](https://github.com/ncbi/datasets) tool. If you encounter unexpected errors during genome download (e.g., `ApiException` or connection drops), please [submit an issue](https://github.com/erinyoung/Spestimator/issues) and use the accessions in the "refseq_accession" column to download with [datasets](https://github.com/ncbi/datasets) separately.

```bash
# download datasets
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets && chmod +x datasets

# get a list of accessions
cut -f 4 -d , results.csv | grep GCF > id_list.txt

# use the list of accessions to download genomes
./datasets download genome accession --inputfile id_list.txt --filename ncbi_dataset.zip

# decompress file and use for additional analysis
unzip ncbi_dataset.zip
```

## AI Usage

Spestimator began as a set of custom shell and Python scripts used for ad-hoc analysis. To make these tools more reliable and accessible, Google's Gemini was used to accelerate the transition into a Python package. Gemini assisted in modularizing the codebase, replacing parsing logic with Pandas operations, and implementing a "mocked" testing suite that allows for safe CI/CD without hitting NCBI servers. Gemini also created the logo.