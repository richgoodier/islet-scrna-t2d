# Files — Sample Metadata and Download Scripts

This directory contains scripts and metadata for identifying and downloading the β-cell samples from ENA.

## Sample Selection

The original ENA project [PRJEB15401](https://www.ebi.ac.uk/ena/browser/view/PRJEB15401) contains mixed islet cell types. We filtered for β-cells only:

- `healthy_beta.txt` — 171 healthy β-cell sample IDs
- `t2d_beta.txt` — 99 T2D β-cell sample IDs
- `donor_IDs.csv` — donor metadata mapping sample IDs to donors

After quality filtering (≥750,000 reads), the final dataset was 97 healthy + 45 T2D = 142 samples.

## Downloading Data

```bash
# Download healthy β-cell FASTQ files
bash download_healthy_beta_samples.sh

# Download T2D β-cell FASTQ files
bash download_t2d_beta_samples.sh
```

These scripts use the ENA file download API. Expect ~500 GB total. Files are downloaded to a directory of your choice (edit paths in the scripts).

## Quality Filtering

After download, count reads per file and remove low-quality samples:

```bash
python count_reads.py
```

This outputs `fastq_readcount.csv`. Samples with fewer than 750,000 reads were excluded.

## Notebooks

- `creating_download_scripts.ipynb` — how the download scripts were generated from the full ENA manifest
- `donor_samples.ipynb` — organizes donor metadata
- `read_counts.ipynb` — visualizes read count distribution across samples
