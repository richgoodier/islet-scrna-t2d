# QC Reports — Quality Control

This directory contains quality control outputs from FastQC and MultiQC.

## Contents

- `multiqc_report.html` — consolidated QC report across all 149 samples (open in browser)
- `multiqc_data/` — raw data underlying the MultiQC report
- `qc_analysis.ipynb` — additional QC analysis and visualization

## Running QC

FastQC and MultiQC were run on all trimmed FASTQ files:

```bash
# Run FastQC on all trimmed files (parallelized)
fastqc trimmed/*.fastq.gz -o fastqc_output/ -t 8

# Aggregate with MultiQC
multiqc fastqc_output/ -o qc_reports/
```

## Key Findings

All 149 samples passed quality thresholds. Per-base quality scores were consistently high (Phred > 30) across the dataset. See `multiqc_report.html` for the full breakdown by sample.
