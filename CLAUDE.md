# CLAUDE.md

## Session Goal (2026-04-02)

Reorganize this project directory into a professional GitHub repository. Work is tracked in [TASKS.md](TASKS.md).

**Rules:**

- Never delete or alter existing project file contents
- Move unwanted/out-of-place files to `_archive/` (gitignored), not deleted
- All task progress tracked in `TASKS.md`
- No personal information should be included in the repo.  For example, a user name or email id.

## Further instructions

- Frame as a reproducibility study.
- Break down goals into discrete tasks and delegate each task when possible.  Document your progress.
- Whenever possible, use subagents.
- Think carefully, take your time.
- Check your work after each task.  Check to see if previous tasks have been broken by your recent work.  Return to these tasks if needed.
- When you think you have finished, walk through iterations of checking everything and making necessary corrections until no more corrections are needed.
- Keep in mind that this is mainly for my portfolio to find a job.

## Project Overview

Bioinformatics course project (BINF 6310) reproducing the study "Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes" (PRJEB15401). The pipeline processes ~3,274 scRNA-seq FASTQ files from human pancreatic islet beta cells, comparing 97 healthy vs. 45 T2D samples across 62,710 genes.

## Key Tools

FastQC, MultiQC, Nextflow, cutadapt, STAR (v2.3.0e used in project; v2.7+ recommended for new runs), featureCounts (subread), edgeR (R)

## Pipeline Architecture

The analysis follows a linear pipeline documented in [pipeline.ipynb](pipeline.ipynb):

1. **Data download** — FASTQ files from ENA via scripts in [samples/](samples/) (`download_healthy_beta_samples.sh`, `download_t2d_beta_samples.sh`)
2. **Read QC filtering** — `samples/count_reads.py` counts reads per file; samples with <750,000 reads were removed
3. **Quality control** — FastQC per sample → MultiQC aggregation → `qc_reports/qc_analysis.ipynb`
4. **Adapter trimming** — Nextflow pipeline `pipeline/fastq_processing.nf` runs cutadapt (quality threshold 20) in parallel
5. **Alignment** — Nextflow pipeline `pipeline/star_align.nf` aligns to GRCh38 (Release 109) using STAR; index built by `pipeline/build_star_index.sh` (SLURM, 16 threads, 32 GB, ~2 hr)
6. **Count matrix** — featureCounts produces `gene_matrix/featureCounts_healthy/` and `gene_matrix/featureCounts_t2d/`; merged in `gene_matrix/combine_tables.ipynb` → `counts_combined.txt`
7. **RPKM normalization** — R/edgeR calculates RPKM; outputs in `rpkm_values/rpkm_combined.txt` (62,711 genes × 142 samples)
8. **Differential expression** — `analysis/differential_expression.ipynb`: filters genes present in <5 cells, computes log2 fold change and t-test p-values, generates volcano plot
9. **ID mapping** — `analysis/ensembl_to_ids.ipynb` uses MyGeneInfo API to convert Ensembl IDs to gene symbols

## Key Files

| File | Purpose |
| ---- | ------- |
| `pipeline.ipynb` | End-to-end pipeline walkthrough (primary documentation) |
| `pipeline/fastq_processing.nf` | Nextflow: parallel cutadapt trimming |
| `pipeline/star_align.nf` | Nextflow: STAR alignment (processes 149 FASTQ in ~3.5 hr) |
| `pipeline/build_star_index.sh` | SLURM job for STAR genome index |
| `samples/create_download_scripts.py` | Filters ENA download scripts by sample type |
| `samples/count_reads.py` | Counts reads per FASTQ file, outputs CSV |
| `gene_matrix/combine_tables.ipynb` | Merges healthy + T2D count matrices |
| `analysis/differential_expression.ipynb` | Main differential expression analysis |
| `rpkm_values/calculate_rpkm.R` | R script for RPKM calculation |

## Data Notes

- Large data files (FASTQ, BAM) are stored locally and are not version-controlled
- Reference genome: `Homo_sapiens.GRCh38.dna.primary_assembly.fa` + `Homo_sapiens.GRCh38.109.gtf`
- Final dataset after filtering: 142 samples (97 healthy, 45 T2D), 16,361 expressed genes
- Genes of biological interest: INS, GCG, SST, PPY (islet cell markers); TCF7L2, CDKAL1, KCNJ11 (known T2D risk loci)
