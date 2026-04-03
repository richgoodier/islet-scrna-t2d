# Islet scRNA-seq in Type 2 Diabetes

A reproduction and extension of the study **"Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes"** (Segerstolpe et al., 2016), produced as a course project for BINF 6310 (Spring 2025).

We reproduced the core RNA-seq pipeline from raw FASTQ files through differential expression analysis, then extended the work by training a logistic regression classifier to identify transcriptomic predictors of T2D — with **IRF2BPL** emerging as the top predictive feature alongside 9 established T2D risk genes (GCK, CDKAL1, KCNJ11, SLC30A8, IRS1, GLIS3, JAZF1, SLC16A11, SREBF1).

## Key Results

| Analysis | Finding |
| --- | --- |
| Marker expression | INS significantly lower in T2D β-cells (p < 0.001) |
| Differential expression | GCG, INS, SST all significantly different (p < 0.05) |
| t-SNE clustering | Healthy and T2D β-cells form distinct clusters |
| GSEA | T2D mellitus and insulin resistance pathways enriched |
| Predictive model | Logistic regression (L1) identified 15 top predictive genes |

Figures are in [`analysis/figures/`](analysis/figures/) and [`presentation/`](presentation/).

## Data

Raw data: ENA accession [PRJEB15401](https://www.ebi.ac.uk/ena/browser/view/PRJEB15401)

- 142 β-cell samples: 97 healthy, 45 T2D (after quality filtering)
- 62,710 genes → 16,361 expressed in ≥5 cells
- Reference genome: Homo sapiens GRCh38, Ensembl Release 109

Large intermediate files (count matrices, RPKM tables) are not tracked in git. See the individual subdirectory READMEs for how to regenerate them.

## Repository Structure

```text
├── pipeline.ipynb          # Master pipeline documentation (start here)
├── pipeline/               # HPC pipeline scripts (Nextflow, SLURM)
├── samples/                # Sample metadata and download scripts
├── gene_matrix/            # Count matrix generation
├── rpkm_values/            # RPKM normalization (R/edgeR)
├── analysis/               # Local Python analysis and figures
├── qc_reports/             # FastQC / MultiQC quality control
├── presentation/           # Final presentation slides
├── docs/                   # Project documentation and writeup
└── _archive/               # Out-of-place files (not tracked in git)
```

## Quickstart

### Environment setup

```bash
conda env create -f environment.yml
conda activate islet-scrna-t2d
```

### HPC pipeline (requires cluster access)

See [`pipeline/README.md`](pipeline/README.md) and [`pipeline.ipynb`](pipeline.ipynb) for full instructions. Requires STAR v2.3.0e, cutadapt, Nextflow, featureCounts, and R/edgeR.

### Local analysis (runs on laptop)

With RPKM data in place, open notebooks in order:

1. [`gene_matrix/combine_tables.ipynb`](gene_matrix/combine_tables.ipynb) — merge count matrices
2. [`analysis/ensembl_to_ids.ipynb`](analysis/ensembl_to_ids.ipynb) — map gene IDs
3. [`analysis/differential_expression.ipynb`](analysis/differential_expression.ipynb) — differential expression
4. [`analysis/expression_analysis.ipynb`](analysis/expression_analysis.ipynb) — expression plots, GSEA
5. [`analysis/ml_classification.ipynb`](analysis/ml_classification.ipynb) — t-SNE, PCA, logistic regression

## References

1. Segerstolpe Å, et al. *Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes.* Cell Metabolism. 2016;24(4):593-607. [doi:10.1016/j.cmet.2016.08.020](https://doi.org/10.1016/j.cmet.2016.08.020)

## Acknowledgments

This project was developed collaboratively by Group 3 for BINF 6310: Dinesh Sambhaji Pradhan, Asmitha Nagajothi Purushotam, Vedant Kulkarni, and Richard Goodier.
