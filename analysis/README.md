# Analysis — Local Python Notebooks

Local analysis of normalized gene expression data. These notebooks run on a laptop once the RPKM matrix is available from the HPC pipeline.

## Prerequisites

```bash
conda activate islet-scrna-t2d
jupyter lab
```

The main input file is `../rpkm_values/rpkm_combined.txt` (62,711 genes × 142 samples). This file is not tracked in git — regenerate it by following the pipeline in `../pipeline/README.md` and `../rpkm_values/calculate_rpkm.R`.

## Notebooks

| Notebook | Description |
| --- | --- |
| `ensembl_to_ids.ipynb` | Converts Ensembl gene IDs to HGNC symbols via MyGeneInfo API |
| `differential_expression.ipynb` | Differential expression: log2 fold change + t-test, volcano plot |
| `expression_analysis.ipynb` | Expression scatter plots, boxplots, t-SNE, GSEA heatmaps |
| `ml_classification.ipynb` | PCA (pre/post normalization), t-SNE by donor, logistic regression classifier |

Run in the order listed above. `ensembl_to_ids.ipynb` must run first to generate the gene symbol mapping used by subsequent notebooks.

## Logistic Regression Model (`ml_classification.ipynb`)

Trains a logistic regression classifier (L1 penalty, liblinear solver) to distinguish healthy vs. T2D β-cells using log2-RPKM gene expression features.

- 70/30 train-test split, stratified by condition
- Features standardized with `StandardScaler`
- Top 15 predictive genes identified by absolute coefficient magnitude
- Top predictor: **IRF2BPL** (note: sex imbalance in dataset — 5M:1F healthy, 2M:2F T2D — may influence this ranking)

## Figures

All output figures are saved as `.png` in the [`figures/`](figures/) subdirectory:

- `figures/volcano.png` — differential expression volcano plot
- `figures/tsne_condition.png` / `tsne_donor.png` / `tsne_donor_top50.png` — t-SNE colored by condition and donor
- `figures/predictive_genes.png` — top 15 logistic regression feature importances
- `figures/reproduced_fig4Bi.png`, `reproduced_fig4Bii.png`, `reproduced_figS2C*.png` — reproduced figures from Segerstolpe et al.
