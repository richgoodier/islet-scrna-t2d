# Gene Matrix — Read Count Generation

This directory contains the raw read count matrices produced by featureCounts after STAR alignment.

## Generating Count Matrices

After alignment (see [`../pipeline/README.md`](../pipeline/README.md)), run featureCounts on the BAM files:

```bash
# Healthy samples
featureCounts -T 8 -a Homo_sapiens.GRCh38.109.gtf \
  -o featureCounts_healthy/counts.txt \
  aligned_healthy/*.bam

# T2D samples
featureCounts -T 8 -a Homo_sapiens.GRCh38.109.gtf \
  -o featureCounts_t2d/counts.txt \
  aligned_t2d/*.bam
```

## Merging

Run `combine_tables.ipynb` to merge the healthy and T2D count matrices into a single file:

```
counts_combined.txt  — 62,711 genes × 142 samples (raw counts + gene length)
```

This file is not tracked in git due to size (~55 MB). It is the input to RPKM normalization.

## Next Step

Pass `counts_combined.txt` to `../rpkm_values/calculate_rpkm.txt` (R/edgeR) to compute RPKM values.
