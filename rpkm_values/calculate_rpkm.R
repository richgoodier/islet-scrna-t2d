# RPKM normalization using edgeR
#
# Run this script on the HPC cluster after featureCounts has produced
# the combined count matrix. Requires an interactive session with R:
#
#   srun --pty --partition=courses --mem=16G -t 4:00:00 bash
#   conda activate islet-scrna-t2d
#   module load R
#   Rscript calculate_rpkm.R

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("edgeR")

library(edgeR)

# Load the combined count matrix produced by gene_matrix/combine_tables.ipynb
counts <- read.table("../gene_matrix/counts_combined.txt", header = TRUE, row.names = 1)

# Extract gene lengths (used for RPKM normalization)
gene_lengths <- counts$Length

# Select only sample count columns (drop featureCounts metadata columns)
count_matrix <- counts[, !(colnames(counts) %in% c("Geneid", "Chr", "Start", "End", "Strand", "Length"))]

# Calculate RPKM
y <- DGEList(counts = count_matrix)
rpkm_values <- rpkm(y, gene.length = gene_lengths)

# Write output
write.table(rpkm_values, "rpkm_combined.txt", sep = "\t", quote = FALSE)
