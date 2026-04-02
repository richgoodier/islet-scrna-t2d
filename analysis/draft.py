import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set plotting style
sns.set(style="whitegrid")

# Placeholder list of beta cell-associated genes (Ensembl IDs)
BETA_CELL_GENES = [
    'ENSG00000254647',  # INS
    'ENSG00000106633',  # GCK
    'ENSG00000139515',  # PDX1
    'ENSG00000163581',  # SLC2A2
    'ENSG00000115565'   # NKX6-1
]

# 1. Load and Prepare the Data
def load_and_prepare_data(file_path, healthy_label='healthy', t2d_label='t2d'):
    """Loads the data, sets index, and preprocesses."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found in {os.getcwd()}")

    print(f"Loading file: {file_path}")
    df = pd.read_csv(file_path, sep='\t', index_col=0, low_memory=False)
    print(f"Data shape: {df.shape}")
    print(f"First few rows:\n{df.head()}")

    conditions = df.iloc[0].to_list()
    print(f"Conditions: {conditions}")
    df = df.iloc[1:]
    df = df.dropna(axis=0, how='any')  # Remove rows with any NaN
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.dropna()
    print(f"After cleaning, data shape: {df.shape}")

    # Ensure the correct labels are used
    if healthy_label not in conditions or t2d_label not in conditions:
        raise ValueError(f"Conditions must contain '{healthy_label}' and '{t2d_label}'. Found: {set(conditions)}")

    return df, conditions, healthy_label, t2d_label

# 2. Differential Expression Analysis (using t-test)
def perform_t_test(df, conditions, healthy_label, t2d_label):
    """Performs t-tests and calculates fold changes."""
    group1_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == healthy_label]]
    group2_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == t2d_label]]

    mean_group1 = group1_data.mean(axis=1)
    mean_group2 = group2_data.mean(axis=1)

    pseudocount = 1e-10
    fold_change = np.log2((mean_group2 + pseudocount) / (mean_group1 + pseudocount))

    t_stats, p_values = ttest_ind(group1_data, group2_data, axis=1)

    results_df = pd.DataFrame({'log2FoldChange': fold_change, 'pvalue': p_values}).dropna()
    return results_df

def apply_fdr_correction(results_df):
    """Applies Benjamini-Hochberg FDR correction."""
    results_df['padj'] = fdrcorrection(results_df['pvalue'])[1]
    results_df['neg_log10_padj'] = -np.log10(results_df['padj'].replace(0, np.nextafter(0, 1)))
    return results_df

def filter_significant_genes(results_df, log2fc_threshold=1, padj_threshold=0.05):
    """Filters for significantly differentially expressed genes."""
    upregulated_genes = results_df[
        (results_df['log2FoldChange'] > log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    downregulated_genes = results_df[
        (results_df['log2FoldChange'] < -log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    return upregulated_genes, downregulated_genes

# 3. Prepare for GSEA
def prepare_for_gsea(results_df, rank_by='log2FoldChange'):
    """Prepares a ranked list for GSEA."""
    ranked_genes = results_df.sort_values(by=rank_by, ascending=False)
    ranked_list = ranked_genes[[rank_by]]
    return ranked_list

# Plotting Functions
def plot_upregulated_genes(upregulated, filename='upregulated_genes.png'):
    """Plots top 10 upregulated genes."""
    if upregulated.empty:
        print("No upregulated genes to plot.")
        return
    top_n = upregulated.head(10)
    plt.figure(figsize=(10, 6))
    sns.barplot(x='log2FoldChange', y=top_n.index, data=top_n, hue='log2FoldChange', palette='Blues_d')
    plt.title('Top 10 Upregulated Genes (log2FC > 1, padj < 0.05)')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")

def plot_downregulated_genes(downregulated, filename='downregulated_genes.png'):
    """Plots top 10 downregulated genes."""
    if downregulated.empty:
        print("No downregulated genes to plot.")
        return
    top_n = downregulated.head(10)
    plt.figure(figsize=(10, 6))
    sns.barplot(x='log2FoldChange', y=top_n.index, data=top_n, hue='log2FoldChange', palette='Reds_d')
    plt.title('Top 10 Downregulated Genes (log2FC < -1, padj < 0.05)')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")

def plot_volcano(results_df, filename='volcano_plot.png'):
    """Plots a volcano plot with upregulated, downregulated, and not significant genes."""
    try:
        plt.figure(figsize=(10, 8))
        # Categorize genes
        results_df['category'] = 'Not Significant'
        results_df.loc[(results_df['log2FoldChange'] > 1) & (results_df['padj'] < 0.05), 'category'] = 'Upregulated'
        results_df.loc[(results_df['log2FoldChange'] < -1) & (results_df['padj'] < 0.05), 'category'] = 'Downregulated'

        # Plot
        sns.scatterplot(x='log2FoldChange', y='neg_log10_padj', data=results_df,
                        hue='category', palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not Significant': 'grey'},
                        alpha=0.6)

        plt.title('Volcano Plot (t2d vs healthy)')
        plt.xlabel('log2 Fold Change (t2d vs healthy)')
        plt.ylabel('-log10 Adjusted P-value')
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=1, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-1, color='black', linestyle='--', alpha=0.5)

        # Label only beta cell genes that are significant
        significant = (results_df['padj'] < 0.05) & (abs(results_df['log2FoldChange']) > 1)
        for i, row in results_df[significant & results_df.index.isin(BETA_CELL_GENES)].iterrows():
            plt.text(row['log2FoldChange'], row['neg_log10_padj'], i, fontsize=8, ha='center', va='bottom')

        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_volcano: {e}")

def plot_pvalue_histogram(results_df, filename='pvalue_histogram.png'):
    """Plots histogram of neg_log10_padj."""
    try:
        plt.figure(figsize=(8, 6))
        sns.histplot(results_df['neg_log10_padj'], bins=50, color='purple')
        plt.title('Distribution of -log10 Adjusted P-values')
        plt.xlabel('-log10 Adjusted P-value')
        plt.ylabel('Count')
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_pvalue_histogram: {e}")

def plot_gsea_heatmap(df, conditions, gene_list, filename='gsea_heatmap.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots a heatmap for the given gene list."""
    try:
        if not gene_list:
            print("No genes provided for GSEA heatmap.")
            return

        valid_genes = [g for g in gene_list if g in df.index]
        if not valid_genes:
            print("No valid genes found in data for GSEA heatmap.")
            return

        subset_df = df.loc[valid_genes, :]

        condition_colors = {healthy_label: 'blue', t2d_label: 'red'}
        col_colors = [condition_colors.get(cond, 'grey') for cond in conditions]

        plt.figure(figsize=(len(conditions) * 0.8, len(valid_genes) * 0.5))
        sns.clustermap(subset_df, cmap='RdBu_r', standard_scale=1, col_colors=col_colors,
                       cbar_kws={'label': 'Z-score'}, xticklabels=True, yticklabels=True)
        abs_path = os.path.abspath(filename)
        plt.savefig(filename, bbox_inches='tight')
        plt.close()
        print(f"Saved heatmap: {abs_path}")
    except Exception as e:
        print(f"Error in plot_gsea_heatmap: {e}")

# Main Execution
if __name__ == '__main__':
    file_path = 'rpkm_combined_with_classes.txt'

    try:
        df, conditions, healthy_label, t2d_label = load_and_prepare_data(file_path, healthy_label='healthy', t2d_label='t2d')

        results_df = perform_t_test(df, conditions, healthy_label, t2d_label)
        results_df = apply_fdr_correction(results_df)

        upregulated, downregulated = filter_significant_genes(results_df)

        print("\nUpregulated Genes (t2d vs healthy):")
        print(upregulated.head())
        plot_upregulated_genes(upregulated)

        print("\nDownregulated Genes (t2d vs healthy):")
        print(downregulated.head())
        plot_downregulated_genes(downregulated)

        ranked_list_fc = prepare_for_gsea(results_df, rank_by='log2FoldChange')
        ranked_list_pvalue = prepare_for_gsea(results_df, rank_by='neg_log10_padj')

        print("\nRanked List for GSEA (Fold Change):")
        print(ranked_list_fc.head())
        plot_volcano(results_df)

        print("\nRanked List for GSEA (P-value):")
        print(ranked_list_pvalue.head())
        plot_pvalue_histogram(results_df)

        ranked_list_fc.to_csv("ranked_genes_fold_change.rnk", sep='\t', header=True)
        ranked_list_pvalue.to_csv("ranked_genes_pvalue.rnk", sep='\t', header=True)
        print("Saved GSEA ranked lists: ranked_genes_fold_change.rnk, ranked_genes_pvalue.rnk")

        example_gene_list = upregulated.index[:20].tolist() if not upregulated.empty else []
        plot_gsea_heatmap(df, conditions, example_gene_list, healthy_label=healthy_label, t2d_label=t2d_label)

    except Exception as e:
        print(f"Error: {e}")


        ###

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set plotting style
sns.set(style="whitegrid")

# 1. Load and Prepare the Data
def load_and_prepare_data(file_path, healthy_label='healthy', t2d_label='t2d'):
    """Loads the data, sets index, and preprocesses."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found in {os.getcwd()}")

    print(f"Loading file: {file_path}")
    df = pd.read_csv(file_path, sep='\t', index_col=0, low_memory=False)
    print(f"Data shape: {df.shape}")
    print(f"First few rows:\n{df.head()}")

    conditions = df.iloc[0].to_list()
    print(f"Conditions: {conditions}")
    df = df.iloc[1:]
    df = df.dropna(axis=0, how='any')  # Remove rows with any NaN
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.dropna()
    print(f"After cleaning, data shape: {df.shape}")

    # Ensure the correct labels are used
    if healthy_label not in conditions or t2d_label not in conditions:
        raise ValueError(f"Conditions must contain '{healthy_label}' and '{t2d_label}'. Found: {set(conditions)}")

    return df, conditions, healthy_label, t2d_label

# 2. Differential Expression Analysis (using t-test)
def perform_t_test(df, conditions, healthy_label, t2d_label):
    """Performs t-tests and calculates fold changes."""
    group1_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == healthy_label]]
    group2_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == t2d_label]]

    mean_group1 = group1_data.mean(axis=1)
    mean_group2 = group2_data.mean(axis=1)

    pseudocount = 1e-10
    fold_change = np.log2((mean_group2 + pseudocount) / (mean_group1 + pseudocount))

    t_stats, p_values = ttest_ind(group1_data, group2_data, axis=1)

    results_df = pd.DataFrame({'log2FoldChange': fold_change, 'pvalue': p_values}).dropna()
    return results_df

def apply_fdr_correction(results_df):
    """Applies Benjamini-Hochberg FDR correction."""
    results_df['padj'] = fdrcorrection(results_df['pvalue'])[1]
    results_df['neg_log10_padj'] = -np.log10(results_df['padj'].replace(0, np.nextafter(0, 1)))
    return results_df

def filter_significant_genes(results_df, log2fc_threshold=1, padj_threshold=0.05):
    """Filters for significantly differentially expressed genes."""
    upregulated_genes = results_df[
        (results_df['log2FoldChange'] > log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    downregulated_genes = results_df[
        (results_df['log2FoldChange'] < -log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    return upregulated_genes, downregulated_genes

# 3. Prepare for GSEA
def prepare_for_gsea(results_df, rank_by='log2FoldChange'):
    """Prepares a ranked list for GSEA."""
    ranked_genes = results_df.sort_values(by=rank_by, ascending=False)
    ranked_list = ranked_genes[[rank_by]]
    return ranked_list

# Plotting Functions
def plot_upregulated_genes(upregulated, filename='upregulated_genes.png'):
    """Plots top 10 upregulated genes."""
    if upregulated.empty:
        print("No upregulated genes to plot.")
        return
    top_n = upregulated.head(10)
    plt.figure(figsize=(10, 6))
    sns.barplot(x='log2FoldChange', y=top_n.index, data=top_n, hue='log2FoldChange', palette='Blues_d')
    plt.title('Top 10 Upregulated Genes (log2FC > 1, padj < 0.05)')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")

def plot_downregulated_genes(downregulated, filename='downregulated_genes.png'):
    """Plots top 10 downregulated genes."""
    if downregulated.empty:
        print("No downregulated genes to plot.")
        return
    top_n = downregulated.head(10)
    plt.figure(figsize=(10, 6))
    sns.barplot(x='log2FoldChange', y=top_n.index, data=top_n, hue='log2FoldChange', palette='Reds_d')
    plt.title('Top 10 Downregulated Genes (log2FC < -1, padj < 0.05)')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")

def plot_volcano(results_df, filename='volcano_plot.png'):
    """Plots a volcano plot."""
    try:
        plt.figure(figsize=(10, 8))
        significant = (results_df['padj'] < 0.05) & (abs(results_df['log2FoldChange']) > 1)
        results_df['significant'] = ['t2d vs healthy' if sig else 'Not Significant' for sig in significant]

        sns.scatterplot(x='log2FoldChange', y='neg_log10_padj', data=results_df,
                        hue='significant', palette={'t2d vs healthy': 'red', 'Not Significant': 'grey'}, alpha=0.6)

        plt.title('Volcano Plot (t2d vs healthy)')
        plt.xlabel('log2 Fold Change (t2d vs healthy)')
        plt.ylabel('-log10 Adjusted P-value')
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=1, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-1, color='black', linestyle='--', alpha=0.5)

        # Add gene names for significant points
        for i, row in results_df[significant].iterrows():
            plt.text(row['log2FoldChange'], row['neg_log10_padj'], i, fontsize=8, ha='center', va='bottom')

        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_volcano: {e}")

def plot_pvalue_histogram(results_df, filename='pvalue_histogram.png'):
    """Plots histogram of neg_log10_padj."""
    try:
        plt.figure(figsize=(8, 6))
        sns.histplot(results_df['neg_log10_padj'], bins=50, color='purple')
        plt.title('Distribution of -log10 Adjusted P-values')
        plt.xlabel('-log10 Adjusted P-value')
        plt.ylabel('Count')
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_pvalue_histogram: {e}")

def plot_gsea_heatmap(df, conditions, gene_list, filename='gsea_heatmap.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots a heatmap for the given gene list."""
    try:
        if not gene_list:
            print("No genes provided for GSEA heatmap.")
            return

        valid_genes = [g for g in gene_list if g in df.index]
        if not valid_genes:
            print("No valid genes found in data for GSEA heatmap.")
            return

        subset_df = df.loc[valid_genes, :]

        condition_colors = {healthy_label: 'blue', t2d_label: 'red'}
        col_colors = [condition_colors.get(cond, 'grey') for cond in conditions]

        plt.figure(figsize=(len(conditions) * 0.8, len(valid_genes) * 0.5))
        sns.clustermap(subset_df, cmap='RdBu_r', standard_scale=1, col_colors=col_colors,
                       cbar_kws={'label': 'Z-score'}, xticklabels=True, yticklabels=True)
        abs_path = os.path.abspath(filename)
        plt.savefig(filename, bbox_inches='tight')
        plt.close()
        print(f"Saved heatmap: {abs_path}")
    except Exception as e:
        print(f"Error in plot_gsea_heatmap: {e}")

# Main Execution
if __name__ == '__main__':
    file_path = 'rpkm_combined_with_classes.txt'

    try:
        df, conditions, healthy_label, t2d_label = load_and_prepare_data(file_path, healthy_label='healthy', t2d_label='t2d')

        results_df = perform_t_test(df, conditions, healthy_label, t2d_label)
        results_df = apply_fdr_correction(results_df)

        upregulated, downregulated = filter_significant_genes(results_df)

        print("\nUpregulated Genes (t2d vs healthy):")
        print(upregulated.head())
        plot_upregulated_genes(upregulated)

        print("\nDownregulated Genes (t2d vs healthy):")
        print(downregulated.head())
        plot_downregulated_genes(downregulated)

        ranked_list_fc = prepare_for_gsea(results_df, rank_by='log2FoldChange')
        ranked_list_pvalue = prepare_for_gsea(results_df, rank_by='neg_log10_padj')

        print("\nRanked List for GSEA (Fold Change):")
        print(ranked_list_fc.head())
        plot_volcano(results_df)

        print("\nRanked List for GSEA (P-value):")
        print(ranked_list_pvalue.head())
        plot_pvalue_histogram(results_df)

        ranked_list_fc.to_csv("ranked_genes_fold_change.rnk", sep='\t', header=True)
        ranked_list_pvalue.to_csv("ranked_genes_pvalue.rnk", sep='\t', header=True)
        print("Saved GSEA ranked lists: ranked_genes_fold_change.rnk, ranked_genes_pvalue.rnk")

        example_gene_list = upregulated.index[:20].tolist() if not upregulated.empty else []
        plot_gsea_heatmap(df, conditions, example_gene_list, healthy_label=healthy_label, t2d_label=t2d_label)

    except Exception as e:
        print(f"Error: {e}")


###
"""With upregulated and down regulated volcano plot"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set plotting style
sns.set(style="whitegrid")

# Placeholder list of beta cell-associated genes (Ensembl IDs)
BETA_CELL_GENES = [
    'ENSG00000254647',  # INS
    'ENSG00000106633',  # GCK
    'ENSG00000139515',  # PDX1
    'ENSG00000163581',  # SLC2A2
    'ENSG00000115565'  # NKX6-1
]


# 1. Load and Prepare the Data
def load_and_prepare_data(file_path, healthy_label='healthy', t2d_label='t2d'):
    """Loads the data, sets index, and preprocesses."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found in {os.getcwd()}")

    print(f"Loading file: {file_path}")
    df = pd.read_csv(file_path, sep='\t', index_col=0, low_memory=False)
    print(f"Data shape: {df.shape}")
    print(f"First few rows:\n{df.head()}")

    conditions = df.iloc[0].to_list()
    print(f"Conditions: {conditions}")
    df = df.iloc[1:]
    df = df.dropna(axis=0, how='any')  # Remove rows with any NaN
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.dropna()
    print(f"After cleaning, data shape: {df.shape}")

    # Ensure the correct labels are used
    if healthy_label not in conditions or t2d_label not in conditions:
        raise ValueError(f"Conditions must contain '{healthy_label}' and '{t2d_label}'. Found: {set(conditions)}")

    return df, conditions, healthy_label, t2d_label


# 2. Differential Expression Analysis (using t-test)
def perform_t_test(df, conditions, healthy_label, t2d_label):
    """Performs t-tests and calculates fold changes."""
    group1_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == healthy_label]]
    group2_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == t2d_label]]

    mean_group1 = group1_data.mean(axis=1)
    mean_group2 = group2_data.mean(axis=1)

    pseudocount = 1e-10
    fold_change = np.log2((mean_group2 + pseudocount) / (mean_group1 + pseudocount))

    t_stats, p_values = ttest_ind(group1_data, group2_data, axis=1)

    results_df = pd.DataFrame({'log2FoldChange': fold_change, 'pvalue': p_values}).dropna()
    return results_df


def apply_fdr_correction(results_df):
    """Applies Benjamini-Hochberg FDR correction."""
    results_df['padj'] = fdrcorrection(results_df['pvalue'])[1]
    results_df['neg_log10_padj'] = -np.log10(results_df['padj'].replace(0, np.nextafter(0, 1)))
    return results_df


def filter_significant_genes(results_df, log2fc_threshold=1, padj_threshold=0.05):
    """Filters for significantly differentially expressed genes."""
    upregulated_genes = results_df[
        (results_df['log2FoldChange'] > log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    downregulated_genes = results_df[
        (results_df['log2FoldChange'] < -log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    return upregulated_genes, downregulated_genes


# 3. Prepare for GSEA
def prepare_for_gsea(results_df, rank_by='log2FoldChange'):
    """Prepares a ranked list for GSEA."""
    ranked_genes = results_df.sort_values(by=rank_by, ascending=False)
    ranked_list = ranked_genes[[rank_by]]
    return ranked_list


# Plotting Functions
def plot_upregulated_genes(upregulated, filename='upregulated_genes.png'):
    """Plots top 10 upregulated genes."""
    if upregulated.empty:
        print("No upregulated genes to plot.")
        return
    top_n = upregulated.head(10)
    plt.figure(figsize=(10, 6))
    sns.barplot(x='log2FoldChange', y=top_n.index, data=top_n, hue='log2FoldChange', palette='Blues_d')
    plt.title('Top 10 Upregulated Genes (log2FC > 1, padj < 0.05)')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")


def plot_downregulated_genes(downregulated, filename='downregulated_genes.png'):
    """Plots top 10 downregulated genes."""
    if downregulated.empty:
        print("No downregulated genes to plot.")
        return
    top_n = downregulated.head(10)
    plt.figure(figsize=(10, 6))
    sns.barplot(x='log2FoldChange', y=top_n.index, data=top_n, hue='log2FoldChange', palette='Reds_d')
    plt.title('Top 10 Downregulated Genes (log2FC < -1, padj < 0.05)')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")


def plot_volcano(results_df, filename='volcano_plot.png'):
    """Plots a volcano plot with upregulated, downregulated, and not significant genes."""
    try:
        plt.figure(figsize=(10, 8))
        # Categorize genes
        results_df['category'] = 'Not Significant'
        results_df.loc[(results_df['log2FoldChange'] > 1) & (results_df['padj'] < 0.05), 'category'] = 'Upregulated'
        results_df.loc[(results_df['log2FoldChange'] < -1) & (results_df['padj'] < 0.05), 'category'] = 'Downregulated'

        # Plot
        sns.scatterplot(x='log2FoldChange', y='neg_log10_padj', data=results_df,
                        hue='category',
                        palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not Significant': 'grey'},
                        alpha=0.6)

        plt.title('Volcano Plot (t2d vs healthy)')
        plt.xlabel('log2 Fold Change (t2d vs healthy)')
        plt.ylabel('-log10 Adjusted P-value')
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=1, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-1, color='black', linestyle='--', alpha=0.5)

        # Label only beta cell genes that are significant
        significant = (results_df['padj'] < 0.05) & (abs(results_df['log2FoldChange']) > 1)
        for i, row in results_df[significant & results_df.index.isin(BETA_CELL_GENES)].iterrows():
            plt.text(row['log2FoldChange'], row['neg_log10_padj'], i, fontsize=8, ha='center', va='bottom')

        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_volcano: {e}")


def plot_pvalue_histogram(results_df, filename='pvalue_histogram.png'):
    """Plots histogram of neg_log10_padj."""
    try:
        plt.figure(figsize=(8, 6))
        sns.histplot(results_df['neg_log10_padj'], bins=50, color='purple')
        plt.title('Distribution of -log10 Adjusted P-values')
        plt.xlabel('-log10 Adjusted P-value')
        plt.ylabel('Count')
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_pvalue_histogram: {e}")


def plot_gsea_heatmap(df, conditions, gene_list, filename='gsea_heatmap.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots a heatmap for the given gene list."""
    try:
        if not gene_list:
            print("No genes provided for GSEA heatmap.")
            return

        valid_genes = [g for g in gene_list if g in df.index]
        if not valid_genes:
            print("No valid genes found in data for GSEA heatmap.")
            return

        subset_df = df.loc[valid_genes, :]

        condition_colors = {healthy_label: 'blue', t2d_label: 'red'}
        col_colors = [condition_colors.get(cond, 'grey') for cond in conditions]

        plt.figure(figsize=(len(conditions) * 0.8, len(valid_genes) * 0.5))
        sns.clustermap(subset_df, cmap='RdBu_r', standard_scale=1, col_colors=col_colors,
                       cbar_kws={'label': 'Z-score'}, xticklabels=True, yticklabels=True)
        abs_path = os.path.abspath(filename)
        plt.savefig(filename, bbox_inches='tight')
        plt.close()
        print(f"Saved heatmap: {abs_path}")
    except Exception as e:
        print(f"Error in plot_gsea_heatmap: {e}")


# Main Execution
if __name__ == '__main__':
    file_path = 'rpkm_combined_with_classes.txt'

    try:
        df, conditions, healthy_label, t2d_label = load_and_prepare_data(file_path, healthy_label='healthy',
                                                                         t2d_label='t2d')

        results_df = perform_t_test(df, conditions, healthy_label, t2d_label)
        results_df = apply_fdr_correction(results_df)

        upregulated, downregulated = filter_significant_genes(results_df)

        print("\nUpregulated Genes (t2d vs healthy):")
        print(upregulated.head())
        plot_upregulated_genes(upregulated)

        print("\nDownregulated Genes (t2d vs healthy):")
        print(downregulated.head())
        plot_downregulated_genes(downregulated)

        ranked_list_fc = prepare_for_gsea(results_df, rank_by='log2FoldChange')
        ranked_list_pvalue = prepare_for_gsea(results_df, rank_by='neg_log10_padj')

        print("\nRanked List for GSEA (Fold Change):")
        print(ranked_list_fc.head())
        plot_volcano(results_df)

        print("\nRanked List for GSEA (P-value):")
        print(ranked_list_pvalue.head())
        plot_pvalue_histogram(results_df)

        ranked_list_fc.to_csv("ranked_genes_fold_change.rnk", sep='\t', header=True)
        ranked_list_pvalue.to_csv("ranked_genes_pvalue.rnk", sep='\t', header=True)
        print("Saved GSEA ranked lists: ranked_genes_fold_change.rnk, ranked_genes_pvalue.rnk")

        example_gene_list = upregulated.index[:20].tolist() if not upregulated.empty else []
        plot_gsea_heatmap(df, conditions, example_gene_list, healthy_label=healthy_label, t2d_label=t2d_label)

    except Exception as e:
        print(f"Error: {e}")


###
"""This code has downregulated genes, upregulated genes, volcano plot, expression boxplot, pca scatter plot, pvalue histogram"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import os

# Set plotting style
sns.set(style="whitegrid")

# Placeholder beta cell genes (Ensembl IDs)
BETA_CELL_GENES = [
    'ENSG00000254647',  # INS
    'ENSG00000106633',  # GCK
    'ENSG00000139515',  # PDX1
    'ENSG00000163581',  # SLC2A2
    'ENSG00000115565'   # NKX6-1
]

# 1. Load and Prepare the Data
def load_and_prepare_data(file_path, healthy_label='healthy', t2d_label='t2d'):
    """Loads the data, sets index, and preprocesses."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found in {os.getcwd()}")

    print(f"Loading file: {file_path}")
    df = pd.read_csv(file_path, sep='\t', index_col=0, low_memory=False)
    print(f"Data shape: {df.shape}")
    print(f"First few rows:\n{df.head()}")

    conditions = df.iloc[0].to_list()
    print(f"Conditions: {conditions}")
    df = df.iloc[1:]
    df = df.dropna(axis=0, how='any')  # Remove rows with any NaN
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.dropna()
    print(f"After cleaning, data shape: {df.shape}")

    if healthy_label not in conditions or t2d_label not in conditions:
        raise ValueError(f"Conditions must contain '{healthy_label}' and '{t2d_label}'. Found: {set(conditions)}")

    return df, conditions, healthy_label, t2d_label

# 2. Differential Expression Analysis
def perform_t_test(df, conditions, healthy_label, t2d_label):
    """Performs t-tests and calculates fold changes."""
    group1_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == healthy_label]]
    group2_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == t2d_label]]

    mean_group1 = group1_data.mean(axis=1)
    mean_group2 = group2_data.mean(axis=1)

    pseudocount = 1e-10
    fold_change = np.log2((mean_group2 + pseudocount) / (mean_group1 + pseudocount))

    t_stats, p_values = ttest_ind(group1_data, group2_data, axis=1)

    results_df = pd.DataFrame({'log2FoldChange': fold_change, 'pvalue': p_values}).dropna()
    return results_df

def apply_fdr_correction(results_df):
    """Applies Benjamini-Hochberg FDR correction."""
    results_df['padj'] = fdrcorrection(results_df['pvalue'])[1]
    results_df['neg_log10_padj'] = -np.log10(results_df['padj'].replace(0, np.nextafter(0, 1)))
    return results_df

def filter_significant_genes(results_df, log2fc_threshold=1, padj_threshold=0.05):
    """Filters for significantly differentially expressed genes."""
    upregulated_genes = results_df[
        (results_df['log2FoldChange'] > log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    downregulated_genes = results_df[
        (results_df['log2FoldChange'] < -log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    return upregulated_genes, downregulated_genes

# 3. Prepare for GSEA
def prepare_for_gsea(results_df, rank_by='log2FoldChange'):
    """Prepares a ranked list for GSEA."""
    ranked_genes = results_df.sort_values(by=rank_by, ascending=False)
    ranked_list = ranked_genes[[rank_by]]
    return ranked_list

# 4. Plotting Functions
def plot_upregulated_genes(upregulated, filename='upregulated_genes.png'):
    """Plots top 10 upregulated genes."""
    if upregulated.empty:
        print("No upregulated genes to plot.")
        return
    top_n = upregulated.head(10)
    plt.figure(figsize=(10, 6))
    sns.barplot(x='log2FoldChange', y=top_n.index, data=top_n, hue='log2FoldChange', palette='Blues_d')
    plt.title('Top 10 Upregulated Genes (log2FC > 1, padj < 0.05)')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")

def plot_downregulated_genes(downregulated, filename='downregulated_genes.png'):
    """Plots top 10 downregulated genes."""
    if downregulated.empty:
        print("No downregulated genes to plot.")
        return
    top_n = downregulated.head(10)
    plt.figure(figsize=(10, 6))
    sns.barplot(x='log2FoldChange', y=top_n.index, data=top_n, hue='log2FoldChange', palette='Reds_d')
    plt.title('Top 10 Downregulated Genes (log2FC < -1, padj < 0.05)')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")

def plot_volcano(results_df, filename='volcano_plot.png'):
    """Plots a volcano plot with upregulated, downregulated, and not significant genes."""
    try:
        plt.figure(figsize=(10, 8))
        results_df['category'] = 'Not Significant'
        results_df.loc[(results_df['log2FoldChange'] > 1) & (results_df['padj'] < 0.05), 'category'] = 'Upregulated'
        results_df.loc[(results_df['log2FoldChange'] < -1) & (results_df['padj'] < 0.05), 'category'] = 'Downregulated'

        sns.scatterplot(x='log2FoldChange', y='neg_log10_padj', data=results_df,
                        hue='category', palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not Significant': 'grey'},
                        alpha=0.6)

        plt.title('Volcano Plot (t2d vs healthy)')
        plt.xlabel('log2 Fold Change (t2d vs healthy)')
        plt.ylabel('-log10 Adjusted P-value')
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=1, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-1, color='black', linestyle='--', alpha=0.5)

        significant = (results_df['padj'] < 0.05) & (abs(results_df['log2FoldChange']) > 1)
        for i, row in results_df[significant & results_df.index.isin(BETA_CELL_GENES)].iterrows():
            plt.text(row['log2FoldChange'], row['neg_log10_padj'], i, fontsize=8, ha='center', va='bottom')

        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_volcano: {e}")

def plot_pvalue_histogram(results_df, filename='pvalue_histogram.png'):
    """Plots histogram of neg_log10_padj."""
    try:
        plt.figure(figsize=(8, 6))
        sns.histplot(results_df['neg_log10_padj'], bins=50, color='purple')
        plt.title('Distribution of -log10 Adjusted P-values')
        plt.xlabel('-log10 Adjusted P-value')
        plt.ylabel('Count')
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_pvalue_histogram: {e}")

def plot_gsea_heatmap(df, conditions, gene_list, filename='gsea_heatmap.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots a heatmap for the given gene list."""
    try:
        if not gene_list:
            print("No genes provided for GSEA heatmap.")
            return

        valid_genes = [g for g in gene_list if g in df.index]
        if not valid_genes:
            print("No valid genes found in data for GSEA heatmap.")
            return

        subset_df = df.loc[valid_genes, :]

        condition_colors = {healthy_label: 'blue', t2d_label: 'red'}
        col_colors = [condition_colors.get(cond, 'grey') for cond in conditions]

        plt.figure(figsize=(len(conditions) * 0.8, len(valid_genes) * 0.5))
        sns.clustermap(subset_df, cmap='RdBu_r', standard_scale=1, col_colors=col_colors,
                       cbar_kws={'label': 'Z-score'}, xticklabels=True, yticklabels=True)
        abs_path = os.path.abspath(filename)
        plt.savefig(filename, bbox_inches='tight')
        plt.close()
        print(f"Saved heatmap: {abs_path}")
    except Exception as e:
        print(f"Error in plot_gsea_heatmap: {e}")

def plot_expression_boxplot(df, conditions, upregulated, downregulated, filename='expression_boxplot.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots boxplot of expression for significant genes."""
    try:
        sig_genes = upregulated.index.tolist() + downregulated.index.tolist()
        if not sig_genes:
            print("No significant genes for boxplot.")
            return

        # Limit to top 10 genes for clarity
        sig_genes = sig_genes[:10]
        subset_df = df.loc[sig_genes, :]

        # Prepare data for boxplot
        data = []
        for gene in subset_df.index:
            for i, sample in enumerate(subset_df.columns):
                data.append({
                    'Gene': gene,
                    'Condition': conditions[i],
                    'Expression': np.log1p(subset_df.loc[gene, sample])  # log1p for RPKM
                })
        plot_df = pd.DataFrame(data)

        plt.figure(figsize=(12, 6))
        sns.boxplot(x='Gene', y='Expression', hue='Condition', data=plot_df,
                    palette={healthy_label: 'blue', t2d_label: 'red'})
        plt.title('Expression of Significant Genes (t2d vs healthy)')
        plt.xlabel('Gene ID')
        plt.ylabel('Log1p(RPKM)')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_expression_boxplot: {e}")

def plot_pca_scatter(df, conditions, filename='pca_scatter.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots PCA scatter of samples."""
    try:
        # Transpose for PCA (samples as rows)
        pca_data = df.T
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(pca_data)
        pca_df = pd.DataFrame({
            'PC1': pca_result[:, 0],
            'PC2': pca_result[:, 1],
            'Condition': conditions
        })

        plt.figure(figsize=(8, 6))
        sns.scatterplot(x='PC1', y='PC2', hue='Condition', data=pca_df,
                        palette={healthy_label: 'blue', t2d_label: 'red'}, alpha=0.8)
        plt.title('PCA of Samples (t2d vs healthy)')
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_pca_scatter: {e}")

# Main Execution
if __name__ == '__main__':
    file_path = 'rpkm_combined_with_classes.txt'

    try:
        df, conditions, healthy_label, t2d_label = load_and_prepare_data(file_path, healthy_label='healthy', t2d_label='t2d')

        results_df = perform_t_test(df, conditions, healthy_label, t2d_label)
        results_df = apply_fdr_correction(results_df)

        upregulated, downregulated = filter_significant_genes(results_df)

        # Save significant genes
        upregulated.to_csv('upregulated_genes.csv')
        downregulated.to_csv('downregulated_genes.csv')
        print("\nUpregulated Genes (t2d vs healthy):")
        print(upregulated.head())
        print(f"Saved upregulated genes to: {os.path.abspath('upregulated_genes.csv')}")
        plot_upregulated_genes(upregulated)

        print("\nDownregulated Genes (t2d vs healthy):")
        print(downregulated.head())
        print(f"Saved downregulated genes to: {os.path.abspath('downregulated_genes.csv')}")
        plot_downregulated_genes(downregulated)

        # Significance summary
        print("\nSignificance Summary:")
        print(f"Upregulated genes ({len(upregulated)}): Higher expression in T2D, potentially linked to stress, inflammation, or compensatory beta cell activity.")
        print(f"Downregulated genes ({len(downregulated)}): Lower expression in T2D, possibly indicating beta cell dysfunction or reduced insulin secretion.")
        print("For detailed pathway analysis, use upregulated_genes.csv and downregulated_genes.csv with tools like DAVID or Enrichr.")

        ranked_list_fc = prepare_for_gsea(results_df, rank_by='log2FoldChange')
        ranked_list_pvalue = prepare_for_gsea(results_df, rank_by='neg_log10_padj')

        print("\nRanked List for GSEA (Fold Change):")
        print(ranked_list_fc.head())
        plot_volcano(results_df)

        print("\nRanked List for GSEA (P-value):")
        print(ranked_list_pvalue.head())
        plot_pvalue_histogram(results_df)

        ranked_list_fc.to_csv("ranked_genes_fold_change.rnk", sep='\t', header=True)
        ranked_list_pvalue.to_csv("ranked_genes_pvalue.rnk", sep='\t', header=True)
        print("Saved GSEA ranked lists: ranked_genes_fold_change.rnk, ranked_genes_pvalue.rnk")

        example_gene_list = upregulated.index[:20].tolist() if not upregulated.empty else []
        plot_gsea_heatmap(df, conditions, example_gene_list, healthy_label=healthy_label, t2d_label=t2d_label)

        # New comparison plots
        plot_expression_boxplot(df, conditions, upregulated, downregulated)
        plot_pca_scatter(df, conditions)

    except Exception as e:
        print(f"Error: {e}")


###
"""Managed to generate some kind of GSEA heatmap."""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import os

# Set plotting style
sns.set(style="whitegrid")

# Placeholder beta cell genes
BETA_CELL_GENES = [
    'ENSG00000254647',  # INS
    'ENSG00000106633',  # GCK
    'ENSG00000139515',  # PDX1
    'ENSG00000163581',  # SLC2A2
    'ENSG00000115565'   # NKX6-1
]

# 1. Load and Prepare the Data
def load_and_prepare_data(file_path, healthy_label='healthy', t2d_label='t2d'):
    """Loads the data, sets index, and preprocesses."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found in {os.getcwd()}")

    print(f"Loading file: {file_path}")
    df = pd.read_csv(file_path, sep='\t', index_col=0, low_memory=False)
    print(f"Data shape: {df.shape}")
    print(f"First few rows:\n{df.head()}")

    conditions = df.iloc[0].to_list()
    print(f"Conditions: {conditions}")
    df = df.iloc[1:]
    df = df.dropna(axis=0, how='any')
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.dropna()
    print(f"After cleaning, data shape: {df.shape}")

    if healthy_label not in conditions or t2d_label not in conditions:
        raise ValueError(f"Conditions must contain '{healthy_label}' and '{t2d_label}'. Found: {set(conditions)}")

    return df, conditions, healthy_label, t2d_label

# 2. Differential Expression Analysis
def perform_t_test(df, conditions, healthy_label, t2d_label):
    """Performs t-tests and calculates fold changes."""
    group1_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == healthy_label]]
    group2_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == t2d_label]]

    mean_group1 = group1_data.mean(axis=1)
    mean_group2 = group2_data.mean(axis=1)

    pseudocount = 1e-10
    fold_change = np.log2((mean_group2 + pseudocount) / (mean_group1 + pseudocount))

    t_stats, p_values = ttest_ind(group1_data, group2_data, axis=1)

    results_df = pd.DataFrame({'log2FoldChange': fold_change, 'pvalue': p_values}).dropna()
    return results_df

def apply_fdr_correction(results_df):
    """Applies Benjamini-Hochberg FDR correction."""
    results_df['padj'] = fdrcorrection(results_df['pvalue'])[1]
    results_df['neg_log10_padj'] = -np.log10(results_df['padj'].replace(0, np.nextafter(0, 1)))
    return results_df

def filter_significant_genes(results_df, log2fc_threshold=1, padj_threshold=0.05):
    """Filters for significantly differentially expressed genes."""
    upregulated_genes = results_df[
        (results_df['log2FoldChange'] > log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    downregulated_genes = results_df[
        (results_df['log2FoldChange'] < -log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    return upregulated_genes, downregulated_genes

# 3. Prepare for GSEA
def prepare_for_gsea(results_df, rank_by='log2FoldChange'):
    """Prepares a ranked list for GSEA."""
    ranked_genes = results_df.sort_values(by=rank_by, ascending=False)
    ranked_list = ranked_genes[[rank_by]]
    return ranked_list

# 4. Plotting Functions
def plot_upregulated_genes(upregulated, filename='upregulated_genes.png'):
    """Plots top 10 upregulated genes."""
    if upregulated.empty:
        print("No upregulated genes to plot.")
        return
    top_n = upregulated.head(10)
    plt.figure(figsize=(10, 6))
    sns.barplot(x='log2FoldChange', y=top_n.index, data=top_n, hue='log2FoldChange', palette='Blues_d')
    plt.title('Top 10 Upregulated Genes (log2FC > 1, padj < 0.05)')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")

def plot_downregulated_genes(downregulated, filename='downregulated_genes.png'):
    """Plots top 10 downregulated genes."""
    if downregulated.empty:
        print("No downregulated genes to plot.")
        return
    top_n = downregulated.head(10)
    plt.figure(figsize=(10, 6))
    sns.barplot(x='log2FoldChange', y=top_n.index, data=top_n, hue='log2FoldChange', palette='Reds_d')
    plt.title('Top 10 Downregulated Genes (log2FC < -1, padj < 0.05)')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")

def plot_volcano(results_df, filename='volcano_plot.png'):
    """Plots a volcano plot with upregulated, downregulated, and not significant genes."""
    try:
        plt.figure(figsize=(10, 8))
        results_df['category'] = 'Not Significant'
        results_df.loc[(results_df['log2FoldChange'] > 1) & (results_df['padj'] < 0.05), 'category'] = 'Upregulated'
        results_df.loc[(results_df['log2FoldChange'] < -1) & (results_df['padj'] < 0.05), 'category'] = 'Downregulated'

        sns.scatterplot(x='log2FoldChange', y='neg_log10_padj', data=results_df,
                        hue='category', palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not Significant': 'grey'},
                        alpha=0.6)

        plt.title('Volcano Plot (t2d vs healthy)')
        plt.xlabel('log2 Fold Change (t2d vs healthy)')
        plt.ylabel('-log10 Adjusted P-value')
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=1, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-1, color='black', linestyle='--', alpha=0.5)

        significant = (results_df['padj'] < 0.05) & (abs(results_df['log2FoldChange']) > 1)
        for i, row in results_df[significant & results_df.index.isin(BETA_CELL_GENES)].iterrows():
            plt.text(row['log2FoldChange'], row['neg_log10_padj'], i, fontsize=8, ha='center', va='bottom')

        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_volcano: {e}")

def plot_pvalue_histogram(results_df, filename='pvalue_histogram.png'):
    """Plots histogram of neg_log10_padj."""
    try:
        plt.figure(figsize=(8, 6))
        sns.histplot(results_df['neg_log10_padj'], bins=50, color='purple')
        plt.title('Distribution of -log10 Adjusted P-values')
        plt.xlabel('-log10 Adjusted P-value')
        plt.ylabel('Count')
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_pvalue_histogram: {e}")

def plot_gsea_heatmap(df, conditions, upregulated, downregulated, filename='gsea_heatmap.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots a readable heatmap for significant genes, comparing healthy vs. T2D."""
    try:
        # Select top 5 upregulated and 5 downregulated genes
        sig_genes = upregulated.head(5).index.tolist() + downregulated.head(5).index.tolist()
        if not sig_genes:
            print("No significant genes for GSEA heatmap.")
            return

        valid_genes = [g for g in sig_genes if g in df.index]
        if not valid_genes:
            print("No valid significant genes found in data for GSEA heatmap.")
            return

        # Log-transform and z-score normalize
        subset_df = df.loc[valid_genes, :].apply(np.log1p)
        subset_df = (subset_df - subset_df.mean()) / subset_df.std()

        # Sort samples: healthy first, then T2D
        sample_order = [s for c, s in sorted(zip(conditions, df.columns), key=lambda x: x[0])]
        subset_df = subset_df[sample_order]
        sorted_conditions = sorted(conditions, key=lambda x: x == healthy_label)

        condition_colors = {healthy_label: 'blue', t2d_label: 'red'}
        col_colors = [condition_colors.get(c, 'grey') for c in sorted_conditions]

        # Plot heatmap
        plt.figure(figsize=(15, max(5, len(valid_genes) * 0.5)))
        sns.heatmap(subset_df, cmap='RdBu_r', center=0, cbar_kws={'label': 'Z-score (log1p)'},
                    xticklabels=sample_order, yticklabels=valid_genes)
        plt.title('GSEA Heatmap: Significant Genes (t2d vs healthy)')
        plt.xlabel('Samples')
        plt.ylabel('Genes')
        plt.xticks(rotation=45, ha='right', fontsize=8)
        plt.yticks(fontsize=10)
        # Add condition colorbar
        for tick, color in zip(plt.gca().get_xticklabels(), col_colors):
            tick.set_color(color)
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename, bbox_inches='tight')
        plt.close()
        print(f"Saved GSEA heatmap: {abs_path}")
    except Exception as e:
        print(f"Error in plot_gsea_heatmap: {e}")

def plot_expression_boxplot(df, conditions, upregulated, downregulated, filename='expression_boxplot.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots boxplot of expression for significant genes."""
    try:
        sig_genes = upregulated.index.tolist() + downregulated.index.tolist()
        if not sig_genes:
            print("No significant genes for boxplot.")
            return

        sig_genes = sig_genes[:10]
        subset_df = df.loc[sig_genes, :]

        data = []
        for gene in subset_df.index:
            for i, sample in enumerate(subset_df.columns):
                data.append({
                    'Gene': gene,
                    'Condition': conditions[i],
                    'Expression': np.log1p(subset_df.loc[gene, sample])
                })
        plot_df = pd.DataFrame(data)

        plt.figure(figsize=(12, 6))
        sns.boxplot(x='Gene', y='Expression', hue='Condition', data=plot_df,
                    palette={healthy_label: 'blue', t2d_label: 'red'})
        plt.title('Expression of Significant Genes (t2d vs healthy)')
        plt.xlabel('Gene ID')
        plt.ylabel('Log1p(RPKM)')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_expression_boxplot: {e}")

def plot_pca_scatter(df, conditions, filename='pca_scatter.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots PCA scatter of samples."""
    try:
        pca_data = df.T
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(pca_data)
        pca_df = pd.DataFrame({
            'PC1': pca_result[:, 0],
            'PC2': pca_result[:, 1],
            'Condition': conditions
        })

        plt.figure(figsize=(8, 6))
        sns.scatterplot(x='PC1', y='PC2', hue='Condition', data=pca_df,
                        palette={healthy_label: 'blue', t2d_label: 'red'}, alpha=0.8)
        plt.title('PCA of Samples (t2d vs healthy)')
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_pca_scatter: {e}")

# Main Execution
if __name__ == '__main__':
    file_path = 'rpkm_combined_with_classes.txt'

    try:
        df, conditions, healthy_label, t2d_label = load_and_prepare_data(file_path, healthy_label='healthy', t2d_label='t2d')

        results_df = perform_t_test(df, conditions, healthy_label, t2d_label)
        results_df = apply_fdr_correction(results_df)

        upregulated, downregulated = filter_significant_genes(results_df)

        upregulated.to_csv('upregulated_genes.csv')
        downregulated.to_csv('downregulated_genes.csv')
        print("\nUpregulated Genes (t2d vs healthy):")
        print(upregulated.head())
        print(f"Saved upregulated genes to: {os.path.abspath('upregulated_genes.csv')}")
        plot_upregulated_genes(upregulated)

        print("\nDownregulated Genes (t2d vs healthy):")
        print(downregulated.head())
        print(f"Saved downregulated genes to: {os.path.abspath('downregulated_genes.csv')}")
        plot_downregulated_genes(downregulated)

        print("\nSignificance Summary:")
        print(f"Upregulated genes ({len(upregulated)}): Higher expression in T2D, potentially linked to stress, inflammation, or compensatory beta cell activity.")
        print(f"Downregulated genes ({len(downregulated)}): Lower expression in T2D, possibly indicating beta cell dysfunction or reduced insulin secretion.")
        print("For detailed pathway analysis, use upregulated_genes.csv and downregulated_genes.csv with tools like DAVID or Enrichr.")

        ranked_list_fc = prepare_for_gsea(results_df, rank_by='log2FoldChange')
        ranked_list_pvalue = prepare_for_gsea(results_df, rank_by='neg_log10_padj')

        print("\nRanked List for GSEA (Fold Change):")
        print(ranked_list_fc.head())
        plot_volcano(results_df)

        print("\nRanked List for GSEA (P-value):")
        print(ranked_list_pvalue.head())
        plot_pvalue_histogram(results_df)

        ranked_list_fc.to_csv("ranked_genes_fold_change.rnk", sep='\t', header=True)
        ranked_list_pvalue.to_csv("ranked_genes_pvalue.rnk", sep='\t', header=True)
        print("Saved GSEA ranked lists: ranked_genes_fold_change.rnk, ranked_genes_pvalue.rnk")

        # GSEA heatmap with significant genes
        plot_gsea_heatmap(df, conditions, upregulated, downregulated)

        plot_expression_boxplot(df, conditions, upregulated, downregulated)
        plot_pca_scatter(df, conditions)

    except Exception as e:
        print(f"Error: {e}")


###
"""This is only going to focus on beta cells along with it's 3 markers (INS, MAFA, PDX1"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import os

# Set plotting style
sns.set(style="whitegrid")

# Beta cell genes (Ensembl IDs for INS, PDX1, MAFA)
BETA_CELL_GENES = [
    'ENSG00000254647',  # INS
    'ENSG00000139515',  # PDX1
    'ENSG00000105697'  # MAFA
]


# 1. Load and Prepare the Data
def load_and_prepare_data(file_path, healthy_label='healthy', t2d_label='t2d'):
    """Loads the data, sets index, and preprocesses."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found in {os.getcwd()}")

    print(f"Loading file: {file_path}")
    df = pd.read_csv(file_path, sep='\t', index_col=0, low_memory=False)
    print(f"Data shape: {df.shape}")
    print(f"First few rows:\n{df.head()}")

    conditions = df.iloc[0].to_list()
    print(f"Conditions: {conditions}")
    df = df.iloc[1:]
    df = df.dropna(axis=0, how='any')
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.dropna()
    print(f"After cleaning, data shape: {df.shape}")

    if healthy_label not in conditions or t2d_label not in conditions:
        raise ValueError(f"Conditions must contain '{healthy_label}' and '{t2d_label}'. Found: {set(conditions)}")

    # Filter for beta cell genes
    valid_beta_genes = [g for g in BETA_CELL_GENES if g in df.index]
    if not valid_beta_genes:
        raise ValueError(f"No beta cell genes {BETA_CELL_GENES} found in data.")
    df = df.loc[valid_beta_genes, :]
    print(f"Filtered for beta cell genes {valid_beta_genes}, data shape: {df.shape}")

    return df, conditions, healthy_label, t2d_label


# 2. Differential Expression Analysis
def perform_t_test(df, conditions, healthy_label, t2d_label):
    """Performs t-tests and calculates fold changes."""
    group1_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == healthy_label]]
    group2_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == t2d_label]]

    mean_group1 = group1_data.mean(axis=1)
    mean_group2 = group2_data.mean(axis=1)

    pseudocount = 1e-10
    fold_change = np.log2((mean_group2 + pseudocount) / (mean_group1 + pseudocount))

    t_stats, p_values = ttest_ind(group1_data, group2_data, axis=1)

    results_df = pd.DataFrame({'log2FoldChange': fold_change, 'pvalue': p_values}, index=df.index).dropna()
    return results_df


def apply_fdr_correction(results_df):
    """Applies Benjamini-Hochberg FDR correction."""
    results_df['padj'] = fdrcorrection(results_df['pvalue'])[1]
    results_df['neg_log10_padj'] = -np.log10(results_df['padj'].replace(0, np.nextafter(0, 1)))
    return results_df


def filter_significant_genes(results_df, log2fc_threshold=1, padj_threshold=0.05):
    """Filters for significantly differentially expressed genes."""
    upregulated_genes = results_df[
        (results_df['log2FoldChange'] > log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    downregulated_genes = results_df[
        (results_df['log2FoldChange'] < -log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    return upregulated_genes, downregulated_genes


# 3. Plotting Functions
def plot_upregulated_genes(upregulated, filename='upregulated_beta_cells.png'):
    """Plots upregulated beta cell genes."""
    if upregulated.empty:
        print("No upregulated beta cell genes to plot.")
        return
    plt.figure(figsize=(8, 4))
    sns.barplot(x='log2FoldChange', y=upregulated.index, data=upregulated, hue='log2FoldChange', palette='Blues_d')
    plt.title('Upregulated Beta Cell Genes (log2FC > 1, padj < 0.05)')
    plt.xlabel('log2 Fold Change (t2d vs healthy)')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")


def plot_downregulated_genes(downregulated, filename='downregulated_beta_cells.png'):
    """Plots downregulated beta cell genes."""
    if downregulated.empty:
        print("No downregulated beta cell genes to plot.")
        return
    plt.figure(figsize=(8, 4))
    sns.barplot(x='log2FoldChange', y=downregulated.index, data=downregulated, hue='log2FoldChange', palette='Reds_d')
    plt.title('Downregulated Beta Cell Genes (log2FC < -1, padj < 0.05)')
    plt.xlabel('log2 Fold Change (t2d vs healthy)')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")


def plot_volcano(results_df, filename='volcano_plot_beta_cells.png'):
    """Plots a volcano plot for beta cell genes."""
    try:
        plt.figure(figsize=(8, 6))
        results_df['category'] = 'Not Significant'
        results_df.loc[(results_df['log2FoldChange'] > 1) & (results_df['padj'] < 0.05), 'category'] = 'Upregulated'
        results_df.loc[(results_df['log2FoldChange'] < -1) & (results_df['padj'] < 0.05), 'category'] = 'Downregulated'

        sns.scatterplot(x='log2FoldChange', y='neg_log10_padj', data=results_df,
                        hue='category',
                        palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not Significant': 'grey'},
                        size=results_df['category'].map(
                            {'Upregulated': 100, 'Downregulated': 100, 'Not Significant': 50}),
                        alpha=0.8)

        plt.title('Volcano Plot: Beta Cell Genes (t2d vs healthy)')
        plt.xlabel('log2 Fold Change (t2d vs healthy)')
        plt.ylabel('-log10 Adjusted P-value')
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=1, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-1, color='black', linestyle='--', alpha=0.5)

        for i, row in results_df.iterrows():
            plt.text(row['log2FoldChange'], row['neg_log10_padj'], i, fontsize=10, ha='center', va='bottom')

        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_volcano: {e}")


def plot_gsea_heatmap(df, conditions, filename='gsea_heatmap_beta_cells.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots a readable heatmap for beta cell genes."""
    try:
        valid_genes = [g for g in df.index]
        if not valid_genes:
            print("No beta cell genes for GSEA heatmap.")
            return

        # Log-transform and z-score normalize
        subset_df = df.apply(np.log1p)
        subset_df = (subset_df - subset_df.mean()) / subset_df.std()

        # Sort samples: healthy first, then T2D
        sample_order = [s for c, s in sorted(zip(conditions, df.columns), key=lambda x: x[0])]
        subset_df = subset_df[sample_order]
        sorted_conditions = sorted(conditions, key=lambda x: x == healthy_label)

        condition_colors = {healthy_label: 'blue', t2d_label: 'red'}
        col_colors = [condition_colors.get(c, 'grey') for c in sorted_conditions]

        # Plot heatmap
        plt.figure(figsize=(15, 4))
        sns.heatmap(subset_df, cmap='RdBu_r', center=0, cbar_kws={'label': 'Z-score (log1p)'},
                    xticklabels=False, yticklabels=valid_genes)
        plt.title('GSEA Heatmap: Beta Cell Genes (t2d vs healthy)')
        plt.xlabel('Samples (Healthy: Blue, T2D: Red)')
        plt.ylabel('Genes')
        plt.yticks(fontsize=12)
        # Add condition colorbar
        ax = plt.gca()
        for i, color in enumerate(col_colors):
            ax.add_patch(
                plt.Rectangle((i, -0.5), 1, 0.2, color=color, transform=ax.get_xaxis_transform(), clip_on=False))
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename, bbox_inches='tight')
        plt.close()
        print(f"Saved GSEA heatmap: {abs_path}")
    except Exception as e:
        print(f"Error in plot_gsea_heatmap: {e}")


def plot_expression_boxplot(df, conditions, filename='expression_boxplot_beta_cells.png', healthy_label='healthy',
                            t2d_label='t2d'):
    """Plots boxplot of expression for beta cell genes."""
    try:
        valid_genes = [g for g in df.index]
        if not valid_genes:
            print("No beta cell genes for boxplot.")
            return

        data = []
        for gene in valid_genes:
            for i, sample in enumerate(df.columns):
                data.append({
                    'Gene': gene,
                    'Condition': conditions[i],
                    'Expression': np.log1p(df.loc[gene, sample])
                })
        plot_df = pd.DataFrame(data)

        plt.figure(figsize=(8, 6))
        sns.boxplot(x='Gene', y='Expression', hue='Condition', data=plot_df,
                    palette={healthy_label: 'blue', t2d_label: 'red'})
        plt.title('Expression of Beta Cell Genes (t2d vs healthy)')
        plt.xlabel('Gene ID')
        plt.ylabel('Log1p(RPKM)')
        plt.xticks(fontsize=10)
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_expression_boxplot: {e}")


def plot_pca_scatter(df, conditions, filename='pca_scatter_beta_cells.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots PCA scatter of samples using beta cell genes."""
    try:
        pca_data = df.T
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(pca_data)
        pca_df = pd.DataFrame({
            'PC1': pca_result[:, 0],
            'PC2': pca_result[:, 1],
            'Condition': conditions
        })

        plt.figure(figsize=(8, 6))
        sns.scatterplot(x='PC1', y='PC2', hue='Condition', data=pca_df,
                        palette={healthy_label: 'blue', t2d_label: 'red'}, alpha=0.8)
        plt.title('PCA of Samples: Beta Cell Genes (t2d vs healthy)')
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_pca_scatter: {e}")


def plot_pvalue_histogram(results_df, filename='pvalue_histogram_beta_cells.png'):
    """Plots histogram of neg_log10_padj for beta cell genes."""
    try:
        plt.figure(figsize=(8, 6))
        sns.histplot(results_df['neg_log10_padj'], bins=10, color='purple')
        plt.title('Distribution of -log10 Adjusted P-values: Beta Cell Genes')
        plt.xlabel('-log10 Adjusted P-value')
        plt.ylabel('Count')
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_pvalue_histogram: {e}")


# Main Execution
if __name__ == '__main__':
    file_path = 'rpkm_combined_with_classes.txt'

    try:
        df, conditions, healthy_label, t2d_label = load_and_prepare_data(file_path, healthy_label='healthy',
                                                                         t2d_label='t2d')

        results_df = perform_t_test(df, conditions, healthy_label, t2d_label)
        results_df = apply_fdr_correction(results_df)

        upregulated, downregulated = filter_significant_genes(results_df)

        results_df.to_csv('beta_cell_genes.csv')
        print("\nBeta Cell Genes Stats (t2d vs healthy):")
        print(results_df)
        print(f"Saved beta cell stats to: {os.path.abspath('beta_cell_genes.csv')}")

        print("\nUpregulated Beta Cell Genes:")
        print(upregulated if not upregulated.empty else "None")
        plot_upregulated_genes(upregulated)

        print("\nDownregulated Beta Cell Genes:")
        print(downregulated if not downregulated.empty else "None")
        plot_downregulated_genes(downregulated)

        print("\nSignificance Summary:")
        for gene in results_df.index:
            padj = results_df.loc[gene, 'padj']
            log2fc = results_df.loc[gene, 'log2FoldChange']
            status = "significant" if padj < 0.05 and abs(log2fc) > 1 else "not significant"
            print(f"{gene}: log2FC={log2fc:.3f}, padj={padj:.3e} ({status})")
        print("Significant genes may indicate altered beta cell function in T2D (e.g., insulin secretion changes).")

        plot_volcano(results_df)
        plot_gsea_heatmap(df, conditions)
        plot_expression_boxplot(df, conditions)
        plot_pca_scatter(df, conditions)
        plot_pvalue_histogram(results_df)

    except Exception as e:
        print(f"Error: {e}")


####
"""This will have everything I hope but specific to markers and including beta cell counts, four markers INS, PPY, SST, GSG"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import os

# Set plotting style
sns.set(style="whitegrid")

# Beta cell genes (Ensembl IDs for INS, PPY, SST, GCG)
BETA_CELL_GENES = [
    'ENSG00000254647',  # INS
    'ENSG00000108849',  # PPY
    'ENSG00000157005',  # SST
    'ENSG00000115263'   # GCG
]

# 1. Load and Prepare the Data
def load_and_prepare_data(file_path, healthy_label='healthy', t2d_label='t2d'):
    """Loads the data, sets index, and preprocesses."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found in {os.getcwd()}")

    print(f"Loading file: {file_path}")
    df = pd.read_csv(file_path, sep='\t', index_col=0, low_memory=False)
    print(f"Data shape: {df.shape}")
    print(f"First few rows:\n{df.head()}")

    conditions = df.iloc[0].to_list()
    print(f"Conditions: {conditions}")
    df = df.iloc[1:]
    df = df.dropna(axis=0, how='any')
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.dropna()
    print(f"After cleaning, data shape: {df.shape}")

    if healthy_label not in conditions or t2d_label not in conditions:
        raise ValueError(f"Conditions must contain '{healthy_label}' and '{t2d_label}'. Found: {set(conditions)}")

    # Filter for beta cell genes
    valid_beta_genes = [g for g in BETA_CELL_GENES if g in df.index]
    if not valid_beta_genes:
        raise ValueError(f"No beta cell genes {BETA_CELL_GENES} found in data.")
    df = df.loc[valid_beta_genes, :]
    print(f"Filtered for beta cell genes {valid_beta_genes}, data shape: {df.shape}")

    return df, conditions, healthy_label, t2d_label

# 2. Differential Expression Analysis
def perform_t_test(df, conditions, healthy_label, t2d_label):
    """Performs t-tests and calculates fold changes."""
    group1_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == healthy_label]]
    group2_data = df.loc[:, [col for i, col in enumerate(df.columns) if conditions[i] == t2d_label]]

    mean_group1 = group1_data.mean(axis=1)
    mean_group2 = group2_data.mean(axis=1)

    pseudocount = 1e-10
    fold_change = np.log2((mean_group2 + pseudocount) / (mean_group1 + pseudocount))

    t_stats, p_values = ttest_ind(group1_data, group2_data, axis=1)

    results_df = pd.DataFrame({'log2FoldChange': fold_change, 'pvalue': p_values}, index=df.index).dropna()
    return results_df

def apply_fdr_correction(results_df):
    """Applies Benjamini-Hochberg FDR correction."""
    results_df['padj'] = fdrcorrection(results_df['pvalue'])[1]
    results_df['neg_log10_padj'] = -np.log10(results_df['padj'].replace(0, np.nextafter(0, 1)))
    return results_df

def filter_significant_genes(results_df, log2fc_threshold=1, padj_threshold=0.05):
    """Filters for significantly differentially expressed genes."""
    upregulated_genes = results_df[
        (results_df['log2FoldChange'] > log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    downregulated_genes = results_df[
        (results_df['log2FoldChange'] < -log2fc_threshold) & (results_df['padj'] < padj_threshold)]
    return upregulated_genes, downregulated_genes

# 3. Plotting Functions
def plot_upregulated_genes(upregulated, filename='upregulated_beta_cells.png'):
    """Plots upregulated beta cell genes."""
    if upregulated.empty:
        print("No upregulated beta cell genes to plot.")
        return
    plt.figure(figsize=(8, 4))
    sns.barplot(x='log2FoldChange', y=upregulated.index, data=upregulated, hue='log2FoldChange', palette='Blues_d')
    plt.title('Upregulated Beta Cell Genes (log2FC > 1, padj < 0.05)')
    plt.xlabel('log2 Fold Change (t2d vs healthy)')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")

def plot_downregulated_genes(downregulated, filename='downregulated_beta_cells.png'):
    """Plots downregulated beta cell genes."""
    if downregulated.empty:
        print("No downregulated beta cell genes to plot.")
        return
    plt.figure(figsize=(8, 4))
    sns.barplot(x='log2FoldChange', y=downregulated.index, data=downregulated, hue='log2FoldChange', palette='Reds_d')
    plt.title('Downregulated Beta Cell Genes (log2FC < -1, padj < 0.05)')
    plt.xlabel('log2 Fold Change (t2d vs healthy)')
    plt.ylabel('Gene ID')
    plt.tight_layout()
    abs_path = os.path.abspath(filename)
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {abs_path}")

def plot_volcano(results_df, filename='volcano_plot_beta_cells.png'):
    """Plots a volcano plot for beta cell genes."""
    try:
        plt.figure(figsize=(8, 6))
        results_df['category'] = 'Not Significant'
        results_df.loc[(results_df['log2FoldChange'] > 1) & (results_df['padj'] < 0.05), 'category'] = 'Upregulated'
        results_df.loc[(results_df['log2FoldChange'] < -1) & (results_df['padj'] < 0.05), 'category'] = 'Downregulated'

        sns.scatterplot(x='log2FoldChange', y='neg_log10_padj', data=results_df,
                        hue='category', palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not Significant': 'grey'},
                        size=results_df['category'].map({'Upregulated': 100, 'Downregulated': 100, 'Not Significant': 50}),
                        alpha=0.8)

        plt.title('Volcano Plot: Beta Cell Genes (t2d vs healthy)')
        plt.xlabel('log2 Fold Change (t2d vs healthy)')
        plt.ylabel('-log10 Adjusted P-value')
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=1, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-1, color='black', linestyle='--', alpha=0.5)

        for i, row in results_df.iterrows():
            plt.text(row['log2FoldChange'], row['neg_log10_padj'], i, fontsize=10, ha='center', va='bottom')

        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_volcano: {e}")

def plot_gsea_heatmap(df, conditions, filename='gsea_heatmap_beta_cells.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots a readable heatmap for beta cell genes."""
    try:
        valid_genes = [g for g in df.index]
        if not valid_genes:
            print("No beta cell genes for GSEA heatmap.")
            return

        # Log-transform and z-score normalize
        subset_df = df.apply(np.log1p)
        subset_df = (subset_df - subset_df.mean()) / subset_df.std()

        # Sort samples: healthy first, then T2D
        sample_order = [s for c, s in sorted(zip(conditions, df.columns), key=lambda x: x[0])]
        subset_df = subset_df[sample_order]
        sorted_conditions = sorted(conditions, key=lambda x: x == healthy_label)

        condition_colors = {healthy_label: 'blue', t2d_label: 'red'}
        col_colors = [condition_colors.get(c, 'grey') for c in sorted_conditions]

        # Plot heatmap
        plt.figure(figsize=(15, 4))
        sns.heatmap(subset_df, cmap='RdBu_r', center=0, cbar_kws={'label': 'Z-score (log1p)'},
                    xticklabels=False, yticklabels=valid_genes)
        plt.title('GSEA Heatmap: Beta Cell Genes (t2d vs healthy)')
        plt.xlabel('Samples (Healthy: Blue, T2D: Red)')
        plt.ylabel('Genes')
        plt.yticks(fontsize=12)
        # Add condition colorbar
        ax = plt.gca()
        for i, color in enumerate(col_colors):
            ax.add_patch(plt.Rectangle((i, -0.5), 1, 0.2, color=color, transform=ax.get_xaxis_transform(), clip_on=False))
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename, bbox_inches='tight')
        plt.close()
        print(f"Saved GSEA heatmap: {abs_path}")
    except Exception as e:
        print(f"Error in plot_gsea_heatmap: {e}")

def plot_expression_boxplot(df, conditions, filename='expression_boxplot_beta_cells.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots boxplot of expression for beta cell genes."""
    try:
        valid_genes = [g for g in df.index]
        if not valid_genes:
            print("No beta cell genes for boxplot.")
            return

        data = []
        for gene in valid_genes:
            for i, sample in enumerate(df.columns):
                data.append({
                    'Gene': gene,
                    'Condition': conditions[i],
                    'Expression': np.log1p(df.loc[gene, sample])
                })
        plot_df = pd.DataFrame(data)

        plt.figure(figsize=(10, 6))
        sns.boxplot(x='Gene', y='Expression', hue='Condition', data=plot_df,
                    palette={healthy_label: 'blue', t2d_label: 'red'})
        plt.title('Expression of Beta Cell Genes (t2d vs healthy)')
        plt.xlabel('Gene ID')
        plt.ylabel('Log1p(RPKM)')
        plt.xticks(fontsize=10)
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_expression_boxplot: {e}")

def plot_pca_scatter(df, conditions, filename='pca_scatter_beta_cells.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots PCA scatter of samples using beta cell genes."""
    try:
        pca_data = df.T
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(pca_data)
        pca_df = pd.DataFrame({
            'PC1': pca_result[:, 0],
            'PC2': pca_result[:, 1],
            'Condition': conditions
        })

        plt.figure(figsize=(8, 6))
        sns.scatterplot(x='PC1', y='PC2', hue='Condition', data=pca_df,
                        palette={healthy_label: 'blue', t2d_label: 'red'}, alpha=0.8)
        plt.title('PCA of Samples: Beta Cell Genes (t2d vs healthy)')
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_pca_scatter: {e}")

def plot_pvalue_histogram(results_df, filename='pvalue_histogram_beta_cells.png'):
    """Plots histogram of neg_log10_padj for beta cell genes."""
    try:
        plt.figure(figsize=(8, 6))
        sns.scatterplot(x=results_df.index, y='neg_log10_padj', data=results_df, color='purple', s=100)
        plt.title('Adjusted P-values: Beta Cell Genes')
        plt.xlabel('Gene ID')
        plt.ylabel('-log10 Adjusted P-value')
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        plt.xticks(rotation=45, fontsize=10)
        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved plot: {abs_path}")
    except Exception as e:
        print(f"Error in plot_pvalue_histogram: {e}")

def plot_beta_cell_summary(df, conditions, results_df, filename='beta_cell_summary.png', healthy_label='healthy', t2d_label='t2d'):
    """Plots beta cell gene expression with sample counts."""
    try:
        valid_genes = [g for g in df.index]
        if not valid_genes:
            print("No beta cell genes for summary plot.")
            return

        # Calculate mean log1p expression
        data = []
        for gene in valid_genes:
            healthy_expr = np.log1p(df.loc[gene, [c for i, c in enumerate(df.columns) if conditions[i] == healthy_label]]).mean()
            t2d_expr = np.log1p(df.loc[gene, [c for i, c in enumerate(df.columns) if conditions[i] == t2d_label]]).mean()
            data.append({'Gene': gene, 'Condition': healthy_label, 'Expression': healthy_expr})
            data.append({'Gene': gene, 'Condition': t2d_label, 'Expression': t2d_expr})
        plot_df = pd.DataFrame(data)

        # Count samples with expression > 1 (proxy for beta cell presence)
        counts = []
        for gene in valid_genes:
            healthy_count = (df.loc[gene, [c for i, c in enumerate(df.columns) if conditions[i] == healthy_label]] > 1).sum()
            t2d_count = (df.loc[gene, [c for i, c in enumerate(df.columns) if conditions[i] == t2d_label]] > 1).sum()
            counts.append({'Gene': gene, 'Condition': healthy_label, 'Count': healthy_count})
            counts.append({'Gene': gene, 'Condition': t2d_label, 'Count': t2d_count})
        count_df = pd.DataFrame(counts)

        # Plot
        plt.figure(figsize=(10, 6))
        bar = sns.barplot(x='Gene', y='Expression', hue='Condition', data=plot_df,
                          palette={healthy_label: 'blue', t2d_label: 'red'})
        plt.title('Beta Cell Gene Expression (t2d vs healthy)')
        plt.xlabel('Gene ID')
        plt.ylabel('Mean Log1p(RPKM)')
        plt.xticks(fontsize=10)

        # Add count annotations
        for i, gene in enumerate(valid_genes):
            healthy_count = count_df[(count_df['Gene'] == gene) & (count_df['Condition'] == healthy_label)]['Count'].iloc[0]
            t2d_count = count_df[(count_df['Gene'] == gene) & (count_df['Condition'] == t2d_label)]['Count'].iloc[0]
            plt.text(i - 0.2, plot_df[(plot_df['Gene'] == gene) & (plot_df['Condition'] == healthy_label)]['Expression'].iloc[0] + 0.1,
                     f'n={healthy_count}', color='blue', ha='center')
            plt.text(i + 0.2, plot_df[(plot_df['Gene'] == gene) & (plot_df['Condition'] == t2d_label)]['Expression'].iloc[0] + 0.1,
                     f'n={t2d_count}', color='red', ha='center')

        plt.tight_layout()
        abs_path = os.path.abspath(filename)
        plt.savefig(filename)
        plt.close()
        print(f"Saved beta cell summary: {abs_path}")
    except Exception as e:
        print(f"Error in plot_beta_cell_summary: {e}")

# Main Execution
if __name__ == '__main__':
    file_path = 'rpkm_combined_with_classes.txt'

    try:
        df, conditions, healthy_label, t2d_label = load_and_prepare_data(file_path, healthy_label='healthy', t2d_label='t2d')

        results_df = perform_t_test(df, conditions, healthy_label, t2d_label)
        results_df = apply_fdr_correction(results_df)

        upregulated, downregulated = filter_significant_genes(results_df)

        results_df.to_csv('beta_cell_genes.csv')
        print("\nBeta Cell Genes Stats (t2d vs healthy):")
        print(results_df)
        print(f"Saved beta cell stats to: {os.path.abspath('beta_cell_genes.csv')}")

        print("\nUpregulated Beta Cell Genes:")
        print(upregulated if not upregulated.empty else "None")
        plot_upregulated_genes(upregulated)

        print("\nDownregulated Beta Cell Genes:")
        print(downregulated if not downregulated.empty else "None")
        plot_downregulated_genes(downregulated)

        print("\nSignificance Summary:")
        for gene in results_df.index:
            padj = results_df.loc[gene, 'padj']
            log2fc = results_df.loc[gene, 'log2FoldChange']
            status = "significant" if padj < 0.05 and abs(log2fc) > 1 else "not significant"
            print(f"{gene}: log2FC={log2fc:.3f}, padj={padj:.3e} ({status})")
        print("Significant genes may indicate altered beta cell function in T2D (e.g., hormone secretion changes).")

        plot_volcano(results_df)
        plot_gsea_heatmap(df, conditions)
        plot_expression_boxplot(df, conditions)
        plot_pca_scatter(df, conditions)
        plot_pvalue_histogram(results_df)
        plot_beta_cell_summary(df, conditions, results_df)

    except Exception as e:
        print(f"Error: {e}")

###
""""""
