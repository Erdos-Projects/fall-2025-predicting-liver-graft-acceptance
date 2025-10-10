#Microarray Processing
import GEOparse
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

# 1 Import annotation data and microarray expression data

# Load the GSE dataset
gse = GEOparse.get_GEO("GSE145780", destdir="C:/Users/nanth/Desktop/Erdo_work")

# Print general info
print(gse.metadata)

# List all samples
print(gse.gsms.keys())

# Access the GPL
gpl = list(gse.gpls.values())[0]
# Define the columns you want to keep
columns_of_interest = [
    'ID', 
    'Gene Symbol', 
    'Gene Ontology Biological Process', 
    'Gene Ontology Cellular Component',
    'Gene Ontology Molecular Function', 
    'Pathway'
]

# Subset the table
gpl_subset = gpl.table[columns_of_interest]

# Check the first few rows
print(gpl_subset.head())

# Combine all samples into one DataFrame
all_expr = pd.concat([gsm.table['VALUE'] for gsm in gse.gsms.values()], axis=1)
all_expr.columns = gse.gsms.keys()
all_expr['ID_REF'] = list(list(gse.gsms.values())[0].table['ID_REF'])

# Merge expression data with gene symbols
expr_annotated = all_expr.merge(
    gpl_subset,
    left_on='ID_REF',   # column in your expression matrix
    right_on='ID',      # column in the GPL annotation
    how='left'          # keep all rows from expression matrix
)

# Set gene symbol as index
expr_annotated.set_index('Gene Symbol', inplace=True)
expr_annotated.drop(columns=['ID_REF'], inplace=True)

# 7. Save to CSV
expr_annotated.to_csv(r"C:/Users/nanth\Desktop/Erdo_work/GSE145780_expression.csv")

print(expr_annotated.head())

# Extract sample and archetype
sample_info = []

for gsm_name, gsm in gse.gsms.items():
    characteristics = gsm.metadata.get('characteristics_ch1', [])
    
    # Find the cluster info
    cluster = None
    for char in characteristics:
        if "cluster" in char.lower():
            parts = char.split(":", 1)
            if len(parts) > 1:
                cluster = parts[1].strip()  # take the part after colon
    
    sample_info.append({'Sample': gsm_name, 'Cluster': cluster})

sample_df = pd.DataFrame(sample_info)
sample_df.to_csv(r"C:/Users/nanth/Desktop/Erdo_work/GSE145780_metadata.csv", index=False)

# Ensure the Gene Symbol column exists
if 'Gene Symbol' not in expr_annotated.columns:
    raise KeyError("Column 'Gene Symbol' not found in expr_annotated")

# Set 'Gene Symbol' as index for grouping
expr_annotated = expr_annotated.set_index('Gene Symbol')

# Separate numeric (expression) and non-numeric (annotation) columns
expr_values = expr_annotated.select_dtypes(include='number')
annotations = expr_annotated.select_dtypes(exclude='number')

# Collapse expression values by taking the median per gene symbol
expr_collapsed = expr_values.groupby(expr_annotated.index).median()

# For annotations, keep the first occurrence for each gene symbol
annotations_unique = annotations[~annotations.index.duplicated(keep='first')]

# 5. Merge the collapsed expression and annotation data
expr_final = expr_collapsed.merge(
    annotations_unique,
    left_index=True,
    right_index=True,
    how='left'
)

print(f"Shape before: {expr_annotated.shape}, after collapsing: {expr_final.shape}")
expr_final.to_csv(r"C:/Users/nanth\Desktop/Erdo_work/GSE145780_collapsed.csv")


# Quality control of RAW data
# Keep only numeric columns (expression values)
expr_values = expr_final.select_dtypes(include='number')

# Perform PCA on transposed data (samples as rows)
pca = PCA()
pca_raw = pca.fit_transform(expr_values.T)

# Calculate variance explained by each PC
percent_var_raw = np.round(100 * pca.explained_variance_ratio_, 1)

# Create a DataFrame for plotting
dataGG_raw = pd.DataFrame({
    "PC1_raw": pca_raw[:, 0],
    "PC2_raw": pca_raw[:, 1],
    # Replace this with your actual cluster info — for example:
    # clusters = sample_metadata["Cluster"]
    # or extract from gse.gsms metadata
    "Cluster": sample_df["Cluster"]  
})

# Plot using seaborn
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=dataGG_raw,
    x="PC1_raw",
    y="PC2_raw",
    hue="Cluster",
    style="Cluster",
    s=80,
    palette=["darkorange", "dodgerblue", "forestgreen", "purple"]
)

plt.title("PCA plot of log-transformed RMA expression data", fontsize=14, ha='center')
plt.xlabel(f"PC1, VarExp: {percent_var_raw[0]}%")
plt.ylabel(f"PC2, VarExp: {percent_var_raw[1]}%")
plt.tight_layout()
plt.savefig(r"C:/Users/nanth/Desktop/Erdo_work/GSE145780_rawPCA.png", dpi=300)
plt.show()

## Box plot of log2 intensities RAW data
# Keep only numeric expression values
expr_values = expr_final.select_dtypes(include='number')

# Compute the median expression for each gene
gene_medians = expr_values.median(axis=1)

# Compute RLE (subtract gene median from each sample)
rle_data = expr_values.subtract(gene_medians, axis=0)

# Convert to long format for plotting
rle_long = rle_data.reset_index().melt(id_vars='Gene Symbol',
                                       var_name='Sample',
                                       value_name='RLE')

# Merge sample metadata (Cluster/Archetype)
rle_long = rle_long.merge(sample_df, left_on='Sample', right_on='Sample', how='left')

# Plot RLE using seaborn
plt.figure(figsize=(16, 6))
sns.boxplot(
    data=rle_long,
    x='Sample',
    y='RLE',
    hue='Cluster',  # or 'Archetype' if you want
    dodge=False,    # one box per sample
    showfliers=False,
    palette="Set2"
)

plt.ylim(-2, 2)
plt.xticks(rotation=90, fontsize=5)
plt.title("RLE Plot of Log2-RMA Expression Data", fontsize=14)
plt.xlabel("Samples")
plt.ylabel("Relative Log Expression (RLE)")
plt.legend(title='Cluster', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(r"C:/Users/nanth/Desktop/Erdo_work/GSE145780_BoxIntensities.png", dpi=300)
plt.show()

## Heatmap clustering, sample-to-sample distances
# Keep only numeric columns (expression values)
expr_numeric = expr_final.select_dtypes(include=np.number)

# Compute Manhattan distances between samples
dists = pd.DataFrame(
    squareform(pdist(expr_numeric.T, metric='cityblock')),
    index=expr_numeric.columns,
    columns=expr_numeric.columns
)

# Optional: mask the diagonal
np.fill_diagonal(dists.values, np.nan)

# Cluster the samples (hierarchical clustering)
linkage_matrix = linkage(pdist(expr_numeric.T, metric='cityblock'), method='average')
dendro = dendrogram(linkage_matrix, labels=expr_numeric.columns, no_plot=True)
ordered_cols = dendro['ivl']
dists_ordered = dists.loc[ordered_cols, ordered_cols]

# Cluster annotation
annotation_for_heatmap = sample_df.set_index('Sample')
cluster_palette = {"1": "darkorange", "2": "dodgerblue", "3": "forestgreen", "4": "purple"}
row_colors = annotation_for_heatmap['Cluster'].map(cluster_palette).loc[ordered_cols]

# Plot heatmap
sns.set(style='white')
plt.figure(figsize=(10, 8))
sns.heatmap(
    dists_ordered,
    cmap="YlOrRd_r",
    linewidths=0.5,
    linecolor='gray',
    square=True,
    cbar_kws={'label': 'Manhattan distance'},
    yticklabels=True,
    xticklabels=False
)

# Add colored bars on the left to indicate cluster
for ytick, color in zip(range(len(row_colors)), row_colors):
    plt.gca().add_patch(plt.Rectangle((-0.5, ytick), 0.5, 1, facecolor=color, edgecolor='none'))
plt.savefig(r"C:/Users/nanth/Desktop/Erdo_work/GSE145780_Heatmap.png", dpi=300)
plt.show()

## Filtering on Intensity
# Ensure numeric values only
expr_numeric = expr_final.select_dtypes(include=np.number)

# Compute row medians correctly
row_medians = expr_numeric.median(axis=1).astype(float)  # <- force float

# Plot histogram
plt.figure(figsize=(8, 6))
plt.hist(
    row_medians,
    bins=100,              # make sure enough bins
    color="peachpuff",
    edgecolor="antiquewhite",
    density=True
)
plt.title("Histogram of the median intensities")
plt.xlabel("Median intensities")
plt.xlim(row_medians.min() - 0.5, row_medians.max() + 0.5)

# Add threshold line
man_threshold = 3.5
plt.axvline(x=man_threshold, color="coral", linewidth=2)
plt.savefig(r"C:/Users/nanth/Desktop/Erdo_work/GSE145780_HistoIntensities.png", dpi=300)
plt.show()

# Filtering by intensity
expr_numeric = expr_final.select_dtypes(include=np.number)

# 1. Count number of samples in each cluster
no_of_samples = sample_df['Cluster'].value_counts()
print(no_of_samples)

# 2. Minimum number of samples across clusters (samples_cutoff)
samples_cutoff = no_of_samples.min()
print(f"Samples cutoff: {samples_cutoff}")

# 3. Filter genes: must have expression > man_threshold in at least samples_cutoff samples
man_threshold = 3.5
idx_man_threshold = (expr_numeric > man_threshold).sum(axis=1) >= samples_cutoff
print(idx_man_threshold.value_counts())  # True = genes kept, False = filtered out

# 4. Subset filtered genes
chang_manfiltered = expr_numeric[idx_man_threshold]
print(f"Number of genes after filtering: {chang_manfiltered.shape[0]}")
