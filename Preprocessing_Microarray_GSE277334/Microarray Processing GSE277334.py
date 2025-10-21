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
gse = GEOparse.get_GEO("GSE277334", destdir="C:/Users/nanth/Desktop/Erdo_work/NewStudy")

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
    left_on='ID_REF',
    right_on='ID',
    how='left'
)

# Clean up the gene symbols BEFORE setting as index
expr_annotated['Gene Symbol'] = expr_annotated['Gene Symbol'].astype(str).str.split(' /// ').str[0]
expr_annotated = expr_annotated[expr_annotated['Gene Symbol'].notna()]
expr_annotated = expr_annotated[expr_annotated['Gene Symbol'] != "---"]

# Now set the cleaned gene symbols as the index
expr_annotated.set_index('Gene Symbol', inplace=True)
expr_annotated.drop(columns=['ID_REF'], inplace=True)

print(expr_annotated.head())

# Extract sample and archetype
sample_info = []

for gsm_name, gsm in gse.gsms.items():
    characteristics = gsm.metadata.get('characteristics_ch1', [])
    
    # Find the cluster info
    cluster = None
    for char in characteristics:
        if "rejection archetype" in char.lower():
            parts = char.split(":", 1)
            if len(parts) > 1:
                cluster = parts[1].strip()  # take the part after colon

    sample_info.append({'Sample': gsm_name, 'Cluster': cluster})

sample_df = pd.DataFrame(sample_info)
sample_df.to_csv(r"C:/Users/nanth/Desktop/Erdo_work/NewStudy/GSE277334_metadata.csv", index=False)

# Ensure the Gene Symbol column exists
if 'Gene Symbol' not in expr_annotated.columns:
    raise KeyError("Column 'Gene Symbol' not found in expr_annotated")

# Set 'Gene Symbol' as index for grouping
expr_annotated = expr_annotated.set_index('Gene Symbol')

# Expression-only matrix (numeric columns)
expr_values = expr_annotated.select_dtypes(include='number')

# Gene-level annotation (non-numeric columns)
gene_annotations = expr_annotated[
    [
        'Gene Ontology Biological Process',
        'Gene Ontology Cellular Component',
        'Gene Ontology Molecular Function',
        'Pathway'
    ]
].copy()

# Remove duplicate annotation rows (keep first)
gene_annotations = gene_annotations[~gene_annotations.index.duplicated(keep='first')]

# Save both separately
expr_values.to_csv(r"C:/Users/nanth/Desktop/Erdo_work/NewStudy/GSE277334_expression_uncollapsed.csv")
gene_annotations.to_csv(r"C:/Users/nanth/Desktop/Erdo_work/NewStudy/GSE277334_gene_annotations.csv")

# Collapse expression values by taking the median per gene symbol
expr_collapsed = expr_values.groupby(expr_annotated.index).median()

expr_collapsed.to_csv(r"C:/Users/nanth/Desktop/Erdo_work/NewStudy/GSE277334_expression_collapsed.csv")


print(f"Shape before: {expr_values.shape}, after collapsing: {expr_collapsed.shape}")

# Quality control of RAW data

# Perform PCA on transposed data (samples as rows)
pca = PCA()
pca_raw = pca.fit_transform(expr_collapsed.T)

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
plt.savefig(r"C:/Users/nanth/Desktop/Erdo_work/NewStudy/GSE277334rawPCA.png", dpi=300)
plt.show()

## Box plot of log2 intensities RAW data

# Compute the median expression for each gene
gene_medians = expr_collapsed.median(axis=1)

# Compute RLE (subtract gene median from each sample)
rle_data = expr_collapsed.subtract(gene_medians, axis=0)

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
plt.savefig(r"C:/Users/nanth/Desktop/Erdo_work/NewStudy/GSE277334_BoxIntensities.png", dpi=300)
plt.show()


## Filtering on Intensity

# Compute row medians correctly
row_medians = expr_collapsed.median(axis=1).astype(float)  # <- force float

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
plt.savefig(r"C:/Users/nanth/Desktop/Erdo_work/NewStudy/GSE277334_HistoIntensities.png", dpi=300)
plt.show()

# Filtering by intensity
# 1. Count number of samples in each cluster
no_of_samples = sample_df['Cluster'].value_counts()
print(no_of_samples)

# 2. Minimum number of samples across clusters (samples_cutoff)
samples_cutoff = no_of_samples.min()
print(f"Samples cutoff: {samples_cutoff}")

# 3. Filter genes: must have expression > man_threshold in at least samples_cutoff samples
man_threshold = 3.5
idx_man_threshold = (expr_collapsed > man_threshold).sum(axis=1) >= samples_cutoff
print(idx_man_threshold.value_counts())  # True = genes kept, False = filtered out

# 4. Subset filtered genes
chang_manfiltered = expr_collapsed[idx_man_threshold]
print(f"Number of genes after filtering: {chang_manfiltered.shape[0]}")

chang_manfiltered.to_csv(r"C:/Users/nanth/Desktop/Erdo_work/NewStudy/GSE277334_expression_collapsed_filtered.csv")
 