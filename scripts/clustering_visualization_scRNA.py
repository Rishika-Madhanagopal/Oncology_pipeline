import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load annotated dataset
adata_combined = sc.read('data/scRNAseq/normalized_annotated.h5ad')

#Set figure parameters
sc.settings.set_figure_params(dpi=300, dpi_save=300)

#clustering
sc.tl.umap(adata_combined)
sc.tl.leiden(adata_combined)

sc.pl.umap(adata_combined, color='leiden', legend_loc=None, show=True, save='leiden.png')
sc.pl.umap(adata_combined, color='sample_id', legend_loc=None, show=True, save='sample_id.png')
#annotate the cluster
sc.tl.rank_genes_groups(adata_combined, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata_combined, n_genes=10, sharey=False, save='_markers.png')


# Number of top genes to extract per cluster
n_top = 3

# Get ranked gene names for each cluster
ranked_genes = adata_combined.uns['rank_genes_groups']['names']

# Collect top N genes per cluster, ensuring uniqueness
top_genes_set = set()
for cluster in adata_combined.obs['leiden'].cat.categories:
    genes = ranked_genes[cluster][:n_top]
    top_genes_set.update(genes)

# Convert to sorted list
top_unique_genes = sorted(top_genes_set)
print(f"Top unique marker genes: {top_unique_genes}")
#Inspect marker gene expression by cluster
sc.pl.dotplot(
    adata_combined,
    var_names = top_unique_genes,  # your marker genes
    groupby = 'leiden',
    standard_scale = 'var',
    dot_max = 0.7,
    dendrogram = False,
    figsize = (18, 10),
    show = True,
    save = 'dot_plot.png'
)


# Get top N genes per cluster
n_top = 5

# Create a dictionary of top genes for each cluster
ranked = adata_combined.uns['rank_genes_groups']
clusters = ranked['names'].dtype.names

# Build a DataFrame
top_genes_df = pd.DataFrame({
    cluster: ranked['names'][cluster][:n_top]
    for cluster in clusters
})

# Transpose for better readability
top_genes_df = top_genes_df.T
top_genes_df.columns = [f"Top{i+1}" for i in range(n_top)]

# Display top genes
top_genes_df.head()

#Annotate These Clusters in Your Mapping
cluster_to_celltype = {
    '0': 'CD4+ T cells',
    '1': 'B cells',
    '2': 'NK cells',
    '3': 'Regulatory/Activated T cells',
    '4': 'CD8+ T cells'
}
adata_combined.obs['cell_type'] = adata_combined.obs['leiden'].map(cluster_to_celltype)
cluster_to_celltype.update({
    '5': 'Monocytes',
    '6': 'Dendritic cells',
    '7': 'T cells',
    '8': 'Plasma cells',
    # continue for remaining clusters
})
# Step 1: Add 'Unknown' to the list of allowed categories
adata_combined.obs['cell_type'] = adata_combined.obs['cell_type'].astype('category')
adata_combined.obs['cell_type'] = adata_combined.obs['cell_type'].cat.add_categories(['Unknown'])

# Step 2: Fill the missing values
adata_combined.obs['cell_type'] = adata_combined.obs['cell_type'].fillna('Unknown')
sc.settings.set_figure_params(dpi=300, dpi_save=300)
sc.pl.umap(adata_combined,color = 'cell_type',legend_loc = 'on data',legend_fontsize = 2,save = 'cell_type_withnames.png')


#Save the final annotated object
adata_combined.write('data/scRNAseq/normalized_annotated.h5ad')

