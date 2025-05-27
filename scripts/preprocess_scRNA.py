import scanpy as sc
import scipy.io
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

base_path = r'C:\Users\rishi\OneDrive\Desktop\Project 1\data\GSE263995_data'
sample_dirs = [os.path.join(base_path, d) for d in os.listdir(base_path) if d.startswith('GSM')]

adatas = []

for sd in sample_dirs:
    print(f"Loading: {sd}")
    
    matrix_file = os.path.join(sd, 'matrix.mtx.gz')
    genes_file = os.path.join(sd, 'features.tsv.gz')  # or 'genes.tsv.gz'
    barcodes_file = os.path.join(sd, 'barcodes.tsv.gz')
    
    # Load matrix
    matrix = scipy.io.mmread(matrix_file).T.tocsc()
    
    # Load genes (first two columns only)
    genes = pd.read_csv(genes_file, header=None, sep='\t')
    if genes.shape[1] < 2:
        raise ValueError(f"Gene file in {sd} has fewer than 2 columns!")
    
    var_names = genes[1].values
    var = pd.DataFrame(index=var_names)
    
    # Load barcodes
    barcodes = pd.read_csv(barcodes_file, header=None)[0].values
    obs = pd.DataFrame(index=barcodes)
    
    # Build AnnData
    adata = sc.AnnData(X=matrix, var=var, obs=obs)
    adata.var_names_make_unique()
    
    adata.obs['sample_id'] = os.path.basename(sd)
    adatas.append(adata)

# Combine
batch_categories = [os.path.basename(d) for d in sample_dirs]
adata_combined = sc.concat(adatas, label='batch', keys=batch_categories)
adata_combined.obs['sample_id'] = adata_combined.obs['batch']
adata_combined.var_names_make_unique()

# Save to file for later steps
save_path = r'data/scRNAseq/normalized_raw.h5ad'
os.makedirs(os.path.dirname(save_path), exist_ok=True)
adata_combined.write(save_path)

print(f"Combined AnnData saved to {save_path}")

#print("Unique sample IDs in AnnData:")
#print(adata_combined.obs['sample_id'].unique().tolist())


# Load combined data
#adata_combined = sc.read('data/scRNAseq/normalized_raw.h5ad')
# Label mitochondrial genes
adata_combined.var['mt'] = adata_combined.var_names.str.startswith('MT-')

# Compute QC metrics: n_genes_by_counts, total_counts, pct_counts_mt
sc.pp.calculate_qc_metrics(adata_combined, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter using the correct variable name
adata_combined = adata_combined[adata_combined.obs['n_genes_by_counts'] > 200, :]
adata_combined = adata_combined[adata_combined.obs['n_genes_by_counts'] < 6000, :]
adata_combined = adata_combined[adata_combined.obs['pct_counts_mt'] < 10, :]

#normalize
sc.pp.normalize_total(adata_combined, target_sum=1e4)
sc.pp.log1p(adata_combined)
#Find highly variable genes
sc.pp.highly_variable_genes(adata_combined, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_combined.var['highly_variable'].sum()  #number of genes marked as highly variable

#dimentionality reduction
adata_combined = adata_combined[:, adata_combined.var.highly_variable]
sc.pp.pca(adata_combined, n_comps=50)
sc.pl.pca_variance_ratio(adata_combined, log=True)
# Find neighbors and cluster:
sc.pp.neighbors(adata_combined, n_neighbors=10, n_pcs=40)

# Convert to DataFrame
df = pd.DataFrame(adata_combined.X.toarray(), columns=adata_combined.var_names)
df['sample_id'] = adata_combined.obs['sample_id'].values

# Map GSM to TCGA
df['sample_id'] = df['sample_id'].map(folder_to_tcga)
df = df[~df['sample_id'].isna()]

# Aggregate to one row per TCGA sample
expr_agg = df.groupby('sample_id').mean()
expr_agg.to_csv('data/expression_matrix_aligned.csv')
