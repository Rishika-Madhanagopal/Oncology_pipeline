import pandas as pd
import gzip
import glob
import os

# Set your base folder
base_path = r'C:\Users\rishi\OneDrive\Desktop\Project 1\data\TGCA_data\gdc_mutation'

# Find all .maf.gz files
maf_files = glob.glob(os.path.join(base_path, '**', '*.maf.gz'), recursive=True)

# Combine all into one DataFrame
dfs = []
for file in maf_files:
    print(f"Loading: {file}")
    with gzip.open(file, 'rt') as f:
        df = pd.read_csv(f, sep='\t', comment='#', low_memory=False)
        df['source_file'] = os.path.basename(file)
        dfs.append(df)

maf_combined = pd.concat(dfs, ignore_index=True)

#print(f"Combined MAF shape: {maf_combined.shape}")
#print(maf_combined.head())
maf_combined.to_csv('data/maf_combined.csv', index=False)

# Preprocess
maf_combined = pd.read_csv('data/maf_combined.csv')
mut_df = maf_combined[['Tumor_Sample_Barcode', 'Hugo_Symbol']]
mut_binary = pd.crosstab(mut_df['Tumor_Sample_Barcode'], mut_df['Hugo_Symbol'])
mut_binary = (mut_binary > 0).astype(int)
mut_binary.to_csv('data/annotated_variants.tsv', sep='\t')

mutation_matrix = pd.read_csv('data/annotated_variants.tsv', index_col=0, sep='\t')
mutation_matrix.index = mutation_matrix.index.str.extract(r'(TCGA-\w\w-\w\w\w\w)')[0]
mutation_matrix = mutation_matrix[~mutation_matrix.index.isna()]
print("Mutation sample IDs:", mutation_matrix.index.unique().tolist())
