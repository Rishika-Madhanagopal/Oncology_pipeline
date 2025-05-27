import pandas as pd
import os
from sklearn.preprocessing import StandardScaler
import pandas as pd
# Step 1: Load mutation matrix
mutation = pd.read_csv('data/annotated_variants.tsv', sep='\t', index_col=0)
mutation.index = mutation.index.str.extract(r'(TCGA-\w\w-\w\w\w\w)')[0]
mutation = mutation[~mutation.index.isna()]
mutation.index = mutation.index.str.upper().str.strip()

# Step 2: Load CNV matrix
cnv = pd.read_csv('data/cnv_matrix.csv', index_col=0)
cnv.index = cnv.index.str.upper().str.strip()

# Step 3: Load expression matrix
expr = pd.read_csv('data/expression_matrix_aligned.csv', index_col=0)
expr.index = expr.index.str.upper().str.strip()

# Step 4: Find common sample IDs
shared_ids = list(set(mutation.index) & set(cnv.index) & set(expr.index))
print(f"Shared TCGA sample IDs: {len(shared_ids)} → {shared_ids}")

#  Step 5: Subset each matrix
mutation_filtered = mutation.loc[shared_ids]
cnv_filtered = cnv.loc[shared_ids]
expr_filtered = expr.loc[shared_ids]

# Step 6: Concatenate to form the master matrix
master = pd.concat([mutation_filtered, cnv_filtered, expr_filtered], axis=1)
print(f"📐 Master matrix shape: {master.shape}")

# Step 7: Save master matrix
os.makedirs('data', exist_ok=True)
master.to_csv('data/master_features_matrix.csv')
print("Saved: data/master_features_matrix.csv")

#Preprocess (Scaling)
# Load master matrix
df = pd.read_csv('data/master_features_matrix.csv', index_col=0)

# Identify expression columns
expression_cols = df.columns[-1588:]

# Scale expression values
scaler = StandardScaler()
df[expression_cols] = scaler.fit_transform(df[expression_cols])

# Save scaled matrix
df.to_csv('data/master_features_scaled.csv')
print("Saved: data/master_features_scaled.csv (expression scaled only)")
