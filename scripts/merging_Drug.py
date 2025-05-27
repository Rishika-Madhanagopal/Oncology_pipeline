import pandas as pd

# Load features and drug data
features = pd.read_csv("data/master_features_scaled.csv", index_col=0)
drug = pd.read_csv("data/drug_data_breast_only.csv", index_col=0)
drug.index = drug.index.str.upper().str.strip()

# Load manual mapping
mapping = pd.read_csv("data/tcga_to_depmap.csv")
mapping["TCGA_ID"] = mapping["TCGA_ID"].str.upper().str.strip()
mapping["DepMap_ID"] = mapping["DepMap_ID"].str.upper().str.strip()

# Merge mapping into features
features = features.reset_index().rename(columns={"index": "TCGA_ID"})
merged = features.merge(mapping, on="TCGA_ID", how="inner").set_index("DepMap_ID")

# Final join with drug data
final = merged.join(drug, how="inner")

print(f"Final shape (features + drug): {final.shape}")
final.to_csv("data/final_features_with_drugs.csv")

