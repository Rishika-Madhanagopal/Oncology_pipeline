import pandas as pd

# Load data
base_path = r'C:\Users\rishi\OneDrive\Desktop\Project 1\data\DepMap_data'
sample_file = f"{base_path}\\sample_info.csv"
drug_file = f"{base_path}\\primary-screen-replicate-collapsed-logfold-change.csv"

sample_info = pd.read_csv(sample_file)
drug_data = pd.read_csv(drug_file)

# Step 1: Filter breast cancer DepMap IDs
breast_ids = sample_info[sample_info['lineage'].str.lower() == 'breast']['DepMap_ID']
breast_ids = breast_ids.str.replace(r'_FAILED.*', '', regex=True).tolist()

# Step 2: Filter rows of drug_data based on breast_ids
drug_data_filtered = drug_data[drug_data['Unnamed: 0'].isin(breast_ids)]

# Step 3: Rename first column to 'DepMap_ID'
drug_data_filtered = drug_data_filtered.rename(columns={'Unnamed: 0': 'DepMap_ID'})

# Preview
print(f"Final filtered shape: {drug_data_filtered.shape}")
print(drug_data_filtered.head())
drug_data_filtered.to_csv(f"data\drug_data_breast_only.csv", index=False)

