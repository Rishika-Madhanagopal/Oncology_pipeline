import pandas as pd
import glob
import os
import re

# Path to SEG files
cnv_path = r'C:\Users\rishi\OneDrive\Desktop\Project 1\data\TGCA_data\gdc_cnv'
cnv_files = glob.glob(os.path.join(cnv_path, '**', '*.seg.txt'), recursive=True)

dfs = []
uuid_to_tcga = {}

for file in cnv_files:
    folder_uuid = os.path.basename(os.path.dirname(file)).lower()
    tcga_match = re.search(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})', file)
    if tcga_match:
        tcga_id = tcga_match.group(1)
        uuid_to_tcga[folder_uuid] = tcga_id

    try:
        with open(file) as f:
            lines = f.readlines()
        header_index = next(i for i, line in enumerate(lines) if 'Segment_Mean' in line)

        df = pd.read_csv(file, sep='\t', skiprows=header_index)
        df['GDC_Aliquot'] = folder_uuid  # Use folder UUID as index key
        dfs.append(df)
        print(f"Loaded: {file}")

    except Exception as e:
        print(f"Failed to read {file}: {e}")

# Combine
cnv_combined = pd.concat(dfs, ignore_index=True)
cnv_combined['copy_number'] = pd.to_numeric(cnv_combined['Segment_Mean'], errors='coerce')
cnv_combined['cnv_status'] = cnv_combined['copy_number'].apply(
    lambda x: 1 if x > 0.3 else -1 if x < -0.3 else 0
)

# Pivot
cnv_matrix = cnv_combined.pivot_table(
    index='GDC_Aliquot',
    columns='Chromosome',
    values='cnv_status',
    aggfunc='median'
).fillna(0).astype(int)

# Map UUIDs to TCGA
cnv_matrix.index = cnv_matrix.index.str.lower()
uuid_to_tcga = {k.lower(): v for k, v in uuid_to_tcga.items()}
cnv_matrix['sample_id'] = cnv_matrix.index.map(uuid_to_tcga)
cnv_matrix = cnv_matrix[~cnv_matrix['sample_id'].isna()]
cnv_matrix = cnv_matrix.set_index('sample_id')

# Save final matrix
cnv_matrix.to_csv('data/cnv_matrix.csv')
print(f"Final CNV matrix shape: {cnv_matrix.shape}")
print(cnv_matrix.head())
