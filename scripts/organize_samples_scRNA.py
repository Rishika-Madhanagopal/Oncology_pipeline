import os
import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help='Path to base folder')
parser.add_argument('--output', required=True, help='Path to output folder')
args = parser.parse_args()

input_dir = args.input
output_dir = args.output

os.makedirs(output_dir, exist_ok=True)

for root, dirs, files in os.walk(input_dir):
    for file in files:
        if any(key in file for key in ['matrix.mtx.gz', 'barcodes.tsv.gz', 'features.tsv.gz', 'genes.tsv.gz']):
            prefix = os.path.basename(root)
            subfolder = os.path.join(output_dir, prefix)
            os.makedirs(subfolder, exist_ok=True)

            suffix_map = {
                'matrix.mtx.gz': 'matrix.mtx.gz',
                'barcodes.tsv.gz': 'barcodes.tsv.gz',
                'features.tsv.gz': 'features.tsv.gz',
                'genes.tsv.gz': 'features.tsv.gz'
            }

            for key, name in suffix_map.items():
                if file.endswith(key):
                    src = os.path.join(root, file)
                    dst = os.path.join(subfolder, name)
                    shutil.copy2(src, dst)
                    print(f"Copied {src} to {dst}")

print('All files organized successfully.')
