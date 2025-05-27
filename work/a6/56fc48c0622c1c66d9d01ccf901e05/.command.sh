#!/bin/bash -ue
mkdir -p organized
python /workspace/scripts/organize_samples_scRNA.py --input GSE263995_data --output organized
