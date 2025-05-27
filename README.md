# 🧬 Oncology Pipeline: Multi-Omics Drug Sensitivity Analysis

This repository contains a complete bioinformatics pipeline for analyzing **scRNA-seq**, **mutation**, **CNV**, and **drug response data** in cancer, specifically breast cancer. The goal is to integrate these data modalities and build predictive models of drug response, followed by interpretability and network analysis.

---

## 📁 Project Structure

```bash
Oncology_pipeline/
├── data/                      # Raw and processed datasets (ignored in repo)
├── scripts/                   # Python and R scripts for each phase
├── models/                    # Saved ML models
├── Result/                    # Figures, SHAP plots, enrichment results
├── Network_analysis/          # Cytoscape files
├── Dockerfile                 # Reproducible environment container
├── environment.yml            # Conda environment specification
├── nextflow.config            # Nextflow configuration
├── main.nf                    # Main Nextflow pipeline
├── .gitignore                 # Large files ignored from Git
└── README.md
```


---

## 🔄 Workflow Overview

### ✅ Phase 1: Data Acquisition
- **scRNA-seq**: Downloaded from GEO (10x Genomics format)
- **Mutations/CNV**: Retrieved from TCGA via GDC API
- **Drug Sensitivity**: Extracted from DepMap (GDSC format)

### ✅ Phase 2: Preprocessing
- `scripts/organize_samples_scRNA.py`: Organize scRNA-seq folders
- `scripts/preprocess_scrna.py`: Quality control, normalization, PCA, clustering using Scanpy
- `scripts/GDC_mutation.py`, `GDC_cnv.py`, `demap.py`: Format mutation, CNV, and drug datasets

### ✅ Phase 3: Integration & Feature Engineering
- Align samples across datasets
- Normalize features, scale expression
- Output: `data/master_features_scaled.csv`

### ✅ Phase 4: Machine Learning
- Train **RandomForestClassifier** on selected drug targets
- Evaluate using cross-validation (AUC, confusion matrix)
- Save model: `models/model_rf_<drug>.pkl`

### ✅ Phase 5: Model Interpretation
- Use **SHAP** to generate global and per-sample explanations
- Save:
  - SHAP summary plot (`Result/shap_summary_bar.png`)
  - Force plot for local explanation
  - Top features: `data/top_predictive_features.csv`

### ✅ Phase 6: Network Analysis
- Functional enrichment with **gProfiler2**
- Visualize pathways and networks in Cytoscape
- Save: enrichment plot + Cytoscape session

---

## ⚙️ Reproducibility

This project is **fully reproducible** using:

- 📦 **Docker**: Environment encapsulated with all dependencies  
- 🧬 **Nextflow**: Pipeline management system for reproducible workflows

### To Run Locally:
```bash
# Install Nextflow and Docker
nextflow run main.nf -profile docker
