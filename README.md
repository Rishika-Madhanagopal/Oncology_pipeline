# 🧬 Oncology Pipeline: Multi-Omics Drug Sensitivity Analysis

This repository contains a complete bioinformatics pipeline for analyzing **scRNA-seq**, **mutation**, **CNV**, and **drug response data** in cancer, specifically breast cancer. The goal is to integrate these data modalities and build predictive models of drug response, followed by interpretability and network analysis.

---

## 📁 Project Structure

Oncology_pipeline/
├── data/ # Raw and processed datasets (ignored in repo)
├── scripts/ # Python and R scripts for each phase
├── models/ # Saved ML models
├── Result/ # Figures, SHAP plots, enrichment results
├── Network_analysis/ # Cytoscape files
├── Dockerfile # Reproducible environment container
├── environment.yml # Conda environment specification
├── nextflow.config # Nextflow configuration
├── main.nf # Main Nextflow pipeline
├── .gitignore # Large files ignored from Git
└── README.md
