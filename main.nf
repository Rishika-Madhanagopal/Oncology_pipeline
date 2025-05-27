workflow {
  organize_scRNA()
  preprocess_scRNA()
  variant_cnv_processing()
  merge_and_train()
  interpret_and_enrich()
}

process organize_scRNA {
  input:
  path 'scripts/organize_samples.py'
  output:
  path "data/raw/*"
  script:
  """
  python scripts/organize_samples.py
  """
}

process preprocess_scRNA {
  input:
  path 'scripts/preprocess_scrna.py'
  output:
  path "data/scRNAseq/normalized_raw.h5ad"
  script:
  """
  python scripts/preprocess_scrna.py
  """
}

process variant_cnv_processing {
  script:
  """
  python scripts/GDC_mutation.py
  python scripts/GDC_cnv.py
  """
}

process merge_and_train {
  script:
  """
  python scripts/feature_engineering.py
  python scripts/merging_drug.py
  python scripts/model_training.py
  """
}

process interpret_and_enrich {
  script:
  """
  python scripts/model_interpretation.py
  Rscript scripts/enrichment_analysis.R
  """
}
