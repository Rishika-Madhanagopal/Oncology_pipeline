import shap
import joblib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Load model and dataset
model = joblib.load("models/model_rf_BRD-A00077.pkl")
df = pd.read_csv("data/final_features_with_drugs.csv", index_col=0)

# Isolate drug column and numeric features
drug_cols = [col for col in df.columns if col.startswith("BRD-")]
X = df.drop(columns=drug_cols).select_dtypes(include="number")

# SHAP TreeExplainer (works with RandomForest)
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X)

# Validate shape
shap_array = shap_values[1]  # class 1
if shap_array.shape[0] != X.shape[0]:
    shap_array = shap_array.T

print("SHAP shape:", shap_array.shape)
print("Feature matrix shape:", X.shape)

# Create output directory
os.makedirs("Result", exist_ok=True)
os.makedirs("data", exist_ok=True)

# Summary Plot (Global Importance)
plt.figure()
shap.summary_plot(shap_array, X, plot_type="bar", max_display=20, show=False)
plt.tight_layout()
plt.savefig("Result/shap_summary_bar.png", bbox_inches="tight", dpi=300)
plt.close()

# Save top features
mean_abs_shap = np.abs(shap_array).mean(axis=0)
top_df = pd.DataFrame({"gene": X.columns, "mean_abs_shap": mean_abs_shap})
top_df.sort_values("mean_abs_shap", ascending=False).head(20).to_csv(
    "data/top_predictive_features.csv", index=False
)

# Force Plot for single sample (local explanation)
shap.initjs()
i = 0
force_plot = shap.force_plot(
    explainer.expected_value[1],
    shap_array[i],
    X.iloc[i],
    matplotlib=True,
    show=False
)
plt.gcf().set_size_inches(10, 3)
plt.savefig("Result/shap_per_sample.png", bbox_inches="tight", dpi=300)
plt.close()

print("All SHAP plots and top features saved.")
