import pandas as pd
import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
import joblib

# Load data
df = pd.read_csv("data/final_features_with_drugs.csv", index_col=0)

# Identify drug and feature columns
drug_cols = [col for col in df.columns if col.startswith("BRD-")]
target_drug = 'BRD-A00077618-236-07-6::2.5::HTS'  # replace with other drug if needed
feature_cols = [col for col in df.columns if col not in drug_cols]

# Prepare features and target
X = df[feature_cols].select_dtypes(include=np.number).fillna(0)
y = (df[target_drug] > df[target_drug].median()).astype(int)

# Scale features
scaler = StandardScaler()
X_scaled = pd.DataFrame(scaler.fit_transform(X), index=X.index, columns=X.columns)

print("Data shape:", X_scaled.shape)
print("Target drug:", target_drug)
print("Class distribution:\n", y.value_counts())

# Model
model = RandomForestClassifier(n_estimators=100, random_state=42)

# Cross-validation (2-fold) since having less no.of.sample
cv = StratifiedKFold(n_splits=2, shuffle=True, random_state=42)
y_pred = cross_val_predict(model, X_scaled, y, cv=cv)
y_proba = cross_val_predict(model, X_scaled, y, cv=cv, method='predict_proba')[:, 1]

# Evaluation
print("\n Classification Report:")
print(classification_report(y, y_pred))

print("Confusion Matrix:")
print(confusion_matrix(y, y_pred))

if len(np.unique(y)) > 1:
    auc = roc_auc_score(y, y_proba)
    print(f"ROC AUC (2-fold CV): {auc:.3f}")
else:
    print("ROC AUC not defined (only one class present).")

# Final model (train on all for deployment)
model.fit(X_scaled, y)

# Save model
os.makedirs("models", exist_ok=True)
model_path = f"models/model_rf_{target_drug[:10]}.pkl"
joblib.dump(model, model_path)
print(f"Model saved to: {model_path}")


# For large data sample use split/train/test/validation
'''
## === Train/test split ===
X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y, test_size=0.2, stratify=y, random_state=42
)

# === Train model ===
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# === Predict and evaluate ===
y_pred = model.predict(X_test)
y_proba = model.predict_proba(X_test)[:, 1]
'''