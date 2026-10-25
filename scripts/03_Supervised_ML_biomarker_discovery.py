# Supervised Machine Learning (Lasso Regression, Xgboost, RF)
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import accuracy_score, confusion_matrix, ConfusionMatrixDisplay
from xgboost import XGBClassifier
import shap
import tqdm as notebook_tqdm

import warnings
warnings.filterwarnings('ignore')

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier

# C=10 kar do to zyada genes select honge, 2 ki jagah 8-10 milenge
lasso = LogisticRegression(penalty='l1', solver='liblinear', C=10.0, random_state=42)
lasso.fit(X, y)

lasso_coefs = pd.DataFrame({'Gene': X.columns, 'Coefficient': lasso.coef_[0]})
lasso_selected = lasso_coefs[lasso_coefs['Coefficient'] != 0].copy()

# Agar abhi bhi 2 hi aaye to ye line use karo - top 10 force select
if len(lasso_selected) < 5:
    # RF se top 10 nikal lo
    temp_rf = RandomForestClassifier(n_estimators=100, random_state=42).fit(X, y)
    importances = pd.Series(temp_rf.feature_importances_, index=X.columns).sort_values(ascending=False)
    selected_genes = importances.head(10).index.tolist()
else:
    selected_genes = lasso_selected['Gene'].tolist()

X_filtered = X[selected_genes]
print(f"Final selected genes for ML: {len(selected_genes)}")

# XGBoost me se use_label_encoder hata diya
rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
xgb_model = XGBClassifier(eval_metric='logloss', random_state=42, verbosity=0)
loo = LeaveOneOut()

# --- For RF ---
y_true_rf, y_pred_rf = [], []
for train_idx, test_idx in loo.split(X_filtered):
    rf_model.fit(X_filtered.iloc[train_idx], y[train_idx])
    y_true_rf.append(y[test_idx][0])
    y_pred_rf.append(rf_model.predict(X_filtered.iloc[test_idx])[0])

print("RF LOOCV Accuracy:", accuracy_score(y_true_rf, y_pred_rf))
cm = confusion_matrix(y_true_rf, y_pred_rf)
ConfusionMatrixDisplay(cm, display_labels=['Control','ACS']).plot(cmap='Blues', values_format='d')
plt.title('RF - LOOCV - ACS')
plt.show()

# --- For XGBoost ---
y_true_xgb, y_pred_xgb = [], []
for train_idx, test_idx in loo.split(X_filtered):
    xgb_model.fit(X_filtered.iloc[train_idx], y[train_idx])
    y_true_xgb.append(y[test_idx][0])
    y_pred_xgb.append(xgb_model.predict(X_filtered.iloc[test_idx])[0])

print("XGBoost LOOCV Accuracy:", accuracy_score(y_true_xgb, y_pred_xgb))
cm2 = confusion_matrix(y_true_xgb, y_pred_xgb)
ConfusionMatrixDisplay(cm2, display_labels=['Control','ACS']).plot(cmap='Oranges', values_format='d')
plt.title('XGBoost - LOOCV - ACS')
plt.show()

# ---------------------------------------------------------
# STEP 4: SHAP Analysis for XGBoost
# ---------------------------------------------------------
print("\n--- 4. SHAP Interpretability ---")
# TreeExplainer XGBoost ke liye optimized hota hai
explainer = shap.TreeExplainer(xgb_model)
shap_values = explainer.shap_values(X_filtered)

# SHAP Summary Plot
shap.summary_plot(shap_values, X_filtered, feature_names=X_filtered.columns, show=False)
plt.title("SHAP Summary Plot - Top PCOS (COC) Biomarkers")
plt.show()

# Downstream pathway analysis ke liye dataframe (rf_df) create kar rahe hain
# taaki aapka aage ka MyGeneInfo ka code bina error ke chal sake.
rf_df = pd.DataFrame({
    'Gene': X_filtered.columns,
    'Importance': np.abs(shap_values).mean(axis=0) 
}).sort_values(by='Importance', ascending=False)

print("\nTop Biomarkers based on SHAP:")
print(rf_df.head(10))


# Final biomarkers saving- for R heatmap plot generation
import pandas as pd

# Selected genes given by lasso
final_df = pd.DataFrame({'gene': selected_genes})

# Symbol adding for genes names
if 'SYMBOL' in deseq_df.columns:
    final_df['symbol'] = final_df['gene'].map(deseq_df['SYMBOL'])
elif 'gene_name' in deseq_df.columns:
    final_df['symbol'] = final_df['gene'].map(deseq_df['gene_name'])

# Saveing file
final_df.to_csv("final_biomarkers.csv", index=False)
print("Saved final_biomarkers.csv")
print(final_df)


import pandas as pd

# 1. File load
df_vst = pd.read_csv("vst_counts_for_ML.csv")
df_bio = pd.read_csv("final_biomarkers.csv")
df_deseq = pd.read_csv("DESeq2_full_results.csv") # for mapping

print("vst ka shape:", df_vst.shape)
print("vst ke first 5 columns:", list(df_vst.columns)[:5])
print("vst ke first column ke first 5 values:", df_vst.iloc[:5,0].tolist())

# 2. Transposing file if file is reverse
first_col = df_vst.columns[0]
if df_vst[first_col].astype(str).str.contains("ENSG").any():
    print("File ulta hai, transpose kar raha hu...")
    df_vst = df_vst.set_index(first_col).T
    df_vst.index.name = "sample"
    df_vst.reset_index(inplace=True)

# 3. ENSG version is removing 
def clean_ensg(x):
    return str(x).split(".")[0]

df_vst.columns = [clean_ensg(c) for c in df_vst.columns]
bio_col = 'gene' if 'gene' in df_bio.columns else df_bio.columns[0]
biomarkers_raw = [clean_ensg(g) for g in df_bio[bio_col].tolist()]

# 4. Now filtering
biomarkers = [g for g in biomarkers_raw if g in df_vst.columns]
print("Now Biomarkers are matched:", biomarkers)
print("How many:", len(biomarkers), "/ 6")

# If stills get 0 then trying with DESeq2 mapping
if len(biomarkers) == 0:
    print("ENSG match nahi hua, symbol se try kar raha hu...")
    
    # DESeq2 me ENSG -> symbol mapping
    mapping = dict(zip(df_deseq['gene'].apply(clean_ensg), df_deseq['symbol']))
    
    # Changing vst columns from ENSG to genes names
    df_vst.rename(columns=mapping, inplace=True)
    biomarkers = [mapping.get(g,g) for g in biomarkers_raw if mapping.get(g,g) in df_vst.columns]
    print("Symbol se mile:", biomarkers)

# 5.X is creating
X = df_vst[biomarkers].values
print("X shape:", X.shape)


from sklearn.model_selection import LeaveOneOut
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_curve, auc, accuracy_score
import matplotlib.pyplot as plt

if "condition" in df_vst.columns:
    y_str = df_vst["condition"].values
elif "Condition" in df_vst.columns:
    y_str = df_vst["Condition"].values
elif "group" in df_vst.columns:
    y_str = df_vst["group"].values
else:
    # if the condition does not exists then giving manually
    print("Condition column not founded , making new columns - checking!")
    print(df_vst['sample'].tolist() if 'sample' in df_vst.columns else "sample column bhi nahi hai")
    y_str = ["Treated","Treated","Treated","Control","Control","Control"]

print("y_str:", y_str)

# 2. y_true level - 0/1
le = LabelEncoder()
y = le.fit_transform(y_str)  # Control=0, Treated=1
print("Label mapping:", dict(zip(le.classes_, le.transform(le.classes_))))
print("y (0/1):", y)

# 3. LOOCV + proba level
loo = LeaveOneOut()
y_true = []
y_proba = []
y_pred = []

model = RandomForestClassifier(n_estimators=100, random_state=42)

for train_idx, test_idx in loo.split(X):
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    
    model.fit(X_train, y_train)
    
    y_true.append(y_test[0])  # true level
    y_pred.append(model.predict(X_test)[0])
    y_proba.append(model.predict_proba(X_test)[:,1][0])  # proba's level of class 1

print("y_true:", y_true)
print("y_proba:", y_proba)
print("Accuracy:", accuracy_score(y_true, y_pred))

# 4. ROC / AUC
fpr, tpr, _ = roc_curve(y_true, y_proba)
roc_auc = auc(fpr, tpr)

plt.figure()
plt.plot(fpr, tpr, label=f'RF ROC (AUC = {roc_auc:.2f})')
plt.plot([0,1], [0,1], 'k--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('CAD (COC) - Final 6 Biomarkers ROC')
plt.legend()
plt.savefig("ROC_AUC_CAD (COC).png", dpi=300)
plt.show()

print(f"FINAL AUC: {roc_auc}")