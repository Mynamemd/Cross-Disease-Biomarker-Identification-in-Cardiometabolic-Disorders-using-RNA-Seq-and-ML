import pandas as pd
import numpy as np 
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Loading Normalized Counts Matrix &  DESeq2 Results
counts_df = pd.read_csv("vst_counts_for_ML.csv", index_col=0)
deseq_df = pd.read_csv("DESeq2_full_results.csv", index_col=0)

print("Vst counts Matrix Shape:", counts_df.shape)
print("DESeq2 Matrix Shape:", deseq_df.shape)

# 2. Filter Top diffrentially expressed genes (DEGs)
deg_mask = (deseq_df['padj'] < 0.07) & (deseq_df['log2FoldChange'].abs() >= 1.0)
top_genes = deseq_df[deg_mask].index.tolist()

print(f"[+] Total raw genes: {counts_df.shape[0]}")
print(f"[+] Total DESeq2 genes: {deseq_df.shape[0]}")
print(f"[+] Seleceted DEGs for ML features: {len(top_genes)}")

# Non-Numeric metadata columns (chr, start, End, Strand, length) is removing
metadata_cols = ['Chr', 'Start', 'End', 'Strand', 'Length']
clean_counts = counts_df.drop(columns=[col for col in metadata_cols if col in counts_df.columns])
 
filtered_counts = clean_counts.loc[clean_counts.index.isin(top_genes)]

# Transposing (sample as row and Genes as columns) due to scikit-learn demand we transpose
X = filtered_counts.T # Shape: (n_samples, n_features)

# Strikly float numeric data type enforce
X = X.apply(pd.to_numeric, errors='coerce').fillna(0).astype(float)

# 4. Target labels (y) making (Strictly for 6 samples)
# for each samples create 1D array for supervised. If CAD found in sample(disease) gives 1 or not gives 0 (Control)
sample_names = X.index.tolist()
treated = ["SRR27462365", "SRR27462366", "SRR27462367"]
y = np.array([1 if any(t in str(s) for t in treated) else 0 for s in sample_names])

print(f"[+] Samples count (X rows): {X.shape[0]}")
print(f"[+] Genes count (X cols): {X.shape[1]}")
print(f"[+] Target labels count (y): {len(y)}")

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# PCA plot Dimensionality Reduction (Scikit-learn uses #2)
pca = PCA(n_components=2) 
pca_coords = pca.fit_transform(X_scaled) 
var_exp = pca.explained_variance_ratio_ * 100 
print(f"[+] PCA Explained Variance: PC1 = {var_exp[0]:.2f}%, PC2 = {var_exp[1]:.2f}%")

y = np.array([0, 0, 0, 1, 1, 1])

# 6. Plotting
pca_df = pd.DataFrame(pca_coords, columns=['PC1', 'PC2'])
pca_df['Condition'] = ['ACS' if label == 1 else 'Control' for label in y]

plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=pca_df, x='PC1', y='PC2', 
    hue='Condition', style='Condition', 
    s=150, palette=['#2b5c8f', '#d95f02']
)
plt.title(f"PCA of Top DEGs (PC1: {var_exp[0]:.2f}%, PC2: {var_exp[1]:.2f}%)", fontsize=12, fontweight='bold')
plt.xlabel(f"Principal Component 1 ({var_exp[0]:.2f}% Variance)")
plt.ylabel(f"Principal Component 2 ({var_exp[1]:.2f}% Variance)")
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("PCA_Sample_Clustering.png", dpi=300)
plt.show()
