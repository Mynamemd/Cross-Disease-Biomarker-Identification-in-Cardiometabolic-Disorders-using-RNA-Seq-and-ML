# Pan-Metabolic Biomarker Discovery

Cross-disease biomarker identification in cardiometabolic disorders using RNA-Seq and machine learning.

## Overview
Identification of minimal predictive biomarkers across five cardiometabolic diseases (ACS, Stable Angina, T2D, PCOS, SA). Final panel: 28 stringent biomarkers (padj<0.05, |log2FC|>1), 18 with high-confidence PPI interactions, 1 pan-biomarker BSG (CD147/EMMPRIN).

## Workflow
Raw Counts → DESeq2 → VST → ML (RF/XGBoost/Lasso) → Merged Panel (33) → Filtered Panel (28) → UpSet/DotPlot/PPI(18)

## Methods
- DESeq2: Differential expression with sizeFactors normalization
- VST: Variance Stabilizing Transformation for ML/PCA/Heatmap
- ML: Random Forest, XGBoost, Lasso with LOOCV
- Pan-analysis: UpSet, DotPlot, STRING PPI

## Results
- Total ML biomarkers: 33
- Final stringent: 28
- PPI network: 18/28
- Pan-biomarker: BSG (ENSG00000162366) - ACS+SA shared, log2FC -1.23, padj 0.045

## Repository Structure

Pan_Metabolics_diseases_Biomarkers/
├── README.md
├── .gitignore
├── data/
│   ├── raw/
│   └── processed/vst_counts_for_ML.csv

├── scripts/
│   ├── 01_DESeq2_analysis.R
│   ├── 02Unsupervised_ML_biomarker_discovery.py
│   ├── 03Supervised_final_table_with_padj.py
│   ├── 04_UpSet_DotPlot_28.R
│   └── 05_PPI_18.R
    └── 06_final_table_with_padj.R

├── results/
│   ├── DESeq2_results/
        └── DESeq2_full_results_*

│   ├── merged_results/
│   │   ├── pan_all_biomarkers_33.csv
│   │   ├── Final_28_with_padj.csv
│   │   └── Borderline_5_Biomarkers.csv
        └── common_pan_biomarkers.csv

│   └── figures/
        ├── Final_figures/
│            ├── UpSet_28.png
│            ├── DotPlot_28.png
│            └── PPI_18.png
             └── Boxplot_Biomarkers_*

        ├── Debug_figures/
             ├──README.md
             └── RF_confusion_matrix_100%.png
             └── XGB_confusion_matrix.png
             └── ROC_AUC_PCOS (ACS).png
             └── SHAP_chart.png

        ├── Enrichment_figures/
             └── KEGG_Pathways_* 
 
## Installation
R: DESeq2, org.Hs.eg.db, ggplot2, UpSetR, STRINGdb
Python: scikit-learn, pandas

## Usage
Rscript scripts/01_DESeq2_analysis.R
Rscript scripts/02_ML_biomarker_discovery.R
Rscript scripts/03_final_table_with_padj.R
Rscript scripts/04_UpSet_DotPlot_28.R
Rscript scripts/05_PPI_18.R

## Figures
- UpSet_28: Unique vs shared
- DotPlot_28: Disease-wise distribution
- PPI_18: 3 modules - HLA immune, Mitochondrial, Inflammatory

## Contact
MD Mobarak Ansari
