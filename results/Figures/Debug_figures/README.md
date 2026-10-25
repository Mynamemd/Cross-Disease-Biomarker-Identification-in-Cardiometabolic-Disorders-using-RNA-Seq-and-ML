# Debug - Initial Incorrect Performance Metrics

This folder contains results with inflated/deflated metrics due to methodological issues.

## Files
- confusion_matrix_RF_100percent.png - 100% accuracy
- confusion_matrix_XGBoost_0percent.png - 0% accuracy
- ROC_AUC_1.0.png - AUC 1.0 (100%)

## Why these incorrect results occurred
1. **RF 100% Accuracy:** Model evaluated on training data itself (no train-test split). Data leakage and overfitting due to small sample size (n=6).
2. **XGBoost 0.0%:** Class imbalance and improper hyperparameter tuning with small sample, model failed to learn.
3. **AUC 1.0:** Perfect separation due to leakage - features included label information.

These debug results are retained only for troubleshooting record.
