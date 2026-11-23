

Heart Failure: 20-Gene Transcriptomic Signature

An end-to-end, pipeline that builds and validates a 20-gene blood/tissue transcriptomic signature for heart failure (HF) from public GEO microarray data. It covers data collection, preprocessing, differential expression, LASSO gene selection, comparison of 12 machine-learning and deep-learning models, statistical validation, external validation, explainability, immune infiltration, enrichment analysis and WGCNA.


 1. Datasets (all public, from NCBI GEO)

Role - GEO accession - Use in this notebook 
Discovery / training + internal test - GSE57338 - Merged with GSE5406, then split 80/20.
Discovery / training + internal test - GSE5406  - Merged with GSE57338, then split 80/20. 
External validation - GSE42955 - Never used for training or gene selection.

Data are downloaded automatically with `GEOparse`, so internet must be ON.


2. Pipeline overview

1. Download & extraction: probe-to-gene mapping from each platform annotation; probes mapped to the same gene are averaged.
2. Cleaning: genes with more than 20% missing values are dropped; remaining missing values are median-imputed; values are log2(x+1) transformed.
3. Merging & normalisation: GSE57338 and GSE5406 are merged on common genes; low-variance genes (bottom 10th percentile) are removed; per-gene z-score normalisation is applied.
   Note: `batch_correct()` supports ComBat when batch labels are passed. In the notebook as written it is called without batch labels, so the z-score path is what runs.
4. Fixed train/test split: 80/20, stratified, `random_state=42`. The same split is reused in every later step (no re-splitting).
5. Leakage-safe feature selection (training set only)**:
   - Differential expression: Welch t-test, Benjamini-Hochberg FDR < 0.05 and |log2FC| > 0.5
   - L1-penalised logistic regression (LASSO, `C=0.05`) to select the 20-gene panel
6. Model comparison (12 models)
   - Classical ML (8): Logistic Regression, Random Forest, SVM, XGBoost, KNN, Naive Bayes, AdaBoost, Gradient Boosting
   - Deep learning (4): MLP, 1D-CNN, ResMLP, FT-Transformer
   - 5-fold stratified CV on the training set; the held-out test set is used for final evaluation
7. Statistical validation
   - Bootstrap 95% CI for AUC (1000 resamples)
   - Friedman test + Nemenyi post-hoc across models
   - Permutation test (200 permutations for classical models, 100 for deep-learning models)
   - LASSO selection stability (100 bootstrap resamples)
   - Learning curves, Brier score, calibration, decision-curve analysis (DCA), threshold tuning
   - Per-dataset test AUC (GSE57338 vs GSE5406)
8. External validation on GSE42955 (several models, with gene coverage reported)
9. Explainability: SHAP for SVM, Naive Bayes and Random Forest, with rank-agreement between models
10. Biology: immune-infiltration scoring (marker-gene scores, Mann-Whitney with FDR), GO-BP / KEGG enrichment (`gseapy`), WGCNA with hub-gene validation
11. Packaging: results and processed data are zipped for download


3. Requirements

### Platform
- Kaggle Notebook
- Accelerator: GPU T4 x2
- Internet: ON
- Output folder used by the notebook: `/kaggle/working`

### Python packages
The first cell installs:

```
GEOparse  gseapy  matplotlib-venn  statsmodels  shap  xgboost
```

Other libraries used (normally already on Kaggle): `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, `seaborn`, `torch`.
WGCNA step auto-installs `dynamicTreeCut` and `fastcluster` if missing; `combat` is only installed if batch labels are supplied.

---

## 4. How to run

1. Upload the notebook to Kaggle (or open it there).
2. In **Settings**, turn **Internet ON** and set **Accelerator = GPU T4 x2**.
3. Run cells **from top to bottom** (`Run All`). Do not skip cells, because later cells reuse variables created earlier (`hf_models`, `hf_genes`, `heart_train_cols`, and so on).
4. When finished, download:
   - `PROCESSED_DATA.zip` (processed matrices, labels, gene list, CSV tables)
   - `HEART_ALL_RESULTS.zip` (all `.svg`, `.png`, `.csv`, `.pdf` outputs)

Expected run time

Setting ( Total time )
Kaggle GPU T4 x2, full notebook ( about 2 h 30 min to 3 h )

The slowest part is the permutation test , which retrains every model many times, deep-learning models especially. It saves progress to `permutation_checkpoint.pkl` and prints `Resuming: finished models = [...]` if you re-run it after an interruption, so already-finished models are not repeated.

---

5. Reproducibility

- Fixed seed `42` for the data split, CV folds, model initialisation and bootstrap resampling.
- A single fixed train/test split is reused throughout.
- DEG analysis, LASSO selection and feature scaling are fitted on the **training set only**.
- The external cohort (GSE42955) is used only for evaluation.

Small numerical differences may still occur between runs because of GPU non-determinism and library versions.


6. Main output files

Category (Files)

Processed data (`heart_corrected.pkl`, `heart_labels.npy`, `hf_genes.npy`)
DEG & gene panel (`heart_DEG_results.csv`, `LASSO_stability_HF.csv`, `Heatmap_20genes_HF.*`, `Panel_overview_HF.*`)
Model performance (`HF_12models_master_table.csv`, `CV_AUC_summary_12models.csv`, `bootstrap_CI_12models.csv`, `ROC_curves_HF_all12.svg`, `PR_AUC_curves_HF_12models.svg`, `Confusion_Matrix_HF.svg`)
Statistics (`statistical_validation_12models.csv`, `Friedman_mean_ranks_12models.csv`, `Nemenyi_*`, `Permutation_test_HF.*`, `Learning_curve_HF*`)
Calibration & utility (`Brier_scores_HF_12models.csv`, `Calibration_HF.svg`, `DCA_HF.svg`, `DCA_HF_values.csv`, `Threshold_Tuning_HF.svg`)
External validation (`external_validation_multimodel.csv`, `external_validation_12models_all.csv`, `external_validation_HF_*.svg`)
Explainability (`SHAP_beeswarm_HF.svg`, `hf_shap_importance.csv`, `SHAP_comparison_3models*.csv`, `SHAP_top20_per_model.csv`)
Biology (`Immune_Infiltration_HF.svg`, `GO_BP_*`, `KEGG_*`, `WGCNA_*`, `Hub_gene_*`)
Archives (`PROCESSED_DATA.zip`, `HEART_ALL_RESULTS.zip`)

7. Notes and limitations

- Sample labels (HF vs Normal) are assigned automatically from GEO sample metadata by keyword matching (`normal`, `control`, `non-fail`, `nf `, `donor`). Anything not matching is labelled   HF. It is worth checking the label counts printed by the notebook.
  
- External-validation cohort GSE42955 is small; interpret its confidence intervals accordingly.
- Not all 20 panel genes may be present on every external platform; the notebook reports gene coverage for each validation run.


8. Notebook and supplementary files

- Notebook: hosted in this GitHub repository.
- Supplementary files (processed data, result tables and figures, i.e. the contents of `PROCESSED_DATA.zip` and `HEART_ALL_RESULTS.zip`): hosted on Figshare.
  Figshare link: *https://doi.org/10.6084/m9.figshare.35037400*

9. Author and contact

- Author: Shadman Sakib
- Contact: sr274276@gmail.com
- ORCID: https://orcid.org/0009-0009-3810-6162
