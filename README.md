# Breast Cancer Gene Expression Clustering, Differential Expression & Survival Analysis (R)

## 📌 Project Overview

This project performs a comprehensive analysis of **breast cancer gene expression data**, integrating transcriptomic profiles with clinical outcomes to:

- Identify molecular subtypes via **unsupervised clustering**
- Detect **differentially expressed genes (DEGs)** between clusters
- Perform **Gene Ontology (GO) enrichment analysis**
- Evaluate **prognostic significance** using survival analysis  
  (Kaplan–Meier curves and Cox proportional hazards models)

---

## 🔬 Analysis Steps

### 1. Preprocessing
- Between-array normalization using `limma`
- Removal of control probes and low-variance genes

### 2. Clustering
- Hierarchical clustering using **Ward’s method**
- Silhouette analysis to determine the optimal number of clusters

### 3. Visualization
- Heatmap of top variable genes
- Principal Component Analysis (PCA)
- Multidimensional Scaling (MDS)

### 4. Differential Expression
- Linear modeling with empirical Bayes moderation
- Multiple testing correction using **FDR**
- Visualization with volcano plots

### 5. Functional Enrichment
- Gene Ontology (GO) enrichment analysis  
  *(Biological Process)*

### 6. Survival Analysis
- Kaplan–Meier survival curves
- Log-rank tests
- Univariate and multivariate Cox regression models

---

## 🛠 Tools & Libraries

- **R / Bioconductor**
- `limma`
- `cluster`
- `gplots`
- `RColorBrewer`
- `matrixStats`
- `clusterProfiler`
- `org.Hs.eg.db`
- `survival`

---

## 📊 Outputs

- Molecular cluster assignments
- Differentially expressed gene lists
- GO enrichment results
- Publication-ready visualizations:
  - Heatmaps
  - PCA / MDS plots
  - Volcano plots
  - Kaplan–Meier survival curves

---


