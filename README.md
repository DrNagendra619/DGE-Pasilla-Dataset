# DGE-Pasilla-Dataset
DGE Pasilla Dataset
# 🔬 DESeq2 DGE Pipeline: Pasilla Knockdown Dataset

This R script automates a full, classic Differential Gene Expression (DGE) analysis workflow using the **`pasilla`** RNA-Seq dataset. This dataset compares gene expression in *Drosophila* cells after a **pasilla gene knockdown** (treated) versus **untreated** controls.

The pipeline covers data loading, meticulous metadata alignment, `DESeq2` statistical testing, and the generation of essential quality control (QC) and results visualizations.

## 🚀 Key Features

* **Classic DGE Workflow:** Demonstrates the gold-standard analysis steps using the core **`DESeq2`** package.
* **Built-in Data & Metadata:** Utilizes the raw count data and sample annotation files provided directly within the `pasilla` R package.
* **Metadata Alignment:** Includes critical steps to ensure the sample metadata (`colData`) is perfectly aligned and ordered to match the count data columns, a vital step for correct statistical modeling.
* **Low-Count Filtering:** Pre-filters genes with consistently low counts to improve the statistical power of the analysis.
* **Comprehensive Visualization:** Generates three essential plots for interpreting the results:
    1.  **PCA Plot:** For quality control and assessment of sample clustering.
    2.  **Sample Distance Heatmap:** To visualize pairwise sample similarity.
    3.  **Volcano Plot:** To visualize the relationship between fold change and statistical significance.
* **Data Persistence:** Saves the full DGE results (including log2FC, p-values, and adjusted p-values) to a CSV file.

---

## 🔬 Analysis Overview

| Component | Method / Test | Purpose |
| :--- | :--- | :--- |
| **Dataset** | `pasilla` (Built-in) | RNA-Seq data from *Drosophila* cells. |
| **DGE Tool** | `DESeq2` | Statistical method optimized for count data. |
| **Comparison** | Treated (Knockdown) vs. Untreated (Control) | Identifies genes regulated by the pasilla gene. |
| **QC Transformation** | **VST (Variance Stabilizing Transformation)** | Stabilizes variance for accurate distance plots. |
| **Significance** | $\text{p-value} < 0.05$ and $|\text{log2FC}| > 1.5$ | Highlighted in the final Volcano Plot. |

---

## 🛠️ Prerequisites and Setup

### 📦 Packages

The script automatically installs and loads the necessary packages:
* `DESeq2`
* `pasilla`
* `pheatmap`
* `EnhancedVolcano`

### ⚙️ Execution

1.  **Download** the `DGE Pasilla Dataset.R` file.
2.  **Ensure all packages are installed** (the script provides the necessary `BiocManager::install` commands).
3.  **Execute** the script in your R environment:
    ```R
    source("DGE Pasilla Dataset.R")
    ```
    *Note: The plots are displayed in the R graphics device and the CSV file is saved to the current working directory.*

---

## 📁 Output Files

The pipeline generates the following key files in the execution directory:

| Filename | Type | Description |
| :--- | :--- | :--- |
| `DESeq2_pasilla_results.csv` | CSV | Full table of DGE results, sorted by adjusted p-value (`padj`). |

### Visualization Plots

| Plot | Analysis Stage | Description |
| :--- | :--- | :--- |
| **PCA Plot** | QC / Results | Shows sample clustering based on condition (Treated vs. Untreated). |
| **Sample Distance Heatmap** | QC | Visualizes the distance matrix calculated from VST-transformed data. |
| **Volcano Plot** | Results | Displays $\log_2 \text{Fold Change}$ versus $P_{\text{value}}$, highlighting significantly regulated genes. |
