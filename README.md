# MycoSifter 🍄🧹

**Automated detection and removal of environmental fungal contaminants from human ITS amplicon sequencing data.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%E2%89%A54.0-blue)](https://www.r-project.org/)

---

## 📖 Overview

Low-biomass microbiome samples are notoriously susceptible to contamination from laboratory reagents, consumables, and environmental sources. In fungal ITS sequencing studies of human mucosal sites (e.g., vaginal, oral, skin), environmental fungi — including soil saprophytes, plant pathogens, and ectomycorrhizal symbionts — can dominate the taxonomic profiles of samples with low biological biomass, masking genuine host-associated signals.

**MycoSifter** is an R-based tool that implements a **data-driven, multi-threshold approach** to identify and remove these environmental contaminants **without requiring negative controls**. It exploits the well-established inverse relationship between total DNA concentration (a proxy for microbial biomass) and contaminant relative abundance (Salter et al., 2014).

---

## 🔬 How It Works

### 1. Multi-Threshold Analysis
Samples are stratified into **Low-DNA** and **High-DNA** groups across 9 DNA concentration thresholds (starting at the 10th percentile, incrementing by 5 ng/µL). At each threshold, species are evaluated using three criteria:

- **Fold-change > 2**: At least twice as abundant in low-DNA vs high-DNA samples
- **Wilcoxon p < 0.05**: Statistically significant enrichment
- **Direction Low > High**: Consistent with contaminant behavior

### 2. Consistency Scoring
Each species receives a **Contamination Score (0–1)** based on the number of thresholds (out of 9) meeting all three criteria. Species are classified as:

| Classification | Thresholds | Evidence |
|:---------------|:----------:|:---------|
| **High** | ≥7/9 | Strong, consistent contaminant signal |
| **Moderate** | 4–6/9 | Intermediate evidence |
| **Low** | 1–3/9 | Weak, sporadic signal |
| **Non-suspect** | 0/9 | No evidence of contamination |

### 3. Ecological Filter (Optional)
Users can provide a list of **obligate environmental taxa** (e.g., mycorrhizal fungi, plant pathogens) that are removed regardless of statistical evidence. This catches contaminants that may not show the classic low-biomass enrichment pattern (e.g., those introduced during sample collection rather than through reagents).

### 4. Sensitivity Analysis
A comprehensive sensitivity analysis compares **9 stringency levels** (≥1/9 to ≥9/9), reporting:
- Species and reads removed at each level
- **Procrustes correlation** (community structure preservation)
- **Spearman's ρ** for alpha diversity metrics

This enables **data-driven threshold selection** rather than arbitrary cutoffs.

### 5. Output
- Cleaned phyloseq object (`.rds`) ready for downstream analysis
- Comprehensive tables of removed and retained species with evidence scores
- Diagnostic visualizations (sensitivity panel, volcano plot, PCoA before/after, taxonomic distribution)

---

## 📦 Installation

```r
# Install required packages
install.packages(c("phyloseq", "vegan", "ggplot2", "dplyr", "tidyr", "RColorBrewer"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")
