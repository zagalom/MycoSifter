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
```
---


git clone https://github.com/[your-username]/FungalContaminantRemover.git
cd FungalContaminantRemover

🚀 Quick Start
1. Prepare your input files

Place these three files in your working directory:
File	Format	Required Columns
otu_table.tsv	Tab-separated, rows=OTUs, columns=samples	—
taxonomy.tsv	Tab-separated, rows=OTUs	Kingdom, Phylum, Class, Order, Family, Genus, Species
metadata.tsv	Tab-separated, rows=samples	Must include a DNA concentration column (e.g., DNA_Load in ng/µL)
2. (Optional) Prepare your ecological filter list

Create known_contaminants.txt with one taxon per line using the format LEVEL_TaxonName:
text

# Obligate environmental fungi
P_Glomeromycota        # Entire phylum (arbuscular mycorrhizal)
F_Erysiphaceae         # Entire family (powdery mildews)
G_Inocybe              # Entire genus (ectomycorrhizal)

Supported levels: P (Phylum), C (Class), O (Order), F (Family), G (Genus).
3. Configure and run

Edit the user parameters at the top of FungalContaminantRemover.R:
r

# --- INPUT FILES ---
otu_file <- "otu_table.tsv"
tax_file <- "taxonomy.tsv"
meta_file <- "metadata.tsv"
contaminant_list_file <- "known_contaminants.txt"

# --- DNA LOAD SETTINGS ---
dna_column <- "DNA_Load"

# --- MULTI-THRESHOLD SETTINGS ---
fc_threshold <- 2           
p_threshold <- 0.05         
start_quantile <- 0.10      # 10th percentile = true low biomass
n_thresholds <- 9           
threshold_step <- 5         # Step size matching Qubit precision

# --- USE ECOLOGICAL FILTER? ---
use_ecological_filter <- TRUE

Then source the script:
r

source("FungalContaminantRemover.R")

📊 Output Files

All results are saved in contaminant_analysis_output/:
File	Description
phyloseq_clean.rds	Cleaned phyloseq object (ready for downstream analysis)
full_contaminant_analysis.tsv	All species with scores, classifications, and per-threshold metrics
sensitivity_analysis.tsv	Comparison of 9 stringency levels
removed_species.tsv	Species removed with removal reasons
retained_species.tsv	Species retained
00_REPORT.txt	Human-readable summary report
01_sensitivity_panel.png	Elbow plot: species removed, Procrustes, Spearman by threshold
02_contaminant_load.png	% contaminant reads vs DNA concentration
03_volcano_plot.png	Fold-change vs consistency evidence
04_taxonomic_distribution.png	Contaminants by phylum
05_pcoa_before_after.png	Community ordination before and after cleaning


Key methodological references:

    Salter, S.J., et al. (2014). Reagent and laboratory contamination can critically impact sequence-based microbiome analyses. BMC Biology, 12, 87.

    Davis, N.M., et al. (2018). Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome, 6, 226.


🛠️ Dependencies

    R ≥ 4.0

    phyloseq (Bioconductor)

    vegan (CRAN)

    ggplot2 (CRAN)

    dplyr (CRAN)

    tidyr (CRAN)

    RColorBrewer (CRAN)

📄 License

This project is licensed under the MIT License — see the LICENSE file for details.
🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request or open an Issue.

    Fork the repository

    Create your feature branch (git checkout -b feature/AmazingFeature)

    Commit your changes (git commit -m 'Add some AmazingFeature')

    Push to the branch (git push origin feature/AmazingFeature)

    Open a Pull Request
