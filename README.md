# 🍄 MycoSifter v1.0

**Automated detection and removal of environmental fungal contaminants
from human ITS amplicon sequencing data**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R ≥ 4.0](https://img.shields.io/badge/R-%E2%89%A54.0-blue)](https://www.r-project.org/)

---

## 📖 Overview

Low‑biomass microbiome samples are notoriously susceptible to contamination
from laboratory reagents, consumables, and environmental sources. In fungal
ITS sequencing studies of human mucosal sites, environmental fungi – including
soil saprophytes, plant pathogens, and ectomycorrhizal symbionts – can dominate
the taxonomic profiles of samples with low biological biomass, masking genuine
host‑associated signals.

**MycoSifter** is an R‑based tool that implements a **data‑driven,
multi‑threshold approach** to identify and remove these environmental
contaminants **without requiring negative controls**. It exploits the
well‑established inverse relationship between total DNA concentration
(a proxy for microbial biomass) and contaminant relative abundance
(Salter *et al.*, 2014).

Optionally, MycoSifter applies an **ecological filter** using a curated list of
obligately environmental fungal taxa (phyla, orders, families, and genera with
no known capacity to colonise human mucosae), capturing contaminants that may
not show the classic low‑biomass enrichment pattern.

---

## 🔬 How It Works

### 1. Multi‑Threshold Contaminant Detection
Samples are stratified into **Low‑DNA** and **High‑DNA** groups across up to
nine DNA concentration thresholds (starting at the 10th percentile, incrementing
by 5 ng µL⁻¹). At each threshold, species are evaluated using three criteria:

- **Fold‑change > 2**: at least twice as abundant in low‑DNA *vs.* high‑DNA samples
- **Wilcoxon *p* < 0.05**: statistically significant enrichment
- **Direction Low > High**: consistent with contaminant behaviour

### 2. Consistency Scoring & Sensitivity Analysis
Each species receives a **Contamination Score (0–1)** based on the proportion
of thresholds (out of 9) meeting all three criteria. A sensitivity analysis
compares **ten stringency levels** (Eco‑only, ≥1/9 … ≥9/9) and reports:

- Species and reads removed at each level
- **Procrustes correlation** (community structure preservation)
- **Spearman's *ρ*** for observed richness and Shannon diversity

This enables **data‑driven threshold selection** rather than arbitrary cutoffs.

### 3. Ecological Filter (Optional)
Users can provide a list of **obligate environmental taxa** (e.g., mycorrhizal
fungi, plant pathogens) that are removed regardless of statistical evidence.
This catches contaminants introduced during sample collection that may not show
the classic reagent‑contamination pattern.

### 4. Comprehensive Output
- Cleaned phyloseq object (`.rds`) ready for downstream analysis
- Complete tables of removed and retained species with evidence scores
- Sensitivity analysis table for threshold selection
- Diagnostic visualisations (sensitivity panel, contaminant load, volcano plot,
  taxonomic distribution, PCoA before/after)
- Human‑readable summary report

---

## 📦 Installation

```r
# Install required CRAN packages
install.packages(c("phyloseq", "vegan", "ggplot2", "dplyr", "tidyr", "RColorBrewer"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")

Then clone this repository:
bash

git clone https://github.com/[your-username]/MycoSifter.git
cd MycoSifter

🚀 Quick Start
1. Prepare your input files

Place these three files in your working directory:
File	Format	Required Columns
otu_table.tsv	Tab‑separated, rows=OTUs, columns=samples	—
taxonomy.tsv	Tab‑separated, rows=OTUs	Kingdom, Phylum, Class, Order, Family, Genus, Species
metadata.tsv	Tab‑separated, rows=samples	Must include a DNA concentration column (e.g., DNA_Load in ng µL⁻¹)
2. (Optional) Prepare your ecological filter list

Create known_contaminants.txt with one taxon per line using the format
LEVEL_TaxonName:
text

# Obligate environmental fungi
P_Glomeromycota        # Entire phylum (arbuscular mycorrhizal)
F_Erysiphaceae         # Entire family (powdery mildews)
G_Inocybe              # Entire genus (ectomycorrhizal)

Supported levels: P (Phylum), C (Class), O (Order), F (Family), G (Genus).
3. Configure and run

Edit the user parameters at the top of MycoSifter.R:
r

# --- INPUT FILES ---
otu_file <- "otu_table.tsv"
tax_file <- "taxonomy.tsv"
meta_file <- "metadata.tsv"
contaminant_list_file <- "known_contaminants.txt"

# --- DNA LOAD SETTINGS ---
dna_column <- "DNA_Load"

# --- MULTI‑THRESHOLD SETTINGS ---
fc_threshold <- 2            # Fold‑change (Low/High DNA)
p_threshold <- 0.05          # Wilcoxon p‑value
start_quantile <- 0.10       # 10th percentile = true low biomass
n_thresholds <- 9            # Number of thresholds
threshold_step <- 5          # Step size matching Qubit precision (ng µL⁻¹)

# --- USE ECOLOGICAL FILTER? ---
use_ecological_filter <- TRUE

# --- OUTPUT ---
output_dir <- "contaminant_analysis_output"

Then source the script:
r

source("MycoSifter.R")

📊 Output Files

All results are saved in contaminant_analysis_output/:
File	Description
phyloseq_clean.rds	Cleaned phyloseq object for downstream analysis
full_contaminant_analysis.tsv	All species with scores, classifications, and per‑threshold metrics
sensitivity_analysis.tsv	Comparison of 10 stringency levels
removed_species.tsv	Species removed with removal reasons
retained_species.tsv	Species retained
00_REPORT.txt	Human‑readable summary report
01_sensitivity_panel.png	Elbow plots for threshold selection
02_contaminant_load.png	% contaminant reads vs. DNA concentration
03_volcano_plot.png	Fold‑change vs. evidence
04_taxonomic_distribution.png	Contaminants by phylum
05_pcoa_before_after.png	Community ordination before and after cleaning

Key methodological references:

    Salter, S.J., et al. (2014). Reagent and laboratory contamination can
    critically impact sequence‑based microbiome analyses. BMC Biology, 12, 87.

    Davis, N.M., et al. (2018). Simple statistical identification and removal
    of contaminant sequences in marker‑gene and metagenomics data. Microbiome,
    6, 226.

    Karstens, L., et al. (2019). Controlling for contaminants in low‑biomass
    microbiome studies. mSystems, 4, e00290‑19.

🛠️ Dependencies

    R ≥ 4.0

    phyloseq (Bioconductor)

    vegan (CRAN)

    ggplot2 (CRAN)

    dplyr (CRAN)

    tidyr (CRAN)

    RColorBrewer (CRAN)
