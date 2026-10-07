# 🍄 MycoSifter v1.0

**Automated detection and removal of environmental fungal contaminants from human ITS amplicon sequencing data**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R ≥ 4.0](https://img.shields.io/badge/R-%E2%89%A54.0-blue)](https://www.r-project.org/)

**Author:** Mariana M. Zagalo Fernandes  
**Repository:** https://github.com/zagalom/MycoSifter  
**License:** MIT

---

## 📖 Overview

Low-biomass microbiome samples are notoriously susceptible to contamination from laboratory reagents, consumables, and environmental sources. In fungal ITS sequencing studies of human mucosal sites, environmental fungi – including soil saprophytes, plant pathogens, and ectomycorrhizal symbionts – can dominate the taxonomic profiles of samples with low biological biomass, masking genuine host-associated signals.

**MycoSifter** is an R-based tool that implements a **data-driven, multi-threshold approach** to identify and remove these environmental contaminants **without requiring negative controls**. It exploits the well-established inverse relationship between total DNA concentration (a proxy for microbial biomass) and contaminant relative abundance (Salter *et al.*, 2014).

Optionally, MycoSifter applies an **ecological filter** using a curated list of obligately environmental fungal taxa (phyla, orders, families, and genera with no known capacity to colonise human mucosae), capturing contaminants that may not show the classic low-biomass enrichment pattern.

---

## 🔬 How It Works

### 1. Multi-Threshold Contaminant Detection

Samples are stratified into **Low-DNA** and **High-DNA** groups across multiple DNA concentration thresholds. At each threshold, species are evaluated using three criteria:

- **Fold-change**: enrichment in low-DNA *vs.* high-DNA samples (default: > 2-fold)
- **Wilcoxon test**: statistical significance of enrichment (default: *p* < 0.05)
- **Direction Low > High**: consistent with contaminant behaviour

The number and spacing of thresholds is **fully configurable**, allowing you to adapt the analysis to your data distribution and study design.

### 2. Consistency Scoring & Sensitivity Analysis

Each species receives a **Contamination Score (0–1)** based on the proportion of thresholds meeting all three criteria. A sensitivity analysis compares multiple stringency levels and reports:

- Species and reads removed at each level
- **Procrustes correlation** (community structure preservation)
- **Spearman's ρ** for observed richness and Shannon diversity

This enables **data-driven threshold selection** rather than arbitrary cutoffs.

### 3. Ecological Filter (Optional)

Users can provide a list of **obligate environmental taxa** (e.g., mycorrhizal fungi, plant pathogens) that are removed regardless of statistical evidence. This catches contaminants introduced during sample collection that may not show the classic reagent-contamination pattern.

### 4. Comprehensive Output

- Cleaned phyloseq object (`.rds`) ready for downstream analysis
- Complete tables of removed and retained species with evidence scores
- Sensitivity analysis table for threshold selection
- Diagnostic visualisations (sensitivity panel, contaminant load, volcano plot, taxonomic distribution, PCoA before/after)
- Human-readable summary report with reproducible parameters

---

## 📦 Installation

### Step 1: Install R packages

```r
# Install required CRAN packages
install.packages(c("phyloseq", "vegan", "ggplot2", "dplyr", "tidyr", "RColorBrewer"))
```

### Step 2: Clone the repository

```bash
git clone https://github.com/zagalom/MycoSifter.git
cd MycoSifter
```

---

## 🚀 Quick Start

### 1. Prepare your input files

Place these three files in your working directory:

| File | Format | Required Columns |
|------|--------|------------------|
| `otu_table.tsv` | Tab-separated, rows=OTUs, columns=samples | — |
| `taxonomy.tsv` | Tab-separated, rows=OTUs | Kingdom, Phylum, Class, Order, Family, Genus, Species |
| `metadata.tsv` | Tab-separated, rows=samples | Must include a DNA concentration column (e.g., `DNA_Load` in ng/µL) |

### 2. (Optional) Prepare your ecological filter list

Create `known_contaminants.txt` with one taxon per line using the format `LEVEL_TaxonName`:

```text
# Obligate environmental fungi (phyla, orders, families, genera)
P_Glomeromycota        # Entire phylum (arbuscular mycorrhizal)
F_Erysiphaceae         # Entire family (powdery mildews)
G_Inocybe              # Entire genus (ectomycorrhizal)
```

**Supported levels:** `P` (Phylum), `C` (Class), `O` (Order), `F` (Family), `G` (Genus)

### 3. Configure the script parameters

Edit the `USER CONFIGURATION` section at the top of `MycoSifter` to set your analysis parameters. See the **[Configuration Guide](#-configuration-guide)** below for detailed descriptions of each option.

### 4. Run the analysis

```bash
Rscript MycoSifter
```

Or in R:

```r
source("MycoSifter")
```

The script will:
1. Load your data and align samples with metadata
2. Agglomerate OTUs to species level (preserving full taxonomy)
3. Apply quality filters
4. Run multi-threshold contaminant detection
5. Perform sensitivity analysis
6. (Optionally) Apply the ecological filter
7. Generate a cleaned phyloseq object and summary report
8. Create diagnostic figures

---

## ⚙️ Configuration Guide

All user-facing options are in the `USER CONFIGURATION` section (section 0) of the `MycoSifter` script. Edit these values to customise the analysis to your data and research questions.

### Input Files

```r
otu_file <- "otu_table.tsv"
tax_file <- "taxonomy.tsv"
meta_file <- "metadata.tsv"
contaminant_list_file <- "known_contaminants.txt"
```

**Description:** Path to your input files.

- `otu_file`: OTU/ASV abundance table (tab-separated, rows=OTUs/ASVs, columns=samples, first column=OTU names)
- `tax_file`: Taxonomy table (tab-separated, rows=OTUs/ASVs, columns=taxonomic ranks)
- `meta_file`: Sample metadata (tab-separated, rows=samples, must include DNA concentration column)
- `contaminant_list_file`: List of known environmental taxa (one per line, optional; can be disabled)

---

### DNA Load Settings

```r
dna_column <- "DNA_Load"
```

**Description:** Name of the DNA concentration column in your metadata file.

- **Default:** `"DNA_Load"`
- **Notes:** If this column is not found, the script will attempt to auto-detect columns containing "dna", "load", "concentr", "qubit", or "quantif" (case-insensitive).
- **Units:** Must be in ng/µL to match the threshold parameters.

---

### Statistical Thresholds

```r
fc_threshold <- 2
p_threshold  <- 0.05
```

**Description:** Criteria for flagging a species as a contaminant candidate at each DNA concentration threshold.

- `fc_threshold`: Fold-change threshold (Low DNA / High DNA). Species with a fold-change **> this value** in low-DNA samples are considered suspect.
  - **Default:** `2` (2-fold enrichment)
  - **Recommendations:** 2–5 (higher = more conservative)

- `p_threshold`: Wilcoxon test p-value threshold.
  - **Default:** `0.05` (0.05 = 5% significance level)
  - **Note:** Adjusted (FDR) p-values can be used instead via the `p_method` parameter below.

---

### Multi-Threshold Settings

```r
threshold_mode <- "quantile"
n_thresholds <- 9
start_quantile <- 0.10
quantile_max <- 0.50
threshold_max_frac <- 0.50
threshold_step <- NULL
custom_thresholds <- NULL
```

**Description:** Control how DNA concentration thresholds are generated.

#### `threshold_mode`

Determines how to generate a series of DNA concentration cutoff values.

- **`"quantile"` (recommended)**
  - Generates evenly-spaced percentiles of DNA load distribution
  - **Pros:** Balanced groups at each threshold; approximately independent tests
  - **Cons:** Uneven DNA spacing in ng/µL
  - **Best for:** exploratory analysis and when sample sizes are balanced
  - **Interpretation:** `N_Suspect` is interpretable as the number of independent thresholds supporting contaminant status

- **`"uniform"` (classic low-biomass)**
  - Generates thresholds at fixed intervals (e.g., 5 ng/µL steps)
  - Based on Salter *et al.* (2014) and Karstens *et al.* (2019)
  - **Pros:** Intuitive interpretation; biologically meaningful spacing
  - **Cons:** Unbalanced groups (few samples in low-DNA bins); correlated tests
  - **Best for:** small sample sizes; specific DNA ranges
  - **Interpretation:** `N_Suspect` reflects **persistence** of signal in the low-biomass range, not independent evidence

- **`"linear"`**
  - Generates evenly-spaced thresholds in ng/µL across the data range
  - **Pros:** Balance between quantile and uniform modes
  - **Cons:** Still somewhat correlated at extremes
  - **Best for:** intermediate scenarios

- **`"custom"`**
  - Manually specify exact thresholds in the `custom_thresholds` parameter
  - **Best for:** hypothesis-driven analysis with specific DNA cutoffs

#### `n_thresholds`

Number of thresholds to generate (for `quantile`, `uniform`, and `linear` modes).

- **Default:** `9`
- **Range:** 3–15 (minimum 3 required for sensitivity analysis)
- **Notes:** More thresholds = finer resolution but longer computation time

#### `start_quantile`

Lower bound for threshold generation (only for `quantile`, `uniform`, and `linear` modes).

- **Default:** `0.10` (10th percentile)
- **Range:** 0–1
- **Interpretation:** Samples below this quantile are considered "low DNA"
- **Recommendation:** 0.10–0.25 to capture true low-biomass samples

#### `quantile_max`

Upper bound for threshold generation (only for `quantile` mode).

- **Default:** `0.50` (median)
- **Range:** 0–1
- **Interpretation:** Defines the highest percentile to use for threshold generation
- **Notes:** Beyond 0.50, the "High DNA" group becomes too small (<10% of samples), reducing statistical power
- **Recommendation:** Keep ≤ 0.50

#### `threshold_max_frac`

Absolute cap on thresholds as a fraction of maximum DNA concentration (applied to **all** modes).

- **Default:** `0.50`
- **Range:** 0–1
- **Example:** If max DNA = 100 ng/µL and `threshold_max_frac = 0.50`, no threshold > 50 ng/µL will be used
- **Purpose:** Prevents creating thresholds above biologically meaningful DNA loads

#### `threshold_step`

(For `uniform` and `linear` modes) Fixed step size in ng/µL between consecutive thresholds.

- **Default:** `NULL` (data-driven: automatically calculated)
- **Example:** `threshold_step <- 5` creates thresholds at 5, 10, 15, 20 ng/µL…
- **Notes:** If `NULL`, the script calculates step size as: `(max_dna - start_dna) / n_thresholds`

#### `custom_thresholds`

(For `threshold_mode = "custom"`) Manually specify exact thresholds as a vector.

- **Default:** `NULL`
- **Example:** `custom_thresholds <- c(5, 10, 20, 40, 80)`
- **Notes:** Thresholds should be in increasing order and in ng/µL

---

### Statistical Method

```r
p_method <- "raw"
```

**Description:** How to handle multiple comparisons when evaluating species across thresholds.

- **`"raw"` (default)**
  - Uses raw Wilcoxon *p*-values (not adjusted)
  - **Pros:** More sensitive; faster computation
  - **Cons:** Higher false positive rate with many thresholds
  - **Best for:** exploratory analysis

- **`"fdr"`**
  - Adjusts *p*-values using Benjamini-Hochberg FDR correction
  - **Pros:** Controls false discovery rate; more conservative
  - **Cons:** Slightly lower sensitivity
  - **Best for:** publication-ready, hypothesis-driven analysis

---

### Ecological Filter

```r
use_ecological_filter <- TRUE
```

**Description:** Enable/disable the ecological filter (removal of known environmental taxa).

- **`TRUE`**: Enables the ecological filter
  - Requires a `known_contaminants.txt` file
  - Taxa listed in the file are removed **in addition to** threshold-based contaminants
  - Useful for catching obligate environmental taxa that don't show low-biomass enrichment

- **`FALSE`**: Disables the ecological filter
  - Only threshold-based detection is used
  - If `known_contaminants.txt` does not exist, the script automatically sets this to `FALSE`

---

### Cutoff Selection

```r
interactive_cutoff <- TRUE
min_thresholds_manual <- NULL
```

**Description:** Control how the final removal threshold is selected.

#### `interactive_cutoff`

- **`TRUE` (default)**
  - After sensitivity analysis, the script prompts you to choose the minimum number of thresholds (k) required to flag a taxon as contaminant
  - The prompt shows all available k values with species/reads removed and statistical preservation metrics
  - **Example prompt:**
    ```
    Please choose the minimum number of thresholds required to classify a taxon as contaminant-like.
    
    0 = ecological filter only (no threshold-based removal)
    1 = ≥1/9 thresholds  (removes 45 species, 8.2% reads)
    2 = ≥2/9 thresholds  (removes 28 species, 4.1% reads)
    ...
    ```
  - **Best for:** interactive sessions (RStudio, command line)

- **`FALSE`**
  - Uses the value in `min_thresholds_manual` (see below)
  - Non-interactive mode suitable for scripted/batch analysis
  - **Best for:** Rscript execution, automated workflows

#### `min_thresholds_manual`

- **Default:** `NULL` (uses default k=1)
- **Range:** 0 to `n_thresholds`
- **Interpretation:**
  - `0`: Remove only taxa flagged by ecological filter (no threshold-based removal)
  - `1`: Remove taxa flagged at ≥1 threshold (most liberal)
  - `n`: Remove taxa flagged at ≥n thresholds (most stringent; n = `n_thresholds`)
- **Used when:** `interactive_cutoff = FALSE`
- **Recommendation:** Start with k=1–2 for exploratory analysis; use higher k for more conservative cleanup

---

### Output

```r
output_dir <- "contaminant_analysis_output_New"
```

**Description:** Directory name for all output files.

- **Default:** `"contaminant_analysis_output_New"`
- **Notes:** Directory is created automatically if it doesn't exist

---

## 📋 Example Configuration Scenarios

### Scenario 1: Conservative threshold-based analysis (publication-ready)

```r
# Statistical thresholds
fc_threshold <- 3               # Stricter fold-change
p_threshold  <- 0.01            # Stricter p-value
p_method <- "fdr"               # FDR correction for multiple comparisons

# Threshold mode
threshold_mode <- "quantile"    # Balanced groups
n_thresholds <- 9
start_quantile <- 0.10
quantile_max <- 0.50

# Cutoff selection
interactive_cutoff <- FALSE
min_thresholds_manual <- 3      # Require ≥3 thresholds (high confidence)

# Ecological filter
use_ecological_filter <- FALSE  # Don't use; rely on statistics

output_dir <- "analysis_conservative"
```

**Result:** Only species with **strong evidence** across multiple independent thresholds are removed.

---

### Scenario 2: Exploratory low-biomass analysis (hypothesis generation)

```r
# Statistical thresholds
fc_threshold <- 2               # Standard fold-change
p_threshold  <- 0.05            # Standard p-value
p_method <- "raw"               # Raw p-values (faster, more sensitive)

# Threshold mode
threshold_mode <- "uniform"     # Classic low-biomass approach
n_thresholds <- 7
start_quantile <- 0.10
threshold_step <- 5             # 5 ng/µL steps

# Cutoff selection
interactive_cutoff <- TRUE      # User selects after seeing results

# Ecological filter
use_ecological_filter <- TRUE   # Also remove known environmental taxa

output_dir <- "analysis_exploratory"
```

**Result:** Sensitive detection with curated ecological knowledge; user chooses stringency interactively.

---

### Scenario 3: Custom thresholds (hypothesis-driven)

```r
# Hypothesis: contaminants enriched at DNA < 2 ng/µL and < 5 ng/µL

threshold_mode <- "custom"
custom_thresholds <- c(2, 5, 10, 20)  # Specific DNA cutoffs

# Other parameters
fc_threshold <- 2.5
p_threshold  <- 0.05
p_method <- "raw"

interactive_cutoff <- FALSE
min_thresholds_manual <- 2      # Must be flagged at ≥2 custom thresholds

output_dir <- "analysis_hypothesis"
```

**Result:** Analysis focused on your specific DNA concentration ranges of interest.

---

## 📊 Output Files

All results are saved in the specified output directory (default: `contaminant_analysis_output_New/`):

### Data Tables

| File | Description |
|------|-------------|
| `phyloseq_clean.rds` | **Main output.** Clean phyloseq object with contaminants removed, ready for downstream analysis (beta diversity, differential abundance, etc.) |
| `full_contaminant_analysis.tsv` | **Complete reference table.** All species with: full taxonomy, read counts, prevalence, per-threshold fold-change & p-values, contamination score, classification, removal flag |
| `sensitivity_analysis.tsv` | Comparison of all stringency levels (Eco-only, ≥1/n, ≥2/n, …, ≥n/n): species/reads removed, Procrustes correlation, Spearman's ρ |
| `removed_species.tsv` | Species flagged as contaminants with removal reason (threshold-based, ecological, or both) |
| `retained_species.tsv` | Species retained in the cleaned dataset |

### Diagnostic Figures

| File | Description |
|------|-------------|
| `01_sensitivity_panel.png` | Three elbow plots showing: (left) % species removed, (middle) Procrustes correlation, (right) Spearman's ρ for observed richness. Used to select an appropriate cutoff. |
| `02_contaminant_load.png` | Scatter plot: DNA concentration (x-axis) vs. % contaminant reads (y-axis). Low-DNA samples should show high contamination. Includes LOESS smoothing curve. |
| `03_volcano_plot.png` | Log₂(fold-change) vs. -log₁₀(N_Suspect): visualise which species show the strongest evidence of contamination. Coloured by classification. |
| `04_taxonomic_distribution.png` | Horizontal bar chart: top 10 phyla among removed species (absolute counts). Identifies which taxonomic groups are most affected. |
| `05_pcoa_before_after.png` | Side-by-side PCoA plots of community composition (Bray-Curtis): left = before cleaning, right = after cleaning. Helps assess if community structure is preserved. |

### Report

| File | Description |
|------|-------------|
| `00_REPORT.txt` | Comprehensive human-readable summary with: input data summary, sample/OTU alignment, species agglomeration, quality filters, threshold parameters, sensitivity analysis, contamination classification, removal impact, top contaminants, taxonomic breakdown, file descriptions, and complete pipeline summary. Reproducible and suitable for inclusion in supplementary materials. |

---

## 🧬 Interpretation Guide

### Understanding the Contamination Score

Each species receives a **Contamination Score (0–1)** calculated as:

```
ContamScore = 0.4 × (N_Suspect / n_thresholds) 
            + 0.3 × (N_FC > fc_threshold / n_thresholds) 
            + 0.3 × (N_p < p_threshold / n_thresholds)
```

Where:
- **N_Suspect:** Number of thresholds where the species meets *all* three criteria (fold-change, p-value, direction)
- **N_FC:** Number of thresholds with fold-change > `fc_threshold`
- **N_p:** Number of thresholds with p-value < `p_threshold`

**Interpretation:**
- **0.7–1.0:** High confidence contaminant (consistent evidence across thresholds)
- **0.4–0.7:** Moderate confidence (patchy evidence)
- **0.1–0.4:** Low confidence (weak evidence)
- **< 0.1:** Non-suspect (little/no evidence)

---

### Understanding Classification

Species are classified into **five categories**:

- **High (≥7/n):** Flagged at ≥7 of 9 thresholds. Very strong evidence.
- **Moderate (4–6/n):** Flagged at 4–6 thresholds. Good evidence.
- **Low (1–3/n):** Flagged at 1–3 thresholds. Weak/variable evidence.
- **Ecological only:** Flagged by the ecological filter but not by statistics. Known environmental taxon.
- **Non-suspect:** No evidence of contamination.

---

### Understanding the Sensitivity Table

The sensitivity analysis compares outcomes at different stringency levels (k):

- **k=0:** Remove only ecological filter hits (most liberal)
- **k=1:** Remove species flagged at ≥1 threshold
- **k=9:** Remove species flagged at ≥9/9 thresholds (most stringent)

**Key metrics:**
- **Procrustes correlation:** How similar the PCoA ordinations are before/after cleaning (0–1, higher is better). ≥0.98 indicates community structure is largely preserved.
- **Spearman's ρ (Observed richness):** Correlation of species counts before/after cleaning (0–1, higher is better). ≥0.95 indicates alpha diversity is preserved.

**Decision principle:** Choose the **highest k** where Procrustes ≥ 0.98 *and* Spearman ≥ 0.95. This balances contaminant removal with preservation of true biological signal.

---

## 📚 Key Methodological References

- **Salter, S.J., et al.** (2014). Reagent and laboratory contamination can critically impact sequence-based microbiome analyses. *BMC Biology*, 12, 87. https://doi.org/10.1186/s12915-014-0087-z

- **Davis, N.M., et al.** (2018). Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. *Microbiome*, 6, 226. https://doi.org/10.1186/s40168-018-0605-2

- **Karstens, L., et al.** (2019). Controlling for contaminants in low-biomass microbiome studies remains challenging. *mSystems*, 4, e00290-19. https://doi.org/10.1128/mSystems.00290-19

---

## 🛠️ System Requirements

- **R:** ≥ 4.0
- **RAM:** ≥ 4 GB (for typical-sized datasets; larger datasets may require more)
- **Disk:** ≥ 500 MB for output (varies by dataset size)

### Dependencies

| Package | Source | Purpose |
|---------|--------|---------|
| `phyloseq` | Bioconductor | Microbiome data structures and operations |
| `vegan` | CRAN | Community ecology and diversity calculations |
| `ggplot2` | CRAN | Visualization |
| `dplyr` | CRAN | Data manipulation |
| `tidyr` | CRAN | Data tidying |
| `RColorBrewer` | CRAN | Colour palettes |

---

## ❓ FAQ

### Q: Do I need negative controls?
**A:** No! MycoSifter works without negative controls, exploiting the inverse relationship between DNA concentration and contamination. However, if you have negative controls, you can use them to validate the ecological filter.

### Q: What if my DNA Load column has a different name?
**A:** The script auto-detects columns containing "dna", "load", "concentr", "qubit", or "quantif" (case-insensitive). Alternatively, set `dna_column` to the exact column name.

### Q: How do I choose between `threshold_mode` options?
**A:** 
- **`"quantile"`** (default): Start here for most analyses. Balanced and approximately independent tests.
- **`"uniform"`**: Use if you want biologically meaningful DNA thresholds (e.g., 5 ng/µL steps) or if your sample sizes are very small.
- **`"custom"`**: Use if you have specific hypotheses about contamination thresholds.

### Q: What's the difference between `p_method = "raw"` and `"fdr"`?
**A:**
- **`"raw"`**: No correction for multiple comparisons. Faster and more sensitive, but higher false positive rate.
- **`"fdr"`**: Benjamini-Hochberg correction. More conservative and suitable for publication.

### Q: How do I know what `min_thresholds` to choose?
**A:** Look at the sensitivity table (`sensitivity_analysis.tsv`) and the elbow plot (`01_sensitivity_panel.png`). Choose the **highest k** where Procrustes ≥ 0.98 and Spearman's ρ ≥ 0.95. These thresholds indicate good balance between removing contaminants and preserving true biology.

### Q: Can I re-run with different parameters?
**A:** Yes! Simply edit the configuration, change `output_dir` to a new name, and run again. All outputs are isolated in their respective directories, making it easy to compare results across parameter sets.

### Q: How long does the analysis take?
**A:** Depends on dataset size. Typical runs (100–500 samples, 200–2000 species) take 5–20 minutes. Larger datasets may take an hour or more.

---

## 🐛 Troubleshooting

### Error: "No DNA concentration column found"
- **Cause:** The column name in `metadata.tsv` doesn't match `dna_column`
- **Solution:** Check the exact name of your DNA concentration column and update `dna_column`, or rename the column in your metadata file

### Error: "Fewer than 3 valid thresholds remain"
- **Cause:** Your threshold generation parameters produced < 3 usable thresholds
- **Solution:** Increase `quantile_max` or `threshold_max_frac`, or reduce `n_thresholds`

### Interactive prompt not appearing (non-interactive mode)
- **Cause:** Running in a non-interactive environment (Rscript, cron job, cluster)
- **Solution:** Set `interactive_cutoff <- FALSE` and choose `min_thresholds_manual`

### Procrustes correlation very low (< 0.9)
- **Cause:** Removing too many species is distorting the community structure
- **Solution:** Use a higher `min_thresholds` value to be more stringent, or relax the statistical thresholds

---

## 📄 License

This project is licensed under the **MIT License** – see the LICENSE file for details.

---

## 👋 Contributing

Found a bug? Have a feature request? Please open an issue or submit a pull request!

---

## 📬 Contact

**Mariana M. Zagalo Fernandes**  
mariana@example.com  
https://github.com/zagalom/MycoSifter

---

## 🙏 Acknowledgements

This tool builds on the foundational work of:
- **Salter *et al.* (2014)** – identifying contamination in low-biomass microbiome studies
- **Davis *et al.* (2018)** – statistical decontamination methods
- **Karstens *et al.* (2019)** – contamination control strategies

Thanks to the R microbiome community and the creators of `phyloseq`, `vegan`, and `dplyr`.

---

**Last updated:** October 2026  
**MycoSifter v1.0**
