![Coral reef fish](images/AdobeStock_466179536.jpeg)

# Widespread Heterogeneity in Density-Dependent Mortality of Nearshore Fishes

**Authors**  
Adrian C. Stier (University of California, Santa Barbara)  
Craig W. Osenberg (University of Georgia)  

**Contact**  
Adrian Stier — astier@ucsb.edu  
Craig Osenberg — osenberg@uga.edu  

**Funding**  
National Science Foundation (NSF) Grant [OCE-1851510](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1851510)

---

## Overview

This repository contains the complete data synthesis and analysis pipeline for a global meta-analysis on density-dependent mortality in reef fishes. We estimated mortality parameters (α and β) from >30 ecological studies and explored how ecological traits, experimental methods, and phylogenetic history explain variation in density dependence.

This repository supports the manuscript:

> **Stier & Osenberg**  
> *Widespread Heterogeneity in Density-Dependent Mortality of Nearshore Fishes*  
> *Ecology Letters* (in review)

---

## Repository Structure

```
code/           # All R scripts for data loading, modeling, and figure generation
data/           # Raw data, phylogenies, covariates, and metadata
output/         # Model outputs, effect sizes, and summary tables
results/        # Intermediate results, AIC tables, phylogenetic diagnostics
figures/        # Final figures used in the manuscript and supplement
images/         # Supplemental images for README (e.g., reef fish photo)
archive/        # Archived old sensitivity analysis files and deprecated code
```

---

## Quick Start

1. Clone the repository:

```bash
git clone https://github.com/stier-lab/stier-osenberg-2025-density-dependence.git
cd stier-osenberg-2025-density-dependence
R -q -e 'source("code/00_run_all.R")'
```

2. Open R or RStudio and run:

```r
source("code/00_run_all.R")
```

3. All outputs will be saved to the `output/`, `results/`, and `figures/` directories.

---

## Code Overview

| Script                                  | Description                                                         |
|-----------------------------------------|---------------------------------------------------------------------|
| `code/00_run_all.R`                     | Master script to reproduce full pipeline                            |
| `code/0_libraries.R`                    | Loads all required R packages                                       |
| `code/0a_individual_study_estimates.R`  | Generates parameter estimates for individual substudies (NLS2 fits) |
| `code/1_data_phylogeny_loading.R`       | Loads raw data and phylogenies                                      |
| `code/2_beta_estimate.R`                | Estimates α and β from nonlinear fits                               |
| `code/3_phylogenetic_analysis.R`        | Quantifies phylogenetic signal in β                                 |
| `code/4_predators_pairedpredators.R`    | Tests for predator effects on β                                     |
| `code/5_model_selection_comparison.R`   | Compares covariate-based and random-effects models                  |
| `code/6_bivariate_plots_predictors.R`   | Generates bivariate plots of predictors vs. β                       |
| `code/7_global_study_distribution_map.R.R`| Creates world map of study locations                              |
| `code/8_sensitivity_supplement.R`       | Generates Appendix D sensitivity analyses for supplement            |

---

## Data Contents

All required data files are in the `data/` directory.

### Required Data Files (7)

**CSV Files:**
- `all_studies_looped-2024-09-11.csv`: Raw digitized mortality time-series data from all studies
- `covariates-2024-09-30.csv`: Study-level traits and methodological variables (predator presence, duration, etc.)
- `range8132024.csv`: Parameter search ranges for NLS2 grid search optimization

**Phylogenetic Trees:**
- `1.newick.tre`: Phylogeny from Open Tree of Life
- `actinopt_12k_treePL.tre`: Time-calibrated actinopterygii tree from Rabosky et al. (2018)
- `Reef_fish_all.tacted.newick.tre`: Species-level reef fish tree from Siqueira et al. (2020)

**Species Metadata:**
- `unique_species_studies.xlsx`: Species-study metadata for phylogenetic matching

### Additional Reference Files

The `data/` directory also contains individual study files used for reference and validation:
- `schmitt_holbrook_2007.csv`, `shima_osenberg2003.csv`, `wilson_osenberg_2002.csv`: Individual study data for validation
- `manual_densities.csv`: Manually extracted densities for quality assurance

### Archived Files

Old analysis files, deprecated sensitivity scripts, and intermediate results have been moved to `archive/`:
- Old `figs/` and `tables/` directories with sensitivity analysis outputs
- Deprecated sensitivity scripts (`8_sensitivity_master.R`, `8b_sensitivity_asinh.R`)
- Old data files (`all_dat2.csv`, `combined_results_2024-09-17.csv`)
- Legacy figure variants and exploratory outputs

The current analysis uses streamlined scripts and updated file naming (2025-10-15 versions).

---

## Output Contents

### `output/`
Primary analysis outputs including:
- `combined_results_2025-10-15.csv`: Final α and β estimates for all substudies with confidence intervals
- `substudy_results_2025-10-15.csv`: Per-substudy effect sizes and variances
- `substudy_plots_2025-10-15.pdf`: Individual fit plots for all substudies
- `meta_analysis_table.html`: Summary statistics and model results
- `merged-covariates-10-21-24.csv`: Combined dataset with covariates for analysis
- Individual study diagnostic plots in `individual_plots/`

### `results/`
Model diagnostics and statistical outputs:
- `model_comparison_aic.csv`: AICc-based model selection results
- `best_model_summary.html`: Final selected model coefficients and fit statistics
- `predator_presence_summary.html`: Predator effects meta-analysis results
- `predator_all_vs_paired.html`: Predator comparison tables
- `species_effects_no_phylo.csv`: Species-level effects without phylogenetic correction
- `species_variance_CI.csv`, `species_variance_LRT.csv`: Species-level variance components
- `session-info.txt`: R session information for reproducibility
- `manifest.json`: Analysis manifest with git SHA and file list

### `figures/`
Publication-ready figures for manuscript and supplement:

**Main Manuscript Figures:**
- `figure1_orchard_bh_plot.pdf`: Main effect size forest plot with orchard plot
- `Figure3_species_by_family_order_beta_clean.png`: Phylogenetic tree with β estimates by family
- `Figure3_species_by_family_order_beta_darjeeling9.png`: Alternative color palette version
- `Fig4_paired_vs_unpaired.png`: Predator presence effects comparison

**Appendix D (Sensitivity Analyses):**
- `Fig_D1A_LOO.png`: Leave-one-out analysis
- `Fig_D1B_Funnel.png`/`.pdf`: Funnel plot for publication bias
- `Fig_D1C_Egger.png`: Egger's regression test
- `Fig_D1D_Models.png`: Null model comparisons
- `Fig_D1_ALL.pdf`/`.png`: Combined 4-panel sensitivity figure

**Supplementary Figures:**
- `bivar_*.png`: Bivariate predictor vs. β relationships (density, duration, size, etc.)
- `species_specific_effects.png`: Species-level effect sizes
- `body_size_by_species.png`, `duration_by_species.png`: Study characteristics by species
- `paired_exp_obs_by_species.pdf`: Expected vs. observed density for paired studies
- `phylo_beta_circular.png`: Circular phylogenetic tree visualization
- `predators_vs_beta_boxplot.png`: Predator effects boxplot

---

## Dependencies

This project was developed using R (≥ 4.2.0). Core packages include:

**Meta-analysis & Modeling:**
- `metafor` – meta-analysis models
- `nlme`, `MCMCglmm` – mixed and phylogenetic models
- `nls2` – nonlinear least squares grid search

**Phylogenetics:**
- `ape`, `phytools` – phylogenetic manipulation and plotting
- `fishtree` – fish phylogeny database access
- `geiger` – phylogenetic comparative methods

**Data Wrangling:**
- `dplyr`, `tibble`, `tidyr` – data manipulation
- `readr`, `readxl` – file I/O

**Visualization:**
- `ggplot2`, `patchwork` – plotting
- `orchaRd` – orchard plots for meta-analysis
- `gt`, `gtsummary` – table generation
- `maps`, `ggmap` – spatial visualization

**Utilities:**
- `here`, `glue`, `stringr` – path and string manipulation

To install all required packages, run:

```r
source("code/0_libraries.R")
```


---

## Workflow Instructions

To fully reproduce the analysis:

```r
# 1. Clone and enter repo
git clone https://github.com/[your-username]/density_dependence_reef_fish_stier_osenberg.git
cd density_dependence_reef_fish_stier_osenberg

# 2. Open R and run:
source("code/00_run_all.R")
```

The pipeline executes in this order:

1. **Script 0**: Load required R packages
2. **Script 0a**: Generate individual substudy parameter estimates (NLS2 grid search)
   - *Note*: This step is computationally expensive (hours). Results are cached in `output/combined_results_2025-10-15.csv`
   - Set `SKIP_IF_EXISTS <- FALSE` in `0a_individual_study_estimates.R` to force regeneration
3. **Script 1**: Load mortality data and phylogenetic trees
4. **Script 2**: Estimate study-level α and β parameters
5. **Script 3**: Phylogenetic analysis and model selection
6. **Script 4**: Test for predator presence effects
7. **Script 5**: AIC-based model comparison and selection
8. **Script 6**: Generate bivariate predictor plots
9. **Script 7**: Create global study distribution map
10. **Script 8**: Generate sensitivity analyses for supplement (Appendix D)

All outputs are saved to `output/`, `results/`, and `figures/` directories.

### Performance Notes

- **First run**: ~2-4 hours (due to NLS2 grid search in script 0a)
- **Subsequent runs**: ~5-10 minutes (cached results are used)
- To force full regeneration: delete `output/combined_results_2025-10-15.csv` or set `SKIP_IF_EXISTS <- FALSE`

---

## Pipeline Features & Optimization

This repository has been optimized for computational efficiency, cross-platform compatibility, and reproducibility:

### Smart Caching
- **NLS2 Results Caching**: Script 0a automatically checks for existing results in `output/combined_results_2025-10-15.csv`
- Only missing substudies are recomputed, saving hours on subsequent runs
- Control via `SKIP_IF_EXISTS` flag in [0a_individual_study_estimates.R](code/0a_individual_study_estimates.R)

### Cross-Platform Compatibility
- **No Browser Dependencies**: All tables saved as HTML (not PNG) for universal compatibility
- Works on any system without Chrome/Chromium installation
- Tested on macOS, Linux, and Windows

### Streamlined Code
- **Consolidated Sensitivity Analyses**: [8_sensitivity_supplement.R](code/8_sensitivity_supplement.R) generates all Appendix D figures in a single script
- Deprecated exploratory scripts moved to `archive/`
- Clean separation between main figures (`figures/`) and intermediate outputs (`output/`)

### Reproducibility Features
- **Deterministic RNG**: Fixed seed (20250708) for reproducible random sampling
- **Session tracking**: Automatic `results/session-info.txt` and `results/manifest.json` generation
- **Git integration**: Analysis manifest includes git SHA for version tracking
- **Complete provenance**: All figures traceable to specific script and data version

---

## Data Provenance

- **Mortality data**: Digitized from published ecological studies (see manuscript Table S1)
- **Phylogenies**:
  - Rabosky et al. (2018): [https://doi.org/10.1038/s41586-018-0273-1](https://doi.org/10.1038/s41586-018-0273-1)
  - Siqueira et al. (2020): [https://doi.org/10.1038/s41467-020-16498-w](https://doi.org/10.1038/s41467-020-16498-w)
  - OpenTree: [https://tree.opentreeoflife.org/curator/study/view/ot_1592](https://tree.opentreeoflife.org/curator/study/view/ot_1592)

---

## License

- **Data**: [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/)
- **Code**: [MIT License](https://opensource.org/licenses/MIT)

---

## Citation

If using this dataset or codebase, please cite:

> Stier, A. C. & Osenberg, C. W. *Widespread heterogeneity in density-dependent mortality of nearshore fishes.* Ecology Letters (in review).

Upon publication, citation DOIs (manuscript + Dryad) will be added.

---

## ORCID IDs

- Adrian C. Stier — [0000-0002-4704-4145](https://orcid.org/0000-0002-4704-4145)

---


## Reproducibility Summary

This repository follows key principles of open, reproducible science. The table below summarizes the key components of reproducibility for this project:

| Component                     | Status                      | Details                                                                 |
|------------------------------|-----------------------------|-------------------------------------------------------------------------|
| **Code availability**        | ✅ Available                 | All scripts included in the `code/` directory                           |
| **Data availability**        | ✅ Available                 | Raw data in `data/`; each file includes a Dryad-style metadata `.txt`  |
| **Environment**              | ✅ Documented                | R ≥ 4.2.0, dependencies listed in `code/0_libraries.R`                  |
| **Full workflow script**     | ✅ Available                 | `code/00_run_all.R` runs full analysis pipeline                         |
| **Intermediate results**     | ✅ Included                  | See `results/` folder                                                   |
| **Final outputs**            | ✅ Included                  | See `output/` and `figures/`                                            |
| **Manuscript linkage**       | ✅ In progress               | *Ecology Letters*, in review                                           |
| **License**                  | ✅ Open                      | Code: MIT License, Data: CC0                                           |
| **Version control**          | ✅ GitHub                    | All development tracked via Git                                         |

**How to reproduce full analysis:**  
1. Clone repo  
2. Run `source("code/00_run_all.R")` in R  

For issues or questions, contact:  
📧 Adrian Stier – astier@ucsb.edu  
📧 Craig Osenberg – osenberg@uga.edu

</pre>


