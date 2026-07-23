# Development, Health-System, and Spatial Determinants of Cervical Cancer Incidence and Mortality

A multi-robustness-checked ecological analysis of global cervical cancer burden (GLOBOCAN 2022/2024), linked to the UNDP Human Development Index, HPV vaccination coverage, national cervical screening programs, female smoking prevalence, HIV prevalence, and death-registration completeness, with paired temporal, mediation, and spatial-lag robustness checks.

This repository contains **all raw data, processed datasets, analysis scripts, figures, a predictive-modeling notebook, and the manuscript** for full reproducibility.

---

## Repository Structure

```
.
├── README.md                          <- you are here
├── requirements.txt                   <- Python package dependencies
├── data/
│   ├── raw/                           <- original, unmodified downloaded source files
│   └── processed/                     <- merged/derived analytic datasets (numbered by pipeline stage)
├── scripts/
│   ├── 01-17_*.py                     <- analysis pipeline, in execution order
│   └── figures/                       <- scripts that generate every figure in the manuscript
├── figures/                           <- final rendered figures (PNG), as they appear in the manuscript
├── notebooks/
│   └── cervical_cancer_prediction_model.ipynb   <- supplementary ML prediction notebook (Random Forest / Gradient Boosting)
└── manuscript/
    ├── Cervical_Cancer_IEEE_Format_Draft.docx      <- IEEE-structured version (sections I-VII, appendices A-F)
    └── Cervical_Cancer_Ecological_Study_Draft.docx <- original journal-article-format version
```

---

## Data Sources — Full Detail

Every dataset below was downloaded as a direct, publicly accessible file (no authentication, no restricted access). Country-level merges use ISO-3 alpha codes throughout.

| # | File (in `data/raw/`) | What it contains | Original source / publisher | How it was obtained | Coverage | Notes / known quirks |
|---|---|---|---|---|---|---|
| 1 | `globocan_2024_incidence.csv` | Cervical cancer (ICD-10 C53) new cases, ASR (World), crude rate, cumulative risk, by country, females, all ages | WHO/IARC **Global Cancer Observatory**, Cancer Today platform, GLOBOCAN 2024 estimates | Manually exported from the GCO Cancer Today interactive tool (Display: Populations; Measure: Incidence; Sex: Females; Cancer: Cervix uteri) — https://gco.iarc.who.int/today | 186 countries/territories + 1 world-aggregate row | "Country" column in the raw export contains numeric ISO/UN population codes, not names — use the "Label" column for country names. Namibia's ISO-3 code "NAM" is fine, but be aware some downstream tools misread the string "NA" (unrelated field) as a missing value. |
| 2 | `globocan_2024_mortality.csv` | Same structure as above, Measure: Mortality | Same as above | Same as above | Same as above | Same as above |
| 3 | `globocan_2022_incidence.csv` | Cervical cancer ASR (World) incidence per country, GLOBOCAN **2022** cycle (used only for the paired temporal-panel comparison) | WHO/IARC GCO, mirrored by **Our World in Data** | OWID grapher "Rate of new cervical cancer cases" (GCO 2022 edition): https://ourworldindata.org/grapher/rate-of-new-cervical-cancer-cases-gco | 185 countries | This is a *separate* GLOBOCAN edition from the 2024 file above; used specifically to test whether the HDI-incidence gap changed between cycles (Section V-A / Fig. 5). Do not merge this with the 2024 file as if they were the same measurement occasion. |
| 4 | `undp_hdi.csv` | Human Development Index, all countries, 1990-2023 (panel) | **UNDP Human Development Report Office** | OWID grapher mirror of the UNDP HDI series: https://ourworldindata.org/grapher/human-development-index.csv | 195 countries × 34 years | Contains aggregate rows (regions, income groups, "World") with a blank `Code` column — these must be dropped before merging (they are not countries). Only the most recent available year (2023) per country was used. |
| 5 | `undp_gni_per_capita.csv` | Gross National Income per capita (PPP), all countries, panel | UNDP Human Development Report (World Bank PPP basis) | OWID grapher mirror: https://ourworldindata.org/grapher/gross-national-income-per-capita-undp.csv | 205 countries × multiple years | **Ultimately excluded from the final multivariable model** — collinear with HDI (VIF ≈ 15-16, since GNI is itself a direct input to the HDI composite index). Retained in this repository for transparency/completeness. |
| 6 | `who_unicef_hpv_vaccination_coverage.csv` | HPV vaccination coverage (% of 15-year-old girls, final recommended dose), by country and year | **WHO/UNICEF Estimates of National Immunization Coverage (WUENIC)** | OWID grapher mirror: https://ourworldindata.org/grapher/coverage-of-the-human-papillomavirus-vaccine.csv | ~126 countries with any reported coverage (of 175 in the base analytic sample) | Countries with no national HPV vaccination program, or that have not reported to WHO/UNICEF, have **no row** for this indicator — these were coded as 0% coverage in the analysis (see `scripts/04_merge_hpv_gni_screening.py`), not treated as missing data. |
| 7 | `who_cervical_screening_program.csv` | Binary indicator: does the country have an existing national cervical cancer screening program (Yes/No), by country and year | WHO Global Health Observatory | OWID grapher mirror: https://ourworldindata.org/grapher/countries-with-national-cervical-cancer-screening-program.csv | 194 countries | A **percentage of women screened** is not available as a clean, complete cross-country indicator (only Europe has that via Eurostat), so this binary "has a program" indicator was used instead as the operationalization of screening infrastructure. |
| 8 | `who_female_smoking_prevalence.csv` | Age-standardized female tobacco use prevalence (%), by country and year | WHO Global Health Observatory | OWID grapher mirror: https://ourworldindata.org/grapher/share-of-women-who-are-smoking.csv | 178 countries with recent data | Most recent available year per country was used (typically 2000-2020s range; exact year varies by country reporting). |
| 9 | `unaids_hiv_prevalence.csv` | HIV prevalence among adults aged 15-49 (%), by country and year | **UNAIDS** Global AIDS Update | OWID grapher mirror: https://ourworldindata.org/grapher/share-of-the-population-infected-with-hiv.csv | 153 countries with any reported data | Most recent available year per country was used. 25 countries in the 156-country full-covariate sample lacked an HIV estimate and were excluded from the HIV/mediation/spatial subsample (final n = 131). |
| 10 | `karlinsky_death_registration_completeness.csv` | Share of all deaths captured by a country's civil/vital registration system (%), 2015-2019 average, by country | Karlinsky, A. (2024). *International completeness of death registration 2015-2019.* Demographic Research, 50, 1151-1170 ("ICDR" dataset) | OWID grapher mirror: https://ourworldindata.org/grapher/share-of-deaths-registered.csv | 195 countries | Used **not** as a direct covariate but as a peer-reviewed proxy for underlying cancer-registry data quality, to test whether the core HDI-incidence relationship is an artifact of GLOBOCAN's imputation for data-poor countries (Section V-B / Fig. 6, right panel). Split at the standard 90% completeness threshold used in the demographic literature. |
| 11 | `country_centroids.json` | Latitude/longitude centroid coordinates for every ISO-3 country code | Public domain country-centroids reference dataset | Fetched from GitHub: https://raw.githubusercontent.com/komsitr/country-centroid/master/country-centroids.json | 245 countries/territories | Used only to build the k-nearest-neighbors spatial weights matrix for the Moran's I / spatial-lag analysis (Section V-D). Not used for any distance or routing calculation beyond adjacency. |

### Not included in `data/raw/` (too large / easily re-downloaded)

- **World country boundary polygons** (GeoJSON, ~14 MB), used only to render the world-map choropleth (`figures/fig03_world_map_incidence.png`). Source: https://raw.githubusercontent.com/datasets/geo-countries/master/data/countries.geojson — re-download this file directly and place it in `data/raw/countries.geojson` if you want to regenerate the world map from `scripts/figures/chart_world_map.py`.

---

## Processed Datasets (`data/processed/`)

Each file represents one stage of the analytic pipeline (see `scripts/` for the code that produces each one). All are plain CSV, one row per country.

| File | n (countries) | Produced by | Key columns added at this stage |
|---|---|---|---|
| `01_globocan2024_merged.csv` | 187 (incl. 1 world-aggregate row) | `scripts/01_merge_globocan_incidence_mortality.py` | `ASIR`, `ASMR`, `MIR` (= ASMR/ASIR), case/death counts |
| `02_with_hdi.csv` | 176 | `scripts/02_merge_hdi.py` | `HDI`, `HDI_tier` (Low/Medium/High/Very High), `HDI_year` |
| `03_full_covariates_n156.csv` | 156 | `scripts/04-07_*.py` | `HPV_vax_coverage`, `Screening_program_bin`, `GNI_per_capita`, `Smoking_prev_female` |
| `04_hiv_mediation_sample_n131.csv` | 131 | `scripts/10-11_*.py` | `HIV_prev` (adds HIV; drops countries missing it) |
| `05_registry_quality_sample_n142.csv` | 142 (subset of the 175-country base sample with registration data) | `scripts/12_registry_quality_stratification.py` | `Death_reg_pct`, `registry_quality` (High ≥90% / Lower <90%) |
| `06_panel_2022_2024.csv` | 331 rows (175 countries × 2022 cycle + 156 countries × 2024 cycle, long format) | `scripts/08_build_temporal_panel.py` | `Year`, time-matched `HDI_val` (2021 HDI for the 2022 cycle, 2023 HDI for the 2024 cycle) |
| `07_paired_country_changes.csv` | 156 (countries present in both GLOBOCAN cycles) | `scripts/09_analyze_temporal_panel.py` | Paired 2022 vs. 2024 ASIR, `pct_change`, `HDI_tier` |
| `08_spatial_analysis_sample.csv` | 131 | `scripts/14_moran_i_spatial_autocorrelation.py` | `latitude`, `longitude`, OLS `residuals` (for Moran's I) |

---

## Scripts (`scripts/`)

Run in numeric order to reproduce the full analysis from raw data to final results. Each script reads from `data/raw/` and/or the previous stage's output in `data/processed/`, and writes its own output CSV.

| Script | Purpose | Manuscript section |
|---|---|---|
| `01_merge_globocan_incidence_mortality.py` | Merge GLOBOCAN 2024 incidence + mortality, compute MIR | Section IV-A |
| `02_merge_hdi.py` | Filter HDI to latest year per country, drop aggregates, merge | Section III-A |
| `03_analyze_hdi_correlations.py` | Spearman correlations, HDI-tier medians, ratio calculations | Section IV-B |
| `04_merge_hpv_gni_screening.py` | Merge HPV vaccination, GNI, screening-program indicator | Section III-A |
| `05_multivariable_model_base.py` | First multivariable OLS (HDI + HPV + screening), VIF diagnostics | Section IV-D |
| `06_add_smoking_prevalence.py` | Merge female smoking prevalence | Section III-A |
| `07_multivariable_model_with_smoking.py` | Multivariable OLS incl. smoking; confounding-by-HDI check | Section IV-D |
| `08_build_temporal_panel.py` | Construct the paired 2022-2024 GLOBOCAN panel | Section V-A |
| `09_analyze_temporal_panel.py` | Wilcoxon signed-rank, Kruskal-Wallis, tier-ratio comparison | Section V-A |
| `10_add_hiv_and_registry_quality.py` | Merge HIV prevalence and death-registration completeness | Section III-A |
| `11_multivariable_model_with_hiv.py` | Full multivariable OLS including HIV; model R² comparison | Section V-B |
| `12_registry_quality_stratification.py` | Split by registration completeness; stratified correlations | Section V-B |
| `13_mediation_analysis.py` | Product-of-coefficients mediation (Sobel test + bootstrap CI) | Section V-C |
| `14_moran_i_spatial_autocorrelation.py` | Build k-NN weights, global Moran's I on ASIR and residuals | Section V-D |
| `15_spatial_lag_error_diagnostics.py` | OLS spatial diagnostics, Lagrange Multiplier tests | Section V-D |
| `16_spatial_lag_model_final.py` | Final maximum-likelihood spatial lag model | Section V-D |
| `17_spatial_k_sensitivity_check.py` | Re-run spatial model at k = 3, 5, 8, 10 for robustness | Section V-D |
| `figures/*.py` | One script per manuscript figure (charts + Graphviz flowcharts) | Figs. 1-10 |

---

## Notebook (`notebooks/`)

`cervical_cancer_prediction_model.ipynb` — a **supplementary, predictive-modeling companion** to the main inferential analysis. Compares Linear/Ridge regression against tuned Random Forest and Gradient Boosting models for predicting log(ASIR) from the same covariate set (best test R² = 0.77 with HIV included), with SHAP interpretability. This notebook is exploratory/methodological and is not part of the main manuscript's statistical inference — see the manuscript for the primary, hypothesis-driven analysis.

---

## Reproducing the Analysis

```bash
# 1. Clone this repository, then install dependencies
pip install -r requirements.txt

# 2. Run the pipeline in order (from the repository root)
for f in scripts/0*.py scripts/1*.py; do python3 "$f"; done

# 3. Regenerate figures
for f in scripts/figures/*.py; do python3 "$f"; done

# 4. (Optional) Run the ML notebook
jupyter nbconvert --to notebook --execute --inplace notebooks/cervical_cancer_prediction_model.ipynb
```

Note: several scripts contain hardcoded file paths from the original development environment (e.g. `/mnt/user-data/outputs/...`) — update these to relative paths (e.g. `../data/processed/...`) before running outside that environment. A fully path-relative refactor is a natural first contribution for anyone extending this repository.

---

## Key Findings Summary

- HDI is strongly, independently associated with cervical cancer incidence (Spearman ρ = -0.63) and mortality (ρ = -0.74), robust across every check applied.
- A naive comparison of GLOBOCAN 2022 vs. 2024 tier-level averages suggested the low-vs-high-HDI gap widened sharply — a **paired, same-country comparison shows this was a compositional artifact**, not a real trend (p = 0.940).
- HIV prevalence explains a modest, statistically significant 7.1% of HDI's total effect on mortality (mediation analysis).
- The HDI-incidence relationship holds equally in countries with high- and low-quality underlying surveillance data, arguing against a GLOBOCAN data-imputation artifact.
- Formal spatial-lag correction confirms HDI's effect is real but **roughly half the magnitude** conventional (spatially-naive) OLS estimates — a caution for the wider ecological-GLOBOCAN literature.

Full results, statistics, and discussion are in `manuscript/`.

---

## Data Licensing / Attribution

All raw data are public and free to reuse for research purposes; consult each original publisher's terms for redistribution or commercial use:
- GLOBOCAN/IARC data: subject to IARC's terms of use (https://gco.iarc.who.int)
- UNDP HDI/GNI data: UNDP Human Development Report Office terms
- WHO/UNICEF, WHO GHO data: WHO terms of use
- UNAIDS data: UNAIDS terms of use
- Our World in Data mirrors: OWID content is available under Creative Commons BY 4.0 (https://ourworldindata.org/about#access-to-owid-content-and-data); underlying data retain their original source's license
- Karlinsky (2024) ICDR dataset: released under the paper's stated open-data terms (Demographic Research is an open-access journal)

## Citation

If you use this repository, please cite the accompanying manuscript (see `manuscript/`) and the original data sources listed above.
