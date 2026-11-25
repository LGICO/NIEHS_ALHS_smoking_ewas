# NIEHS ALHS smoking EWAS

Scripts for quality control and analysis for *Smoking drives an epigenetic memory of aberrant hematopoiesis*: [medRxiv preprint](https://www.medrxiv.org/content/10.64898/2026.05.14.26353250v1).

Run the scripts from the repository root. The cleaned analytical dataset, `751_clean_pc.rdata`, must be located in the repository root. Input IDATs, reference files, and sample sheets belong in `input_idats/`, `references/`, and `samplesheets/`, respectively. Scripts create their required output directories when they are run.

## Scripts

- `00_sesame_clean_probe_filt.R` — Processes MSA IDATs with SeSAMe, applies sample and probe quality control, performs batch correction, and creates `751_clean_pc.rdata`.

- `01_celldmc.R` — Runs CellDMC smoking contrasts and ordinal-trend analyses at the probe-design level and after mean or median CpG collapsing; also creates CellDMC tables and plots.

- `02_main_models_final.R` — Runs the primary weighted and unweighted current- and former-smoking EWAS models in the European and transancestry populations, then annotates and filters the results.

- `03_sensitivity_models.R` — Runs EWAS sensitivity models, including case-control, sex-stratified, state-stratified, and interaction analyses, and writes the raw results under `output/model_output/`.

- `04_annotate_filter_output.R` — Annotates and significance-filters the sensitivity-model results in `output/model_output/` and creates selected methylation boxplots.

- `05_manhattan.R` — Creates Miami and stacked Manhattan plots, compares MSA results with earlier arrays and replication results, and identifies MSA-specific significant CpGs.

- `06_fig1_numbers.R` — Summarizes the numbers of significant probe designs, base CpGs, and MSA-specific CpGs within and across the main current- and former-smoking models.

- `07_plotDMR.R` — Generates ECEL1P1 regional methylation and DMR plots, including participant-level values and CpG-level summaries.

- `08_generate_boxplots.R` — Creates smoking-status plots for significant CpG probe designs from the main EWAS models using the revised and original color schemes.

- `09_former_smoker_stratified_analysis.R` — Runs former-smoker analyses using categorical cessation intervals, five- and ten-year groupings, and continuous years since cessation; also creates exploratory AHRR and ECEL1P1 summaries.

- `10_former_strat_loess_scatter.R` — Creates continuous and grouped AHRR and ECEL1P1 methylation trajectories by years since smoking cessation.

- `11_table1.R` — Creates participant-characteristic tables and overall and pairwise statistical comparisons by smoking status.

- `12_former_vs_current_ewas_restricted.R` — Runs the direct genome-wide former-versus-current EWAS in weighted and unweighted European and transancestry models, then annotates and filters the output.

- `13_validate_dmr.R` — Evaluates age, smoking status, and age-by-smoking interactions at cg04284215 and across the ECEL1P1 regional mean and creates the age scatterplots and summary tables.

- `14_neutrophil_cessation.R` — Tests associations between neutrophil proportions and smoking-cessation time and creates the corresponding tables and figure.

- `15_former_reversion.R` — Performs the targeted smoking-cessation reversion analysis, including attenuation estimates, former-versus-current and former-versus-never contrasts, and equivalence testing.

- `16_mean_median_collapsed_models.R` — Re-runs the four primary weighted EWAS models after collapsing multiple MSA probe designs to base CpGs using either the mean or median.

- `17_mean_median_manhattan.R` — Creates stacked Manhattan plots from the mean- and median-collapsed CpG analyses.

- `18_genomic_inflation.R` — Calculates genomic-inflation statistics, performs BACON diagnostic analyses, and runs the European and transancestry permutation analyses.

- `19_create_supplementary_tables.R` — Compiles analysis outputs into the supplementary-table Excel workbooks, with five numbered tables per workbook.

- `20_extract_msa_idat_calls.R` — Command-line utility that extracts MSA SNP genotypes and inferred sex from paired IDAT files for sample-identity quality control.

- `21_compare_msa_axiom.py` — Cross-matches MSA SNP calls against Axiom genotypes, evaluates expected and best genetic matches, and identifies possible sample mismatches or swaps.

- `22_check_snps.txt` — Example `bcftools`, R, and Python commands used to inspect the genotype VCF and run the MSA-Axiom identity workflow.

- `99_basic_stats.R` — Performs quick descriptive checks of cohort counts, smoking variables, and participant characteristics.

- `99_cpg_overlap.R` — Performs supporting CpG overlap comparisons among MSA, EPIC, and the principal smoking models.

- `99_former_strat_loess_sensitivity.R` — Contains supporting and legacy former-smoker, CellDMC, and ECEL1P1 sensitivity analyses and plots.

## Directory structure

```text
NIEHS_ALHS_smoking_ewas/
├── 751_clean_pc.rdata
├── README.md
├── input_idats/
├── references/
├── samplesheets/
├── scripts/
│   ├── 00_sesame_clean_probe_filt.R
│   ├── 01_celldmc.R
│   ├── ...
│   └── 99_former_strat_loess_sensitivity.R
└── output/
    ├── age_ecel1p1/
    ├── boxplots/
    ├── case_control_sensitivity/
    ├── cdmc/
    ├── cdmc_mean/
    ├── cdmc_median/
    ├── continuous/
    ├── cpg_collapsed_mean/
    ├── cpg_collapsed_median/
    ├── ecel1p1/
    ├── five_year/
    ├── former_stratified/
    ├── former_vs_current/
    ├── genomic_inflation/
    │   ├── full/
    │   └── significant/
    ├── merged_ewas_comparison/
    ├── miami/
    ├── model_output/
    │   ├── svyglm/
    │   ├── svyglm_anno/
    │   └── svyglm_significant_cpgs/
    ├── neutrophil_cessation/
    ├── perm/
    ├── reversion/
    ├── significant_models/
    ├── supplementary_figures/
    ├── supplementary_tables/
    ├── svyglm/
    ├── svyglm_anno/
    ├── table1/
    └── ten_year/
```
