################################################################################
## Create supplementary tables workbooks with five tables per workbook
##
## Current analysis only: output
## This script does not use the pre-ComBat analysis outputs.
################################################################################

library(data.table)
library(dplyr)
library(openxlsx)

project_dir <- ""
analysis_output_dir <- file.path(project_dir, "output")
reference_dir <- file.path(project_dir, "references")
output_dir <- file.path(analysis_output_dir, "supplementary_tables")
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

output_file_1 <- file.path(
  output_dir,
  "supplementary_tables_S1-S5.xlsx"
)

output_file_2 <- file.path(
  output_dir,
  "supplementary_tables_S6-S10.xlsx"
)

output_file_3 <- file.path(
  output_dir,
  "supplementary_tables_S11-S15.xlsx"
)

output_file_4 <- file.path(
  output_dir,
  "supplementary_tables_S16-S20.xlsx"
)

################################################################################
## Input paths
################################################################################

paths <- list(
  table1=file.path(analysis_output_dir, "table1", "Table1_final.tsv"),
  pairwise=file.path(analysis_output_dir, "table1", "Table1_pairwise_comparisons.tsv"),

  s3=file.path(
    analysis_output_dir, "svyglm_anno",
    "unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv"
  ),
  s4=file.path(
    analysis_output_dir, "svyglm_anno",
    "unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv"
  ),
  s5=file.path(
    analysis_output_dir, "svyglm_anno",
    "unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv"
  ),
  s6=file.path(
    analysis_output_dir, "svyglm_anno",
    "unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv"
  ),
  s7=file.path(
    analysis_output_dir, "former_vs_current", "annotated",
    "unimputed_weighted_restricted_former_vs_current_euro_annotated.tsv"
  ),
  s8=file.path(
    analysis_output_dir, "former_vs_current", "annotated",
    "unimputed_weighted_restricted_former_vs_current_transancestry_annotated.tsv"
  ),

  epic_current=file.path(analysis_output_dir, "EWAS_ALHS_EU_ALL_smok1vs3_e.cov11_817235x1619.csv.s"),
  epic_former=file.path(analysis_output_dir, "EWAS_ALHS_EU_ALL_smok1vs2_e.cov11_2190_2025-04-23.csv"),

  cdmc=file.path(analysis_output_dir, "cdmc", "cdmc_complete_results.rdata"),

  age_primary=file.path(
    analysis_output_dir, "age_ecel1p1", "cg04284215_primary_EWAS_coefficient.tsv"
  ),
  age_summary=file.path(
    analysis_output_dir, "age_ecel1p1", "age_smoking_interaction_summary.tsv"
  ),
  age_contrasts=file.path(
    analysis_output_dir, "age_ecel1p1", "age_smoking_interaction_contrasts.tsv"
  ),
  neutrophil_interval=file.path(
    analysis_output_dir, "neutrophil_cessation", "neutrophil_by_cessation_interval.tsv"
  ),
  neutrophil_interval_test=file.path(
    analysis_output_dir, "neutrophil_cessation", "neutrophil_cessation_interval_test.tsv"
  ),
  neutrophil_spearman=file.path(
    analysis_output_dir, "neutrophil_cessation", "neutrophil_cessation_spearman.tsv"
  ),
  neutrophil_adjusted=file.path(
    analysis_output_dir, "neutrophil_cessation", "neutrophil_cessation_adjusted_models.tsv"
  ),

  permutation_summary=file.path(
    analysis_output_dir, "perm", "permutation_summary.tsv"
  ),

  reversion_summary=file.path(
    analysis_output_dir, "reversion", "reversion_classification_summary.csv"
  ),
  reversion_overlap=file.path(
    analysis_output_dir, "reversion", "weighted_residual_unique_cpg_overlap_summary.csv"
  ),
  reversion_euro=file.path(
    analysis_output_dir, "reversion", "model3_weighted_euro_reversion_long_annotated.csv"
  ),
  reversion_trans=file.path(
    analysis_output_dir, "reversion", "model12_weighted_trans_reversion_long_annotated.csv"
  ),

  msa_manifest=file.path(
    reference_dir, "MSA.hg38.manifest.gencode.v41.tsv.gz"
  ),
  epic_manifest=file.path(
    reference_dir, "EPIC.hg38.manifest.gencode.v36.tsv.gz"
  ),
  epicv2_manifest=file.path(
    reference_dir, "EPICv2.hg38.manifest.gencode.v41.tsv.gz"
  ),
  hm450_manifest=file.path(
    reference_dir, "HM450.hg38.manifest.tsv.gz"
  )
)

################################################################################
## General helpers
################################################################################

must_exist <- function(path, label=basename(path)) {
  if (!file.exists(path)) {
    stop("Missing required file for ", label, ":\n", path, call.=FALSE)
  }
  path
}

read_any <- function(path) {
  must_exist(path)
  as.data.frame(data.table::fread(path, data.table=FALSE, check.names=FALSE))
}

first_existing <- function(candidates, label) {
  hit <- candidates[file.exists(candidates)]
  if (length(hit)==0) {
    stop(
      "Could not find ", label, ". Checked:\n",
      paste(candidates, collapse="\n"),
      call.=FALSE
    )
  }
  hit[[1]]
}

get_column <- function(x, candidates, default=NA) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit)==0) return(rep(default, nrow(x)))
  x[[hit[[1]]]]
}

numeric_column <- function(x, candidates, default=NA_real_) {
  suppressWarnings(as.numeric(get_column(x, candidates, default)))
}

character_column <- function(x, candidates, default=NA_character_) {
  as.character(get_column(x, candidates, default))
}

base_cpg <- function(x) sub("_.*$", "", as.character(x))

safe_left_join <- function(x, y, by) {
  x %>% dplyr::left_join(y, by=by)
}

################################################################################
## Manifest annotation and array overlap
################################################################################

msa_raw <- read_any(paths$msa_manifest)
epic_raw <- read_any(paths$epic_manifest)
epicv2_raw <- read_any(paths$epicv2_manifest)
hm450_raw <- read_any(paths$hm450_manifest)

msa_annotation <- data.frame(
  Probe_ID=character_column(msa_raw, c("probeID", "Probe_ID", "ProbeID")),
  Chromosome=character_column(msa_raw, c("CpG_chrm", "Chromosome", "chr")),
  CpG_beginning=numeric_column(msa_raw, c("CpG_beg", "CpG_beginning", "start")),
  CpG_end=numeric_column(msa_raw, c("CpG_end", "end")),
  Probe_strand=character_column(msa_raw, c("probe_strand", "Probe_strand")),
  Gene_name_unique=character_column(msa_raw, c("genesUniq", "Gene_name_unique")),
  Gene_name=character_column(msa_raw, c("geneNames", "Gene_name")),
  Transcript_type=character_column(msa_raw, c("transcriptTypes", "Transcript_type")),
  Transcript_ID=character_column(msa_raw, c("transcriptIDs", "Transcript_ID")),
  Distance_to_transcription_start_site=numeric_column(
    msa_raw,
    c("distToTSS", "Distance_to_transcription_start_site")
  ),
  stringsAsFactors=FALSE
) %>%
  mutate(CpG=base_cpg(Probe_ID)) %>%
  distinct(Probe_ID, .keep_all=TRUE)

epic_annotation <- data.frame(
  Probe_ID=character_column(epic_raw, c("probeID", "Probe_ID", "ProbeID")),
  Chromosome=character_column(epic_raw, c("CpG_chrm", "Chromosome", "chr")),
  CpG_beginning=numeric_column(epic_raw, c("CpG_beg", "CpG_beginning", "start")),
  CpG_end=numeric_column(epic_raw, c("CpG_end", "end")),
  Probe_strand=character_column(epic_raw, c("probe_strand", "Probe_strand")),
  Gene_name_unique=character_column(epic_raw, c("genesUniq", "Gene_name_unique")),
  Gene_name=character_column(epic_raw, c("geneNames", "Gene_name")),
  Transcript_type=character_column(epic_raw, c("transcriptTypes", "Transcript_type")),
  Transcript_ID=character_column(epic_raw, c("transcriptIDs", "Transcript_ID")),
  Distance_to_transcription_start_site=numeric_column(
    epic_raw,
    c("distToTSS", "Distance_to_transcription_start_site")
  ),
  stringsAsFactors=FALSE
) %>%
  mutate(
    Probe_ID=base_cpg(Probe_ID),
    CpG=Probe_ID
  ) %>%
  distinct(Probe_ID, .keep_all=TRUE)

epic_ids <- unique(base_cpg(character_column(
  epic_raw, c("probeID", "Probe_ID", "ProbeID", "CpG")
)))
epicv2_ids <- unique(base_cpg(character_column(
  epicv2_raw, c("probeID", "Probe_ID", "ProbeID", "CpG")
)))
hm450_ids <- unique(base_cpg(character_column(
  hm450_raw, c("probeID", "Probe_ID", "ProbeID", "CpG")
)))

add_array_overlap <- function(x) {
  if (!"CpG" %in% names(x)) x$CpG <- base_cpg(x$Probe_ID)

  x %>%
    mutate(
      on_EPIC=CpG %in% epic_ids,
      on_EPICv2=CpG %in% epicv2_ids,
      on_450K=CpG %in% hm450_ids,
      array_overlap=case_when(
        on_EPIC & on_EPICv2 & on_450K ~ "EPIC_EPICv2_450K",
        on_EPIC & on_EPICv2 ~ "EPIC_EPICv2",
        on_EPIC & on_450K ~ "EPIC_450K",
        on_EPICv2 & on_450K ~ "EPICv2_450K",
        on_EPIC ~ "EPIC",
        on_EPICv2 ~ "EPICv2",
        on_450K ~ "450K",
        TRUE ~ "MSA_ONLY"
      )
    )
}

annotation_columns <- c(
  "Probe_ID", "CpG", "Beta", "Standard_error", "P_value", "N",
  "n_smoker", "n_never_smoker", "Chromosome", "CpG_beginning", "CpG_end",
  "Probe_strand", "Gene_name_unique", "Gene_name", "Transcript_type",
  "Transcript_ID", "Distance_to_transcription_start_site",
  "on_EPIC", "on_EPICv2", "on_450K", "array_overlap"
)

standardize_ewas <- function(
  x,
  exposed_label=c("smoker", "former", "current"),
  annotation=msa_annotation
) {
  exposed_label <- match.arg(exposed_label)

  probe <- character_column(x, c("Probe_ID", "V1", "CpG", "probeID"))
  output <- data.frame(
    Probe_ID=probe,
    CpG=base_cpg(probe),
    Beta=numeric_column(x, c("Beta", "Estimate", "beta")),
    Standard_error=numeric_column(x, c("SE", "Standard_error", "Std.Error")),
    P_value=numeric_column(x, c("P", "P_value", "p", "p.value")),
    N=numeric_column(x, c("N", "n")),
    n_exposed=numeric_column(x, c("N.exp", "n_smoker", "n_exposed")),
    n_reference=numeric_column(x, c("N.unexp", "n_never_smoker", "n_reference")),
    stringsAsFactors=FALSE
  )

  output <- output %>%
    safe_left_join(annotation, by=c("Probe_ID", "CpG")) %>%
    add_array_overlap()

  ## n_smoker is the exposed smoking group for the sheet: current smokers in
  ## current-versus-never models and former smokers in former-versus-never models.
  names(output)[names(output)=="n_exposed"] <- "n_smoker"
  names(output)[names(output)=="n_reference"] <- "n_never_smoker"

  first_columns <- c(
    "Probe_ID", "CpG", "Beta", "Standard_error", "P_value", "N",
    "n_smoker", "n_never_smoker"
  )

  output %>%
    select(
      all_of(first_columns),
      Chromosome, CpG_beginning, CpG_end, Probe_strand,
      Gene_name_unique, Gene_name, Transcript_type, Transcript_ID,
      Distance_to_transcription_start_site,
      on_EPIC, on_EPICv2, on_450K, array_overlap
    ) %>%
    arrange(P_value)
}

standardize_former_current <- function(x) {
  probe <- character_column(x, c("Probe_ID", "V1", "CpG", "probeID"))

  data.frame(
    Probe_ID=probe,
    CpG=base_cpg(probe),
    Beta=numeric_column(x, c("Beta", "Estimate", "beta")),
    Standard_error=numeric_column(x, c("SE", "Standard_error", "Std.Error")),
    P_value=numeric_column(x, c("P", "P_value", "p", "p.value")),
    N=numeric_column(x, c("N", "n")),
    n_former_smoker=numeric_column(x, c("N.exp", "n_former_smoker", "n_exposed")),
    n_current_smoker=numeric_column(x, c("N.unexp", "n_current_smoker", "n_reference")),
    stringsAsFactors=FALSE
  ) %>%
    safe_left_join(msa_annotation, by=c("Probe_ID", "CpG")) %>%
    add_array_overlap() %>%
    select(
      Probe_ID, CpG, Beta, Standard_error, P_value, N,
      n_former_smoker, n_current_smoker,
      Chromosome, CpG_beginning, CpG_end, Probe_strand,
      Gene_name_unique, Gene_name, Transcript_type, Transcript_ID,
      Distance_to_transcription_start_site,
      on_EPIC, on_EPICv2, on_450K, array_overlap
    ) %>%
    arrange(P_value)
}

standardize_annotated_detail <- function(x) {
  probe <- character_column(x, c("Probe_ID", "V1", "CpG", "probeID"))

  annotation_aliases <- c(
    "Probe_ID", "CpG", "CpG_ID", "V1", "probeID",
    "CpG_chrm", "Chromosome", "CpG_beg", "CpG_beginning", "CpG_end",
    "probe_strand", "Probe_strand", "genesUniq", "Gene_name_unique",
    "geneNames", "Gene_name", "transcriptTypes", "Transcript_type",
    "transcriptIDs", "Transcript_ID", "distToTSS",
    "Distance_to_transcription_start_site",
    "on_EPIC", "on_EPICv2", "on_450K", "array_overlap"
  )

  analysis_columns <- x %>% select(-any_of(annotation_aliases))

  data.frame(
    Probe_ID=probe,
    CpG=base_cpg(probe),
    stringsAsFactors=FALSE
  ) %>%
    bind_cols(analysis_columns) %>%
    safe_left_join(msa_annotation, by=c("Probe_ID", "CpG")) %>%
    add_array_overlap() %>%
    select(
      Probe_ID, CpG,
      Chromosome, CpG_beginning, CpG_end, Probe_strand,
      Gene_name_unique, Gene_name, Transcript_type, Transcript_ID,
      Distance_to_transcription_start_site,
      on_EPIC, on_EPICv2, on_450K, array_overlap,
      everything()
    )
}

################################################################################
## S1-S10: participant and EWAS tables
################################################################################

s1 <- read_any(paths$table1)
s2 <- read_any(paths$pairwise)

s3 <- standardize_ewas(read_any(paths$s3), exposed_label="current")
s4 <- standardize_ewas(read_any(paths$s4), exposed_label="current")
s5 <- standardize_ewas(read_any(paths$s5), exposed_label="former")
s6 <- standardize_ewas(read_any(paths$s6), exposed_label="former")
s7 <- standardize_former_current(read_any(paths$s7))
s8 <- standardize_former_current(read_any(paths$s8))

## EPIC replication files use the base CpG identifier as the probe identifier.
s9 <- standardize_ewas(
  read_any(paths$epic_current),
  exposed_label="current",
  annotation=epic_annotation
)
s10 <- standardize_ewas(
  read_any(paths$epic_former),
  exposed_label="former",
  annotation=epic_annotation
)

################################################################################
## S11-S17: CellDMC results
################################################################################

cdmc_environment <- new.env()
load(must_exist(paths$cdmc, "CellDMC results"), envir=cdmc_environment)

required_cdmc_objects <- c("cdmc_nf", "cdmc_nc", "cdmc_fc", "cdmc_trend")
missing_cdmc_objects <- required_cdmc_objects[
  !vapply(required_cdmc_objects, exists, logical(1), envir=cdmc_environment)
]
if (length(missing_cdmc_objects)>0) {
  stop(
    "The CellDMC RData file is missing: ",
    paste(missing_cdmc_objects, collapse=", "),
    call.=FALSE
  )
}

cell_types <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")

cdmc_contrasts <- list(
  list(
    Comparison="Never vs Former",
    object=get("cdmc_nf", envir=cdmc_environment),
    N=651, Group_1="Former", n_group_1=202,
    Group_2="Never", n_group_2=449
  ),
  list(
    Comparison="Never vs Current",
    object=get("cdmc_nc", envir=cdmc_environment),
    N=515, Group_1="Current", n_group_1=66,
    Group_2="Never", n_group_2=449
  ),
  list(
    Comparison="Former vs Current",
    object=get("cdmc_fc", envir=cdmc_environment),
    N=268, Group_1="Current", n_group_1=66,
    Group_2="Former", n_group_2=202
  ),
  list(
    Comparison="Ordinal trend",
    object=get("cdmc_trend", envir=cdmc_environment),
    N=717, Group_1="Ordinal smoking status", n_group_1=717,
    Group_2=NA_character_, n_group_2=NA_real_
  )
)

extract_cdmc_cell <- function(contrast, cell) {
  result_object <- contrast$object

  if (is.null(result_object$coe[[cell]])) {
    stop("CellDMC coefficient output is missing cell type: ", cell, call.=FALSE)
  }

  coefficients <- as.data.frame(result_object$coe[[cell]], check.names=FALSE)
  probe <- rownames(coefficients)
  if (is.null(probe)) {
    probe <- character_column(coefficients, c("Probe_ID", "CpG", "V1"))
  }

  direction <- rep(0, length(probe))
  if (!is.null(result_object$dmct) && cell %in% colnames(result_object$dmct)) {
    direction <- as.numeric(result_object$dmct[probe, cell])
    direction[is.na(direction)] <- 0
  }

  data.frame(
    Comparison=contrast$Comparison,
    Probe_ID=probe,
    CpG=base_cpg(probe),
    Cell_type=cell,
    Direction=direction,
    Beta=numeric_column(
      coefficients,
      c("Estimate", "estimate", "Beta", "beta", "coef", "Coefficient")
    ),
    Standard_error=numeric_column(
      coefficients,
      c("SE", "Std. Error", "Std.Error", "Standard_error", "se")
    ),
    T_statistic=numeric_column(
      coefficients,
      c("T-statistic", "T.statistic", "t", "tstat", "Statistic")
    ),
    P_value=numeric_column(
      coefficients,
      c("p", "P", "P_value", "p.value", "Pr(>|t|)")
    ),
    FDR=numeric_column(
      coefficients,
      c("adjP", "FDR", "fdr", "adjusted_p", "adj.P.Val")
    ),
    N=contrast$N,
    Group_1=contrast$Group_1,
    n_group_1=contrast$n_group_1,
    Group_2=contrast$Group_2,
    n_group_2=contrast$n_group_2,
    stringsAsFactors=FALSE
  ) %>%
    safe_left_join(msa_annotation, by=c("Probe_ID", "CpG")) %>%
    add_array_overlap() %>%
    select(
      Comparison, Probe_ID, CpG, Cell_type, Direction,
      Beta, Standard_error, T_statistic, P_value, FDR,
      N, Group_1, n_group_1, Group_2, n_group_2,
      Chromosome, CpG_beginning, CpG_end, Probe_strand,
      Gene_name_unique, Gene_name, Transcript_type, Transcript_ID,
      Distance_to_transcription_start_site,
      on_EPIC, on_EPICv2, on_450K, array_overlap
    )
}

cell_results <- setNames(vector("list", length(cell_types)), cell_types)
for (cell in cell_types) {
  cell_results[[cell]] <- bind_rows(lapply(
    cdmc_contrasts,
    extract_cdmc_cell,
    cell=cell
  ))
}

## S11 contains every significant probe-design-level DMCT across all contrasts.
## S12-S17 contain the complete CellDMC coefficient output for each cell type.
s11 <- bind_rows(cell_results) %>%
  filter(Direction != 0) %>%
  arrange(Comparison, Cell_type, P_value)

s11_summary_observed <- s11 %>%
  group_by(Comparison, Cell_type, Direction) %>%
  summarise(
    Significant_probe_designs=n(),
    Significant_unique_CpGs=n_distinct(CpG),
    .groups="drop"
  )

s11_summary <- expand.grid(
  Comparison=vapply(cdmc_contrasts, `[[`, character(1), "Comparison"),
  Cell_type=cell_types,
  Direction=c(-1, 1),
  stringsAsFactors=FALSE
) %>%
  left_join(
    s11_summary_observed,
    by=c("Comparison", "Cell_type", "Direction")
  ) %>%
  mutate(
    Significant_probe_designs=ifelse(
      is.na(Significant_probe_designs), 0L, Significant_probe_designs
    ),
    Significant_unique_CpGs=ifelse(
      is.na(Significant_unique_CpGs), 0L, Significant_unique_CpGs
    )
  ) %>%
  arrange(Comparison, Cell_type, Direction)

s12 <- cell_results$B
s13 <- cell_results$NK
s14 <- cell_results$CD4T
s15 <- cell_results$CD8T
s16 <- cell_results$Mono
s17 <- cell_results$Neutro

################################################################################
## S18: cessation-time, cg05575921, neutrophil, and age analyses
################################################################################

cessation_distribution <- data.frame(
  Population=c(
    rep("Transancestry", 4),
    rep("European", 4)
  ),
  Cessation_interval=rep(
    c("<=10 years", "11-20 years", "21-30 years", ">30 years"),
    2
  ),
  N=c(33, 84, 30, 55, 29, 70, 27, 49),
  Percent=c(16.3, 41.6, 14.9, 27.2, 16.6, 40.0, 15.4, 28.0),
  stringsAsFactors=FALSE
)

cessation_overall <- data.frame(
  Population="Transancestry former smokers",
  N=202,
  N_with_cessation_time=202,
  Missing=0,
  Median_years=18,
  Q1=14.25,
  Q3=32,
  Minimum=0,
  Maximum=66
)

cg05575921_cessation <- data.frame(
  Population="Transancestry former smokers",
  Probe_ID=c(
    "cg05575921_TC11", "cg05575921_TC12",
    "cg05575921_TC13", "cg05575921_TC14"
  ),
  CpG="cg05575921",
  N=202,
  Beta_per_year=c(0.001773, 0.001717, 0.001698, 0.001712),
  Standard_error=c(0.000328, 0.000315, 0.000347, 0.000329),
  Beta_per_10_years=c(0.01773, 0.01717, 0.01698, 0.01712),
  Percentage_points_per_10_years=c(1.773, 1.717, 1.698, 1.712),
  P_value=c(2.01e-7, 1.63e-7, 2.21e-6, 5.52e-7),
  stringsAsFactors=FALSE
) %>%
  safe_left_join(msa_annotation, by=c("Probe_ID", "CpG")) %>%
  add_array_overlap()

neutrophil_interval <- read_any(paths$neutrophil_interval)
neutrophil_interval_test <- read_any(paths$neutrophil_interval_test)
neutrophil_spearman <- read_any(paths$neutrophil_spearman)
neutrophil_adjusted <- read_any(paths$neutrophil_adjusted)

age_analysis <- read_any(paths$age_summary)
age_contrasts <- read_any(paths$age_contrasts)

primary_cg04284215_raw <- read_any(paths$age_primary)
primary_cg04284215 <- primary_cg04284215_raw %>%
  transmute(
    Population,
    Probe_ID=as.character(Probe_ID),
    CpG=base_cpg(Probe_ID),
    Beta=as.numeric(Beta),
    Standard_error=as.numeric(SE),
    P_value=as.numeric(P),
    Difference_percentage_points,
    Absolute_difference_percentage_points,
    Greater_than_10_percentage_points
  ) %>%
  safe_left_join(msa_annotation, by=c("Probe_ID", "CpG")) %>%
  add_array_overlap()

################################################################################
## S19: reversion analyses
################################################################################

reversion_summary <- read_any(paths$reversion_summary)
reversion_overlap <- read_any(paths$reversion_overlap)

reversion_euro <- standardize_annotated_detail(read_any(paths$reversion_euro)) %>%
  mutate(Population="Weighted European", .before=1)
reversion_trans <- standardize_annotated_detail(read_any(paths$reversion_trans)) %>%
  mutate(Population="Weighted transancestry", .before=1)
weighted_reversion_detail <- bind_rows(reversion_euro, reversion_trans)

################################################################################
## S20: raw genomic inflation and permutation analyses
################################################################################

lambda_gc <- function(beta, standard_error) {
  z <- as.numeric(beta) / as.numeric(standard_error)
  z <- z[is.finite(z)]
  median(z^2, na.rm=TRUE) / qchisq(0.5, df=1)
}

raw_lambda <- data.frame(
  Model=c(
    "Current vs never, European",
    "Former vs never, European",
    "Current vs never, transancestry",
    "Former vs never, transancestry"
  ),
  Lambda_GC=c(
    lambda_gc(s3$Beta, s3$Standard_error),
    lambda_gc(s5$Beta, s5$Standard_error),
    lambda_gc(s4$Beta, s4$Standard_error),
    lambda_gc(s6$Beta, s6$Standard_error)
  )
)

euro_permutation_file <- first_existing(
  c(
    file.path(analysis_output_dir, "perm", "euro", "permutation_lambdas.tsv"),
    file.path(analysis_output_dir, "perm", "euro", "permutation_lambdas_100.tsv")
  ),
  "European permutation lambda output"
)

trans_permutation_file <- first_existing(
  c(
    file.path(analysis_output_dir, "perm", "trans", "permutation_lambdas.tsv"),
    file.path(analysis_output_dir, "perm", "trans", "permutation_lambdas_100.tsv")
  ),
  "transancestry permutation lambda output"
)

read_permutation <- function(path, population) {
  x <- read_any(path)
  lambda <- numeric_column(
    x,
    c("Lambda_GC", "lambda", "Lambda", "lambda_gc", "V1")
  )
  data.frame(
    Population=population,
    Permutation=seq_along(lambda),
    Lambda_GC=lambda
  ) %>%
    filter(is.finite(Lambda_GC))
}

permutation_values <- bind_rows(
  read_permutation(euro_permutation_file, "European"),
  read_permutation(trans_permutation_file, "Transancestry")
)

permutation_summary <- read_any(paths$permutation_summary)

################################################################################
## Workbook-writing helpers
##
## Keep Excel's default appearance: standard gridlines, no filters/dropdowns,
## no frozen panes, and no color formatting. Sheet titles are bold. Column
## headers are bold with a thick bottom border. Columns are moderately widened.
################################################################################

wb <- createWorkbook(creator="ALHS smoking EWAS")
title_style <- createStyle(textDecoration="bold")
header_style <- createStyle(
  textDecoration="bold",
  border="bottom",
  borderStyle="thick",
  borderColour="#000000"
)

set_readable_widths <- function(sheet, number_of_columns) {
  if (number_of_columns > 0) {
    setColWidths(
      wb,
      sheet,
      cols=seq_len(number_of_columns),
      widths=18
    )
  }
}

write_single_table <- function(sheet, title, x, note=NULL) {
  if (nrow(x) + 4 > 1048576) {
    stop(sheet, " exceeds the Excel row limit.", call.=FALSE)
  }

  addWorksheet(wb, sheet)
  writeData(wb, sheet, title, startRow=1, startCol=1, colNames=FALSE)
  addStyle(wb, sheet, title_style, rows=1, cols=1, stack=TRUE)

  header_row <- 3
  if (!is.null(note)) {
    writeData(wb, sheet, note, startRow=2, startCol=1, colNames=FALSE)
    header_row <- 4
  }

  writeData(
    wb, sheet, x,
    startRow=header_row, startCol=1,
    colNames=TRUE, rowNames=FALSE, withFilter=FALSE
  )
  addStyle(
    wb, sheet, header_style,
    rows=header_row, cols=seq_len(ncol(x)),
    gridExpand=TRUE, stack=TRUE
  )
  set_readable_widths(sheet, ncol(x))
}

write_sections <- function(sheet, title, sections) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, title, startRow=1, startCol=1, colNames=FALSE)
  addStyle(wb, sheet, title_style, rows=1, cols=1, stack=TRUE)

  current_row <- 3
  maximum_columns <- 0
  for (section in sections) {
    section_title <- section$title
    x <- section$data

    writeData(
      wb, sheet, section_title,
      startRow=current_row, startCol=1, colNames=FALSE
    )
    current_row <- current_row + 1

    writeData(
      wb, sheet, x,
      startRow=current_row, startCol=1,
      colNames=TRUE, rowNames=FALSE, withFilter=FALSE
    )
    addStyle(
      wb, sheet, header_style,
      rows=current_row, cols=seq_len(ncol(x)),
      gridExpand=TRUE, stack=TRUE
    )
    maximum_columns <- max(maximum_columns, ncol(x))
    current_row <- current_row + nrow(x) + 2
  }
  set_readable_widths(sheet, maximum_columns)
}

################################################################################
## Table of contents and final workbook
################################################################################

sheet_names <- c(
  "Contents",
  "S1 Table 1",
  "S2 Pairwise comparisons",
  "S3 EUR current EWAS",
  "S4 TRANS current EWAS",
  "S5 EUR former EWAS",
  "S6 TRANS former EWAS",
  "S7 EUR former-current",
  "S8 TRANS former-current",
  "S9 EPIC current replication",
  "S10 EPIC former replication",
  "S11 Significant DMCTs",
  "S12 CellDMC B cells",
  "S13 CellDMC NK cells",
  "S14 CellDMC CD4 T cells",
  "S15 CellDMC CD8 T cells",
  "S16 CellDMC monocytes",
  "S17 CellDMC neutrophils",
  "S18 Cessation and age",
  "S19 Reversion analyses",
  "S20 Inflation permutations"
)

contents <- data.frame(
  Table=c("S1", "S2", paste0("S", 3:20)),
  Sheet=sheet_names[-1],
  Description=c(
    "Participant characteristics overall and by smoking status",
    "Overall and pairwise comparisons of participant characteristics",
    "Unimputed weighted European current-versus-never EWAS with array overlap",
    "Unimputed weighted transancestry current-versus-never EWAS with array overlap",
    "Unimputed weighted European former-versus-never EWAS with array overlap",
    "Unimputed weighted transancestry former-versus-never EWAS with array overlap",
    "Unimputed weighted European former-versus-current EWAS with array overlap",
    "Unimputed weighted transancestry former-versus-current EWAS with array overlap",
    "EPIC European current-versus-never replication EWAS",
    "EPIC European former-versus-never replication EWAS",
    "Significant probe-design-level DMCTs across all CellDMC contrasts",
    "Complete B-cell CellDMC output across all contrasts",
    "Complete NK-cell CellDMC output across all contrasts",
    "Complete CD4 T-cell CellDMC output across all contrasts",
    "Complete CD8 T-cell CellDMC output across all contrasts",
    "Complete monocyte CellDMC output across all contrasts",
    "Complete neutrophil CellDMC output across all contrasts",
    "Cessation-time, cg05575921, neutrophil, and age interaction analyses",
    "Smoking-cessation reversion analyses",
    "Raw genomic inflation and smoking-status permutation analyses"
  ),
  stringsAsFactors=FALSE
)

start_workbook <- function(first_table, last_table) {
  wb <<- createWorkbook(creator="ALHS smoking EWAS")
  write_single_table(
    sheet_names[1],
    paste0(
      "Supplementary tables S", first_table, "-S", last_table,
      ": table of contents"
    ),
    contents %>%
      filter(Table %in% paste0("S", first_table:last_table))
  )
}

################################################################################
## Workbook 1: S1-S5
################################################################################

start_workbook(1, 5)

write_single_table(
  sheet_names[2],
  "Table S1. Participant characteristics by smoking status",
  s1
)
write_single_table(
  sheet_names[3],
  "Table S2. Overall and pairwise comparisons by smoking status",
  s2
)
write_single_table(
  sheet_names[4],
  "Table S3. Unimputed weighted European current-versus-never EWAS",
  s3
)
write_single_table(
  sheet_names[5],
  "Table S4. Unimputed weighted transancestry current-versus-never EWAS",
  s4
)
write_single_table(
  sheet_names[6],
  "Table S5. Unimputed weighted European former-versus-never EWAS",
  s5
)

saveWorkbook(wb, output_file_1, overwrite=TRUE)

################################################################################
## Workbook 2: S6-S10
################################################################################

start_workbook(6, 10)

write_single_table(
  sheet_names[7],
  "Table S6. Unimputed weighted transancestry former-versus-never EWAS",
  s6
)
write_single_table(
  sheet_names[8],
  "Table S7. Unimputed weighted European former-versus-current EWAS",
  s7,
  note="Positive beta values indicate higher methylation in former smokers relative to current smokers."
)
write_single_table(
  sheet_names[9],
  "Table S8. Unimputed weighted transancestry former-versus-current EWAS",
  s8,
  note="Positive beta values indicate higher methylation in former smokers relative to current smokers."
)
write_single_table(
  sheet_names[10],
  "Table S9. EPIC European current-versus-never replication EWAS",
  s9
)
write_single_table(
  sheet_names[11],
  "Table S10. EPIC European former-versus-never replication EWAS",
  s10
)

saveWorkbook(wb, output_file_2, overwrite=TRUE)

################################################################################
## Workbook 3: S11-S15
################################################################################

start_workbook(11, 15)

write_sections(
  sheet_names[12],
  "Table S11. Significant probe-design-level DMCTs across CellDMC contrasts",
  list(
    list(
      title=paste(
        "DMCT counts by contrast, cell type, and direction; significance was defined",
        "using the CellDMC Benjamini-Hochberg FDR threshold of 0.05"
      ),
      data=s11_summary
    ),
    list(
      title=paste(
        "Significant DMCT details; multiple MSA probe designs can target the same base CpG"
      ),
      data=s11
    )
  )
)
write_single_table(
  sheet_names[13],
  "Table S12. Complete CellDMC results for B cells",
  s12
)
write_single_table(
  sheet_names[14],
  "Table S13. Complete CellDMC results for NK cells",
  s13
)
write_single_table(
  sheet_names[15],
  "Table S14. Complete CellDMC results for CD4 T cells",
  s14
)
write_single_table(
  sheet_names[16],
  "Table S15. Complete CellDMC results for CD8 T cells",
  s15
)

saveWorkbook(wb, output_file_3, overwrite=TRUE)

################################################################################
## Workbook 4: S16-S20
################################################################################

start_workbook(16, 20)

write_single_table(
  sheet_names[17],
  "Table S16. Complete CellDMC results for monocytes",
  s16
)
write_single_table(
  sheet_names[18],
  "Table S17. Complete CellDMC results for neutrophils",
  s17
)

write_sections(
  sheet_names[19],
  "Table S18. Smoking cessation, cg05575921, neutrophil, and age analyses",
  list(
    list(
      title="Overall distribution of years since cessation",
      data=cessation_overall
    ),
    list(
      title="Cessation-interval distribution",
      data=cessation_distribution
    ),
    list(
      title="cg05575921 methylation versus years since cessation among former smokers",
      data=cg05575921_cessation
    ),
    list(
      title="Neutrophil proportions by cessation interval",
      data=neutrophil_interval
    ),
    list(
      title="Kruskal-Wallis test of neutrophil proportions across cessation intervals",
      data=neutrophil_interval_test
    ),
    list(
      title="Spearman correlation of neutrophil proportion with years since cessation",
      data=neutrophil_spearman
    ),
    list(
      title="Adjusted association of neutrophil proportion with years since cessation",
      data=neutrophil_adjusted
    ),
    list(
      title="Primary adjusted cg04284215 current-versus-never estimates",
      data=primary_cg04284215
    ),
    list(
      title="Age and smoking-status interaction analyses",
      data=age_analysis
    ),
    list(
      title="Age-model smoking-status contrasts",
      data=age_contrasts
    )
  )
)

write_sections(
  sheet_names[20],
  "Table S19. Smoking-cessation reversion analyses",
  list(
    list(
      title="Reversion-classification summary",
      data=reversion_summary
    ),
    list(
      title="Weighted residual unique-CpG overlap summary",
      data=reversion_overlap
    ),
    list(
      title="Weighted European and transancestry probe-design-level results",
      data=weighted_reversion_detail
    )
  )
)

write_sections(
  sheet_names[21],
  "Table S20. Raw genomic inflation and permutation analyses",
  list(
    list(
      title="Raw genomic inflation factors",
      data=raw_lambda
    ),
    list(
      title="Permutation summary",
      data=permutation_summary
    ),
    list(
      title="Individual permutation genomic inflation factors",
      data=permutation_values
    )
  )
)

saveWorkbook(wb, output_file_4, overwrite=TRUE)

message("Supplementary workbooks written to:")
message(output_file_1)
message(output_file_2)
message(output_file_3)
message(output_file_4)

