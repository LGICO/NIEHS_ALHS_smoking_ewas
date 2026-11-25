library(dplyr)
library(data.table)
library(survey)
library(parallel)

project_dir <- ""
input_file <- "751_clean_pc.rdata"
output_dir <- file.path(
  project_dir,
  "output",
  "former_vs_current"
)

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

load(input_file)

## Use the unimputed methylation values, matching the primary reported analyses.
betap <- t(betas_miss_masked)

## Recode smoking status using Current as the reference category.
smoke_raw <- as.character(pheno$LH_smoke_status)
pheno$LH_smoke_status <- dplyr::case_when(
  smoke_raw %in% c("Never", "0") ~ "Never",
  smoke_raw %in% c("Former", "1") ~ "Former",
  smoke_raw %in% c("Current", "2") ~ "Current",
  smoke_raw %in% c("M", "Missing", "Unknown", "99") ~ "Missing",
  TRUE ~ NA_character_
)

## Retain the same coding used in the primary EWAS models.
anyets_raw <- as.character(pheno$ANYETS)
pheno$ANYETS <- dplyr::case_when(
  anyets_raw %in% c("No", "0") ~ "No",
  anyets_raw %in% c("Yes", "1") ~ "Yes",
  anyets_raw %in% c("M", "Missing", "Unknown", "99") ~ "Missing",
  TRUE ~ NA_character_
)
pheno$ANYETS <- factor(pheno$ANYETS, levels=c("No", "Yes", "Missing"))

screen_raw <- as.character(pheno$LHSSCREEN_12142015)
pheno$LHSSCREEN_12142015 <- dplyr::case_when(
  screen_raw %in% c("Control", "0") ~ "Control",
  screen_raw %in% c("Case", "1") ~ "Case",
  TRUE ~ NA_character_
)
pheno$LHSSCREEN_12142015 <- factor(
  pheno$LHSSCREEN_12142015,
  levels=c("Control", "Case")
)

gender_raw <- tolower(as.character(pheno$Gender))
pheno$Gender <- dplyr::case_when(
  gender_raw %in% c("male", "0") ~ "male",
  gender_raw %in% c("female", "1") ~ "female",
  TRUE ~ NA_character_
)
pheno$Gender <- factor(pheno$Gender, levels=c("male", "female"))

state_raw <- as.character(pheno$STATE)
pheno$STATE <- dplyr::case_when(
  state_raw %in% c("IA", "0") ~ "IA",
  state_raw %in% c("NC", "1") ~ "NC",
  TRUE ~ NA_character_
)
pheno$STATE <- factor(pheno$STATE, levels=c("IA", "NC"))

european_raw <- as.character(pheno$European)
pheno$European <- dplyr::case_when(
  european_raw %in% c("European", "1") ~ "European",
  european_raw %in% c("Non_European", "Non-European", "0") ~ "Non_European",
  european_raw %in% c("GWAS sample failure", "99") ~ "GWAS_sample_failure",
  TRUE ~ NA_character_
)
pheno$European <- factor(
  pheno$European,
  levels=c("Non_European", "European", "GWAS_sample_failure")
)

## Match methylation rows explicitly to phenotype order.
stopifnot(
  all(pheno$Sample_ID %in% rownames(betap)),
  all(c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino") %in% names(pheno))
)
betap <- betap[match(pheno$Sample_ID, rownames(betap)), , drop=FALSE]
stopifnot(identical(rownames(betap), pheno$Sample_ID))

## Calculate the same asthma sampling weights before restricting smoking groups.
p_asthma_overall <- 0.06852481

add_asthma_weights <- function(dat) {
  n_asthma <- sum(dat$LHSSCREEN_12142015 == "Case", na.rm=TRUE)
  n_nonasthma <- sum(dat$LHSSCREEN_12142015 == "Control", na.rm=TRUE)
  p_asthma_nested <- n_asthma/(n_asthma+n_nonasthma)
  p_nonasthma_overall <- 1-p_asthma_overall
  p_nonasthma_nested <- 1-p_asthma_nested
  w_asthma <- 1/(p_asthma_nested/p_asthma_overall)
  w_nonasthma <- 1/(p_nonasthma_nested/p_nonasthma_overall)
  dat$asthma_weights <- ifelse(
    dat$LHSSCREEN_12142015 == "Case",
    w_asthma,
    w_nonasthma
  )
  dat
}

pheno_pc_all <- pheno %>%
  dplyr::filter(European == "European", !is.na(PC1)) %>%
  droplevels()

pheno <- add_asthma_weights(pheno)
pheno_pc_all <- add_asthma_weights(pheno_pc_all)

## Restrict the EWAS to former and current smokers only.
## Current is the reference, so LH_smoke_statusFormer is Former minus Current.
pheno_fc <- pheno %>%
  dplyr::filter(LH_smoke_status %in% c("Current", "Former")) %>%
  dplyr::mutate(
    LH_smoke_status=factor(
      LH_smoke_status,
      levels=c("Current", "Former")
    )
  ) %>%
  droplevels()

pheno_pc_fc <- pheno_pc_all %>%
  dplyr::filter(LH_smoke_status %in% c("Current", "Former")) %>%
  dplyr::mutate(
    LH_smoke_status=factor(
      LH_smoke_status,
      levels=c("Current", "Former")
    )
  ) %>%
  droplevels()

betap_fc <- betap[match(pheno_fc$Sample_ID, rownames(betap)), , drop=FALSE]
betap_pc_fc <- betap[match(pheno_pc_fc$Sample_ID, rownames(betap)), , drop=FALSE]

stopifnot(
  identical(rownames(betap_fc), pheno_fc$Sample_ID),
  identical(rownames(betap_pc_fc), pheno_pc_fc$Sample_ID),
  identical(levels(pheno_fc$LH_smoke_status), c("Current", "Former")),
  identical(levels(pheno_pc_fc$LH_smoke_status), c("Current", "Former"))
)

## Extract the direct Former minus Current coefficient.
extract_former_current <- function(mdl, meth, dat1, use_t=FALSE) {
  b <- coef(mdl)
  V <- vcov(mdl)
  coefficient_name <- "LH_smoke_statusFormer"

  if (!coefficient_name %in% names(b)) {
    stop("The direct Former minus Current coefficient was not estimable")
  }

  estimate <- as.numeric(b[coefficient_name])
  se <- as.numeric(sqrt(V[coefficient_name, coefficient_name]))
  statistic <- estimate/se

  if (use_t) {
    p_value <- 2*pt(abs(statistic), df=mdl$df.residual, lower.tail=FALSE)
  } else {
    p_value <- 2*pnorm(abs(statistic), lower.tail=FALSE)
  }

  N <- sum(!is.na(meth))
  n_former <- sum(!is.na(meth) & dat1$LH_smoke_status == "Former")
  n_current <- sum(!is.na(meth) & dat1$LH_smoke_status == "Current")

  out <- c(estimate, se, statistic, p_value, N, n_former, n_current)
  names(out) <- c("Beta", "SE", "tstat", "P", "N", "N.exp", "N.unexp")
  out
}

fit_unweighted_euro <- function(meth, dat1) {
  dat <- svydesign(id=~Sample_ID, data=data.frame(y=meth, dat1))
  mdl <- svyglm(
    y ~ LH_smoke_status + ANYETS + LHSSCREEN_12142015 + Gender +
      LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono +
      Neutro + Eosino + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
      PC8 + PC9 + PC10,
    design=dat
  )
  extract_former_current(mdl, meth, dat1, use_t=TRUE)
}

fit_unweighted_trans <- function(meth, dat1) {
  dat <- svydesign(id=~Sample_ID, data=data.frame(y=meth, dat1))
  mdl <- svyglm(
    y ~ LH_smoke_status + ANYETS + European + LHSSCREEN_12142015 +
      Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T +
      Mono + Neutro + Eosino + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +
      PC7 + PC8 + PC9 + PC10,
    design=dat
  )
  extract_former_current(mdl, meth, dat1, use_t=FALSE)
}

fit_weighted_euro <- function(meth, dat1) {
  dat <- svydesign(
    id=~Sample_ID,
    weights=~asthma_weights,
    data=data.frame(y=meth, dat1)
  )
  mdl <- svyglm(
    y ~ LH_smoke_status + ANYETS + Gender + LH_BMI + AgeAtHomeVisit +
      STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino + PC1 +
      PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
    design=dat
  )
  extract_former_current(mdl, meth, dat1, use_t=FALSE)
}

fit_weighted_trans <- function(meth, dat1) {
  dat <- svydesign(
    id=~Sample_ID,
    weights=~asthma_weights,
    data=data.frame(y=meth, dat1)
  )
  mdl <- svyglm(
    y ~ LH_smoke_status + ANYETS + European + Gender + LH_BMI +
      AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro +
      Eosino + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +
      PC10,
    design=dat
  )
  extract_former_current(mdl, meth, dat1, use_t=FALSE)
}

cl <- makeCluster(getOption("cl.cores", 32))
clusterEvalQ(cl, {
  library(survey)
  NULL
})
clusterExport(cl, "extract_former_current", envir=environment())

run_model <- function(beta_data, phenotype_data, model_function, output_name) {
  results <- parApply(
    cl,
    beta_data,
    2,
    model_function,
    dat1=phenotype_data
  )
  results <- as.data.frame(t(results))
  write.csv(
    results,
    file=file.path(output_dir, output_name),
    row.names=TRUE
  )
  results
}

## N.exp is former smokers and N.unexp is current smokers.
results_unweighted_euro <- run_model(
  betap_pc_fc,
  pheno_pc_fc,
  fit_unweighted_euro,
  "unimputed_unweighted_restricted_former_vs_current_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv"
)

results_unweighted_trans <- run_model(
  betap_fc,
  pheno_fc,
  fit_unweighted_trans,
  "unimputed_unweighted_restricted_former_vs_current_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv"
)

results_weighted_euro <- run_model(
  betap_pc_fc,
  pheno_pc_fc,
  fit_weighted_euro,
  "unimputed_weighted_restricted_former_vs_current_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv"
)

results_weighted_trans <- run_model(
  betap_fc,
  pheno_fc,
  fit_weighted_trans,
  "unimputed_weighted_restricted_former_vs_current_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv"
)

stopCluster(cl)

## Quick genome-wide counts using the same probe-design Bonferroni denominator.
bonf <- 0.05/230325

significant_counts <- data.frame(
  Model=c(
    "Unweighted European",
    "Unweighted transancestry",
    "Weighted European",
    "Weighted transancestry"
  ),
  Significant_probe_designs=c(
    sum(results_unweighted_euro$P <= bonf, na.rm=TRUE),
    sum(results_unweighted_trans$P <= bonf, na.rm=TRUE),
    sum(results_weighted_euro$P <= bonf, na.rm=TRUE),
    sum(results_weighted_trans$P <= bonf, na.rm=TRUE)
  )
)

write.table(
  significant_counts,
  file=file.path(
    output_dir,
    "restricted_former_vs_current_significant_counts.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

print(significant_counts)


###############################################################################
## Annotate all results and save Bonferroni-significant results separately
###############################################################################

annotated_dir <- file.path(output_dir, "annotated")
annotated_filtered_dir <- file.path(output_dir, "annotated_filtered")

dir.create(annotated_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(annotated_filtered_dir, recursive=TRUE, showWarnings=FALSE)

reference_dir <- file.path(project_dir, "references")

msa_annotation <- data.table::fread(
  file.path(reference_dir, "MSA.hg38.manifest.gencode.v41.tsv.gz")
)

epic_cpgs <- data.table::fread(
  file.path(reference_dir, "EPIC.hg38.manifest.gencode.v36.tsv.gz"),
  select="probeID"
)$probeID

epicv2_cpgs <- data.table::fread(
  file.path(reference_dir, "EPICv2.hg38.manifest.gencode.v41.tsv.gz"),
  select="probeID"
)$probeID

hm450_cpgs <- data.table::fread(
  file.path(reference_dir, "HM450.hg38.manifest.tsv.gz"),
  select="Probe_ID"
)$Probe_ID

epic_cpgs <- unique(sub("_.*$", "", epic_cpgs))
epicv2_cpgs <- unique(sub("_.*$", "", epicv2_cpgs))
hm450_cpgs <- unique(sub("_.*$", "", hm450_cpgs))

annotate_ewas_results <- function(results) {
  annotated <- data.frame(
    Probe_ID=rownames(results),
    results,
    row.names=NULL,
    check.names=FALSE
  ) %>%
    dplyr::mutate(
      CpG_ID=sub("_.*$", "", Probe_ID),
      on_EPIC=CpG_ID %in% epic_cpgs,
      on_EPICv2=CpG_ID %in% epicv2_cpgs,
      on_450K=CpG_ID %in% hm450_cpgs,
      array_overlap=dplyr::case_when(
        on_EPIC & on_EPICv2 & on_450K ~ "EPIC_EPICv2_450K",
        on_EPIC & on_EPICv2 ~ "EPIC_EPICv2",
        on_EPIC & on_450K ~ "EPIC_450K",
        on_EPICv2 & on_450K ~ "EPICv2_450K",
        on_EPIC ~ "EPIC",
        on_EPICv2 ~ "EPICv2",
        on_450K ~ "450K",
        TRUE ~ "MSA_ONLY"
      )
    ) %>%
    dplyr::left_join(msa_annotation, by=c("Probe_ID"="probeID")) %>%
    dplyr::relocate(
      dplyr::any_of(c(
        "Probe_ID", "CpG_ID", "CpG_chrm", "CpG_beg", "CpG_end",
        "probe_strand", "genesUniq", "on_EPIC", "on_EPICv2", "on_450K",
        "array_overlap"
      ))
    )
  
  annotated
}

model_results <- list(
  unimputed_unweighted_restricted_former_vs_current_euro=results_unweighted_euro,
  unimputed_unweighted_restricted_former_vs_current_transancestry=results_unweighted_trans,
  unimputed_weighted_restricted_former_vs_current_euro=results_weighted_euro,
  unimputed_weighted_restricted_former_vs_current_transancestry=results_weighted_trans
)

for (model_name in names(model_results)) {
  annotated_results <- annotate_ewas_results(model_results[[model_name]])
  significant_annotated_results <- annotated_results %>%
    dplyr::filter(P <= bonf)
  
  data.table::fwrite(
    annotated_results,
    file=file.path(
      annotated_dir,
      paste0(model_name, "_annotated.tsv")
    ),
    sep="\t",
    na="NA"
  )
  
  data.table::fwrite(
    significant_annotated_results,
    file=file.path(
      annotated_filtered_dir,
      paste0(model_name, "_significant_annotated.tsv")
    ),
    sep="\t",
    na="NA"
  )
}
