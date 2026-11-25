library(survey)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(ggplot2)

project_dir <- ""
output_dir <- file.path(project_dir, "output", "age_ecel1p1")
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

load('751_clean_pc.rdata')

## Use the unimputed/missingness-masked methylation values
if (!exists("betap")) betap <- t(betas_miss_masked)
betap <- as.matrix(betap)

stopifnot(exists("pheno"))
stopifnot(!is.null(rownames(betap)))
stopifnot(!anyDuplicated(pheno$Sample_ID))
stopifnot(!anyDuplicated(rownames(betap)))
stopifnot(all(pheno$Sample_ID %in% rownames(betap)))

################################################################################
## Confirm that the reported >10-percentage-point result is from the adjusted EWAS

target_probe <- "cg04284215_BC11"

primary_files <- c(
  European=file.path(project_dir, "output", "svyglm_anno", "unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv"),
  Transancestry=file.path(project_dir, "output", "svyglm_anno", "unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv")
)

extract_primary_result <- function(input_file, population) {
  x <- data.table::fread(input_file)
  id_column <- intersect(c("Probe_ID", "V1", "probeID"), names(x))[1]
  if (is.na(id_column)) stop("No probe-ID column found in ", input_file)
  x <- x[x[[id_column]] == target_probe, ]
  if (nrow(x) != 1) stop("Expected one row for ", target_probe, " in ", input_file)
  
  data.frame(
    Population=population,
    Probe_ID=target_probe,
    Beta=x$Beta,
    SE=x$SE,
    P=x$P,
    Difference_percentage_points=100*x$Beta,
    Absolute_difference_percentage_points=100*abs(x$Beta),
    Greater_than_10_percentage_points=abs(x$Beta) > 0.10
  )
}

primary_result <- bind_rows(
  extract_primary_result(primary_files["European"], "European"),
  extract_primary_result(primary_files["Transancestry"], "Transancestry")
)

data.table::fwrite(primary_result, file.path(output_dir, "cg04284215_primary_EWAS_coefficient.tsv"), sep="\t")

################################################################################
## Standardize phenotype variables without changing the underlying data

clean_smoking <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  out[x %in% c("Never", "0")] <- "Never"
  out[x %in% c("Former", "1")] <- "Former"
  out[x %in% c("Current", "2")] <- "Current"
  factor(out, levels=c("Never", "Former", "Current"))
}

clean_anyets <- function(x) {
  x <- trimws(as.character(x))
  out <- rep("Unknown", length(x))
  out[x %in% c("No", "0")] <- "No"
  out[x %in% c("Yes", "1")] <- "Yes"
  factor(out, levels=c("No", "Yes", "Unknown"))
}

clean_gender <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[x %in% c("male", "m", "0")] <- "Male"
  out[x %in% c("female", "f", "1")] <- "Female"
  factor(out, levels=c("Male", "Female"))
}

clean_state <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  out[x %in% c("IA", "0")] <- "IA"
  out[x %in% c("NC", "1")] <- "NC"
  factor(out, levels=c("IA", "NC"))
}

clean_ancestry <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  out[x %in% c("European", "1")] <- "European"
  out[x %in% c("Non_European", "Non-European", "0")] <- "Non-European"
  factor(out, levels=c("European", "Non-European"))
}

clean_screening <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  out[x %in% c("Control", "0")] <- "Control"
  out[x %in% c("Case", "1")] <- "Case"
  factor(out, levels=c("Control", "Case"))
}

pheno_age <- pheno %>%
  mutate(
    smoke_status=clean_smoking(LH_smoke_status),
    ANYETS_model=clean_anyets(ANYETS),
    Gender_model=clean_gender(Gender),
    STATE_model=clean_state(STATE),
    European_model=clean_ancestry(European),
    screening_model=clean_screening(LHSSCREEN_12142015)
  )

################################################################################
## Create cg04284215 and the ECEL1P1 regional mean

ecel1p1_displayed_cpgs <- c(
  "cg03329539",
  "cg04284170",
  "cg04284171",
  "cg06644428",
  "cg04284185",
  "cg04284186",
  "cg04284205",
  "cg04284207",
  "cg04284209",
  "cg04284215",
  "cg21566642",
  "cg13193840"
)

## cg13193840 was added as the right-side boundary marker, not as part of the DMR
ecel1p1_core_cpgs <- setdiff(ecel1p1_displayed_cpgs, "cg13193840")

probe_base_cpg <- sub("_.*$", "", colnames(betap))

probe_mapping <- data.frame(
  Probe_ID=colnames(betap)[probe_base_cpg %in% ecel1p1_displayed_cpgs],
  CpG_ID=probe_base_cpg[probe_base_cpg %in% ecel1p1_displayed_cpgs]
) %>%
  mutate(DMR_role=ifelse(CpG_ID == "cg13193840", "Flanking boundary marker", "Core DMR"))

missing_core_cpgs <- setdiff(ecel1p1_core_cpgs, unique(probe_mapping$CpG_ID))
if (length(missing_core_cpgs) > 0) stop("Missing ECEL1P1 CpGs: ", paste(missing_core_cpgs, collapse=", "))

data.table::fwrite(probe_mapping, file.path(output_dir, "ecel1p1_probe_mapping.tsv"), sep="\t")

## Average multiple MSA probe designs within each base CpG
collapsed_cpgs <- do.call(cbind, lapply(ecel1p1_core_cpgs, function(cpg) {
  designs <- which(probe_base_cpg == cpg)
  value <- rowMeans(betap[, designs, drop=FALSE], na.rm=TRUE)
  value[!is.finite(value)] <- NA_real_
  value
}))

colnames(collapsed_cpgs) <- ecel1p1_core_cpgs
rownames(collapsed_cpgs) <- rownames(betap)

## Require at least 80% of the core CpGs for the regional mean
minimum_cpgs <- ceiling(0.80*ncol(collapsed_cpgs))
dmr_n_cpgs <- rowSums(!is.na(collapsed_cpgs))
dmr_mean <- rowMeans(collapsed_cpgs, na.rm=TRUE)
dmr_mean[dmr_n_cpgs < minimum_cpgs | !is.finite(dmr_mean)] <- NA_real_

methylation_index <- match(pheno_age$Sample_ID, rownames(betap))
pheno_age$cg04284215 <- betap[methylation_index, target_probe]
pheno_age$ECEL1P1_DMR_mean <- dmr_mean[methylation_index]
pheno_age$ECEL1P1_DMR_n_CpGs <- dmr_n_cpgs[methylation_index]

################################################################################
## Survey weights and population datasets

p_asthma_overall <- 0.06852481

prepare_population <- function(dat, population) {
  if (population == "European") dat <- dat %>% filter(European_model == "European")
  
  p_asthma_nested <- mean(dat$screening_model == "Case", na.rm=TRUE)
  w_asthma <- 1/(p_asthma_nested/p_asthma_overall)
  w_nonasthma <- 1/((1-p_asthma_nested)/(1-p_asthma_overall))
  
  dat <- dat %>%
    mutate(asthma_weights=ifelse(screening_model == "Case", w_asthma, w_nonasthma)) %>%
    filter(!is.na(smoke_status)) %>%
    mutate(Age_c=AgeAtHomeVisit-mean(AgeAtHomeVisit, na.rm=TRUE))
  
  dat
}

################################################################################
## Fit adjusted smoking, age, and smoking-by-age models

cell_covariates <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino")
pc_covariates <- paste0("PC", 1:10)

contrast_result <- function(model, weights, label, population, outcome) {
  coefficient_names <- names(coef(model))
  contrast_vector <- setNames(rep(0, length(coefficient_names)), coefficient_names)
  contrast_vector[names(weights)] <- weights
  
  test <- survey::svycontrast(model, contrast_vector)
  estimate <- as.numeric(coef(test))
  standard_error <- as.numeric(survey::SE(test))
  degrees_freedom <- df.residual(model)
  statistic <- estimate/standard_error
  p_value <- 2*pt(abs(statistic), df=degrees_freedom, lower.tail=FALSE)
  critical_value <- qt(0.975, df=degrees_freedom)
  
  data.frame(
    Population=population,
    Outcome=outcome,
    Contrast=label,
    Estimate=estimate,
    SE=standard_error,
    Lower_95_CI=estimate-critical_value*standard_error,
    Upper_95_CI=estimate+critical_value*standard_error,
    P=p_value
  )
}

fit_age_interaction_model <- function(dat, outcome, population) {
  covariates <- c("ANYETS_model", "Gender_model", "LH_BMI", "STATE_model", cell_covariates, pc_covariates)
  if (population == "Transancestry") covariates <- c(covariates, "European_model")
  
  dat$.outcome <- dat[[outcome]]
  required_variables <- unique(c("Sample_ID", "asthma_weights", ".outcome", "smoke_status", "Age_c", covariates))
  dat <- dat[complete.cases(dat[, required_variables]), , drop=FALSE]
  
  model_formula <- as.formula(paste(".outcome ~ smoke_status*Age_c +", paste(covariates, collapse=" + ")))
  design <- survey::svydesign(id=~Sample_ID, weights=~asthma_weights, data=dat)
  model <- survey::svyglm(model_formula, design=design)
  
  coefficient_names <- names(coef(model))
  former_interaction <- intersect(c("smoke_statusFormer:Age_c", "Age_c:smoke_statusFormer"), coefficient_names)
  current_interaction <- intersect(c("smoke_statusCurrent:Age_c", "Age_c:smoke_statusCurrent"), coefficient_names)
  
  if (length(former_interaction) != 1 || length(current_interaction) != 1) stop("Could not identify smoking-by-age coefficients")
  
  contrasts <- bind_rows(
    contrast_result(model, c(smoke_statusCurrent=1), "Current vs Never at mean age", population, outcome),
    contrast_result(model, c(smoke_statusFormer=1), "Former vs Never at mean age", population, outcome),
    contrast_result(model, c(smoke_statusCurrent=1, smoke_statusFormer=-1), "Current vs Former at mean age", population, outcome),
    contrast_result(model, c(Age_c=1), "Age slope among Never smokers", population, outcome),
    contrast_result(model, setNames(c(1, 1), c("Age_c", former_interaction)), "Age slope among Former smokers", population, outcome),
    contrast_result(model, setNames(c(1, 1), c("Age_c", current_interaction)), "Age slope among Current smokers", population, outcome)
  )
  
  smoking_test <- survey::regTermTest(model, ~smoke_status)
  interaction_test <- survey::regTermTest(model, ~smoke_status:Age_c)
  
  relative_estimates <- c(
    Never=0,
    Former=unname(coef(model)["smoke_statusFormer"]),
    Current=unname(coef(model)["smoke_statusCurrent"])
  )
  
  current_result <- contrasts %>% filter(Contrast == "Current vs Never at mean age")
  never_age_result <- contrasts %>% filter(Contrast == "Age slope among Never smokers")
  
  summary_result <- data.frame(
    Population=population,
    Outcome=outcome,
    N=nrow(dat),
    Mean_age=mean(dat$AgeAtHomeVisit),
    Current_vs_never_estimate=current_result$Estimate,
    Current_vs_never_SE=current_result$SE,
    Current_vs_never_P=current_result$P,
    Current_vs_never_percentage_points=100*current_result$Estimate,
    Age_coefficient_among_never=never_age_result$Estimate,
    Age_coefficient_SE=never_age_result$SE,
    Age_coefficient_P=never_age_result$P,
    Overall_smoking_P=as.numeric(smoking_test$p),
    Smoking_by_age_interaction_P=as.numeric(interaction_test$p),
    Adjusted_order_at_mean_age=paste(names(sort(relative_estimates)), collapse=" < ")
  )
  
  list(model=model, contrasts=contrasts, summary=summary_result, analysis_data=dat)
}

################################################################################
## Run transancestry and European analyses

all_summaries <- list()
all_contrasts <- list()
saved_models <- list()
analysis_counter <- 1

for (population in c("Transancestry", "European")) {
  population_data <- prepare_population(pheno_age, population)
  
  for (outcome in c("cg04284215", "ECEL1P1_DMR_mean")) {
    result <- fit_age_interaction_model(population_data, outcome, population)
    
    all_summaries[[analysis_counter]] <- result$summary
    all_contrasts[[analysis_counter]] <- result$contrasts
    saved_models[[paste(population, outcome, sep="_")]] <- result$model
    analysis_counter <- analysis_counter+1
  }
}

age_model_summary <- bind_rows(all_summaries)
age_model_contrasts <- bind_rows(all_contrasts)

data.table::fwrite(age_model_summary, file.path(output_dir, "age_smoking_interaction_summary.tsv"), sep="\t")
data.table::fwrite(age_model_contrasts, file.path(output_dir, "age_smoking_interaction_contrasts.tsv"), sep="\t")
saveRDS(saved_models, file.path(output_dir, "age_smoking_interaction_models.rds"))

################################################################################
## Unadjusted descriptive age scatterplots

plot_data <- prepare_population(pheno_age, "Transancestry") %>%
  dplyr::select(Sample_ID, smoke_status, AgeAtHomeVisit, cg04284215, ECEL1P1_DMR_mean) %>%
  pivot_longer(cols=c(cg04284215, ECEL1P1_DMR_mean), names_to="Outcome", values_to="Beta") %>%
  mutate(
    Outcome=recode(
      Outcome,
      cg04284215="cg04284215",
      ECEL1P1_DMR_mean="ECEL1P1 regional mean"
    )
  )

smoking_counts <- plot_data %>%
  dplyr::distinct(Sample_ID, smoke_status) %>%
  dplyr::count(smoke_status)

legend_labels <- setNames(
  paste0(smoking_counts$smoke_status, " (n=", smoking_counts$n, ")"),
  smoking_counts$smoke_status
)

smoking_colors <- c(
  Never="#D99A9A",
  Former="#AA3732",
  Current="#32100E"
)

age_plot <- ggplot(plot_data, aes(x=AgeAtHomeVisit, y=Beta, color=smoke_status, shape=smoke_status)) +
  geom_point(alpha=0.45, size=1.5) +
  geom_smooth(method="lm", formula=y~x, se=TRUE, linewidth=0.9) +
  facet_wrap(~Outcome, ncol=2, scales="free_y") +
  scale_color_manual(values=smoking_colors, breaks=names(smoking_colors), labels=legend_labels) +
  scale_shape_manual(values=c(Never=16, Former=17, Current=15), breaks=names(smoking_colors), labels=legend_labels) +
  labs(
    title="Age and ECEL1P1 methylation by smoking status",
    subtitle="Points and fitted lines are unadjusted; statistical results are from covariate-adjusted survey-weighted models",
    x="Age at home visit, years",
    y="Unimputed methylation beta value",
    color="Smoking status",
    shape="Smoking status"
  ) +
  theme_classic(base_size=13) +
  theme(
    strip.text=element_text(face="bold"),
    legend.position="bottom",
    plot.title=element_text(face="bold")
  )

ggsave(file.path(output_dir, "cg04284215_ECEL1P1_age_scatter.pdf"), age_plot, width=11, height=5.5)
ggsave(file.path(output_dir, "cg04284215_ECEL1P1_age_scatter.png"), age_plot, width=11, height=5.5, dpi=300)

data.table::fwrite(plot_data, file.path(output_dir, "cg04284215_ECEL1P1_age_scatter_data.tsv"), sep="\t")

################################################################################
## Raw descriptive summaries

raw_summary <- prepare_population(pheno_age, "Transancestry") %>%
  group_by(smoke_status) %>%
  summarise(
    N=dplyr::n(),
    Age_median=median(AgeAtHomeVisit, na.rm=TRUE),
    Age_Q1=quantile(AgeAtHomeVisit, 0.25, na.rm=TRUE),
    Age_Q3=quantile(AgeAtHomeVisit, 0.75, na.rm=TRUE),
    cg04284215_mean=mean(cg04284215, na.rm=TRUE),
    cg04284215_median=median(cg04284215, na.rm=TRUE),
    ECEL1P1_DMR_mean=mean(ECEL1P1_DMR_mean, na.rm=TRUE),
    .groups="drop"
  )

data.table::fwrite(raw_summary, file.path(output_dir, "age_methylation_descriptive_summary.tsv"), sep="\t")

print(primary_result)
print(age_model_summary)
print(age_model_contrasts)
print(raw_summary)
