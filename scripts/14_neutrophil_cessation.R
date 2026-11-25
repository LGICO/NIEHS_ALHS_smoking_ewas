library(survey)
library(dplyr)
library(data.table)
library(ggplot2)
library(scales)

project_dir <- ""
output_dir <- file.path(project_dir, "output", "neutrophil_cessation")
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

load('751_clean_pc.rdata')

new_pheno <- data.table::fread(file.path(project_dir,'samplesheets/lhs28901_Methylation_Phenotype_2025APR29.txt'),sep=',')
np <- new_pheno %>% filter(LHID %in% pheno$LHID, COMMENT != 'Replicate')
x <- colnames(np)[!colnames(np) %in% colnames(pheno)]
np <- np %>% dplyr::select(LHID, partid, all_of(x))
pheno <- left_join(pheno, np, by=c('LHID','partid'))


clean_smoking <- function(x) {
  x <- as.character(x)
  out <- ifelse(x %in% c("Never", "0"), "Never", ifelse(x %in% c("Former", "1"), "Former", ifelse(x %in% c("Current", "2"), "Current", NA)))
  factor(out, levels=c("Never", "Former", "Current"))
}

clean_ancestry <- function(x) {
  x <- as.character(x)
  out <- ifelse(x %in% c("European", "1"), "European", ifelse(x %in% c("Non_European", "Non-European", "0"), "Non-European", NA))
  factor(out, levels=c("European", "Non-European"))
}

clean_case <- function(x) {
  x <- as.character(x)
  out <- ifelse(x %in% c("Case", "1"), "Case", ifelse(x %in% c("Control", "0"), "Control", NA))
  factor(out, levels=c("Control", "Case"))
}

clean_gender <- function(x) {
  x <- tolower(as.character(x))
  out <- ifelse(x %in% c("male", "m", "0"), "Male", ifelse(x %in% c("female", "f", "1"), "Female", NA))
  factor(out, levels=c("Male", "Female"))
}

clean_state <- function(x) {
  x <- as.character(x)
  out <- ifelse(x %in% c("IA", "0"), "IA", ifelse(x %in% c("NC", "1"), "NC", NA))
  factor(out, levels=c("IA", "NC"))
}

clean_anyets <- function(x) {
  x <- as.character(x)
  out <- ifelse(x %in% c("No", "0"), "No", ifelse(x %in% c("Yes", "1"), "Yes", "Unknown"))
  factor(out, levels=c("No", "Yes", "Unknown"))
}

pheno_neutro <- pheno %>%
  mutate(
    smoke_status=clean_smoking(LH_smoke_status),
    ancestry=clean_ancestry(European),
    screening_status=clean_case(LHSSCREEN_12142015),
    gender_model=clean_gender(Gender),
    state_model=clean_state(STATE),
    anyets_model=clean_anyets(ANYETS),
    years_since_quit=as.numeric(LH_TimeSinceSmked_v1)
  )

################################################################################
## Add asthma sampling weights before restricting to former smokers

add_weights <- function(dat) {
  p_asthma_overall <- 0.06852481
  p_asthma_nested <- mean(dat$screening_status == "Case", na.rm=TRUE)
  w_asthma <- 1/(p_asthma_nested/p_asthma_overall)
  w_nonasthma <- 1/((1-p_asthma_nested)/(1-p_asthma_overall))
  dat %>% mutate(asthma_weights=ifelse(screening_status == "Case", w_asthma, w_nonasthma))
}

trans_data <- add_weights(pheno_neutro) %>%
  filter(smoke_status == "Former", !is.na(years_since_quit), !is.na(Neutro))

euro_data <- add_weights(pheno_neutro %>% filter(ancestry == "European")) %>%
  filter(smoke_status == "Former", !is.na(years_since_quit), !is.na(Neutro))

################################################################################
## Cessation-interval summaries

trans_data <- trans_data %>%
  mutate(cessation_interval=case_when(
    years_since_quit <= 10 ~ "≤10 years",
    years_since_quit <= 20 ~ "11–20 years",
    years_since_quit <= 30 ~ "21–30 years",
    years_since_quit > 30 ~ ">30 years"
  ))

trans_data$cessation_interval <- factor(trans_data$cessation_interval, levels=c("≤10 years", "11–20 years", "21–30 years", ">30 years"))

neutrophil_interval_summary <- trans_data %>%
  group_by(cessation_interval) %>%
  summarise(
    N=dplyr::n(),
    Median_neutrophil_fraction=median(Neutro, na.rm=TRUE),
    Q1=quantile(Neutro, 0.25, na.rm=TRUE),
    Q3=quantile(Neutro, 0.75, na.rm=TRUE),
    Mean_neutrophil_fraction=mean(Neutro, na.rm=TRUE),
    SD=sd(Neutro, na.rm=TRUE),
    .groups="drop"
  )

interval_test <- kruskal.test(Neutro ~ cessation_interval, data=trans_data)

################################################################################
## Unadjusted Spearman correlations

get_spearman <- function(dat, population) {
  test <- cor.test(dat$years_since_quit, dat$Neutro, method="spearman", exact=FALSE)
  
  data.frame(
    Population=population,
    N=nrow(dat),
    Spearman_rho=unname(test$estimate),
    Spearman_P=test$p.value
  )
}

spearman_results <- bind_rows(
  get_spearman(trans_data, "Transancestry"),
  get_spearman(euro_data, "European")
)

################################################################################
## Adjusted survey-weighted models

fit_neutrophil_model <- function(dat, population) {
  covariates <- c("years_since_quit", "AgeAtHomeVisit", "gender_model", "LH_BMI", "state_model", "anyets_model", paste0("PC", 1:10))
  if (population == "Transancestry") covariates <- c(covariates, "ancestry")
  
  required_variables <- c("Sample_ID", "Neutro", "asthma_weights", covariates)
  dat <- dat[complete.cases(dat[, required_variables]), , drop=FALSE]
  
  model_formula <- as.formula(paste("Neutro ~", paste(covariates, collapse=" + ")))
  design <- svydesign(id=~Sample_ID, weights=~asthma_weights, data=dat)
  model <- svyglm(model_formula, design=design)
  
  coefficient_table <- summary(model)$coefficients
  confidence_interval <- confint(model)["years_since_quit", ]
  
  data.frame(
    Population=population,
    N=nrow(dat),
    Beta_per_year=coefficient_table["years_since_quit", "Estimate"],
    SE=coefficient_table["years_since_quit", "Std. Error"],
    Lower_95_CI=confidence_interval[1],
    Upper_95_CI=confidence_interval[2],
    P=coefficient_table["years_since_quit", "Pr(>|t|)"],
    Beta_per_10_years=10*coefficient_table["years_since_quit", "Estimate"],
    Lower_95_CI_per_10_years=10*confidence_interval[1],
    Upper_95_CI_per_10_years=10*confidence_interval[2]
  )
}

adjusted_results <- bind_rows(
  fit_neutrophil_model(trans_data, "Transancestry"),
  fit_neutrophil_model(euro_data, "European")
)

################################################################################
## Descriptive scatterplot

ancestry_colors <- c("European"="#4C78A8", "Non-European"="#E45756")

neutrophil_plot <- ggplot(trans_data, aes(x=years_since_quit, y=Neutro, color=ancestry)) +
  geom_point(alpha=0.55, size=1.7) +
  geom_smooth(aes(group=1), method="lm", formula=y~x, se=TRUE, color="black", linewidth=0.9) +
  scale_color_manual(values=ancestry_colors) +
  scale_y_continuous(labels=label_percent(accuracy=1)) +
  labs(
    title="Neutrophil proportion and time since smoking cessation",
    subtitle=paste0("Former smokers, n=", nrow(trans_data), "; fitted line is unadjusted"),
    x="Years since smoking cessation",
    y="Estimated neutrophil proportion",
    color="Ancestry group"
  ) +
  theme_classic(base_size=13) +
  theme(
    plot.title=element_text(face="bold"),
    legend.position="bottom"
  )

ggsave(file.path(output_dir, "neutrophil_fraction_years_since_cessation.pdf"), neutrophil_plot, width=8, height=6)
ggsave(file.path(output_dir, "neutrophil_fraction_years_since_cessation.png"), neutrophil_plot, width=8, height=6, dpi=300)

################################################################################
## Save and print results

data.table::fwrite(neutrophil_interval_summary, file.path(output_dir, "neutrophil_by_cessation_interval.tsv"), sep="\t")
data.table::fwrite(spearman_results, file.path(output_dir, "neutrophil_cessation_spearman.tsv"), sep="\t")
data.table::fwrite(adjusted_results, file.path(output_dir, "neutrophil_cessation_adjusted_models.tsv"), sep="\t")

interval_test_result <- data.frame(
  Test="Kruskal-Wallis across cessation intervals",
  Statistic=unname(interval_test$statistic),
  DF=unname(interval_test$parameter),
  P=interval_test$p.value
)

data.table::fwrite(interval_test_result, file.path(output_dir, "neutrophil_cessation_interval_test.tsv"), sep="\t")

print(neutrophil_interval_summary)
print(interval_test_result)
print(spearman_results)
print(adjusted_results)
