library(dplyr)
library(data.table)
library(survey)
library(lmtest)
library(parallel)

project_dir <- ""
load('751_clean_pc.rdata')

betap <- t(betas_miss_masked)
storage.mode(betap) <- "double"

requested_method <- commandArgs(trailingOnly=TRUE)
if (length(requested_method) == 0) requested_method <- "both"
methods <- if (requested_method[1] == "both") c("mean", "median") else requested_method[1]
if (!all(methods %in% c("mean", "median"))) stop("Argument must be mean, median, or both")

################################################################################
## Preserve existing phenotype coding, including unknown ANYETS

code_smoking <- function(x) {
  x <- as.character(x)
  out <- ifelse(x %in% c("Never", "0"), 0,
                ifelse(x %in% c("Former", "1"), 1,
                       ifelse(x %in% c("Current", "2"), 2, 99)))
  factor(out, levels=c(0, 1, 2, 99))
}

code_ancestry <- function(x) {
  x <- as.character(x)
  out <- ifelse(x %in% c("European", "1"), 1,
                ifelse(x %in% c("Non_European", "Non-European", "0"), 0, 99))
  factor(out, levels=c(0, 1, 99))
}

code_binary <- function(x, negative, positive) {
  x <- as.character(x)
  out <- ifelse(x %in% c(negative, "0"), 0,
                ifelse(x %in% c(positive, "1"), 1, 99))
  factor(out, levels=c(0, 1, 99))
}

pheno <- pheno %>%
  mutate(
    LH_smoke_status=code_smoking(LH_smoke_status),
    European=code_ancestry(European),
    ANYETS=code_binary(ANYETS, "No", "Yes"),
    LHSSCREEN_12142015=code_binary(LHSSCREEN_12142015, "Control", "Case"),
    Gender=code_binary(tolower(as.character(Gender)), "male", "female"),
    STATE=code_binary(STATE, "IA", "NC")
  )

################################################################################
## Add asthma sampling weights exactly as in the weighted models

add_asthma_weights <- function(dat) {
  p_asthma_overall <- 0.06852481
  n_asthma <- sum(dat$LHSSCREEN_12142015 == 1, na.rm=TRUE)
  n_nonasthma <- sum(dat$LHSSCREEN_12142015 == 0, na.rm=TRUE)
  p_asthma_nested <- n_asthma/(n_asthma+n_nonasthma)
  w_asthma <- 1/(p_asthma_nested/p_asthma_overall)
  w_nonasthma <- 1/((1-p_asthma_nested)/(1-p_asthma_overall))
  
  dat$asthma_weights <- ifelse(
    dat$LHSSCREEN_12142015 == 1,
    w_asthma,
    ifelse(dat$LHSSCREEN_12142015 == 0, w_nonasthma, NA)
  )
  
  dat
}

pheno_trans <- add_asthma_weights(pheno %>% filter(!is.na(PC1)))
pheno_euro <- add_asthma_weights(pheno %>% filter(European == 1, !is.na(PC1)))

################################################################################
## Collapse MSA probe designs to base CpGs

probe_ids <- colnames(betap)
base_cpg_ids <- sub("_.*$", "", probe_ids)
unique_cpg_ids <- unique(base_cpg_ids)
first_design_index <- match(unique_cpg_ids, base_cpg_ids)
design_counts <- table(base_cpg_ids)
duplicated_cpg_ids <- names(design_counts[design_counts > 1])

collapse_cpgs <- function(x, method) {
  message("Creating ", method, "-collapsed CpG matrix")
  
  collapsed <- x[, first_design_index, drop=FALSE]
  colnames(collapsed) <- unique_cpg_ids
  
  for (i in seq_along(duplicated_cpg_ids)) {
    cpg <- duplicated_cpg_ids[i]
    design_index <- which(base_cpg_ids == cpg)
    
    if (method == "mean") {
      value <- rowMeans(x[, design_index, drop=FALSE], na.rm=TRUE)
    } else {
      value <- apply(x[, design_index, drop=FALSE], 1, median, na.rm=TRUE)
    }
    
    value[!is.finite(value)] <- NA_real_
    collapsed[, match(cpg, unique_cpg_ids)] <- value
    
    if (i %% 1000 == 0) message(method, ": collapsed ", i, " of ", length(duplicated_cpg_ids), " duplicated CpGs")
  }
  
  collapsed
}

################################################################################
## Existing weighted-model formulas

formula_euro <- paste(
  "y ~ LH_smoke_status + ANYETS + Gender + LH_BMI + AgeAtHomeVisit + STATE +",
  "B + NK + CD4T + CD8T + Mono + Neutro + Eosino +",
  "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
)

formula_trans <- paste(
  "y ~ LH_smoke_status + ANYETS + European + Gender + LH_BMI + AgeAtHomeVisit + STATE +",
  "B + NK + CD4T + CD8T + Mono + Neutro + Eosino +",
  "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
)

fit_weighted_cpg <- function(meth, dat1, model_formula, coefficient_name, exposed_code) {
  tryCatch({
    design <- survey::svydesign(id=~Sample_ID, weights=~asthma_weights, data=data.frame(y=meth, dat1))
    model <- survey::svyglm(as.formula(model_formula), design=design)
    coefficient_table <- lmtest::coeftest(model, vcov=vcov(model))
    p_column <- grep("^Pr", colnames(coefficient_table), value=TRUE)[1]
    
    out <- c(
      Beta=coefficient_table[coefficient_name, "Estimate"],
      SE=coefficient_table[coefficient_name, "Std. Error"],
      P=coefficient_table[coefficient_name, p_column],
      N=sum(!is.na(meth)),
      N.exp=sum(!is.na(meth) & as.character(dat1$LH_smoke_status) == exposed_code),
      N.unexp=sum(!is.na(meth) & as.character(dat1$LH_smoke_status) == "0")
    )
    
    out
  }, error=function(e) {
    c(Beta=NA_real_, SE=NA_real_, P=NA_real_, N=sum(!is.na(meth)), N.exp=NA_real_, N.unexp=NA_real_)
  })
}

################################################################################
## Manifest collapsed to one annotation row per base CpG

manifest <- data.table::fread(file.path(project_dir, "references", "MSA.hg38.manifest.gencode.v41.tsv.gz"))
manifest <- manifest[probeID %in% probe_ids]
manifest[, CpG_ID := sub("_.*$", "", probeID)]

manifest_base <- manifest[, .(
  CpG_chrm=first(CpG_chrm),
  CpG_beg=first(CpG_beg),
  CpG_end=first(CpG_end),
  probe_strand=first(probe_strand),
  Number_MSA_designs=.N,
  MSA_probe_designs=paste(probeID, collapse=";"),
  genesUniq=first(genesUniq),
  geneNames=first(geneNames),
  transcriptTypes=first(transcriptTypes),
  transcriptIDs=first(transcriptIDs),
  distToTSS=first(distToTSS)
), by=CpG_ID]

################################################################################
## Parallel cluster

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset=getOption("cl.cores", 60)))
n_cores <- max(1, min(n_cores, parallel::detectCores()))

cl <- parallel::makeCluster(n_cores)
clusterEvalQ(cl, {
  library(survey)
  library(lmtest)
})

################################################################################
## Run one model and save full and significant annotated results

run_model <- function(meth_matrix, phenotype, model_formula, coefficient_name, exposed_code, model_name, output_dir) {
  phenotype <- phenotype[match(phenotype$Sample_ID, rownames(meth_matrix)), , drop=FALSE]
  stopifnot(identical(phenotype$Sample_ID, rownames(meth_matrix)))
  
  message("Running ", model_name, " with ", ncol(meth_matrix), " base CpGs")
  
  result <- parallel::parApply(
    cl,
    meth_matrix,
    2,
    fit_weighted_cpg,
    dat1=phenotype,
    model_formula=model_formula,
    coefficient_name=coefficient_name,
    exposed_code=exposed_code
  )
  
  result <- as.data.frame(t(result))
  result <- data.table::as.data.table(result, keep.rownames="CpG_ID")
  
  number_tested <- sum(!is.na(result$P))
  bonferroni_threshold <- 0.05/number_tested
  
  result[, Bonferroni_threshold := bonferroni_threshold]
  result[, Significant_Bonferroni := !is.na(P) & P <= bonferroni_threshold]
  
  result_annotated <- merge(result, manifest_base, by="CpG_ID", all.x=TRUE, sort=FALSE)
  significant_annotated <- result_annotated[Significant_Bonferroni == TRUE]
  
  data.table::fwrite(result_annotated, file.path(output_dir, paste0(model_name, ".tsv")), sep="\t")
  data.table::fwrite(significant_annotated, file.path(output_dir, paste0(model_name, "_significant.tsv")), sep="\t")
  
  data.frame(
    Model=model_name,
    Number_base_CpGs=ncol(meth_matrix),
    Number_tested=number_tested,
    Bonferroni_threshold=bonferroni_threshold,
    Significant_CpGs=nrow(significant_annotated)
  )
}

################################################################################
## Run mean and median analyses

all_model_summaries <- list()
summary_index <- 1

for (method in methods) {
  output_dir <- file.path(project_dir, "output", paste0("cpg_collapsed_", method))
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  
  collapsed <- collapse_cpgs(betap, method)
  rownames(collapsed) <- rownames(betap)
  
  collapse_summary <- data.frame(
    Collapse_method=method,
    Original_probe_designs=ncol(betap),
    Unique_base_CpGs=ncol(collapsed),
    CpGs_with_multiple_designs=length(duplicated_cpg_ids)
  )
  
  data.table::fwrite(collapse_summary, file.path(output_dir, "collapse_summary.tsv"), sep="\t")
  
  collapsed_euro <- collapsed[match(pheno_euro$Sample_ID, rownames(collapsed)), , drop=FALSE]
  collapsed_trans <- collapsed[match(pheno_trans$Sample_ID, rownames(collapsed)), , drop=FALSE]
  
  stopifnot(identical(rownames(collapsed_euro), pheno_euro$Sample_ID))
  stopifnot(identical(rownames(collapsed_trans), pheno_trans$Sample_ID))
  
  all_model_summaries[[summary_index]] <- run_model(
    collapsed_euro,
    pheno_euro,
    formula_euro,
    "LH_smoke_status2",
    "2",
    "weighted_euro_current_vs_never",
    output_dir
  )
  all_model_summaries[[summary_index]]$Collapse_method <- method
  summary_index <- summary_index+1
  
  all_model_summaries[[summary_index]] <- run_model(
    collapsed_euro,
    pheno_euro,
    formula_euro,
    "LH_smoke_status1",
    "1",
    "weighted_euro_former_vs_never",
    output_dir
  )
  all_model_summaries[[summary_index]]$Collapse_method <- method
  summary_index <- summary_index+1
  
  all_model_summaries[[summary_index]] <- run_model(
    collapsed_trans,
    pheno_trans,
    formula_trans,
    "LH_smoke_status2",
    "2",
    "weighted_trans_current_vs_never",
    output_dir
  )
  all_model_summaries[[summary_index]]$Collapse_method <- method
  summary_index <- summary_index+1
  
  all_model_summaries[[summary_index]] <- run_model(
    collapsed_trans,
    pheno_trans,
    formula_trans,
    "LH_smoke_status1",
    "1",
    "weighted_trans_former_vs_never",
    output_dir
  )
  all_model_summaries[[summary_index]]$Collapse_method <- method
  summary_index <- summary_index+1
  
  rm(collapsed, collapsed_euro, collapsed_trans)
  gc()
}

stopCluster(cl)

model_summary <- bind_rows(all_model_summaries) %>%
  select(Collapse_method, everything())

data.table::fwrite(
  model_summary,
  file.path(project_dir, "output", "cpg_collapsed_mean_median_summary.tsv"),
  sep="\t"
)

print(model_summary)
