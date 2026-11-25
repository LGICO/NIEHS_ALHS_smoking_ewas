library(data.table)

project_dir <- ""
sig_dir <- file.path(project_dir, "output/significant_models")
dir.create(sig_dir, recursive=TRUE, showWarnings=FALSE)

## Array manifests
epic <- fread(file.path(project_dir, "references/EPIC.hg38.manifest.gencode.v36.tsv.gz"), select="probeID")[[1]]
epic2 <- fread(file.path(project_dir, "references/EPICv2.hg38.manifest.gencode.v41.tsv.gz"), select="probeID")[[1]]
hm450 <- fread(file.path(project_dir, "references/HM450.hg38.manifest.tsv.gz"), select="Probe_ID")[[1]]

## CpGs present on EPIC, EPICv2, or 450K
known_cpgs <- unique(sub("_.*$", "", c(epic, epic2, hm450)))

## Filenames for the four principal models
model_files <- c(
  "Unweighted Euro"="unimputed_unweighted_%s_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv",
  "Weighted Euro"="unimputed_weighted_%s_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv",
  "Unweighted Transancestry"="unimputed_unweighted_%s_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv",
  "Weighted Transancestry"="unimputed_weighted_%s_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv"
)

summarize_models <- function(group, group_label) {
  probe_lists <- lapply(model_files, function(x) {
    fread(file.path(sig_dir, sprintf(x, group)), select="Probe_ID")[[1]]
  })
  
  combined <- rbindlist(Map(function(model, probes) {
    data.table(Model=model, Probe_ID=probes)
  }, names(probe_lists), probe_lists))
  
  combined[, CpG_ID := sub("_.*$", "", Probe_ID)]
  
  model_summary <- combined[, .(
    Unique_probes=uniqueN(Probe_ID),
    Unique_CpGs=uniqueN(CpG_ID),
    MSA_only_CpGs=uniqueN(CpG_ID[!CpG_ID %chin% known_cpgs])
  ), by=Model]
  
  pooled_summary <- combined[, .(
    Model="Across all 4 models",
    Unique_probes=uniqueN(Probe_ID),
    Unique_CpGs=uniqueN(CpG_ID),
    MSA_only_CpGs=uniqueN(CpG_ID[!CpG_ID %chin% known_cpgs])
  )]
  
  rbind(model_summary, pooled_summary)[, Smoking_group := group_label]
}

current_summary <- summarize_models("smokers", "Current")
former_summary <- summarize_models("formersmokers", "Former")

probe_summary <- rbind(current_summary, former_summary)
setcolorder(probe_summary, c("Smoking_group", "Model", "Unique_probes", "Unique_CpGs", "MSA_only_CpGs"))
all_results <- rbindlist(lapply(c("smokers", "formersmokers"), function(group) {
  rbindlist(lapply(model_files, function(x) {
    fread(file.path(sig_dir, sprintf(x, group)), select="Probe_ID")
  }))
}))

all_results[, CpG_ID := sub("_.*$", "", Probe_ID)]

total_unique <- all_results[, .(
  Smoking_group="Current and Former",
  Model="Across all 8 models",
  Unique_probes=uniqueN(Probe_ID),
  Unique_CpGs=uniqueN(CpG_ID),
  MSA_only_CpGs=uniqueN(CpG_ID[!CpG_ID %chin% known_cpgs])
)]

probe_summary <- rbind(probe_summary, total_unique)
print(probe_summary)

fwrite(
  probe_summary,
  file.path(sig_dir, "current_former_probe_cpg_msa_summary.tsv"),
  sep="\t"
)
