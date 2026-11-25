library(data.table)
library(sesame)
library(wheatmap)
library(multtest)
library(limma)
library(RColorBrewer)
library(ggplot2)
library(minfi)
library(ChAMP)
library(impute)
library(annotatr)
library(tidyr)
library(ggpubr)
library(dplyr)

bonf <- 0.05/230325

path <- ""
boxplot_dir <- file.path(path, "output", "boxplots")
dir.create(boxplot_dir, recursive=TRUE, showWarnings=FALSE)

load('751_clean_pc.rdata')

## Compatibility aliases for the current RData file
clean_betas <- betas
clean_betas_pval_masked <- betas_miss_masked
clean_sampsheet <- pheno

beta2 <- t(clean_betas)
betap <- t(clean_betas_pval_masked)
pheno <- clean_sampsheet

pheno_pc <- pheno %>%
  filter(European == "European", !is.na(PC1))

beta_pc <- beta2[
  rownames(beta2) %in% pheno_pc$Sample_ID,
  ,
  drop=FALSE
]

betap_pc <- betap[
  rownames(betap) %in% pheno_pc$Sample_ID,
  ,
  drop=FALSE
]

smoking_colors <- c(
  Never="#db9a9a",
  Former="#a33636",
  Current="#290101"
)

###### Model 1 #####################

siguws <- fread(
  file.path(
    path,
    "output",
    "svyglm_anno",
    "unimputed_unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv"
  )
)

siguws %>% arrange(P)

probes <- siguws %>%
  arrange(P) %>%
  filter(P <= bonf) %>%
  dplyr::select(V1)

pdf(
  file.path(
    boxplot_dir,
    "model1_significant_smoking_probes_europc_only.pdf"
  )
)

for (p in as.list(probes$V1)) {
  ahrr <- data.frame(beta_pc[, p])
  
  matching_rows <- grep(
    p,
    rownames(t(beta_pc)),
    value=TRUE
  )
  
  cpg_data <- t(beta_pc)[
    matching_rows,
    ,
    drop=FALSE
  ]
  
  cpg_long <- pivot_longer(
    as.data.frame(cpg_data),
    cols=everything(),
    names_to="Sample_ID",
    values_to="Beta_Value"
  )
  
  unique_samples <- unique(cpg_long$Sample_ID)
  
  ahrr <- left_join(
    cpg_long,
    pheno_pc,
    by="Sample_ID"
  )
  
  ahrr$LH_smoke_status <- factor(
    ahrr$LH_smoke_status,
    levels=c("Never", "Former", "Current", "M")
  )
  
  ahrr <- filter(
    ahrr,
    LH_smoke_status != "M"
  )
  
  gene <- ungroup(siguws) %>%
    dplyr::select(V1, genesUniq) %>%
    filter(V1 == p) %>%
    dplyr::select(genesUniq)
  
  print(
    ggboxplot(
      ahrr,
      x="LH_smoke_status",
      y="Beta_Value",
      color="LH_smoke_status",
      palette=smoking_colors,
      order=c("Never", "Former", "Current"),
      add="jitter",
      xlab="Smoking Status",
      ylab="Beta",
      title=paste0(
        "Probe: ",
        p,
        " Gene: ",
        gene,
        "\n by Smoking Status"
      )
    ) +
      theme(
        plot.title=element_text(hjust=0.5)
      )
  )
}

dev.off()

###### Model 3 #####################

siguws <- fread(
  file.path(
    path,
    "output",
    "svyglm_anno",
    "unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv"
  )
)

siguws %>% arrange(P)

probes <- siguws %>%
  arrange(P) %>%
  filter(P <= bonf) %>%
  dplyr::select(V1)

pdf(
  file.path(
    boxplot_dir,
    "model3_significant_smoking_probes_europc_only.pdf"
  )
)

for (p in as.list(probes$V1)) {
  ahrr <- data.frame(beta_pc[, p])
  
  matching_rows <- grep(
    p,
    rownames(t(beta_pc)),
    value=TRUE
  )
  
  cpg_data <- t(beta_pc)[
    matching_rows,
    ,
    drop=FALSE
  ]
  
  cpg_long <- pivot_longer(
    as.data.frame(cpg_data),
    cols=everything(),
    names_to="Sample_ID",
    values_to="Beta_Value"
  )
  
  unique_samples <- unique(cpg_long$Sample_ID)
  
  ahrr <- left_join(
    cpg_long,
    pheno_pc,
    by="Sample_ID"
  )
  
  ahrr$LH_smoke_status <- factor(
    ahrr$LH_smoke_status,
    levels=c("Never", "Former", "Current", "M")
  )
  
  ahrr <- filter(
    ahrr,
    LH_smoke_status != "M"
  )
  
  gene <- ungroup(siguws) %>%
    dplyr::select(V1, genesUniq) %>%
    filter(V1 == p) %>%
    dplyr::select(genesUniq)
  
  print(
    ggboxplot(
      ahrr,
      x="LH_smoke_status",
      y="Beta_Value",
      color="LH_smoke_status",
      palette=smoking_colors,
      order=c("Never", "Former", "Current"),
      add="jitter",
      xlab="Smoking Status",
      ylab="Beta",
      title=paste0(
        "Probe: ",
        p,
        " Gene: ",
        gene,
        "\n by Smoking Status"
      )
    ) +
      theme(
        plot.title=element_text(hjust=0.5)
      )
  )
}

dev.off()

###### Model 10 #####################

siguws <- fread(
  file.path(
    path,
    "output",
    "svyglm_anno",
    "unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv"
  )
)

siguws %>% arrange(P)

probes <- siguws %>%
  arrange(P) %>%
  filter(P <= bonf) %>%
  dplyr::select(V1)

pdf(
  file.path(
    boxplot_dir,
    "model10_significant_smoking_probes_transanc.pdf"
  )
)

for (p in as.list(probes$V1)) {
  ahrr <- data.frame(beta2[, p])
  
  matching_rows <- grep(
    p,
    rownames(t(beta2)),
    value=TRUE
  )
  
  cpg_data <- t(beta2)[
    matching_rows,
    ,
    drop=FALSE
  ]
  
  cpg_long <- pivot_longer(
    as.data.frame(cpg_data),
    cols=everything(),
    names_to="Sample_ID",
    values_to="Beta_Value"
  )
  
  unique_samples <- unique(cpg_long$Sample_ID)
  
  ahrr <- left_join(
    cpg_long,
    pheno,
    by="Sample_ID"
  )
  
  ahrr$LH_smoke_status <- factor(
    ahrr$LH_smoke_status,
    levels=c("Never", "Former", "Current", "M")
  )
  
  ahrr <- filter(
    ahrr,
    LH_smoke_status != "M"
  )
  
  gene <- ungroup(siguws) %>%
    dplyr::select(V1, genesUniq) %>%
    filter(V1 == p) %>%
    dplyr::select(genesUniq)
  
  print(
    ggboxplot(
      ahrr,
      x="LH_smoke_status",
      y="Beta_Value",
      color="LH_smoke_status",
      palette=smoking_colors,
      order=c("Never", "Former", "Current"),
      add="jitter",
      xlab="Smoking Status",
      ylab="Beta",
      title=paste0(
        "Probe: ",
        p,
        " Gene: ",
        gene,
        "\n by Smoking Status"
      )
    ) +
      theme(
        plot.title=element_text(hjust=0.5)
      )
  )
}

dev.off()

###### Model 12 #####################

siguws <- fread(
  file.path(
    path,
    "output",
    "svyglm_anno",
    "unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv"
  )
)

siguws %>% arrange(P)

probes <- siguws %>%
  arrange(P) %>%
  filter(P <= bonf) %>%
  dplyr::select(V1)

pdf(
  file.path(
    boxplot_dir,
    "model12_significant_smoking_probes_transanc.pdf"
  )
)

for (p in as.list(probes$V1)) {
  ahrr <- data.frame(beta2[, p])
  
  matching_rows <- grep(
    p,
    rownames(t(beta2)),
    value=TRUE
  )
  
  cpg_data <- t(beta2)[
    matching_rows,
    ,
    drop=FALSE
  ]
  
  cpg_long <- pivot_longer(
    as.data.frame(cpg_data),
    cols=everything(),
    names_to="Sample_ID",
    values_to="Beta_Value"
  )
  
  unique_samples <- unique(cpg_long$Sample_ID)
  
  ahrr <- left_join(
    cpg_long,
    pheno,
    by="Sample_ID"
  )
  
  ahrr$LH_smoke_status <- factor(
    ahrr$LH_smoke_status,
    levels=c("Never", "Former", "Current", "M")
  )
  
  ahrr <- filter(
    ahrr,
    LH_smoke_status != "M"
  )
  
  gene <- ungroup(siguws) %>%
    dplyr::select(V1, genesUniq) %>%
    filter(V1 == p) %>%
    dplyr::select(genesUniq)
  
  print(
    ggboxplot(
      ahrr,
      x="LH_smoke_status",
      y="Beta_Value",
      color="LH_smoke_status",
      palette=smoking_colors,
      order=c("Never", "Former", "Current"),
      add="jitter",
      xlab="Smoking Status",
      ylab="Beta",
      title=paste0(
        "Probe: ",
        p,
        " Gene: ",
        gene,
        "\n by Smoking Status"
      )
    ) +
      theme(
        plot.title=element_text(hjust=0.5)
      )
  )
}

dev.off()

###############################################################################






###### Model 1 #####################

siguws <- fread(
  file.path(
    path,
    "output",
    "svyglm_anno",
    "unimputed_unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv"
  )
)

siguws %>% arrange(P)

probes <- siguws %>%
  arrange(P) %>%
  filter(P <= bonf) %>%
  dplyr::select(V1)

pdf(
  file.path(
    boxplot_dir,
    "model1_significant_smoking_probes_europc_only_original_colors.pdf"
  )
)

for (p in as.list(probes$V1)) {
  ahrr <- data.frame(beta_pc[, p])
  
  matching_rows <- grep(
    p,
    rownames(t(beta_pc)),
    value=TRUE
  )
  
  cpg_data <- t(beta_pc)[
    matching_rows,
    ,
    drop=FALSE
  ]
  
  cpg_long <- pivot_longer(
    as.data.frame(cpg_data),
    cols=everything(),
    names_to="Sample_ID",
    values_to="Beta_Value"
  )
  
  unique_samples <- unique(cpg_long$Sample_ID)
  
  ahrr <- left_join(
    cpg_long,
    pheno_pc,
    by="Sample_ID"
  )
  
  ahrr$LH_smoke_status <- factor(
    ahrr$LH_smoke_status,
    levels=c("Never", "Former", "Current", "M")
  )
  
  ahrr <- filter(
    ahrr,
    LH_smoke_status != "M"
  )
  
  gene <- ungroup(siguws) %>%
    dplyr::select(V1, genesUniq) %>%
    filter(V1 == p) %>%
    dplyr::select(genesUniq)
  
  print(
    ggboxplot(
      ahrr,
      x="LH_smoke_status",
      y="Beta_Value",
      color="LH_smoke_status",
      palette="jama",
      order=c("Never", "Former", "Current"),
      add="jitter",
      xlab="Smoking Status",
      ylab="Beta",
      title=paste0(
        "Probe: ",
        p,
        " Gene: ",
        gene,
        "\n by Smoking Status"
      )
    ) +
      theme(
        plot.title=element_text(hjust=0.5)
      )
  )
}

dev.off()

###### Model 3 #####################

siguws <- fread(
  file.path(
    path,
    "output",
    "svyglm_anno",
    "unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv"
  )
)

siguws %>% arrange(P)

probes <- siguws %>%
  arrange(P) %>%
  filter(P <= bonf) %>%
  dplyr::select(V1)

pdf(
  file.path(
    boxplot_dir,
    "model3_significant_smoking_probes_europc_only_original_colors.pdf"
  )
)

for (p in as.list(probes$V1)) {
  ahrr <- data.frame(beta_pc[, p])
  
  matching_rows <- grep(
    p,
    rownames(t(beta_pc)),
    value=TRUE
  )
  
  cpg_data <- t(beta_pc)[
    matching_rows,
    ,
    drop=FALSE
  ]
  
  cpg_long <- pivot_longer(
    as.data.frame(cpg_data),
    cols=everything(),
    names_to="Sample_ID",
    values_to="Beta_Value"
  )
  
  unique_samples <- unique(cpg_long$Sample_ID)
  
  ahrr <- left_join(
    cpg_long,
    pheno_pc,
    by="Sample_ID"
  )
  
  ahrr$LH_smoke_status <- factor(
    ahrr$LH_smoke_status,
    levels=c("Never", "Former", "Current", "M")
  )
  
  ahrr <- filter(
    ahrr,
    LH_smoke_status != "M"
  )
  
  gene <- ungroup(siguws) %>%
    dplyr::select(V1, genesUniq) %>%
    filter(V1 == p) %>%
    dplyr::select(genesUniq)
  
  print(
    ggboxplot(
      ahrr,
      x="LH_smoke_status",
      y="Beta_Value",
      color="LH_smoke_status",
      palette="jama",
      order=c("Never", "Former", "Current"),
      add="jitter",
      xlab="Smoking Status",
      ylab="Beta",
      title=paste0(
        "Probe: ",
        p,
        " Gene: ",
        gene,
        "\n by Smoking Status"
      )
    ) +
      theme(
        plot.title=element_text(hjust=0.5)
      )
  )
}

dev.off()

###### Model 10 #####################

siguws <- fread(
  file.path(
    path,
    "output",
    "svyglm_anno",
    "unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv"
  )
)

siguws %>% arrange(P)

probes <- siguws %>%
  arrange(P) %>%
  filter(P <= bonf) %>%
  dplyr::select(V1)

pdf(
  file.path(
    boxplot_dir,
    "model10_significant_smoking_probes_transanc_original_colors.pdf"
  )
)

for (p in as.list(probes$V1)) {
  ahrr <- data.frame(beta2[, p])
  
  matching_rows <- grep(
    p,
    rownames(t(beta2)),
    value=TRUE
  )
  
  cpg_data <- t(beta2)[
    matching_rows,
    ,
    drop=FALSE
  ]
  
  cpg_long <- pivot_longer(
    as.data.frame(cpg_data),
    cols=everything(),
    names_to="Sample_ID",
    values_to="Beta_Value"
  )
  
  unique_samples <- unique(cpg_long$Sample_ID)
  
  ahrr <- left_join(
    cpg_long,
    pheno,
    by="Sample_ID"
  )
  
  ahrr$LH_smoke_status <- factor(
    ahrr$LH_smoke_status,
    levels=c("Never", "Former", "Current", "M")
  )
  
  ahrr <- filter(
    ahrr,
    LH_smoke_status != "M"
  )
  
  gene <- ungroup(siguws) %>%
    dplyr::select(V1, genesUniq) %>%
    filter(V1 == p) %>%
    dplyr::select(genesUniq)
  
  print(
    ggboxplot(
      ahrr,
      x="LH_smoke_status",
      y="Beta_Value",
      color="LH_smoke_status",
      palette="jama",
      order=c("Never", "Former", "Current"),
      add="jitter",
      xlab="Smoking Status",
      ylab="Beta",
      title=paste0(
        "Probe: ",
        p,
        " Gene: ",
        gene,
        "\n by Smoking Status"
      )
    ) +
      theme(
        plot.title=element_text(hjust=0.5)
      )
  )
}

dev.off()

###### Model 12 #####################

siguws <- fread(
  file.path(
    path,
    "output",
    "svyglm_anno",
    "unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv"
  )
)

siguws %>% arrange(P)

probes <- siguws %>%
  arrange(P) %>%
  filter(P <= bonf) %>%
  dplyr::select(V1)

pdf(
  file.path(
    boxplot_dir,
    "model12_significant_smoking_probes_transanc_original_colors.pdf"
  )
)

for (p in as.list(probes$V1)) {
  ahrr <- data.frame(beta2[, p])
  
  matching_rows <- grep(
    p,
    rownames(t(beta2)),
    value=TRUE
  )
  
  cpg_data <- t(beta2)[
    matching_rows,
    ,
    drop=FALSE
  ]
  
  cpg_long <- pivot_longer(
    as.data.frame(cpg_data),
    cols=everything(),
    names_to="Sample_ID",
    values_to="Beta_Value"
  )
  
  unique_samples <- unique(cpg_long$Sample_ID)
  
  ahrr <- left_join(
    cpg_long,
    pheno,
    by="Sample_ID"
  )
  
  ahrr$LH_smoke_status <- factor(
    ahrr$LH_smoke_status,
    levels=c("Never", "Former", "Current", "M")
  )
  
  ahrr <- filter(
    ahrr,
    LH_smoke_status != "M"
  )
  
  gene <- ungroup(siguws) %>%
    dplyr::select(V1, genesUniq) %>%
    filter(V1 == p) %>%
    dplyr::select(genesUniq)
  
  print(
    ggboxplot(
      ahrr,
      x="LH_smoke_status",
      y="Beta_Value",
      color="LH_smoke_status",
      palette="jama",
      order=c("Never", "Former", "Current"),
      add="jitter",
      xlab="Smoking Status",
      ylab="Beta",
      title=paste0(
        "Probe: ",
        p,
        " Gene: ",
        gene,
        "\n by Smoking Status"
      )
    ) +
      theme(
        plot.title=element_text(hjust=0.5)
      )
  )
}

dev.off()
