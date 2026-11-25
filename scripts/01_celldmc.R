library(EpiDISH)
library(data.table)
library(tidyverse)
library(readxl)
library(knitr)
library(sesame)
library(wheatmap)
library(multtest)
library(limma)
library(RColorBrewer)
library(ggplot2)
library(minfi)
library(ChAMP)
load('751_clean_pc.rdata')
betap <- t(betas_miss_masked)
mvalp <- t(mvals_miss_masked)

#pheno=pheno%>%mutate(LEL19=ifelse(LEL19=="Yes",1,ifelse(LEL19=="No",0,99)))
#pheno$LEL19=as.factor(pheno$LEL19)


pheno$LH_smoke_status=factor(pheno$LH_smoke_status, levels = c("Never","Former","Current","M"))
pheno$European=factor(pheno$European,levels = c("European","Non_European"))
pheno <- pheno %>% mutate( ANYETS = case_when(ANYETS =="No"~"No",ANYETS=="Yes"~"Yes",.default = "M"))
pheno$ANYETS <- factor(pheno$ANYETS,levels=c("No","Yes","M"))
pheno$LHSSCREEN_12142015 = factor(pheno$LHSSCREEN_12142015,levels=c("Control","Case"))
pheno$Gender <- factor(pheno$Gender,levels = c("male","female"))
pheno$STATE <- factor(pheno$STATE, levels=c('IA','NC'))

## Split the data into never former and never current
nf = pheno %>% filter(LH_smoke_status == "Never" | LH_smoke_status == "Former")
nc = pheno %>% filter(LH_smoke_status == "Never" | LH_smoke_status == "Current")
fc = pheno %>% filter(LH_smoke_status == "Former" | LH_smoke_status == "Current")

betas_nf = betas[,as.vector(nf$Sample_ID)]
betas_nc = betas[,as.vector(nc$Sample_ID)]
betas_fc = betas[,as.vector(fc$Sample_ID)]

ed_nf = nf %>% select(B, NK, CD4T, CD8T, Mono, Neutro)
ed_nc = nc %>% select(B, NK, CD4T, CD8T, Mono, Neutro)
ed_fc = fc %>% select(B, NK, CD4T, CD8T, Mono, Neutro)

cdmc_nf = CellDMC(beta.m=betas_nf, pheno.v=nf$LH_smoke_status,frac.m=as.matrix(ed_nf),
               cov.mod = model.matrix(~ Gender + AgeAtHomeVisit + STATE + ANYETS +
                                        LH_BMI + European + LHSSCREEN_12142015 +
                                        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
                                        PC8 + PC9 + PC10,data = nf),mc.cores=16)

cdmc_nc = CellDMC(beta.m=betas_nc, pheno.v=nc$LH_smoke_status,frac.m=as.matrix(ed_nc),
               cov.mod = model.matrix(~ Gender + AgeAtHomeVisit + STATE + ANYETS + 
                                        LH_BMI + European + LHSSCREEN_12142015 +
                                        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
                                        PC8 + PC9 + PC10,data = nc),mc.cores=16)

cdmc_fc = CellDMC(beta.m=betas_fc, pheno.v=fc$LH_smoke_status,frac.m=as.matrix(ed_fc),
                  cov.mod = model.matrix(~ Gender + AgeAtHomeVisit + STATE + ANYETS + 
                                           LH_BMI + European + LHSSCREEN_12142015 +
                                           PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
                                           PC8 + PC9 + PC10,data = fc),mc.cores=16)

#cdmc_nf$dmct %>% as_tibble(.,rownames="CpG") %>% filter(DMC == 1) 
#nccpg <- cdmc_nc$dmct %>% as_tibble(.,rownames="CpG") %>% filter(DMC == 1) %>% select(CpG)
#fccpg <- cdmc_fc$dmct %>% as_tibble(.,rownames="CpG") %>% filter(DMC == 1) %>% select(CpG)
#library(janitor)
# ed_nc = nc %>% select(Sample_ID,B, NK, CD4T, CD8T, Mono, Neutro, Eosino,LH_smoke_status)
# ed_fc = fc %>% select(Sample_ID,B, NK, CD4T, CD8T, Mono, Neutro, Eosino,LH_smoke_status)
# ncbetas <- as_tibble(subset(betas_nc,rownames(betas) %in% nccpg$CpG),rownames="CpG") %>% t %>%
#   as.data.frame %>% row_to_names(row_number = 1) 
# ncbetas <- left_join(rownames_to_column(ncbetas),ed_nc,by =c("rowname"= 'Sample_ID'))
# fcbetas <- as_tibble(subset(betas_fc,rownames(betas) %in% fccpg$CpG),rownames="CpG")%>% t %>%
#   as.data.frame %>% row_to_names(row_number = 1) 
# fcbetas <- left_join(rownames_to_column(fcbetas),ed_fc,by =c("rowname"= 'Sample_ID'))


## Ordinal Never=0, Former=1, Current=2 analysis
cell_types <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")

trend <- pheno %>%
  filter(LH_smoke_status %in% c("Never", "Former", "Current")) %>%
  mutate(smoke_trend=case_when(LH_smoke_status=="Never" ~ 0,
                               LH_smoke_status=="Former" ~ 1,
                               LH_smoke_status=="Current" ~ 2))
cov_trend <- model.matrix(
  ~ Gender + AgeAtHomeVisit + STATE + ANYETS + LH_BMI +
    European + LHSSCREEN_12142015 +
    PC1 + PC2 + PC3 + PC4 + PC5 +
    PC6 + PC7 + PC8 + PC9 + PC10,
  data=trend
)
cov_formula <- ~ Gender + AgeAtHomeVisit + STATE + ANYETS + LH_BMI + European +
  LHSSCREEN_12142015 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
betas_trend <- betas[, trend$Sample_ID, drop=FALSE]
ed_trend <- trend %>% dplyr::select(all_of(cell_types))
cov_trend <- model.matrix(cov_formula, data=trend)

cdmc_trend <- CellDMC(beta.m=betas_trend, pheno.v=trend$smoke_trend,
                           frac.m=as.matrix(ed_trend), cov.mod=cov_trend,
                           mc.cores=16)

trend_sig <- cdmc_trend$dmct %>%
  as_tibble(rownames="CpG") %>%
  filter(DMC == 1)

## Check monotonicity using the three post-ComBat pairwise models
trend_long <- trend_sig %>%
  pivot_longer(all_of(cell_types), names_to="Cell_type", values_to="Trend_direction") %>%
  filter(Trend_direction != 0)

get_cell_coef <- function(object, cpg, cell) {
  x <- object$coe[[cell]][cpg, , drop=FALSE]
  c(Estimate=unname(x$Estimate), P=unname(x$p))
}

dose_check <- bind_rows(lapply(seq_len(nrow(trend_long)), function(i) {
  cpg <- trend_long$CpG[i]
  cell <- trend_long$Cell_type[i]
  nf_coef <- get_cell_coef(cdmc_nf, cpg, cell)
  nc_coef <- get_cell_coef(cdmc_nc, cpg, cell)
  fc_coef <- get_cell_coef(cdmc_fc, cpg, cell)
  
  data.frame(CpG=cpg, Cell_type=cell,
             Trend_direction=trend_long$Trend_direction[i],
             Former_vs_Never=nf_coef["Estimate"], Former_vs_Never_P=nf_coef["P"],
             Current_vs_Never=nc_coef["Estimate"], Current_vs_Never_P=nc_coef["P"],
             Current_vs_Former=fc_coef["Estimate"], Current_vs_Former_P=fc_coef["P"])
}))

dose_check <- dose_check %>%
  mutate(monotonic_direction=case_when(
    Trend_direction == -1 & Former_vs_Never < 0 & Current_vs_Never < 0 & Current_vs_Former < 0 ~ "Never > Former > Current",
    Trend_direction == 1 & Former_vs_Never > 0 & Current_vs_Never > 0 & Current_vs_Former > 0 ~ "Never < Former < Current",
    TRUE ~ "Not consistently monotonic"
  ))

table(dose_check$monotonic_direction)
nrow(trend_sig)


###############################################################################
## Annotate and save CellDMC results
###############################################################################

library(data.table)
library(dplyr)
library(tidyr)
library(tibble)

cdmc_dir <- "output/cdmc"
dir.create(cdmc_dir, recursive=TRUE, showWarnings=FALSE)

annotation_file <- paste0(
  "references/",
  "MSA.hg38.manifest.gencode.v41.tsv.gz"
)

msa_annotation <- fread(annotation_file)

## Create and annotate significant-result tables
annotate_cdmc <- function(object) {
  object$dmct %>%
    as_tibble(rownames="CpG") %>%
    dplyr::filter(DMC == 1) %>%
    dplyr::left_join(msa_annotation, by=c("CpG"="probeID"))
}

nf_sig_anno <- annotate_cdmc(cdmc_nf)
nc_sig_anno <- annotate_cdmc(cdmc_nc)
fc_sig_anno <- annotate_cdmc(cdmc_fc)
trend_sig_anno <- annotate_cdmc(cdmc_trend)

## Annotate the monotonicity results
dose_check_anno <- dose_check %>%
  dplyr::left_join(msa_annotation, by=c("CpG"="probeID"))

## Summary by comparison
cdmc_summary <- data.frame(
  Comparison=c(
    "Never vs Former",
    "Never vs Current",
    "Former vs Current",
    "Ordinal trend"
  ),
  Probe_designs=c(
    nrow(nf_sig_anno),
    nrow(nc_sig_anno),
    nrow(fc_sig_anno),
    nrow(trend_sig_anno)
  ),
  Unique_CpGs=c(
    n_distinct(sub("_.*$", "", nf_sig_anno$CpG)),
    n_distinct(sub("_.*$", "", nc_sig_anno$CpG)),
    n_distinct(sub("_.*$", "", fc_sig_anno$CpG)),
    n_distinct(sub("_.*$", "", trend_sig_anno$CpG))
  )
)

## Summary of significant cell types and directions
cell_types_all <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")

summarize_cell_signals <- function(x, comparison) {
  available_cells <- intersect(cell_types_all, names(x))
  
  x %>%
    tidyr::pivot_longer(
      cols=all_of(available_cells),
      names_to="Cell_type",
      values_to="Direction"
    ) %>%
    dplyr::filter(Direction != 0) %>%
    dplyr::group_by(Cell_type, Direction) %>%
    dplyr::summarise(N=dplyr::n(), .groups="drop") %>%
    dplyr::mutate(Comparison=comparison, .before=1)
}

cell_signal_summary <- bind_rows(
  summarize_cell_signals(nf_sig_anno, "Never vs Former"),
  summarize_cell_signals(nc_sig_anno, "Never vs Current"),
  summarize_cell_signals(fc_sig_anno, "Former vs Current"),
  summarize_cell_signals(trend_sig_anno, "Ordinal trend")
)

## Save annotated significant tables
fwrite(
  nf_sig_anno,
  file.path(cdmc_dir, "never_vs_former_significant_annotated.tsv"),
  sep="\t",
  na=""
)

fwrite(
  nc_sig_anno,
  file.path(cdmc_dir, "never_vs_current_significant_annotated.tsv"),
  sep="\t",
  na=""
)

fwrite(
  fc_sig_anno,
  file.path(cdmc_dir, "former_vs_current_significant_annotated.tsv"),
  sep="\t",
  na=""
)

fwrite(
  trend_sig_anno,
  file.path(cdmc_dir, "ordinal_trend_significant_annotated.tsv"),
  sep="\t",
  na=""
)

fwrite(
  dose_check_anno,
  file.path(cdmc_dir, "ordinal_trend_monotonicity_annotated.tsv"),
  sep="\t",
  na=""
)

fwrite(
  cdmc_summary,
  file.path(cdmc_dir, "cdmc_significant_summary.tsv"),
  sep="\t",
  na=""
)

fwrite(
  cell_signal_summary,
  file.path(cdmc_dir, "cdmc_cell_type_summary.tsv"),
  sep="\t",
  na=""
)

## Save methylation tables used for downstream plots
fwrite(
  ncbetas,
  file.path(cdmc_dir, "never_vs_current_significant_betas.tsv"),
  sep="\t",
  na=""
)

fwrite(
  fcbetas,
  file.path(cdmc_dir, "former_vs_current_significant_betas.tsv"),
  sep="\t",
  na=""
)

## Save complete CellDMC objects and derived tables
save(
  cdmc_nf,
  cdmc_nc,
  cdmc_fc,
  cdmc_trend,
  nf_sig_anno,
  nc_sig_anno,
  fc_sig_anno,
  trend_sig_anno,
  dose_check_anno,
  cdmc_summary,
  cell_signal_summary,
  file=file.path(cdmc_dir, "cdmc_complete_results.rdata")
)

print(cdmc_summary)
print(cell_signal_summary)
table(dose_check$monotonic_direction)








library(dplyr)
library(tibble)
library(ggplot2)
library(cowplot)

cdmc_dir <- "output/cdmc"
dir.create(cdmc_dir, recursive=TRUE, showWarnings=FALSE)

## Exact significant DMC probe IDs
nccpg <- cdmc_nc$dmct %>%
  as_tibble(rownames="CpG") %>%
  filter(DMC == 1) %>%
  pull(CpG)

fccpg <- cdmc_fc$dmct %>%
  as_tibble(rownames="CpG") %>%
  filter(DMC == 1) %>%
  pull(CpG)

## Keep only probes present in the corresponding beta matrices
nccpg <- intersect(nccpg, rownames(betas_nc))
fccpg <- intersect(fccpg, rownames(betas_fc))

## Create plotting data
ed_nc <- nc %>%
  dplyr::select(Sample_ID, B, NK, CD4T, CD8T, Mono, Neutro, LH_smoke_status)

ed_fc <- fc %>%
  dplyr::select(Sample_ID, B, NK, CD4T, CD8T, Mono, Neutro, LH_smoke_status)

ncbetas <- as.data.frame(t(betas_nc[nccpg, , drop=FALSE])) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(ed_nc, by="Sample_ID")

fcbetas <- as.data.frame(t(betas_fc[fccpg, , drop=FALSE])) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(ed_fc, by="Sample_ID")

## Plot only CpGs identified as significant DMCs
plot_dmc_pdf <- function(plot_data, dmc_probes, output_file) {
  cell_types <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")
  
  pdf(output_file, width=12, height=7)
  
  for (probe in dmc_probes) {
    plots <- lapply(cell_types, function(cell) {
      ggplot(
        plot_data,
        aes(
          x=as.numeric(.data[[cell]]),
          y=as.numeric(.data[[probe]]),
          shape=LH_smoke_status,
          color=LH_smoke_status
        )
      ) +
        geom_point(alpha=0.8) +
        geom_smooth(
          aes(group=LH_smoke_status),
          method="loess",
          se=TRUE
        ) +
        labs(
          title=cell,
          x=paste(cell, "fraction"),
          y="Beta value",
          color="Smoking status",
          shape="Smoking status"
        ) +
        theme_classic(base_size=13) +
        theme(
          plot.title=element_text(
            face="bold",
            size=15,
            hjust=0.5
          ),
          legend.position="bottom"
        )
    })
    
    shared_legend <- cowplot::get_legend(
      plots[[1]] +
        theme(
          legend.position="bottom",
          legend.direction="horizontal"
        )
    )
    
    plots <- lapply(
      plots,
      function(p) p + theme(legend.position="none")
    )
    
    plot_grid_body <- cowplot::plot_grid(
      plotlist=plots,
      ncol=3,
      align="hv"
    )
    
    page_title <- cowplot::ggdraw() +
      cowplot::draw_label(
        probe,
        x=0.5,
        hjust=0.5,
        fontface="bold",
        size=18
      )
    
    complete_page <- cowplot::plot_grid(
      page_title,
      plot_grid_body,
      shared_legend,
      ncol=1,
      rel_heights=c(0.07, 1, 0.10)
    )
    
    print(complete_page)
  }
  
  dev.off()
}

plot_dmc_pdf(
  plot_data=ncbetas,
  dmc_probes=nccpg,
  output_file=file.path(cdmc_dir, "never_current_beta_dmc.pdf")
)

plot_dmc_pdf(
  plot_data=fcbetas,
  dmc_probes=fccpg,
  output_file=file.path(cdmc_dir, "former_current_beta_dmc.pdf")
)



###############################################################################


plot_dmc_pdf <- function(plot_data, dmc_probes, celldmc_result, output_file) {
  cell_types <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")
  
  cell_labels <- c(
    B="% B cells",
    NK="% NK cells",
    CD4T="% CD4 T cells",
    CD8T="% CD8 T cells",
    Mono="% monocytes",
    Neutro="% neutrophils"
  )
  
  if (is.data.frame(dmc_probes)) dmc_probes <- dmc_probes[[1]]
  dmc_probes <- unique(as.character(dmc_probes))
  
  pdf(output_file, width=12, height=7)
  
  for (probe in dmc_probes) {
    plots <- lapply(cell_types, function(cell) {
      cell_result <- celldmc_result$coe[[cell]]
      
      raw_p <- cell_result[probe, "p"]
      adjusted_p <- cell_result[probe, "adjP"]
      
      p_label <- paste0(
        "Smoking status × cell type: P=",
        format.pval(raw_p, digits=2, eps=1e-300),
        "; FDR=",
        format.pval(adjusted_p, digits=2, eps=1e-300)
      )
      
      ggplot(
        plot_data,
        aes(
          x=.data[[cell]],
          y=.data[[probe]],
          shape=LH_smoke_status,
          color=LH_smoke_status
        )
      ) +
        geom_point(alpha=0.8) +
        geom_smooth(
          aes(group=LH_smoke_status),
          method="loess",
          se=TRUE
        ) +
        scale_x_continuous(labels=scales::label_percent(accuracy=1)) +
        labs(
          title=cell_labels[[cell]],
          subtitle=p_label,
          x=cell_labels[[cell]],
          y="Beta value",
          color="Smoking status",
          shape="Smoking status"
        ) +
        theme_classic(base_size=13) +
        theme(
          plot.title=element_text(face="bold", size=15, hjust=0.5),
          plot.subtitle=element_text(size=9, hjust=0.5),
          legend.position="bottom"
        )
    })
    
    shared_legend <- cowplot::get_legend(
      plots[[1]] +
        theme(
          legend.position="bottom",
          legend.direction="horizontal"
        )
    )
    
    plots <- lapply(
      plots,
      function(p) p + theme(legend.position="none")
    )
    
    plot_grid_body <- cowplot::plot_grid(
      plotlist=plots,
      ncol=3,
      align="hv"
    )
    
    page_title <- cowplot::ggdraw() +
      cowplot::draw_label(
        probe,
        x=0.5,
        hjust=0.5,
        fontface="bold",
        size=18
      )
    
    complete_page <- cowplot::plot_grid(
      page_title,
      plot_grid_body,
      shared_legend,
      ncol=1,
      rel_heights=c(0.07, 1, 0.10)
    )
    
    print(complete_page)
  }
  
  dev.off()
}

cdmc_dir = 'output/cdmc/'
plot_dmc_pdf(
  plot_data=ncbetas,
  dmc_probes=nccpg,
  celldmc_result=cdmc_nc,
  output_file=file.path(cdmc_dir, "never_current_beta_dmc2.pdf")
)

plot_dmc_pdf(
  plot_data=fcbetas,
  dmc_probes=fccpg,
  celldmc_result=cdmc_fc,
  output_file=file.path(cdmc_dir, "former_current_beta_dmc2.pdf")
)


########################################################################
## Mean analysis
################################################################################
## CellDMC using mean methylation across probe designs for the same CpG

library(EpiDISH)
library(dplyr)
library(data.table)
library(tibble)
library(ggplot2)
library(cowplot)

project_dir <- ""
cdmc_mean_dir <- file.path(project_dir, "output", "cdmc_mean")
dir.create(cdmc_mean_dir, recursive=TRUE, showWarnings=FALSE)

load('751_clean_pc.rdata')



################################################################################
## Prepare phenotype variables exactly as in the original analysis

pheno$LH_smoke_status <- factor(
  pheno$LH_smoke_status,
  levels=c("Never", "Former", "Current", "M")
)

pheno$European <- factor(
  pheno$European,
  levels=c("European", "Non_European")
)

pheno <- pheno %>%
  mutate(
    ANYETS=case_when(
      ANYETS == "No" ~ "No",
      ANYETS == "Yes" ~ "Yes",
      TRUE ~ "M"
    )
  )

pheno$ANYETS <- factor(pheno$ANYETS, levels=c("No", "Yes", "M"))
pheno$LHSSCREEN_12142015 <- factor(pheno$LHSSCREEN_12142015, levels=c("Control", "Case"))
pheno$Gender <- factor(pheno$Gender, levels=c("male", "female"))
pheno$STATE <- factor(pheno$STATE, levels=c("IA", "NC"))

stopifnot(identical(colnames(betas), pheno$Sample_ID))

################################################################################
## Average probe designs belonging to the same CpG

average_msa_designs <- function(beta_matrix) {
  cpg_id <- sub("_.*$", "", rownames(beta_matrix))
  cpg_order <- unique(cpg_id)
  
  beta_sum <- rowsum(
    beta_matrix,
    group=factor(cpg_id, levels=cpg_order),
    reorder=FALSE,
    na.rm=TRUE
  )
  
  n_observed <- rowsum(
    1L * !is.na(beta_matrix),
    group=factor(cpg_id, levels=cpg_order),
    reorder=FALSE
  )
  
  beta_mean <- beta_sum/n_observed
  beta_mean[n_observed == 0] <- NA_real_
  rownames(beta_mean) <- cpg_order
  
  return(beta_mean)
}

betas_mean <- average_msa_designs(betas)

dim(betas_mean)
################################################################################
## Split into Never–Former, Never–Current, and Former–Current

nf_mean <- pheno %>%
  filter(LH_smoke_status %in% c("Never", "Former")) %>%
  droplevels()

nc_mean <- pheno %>%
  filter(LH_smoke_status %in% c("Never", "Current")) %>%
  droplevels()

fc_mean <- pheno %>%
  filter(LH_smoke_status %in% c("Former", "Current")) %>%
  droplevels()

nf_mean$LH_smoke_status <- factor(nf_mean$LH_smoke_status, levels=c("Never", "Former"))
nc_mean$LH_smoke_status <- factor(nc_mean$LH_smoke_status, levels=c("Never", "Current"))
fc_mean$LH_smoke_status <- factor(fc_mean$LH_smoke_status, levels=c("Former", "Current"))

betas_nf_mean <- betas_mean[, nf_mean$Sample_ID, drop=FALSE]
betas_nc_mean <- betas_mean[, nc_mean$Sample_ID, drop=FALSE]
betas_fc_mean <- betas_mean[, fc_mean$Sample_ID, drop=FALSE]

cell_types <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")

frac_nf_mean <- nf_mean %>% dplyr::select(all_of(cell_types))
frac_nc_mean <- nc_mean %>% dplyr::select(all_of(cell_types))
frac_fc_mean <- fc_mean %>% dplyr::select(all_of(cell_types))

cov_formula <- ~ Gender + AgeAtHomeVisit + STATE + ANYETS +
  LH_BMI + European + LHSSCREEN_12142015 +
  PC1 + PC2 + PC3 + PC4 + PC5 +
  PC6 + PC7 + PC8 + PC9 + PC10

cov_nf_mean <- model.matrix(cov_formula, data=nf_mean)
cov_nc_mean <- model.matrix(cov_formula, data=nc_mean)
cov_fc_mean <- model.matrix(cov_formula, data=fc_mean)

stopifnot(
  identical(colnames(betas_nf_mean), nf_mean$Sample_ID),
  identical(colnames(betas_nc_mean), nc_mean$Sample_ID),
  identical(colnames(betas_fc_mean), fc_mean$Sample_ID),
  ncol(betas_nf_mean) == nrow(frac_nf_mean),
  ncol(betas_nc_mean) == nrow(frac_nc_mean),
  ncol(betas_fc_mean) == nrow(frac_fc_mean),
  nrow(cov_nf_mean) == ncol(betas_nf_mean),
  nrow(cov_nc_mean) == ncol(betas_nc_mean),
  nrow(cov_fc_mean) == ncol(betas_fc_mean)
)

################################################################################
## Run CellDMC

cdmc_nf_mean <- CellDMC(
  beta.m=betas_nf_mean,
  pheno.v=nf_mean$LH_smoke_status,
  frac.m=as.matrix(frac_nf_mean),
  cov.mod=cov_nf_mean,
  adjPMethod="fdr",
  adjPThresh=0.05,
  mc.cores=16
)

cdmc_nc_mean <- CellDMC(
  beta.m=betas_nc_mean,
  pheno.v=nc_mean$LH_smoke_status,
  frac.m=as.matrix(frac_nc_mean),
  cov.mod=cov_nc_mean,
  adjPMethod="fdr",
  adjPThresh=0.05,
  mc.cores=16
)

cdmc_fc_mean <- CellDMC(
  beta.m=betas_fc_mean,
  pheno.v=fc_mean$LH_smoke_status,
  frac.m=as.matrix(frac_fc_mean),
  cov.mod=cov_fc_mean,
  adjPMethod="fdr",
  adjPThresh=0.05,
  mc.cores=16
)

################################################################################
## Extract significant CpGs

nf_sig_mean <- cdmc_nf_mean$dmct %>%
  as.data.frame() %>%
  rownames_to_column("CpG") %>%
  filter(DMC == 1)

nc_sig_mean <- cdmc_nc_mean$dmct %>%
  as.data.frame() %>%
  rownames_to_column("CpG") %>%
  filter(DMC == 1)

fc_sig_mean <- cdmc_fc_mean$dmct %>%
  as.data.frame() %>%
  rownames_to_column("CpG") %>%
  filter(DMC == 1)

nfcpg_mean <- nf_sig_mean$CpG
nccpg_mean <- nc_sig_mean$CpG
fccpg_mean <- fc_sig_mean$CpG

data.frame(
  Comparison=c("Never vs Former", "Never vs Current", "Former vs Current"),
  Significant_CpGs=c(
    length(nfcpg_mean),
    length(nccpg_mean),
    length(fccpg_mean)
  )
)

fwrite(nf_sig_mean, file.path(cdmc_mean_dir, "never_vs_former_mean_significant.tsv"), sep="\t")
fwrite(nc_sig_mean, file.path(cdmc_mean_dir, "never_vs_current_mean_significant.tsv"), sep="\t")
fwrite(fc_sig_mean, file.path(cdmc_mean_dir, "former_vs_current_mean_significant.tsv"), sep="\t")

################################################################################
## Create plotting datasets
## These contain methylation, cell fractions, and smoking status

make_plot_data <- function(beta_matrix, pheno_data, probes) {
  beta_data <- as.data.frame(
    t(beta_matrix[probes, , drop=FALSE]),
    check.names=FALSE
  )
  
  beta_data <- rownames_to_column(beta_data, var="Sample_ID")
  
  plot_metadata <- pheno_data %>%
    dplyr::select(
      Sample_ID,
      all_of(cell_types),
      LH_smoke_status
    )
  
  left_join(beta_data, plot_metadata, by="Sample_ID")
}

nfbetas_mean <- make_plot_data(betas_nf_mean, nf_mean, nfcpg_mean)
ncbetas_mean <- make_plot_data(betas_nc_mean, nc_mean, nccpg_mean)
fcbetas_mean <- make_plot_data(betas_fc_mean, fc_mean, fccpg_mean)

stopifnot(
  all(cell_types %in% names(nfbetas_mean)),
  all(cell_types %in% names(ncbetas_mean)),
  all(cell_types %in% names(fcbetas_mean)),
  all(nfcpg_mean %in% names(nfbetas_mean)),
  all(nccpg_mean %in% names(ncbetas_mean)),
  all(fccpg_mean %in% names(fcbetas_mean))
)

################################################################################
## Plot significant mean CpGs and report interaction FDR and raw P value

plot_mean_dmc_pdf <- function(plot_data, dmc_probes, celldmc_result, output_file) {
  cell_types <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")
  
  cell_labels <- c(
    B="% B cells",
    NK="% NK cells",
    CD4T="% CD4 T cells",
    CD8T="% CD8 T cells",
    Mono="% monocytes",
    Neutro="% neutrophils"
  )
  
  if (is.data.frame(dmc_probes)) dmc_probes <- dmc_probes[[1]]
  dmc_probes <- unique(as.character(dmc_probes))
  
  pdf(output_file, width=12, height=7)
  
  for (probe in dmc_probes) {
    plots <- lapply(cell_types, function(cell) {
      cell_result <- celldmc_result$coe[[cell]]
      raw_p <- cell_result[probe, "p"]
      adjusted_p <- cell_result[probe, "adjP"]
      
      p_label <- paste0(
        "Smoking status × cell type: P=",
        format.pval(raw_p, digits=2, eps=1e-300),
        "; FDR=",
        format.pval(adjusted_p, digits=2, eps=1e-300)
      )
      
      ggplot(
        plot_data,
        aes(
          x=.data[[cell]],
          y=.data[[probe]],
          shape=LH_smoke_status,
          color=LH_smoke_status
        )
      ) +
        geom_point(alpha=0.8) +
        geom_smooth(
          aes(group=LH_smoke_status),
          method="loess",
          se=TRUE
        ) +
        scale_x_continuous(labels=scales::label_percent(accuracy=1)) +
        labs(
          title=cell_labels[[cell]],
          subtitle=p_label,
          x=cell_labels[[cell]],
          y="Beta value",
          color="Smoking status",
          shape="Smoking status"
        ) +
        theme_classic(base_size=13) +
        theme(
          plot.title=element_text(face="bold", size=15, hjust=0.5),
          plot.subtitle=element_text(size=9, hjust=0.5),
          legend.position="bottom"
        )
    })
    
    shared_legend <- cowplot::get_legend(
      plots[[1]] +
        theme(
          legend.position="bottom",
          legend.direction="horizontal"
        )
    )
    
    plots <- lapply(
      plots,
      function(p) p + theme(legend.position="none")
    )
    
    plot_grid_body <- cowplot::plot_grid(
      plotlist=plots,
      ncol=3,
      align="hv"
    )
    
    page_title <- cowplot::ggdraw() +
      cowplot::draw_label(
        probe,
        x=0.5,
        hjust=0.5,
        fontface="bold",
        size=18
      )
    
    complete_page <- cowplot::plot_grid(
      page_title,
      plot_grid_body,
      shared_legend,
      ncol=1,
      rel_heights=c(0.07, 1, 0.10)
    )
    
    print(complete_page)
  }
  
  dev.off()
}
################################################################################
## Generate PDFs and interaction-statistic tables

plot_mean_dmc_pdf(
  plot_data=nfbetas_mean,
  dmc_probes=nfcpg_mean,
  celldmc_result=cdmc_nf_mean,
  output_file=file.path(cdmc_mean_dir, "never_former_mean_cpg_beta_dmc.pdf")
)

plot_mean_dmc_pdf(
  plot_data=ncbetas_mean,
  dmc_probes=nccpg_mean,
  celldmc_result=cdmc_nc_mean,
  output_file=file.path(cdmc_mean_dir, "never_current_mean_cpg_beta_dmc.pdf")
)

plot_mean_dmc_pdf(
  plot_data=fcbetas_mean,
  dmc_probes=fccpg_mean,
  celldmc_result=cdmc_fc_mean,
  output_file=file.path(cdmc_mean_dir, "former_current_mean_cpg_beta_dmc.pdf")
)

################################################################################
## Save the complete mean-CpG analysis

save(
  betas_mean,
  cdmc_nf_mean,
  cdmc_nc_mean,
  cdmc_fc_mean,
  nf_sig_mean,
  nc_sig_mean,
  fc_sig_mean,
  nf_interaction_mean,
  nc_interaction_mean,
  fc_interaction_mean,
  file=file.path(cdmc_mean_dir, "cdmc_mean_cpg_complete_results.rdata")
)




########################################################################
## Median analysis
################################################################################
## CellDMC using median methylation across probe designs for the same CpG

library(EpiDISH)
library(dplyr)
library(data.table)
library(tibble)
library(ggplot2)
library(cowplot)

project_dir <- ""
cdmc_median_dir <- file.path(project_dir, "output", "cdmc_median")
dir.create(cdmc_median_dir, recursive=TRUE, showWarnings=FALSE)

load('751_clean_pc.rdata')

################################################################################
## Prepare phenotype variables exactly as in the original analysis

pheno$LH_smoke_status <- factor(
  pheno$LH_smoke_status,
  levels=c("Never", "Former", "Current", "M")
)

pheno$European <- factor(
  pheno$European,
  levels=c("European", "Non_European")
)

pheno <- pheno %>%
  mutate(
    ANYETS=case_when(
      ANYETS == "No" ~ "No",
      ANYETS == "Yes" ~ "Yes",
      TRUE ~ "M"
    )
  )

pheno$ANYETS <- factor(pheno$ANYETS, levels=c("No", "Yes", "M"))
pheno$LHSSCREEN_12142015 <- factor(pheno$LHSSCREEN_12142015, levels=c("Control", "Case"))
pheno$Gender <- factor(pheno$Gender, levels=c("male", "female"))
pheno$STATE <- factor(pheno$STATE, levels=c("IA", "NC"))

stopifnot(identical(colnames(betas), pheno$Sample_ID))

################################################################################
## Calculate median across probe designs belonging to the same CpG

median_msa_designs <- function(beta_matrix) {
  cpg_id <- sub("_.*$", "", rownames(beta_matrix))
  cpg_order <- unique(cpg_id)
  
  first_design <- match(cpg_order, cpg_id)
  beta_median <- beta_matrix[first_design, , drop=FALSE]
  rownames(beta_median) <- cpg_order
  
  n_designs <- tabulate(
    match(cpg_id, cpg_order),
    nbins=length(cpg_order)
  )
  
  multiple_design_cpgs <- cpg_order[n_designs > 1]
  multiple_rows <- which(cpg_id %in% multiple_design_cpgs)
  rows_by_cpg <- split(multiple_rows, cpg_id[multiple_rows])
  
  for (cpg in names(rows_by_cpg)) {
    beta_median[cpg, ] <- apply(
      beta_matrix[rows_by_cpg[[cpg]], , drop=FALSE],
      2,
      median,
      na.rm=TRUE
    )
  }
  
  return(beta_median)
}

betas_median <- median_msa_designs(betas)

dim(betas_median)

################################################################################
## Split into Never–Former, Never–Current, and Former–Current

nf_median <- pheno %>%
  filter(LH_smoke_status %in% c("Never", "Former")) %>%
  droplevels()

nc_median <- pheno %>%
  filter(LH_smoke_status %in% c("Never", "Current")) %>%
  droplevels()

fc_median <- pheno %>%
  filter(LH_smoke_status %in% c("Former", "Current")) %>%
  droplevels()

nf_median$LH_smoke_status <- factor(nf_median$LH_smoke_status, levels=c("Never", "Former"))
nc_median$LH_smoke_status <- factor(nc_median$LH_smoke_status, levels=c("Never", "Current"))
fc_median$LH_smoke_status <- factor(fc_median$LH_smoke_status, levels=c("Former", "Current"))

betas_nf_median <- betas_median[, nf_median$Sample_ID, drop=FALSE]
betas_nc_median <- betas_median[, nc_median$Sample_ID, drop=FALSE]
betas_fc_median <- betas_median[, fc_median$Sample_ID, drop=FALSE]

cell_types <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")

frac_nf_median <- nf_median %>% dplyr::select(all_of(cell_types))
frac_nc_median <- nc_median %>% dplyr::select(all_of(cell_types))
frac_fc_median <- fc_median %>% dplyr::select(all_of(cell_types))

cov_formula <- ~ Gender + AgeAtHomeVisit + STATE + ANYETS +
  LH_BMI + European + LHSSCREEN_12142015 +
  PC1 + PC2 + PC3 + PC4 + PC5 +
  PC6 + PC7 + PC8 + PC9 + PC10

cov_nf_median <- model.matrix(cov_formula, data=nf_median)
cov_nc_median <- model.matrix(cov_formula, data=nc_median)
cov_fc_median <- model.matrix(cov_formula, data=fc_median)

stopifnot(
  identical(colnames(betas_nf_median), nf_median$Sample_ID),
  identical(colnames(betas_nc_median), nc_median$Sample_ID),
  identical(colnames(betas_fc_median), fc_median$Sample_ID),
  ncol(betas_nf_median) == nrow(frac_nf_median),
  ncol(betas_nc_median) == nrow(frac_nc_median),
  ncol(betas_fc_median) == nrow(frac_fc_median),
  nrow(cov_nf_median) == ncol(betas_nf_median),
  nrow(cov_nc_median) == ncol(betas_nc_median),
  nrow(cov_fc_median) == ncol(betas_fc_median)
)

################################################################################
## Run CellDMC

cdmc_nf_median <- CellDMC(
  beta.m=betas_nf_median,
  pheno.v=nf_median$LH_smoke_status,
  frac.m=as.matrix(frac_nf_median),
  cov.mod=cov_nf_median,
  adjPMethod="fdr",
  adjPThresh=0.05,
  mc.cores=16
)

cdmc_nc_median <- CellDMC(
  beta.m=betas_nc_median,
  pheno.v=nc_median$LH_smoke_status,
  frac.m=as.matrix(frac_nc_median),
  cov.mod=cov_nc_median,
  adjPMethod="fdr",
  adjPThresh=0.05,
  mc.cores=16
)

cdmc_fc_median <- CellDMC(
  beta.m=betas_fc_median,
  pheno.v=fc_median$LH_smoke_status,
  frac.m=as.matrix(frac_fc_median),
  cov.mod=cov_fc_median,
  adjPMethod="fdr",
  adjPThresh=0.05,
  mc.cores=16
)

################################################################################
## Extract significant CpGs

nf_sig_median <- cdmc_nf_median$dmct %>%
  as.data.frame() %>%
  rownames_to_column("CpG") %>%
  filter(DMC == 1)

nc_sig_median <- cdmc_nc_median$dmct %>%
  as.data.frame() %>%
  rownames_to_column("CpG") %>%
  filter(DMC == 1)

fc_sig_median <- cdmc_fc_median$dmct %>%
  as.data.frame() %>%
  rownames_to_column("CpG") %>%
  filter(DMC == 1)

nfcpg_median <- nf_sig_median$CpG
nccpg_median <- nc_sig_median$CpG
fccpg_median <- fc_sig_median$CpG

data.frame(
  Comparison=c("Never vs Former", "Never vs Current", "Former vs Current"),
  Significant_CpGs=c(
    length(nfcpg_median),
    length(nccpg_median),
    length(fccpg_median)
  )
)

fwrite(nf_sig_median, file.path(cdmc_median_dir, "never_vs_former_median_significant.tsv"), sep="\t")
fwrite(nc_sig_median, file.path(cdmc_median_dir, "never_vs_current_median_significant.tsv"), sep="\t")
fwrite(fc_sig_median, file.path(cdmc_median_dir, "former_vs_current_median_significant.tsv"), sep="\t")

################################################################################
## Create plotting datasets
## These contain methylation, cell fractions, and smoking status

make_plot_data <- function(beta_matrix, pheno_data, probes) {
  beta_data <- as.data.frame(
    t(beta_matrix[probes, , drop=FALSE]),
    check.names=FALSE
  )
  
  beta_data <- rownames_to_column(beta_data, var="Sample_ID")
  
  plot_metadata <- pheno_data %>%
    dplyr::select(
      Sample_ID,
      all_of(cell_types),
      LH_smoke_status
    )
  
  left_join(beta_data, plot_metadata, by="Sample_ID")
}

nfbetas_median <- make_plot_data(betas_nf_median, nf_median, nfcpg_median)
ncbetas_median <- make_plot_data(betas_nc_median, nc_median, nccpg_median)
fcbetas_median <- make_plot_data(betas_fc_median, fc_median, fccpg_median)

stopifnot(
  all(cell_types %in% names(nfbetas_median)),
  all(cell_types %in% names(ncbetas_median)),
  all(cell_types %in% names(fcbetas_median)),
  all(nfcpg_median %in% names(nfbetas_median)),
  all(nccpg_median %in% names(ncbetas_median)),
  all(fccpg_median %in% names(fcbetas_median))
)

################################################################################
## Plot significant median CpGs and report interaction FDR and raw P value

plot_median_dmc_pdf <- function(plot_data, dmc_probes, celldmc_result, output_file) {
  cell_types <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")
  
  cell_labels <- c(
    B="% B cells",
    NK="% NK cells",
    CD4T="% CD4 T cells",
    CD8T="% CD8 T cells",
    Mono="% monocytes",
    Neutro="% neutrophils"
  )
  
  if (is.data.frame(dmc_probes)) dmc_probes <- dmc_probes[[1]]
  dmc_probes <- unique(as.character(dmc_probes))
  
  pdf(output_file, width=12, height=7)
  
  for (probe in dmc_probes) {
    plots <- lapply(cell_types, function(cell) {
      cell_result <- celldmc_result$coe[[cell]]
      raw_p <- cell_result[probe, "p"]
      adjusted_p <- cell_result[probe, "adjP"]
      
      p_label <- paste0(
        "Smoking status × cell type: P=",
        format.pval(raw_p, digits=2, eps=1e-300),
        "; FDR=",
        format.pval(adjusted_p, digits=2, eps=1e-300)
      )
      
      ggplot(
        plot_data,
        aes(
          x=.data[[cell]],
          y=.data[[probe]],
          shape=LH_smoke_status,
          color=LH_smoke_status
        )
      ) +
        geom_point(alpha=0.8) +
        geom_smooth(
          aes(group=LH_smoke_status),
          method="loess",
          se=TRUE
        ) +
        scale_x_continuous(labels=scales::label_percent(accuracy=1)) +
        labs(
          title=cell_labels[[cell]],
          subtitle=p_label,
          x=cell_labels[[cell]],
          y="Beta value",
          color="Smoking status",
          shape="Smoking status"
        ) +
        theme_classic(base_size=13) +
        theme(
          plot.title=element_text(face="bold", size=15, hjust=0.5),
          plot.subtitle=element_text(size=9, hjust=0.5),
          legend.position="bottom"
        )
    })
    
    shared_legend <- cowplot::get_legend(
      plots[[1]] +
        theme(
          legend.position="bottom",
          legend.direction="horizontal"
        )
    )
    
    plots <- lapply(
      plots,
      function(p) p + theme(legend.position="none")
    )
    
    plot_grid_body <- cowplot::plot_grid(
      plotlist=plots,
      ncol=3,
      align="hv"
    )
    
    page_title <- cowplot::ggdraw() +
      cowplot::draw_label(
        probe,
        x=0.5,
        hjust=0.5,
        fontface="bold",
        size=18
      )
    
    complete_page <- cowplot::plot_grid(
      page_title,
      plot_grid_body,
      shared_legend,
      ncol=1,
      rel_heights=c(0.07, 1, 0.10)
    )
    
    print(complete_page)
  }
  
  dev.off()
}

################################################################################
## Generate PDFs

plot_median_dmc_pdf(
  plot_data=nfbetas_median,
  dmc_probes=nfcpg_median,
  celldmc_result=cdmc_nf_median,
  output_file=file.path(cdmc_median_dir, "never_former_median_cpg_beta_dmc.pdf")
)

plot_median_dmc_pdf(
  plot_data=ncbetas_median,
  dmc_probes=nccpg_median,
  celldmc_result=cdmc_nc_median,
  output_file=file.path(cdmc_median_dir, "never_current_median_cpg_beta_dmc.pdf")
)

plot_median_dmc_pdf(
  plot_data=fcbetas_median,
  dmc_probes=fccpg_median,
  celldmc_result=cdmc_fc_median,
  output_file=file.path(cdmc_median_dir, "former_current_median_cpg_beta_dmc.pdf")
)



###############################################################################
## Barchart
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

root <- "output"

cell_levels <- c(
  "B",
  "CD4T",
  "CD8T",
  "Mono",
  "Neutro",
  "NK",
  "Eosino"
)

cell_labels <- c(
  B="B",
  CD4T="CD4T",
  CD8T="CD8T",
  Mono="Monocytes",
  Neutro="Neutrophils",
  NK="NK",
  Eosino="Eosinophils"
)

fill_colors <- c(
  "Negative current vs. former"="#F8766D",
  "Negative current vs. never"="#7CAE00",
  "Positive current vs. former"="#00BFC4",
  "Positive current vs. never"="#C77CFF"
)

read_dmct_counts <- function(current_never_file, current_former_file, analysis) {
  current_never <- fread(current_never_file)
  current_former <- fread(current_former_file)
  
  current_never_long <- current_never %>%
    dplyr::select(any_of(cell_levels)) %>%
    tidyr::pivot_longer(
      cols=everything(),
      names_to="Cell_type",
      values_to="Direction"
    ) %>%
    filter(Direction != 0) %>%
    mutate(
      Direction_contrast=ifelse(
        Direction == -1,
        "Negative current vs. never",
        "Positive current vs. never"
      )
    )
  
  current_former_long <- current_former %>%
    dplyr::select(any_of(cell_levels)) %>%
    tidyr::pivot_longer(
      cols=everything(),
      names_to="Cell_type",
      values_to="Direction"
    ) %>%
    filter(Direction != 0) %>%
    mutate(
      Direction_contrast=ifelse(
        Direction == -1,
        "Negative current vs. former",
        "Positive current vs. former"
      )
    )
  
  bind_rows(current_never_long, current_former_long) %>%
    dplyr::count(
      Cell_type,
      Direction,
      Direction_contrast,
      name="n"
    ) %>%
    mutate(
      Analysis=analysis,
      Cell_type=factor(Cell_type, levels=cell_levels)
    )
}

plot_dmct_counts <- function(counts, title, output_file) {
  p <- ggplot(
    counts,
    aes(
      x=Cell_type,
      y=n,
      fill=Direction_contrast
    )
  ) +
    geom_col(
      position=position_dodge(width=0.8),
      width=0.7
    ) +
    scale_x_discrete(
      limits=cell_levels,
      labels=cell_labels,
      drop=FALSE
    ) +
    scale_fill_manual(
      values=fill_colors,
      drop=FALSE
    ) +
    labs(
      title=title,
      x="Cell Type",
      y="Number of DMCTs",
      fill="Direction and contrast"
    ) +
    theme_minimal() +
    theme(
      plot.title=element_text(hjust=0.5),
      axis.text.x=element_text(angle=45, hjust=1),
      panel.grid.minor=element_blank()
    )
  
  ggsave(
    output_file,
    plot=p,
    width=8,
    height=5
  )
  
  return(p)
}

################################################################################
## Probe-design results

probe_counts <- read_dmct_counts(
  current_never_file=file.path(
    root,
    "cdmc",
    "never_vs_current_significant_annotated.tsv"
  ),
  current_former_file=file.path(
    root,
    "cdmc",
    "former_vs_current_significant_annotated.tsv"
  ),
  analysis="Probe design"
)

probe_plot <- plot_dmct_counts(
  counts=probe_counts,
  title="Barplot of DMCTs by Cell Type — Probe Design",
  output_file=file.path(
    root,
    "cdmc",
    "dmct_barplot_by_probe.pdf"
  )
)

################################################################################
## Mean-CpG results

mean_counts <- read_dmct_counts(
  current_never_file=file.path(
    root,
    "cdmc_mean",
    "never_vs_current_mean_significant.tsv"
  ),
  current_former_file=file.path(
    root,
    "cdmc_mean",
    "former_vs_current_mean_significant.tsv"
  ),
  analysis="Mean CpG"
)

mean_plot <- plot_dmct_counts(
  counts=mean_counts,
  title="Barplot of DMCTs by Cell Type — Mean CpG",
  output_file=file.path(
    root,
    "cdmc_mean",
    "dmct_barplot_by_mean_cpg.pdf"
  )
)

################################################################################
## Median-CpG results

median_counts <- read_dmct_counts(
  current_never_file=file.path(
    root,
    "cdmc_median",
    "never_vs_current_median_significant.tsv"
  ),
  current_former_file=file.path(
    root,
    "cdmc_median",
    "former_vs_current_median_significant.tsv"
  ),
  analysis="Median CpG"
)

median_plot <- plot_dmct_counts(
  counts=median_counts,
  title="Barplot of DMCTs by Cell Type — Median CpG",
  output_file=file.path(
    root,
    "cdmc_median",
    "dmct_barplot_by_median_cpg.pdf"
  )
)

################################################################################
## Save and display counts

all_dmct_counts <- bind_rows(
  probe_counts,
  mean_counts,
  median_counts
) %>%
  arrange(Analysis, Direction_contrast, Cell_type)

fwrite(
  all_dmct_counts,
  file.path(root, "cdmc", "dmct_counts_probe_mean_median.tsv"),
  sep="\t"
)

print(all_dmct_counts, n=Inf)
















#####################
# Adjusting beta values for each CpG by cell proportions


betas_nc <- as_tibble(betas_nc,rownames="CpG")
# Load necessary libraries
library(dplyr)

# Create the beta values data frame
beta_values <- data.frame(
  CpG = c("cg00000109_BC21", "cg00000658_TC21", "cg00000714_TC21", 
          "cg00000721_TC21", "cg00001099_BC21"),
  A1 = c(0.8435673, 0.8267204, 0.2360347, 0.9349584, 0.6986932),
  A2 = c(0.8497762, 0.8516334, 0.2540831, 0.9464118, 0.8425044)
)

# Create the cell proportions data frame
cell_proportions <- data.frame(
  Sample_ID = c("A1", "A2"),
  B = c(0.0329, 0.0396),
  NK = c(0.0248, 0.0215),
  CD4T = c(0.164, 0.0723),
  CD8T = c(0.0779, 0.0380),
  Mono = c(0.0777, 0.0887),
  Neutro = c(0.622, 0.740)
)

# Function to adjust beta values for a specific sample
adjust_beta_values <- function(sample_name) {
  # Extract the beta values for the specified sample
  beta_values_sample <- betas_nc[[sample_name]]
  
  # Get the corresponding cell proportions for the specified sample
  cell_props_sample <- ed_nc %>%
    filter(Sample_ID == sample_name) %>%
    select(-Sample_ID,-LH_smoke_status)  # Remove Sample_ID for modeling
  
  # Create a data frame for modeling
  model_data <- data.frame(beta = beta_values_sample, cell_props_sample)
  
  # Fit a linear model for the current sample
  model <- lm(beta ~ ., data = model_data)
  
  # Calculate adjusted values
  adjusted_values <- predict(model)  # Use the model to predict adjusted values
  
  # Return adjusted values
  return(adjusted_values)
}

# Adjust beta values for sample A1
adjusted_A1 <- adjust_beta_values("A1")

# Adjust beta values for sample A2
adjusted_A2 <- adjust_beta_values("A2")

# Combine adjusted values into a data frame
adjusted_betas <- data.frame(
  CpG = beta_values$CpG,
  adjusted_A1 = adjusted_A1,
  adjusted_A2 = adjusted_A2
)

# View the adjusted beta values
print(adjusted_betas)

pdf(file = "output/cdmc/never_current_estimate_dmc.pdf", width=12, height=4)
for(i in 2:(length(colnames(ncbetas))-7)){
  B <- ggplot(ncbetas, aes(x=as.numeric(B), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "B")
  NK <- ggplot(ncbetas, aes(x=as.numeric(NK), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "NK")
  CD4T <- ggplot(ncbetas, aes(x=as.numeric(CD4T), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "CD4T")
  CD8T<- ggplot(ncbetas, aes(x=as.numeric(CD8T), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "CD8T")
  Mono <- ggplot(ncbetas, aes(x=as.numeric(Mono), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "Mono")
  Neutro<- ggplot(ncbetas, aes(x=as.numeric(Neutro), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "Neutro")
  p <- plot_grid(B,NK,CD4T,CD8T,Mono,Neutro, labels = "AUTO",scale = 1)
  print(p)
}
dev.off()





save(cdmc_nc,cdmc_nf,cdmc_fc, file='output/cdmc/cell_dmc_sensitvity.rdata')
write.table(cdmc_nc$dmct,sep='\t',file= "output/cdmc/cell_dmc_never_current.tsv",quote = F)
write.table(cdmc_fc$dmct,sep='\t',file= "output/cdmc/cell_dmc_former_current.tsv",quote = F)
write.table(cdmc_nf$dcmt,sep='\t',file= "output/cdmc/cell_dmc_never_former.tsv",quote = F)


