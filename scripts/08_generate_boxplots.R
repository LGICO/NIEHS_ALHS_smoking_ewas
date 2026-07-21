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
bonf<- 0.05/230325

load('/Users/goodneyga/OneDrive - National Institutes of Health/NHLBI/NIEHS_ahs/output/cleaned_MSA_niehs_alhs.rdata')
beta2 = t(clean_betas) # methylation data
betap <- t(clean_betas_pval_masked)
pheno = clean_sampsheet# phenotype data
pheno_pc <- pheno %>% filter(European == 'European',!is.na(PC1))
beta_pc <- beta2[row.names(beta2) %in% pheno_pc$Sample_ID,]
betap_pc <- betap[row.names(betap) %in% pheno_pc$Sample_ID,]

path = "/Users/goodneyga/OneDrive - National Institutes of Health/NHLBI/NIEHS_ahs/"

###### Model 1 #####################
siguws <- fread(file=paste0(path,'output/svyglm/unimputed_unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))
siguws %>% arrange(P)
probes<-siguws %>% arrange(P) %>% filter(P<=bonf) %>% dplyr::select(V1)
pdf(paste0(path,"model1_significant_smoking_probes_europc_only.pdf"))
for(p in as.list(probes$V1)){
  ahrr <- data.frame(beta_pc[,p])
  matching_rows <- grep(p, rownames(t(beta_pc)), value = TRUE)
  cpg_data <- t(beta_pc)[matching_rows, , drop = FALSE]
  cpg_long <- pivot_longer(as.data.frame(cpg_data),
                           cols = everything(),
                           names_to = "Sample_ID",
                           values_to = "Beta_Value")
  unique_samples <- unique(cpg_long$Sample_ID)
  ahrr<-left_join(cpg_long,pheno_pc,by="Sample_ID")
  
  #ahrr=ahrr%>%mutate(LH_smoke_status=case_when(LH_smoke_status==0~"Never",
  #                                             LH_smoke_status==1~"Former",
  #                                             LH_smoke_status==2~"Current",
  #                                             LH_smoke_status==99~"Missing"))
  ahrr$LH_smoke_status=as.factor(ahrr$LH_smoke_status)
  ahrr = filter(ahrr,LH_smoke_status!="M")
  gene <- ungroup(siguws) %>% dplyr::select(V1,genesUniq) %>% filter(V1 == p) %>%
    dplyr::select(genesUniq) #%>% separate(annot.symbol,into = 'x',extra='drop')
  #gene="AHRR;PDCD6-AHRR"
  print(ggboxplot(ahrr, x = "LH_smoke_status", y = "Beta_Value",
                  color = "LH_smoke_status", palette ="jama",order=c('Never','Former','Current'),
                  add = "jitter",xlab="Smoking Status", ylab = "Beta", title=paste0("Probe: ",p," Gene: ",gene,"\n by Smoking Status") )+
          theme(plot.title = element_text(hjust = 0.5)))
  
}
dev.off()

###### Model 3 #####################
siguws <- fread(file=paste0(path,'output/svyglm/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))
siguws %>% arrange(P)
probes<-siguws %>% arrange(P) %>% filter(P<=bonf) %>% dplyr::select(V1)
pdf(paste0(path,"model3_significant_smoking_probes_europc_only.pdf"))
for(p in as.list(probes$V1)){
  ahrr <- data.frame(beta_pc[,p])
  matching_rows <- grep(p, rownames(t(beta_pc)), value = TRUE)
  cpg_data <- t(beta_pc)[matching_rows, , drop = FALSE]
  cpg_long <- pivot_longer(as.data.frame(cpg_data),
                           cols = everything(),
                           names_to = "Sample_ID",
                           values_to = "Beta_Value")
  unique_samples <- unique(cpg_long$Sample_ID)
  ahrr<-left_join(cpg_long,pheno_pc,by="Sample_ID")
  
  #ahrr=ahrr%>%mutate(LH_smoke_status=case_when(LH_smoke_status==0~"Never",
  #                                             LH_smoke_status==1~"Former",
  #                                             LH_smoke_status==2~"Current",
  #                                             LH_smoke_status==99~"Missing"))
  ahrr$LH_smoke_status=as.factor(ahrr$LH_smoke_status)
  ahrr = filter(ahrr,LH_smoke_status!="M")
  gene <- ungroup(siguws) %>% dplyr::select(V1,genesUniq) %>% filter(V1 == p) %>%
    dplyr::select(genesUniq) #%>% separate(annot.symbol,into = 'x',extra='drop')
  #gene="AHRR;PDCD6-AHRR"
  print(ggboxplot(ahrr, x = "LH_smoke_status", y = "Beta_Value",
                  color = "LH_smoke_status", palette ="jama",order=c('Never','Former','Current'),
                  add = "jitter",xlab="Smoking Status", ylab = "Beta", title=paste0("Probe: ",p," Gene: ",gene,"\n by Smoking Status") )+
          theme(plot.title = element_text(hjust = 0.5)))
  
}
dev.off()

###### Model 10 #####################
siguws <- fread(file=paste0(path,'output/svyglm/unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))
siguws %>% arrange(P)
probes<-siguws %>% arrange(P) %>% filter(P<=bonf) %>% dplyr::select(V1)
pdf(paste0(path,"model10_significant_smoking_probes_europc_only.pdf"))
for(p in as.list(probes$V1)){
  ahrr <- data.frame(beta_pc[,p])
  matching_rows <- grep(p, rownames(t(beta_pc)), value = TRUE)
  cpg_data <- t(beta_pc)[matching_rows, , drop = FALSE]
  cpg_long <- pivot_longer(as.data.frame(cpg_data),
                           cols = everything(),
                           names_to = "Sample_ID",
                           values_to = "Beta_Value")
  unique_samples <- unique(cpg_long$Sample_ID)
  ahrr<-left_join(cpg_long,pheno_pc,by="Sample_ID")
  
  #ahrr=ahrr%>%mutate(LH_smoke_status=case_when(LH_smoke_status==0~"Never",
  #                                             LH_smoke_status==1~"Former",
  #                                             LH_smoke_status==2~"Current",
  #                                             LH_smoke_status==99~"Missing"))
  ahrr$LH_smoke_status=as.factor(ahrr$LH_smoke_status)
  ahrr = filter(ahrr,LH_smoke_status!="M")
  gene <- ungroup(siguws) %>% dplyr::select(V1,genesUniq) %>% filter(V1 == p) %>%
    dplyr::select(genesUniq) #%>% separate(annot.symbol,into = 'x',extra='drop')
  #gene="AHRR;PDCD6-AHRR"
  print(ggboxplot(ahrr, x = "LH_smoke_status", y = "Beta_Value",
                  color = "LH_smoke_status", palette ="jama",order=c('Never','Former','Current'),
                  add = "jitter",xlab="Smoking Status", ylab = "Beta", title=paste0("Probe: ",p," Gene: ",gene,"\n by Smoking Status") )+
          theme(plot.title = element_text(hjust = 0.5)))
  
}
dev.off()

###### Model 12 #####################
siguws <- fread(file=paste0(path,'output/svyglm/unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))
siguws %>% arrange(P)
probes<-siguws %>% arrange(P) %>% filter(P<=bonf) %>% dplyr::select(V1)
pdf(paste0(path,"model12_significant_smoking_probes_europc_only.pdf"))
for(p in as.list(probes$V1)){
  ahrr <- data.frame(beta_pc[,p])
  matching_rows <- grep(p, rownames(t(beta_pc)), value = TRUE)
  cpg_data <- t(beta_pc)[matching_rows, , drop = FALSE]
  cpg_long <- pivot_longer(as.data.frame(cpg_data),
                           cols = everything(),
                           names_to = "Sample_ID",
                           values_to = "Beta_Value")
  unique_samples <- unique(cpg_long$Sample_ID)
  ahrr<-left_join(cpg_long,pheno_pc,by="Sample_ID")
  
  #ahrr=ahrr%>%mutate(LH_smoke_status=case_when(LH_smoke_status==0~"Never",
  #                                             LH_smoke_status==1~"Former",
  #                                             LH_smoke_status==2~"Current",
  #                                             LH_smoke_status==99~"Missing"))
  ahrr$LH_smoke_status=as.factor(ahrr$LH_smoke_status)
  ahrr = filter(ahrr,LH_smoke_status!="M")
  gene <- ungroup(siguws) %>% dplyr::select(V1,genesUniq) %>% filter(V1 == p) %>%
    dplyr::select(genesUniq) #%>% separate(annot.symbol,into = 'x',extra='drop')
  #gene="AHRR;PDCD6-AHRR"
  print(ggboxplot(ahrr, x = "LH_smoke_status", y = "Beta_Value",
                  color = "LH_smoke_status", palette ="jama",order=c('Never','Former','Current'),
                  add = "jitter",xlab="Smoking Status", ylab = "Beta", title=paste0("Probe: ",p," Gene: ",gene,"\n by Smoking Status") )+
          theme(plot.title = element_text(hjust = 0.5)))
  
}
dev.off()

