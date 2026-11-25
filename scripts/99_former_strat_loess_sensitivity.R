
### methylation & phenotype data
library(dplyr)

load('751_clean_pc.rdata')
dir.create("output/model_output/svyglm", recursive=TRUE, showWarnings=FALSE)
dir.create("output/model_output/svyglm_anno", recursive=TRUE, showWarnings=FALSE)
dir.create("output/model_output/svyglm_significant_cpgs", recursive=TRUE, showWarnings=FALSE)
dir.create("output/cdmc", recursive=TRUE, showWarnings=FALSE)
dir.create("output/former_stratified", recursive=TRUE, showWarnings=FALSE)
dir.create("output/ecel1p1", recursive=TRUE, showWarnings=FALSE)
beta2 = t(clean_betas) # methylation data
betap <- t(clean_betas_pval_masked)
pheno = clean_sampsheet# phenotype data

### EWAS

library(survey)
library(MASS)
library(lmtest)
library(sandwich) 

library(parallel) 
#cl <- makeCluster(1, type="FORK")
cl <- makeCluster(getOption("cl.cores", 50))
clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(sandwich))
clusterEvalQ(cl,library(lmtest))
clusterEvalQ(cl,library(survey))
clusterEvalQ(cl,library(estimatr))
clusterEvalQ(cl,library(dplyr))
clusterEvalQ(cl,library(robustbase))


################################

pheno=pheno%>%mutate(LEL19=ifelse(LEL19=="Yes",1,ifelse(LEL19=="No",0,99)))
pheno$LEL19=as.factor(pheno$LEL19)

pheno=pheno%>%mutate(LH_smoke_status=case_when(LH_smoke_status=="Never"~0,
                                               LH_smoke_status=="Former"~1,
                                               LH_smoke_status=="Current"~2,
                                               LH_smoke_status=="M"~99))
pheno$LH_smoke_status=as.factor(pheno$LH_smoke_status)

pheno=pheno%>%mutate(European=case_when(European=="European"~1,
                                        European=="Non_European"~0,
                                        European=="GWAS sample failure"~99))
pheno$European=as.factor(pheno$European)

pheno$Gender <- as.factor(pheno$Gender)
pheno$inrdsx <-as.factor(pheno$inrdsx)
pheno <- pheno %>% mutate( ANYETS = case_when(ANYETS =="No"~0,ANYETS=="Yes"~1))
pheno <- pheno %>% mutate(LHSSCREEN_12142015 = case_when(LHSSCREEN_12142015 == "Control"~0,LHSSCREEN_12142015 == "Case"~1))
pheno <- pheno %>% mutate(Gender = case_when(Gender == 'male'~0,Gender=='female'~1))
pheno <- pheno %>% mutate(STATE = case_when(STATE == 'IA'~0,STATE=='NC'~1))
pheno$LHSSCREEN_12142015 = as.factor(pheno$LHSSCREEN_12142015)
pheno$ANYETS = as.factor(pheno$ANYETS)
pheno$STATE = as.factor(pheno$STATE)
pheno <- pheno %>% mutate(smokeXasthma = paste0(LH_smoke_status,'_',LHSSCREEN_12142015))
pheno$smokeXasthma <- as.factor(pheno$smokeXasthma)
dim(pheno)
pheno <- pheno %>% filter(Sample_ID != "WG0342743_DNA_D07", Sample_ID != "WG0342734_DNA_H09")
beta2 = beta2[row.names(beta2) %in% pheno$Sample_ID,]
betap = betap[row.names(betap) %in% pheno$Sample_ID,]
dim(pheno)
pheno_pc <- pheno %>% filter(European == 1,!is.na(PC1))
beta_pc <- beta2[row.names(beta2) %in% pheno_pc$Sample_ID,]
betap_pc <- betap[row.names(betap) %in% pheno_pc$Sample_ID,]


pheno_pc_case <- pheno_pc %>% filter(LHSSCREEN_12142015 == 1)
beta_pc_case <- beta_pc[row.names(beta_pc) %in% pheno_pc_case$Sample_ID,]
betap_pc_case <- betap_pc[row.names(betap_pc) %in% pheno_pc_case$Sample_ID,]
pheno_pc_control <- pheno_pc %>% filter(LHSSCREEN_12142015 == 0)
beta_pc_control <- beta_pc[row.names(beta_pc) %in% pheno_pc_control$Sample_ID,]
betap_pc_control <- betap_pc[row.names(betap_pc) %in% pheno_pc_control$Sample_ID,]


##########################################################################################################################################
## former-smokers
##########
##########################################################################################################################################

##Model 1 unweighted euro
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf = summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 1)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 0)
  coef = cf$coefficients[2,"Estimate"]
  se = cf$coefficients[2,"Std. Error"]
  p.value = cf$coefficients[2,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

##Model 10
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + European + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf = coeftest(mdl, vcov=vcov(mdl))
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 1)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 0)
  coef = cf[3,"Estimate"]
  se = cf[3,"Std. Error"]
  p.value = cf[3,"Pr(>|z|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta2, 2, f.RLM.Adj.Robust.unwgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.unwgt, pheno)
results <- as.data.frame(t(results))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")




##weighted
################
##formersmokers
p_asthma_overall = 0.06852481 # n_asthma_overall/(n_asthma_overall+n_nonasthma_overall) 

n_asthma = sum(pheno_pc$LHSSCREEN_12142015 == 1)
n_nonasthma = sum(pheno_pc$LHSSCREEN_12142015 == 0)
p_asthma_nested = n_asthma/(n_asthma+n_nonasthma) # n_asthma_nested/(n_asthma_nested+n_nonasthma_nested)   
p_nonasthma_overall = 1- p_asthma_overall # n_nonasthma_overall/(n_asthma_overall+n_nonasthma_overall) 
p_nonasthma_nested = 1- p_asthma_nested # n_nonasthma_nested/(n_asthma_nested+n_nonasthma_nested) 
w_asthma <- 1/(p_asthma_nested/p_asthma_overall) 
w_nonasthma <- 1/(p_nonasthma_nested/p_nonasthma_overall) 
pheno_pc$asthma_weights [pheno_pc$LHSSCREEN_12142015 == 0] <- w_nonasthma 
pheno_pc$asthma_weights [pheno_pc$LHSSCREEN_12142015 == 1] <- w_asthma 


###########
## Model 3
f.RLM.Adj.Robust.wgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, weights = ~asthma_weights, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS  + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf = coeftest(mdl, vcov=vcov(mdl))
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 1)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 0)
  coef = cf[2,"Estimate"]
  se = cf[2,"Std. Error"]
  p.value = cf[2,"Pr(>|z|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc)
results <- as.data.frame(t(results))
write.csv(results, "output/model_output/svyglm/unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

##############################
##Model 12
n_asthma = sum(pheno$LHSSCREEN_12142015 == 1)
n_nonasthma = sum(pheno$LHSSCREEN_12142015 == 0)
p_asthma_nested = n_asthma/(n_asthma+n_nonasthma) # n_asthma_nested/(n_asthma_nested+n_nonasthma_nested)   
p_nonasthma_overall = 1- p_asthma_overall # n_nonasthma_overall/(n_asthma_overall+n_nonasthma_overall) 
p_nonasthma_nested = 1- p_asthma_nested # n_nonasthma_nested/(n_asthma_nested+n_nonasthma_nested) 
w_asthma <- 1/(p_asthma_nested/p_asthma_overall) 
w_nonasthma <- 1/(p_nonasthma_nested/p_nonasthma_overall) 
pheno$asthma_weights [pheno$LHSSCREEN_12142015 == 0] <- w_nonasthma 
pheno$asthma_weights [pheno$LHSSCREEN_12142015 == 1] <- w_asthma 
f.RLM.Adj.Robust.wgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, weights = ~asthma_weights, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + European + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf = coeftest(mdl, vcov=vcov(mdl))
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 1)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 0)
  coef = cf[2,"Estimate"]
  se = cf[2,"Std. Error"]
  p.value = cf[2,"Pr(>|z|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta2, 2, f.RLM.Adj.Robust.wgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.wgt, pheno)
results <- as.data.frame(t(results))

write.csv(results, "output/model_output/svyglm/unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")



#########################################################
########################################################
## CellDMC
library(EpiDISH)
## Split the data into never former and never current
pheno <- pheno %>% mutate( ANYETS = case_when(ANYETS ==0~"No",ANYETS==1~"Yes",.default = "M"))
pheno$ANYETS <- factor(pheno$ANYETS,levels=c("No","Yes","M"))
nf = pheno %>% filter(LH_smoke_status == 0 | LH_smoke_status == 1)
nc = pheno %>% filter(LH_smoke_status == 0 | LH_smoke_status == 2)
fc = pheno %>% filter(LH_smoke_status == 1 | LH_smoke_status == 2)

betas_nf = clean_betas[,as.vector(nf$Sample_ID)]
betas_nc = clean_betas[,as.vector(nc$Sample_ID)]
betas_fc = clean_betas[,as.vector(fc$Sample_ID)]

ed_nf = nf %>% dplyr::select(B, NK, CD4T, CD8T, Mono, Neutro)
ed_nc = nc %>% dplyr::select(B, NK, CD4T, CD8T, Mono, Neutro)
ed_fc = fc %>% dplyr::select(B, NK, CD4T, CD8T, Mono, Neutro)

cdmc_nf = CellDMC(beta.m=betas_nf, pheno.v=nf$LH_smoke_status,frac.m=as.matrix(ed_nf),
                  cov.mod = model.matrix(~ Gender + AgeAtHomeVisit + STATE + ANYETS +
                                           LH_BMI + European + LHSSCREEN_12142015 +
                                           PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
                                           PC8 + PC9 + PC10,data = nf),mc.cores=50)

cdmc_nc = CellDMC(beta.m=betas_nc, pheno.v=nc$LH_smoke_status,frac.m=as.matrix(ed_nc),
                  cov.mod = model.matrix(~ Gender + AgeAtHomeVisit + STATE + ANYETS + 
                                           LH_BMI + European + LHSSCREEN_12142015 +
                                           PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
                                           PC8 + PC9 + PC10,data = nc),mc.cores=50)

cdmc_fc = CellDMC(beta.m=betas_fc, pheno.v=fc$LH_smoke_status,frac.m=as.matrix(ed_fc),
                  cov.mod = model.matrix(~ Gender + AgeAtHomeVisit + STATE + ANYETS + 
                                           LH_BMI + European + LHSSCREEN_12142015 +
                                           PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
                                           PC8 + PC9 + PC10,data = fc),mc.cores=50)

nfcpg <- cdmc_nf$dmct %>% as_tibble(.,rownames="CpG") %>% filter(DMC == 1) %>% dplyr::select(CpG)
nccpg <- cdmc_nc$dmct %>% as_tibble(.,rownames="CpG") %>% filter(DMC == 1) %>% dplyr::select(CpG)
fccpg <- cdmc_fc$dmct %>% as_tibble(.,rownames="CpG") %>% filter(DMC == 1) %>% dplyr::select(CpG)
library(janitor)
library(tibble)
ed_nc = nc %>% dplyr::select(Sample_ID,B, NK, CD4T, CD8T, Mono, Neutro,LH_smoke_status)
ed_fc = fc %>% dplyr::select(Sample_ID,B, NK, CD4T, CD8T, Mono, Neutro,LH_smoke_status)
ncbetas <- as_tibble(subset(betas_nc,rownames(clean_betas) %in% nccpg$CpG),rownames="CpG") %>% t %>%
  as.data.frame %>% row_to_names(row_number = 1) 
ncbetas <- left_join(rownames_to_column(ncbetas),ed_nc,by =c("rowname"= 'Sample_ID'))
fcbetas <- as_tibble(subset(betas_fc,rownames(clean_betas) %in% fccpg$CpG),rownames="CpG")%>% t %>%
  as.data.frame %>% row_to_names(row_number = 1) 
fcbetas <- left_join(rownames_to_column(fcbetas),ed_fc,by =c("rowname"= 'Sample_ID'))

library(cowplot)
library(egg)
pdf(file = "output/cdmc/never_current_beta_dmc.pdf", width=12, height=4)
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

pdf(file = "output/cdmc/former_current_beta_dmc.pdf", width=12, height=4)
for(i in 2:(length(colnames(fcbetas))-7)){
  B <- ggplot(fcbetas, aes(x=as.numeric(B), y=as.numeric(!!sym(colnames(fcbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(fcbetas)[i],y = "B")
  NK <- ggplot(fcbetas, aes(x=as.numeric(NK), y=as.numeric(!!sym(colnames(fcbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(fcbetas)[i],y = "NK")
  CD4T <- ggplot(fcbetas, aes(x=as.numeric(CD4T), y=as.numeric(!!sym(colnames(fcbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(fcbetas)[i],y = "CD4T")
  CD8T<- ggplot(fcbetas, aes(x=as.numeric(CD8T), y=as.numeric(!!sym(colnames(fcbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(fcbetas)[i],y = "CD8T")
  Mono <- ggplot(fcbetas, aes(x=as.numeric(Mono), y=as.numeric(!!sym(colnames(fcbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(fcbetas)[i],y = "Mono")
  Neutro<- ggplot(fcbetas, aes(x=as.numeric(Neutro), y=as.numeric(!!sym(colnames(fcbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(fcbetas)[i],y = "Neutro")
  p <- plot_grid(B,NK,CD4T,CD8T,Mono,Neutro, labels = "AUTO",scale=1)
  print(p)
}
dev.off()

save(cdmc_nc,cdmc_nf,cdmc_fc, file='output/cdmc/cell_dmc_sensitvity.rdata')
write.table(cdmc_nc$dmct,sep='\t',file= "output/cdmc/cell_dmc_never_current.tsv",quote = F)
write.table(cdmc_fc$dmct,sep='\t',file= "output/cdmc/cell_dmc_former_current.tsv",quote = F)
write.table(cdmc_nf$dcmt,sep='\t',file= "output/cdmc/cell_dmc_never_former.tsv",quote = F)


stopCluster(cl)


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
annots<- c("hg38_genes_1to5kb","hg38_genes_promoters", 
           "hg38_genes_cds", 
           "hg38_genes_5UTRs", 
           "hg38_genes_exons",
           "hg38_basicgenes", 
           "hg38_cpgs",
           "hg38_genes_firstexons","hg38_genes_introns", 
           "hg38_genes_intronexonboundaries", 
           "hg38_genes_exonintronboundaries","hg38_genes_3UTRs", 
           "hg38_genes_intergenic", "hg38_cpg_islands", "hg38_cpg_shores", 
           "hg38_cpg_shelves")#, 
#"hg38_cpg_inter", "hg38_lncrna_gencode")
annotations = build_annotations(genome = 'hg38', annotations = annots)



################
## SVYGLM
###############
fin<- system('ls output/model_output/svyglm/*.csv',intern=T)
bonf<- 0.05/230325 
tsv_path <- "references/MSA.hg38.manifest.gencode.v41.tsv.gz"
manifest = fread(tsv_path)

for(f in fin){
  uws <- fread(f)
  siguws <- uws %>% arrange(P)#%>% dplyr::filter(uws$P<=bonf)
  siguw_anno <- left_join(siguws,manifest,by=c("V1"="probeID"))
  f <- substring(f,74) %>% substring(.,1,nchar(.)-4)
  write.table(siguw_anno,file=paste0("output/model_output/svyglm_anno/",f,".tsv"),sep='\t',quote = F,na = "",row.names = F)
}     
fin<- system('ls output/model_output/svyglm_anno/*.tsv',intern=T)

for(f in fin){
  uws <- fread(f,header = T)
  print(f)
  colnames(uws)<- c('Probe_ID','Beta','SE','P','N','N.exp','N.unexp','CpG_chrm','CpG_beg','CpG_end','probe_strand','genesUniq','geneNames','transcriptTypes','transcriptIDs','distToTSS')
  siguws <- uws %>% dplyr::filter(uws$P<=bonf) %>% arrange(P) 
  f <- substring(f,79)
  write.table(siguws,file=paste0("output/model_output/svyglm_significant_cpgs/",f),sep='\t',quote = F,na = "",row.names = F)
}    



#####################################
####################################
##plot loess
fm1 <- fread('output/model_output/svyglm_significant_cpgs/unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv') %>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1|AHRR|PDCD6-AHRR')) %>%
  dplyr::select(Probe_ID,genesUniq)
fm3 <- fread('output/model_output/svyglm_significant_cpgs/unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1|AHRR|PDCD6-AHRR')) %>%
  dplyr::select(Probe_ID,genesUniq)
fm10 <- fread('output/model_output/svyglm_significant_cpgs/weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1|AHRR|PDCD6-AHRR')) %>%
  dplyr::select(Probe_ID,genesUniq)
fm12 <- fread('output/model_output/svyglm_significant_cpgs/weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1|AHRR|PDCD6-AHRR')) %>%
  dplyr::select(Probe_ID,genesUniq)



fm1 <- fread('output/model_output/svyglm_significant_cpgs/unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv') %>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1')) %>%
  dplyr::select(Probe_ID,genesUniq)
fm3 <- fread('output/model_output/svyglm_significant_cpgs/unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1')) %>%
  dplyr::select(Probe_ID,genesUniq)
fm10 <- fread('output/model_output/svyglm_significant_cpgs/weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1')) %>%
  dplyr::select(Probe_ID,genesUniq)
fm12 <- fread('output/model_output/svyglm_significant_cpgs/weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1')) %>%
  dplyr::select(Probe_ID,genesUniq)


fm <- bind_rows(list(fm1,fm3,fm10,fm12)) %>% distinct()
probes <- fm$Probe_ID
probes <- c(probes,'cg17087741_TC21','cg21566642_TC11','cg21566642_TC12','cg21566642_TC13')

new_pheno <- data.table::fread('samplesheets/lhs28901_Methylation_Phenotype_2025APR29.txt',sep=',')
np<-new_pheno%>% filter(LHID %in% pheno$LHID,COMMENT!='Replicate')
x<-colnames(np)[!colnames(np) %in% colnames(pheno)]
np <- np %>% dplyr::select(LHID,partid, x)
pheno <- left_join(pheno,np,by=c('LHID','partid'))
pheno=pheno%>%mutate(LEL19=ifelse(LEL19=="Yes",1,ifelse(LEL19=="No",0,99)))
pheno$LEL19=as.factor(pheno$LEL19)

pheno=pheno%>%mutate(LH_smoke_status=case_when(LH_smoke_status=="Never"~0,
                                               LH_smoke_status=="Former"~1,
                                               LH_smoke_status=="Current"~2,
                                               LH_smoke_status=="M"~99))
pheno$LH_smoke_status=as.factor(pheno$LH_smoke_status)

pheno=pheno%>%mutate(European=case_when(European=="European"~1,
                                        European=="Non_European"~0,
                                        European=="GWAS sample failure"~99))
pheno$European=as.factor(pheno$European)
#inrdsx<- tibble::rownames_to_column(data.frame(inrdsx),var="Sample_ID")
#pheno <- left_join(pheno,inrdsx,by="Sample_ID")
#pheno <- pheno %>% mutate(inrdsx = tolower(inrdsx))
pheno$Gender <- as.factor(pheno$Gender)
pheno$inrdsx <-as.factor(pheno$inrdsx)
pheno <- pheno %>% mutate( ANYETS = case_when(ANYETS =="No"~0,ANYETS=="Yes"~1,.default=99))
pheno <- pheno %>% mutate(LHSSCREEN_12142015 = case_when(LHSSCREEN_12142015 == "Control"~0,LHSSCREEN_12142015 == "Case"~1))
pheno <- pheno %>% mutate(Gender = case_when(Gender == 'male'~0,Gender=='female'~1))
pheno <- pheno %>% mutate(STATE = case_when(STATE == 'IA'~0,STATE=='NC'~1))
pheno$LHSSCREEN_12142015 = as.factor(pheno$LHSSCREEN_12142015)
pheno$ANYETS = as.factor(pheno$ANYETS)
pheno$STATE = as.factor(pheno$STATE)
pheno <- pheno %>% mutate(smokeXasthma = paste0(LH_smoke_status,'_',LHSSCREEN_12142015))
pheno$smokeXasthma <- as.factor(pheno$smokeXasthma)
pheno$LH_LastSmkedLE5YrAgo_v1 <-factor(pheno$LH_LastSmkedLE5YrAgo_v1, levels = c("0) LE 5 yrs ago","1) GT 5 yrs ago"))

pheno = pheno %>% mutate(smoke_statusF = case_when(LH_smoke_status == 1 & LH_LastSmkedLE5YrAgo_v1=="0) LE 5 yrs ago" ~ factor(3),
                                                   LH_smoke_status == 1 & LH_LastSmkedLE5YrAgo_v1=="1) GT 5 yrs ago" ~ factor(4),
                                                   .default = LH_smoke_status))
pheno$smoke_statusF <- as.factor(pheno$smoke_statusF) %>% droplevels(.,exclude=1)
pheno$smoke_statusF <- factor(pheno$smoke_statusF, levels = c(0,2,3,4,99))

pheno <- pheno %>% mutate(former_5yrs_split = case_when(LH_smoke_status == 0 ~ factor('Never Smoker'),
                                                        LH_smoke_status == 2 ~ factor('Current Smoker'),
                                                        LH_smoke_status == 1 & LH_TimeSinceSmked_v2 <= 5 ~ factor('<=5yrs'),
                                                        LH_smoke_status == 1 & LH_TimeSinceSmked_v2 > 5 & LH_TimeSinceSmked_v2 <= 10 ~ factor('6-10yrs'),
                                                        LH_smoke_status == 1 & LH_TimeSinceSmked_v2 > 10 & LH_TimeSinceSmked_v2 <= 15 ~ factor('11-15yrs'),
                                                        LH_smoke_status == 1 & LH_TimeSinceSmked_v2 > 15 ~ factor('>15yrs'),
                                                        LH_smoke_status == 99 ~ factor('Missing')))
pheno <- pheno %>% mutate(former_10yrs_split = case_when(LH_smoke_status == 0 ~ factor('Never Smoker'),
                                                         LH_smoke_status == 2 ~ factor('Current Smoker'),
                                                         LH_smoke_status == 1 & LH_TimeSinceSmked_v2 <= 10 ~ factor('<=10yrs'),
                                                         LH_smoke_status == 1 & LH_TimeSinceSmked_v2 > 10 & LH_TimeSinceSmked_v2 <= 20 ~ factor('11-20yrs'),
                                                         LH_smoke_status == 1 & LH_TimeSinceSmked_v2 > 20 & LH_TimeSinceSmked_v2 <= 30 ~ factor('21-30yrs'),
                                                         LH_smoke_status == 1 & LH_TimeSinceSmked_v2 > 30 ~ factor('>30yrs'),
                                                         LH_smoke_status == 99 ~ factor('Missing')))
epic <- fread('references/EPIC.hg38.manifest.gencode.v36.tsv.gz') %>% mutate(epicT =T) %>% dplyr::select(probeID,epicT)
fourfity<- fread('references/HM450.hg38.manifest.tsv.gz')%>% mutate(ffiftyT =T) %>% dplyr::select(Probe_ID,ffiftyT)
##reduce probes to cpgs
probies <- manifest %>% filter(probeID %in% probes) %>% dplyr::select(probeID,genesUniq) %>%
  mutate(cpg = gsub(x = probeID, '\\_.*','')) %>% 
  left_join(epic,by=c('cpg'='probeID')) %>% left_join(fourfity,by=c('cpg'='Probe_ID')) %>% replace_na(list(epicT=FALSE,ffiftyT=FALSE)) %>%
  mutate(msa_only=case_when(epicT==TRUE | ffiftyT==TRUE ~ FALSE,.default=TRUE))

probes<-c(
  "cg04284215_BC11", "cg21566642_TC13", "cg21566642_TC11",
  "cg04284185_TC21", "cg21566642_TC12", "cg04284186_TC11",
  "cg04284209_BC11", "cg04284216_BC11", "cg04284207_BC11",
  "cg04284205_TC11", "cg06644428_BC13", "cg06644428_BC11",
  "cg06644428_BC12", "cg06644428_BC14", "cg04284171_TC11",
  "cg03329539_BC21", "cg04284170_TC11", "cg17087741_TC21",
  "cg21566642_TC11", "cg21566642_TC12", "cg21566642_TC13"
)
####################
## With avg lines
pdf("output/former_stratified/ahrr_ecel1p1_mean_scatter_plot.pdf",width=11)
for(p in probes){
  ##gsub to get correct names
  ##we need to grab current and never betas to get mean and %
  ahrr <- data.frame(beta2[,p])
  matching_rows <- grep(p, rownames(t(beta2)), value = TRUE)
  cpg_data <- t(beta2)[matching_rows, , drop = FALSE]
  cpg_long <- pivot_longer(as.data.frame(cpg_data),
                           cols = everything(),
                           names_to = "Sample_ID",
                           values_to = "Beta_Value")
  unique_samples <- unique(cpg_long$Sample_ID)
  ahrr<-left_join(cpg_long,pheno,by="Sample_ID")
  
  ahrr = filter(ahrr,former_5yrs_split!="Missing")
  ahrr = ahrr %>% dplyr::mutate(LH_TimeSinceSmked_v1.2 = case_when(LH_smoke_status == 0 ~ 67,
                                                                   LH_smoke_status == 2 ~ 0,
                                                                   .default = LH_TimeSinceSmked_v1),
                                smoke_status = case_when(LH_smoke_status == 0 ~ 'Never Smoker',
                                                         LH_smoke_status == 2 ~ 'Current Smoker',
                                                         LH_smoke_status == 1 ~ 'Former Smoker')) 
  
  curr = ahrr %>% filter(smoke_status == 'Current Smoker')
  currb = mean(curr$Beta_Value)
  curr95 <- quantile(curr$Beta_Value, 0.75) 
  curr5 <- quantile(curr$Beta_Value,0.25)
  nev = ahrr %>% filter(smoke_status == 'Never Smoker')
  nevb = mean(nev$Beta_Value)
  nev95 = quantile(nev$Beta_Value, 0.75) 
  nev5 <- quantile(nev$Beta_Value,0.25)
  ahrrf = ahrr %>% filter(smoke_status=='Former Smoker')
  gene <- manifest %>% dplyr::filter(probeID == p) %>% dplyr::select(genesUniq)
  print(ggplot(ahrrf, aes(x = LH_TimeSinceSmked_v1.2, y = Beta_Value)) + geom_point(aes(color=smoke_status))  + scale_color_manual(values = c('Former Smoker'='#DF8F44FF')) + labs(color='Smoke Status') + 
          geom_smooth(method='loess', formula= y~x,color = "#374E55FF") + geom_hline(aes(yintercept=currb,linetype="Current Smokers"),color = "#B24745FF") +  ylim(0,1) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = curr5, ymax = curr95, fill = "#B24745FF", alpha = .1, color = NA) +
          geom_hline(aes(yintercept=nevb,linetype="Never Smokers"),color = "#00A1D5FF") + annotate("rect", xmin = -Inf, xmax = Inf, ymin = nev5, ymax = nev95, fill = "#00A1D5FF", alpha = .1, color = NA) +
          scale_linetype_manual(name = "Mean Beta", 
                                values = c(2, 2)) +
          theme(plot.title = element_text(hjust = 0.5)) + theme_classic() + ggtitle(paste0('Probe: ',p,' Gene: ',gene)))
}
dev.off()

####################
## With loess degree 0 defaults to moving average
pdf("output/former_stratified/ahrr_ecel1p1_moving_avg_scatter_plot.pdf",width=11)
for(p in probes){
  ##gsub to get correct names
  ##we need to grab current and never betas to get mean and %
  ahrr <- data.frame(beta2[,p])
  matching_rows <- grep(p, rownames(t(beta2)), value = TRUE)
  cpg_data <- t(beta2)[matching_rows, , drop = FALSE]
  cpg_long <- pivot_longer(as.data.frame(cpg_data),
                           cols = everything(),
                           names_to = "Sample_ID",
                           values_to = "Beta_Value")
  unique_samples <- unique(cpg_long$Sample_ID)
  ahrr<-left_join(cpg_long,pheno,by="Sample_ID")
  
  ahrr = filter(ahrr,former_5yrs_split!="Missing")
  ahrr = ahrr %>% dplyr::mutate(LH_TimeSinceSmked_v1.2 = case_when(LH_smoke_status == 0 ~ 67,
                                                                   LH_smoke_status == 2 ~ 0,
                                                                   .default = LH_TimeSinceSmked_v1),
                                smoke_status = case_when(LH_smoke_status == 0 ~ 'Never Smoker',
                                                         LH_smoke_status == 2 ~ 'Current Smoker',
                                                         LH_smoke_status == 1 ~ 'Former Smoker')) 
  
  curr = ahrr %>% filter(smoke_status == 'Current Smoker')
  currb = mean(curr$Beta_Value)
  curr95 <- quantile(curr$Beta_Value, 0.75) 
  curr5 <- quantile(curr$Beta_Value,0.25)
  nev = ahrr %>% filter(smoke_status == 'Never Smoker')
  nevb = mean(nev$Beta_Value)
  nev95 = quantile(nev$Beta_Value, 0.75) 
  nev5 <- quantile(nev$Beta_Value,0.25)
  ahrrf = ahrr %>% filter(smoke_status=='Former Smoker')
  gene <- manifest %>% dplyr::filter(probeID == p) %>% dplyr::select(genesUniq)
  print(ggplot(ahrrf, aes(x = LH_TimeSinceSmked_v1.2, y = Beta_Value)) + geom_point(aes(color=smoke_status))  + scale_color_manual(values = c('Former Smoker'='#DF8F44FF')) + labs(color='Smoke Status') + 
          geom_smooth(method=stats::loess, formula= y~x, method.args=list(degree=0),color = "#374E55FF") + geom_hline(aes(yintercept=currb,linetype="Current Smokers"),color = "#B24745FF") +  ylim(0,1) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = curr5, ymax = curr95, fill = "#B24745FF", alpha = .1, color = NA) +
          geom_hline(aes(yintercept=nevb,linetype="Never Smokers"),color = "#00A1D5FF") + annotate("rect", xmin = -Inf, xmax = Inf, ymin = nev5, ymax = nev95, fill = "#00A1D5FF", alpha = .1, color = NA) +
          scale_linetype_manual(name = "Mean Beta", 
                                values = c(2, 2)) +
          theme(plot.title = element_text(hjust = 0.5)) + theme_classic() + ggtitle(paste0('Probe: ',p,' Gene: ',gene)))
}
dev.off()

pdf("output/former_stratified/ahrr_ecel1p1_median_scatter_plot.pdf",width=11)
for(p in probes){
  ##gsub to get correct names
  ##we need to grab current and never betas to get mean and %
  ahrr <- data.frame(beta2[,p])
  matching_rows <- grep(p, rownames(t(beta2)), value = TRUE)
  cpg_data <- t(beta2)[matching_rows, , drop = FALSE]
  cpg_long <- pivot_longer(as.data.frame(cpg_data),
                           cols = everything(),
                           names_to = "Sample_ID",
                           values_to = "Beta_Value")
  unique_samples <- unique(cpg_long$Sample_ID)
  ahrr<-left_join(cpg_long,pheno,by="Sample_ID")
  
  ahrr = filter(ahrr,former_5yrs_split!="Missing")
  ahrr = ahrr %>% dplyr::mutate(LH_TimeSinceSmked_v1.2 = case_when(LH_smoke_status == 0 ~ 67,
                                                                   LH_smoke_status == 2 ~ 0,
                                                                   .default = LH_TimeSinceSmked_v1),
                                smoke_status = case_when(LH_smoke_status == 0 ~ 'Never Smoker',
                                                         LH_smoke_status == 2 ~ 'Current Smoker',
                                                         LH_smoke_status == 1 ~ 'Former Smoker')) 
  
  curr = ahrr %>% filter(smoke_status == 'Current Smoker')
  currb = median(curr$Beta_Value)
  curr95 <- quantile(curr$Beta_Value, 0.75) 
  curr5 <- quantile(curr$Beta_Value,0.25)
  nev = ahrr %>% filter(smoke_status == 'Never Smoker')
  nevb = median(nev$Beta_Value)
  nev95 = quantile(nev$Beta_Value, 0.75) 
  nev5 <- quantile(nev$Beta_Value,0.25)
  ahrrf = ahrr %>% filter(smoke_status=='Former Smoker')
  gene <- manifest %>% dplyr::filter(probeID == p) %>% dplyr::select(genesUniq)
  print(ggplot(ahrrf, aes(x = LH_TimeSinceSmked_v1.2, y = Beta_Value)) + geom_point(aes(color=smoke_status))  + scale_color_manual(values = c('Former Smoker'='#DF8F44FF')) + labs(color='Smoke Status') + 
          geom_quantile(method='rqss',quantiles=0.5,lambda=0.75, color = "#374E55FF") + geom_hline(aes(yintercept=currb,linetype="Current Smokers"),color = "#B24745FF") +  
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = curr5, ymax = curr95, fill = "#B24745FF", alpha = .1, color = NA) +
          geom_hline(aes(yintercept=nevb,linetype="Never Smokers"),color = "#00A1D5FF") + annotate("rect", xmin = -Inf, xmax = Inf, ymin = nev5, ymax = nev95, fill = "#00A1D5FF", alpha = .1, color = NA) +
          scale_linetype_manual(name = "Median Beta", 
                                values = c(2, 2)) +
          theme(plot.title = element_text(hjust = 0.5)) + theme_classic() + ggtitle(paste0('Probe: ',p,' Gene: ',gene)))
}
dev.off()


########################
## betas instead of t-statistic
#####################################################################################################################################
genes <- manifest %>% filter(grepl('ECEL1P1|AHRR',genesUniq)) %>% dplyr::select(probeID,genesUniq) %>%
  mutate(cpg = gsub(x = probeID, '\\_.*','')) %>% 
  left_join(epic,by=c('cpg'='probeID')) %>% left_join(fourfity,by=c('cpg'='Probe_ID')) %>% replace_na(list(epicT=FALSE,ffiftyT=FALSE)) %>%
  mutate(msa_only=case_when(epicT==TRUE | ffiftyT==TRUE ~ FALSE,.default=TRUE))
ahec = beta2 %>% data.frame(sample_id = row.names(beta2),beta2) %>% dplyr::select(sample_id,genes$probeID) 
ahec <- ahec %>% pivot_longer(!sample_id,names_to = 'probeID', values_to = 'Beta' ) 
ahec <- ahec %>% left_join(pheno,by=c("sample_id" = "Sample_ID"))
ahec <- ahec %>% left_join(genes,by = 'probeID')
ahec5 <- ahec %>% dplyr::select(sample_id,Beta,genesUniq,former_5yrs_split) %>% dplyr::group_by(genesUniq,former_5yrs_split) %>% summarise(mean_beta = mean(Beta,na.rm=T))%>%ungroup()
ahec10 <- ahec %>% dplyr::select(sample_id,Beta,genesUniq,former_10yrs_split) %>% dplyr::group_by(genesUniq,former_10yrs_split) %>% summarise(mean_beta = mean(Beta,na.rm=T))%>%ungroup()
pdf("output/former_stratified/ahrr_ecel1p1_former_5yr_mean_beta.pdf",width=10)
ggplot(data=ahec5, aes(x=former_5yrs_split, y=mean_beta,group = genesUniq)) + geom_line(aes(color=genesUniq)) +
  ggtitle(paste0('Model ',i,'long-term')) + xlab('Years after quitting') + ylab('T-statistics') + theme_classic()

dev.off()



###########
## check DMR
ecels <- manifest %>% filter(grepl('ECEL1P1',genesUniq)) %>% dplyr::select(probeID) 
ec = beta2 %>% data.frame(sample_id = row.names(beta2),beta2) %>% dplyr::select(sample_id,any_of(ecels$probeID))
lec <- ec %>% pivot_longer(!sample_id,names_to = 'probeID', values_to = 'Beta' ) %>% mutate(cpg = gsub(x = probeID, '\\_.*',''))
lec <- lec %>% left_join(pheno,by=c("sample_id" = "Sample_ID")) %>% dplyr::filter(former_10yrs_split == '>30yrs') %>%
  dplyr::select(sample_id,probeID,Beta,cpg) %>% 
  dplyr::mutate(dmr = case_when(cpg == 'cg04284171' | cpg == 'cg04284185' | cpg == 'cg04284186' |
                                  cpg == 'cg04284205' | cpg == 'cg04284207' | cpg == 'cg04284209' |
                                  cpg == 'cg04284215' | cpg == 'cg04284216' ~ 1,.default=0))


t.test(Beta~dmr,data=lec)
wilcox.test(Beta~dmr,data=lec)

library(lme4)
library(lmerTest)
summary(lmer(Beta~dmr + (1 | sample_id), data=lec))


###############################################
##DMR
library(karyoploteR)
library(AnnotationHub)
library(data.table)
#library(sesame)
library(dplyr)
library(gtools)
library(ggplot2)
library(ggrepel)
library(stringr)
library(gdata)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(ChIPseeker)
library(plyranges)
#ahub <- AnnotationHub()
#ahub["AH5086"]
bonf<- 0.05/230325
bonf_epic = 0.05/817235
#sesameDataCache()
manifest = fread(tsv_path)
load('751_clean_pc.rdata')
mani = fread("references/MSA.hg38.manifest.gencode.v41.tsv.gz") %>% dplyr::select(CpG_chrm, CpG_beg, CpG_end, probe_strand,probeID) %>%
  dplyr::filter(probeID %in% rownames(clean_betas))
smokers=fread('output/model_output/svyglm_anno/unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)
names(smokers)<- c("probeID","chr","pos","smokers_pvalue","gene")
smokers <- smokers %>% dplyr::select(probeID,smokers_pvalue)

formers=fread('output/model_output/svyglm_anno/unimputed_unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)
names(formers)<- c("probeID","chr","pos","formers_pvalue","gene")
formers <- formers %>% dplyr::select(probeID,formers_pvalue)
pvals <- dplyr::left_join(smokers,formers)
clean_betas <- data.frame(probeID = row.names(clean_betas), clean_betas)
cb <- dplyr::left_join(pvals,clean_betas,by = 'probeID')



mani_cpg <- mani %>% mutate(probeID=gsub("\\_.*","",probeID)) %>% distinct()
## We need to make a granges object
mani_gr <- makeGRangesFromDataFrame(mani,keep.extra.columns = T,ignore.strand = F,
                                    seqnames.field = "CpG_chrm",start.field = "CpG_beg",
                                    end.field = "CpG_end",strand.field = "probe_strand",na.rm=T)
names(mani_gr)<-mani_gr$probeID
rownames(cb) <- cb$probeID
mcols(mani_gr)<-cb[names(mani_gr), ,on=.(probeID)]
mani_cpg <- makeGRangesFromDataFrame(mani,keep.extra.columns = T,ignore.strand = F,
                                     seqnames.field = "CpG_chrm",start.field = "CpG_beg",
                                     end.field = "CpG_end",strand.field = "probe_strand",na.rm=T)
names(mani_cpg)<-mani_cpg$probeID

mani_ug = fread(paste0("references/MSA.hg38.manifest.gencode.v41.tsv.gz")) %>% dplyr::select(CpG_chrm, CpG_beg, CpG_end, probe_strand,probeID,genesUniq)
#chr2:232418299-232418301 --- chr2:232420896-232420898  ECEL1P1
library(GenomicState)
txdb <- gencode_txdb("47", "hg38", chrs = "chr2")

ecel.region<-toGRanges("2:232416025-232419910",genome = "hg38")
ecel1p1 <- subsetByOverlaps(mani_gr, ecel.region)
ecel_probe <- ecel1p1$probeID

## We are going to collapse everything into CpGs from probeIDs
cb_ecel <- cb %>% dplyr::filter(probeID %in% ecel_probe) %>% dplyr::mutate(probeID=gsub("\\_.*","",probeID)) %>% group_by(probeID) %>% summarise(across(everything(), list(mean)))
rownames(cb_ecel) <- cb_ecel$probeID
mani_ecel <-mani_cpg
mcols(mani_ecel)<-cb[names(mani_ecel), ,on=.(probeID)]
mani_ecel<- subsetByOverlaps(mani_ecel, ecel.region)
mani_cpg<- subsetByOverlaps(mani_ecel, ecel.region)
## Annotate
mani_cpg<- annotatePeak(mani_cpg, tssRegion=c(-3000, 3000),
                        TxDb=txdb, annoDb="org.Hs.eg.db",level = 'gene') %>% as_tibble()


######
ss<- pheno %>% dplyr::select(Sample_ID,LH_smoke_status)
mcpg <- manifest %>% dplyr::mutate(cpg=gsub("\\_.*","",Probe_ID)) %>% dplyr::select(-Probe_ID)%>%
  dplyr::select(cpg,CpG_chrm,CpG_beg) %>% dplyr::distinct()
cb_ecel_pval <- cb_ecel %>% select(probeID,smokers_pvalue_1,formers_pvalue_1)
cb_ecel_betas <- cb_ecel %>% select(-smokers_pvalue_1,-formers_pvalue_1)
cb_ecel_betas_long <- cb_ecel_betas %>% tidyr::pivot_longer(cols=-probeID) %>%
  dplyr::mutate(name = substr(name,1,nchar(name)-2)) %>%
  dplyr::left_join(ss, by = c('name'='Sample_ID'))  %>% dplyr::filter(LH_smoke_status!='M') %>% 
  dplyr::left_join(mcpg,by=c('probeID'='cpg')) %>% dplyr::group_by(probeID) %>% arrange(CpG_beg)

names(cb_ecel_betas_long) <- c('cpg','sample','beta','smoke_status','chr','start')
custom.col <- c("","#7fcdbb","")
custom.col <- c("#db9a9a", "#a33636","#290101")
#cb_ecel_betas_long$smoke_status <- factor(cb_ecel_betas_long$smoke_status,levels=c('Never','Former','Current'))
#group.colors = c(Current = "#3261cf", Former = "#a33636",Never)
pdf('output/ecel1p1/dmr.pdf',width=12)
ggplot(cb_ecel_betas_long,aes(x=cpg,y=beta,dot.size=0.03,color=smoke_status,fill=smoke_status)) +
  geom_dotplot(binwidth=0.006,stackratio =0.5, stackdir='center',position='identity',binaxis = "y",stroke=NA, alpha=0.9) +
  geom_smooth(aes(group=smoke_status),method='loess') +
  labs(x = "CpG", y = bquote(atop('Beta Value'))) + 
  scale_fill_manual(values=custom.col,guide = "legend",name = "Smoking Status",aesthetics = c("colour", "fill"))+
  theme_classic() + ggtitle("Methylation beta values across ECEL1P1 on the MSA",)



ggplot(cb_ecel_betas_long,aes(x=cpg,y=beta,dot.size=0.03,color=smoke_status,fill=smoke_status)) +
  ggviolin(binwidth=0.006,stackratio =0.1, stackdir='centerwhole',position='dodge',binaxis = "y",stroke=NA, alpha=0.5) +
  geom_smooth(aes(group=smoke_status),method='loess') +
  labs(x = "CpG", y = bquote(atop('Beta Value'))) + 
  scale_color_discrete(guide = "legend",name = "Smoking Status",aesthetics = c("colour", "fill"))+
  theme_classic() + ggtitle("Methylation beta values across ECEL1P1 on the MSA")

dev.off()

save(cb_ecel_betas_long,file='output/ecel1p1/ecel1p1_beta_long.rdata')


