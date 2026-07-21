### run an epigenome-wide association analysis with and without weights

### methylation & phenotype data
library(dplyr)
load('output/rdata/cleaned_MSA_niehs_alhs.rdata')
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
## Prepare phenotype data for models

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

## Filter data
pheno_pc <- pheno %>% filter(European == 1,!is.na(PC1))
beta_pc <- beta2[row.names(beta2) %in% pheno_pc$Sample_ID,]
betap_pc <- betap[row.names(betap) %in% pheno_pc$Sample_ID,]


pheno_pc_case <- pheno_pc %>% filter(LHSSCREEN_12142015 == 1)
beta_pc_case <- beta_pc[row.names(beta_pc) %in% pheno_pc_case$Sample_ID,]
betap_pc_case <- betap_pc[row.names(betap_pc) %in% pheno_pc_case$Sample_ID,]
pheno_pc_control <- pheno_pc %>% filter(LHSSCREEN_12142015 == 0)
beta_pc_control <- beta_pc[row.names(beta_pc) %in% pheno_pc_control$Sample_ID,]
betap_pc_control <- betap_pc[row.names(betap_pc) %in% pheno_pc_control$Sample_ID,]

###############################################################################
## Run Main Models
system("mkdir -p output/model_output")
system("mkdir -p output/model_output/svyglm")
system("mkdir -p output/model_output/svyglm_anno")
system("mkdir -p output/model_output/svyglm_significant_cpgs")


##########
##Model 1 unweighted europeans
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf = summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 0)
  coef = cf$coefficients[3,"Estimate"]
  se = cf$coefficients[3,"Std. Error"]
  p.value = cf$coefficients[3,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")




##Model 10 unweighted transancestry
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + European + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf = coeftest(mdl, vcov=vcov(mdl))
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 2)
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
write.csv(results, "output/model_output/svyglm/unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.unwgt, pheno)
results <- as.data.frame(t(results))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")

##########
##
##########################################################################################################################
##weighted
################
##smokers

## set weights
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
## Model 3 weighted european
f.RLM.Adj.Robust.wgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, weights = ~asthma_weights, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS  + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf = coeftest(mdl, vcov=vcov(mdl))
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 0)
  coef = cf[3,"Estimate"]
  se = cf[3,"Std. Error"]
  p.value = cf[3,"Pr(>|z|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc)
results <- as.data.frame(t(results))
write.csv(results, "output/model_output/svyglm/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

##############################
##Model 12 weighted transancestry
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
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 0)
  coef = cf[3,"Estimate"]
  se = cf[3,"Std. Error"]
  p.value = cf[3,"Pr(>|z|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta2, 2, f.RLM.Adj.Robust.wgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.wgt, pheno)
results <- as.data.frame(t(results))

write.csv(results, "output/model_output/svyglm/unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")



##########################################################################################################################################
## former-smokers
##########
##########################################################################################################################################

##Model 1 unweighted european
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

##Model 10 unweighted transancestry
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
## Model 3 weighted european
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
##Model 12 weighted transancestry 
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

stopCluster(cl)




