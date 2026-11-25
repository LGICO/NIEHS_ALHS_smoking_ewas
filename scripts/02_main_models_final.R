### Example script - run an epigenome-wide association analysis with weights
### Oct 2024
### Mikyeong Lee (NIEHS)


### methylation & phenotype data
library(dplyr)

load('751_clean_pc.rdata')
dir.create("output/svyglm", recursive=TRUE, showWarnings=FALSE)
dir.create("output/svyglm_anno", recursive=TRUE, showWarnings=FALSE)
dir.create("output/significant_models", recursive=TRUE, showWarnings=FALSE)
betap <- t(betas_miss_masked)
mvalp <- t(mvals_miss_masked)



### EWAS

library(survey)
library(MASS)
library(lmtest)
library(sandwich) 

library(parallel) 
#cl <- makeCluster(1, type="FORK")
cl <- makeCluster(getOption("cl.cores", 60))
clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(sandwich))
clusterEvalQ(cl,library(lmtest))
clusterEvalQ(cl,library(survey))
clusterEvalQ(cl,library(estimatr))
clusterEvalQ(cl,library(dplyr))
clusterEvalQ(cl,library(robustbase))


################################
##pc only
#pheno=pheno%>%mutate(asthma_weights=ifelse(LHSSCREEN_12142015=="Case",w_asthma,w_nonasthma))
#pheno = pheno %>% filter(Gender == inrdsx)
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
pheno <- pheno %>% mutate( ANYETS = case_when(ANYETS =="No"~0,ANYETS=="Yes"~1))
pheno <- pheno %>% mutate(LHSSCREEN_12142015 = case_when(LHSSCREEN_12142015 == "Control"~0,LHSSCREEN_12142015 == "Case"~1))
pheno <- pheno %>% mutate(Gender = case_when(Gender == 'male'~0,Gender=='female'~1))
pheno <- pheno %>% mutate(STATE = case_when(STATE == 'IA'~0,STATE=='NC'~1))
pheno$LHSSCREEN_12142015 = as.factor(pheno$LHSSCREEN_12142015)
pheno$ANYETS = as.factor(pheno$ANYETS)
pheno$STATE = as.factor(pheno$STATE)
pheno <- pheno %>% mutate(smokeXasthma = paste0(LH_smoke_status,'_',LHSSCREEN_12142015))
pheno$smokeXasthma <- as.factor(pheno$smokeXasthma)


pheno_pc <- pheno %>% filter(European == 1,!is.na(PC1))
beta_pc <- beta2[row.names(beta2) %in% pheno_pc$Sample_ID,]
betap_pc <- betap[row.names(betap) %in% pheno_pc$Sample_ID,]


pheno_pc_case <- pheno_pc %>% filter(LHSSCREEN_12142015 == 1)
beta_pc_case <- beta_pc[row.names(beta_pc) %in% pheno_pc_case$Sample_ID,]
betap_pc_case <- betap_pc[row.names(betap_pc) %in% pheno_pc_case$Sample_ID,]
pheno_pc_control <- pheno_pc %>% filter(LHSSCREEN_12142015 == 0)
beta_pc_control <- beta_pc[row.names(beta_pc) %in% pheno_pc_control$Sample_ID,]
betap_pc_control <- betap_pc[row.names(betap_pc) %in% pheno_pc_control$Sample_ID,]

pp <- pheno_pc %>% filter(LH_smoke_status != 'M')
bb <- beta_pc[row.names(beta_pc) %in% pp$Sample_ID,]
##########
##Model 1 unweighted euro
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
write.csv(results, "output/svyglm/unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
write.csv(results, "output/svyglm/unimputed_unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")




##Model 10
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
write.csv(results, "output/svyglm/unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.unwgt, pheno)
results <- as.data.frame(t(results))
write.csv(results, "output/svyglm/unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")

##########
##
##########################################################################################################################
##weighted
################
##smokers
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
pp$asthma_weights [pp$LHSSCREEN_12142015 == 0] <- w_nonasthma 
pp$asthma_weights [pp$LHSSCREEN_12142015 == 1] <- w_asthma 
###########
## Model 3
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
write.csv(results, "output/svyglm/weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc)
results <- as.data.frame(t(results))
write.csv(results, "output/svyglm/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

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
write.csv(results, "output/svyglm/weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.wgt, pheno)
results <- as.data.frame(t(results))

write.csv(results, "output/svyglm/unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")



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
write.csv(results, "output/svyglm/unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
write.csv(results, "output/svyglm/unimputed_unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

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
write.csv(results, "output/svyglm/unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.unwgt, pheno)
results <- as.data.frame(t(results))
write.csv(results, "output/svyglm/unimputed_unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")




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
write.csv(results, "output/svyglm/weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc)
results <- as.data.frame(t(results))
write.csv(results, "output/svyglm/unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

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
write.csv(results, "output/svyglm/weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.wgt, pheno)
results <- as.data.frame(t(results))

write.csv(results, "output/svyglm/unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv")



########################################################################################
#######################################################################################
## CASE CONTROL SENSITIVITY Models 4, 5, 13, 14

#########
##Model 4 controls euro
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 0)
  coef = cf$coefficients[3,"Estimate"]
  se = cf$coefficients[3,"Std. Error"]
  p.value = cf$coefficients[3,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta_pc_control, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_control)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/svyglm/unweighted_smokers_controls_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

results <- parApply(cl, betap_pc_control, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_control)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/svyglm/unimputed_unweighted_smokers_controls_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

##########
##Model 5 cases euro

results <- parApply(cl, beta_pc_case, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_case)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/svyglm/unweighted_smokers_cases_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

results <- parApply(cl, betap_pc_case, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_case)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/svyglm/unimputed_unweighted_smokers_cases_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

#########
##Model 13 controls trans
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + European  + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$LH_smoke_status == 0)
  coef = cf$coefficients[3,"Estimate"]
  se = cf$coefficients[3,"Std. Error"]
  p.value = cf$coefficients[3,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}
pheno_control <- pheno %>% filter(LHSSCREEN_12142015 == 0)
beta_control <- beta2[row.names(beta2) %in% pheno_control$Sample_ID,]
betap_control <- betap[row.names(betap) %in% pheno_control$Sample_ID,]
results <- parApply(cl, beta_control, 2, f.RLM.Adj.Robust.unwgt, pheno_control)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/svyglm/unweighted_smokers_controls_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_control, 2, f.RLM.Adj.Robust.unwgt, pheno_control)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/svyglm/unimputed_unweighted_smokers_controls_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
##########
##Model 14 cases transanc
pheno_case <- pheno %>% filter(LHSSCREEN_12142015== 1)
beta_case <- beta2[row.names(beta2) %in% pheno_case$Sample_ID,]
betap_case <- betap[row.names(betap) %in% pheno_case$Sample_ID,]
results <- parApply(cl, beta_case, 2, f.RLM.Adj.Robust.unwgt, pheno_case)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/svyglm/unweighted_smokers_cases_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_case, 2, f.RLM.Adj.Robust.unwgt, pheno_case)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/svyglm/unimputed_unweighted_smokers_cases_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
#########






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
fin<- system('ls output/svyglm/*.csv',intern=T)
bonf<- 0.05/230325 
tsv_path <- "references/MSA.hg38.manifest.gencode.v41.tsv.gz"
manifest = fread(tsv_path)

for(f in fin){
  uws <- fread(f)
  siguws <- uws %>% arrange(P)#%>% dplyr::filter(uws$P<=bonf)
  siguw_anno <- left_join(siguws,manifest,by=c("V1"="probeID"))
  f <- substring(f,40) %>% substring(.,1,nchar(.)-4)
  write.table(siguw_anno,file=paste0("output/svyglm_anno/",f,".tsv"),sep='\t',quote = F,na = "",row.names = F)
}     
fin<- system('ls output/svyglm_anno/*.tsv',intern=T)

for(f in fin){
  uws <- fread(f,header = T)
  print(f)
  colnames(uws)<- c('Probe_ID','Beta','SE','P','N','N.exp','N.unexp','CpG_chrm','CpG_beg','CpG_end','probe_strand','genesUniq','geneNames','transcriptTypes','transcriptIDs','distToTSS')
  siguws <- uws %>% dplyr::filter(uws$P<=bonf) %>% arrange(P) 
  f <- substring(f,45)
  write.table(siguws,file=paste0("output/significant_models/",f),sep='\t',quote = F,na = "",row.names = F)
}    



