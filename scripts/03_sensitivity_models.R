### methylation & phenotype data
library(dplyr)
load('751_clean_pc.rdata')
dir.create("output/model_output/svyglm", recursive=TRUE, showWarnings=FALSE)
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

###############################################################################
## Smokers
##########

#########
##Model 2 cross euro
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + LHSSCREEN_12142015 + LH_smoke_status:LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[3,"Estimate"]
  se = cf$coefficients[3,"Std. Error"]
  p.value = cf$coefficients[3,"Pr(>|t|)"]
  LH_smoke_status1_LHSSCREEN_121420151_pval = cf$coefficients[27,"Pr(>|t|)"]
  LH_smoke_status2_LHSSCREEN_121420151_pval = cf$coefficients[28,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp,LH_smoke_status1_LHSSCREEN_121420151_pval,LH_smoke_status2_LHSSCREEN_121420151_pval) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp","LH_smoke_status1_LHSSCREEN_121420151_pval","LH_smoke_status2_LHSSCREEN_121420151_pval"); return(out) 
}

results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_smokers_EWAS_Results_europc_only_smoke_gender_crossprod_BMI_age_state_cells_pcs.csv")

results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_EWAS_Results_europc_only_smoke_gender_crossprod_BMI_age_state_cells_pcs.csv")

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
write.csv(results, "output/model_output/svyglm/unweighted_smokers_controls_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

results <- parApply(cl, betap_pc_control, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_control)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_controls_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

##########
##Model 5 cases euro

results <- parApply(cl, beta_pc_case, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_case)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_smokers_cases_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

results <- parApply(cl, betap_pc_case, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_case)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_cases_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

#########
##Model 6 men euro
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + LHSSCREEN_12142015 + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[3,"Estimate"]
  se = cf$coefficients[3,"Std. Error"]
  p.value = cf$coefficients[3,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}

pheno_pc_men <- pheno_pc %>% filter(Gender == 0)
beta_pc_men <- beta_pc[row.names(beta_pc) %in% pheno_pc_men$Sample_ID,]
betap_pc_men <- betap_pc[row.names(betap_pc) %in% pheno_pc_men$Sample_ID,]
results <- parApply(cl, beta_pc_men, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_men)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_smokers_men_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_pc_men, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_men)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_men_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
########
##Model 7 women euro
pheno_pc_women <- pheno_pc %>% filter(Gender == 1)
beta_pc_women <- beta_pc[row.names(beta_pc) %in% pheno_pc_women$Sample_ID,]
betap_pc_women <- betap_pc[row.names(betap_pc) %in% pheno_pc_women$Sample_ID,]
results <- parApply(cl, beta_pc_women, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_women)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_smokers_women_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_pc_women, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_women)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_women_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

##########
##Model 8 Iowa euro
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[3,"Estimate"]
  se = cf$coefficients[3,"Std. Error"]
  p.value = cf$coefficients[3,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}
pheno_pc_state <- pheno_pc %>% filter(STATE == 0)
beta_pc_state <- beta_pc[row.names(beta_pc) %in% pheno_pc_state$Sample_ID,]
betap_pc_state <- betap_pc[row.names(betap_pc) %in% pheno_pc_state$Sample_ID,]
results <- parApply(cl, beta_pc_state, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_smokers_iowa_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_pc_state, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_iowa_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

##########
##Model 9 North carolina euro
pheno_pc_state <- pheno_pc %>% filter(STATE == 1)
beta_pc_state <- beta_pc[row.names(beta_pc) %in% pheno_pc_state$Sample_ID,]
betap_pc_state <- betap_pc[row.names(betap_pc) %in% pheno_pc_state$Sample_ID,]
results <- parApply(cl, beta_pc_state, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_smokers_nocarolina_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_pc_state, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_nocarolina_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

##########
##Model 11 cross-prod
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + European + LHSSCREEN_12142015 + LH_smoke_status:LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[3,"Estimate"]
  se = cf$coefficients[3,"Std. Error"]
  p.value = cf$coefficients[3,"Pr(>|t|)"]
  LH_smoke_status1_LHSSCREEN_121420151_pval = cf$coefficients[27,"Pr(>|t|)"]
  LH_smoke_status2_LHSSCREEN_121420151_pval = cf$coefficients[28,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp,LH_smoke_status1_LHSSCREEN_121420151_pval,LH_smoke_status2_LHSSCREEN_121420151_pval) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp","LH_smoke_status1_LHSSCREEN_121420151_pval","LH_smoke_status2_LHSSCREEN_121420151_pval"); return(out) 
}

results <- parApply(cl, beta2, 2, f.RLM.Adj.Robust.unwgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_smokers_EWAS_Results_transanc_smoke_gender_crossprod_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap2, 2, f.RLM.Adj.Robust.unwgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_crossprod_BMI_age_state_cells_pcs.csv")
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
write.csv(results, "output/model_output/svyglm/unweighted_smokers_controls_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_control, 2, f.RLM.Adj.Robust.unwgt, pheno_control)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_controls_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
##########
##Model 14 cases transanc
pheno_case <- pheno %>% filter(LHSSCREEN_12142015== 1)
beta_case <- beta2[row.names(beta2) %in% pheno_case$Sample_ID,]
betap_case <- betap[row.names(betap) %in% pheno_case$Sample_ID,]
results <- parApply(cl, beta_case, 2, f.RLM.Adj.Robust.unwgt, pheno_case)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_smokers_cases_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_case, 2, f.RLM.Adj.Robust.unwgt, pheno_case)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_cases_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
#########
##Model 15 men transanc
pheno_men <- pheno %>% filter(Gender == 0)
beta_men <- beta2[row.names(beta2) %in% pheno_men$Sample_ID,]
betap_men <- betap[row.names(betap) %in% pheno_men$Sample_ID,]

f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + European  + LHSSCREEN_12142015 + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[3,"Estimate"]
  se = cf$coefficients[3,"Std. Error"]
  p.value = cf$coefficients[3,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}


results <- parApply(cl, beta_men, 2, f.RLM.Adj.Robust.unwgt, pheno_men)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_smokers_men_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_men, 2, f.RLM.Adj.Robust.unwgt, pheno_men)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_men_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
########
##Model 16 women transanc
pheno_women <- pheno %>% filter(Gender == 1)
beta_women <- beta2[row.names(beta2) %in% pheno_women$Sample_ID,]
betap_women <- betap[row.names(betap) %in% pheno_women$Sample_ID,]
results <- parApply(cl, beta_women, 2, f.RLM.Adj.Robust.unwgt, pheno_women)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_smokers_women_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_women, 2, f.RLM.Adj.Robust.unwgt, pheno_women)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_women_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

##########
##Model 17 Iowa transanc
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + European  + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[3,"Estimate"]
  se = cf$coefficients[3,"Std. Error"]
  p.value = cf$coefficients[3,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}
pheno_state <- pheno %>% filter(STATE == 0)
beta_state <- beta2[row.names(beta2) %in% pheno_state$Sample_ID,]
betap_state <- betap[row.names(betap) %in% pheno_state$Sample_ID,]
results <- parApply(cl, beta_state, 2, f.RLM.Adj.Robust.unwgt, pheno_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_smokers_iowa_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_state, 2, f.RLM.Adj.Robust.unwgt, pheno_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_iowa_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
##########
##Model 18 North carolina transanc
pheno_state <- pheno %>% filter(STATE == 1)
beta_state <- beta2[row.names(beta2) %in% pheno_state$Sample_ID,]
betap_state <- betap[row.names(betap) %in% pheno_state$Sample_ID,]
results <- parApply(cl, beta_state, 2, f.RLM.Adj.Robust.unwgt, pheno_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_smokers_nocarolina_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_state, 2, f.RLM.Adj.Robust.unwgt, pheno_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_smokers_nocarolina_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
##########################################################################################################################
##weights
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



##########################################################################################################################################
## former-smokers
##########
##########################################################################################################################################

#########
##Model 2 cross euro
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + LHSSCREEN_12142015 + LH_smoke_status:LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[2,"Estimate"]
  se = cf$coefficients[2,"Std. Error"]
  p.value = cf$coefficients[2,"Pr(>|t|)"]
  LH_smoke_status1_LHSSCREEN_121420151_pval = cf$coefficients[27,"Pr(>|t|)"]
  LH_smoke_status2_LHSSCREEN_121420151_pval = cf$coefficients[28,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp,LH_smoke_status1_LHSSCREEN_121420151_pval,LH_smoke_status2_LHSSCREEN_121420151_pval) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp","LH_smoke_status1_LHSSCREEN_121420151_pval","LH_smoke_status2_LHSSCREEN_121420151_pval"); return(out) 
}

results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_crossprod_BMI_age_state_cells_pcs.csv")

results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_crossprod_BMI_age_state_cells_pcs.csv")

#########
##Model 4 controls euro
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[2,"Estimate"]
  se = cf$coefficients[2,"Std. Error"]
  p.value = cf$coefficients[2,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta_pc_control, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_control)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_controls_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

results <- parApply(cl, betap_pc_control, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_control)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_controls_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

##########
##Model 5 cases euro

results <- parApply(cl, beta_pc_case, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_case)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_cases_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

results <- parApply(cl, betap_pc_case, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_case)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_cases_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

#########
##Model 6 men euro
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + LHSSCREEN_12142015 + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[2,"Estimate"]
  se = cf$coefficients[2,"Std. Error"]
  p.value = cf$coefficients[2,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}

pheno_pc_men <- pheno_pc %>% filter(Gender == 0)
beta_pc_men <- beta_pc[row.names(beta_pc) %in% pheno_pc_men$Sample_ID,]
betap_pc_men <- betap_pc[row.names(betap_pc) %in% pheno_pc_men$Sample_ID,]
results <- parApply(cl, beta_pc_men, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_men)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_men_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_pc_men, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_men)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_men_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
########
##Model 7 women euro
pheno_pc_women <- pheno_pc %>% filter(Gender == 1)
beta_pc_women <- beta_pc[row.names(beta_pc) %in% pheno_pc_women$Sample_ID,]
betap_pc_women <- betap_pc[row.names(betap_pc) %in% pheno_pc_women$Sample_ID,]
results <- parApply(cl, beta_pc_women, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_women)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_women_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_pc_women, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_women)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_women_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

##########
##Model 8 Iowa euro
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[2,"Estimate"]
  se = cf$coefficients[2,"Std. Error"]
  p.value = cf$coefficients[2,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}
pheno_pc_state <- pheno_pc %>% filter(STATE == 0)
beta_pc_state <- beta_pc[row.names(beta_pc) %in% pheno_pc_state$Sample_ID,]
betap_pc_state <- betap_pc[row.names(betap_pc) %in% pheno_pc_state$Sample_ID,]
results <- parApply(cl, beta_pc_state, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_iowa_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_pc_state, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_iowa_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

##########
##Model 9 North carolina euro
pheno_pc_state <- pheno_pc %>% filter(STATE == 1)
beta_pc_state <- beta_pc[row.names(beta_pc) %in% pheno_pc_state$Sample_ID,]
betap_pc_state <- betap_pc[row.names(betap_pc) %in% pheno_pc_state$Sample_ID,]
results <- parApply(cl, beta_pc_state, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_nocarolina_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_pc_state, 2, f.RLM.Adj.Robust.unwgt, pheno_pc_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_nocarolina_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv")


##########
##Model 11 cross-prod
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + European + LHSSCREEN_12142015 + LH_smoke_status:LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[2,"Estimate"]
  se = cf$coefficients[2,"Std. Error"]
  p.value = cf$coefficients[2,"Pr(>|t|)"]
  LH_smoke_status1_LHSSCREEN_121420151_pval = cf$coefficients[27,"Pr(>|t|)"]
  LH_smoke_status2_LHSSCREEN_121420151_pval = cf$coefficients[28,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp,LH_smoke_status1_LHSSCREEN_121420151_pval,LH_smoke_status2_LHSSCREEN_121420151_pval) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp","LH_smoke_status1_LHSSCREEN_121420151_pval","LH_smoke_status2_LHSSCREEN_121420151_pval"); return(out) 
}

results <- parApply(cl, beta2, 2, f.RLM.Adj.Robust.unwgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_crossprod_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap2, 2, f.RLM.Adj.Robust.unwgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_crossprod_BMI_age_state_cells_pcs.csv")
#########
##Model 13 controls trans
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + European  + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[2,"Estimate"]
  se = cf$coefficients[2,"Std. Error"]
  p.value = cf$coefficients[2,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}
pheno_control <- pheno %>% filter(LHSSCREEN_12142015 == 0)
beta_control <- beta2[row.names(beta2) %in% pheno_control$Sample_ID,]
betap_control <- betap[row.names(betap) %in% pheno_control$Sample_ID,]
results <- parApply(cl, beta_control, 2, f.RLM.Adj.Robust.unwgt, pheno_control)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_controls_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_control, 2, f.RLM.Adj.Robust.unwgt, pheno_control)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_controls_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
##########
##Model 14 cases transanc
pheno_case <- pheno %>% filter(LHSSCREEN_12142015== 1)
beta_case <- beta2[row.names(beta2) %in% pheno_case$Sample_ID,]
betap_case <- betap[row.names(betap) %in% pheno_case$Sample_ID,]
results <- parApply(cl, beta_case, 2, f.RLM.Adj.Robust.unwgt, pheno_case)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_cases_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_case, 2, f.RLM.Adj.Robust.unwgt, pheno_case)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_cases_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
#########
##Model 15 men transanc
pheno_men <- pheno %>% filter(Gender == 0)
beta_men <- beta2[row.names(beta2) %in% pheno_men$Sample_ID,]
betap_men <- betap[row.names(betap) %in% pheno_men$Sample_ID,]

f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + European  + LHSSCREEN_12142015 + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[2,"Estimate"]
  se = cf$coefficients[2,"Std. Error"]
  p.value = cf$coefficients[2,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}


results <- parApply(cl, beta_men, 2, f.RLM.Adj.Robust.unwgt, pheno_men)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_men_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_men, 2, f.RLM.Adj.Robust.unwgt, pheno_men)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_men_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
########
##Model 16 women transanc
pheno_women <- pheno %>% filter(Gender == 1)
beta_women <- beta2[row.names(beta2) %in% pheno_women$Sample_ID,]
betap_women <- betap[row.names(betap) %in% pheno_women$Sample_ID,]
results <- parApply(cl, beta_women, 2, f.RLM.Adj.Robust.unwgt, pheno_women)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_women_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_women, 2, f.RLM.Adj.Robust.unwgt, pheno_women)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_women_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

##########
##Model 17 Iowa transanc
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_smoke_status + ANYETS + European  + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 2)
  n.unexp <- sum(!is.na(dat$y) & dat$LH_smoke_status == 0)
  coef = cf$coefficients[2,"Estimate"]
  se = cf$coefficients[2,"Std. Error"]
  p.value = cf$coefficients[2,"Pr(>|t|)"]
  out = c(coef, se, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "P", "N", "N.exp", "N.unexp"); return(out) 
}
pheno_state <- pheno %>% filter(STATE == 0)
beta_state <- beta2[row.names(beta2) %in% pheno_state$Sample_ID,]
betap_state <- betap[row.names(betap) %in% pheno_state$Sample_ID,]
results <- parApply(cl, beta_state, 2, f.RLM.Adj.Robust.unwgt, pheno_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_iowa_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_state, 2, f.RLM.Adj.Robust.unwgt, pheno_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_iowa_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
##########
##Model 18 North carolina transanc
pheno_state <- pheno %>% filter(STATE == 1)
beta_state <- beta2[row.names(beta2) %in% pheno_state$Sample_ID,]
betap_state <- betap[row.names(betap) %in% pheno_state$Sample_ID,]
results <- parApply(cl, beta_state, 2, f.RLM.Adj.Robust.unwgt, pheno_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unweighted_formersmokers_nocarolina_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")
results <- parApply(cl, betap_state, 2, f.RLM.Adj.Robust.unwgt, pheno_state)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/model_output/svyglm/unimputed_unweighted_formersmokers_nocarolina_EWAS_Results_transanc_only_smoke_gender_BMI_age_state_cells_pcs.csv")

stopCluster(cl)



