library(dplyr)
load('/output/rdata/cleaned_MSA_niehs_alhs.rdata')
beta2 = t(clean_betas) # methylation data
pheno = clean_sampsheet
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



################################
##prepare data

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
pheno$LH_LastSmkedLE5YrAgo_v1 <-factor(pheno$LH_LastSmkedLE5YrAgo_v1, levels = c("0) LE 5 yrs ago","1) GT 5 yrs ago"))

pheno = pheno %>% mutate(smoke_statusF = case_when(LH_smoke_status == 1 & LH_LastSmkedLE5YrAgo_v1=="0) LE 5 yrs ago" ~ factor(3),
                                                   LH_smoke_status == 1 & LH_LastSmkedLE5YrAgo_v1=="1) GT 5 yrs ago" ~ factor(4),
                                                   .default = LH_smoke_status))
pheno$smoke_statusF <- as.factor(pheno$smoke_statusF) %>% droplevels(.,exclude=1)
pheno$smoke_statusF <- factor(pheno$smoke_statusF, levels = c(0,2,3,4,99))

pheno <- pheno %>% mutate(former_5yrs_split = case_when(LH_smoke_status == 0 ~ factor(0),
                                               LH_smoke_status == 2 ~ factor(1),
                                               LH_smoke_status == 1 & LH_TimeSinceSmked_v1 <= 5 ~ factor(2),
                                             LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 5 & LH_TimeSinceSmked_v1 <= 10 ~ factor(3),
                                             LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 10 & LH_TimeSinceSmked_v1 <= 15 ~ factor(4),
                                             LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 15 ~ factor(5),
                                             LH_smoke_status == 99 ~ factor(99)))
pheno <- pheno %>% mutate(former_10yrs_split = case_when(LH_smoke_status == 0 ~ factor(0),
                                                       LH_smoke_status == 2 ~ factor(1),
                                                       LH_smoke_status == 1 & LH_TimeSinceSmked_v1 <= 10 ~ factor(2),
                                                       LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 10 & LH_TimeSinceSmked_v1 <= 20 ~ factor(3),
                                                       LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 20 & LH_TimeSinceSmked_v1 <= 30 ~ factor(4),
                                                       LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 30 ~ factor(5),
                                                       LH_smoke_status == 99 ~ factor(99)))
##need to keep never smokers

pheno_pc <- pheno %>% filter(European == 1,!is.na(PC1))
beta_pc <- beta2[row.names(beta2) %in% pheno_pc$Sample_ID,]
betap_pc <- betap[row.names(betap) %in% pheno_pc$Sample_ID,]


pheno_pc_case <- pheno_pc %>% filter(LHSSCREEN_12142015 == 1)
beta_pc_case <- beta_pc[row.names(beta_pc) %in% pheno_pc_case$Sample_ID,]
betap_pc_case <- betap_pc[row.names(betap_pc) %in% pheno_pc_case$Sample_ID,]
pheno_pc_control <- pheno_pc %>% filter(LHSSCREEN_12142015 == 0)
beta_pc_control <- beta_pc[row.names(beta_pc) %in% pheno_pc_control$Sample_ID,]
betap_pc_control <- betap_pc[row.names(betap_pc) %in% pheno_pc_control$Sample_ID,]
pheno_former <- pheno %>% filter(LH_smoke_status==1)
beta_former <- beta2[row.names(beta2) %in% pheno_former$Sample_ID,]
pheno_formerpc <- pheno %>% filter(LH_smoke_status==1,European == 1,!is.na(PC1))
beta_formerpc <- beta2[row.names(beta2) %in% pheno_formerpc$Sample_ID,]
#######################################################################################
##########
##Model 1 unweighted euro
## intercept 1
## smoke_2 2
## smoke_3 3
## smoke_4 4
## smoke_99 5
system('mkdir -p output/former_strat')
##########
##Model 1 unweighted euro
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ smoke_statusF + ANYETS + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 3)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 0)
  coef = cf$coefficients[3,"Estimate"]
  se = cf$coefficients[3,"Std. Error"]
  ts = cf$coefficients[3,"t value"]
  p.value = cf$coefficients[3,"Pr(>|t|)"]
  out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_strat/model1_unweighted.3.csv")

f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ smoke_statusF + ANYETS + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 4)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 0)
  coef = cf$coefficients[4,"Estimate"]
  se = cf$coefficients[4,"Std. Error"]
  ts = cf$coefficients[4,"t value"]
  p.value = cf$coefficients[4,"Pr(>|t|)"]
  out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_strat/model1_unweighted.4.csv")
##########
##Model 10
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ smoke_statusF + ANYETS + European + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 3)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 0)
  coef = cf$coefficients[3,"Estimate"]
  se = cf$coefficients[3,"Std. Error"]
  ts = cf$coefficients[3,"t value"]
  p.value = cf$coefficients[3,"Pr(>|t|)"]
  out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta2, 2, f.RLM.Adj.Robust.unwgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_strat/model10_unweighted.3.csv")
##########
##Model 10
f.RLM.Adj.Robust.unwgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ smoke_statusF + ANYETS + European + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 4)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 0)
  coef = cf$coefficients[4,"Estimate"]
  se = cf$coefficients[4,"Std. Error"]
  ts = cf$coefficients[4,"t value"]
  p.value = cf$coefficients[4,"Pr(>|t|)"]
  out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta2, 2, f.RLM.Adj.Robust.unwgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_strat/model10_unweighted.4.csv")


##########################################################################################################################
##weighted
################
##smokers
p_asthma_overall = 0.06852481 # n_asthma_overall/(n_asthma_overall+n_nonasthma_overall) 
n_asthma = 186
n_nonasthma = 4579
p_asthma_nested = n_asthma/(n_asthma+n_nonasthma) # n_asthma_nested/(n_asthma_nested+n_nonasthma_nested)   
p_nonasthma_overall = 1- p_asthma_overall # n_nonasthma_overall/(n_asthma_overall+n_nonasthma_overall) 
p_nonasthma_nested = 1- p_asthma_nested # n_nonasthma_nested/(n_asthma_nested+n_nonasthma_nested) 
w_asthma <- 1/(p_asthma_nested/p_asthma_overall) 
w_nonasthma <- 1/(p_nonasthma_nested/p_nonasthma_overall) 
pheno$asthma_weights [pheno$LHSSCREEN_12142015 == 0] <- w_nonasthma 
pheno$asthma_weights [pheno$LHSSCREEN_12142015 == 1] <- w_asthma 
pheno_pc$asthma_weights [pheno_pc$LHSSCREEN_12142015 == 0] <- w_nonasthma 
pheno_pc$asthma_weights [pheno_pc$LHSSCREEN_12142015 == 1] <- w_asthma 


###########
## Model 3
f.RLM.Adj.Robust.wgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, weights = ~asthma_weights, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ smoke_statusF + ANYETS + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 3)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 0)
  coef = cf$coefficients[3,"Estimate"]
  se = cf$coefficients[3,"Std. Error"]
  ts = cf$coefficients[3,"t value"]
  p.value = cf$coefficients[3,"Pr(>|t|)"]
  out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_strat/model3_weighted.3.csv")

###########
## Model 3
f.RLM.Adj.Robust.wgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, weights = ~asthma_weights, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ smoke_statusF + ANYETS + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 4)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 0)
  coef = cf$coefficients[4,"Estimate"]
  se = cf$coefficients[4,"Std. Error"]
  ts = cf$coefficients[4,"t value"]
  p.value = cf$coefficients[4,"Pr(>|t|)"]
  out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_strat/model3_weighted.4.csv")

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
  mdl.formula <- formula(paste("y ~ smoke_statusF + ANYETS + European + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 3)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 0)
  coef = cf$coefficients[3,"Estimate"]
  se = cf$coefficients[3,"Std. Error"]
  ts = cf$coefficients[3,"t value"]
  p.value = cf$coefficients[3,"Pr(>|t|)"]
  out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta2, 2, f.RLM.Adj.Robust.wgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_strat/model12_weighted.3.csv")

f.RLM.Adj.Robust.wgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, weights = ~asthma_weights, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ smoke_statusF + ANYETS + European + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  n.exp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 4)
  n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$smoke_statusF == 0)
  coef = cf$coefficients[4,"Estimate"]
  se = cf$coefficients[4,"Std. Error"]
  ts = cf$coefficients[4,"t value"]
  p.value = cf$coefficients[4,"Pr(>|t|)"]
  out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
  names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
}

results <- parApply(cl, beta2, 2, f.RLM.Adj.Robust.wgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_strat/model12_weighted.4.csv")


################################################################################################################
## Split by 0-5,6-10,11-15,>15

##Model 1 unweighted euro
## intercept 1
## current 2
## 5yr 3
## 6-10 4
## 11-15 5
## >15 6


##Model 1 unweighted euro
system('mkdir -p output/five_year')
for (i in 3:6) {
  f.RLM.Adj.Robust.unwgt <- function(meth, dat1,j) {
    dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
    mdl.formula <- formula(paste("y ~ former_5yrs_split + ANYETS + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
    mdl <- svyglm(formula = mdl.formula, design = dat)
    cf <- summary(mdl)
    N <- sum(!is.na(meth))
    n.exp <- sum(!is.na(dat$variables$y) & dat$variables$former_5yrs_split == j)
    n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$former_5yrs_split == 0)
    coef = cf$coefficients[j,"Estimate"]
    se = cf$coefficients[j,"Std. Error"]
    ts = cf$coefficients[j,"t value"]
    p.value = cf$coefficients[j,"Pr(>|t|)"]
    out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
    names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
  }
  
  results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc,i)
  results <- as.data.frame(t(results))
  Nbeta <- sum(!is.na(results$P))
  write.csv(results, paste0("output/five_year/model1_unweighted.",i,".csv"))
}

##########
##Model 10
for (i in 3:6) {
  f.RLM.Adj.Robust.unwgt <- function(meth, dat1,j) {
    dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
    mdl.formula <- formula(paste("y ~ former_5yrs_split + ANYETS + European + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
    mdl <- svyglm(formula = mdl.formula, design = dat)
    cf <- summary(mdl)
    N <- sum(!is.na(meth))
    n.exp <- sum(!is.na(dat$variables$y) & dat$variables$former_5yrs_split == j)
    n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$former_5yrs_split == 0)
    coef = cf$coefficients[j,"Estimate"]
    se = cf$coefficients[j,"Std. Error"]
    ts = cf$coefficients[j,"t value"]
    p.value = cf$coefficients[j,"Pr(>|t|)"]
    out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
    names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
  }
  
  results <- parApply(cl, beta2, 2, f.RLM.Adj.Robust.unwgt, pheno,i)
  results <- as.data.frame(t(results))
  Nbeta <- sum(!is.na(results$P))
  write.csv(results, paste0("output/five_year/model10_unweighted.",i,".csv"))
}



##########################################################################################################################
##weighted
################
##smokers
p_asthma_overall = 0.06852481 # n_asthma_overall/(n_asthma_overall+n_nonasthma_overall) 
n_asthma = 186
n_nonasthma = 4579
p_asthma_nested = n_asthma/(n_asthma+n_nonasthma) # n_asthma_nested/(n_asthma_nested+n_nonasthma_nested)   
p_nonasthma_overall = 1- p_asthma_overall # n_nonasthma_overall/(n_asthma_overall+n_nonasthma_overall) 
p_nonasthma_nested = 1- p_asthma_nested # n_nonasthma_nested/(n_asthma_nested+n_nonasthma_nested) 
w_asthma <- 1/(p_asthma_nested/p_asthma_overall) 
w_nonasthma <- 1/(p_nonasthma_nested/p_nonasthma_overall) 
pheno$asthma_weights [pheno$LHSSCREEN_12142015 == 0] <- w_nonasthma 
pheno$asthma_weights [pheno$LHSSCREEN_12142015 == 1] <- w_asthma 
pheno_pc$asthma_weights [pheno_pc$LHSSCREEN_12142015 == 0] <- w_nonasthma 
pheno_pc$asthma_weights [pheno_pc$LHSSCREEN_12142015 == 1] <- w_asthma 
pheno_former$asthma_weights [pheno_former$LHSSCREEN_12142015 == 0] <- w_nonasthma 
pheno_former$asthma_weights [pheno_former$LHSSCREEN_12142015 == 1] <- w_asthma 
pheno_formerpc$asthma_weights [pheno_formerpc$LHSSCREEN_12142015 == 0] <- w_nonasthma 
pheno_formerpc$asthma_weights [pheno_formerpc$LHSSCREEN_12142015 == 1] <- w_asthma 


###########
## Model 3
for (i in 3:6) {
  f.RLM.Adj.Robust.wgt <- function(meth, dat1,j) {
    dat <- svydesign(id = ~Sample_ID, weights = ~asthma_weights, data = data.frame(y = meth, dat1))
    mdl.formula <- formula(paste("y ~ former_5yrs_split + ANYETS + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
    mdl <- svyglm(formula = mdl.formula, design = dat)
    cf <- summary(mdl)
    N <- sum(!is.na(meth))
    n.exp <- sum(!is.na(dat$variables$y) & dat$variables$former_5yrs_split == j)
    n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$former_5yrs_split == 0)
    coef = cf$coefficients[j,"Estimate"]
    se = cf$coefficients[j,"Std. Error"]
    ts = cf$coefficients[j,"t value"]
    p.value = cf$coefficients[j,"Pr(>|t|)"]
    out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
    names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
  }
  
  results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc,i)
  results <- as.data.frame(t(results))
  Nbeta <- sum(!is.na(results$P))
  write.csv(results, paste0("output/five_year/model3_weighted.",i,".csv"))
}

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
for (i in 3:6) {
  f.RLM.Adj.Robust.wgt <- function(meth, dat1,j) {
    dat <- svydesign(id = ~Sample_ID, weights = ~asthma_weights, data = data.frame(y = meth, dat1))
    mdl.formula <- formula(paste("y ~ former_5yrs_split + ANYETS + European + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
    mdl <- svyglm(formula = mdl.formula, design = dat)
    cf <- summary(mdl)
    N <- sum(!is.na(meth))
    n.exp <- sum(!is.na(dat$variables$y) & dat$variables$former_5yrs_split == j)
    n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$former_5yrs_split == 0)
    coef = cf$coefficients[j,"Estimate"]
    se = cf$coefficients[j,"Std. Error"]
    ts = cf$coefficients[j,"t value"]
    p.value = cf$coefficients[j,"Pr(>|t|)"]
    out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
    names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
  }
  
  results <- parApply(cl, beta2, 2, f.RLM.Adj.Robust.wgt, pheno,i)
  results <- as.data.frame(t(results))
  Nbeta <- sum(!is.na(results$P))
  write.csv(results, paste0("output/five_year/model12_weighted.",i,".csv"))
}

################################################################################################################
## Split by 0-10,11-20,21-30,>30

##Model 1 unweighted euro
## intercept 1
## current 2
## 10yr 3
## 11-20 4
## 21-30 5
## >30 6


##Model 1 unweighted euro
system('mkdir -p ten_year')
for (i in 3:6) {
  f.RLM.Adj.Robust.unwgt <- function(meth, dat1,j) {
    dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
    mdl.formula <- formula(paste("y ~ former_10yrs_split + ANYETS + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
    mdl <- svyglm(formula = mdl.formula, design = dat)
    cf <- summary(mdl)
    N <- sum(!is.na(meth))
    n.exp <- sum(!is.na(dat$variables$y) & dat$variables$former_10yrs_split == j)
    n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$former_10yrs_split == 0)
    coef = cf$coefficients[j,"Estimate"]
    se = cf$coefficients[j,"Std. Error"]
    ts = cf$coefficients[j,"t value"]
    p.value = cf$coefficients[j,"Pr(>|t|)"]
    out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
    names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
  }
  
  results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc,i)
  results <- as.data.frame(t(results))
  Nbeta <- sum(!is.na(results$P))
  write.csv(results, paste0("output/ten_year/model1_unweighted.",i,".csv"))
}

##########
##Model 10
for (i in 3:6) {
  f.RLM.Adj.Robust.unwgt <- function(meth, dat1,j) {
    dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
    mdl.formula <- formula(paste("y ~ former_10yrs_split + ANYETS + European + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
    mdl <- svyglm(formula = mdl.formula, design = dat)
    cf <- summary(mdl)
    N <- sum(!is.na(meth))
    n.exp <- sum(!is.na(dat$variables$y) & dat$variables$former_10yrs_split == j)
    n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$former_10yrs_split == 0)
    coef = cf$coefficients[j,"Estimate"]
    se = cf$coefficients[j,"Std. Error"]
    ts = cf$coefficients[j,"t value"]
    p.value = cf$coefficients[j,"Pr(>|t|)"]
    out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
    names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
  }
  
  results <- parApply(cl, beta2, 2, f.RLM.Adj.Robust.unwgt, pheno,i)
  results <- as.data.frame(t(results))
  Nbeta <- sum(!is.na(results$P))
  write.csv(results, paste0("output/ten_year/model10_unweighted.",i,".csv"))
}



##########################################################################################################################
##weighted
################
##smokers
p_asthma_overall = 0.06852481 # n_asthma_overall/(n_asthma_overall+n_nonasthma_overall) 
n_asthma = 186
n_nonasthma = 4579
p_asthma_nested = n_asthma/(n_asthma+n_nonasthma) # n_asthma_nested/(n_asthma_nested+n_nonasthma_nested)   
p_nonasthma_overall = 1- p_asthma_overall # n_nonasthma_overall/(n_asthma_overall+n_nonasthma_overall) 
p_nonasthma_nested = 1- p_asthma_nested # n_nonasthma_nested/(n_asthma_nested+n_nonasthma_nested) 
w_asthma <- 1/(p_asthma_nested/p_asthma_overall) 
w_nonasthma <- 1/(p_nonasthma_nested/p_nonasthma_overall) 
pheno$asthma_weights [pheno$LHSSCREEN_12142015 == 0] <- w_nonasthma 
pheno$asthma_weights [pheno$LHSSCREEN_12142015 == 1] <- w_asthma 
pheno_pc$asthma_weights [pheno_pc$LHSSCREEN_12142015 == 0] <- w_nonasthma 
pheno_pc$asthma_weights [pheno_pc$LHSSCREEN_12142015 == 1] <- w_asthma 
pheno_former$asthma_weights [pheno_former$LHSSCREEN_12142015 == 0] <- w_nonasthma 
pheno_former$asthma_weights [pheno_former$LHSSCREEN_12142015 == 1] <- w_asthma 
pheno_formerpc$asthma_weights [pheno_formerpc$LHSSCREEN_12142015 == 0] <- w_nonasthma 
pheno_formerpc$asthma_weights [pheno_formerpc$LHSSCREEN_12142015 == 1] <- w_asthma 


###########
## Model 3
for (i in 3:6) {
  f.RLM.Adj.Robust.wgt <- function(meth, dat1,j) {
    dat <- svydesign(id = ~Sample_ID, weights = ~asthma_weights, data = data.frame(y = meth, dat1))
    mdl.formula <- formula(paste("y ~ former_10yrs_split + ANYETS + LHSSCREEN_12142015 + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
    mdl <- svyglm(formula = mdl.formula, design = dat)
    cf <- summary(mdl)
    N <- sum(!is.na(meth))
    n.exp <- sum(!is.na(dat$variables$y) & dat$variables$former_10yrs_split == j)
    n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$former_10yrs_split == 0)
    coef = cf$coefficients[j,"Estimate"]
    se = cf$coefficients[j,"Std. Error"]
    ts = cf$coefficients[j,"t value"]
    p.value = cf$coefficients[j,"Pr(>|t|)"]
    out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
    names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
  }
  
  results <- parApply(cl, beta_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc,i)
  results <- as.data.frame(t(results))
  Nbeta <- sum(!is.na(results$P))
  write.csv(results, paste0("output/ten_year/model3_weighted.",i,".csv"))
}

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
for (i in 3:6) {
  f.RLM.Adj.Robust.wgt <- function(meth, dat1,j) {
    dat <- svydesign(id = ~Sample_ID, weights = ~asthma_weights, data = data.frame(y = meth, dat1))
    mdl.formula <- formula(paste("y ~ former_10yrs_split + ANYETS + European + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
    mdl <- svyglm(formula = mdl.formula, design = dat)
    cf <- summary(mdl)
    N <- sum(!is.na(meth))
    n.exp <- sum(!is.na(dat$variables$y) & dat$variables$former_10yrs_split == j)
    n.unexp <- sum(!is.na(dat$variables$y) & dat$variables$former_10yrs_split == 0)
    coef = cf$coefficients[j,"Estimate"]
    se = cf$coefficients[j,"Std. Error"]
    ts = cf$coefficients[j,"t value"]
    p.value = cf$coefficients[j,"Pr(>|t|)"]
    out = c(coef, se, ts, p.value, N, n.exp, n.unexp) 
    names(out) = c("Beta","SE", "tstat", "P", "N", "N.exp", "N.unexp"); return(out) 
  }
  
  results <- parApply(cl, beta2, 2, f.RLM.Adj.Robust.wgt, pheno,i)
  results <- as.data.frame(t(results))
  Nbeta <- sum(!is.na(results$P))
  write.csv(results, paste0("output/ten_year/model12_weighted.",i,".csv"))
}


############################################################################################################
##Formers only with cont years since quit
system('mkdir -p output/continuous')
f.RLM.Adj.Robust.wgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_TimeSinceSmked_v1 + ANYETS + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  coef = cf$coefficients[2,"Estimate"]
  se = cf$coefficients[2,"Std. Error"]
  ts = cf$coefficients[2,"t value"]
  p.value = cf$coefficients[2,"Pr(>|t|)"]
  out = c(coef, se, ts, p.value, N) 
  names(out) = c("Beta","SE", "tstat", "P", "N"); return(out) 
}
results <- parApply(cl, beta_formerpc, 2, f.RLM.Adj.Robust.wgt, pheno_formerpc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/continuous/former_cont_EU_unweighted.csv")


f.RLM.Adj.Robust.wgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, weights = ~asthma_weights, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_TimeSinceSmked_v1 + ANYETS + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  coef = cf$coefficients[2,"Estimate"]
  se = cf$coefficients[2,"Std. Error"]
  ts = cf$coefficients[2,"t value"]
  p.value = cf$coefficients[2,"Pr(>|t|)"]
  out = c(coef, se, ts, p.value, N) 
  names(out) = c("Beta","SE", "tstat", "P", "N"); return(out) 
}

results <- parApply(cl, beta_formerpc, 2, f.RLM.Adj.Robust.wgt, pheno_formerpc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/continuous/former_cont_EU_weighted.csv")




f.RLM.Adj.Robust.wgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_TimeSinceSmked_v1 + ANYETS + European + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  coef = cf$coefficients[2,"Estimate"]
  se = cf$coefficients[2,"Std. Error"]
  ts = cf$coefficients[2,"t value"]
  p.value = cf$coefficients[2,"Pr(>|t|)"]
  out = c(coef, se, ts, p.value, N) 
  names(out) = c("Beta","SE", "tstat", "P", "N"); return(out) 
}
results <- parApply(cl, beta_former, 2, f.RLM.Adj.Robust.wgt, pheno_former)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/continuous/former_cont_TA_unweighted.csv")


f.RLM.Adj.Robust.wgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, weights = ~asthma_weights, data = data.frame(y = meth, dat1))
  mdl.formula <- formula(paste("y ~ LH_TimeSinceSmked_v1 + ANYETS + European + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
                                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  cf <- summary(mdl)
  N <- sum(!is.na(meth))
  coef = cf$coefficients[2,"Estimate"]
  se = cf$coefficients[2,"Std. Error"]
  ts = cf$coefficients[2,"t value"]
  p.value = cf$coefficients[2,"Pr(>|t|)"]
  out = c(coef, se, ts, p.value, N) 
  names(out) = c("Beta","SE", "tstat", "P", "N"); return(out) 
}

results <- parApply(cl, beta_former, 2, f.RLM.Adj.Robust.wgt, pheno_former)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/continuous/former_cont_TA_weighted.csv")







################################################################################################################

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
##unweighted
#uws <- fread('unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.csv')
#uwf <- fread('unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.csv')
fin<- c('model1_unweighted.3.csv','model1_unweighted.4.csv','model1_unweighted.5.csv','model1_unweighted.6.csv',
        'model10_unweighted.3.csv','model10_unweighted.4.csv','model10_unweighted.5.csv','model10_unweighted.6.csv',
        'model3_weighted.3.csv','model3_weighted.4.csv','model3_weighted.5.csv','model3_weighted.6.csv',
        'model12_weighted.3.csv','model12_weighted.4.csv','model12_weighted.5.csv','model12_weighted.6.csv')
bonf<- 0.05/230329 
tsv_path <- "references/MSA.hg38.manifest.gencode.v41.tsv.gz"
manifest = fread(tsv_path)
for(f in fin){
  uws <- fread(paste0('output/five_year/',f))
  siguws <- uws %>% dplyr::filter(uws$P<=bonf)
  siguw_anno <- left_join(siguws,manifest,by=c("V1"="probeID"))
  f1 <- substring(f,1,nchar(f)-4)
  write.table(siguw_anno,file=paste0("output/five_year/signifcant_",f1,".tsv"),sep='\t',quote = F,na = "",row.names = F)
  uws <- fread(paste0('output/ten_year/',f))
  siguws <- uws %>% dplyr::filter(uws$P<=bonf)
  siguw_anno <- left_join(siguws,manifest,by=c("V1"="probeID"))
  f2 <- substring(f,1,nchar(f)-4)
  write.table(siguw_anno,file=paste0("output/ten_year/signifcant_",f2,".tsv"),sep='\t',quote = F,na = "",row.names = F)
}     
fin=c('former_cont_EU_unweighted.csv','former_cont_EU_weighted.csv',
      'former_cont_TA_unweighted.csv','former_cont_TA_weighted.csv')
for(f in fin){
  uws <- fread(paste0('output/continuous/',f))
  siguws <- uws %>% dplyr::filter(uws$P<=bonf)
  siguw_anno <- left_join(siguws,manifest,by=c("V1"="probeID"))
  f1 <- substring(f,1,nchar(f)-4)
  write.table(siguw_anno,file=paste0("output/continuous/signifcant_",f1,".tsv"),sep='\t',quote = F,na = "",row.names = F)
}     
  ####################################################
m1 <- fread('output/model_output/svyglm_sginificant_cpgs/unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
m1le <- fread('signifcant_model1_unweighted.3.tsv') %>% select(V1)
colnames(m1le)<-'Probe_ID'
m1gt <- fread('signifcant_model1_unweighted.4.tsv') %>% select(V1)
colnames(m1gt)<-'Probe_ID'
m1gtnole<- m1 %>% select(Probe_ID,CpG_chrm,CpG_beg,probe_strand,genesUniq) %>% inner_join(.,m1gt) 


m10 <- fread('output/model_output/svyglm_sginificant_cpgs/unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')
m10le <- fread('signifcant_model10_unweighted.3.tsv') %>% select(V1)
colnames(m10le)<-'Probe_ID'
m10gt <- fread('signifcant_model10_unweighted.4.tsv') %>% select(V1)
colnames(m10gt)<-'Probe_ID'
m10gtnole<- m10 %>% select(Probe_ID,CpG_chrm,CpG_beg,probe_strand,genesUniq) %>% inner_join(.,m10gt)


m3 <- fread('output/model_output/svyglm_sginificant_cpgs/weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
m3le <- fread('signifcant_model3_weighted.3.tsv') %>% select(V1)
colnames(m3le)<-'Probe_ID'
m3gt <- fread('signifcant_model3_weighted.4.tsv') %>% select(V1)
colnames(m3gt)<-'Probe_ID'
m3gtnole<- m3 %>% select(Probe_ID,CpG_chrm,CpG_beg,probe_strand,genesUniq) %>% inner_join(.,m3gt) 


m12 <- fread('output/model_output/svyglm_sginificant_cpgs/weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')
m12le <- fread('signifcant_model12_weighted.3.tsv') %>% select(V1)
colnames(m12le)<-'Probe_ID'
m12gt <- fread('signifcant_model12_weighted.4.tsv') %>% select(V1)
colnames(m12gt)<-'Probe_ID'
m12gtnole<- m12 %>% select(Probe_ID,CpG_chrm,CpG_beg,probe_strand,genesUniq) %>% inner_join(.,m12gt)


write.table(m1gtnole,file='output/former_strat/model1_overlap.tsv',sep='\t',quote=F,row.names = F)
write.table(m10gtnole,file='output/former_strat/model10_overlap.tsv',sep='\t',quote=F,row.names = F)
write.table(m3gtnole,file='output/former_strat/model3_overlap.tsv',sep='\t',quote=F,row.names = F)
write.table(m12gtnole,file='output/former_strat/model12_overlap.tsv',sep='\t',quote=F,row.names = F)

library(openxlsx2)
fy1 <- wb_workbook()
fin = c("output/five_year/signifcant_model10_unweighted.3.tsv", "output/five_year/signifcant_model10_unweighted.4.tsv", 
        "output/five_year/signifcant_model10_unweighted.5.tsv", "output/five_year/signifcant_model10_unweighted.6.tsv", 
        "output/five_year/signifcant_model12_weighted.3.tsv", "output/five_year/signifcant_model12_weighted.4.tsv", 
        "output/five_year/signifcant_model12_weighted.5.tsv", "output/five_year/signifcant_model12_weighted.6.tsv", 
        "output/five_year/signifcant_model1_unweighted.3.tsv", "output/five_year/signifcant_model1_unweighted.4.tsv", 
        "output/five_year/signifcant_model1_unweighted.5.tsv", "output/five_year/signifcant_model1_unweighted.6.tsv", 
        "output/five_year/signifcant_model3_weighted.3.tsv", "output/five_year/signifcant_model3_weighted.4.tsv", 
        "output/five_year/signifcant_model3_weighted.5.tsv", "output/five_year/signifcant_model3_weighted.6.tsv"
)
for(f in fin){
  uws <- fread(f)
  f1 = substring(f,1,nchar(f)-4) %>% substring(.,23)
  fy1 <- wb_add_worksheet(fy1, sheet = f1)
  fy1 <- wb_add_data(fy1, f1, uws)
} 
fy1$save("output/five_year_split_05.22.2025.xlsx")

fin = c("output/ten_year/signifcant_model10_unweighted.3.tsv", "output/ten_year/signifcant_model10_unweighted.4.tsv", 
  "output/ten_year/signifcant_model10_unweighted.5.tsv", "output/ten_year/signifcant_model10_unweighted.6.tsv", 
  "output/ten_year/signifcant_model12_weighted.3.tsv", "output/ten_year/signifcant_model12_weighted.4.tsv", 
  "output/ten_year/signifcant_model12_weighted.5.tsv", "output/ten_year/signifcant_model12_weighted.6.tsv", 
  "output/ten_year/signifcant_model1_unweighted.3.tsv", "output/ten_year/signifcant_model1_unweighted.4.tsv", 
  "output/ten_year/signifcant_model1_unweighted.5.tsv", "output/ten_year/signifcant_model1_unweighted.6.tsv", 
  "output/ten_year/signifcant_model3_weighted.3.tsv", "output/ten_year/signifcant_model3_weighted.4.tsv", 
  "output/ten_year/signifcant_model3_weighted.5.tsv", "output/ten_year/signifcant_model3_weighted.6.tsv"
)
ty1 <- wb_workbook()
for(f in fin){
  uws <- fread(f)
  f1 = substring(f,1,nchar(f)-4) %>% substring(.,22)
  ty1 <- wb_add_worksheet(ty1, sheet = f1)
  ty1 <- wb_add_data(ty1, f1, uws)
} 
ty1$save("output/ten_year_split_05.22.2025.xlsx")

fin = c("output/continuous/signifcant_former_cont_EU_unweighted.tsv", "output/output/continuous/signifcant_former_cont_EU_weighted.tsv", 
        "output/continuous/signifcant_former_cont_TA_unweighted.tsv", "output/output/continuous/signifcant_former_cont_TA_weighted.tsv"
)
cont <- wb_workbook()
for(f in fin){
  uws <- fread(f)
  f1 = substring(f,1,nchar(f)-4) %>% substring(.,23)
  cont <- wb_add_worksheet(cont, sheet = f1)
  cont <- wb_add_data(cont, f1, uws)
} 
cont$save("output/cont_yrs_since_smkd_05.22.2025.xlsx")





###############################################################################################################
## Plot AHRR and ECEL1P1
## We need to get the AHRR and ECEL1P1 probes that were significant in former smokers in the original models
## Subset beta df to former smokers and add new groups
## First we can generate boxplots by groups to see if there is a trend
## Next we need to generate CB's plot with x-axis being years since 
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggpubr)
load('../mini_models_probe_filt.rdata')
beta2 = t(betas) # methylation data
pheno = sampsheetGF

fm1 <- fread('output/model_output/svyglm_sginificant_cpgs/unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv') %>%
            filter(stringr::str_detect(genesUniq,'ECEL1P1|AHRR|PDCD6-AHRR')) %>%
            dplyr::select(Probe_ID,genesUniq)
fm3 <- fread('output/model_output/svyglm_sginificant_cpgs/unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1|AHRR|PDCD6-AHRR')) %>%
  dplyr::select(Probe_ID,genesUniq)
fm10 <- fread('output/model_output/svyglm_sginificant_cpgs/weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1|AHRR|PDCD6-AHRR')) %>%
  dplyr::select(Probe_ID,genesUniq)
fm12 <- fread('output/model_output/svyglm_sginificant_cpgs/weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1|AHRR|PDCD6-AHRR')) %>%
  dplyr::select(Probe_ID,genesUniq)
fm <- bind_rows(list(fm1,fm3,fm10,fm12)) %>% distinct()
probes <- fm$Probe_ID


## plot the probes by split
epic <- fread('references/EPIC.hg38.manifest.gencode.v36.tsv.gz') %>% mutate(epicT =T) %>% dplyr::select(probeID,epicT)
fourfity<- fread('references/HM450.hg38.manifest.tsv.gz')%>% mutate(ffiftyT =T) %>% dplyr::select(Probe_ID,ffiftyT)
##reduce probes to cpgs
probies <- manifest %>% filter(probeID %in% probes) %>% dplyr::select(probeID,genesUniq) %>%
  mutate(cpg = gsub(x = probeID, '\\_.*','')) %>% 
  left_join(epic,by=c('cpg'='probeID')) %>% left_join(fourfity,by=c('cpg'='Probe_ID')) %>% replace_na(list(epicT=FALSE,ffiftyT=FALSE)) %>%
  mutate(msa_only=case_when(epicT==TRUE | ffiftyT==TRUE ~ FALSE,.default=TRUE))
pdf("output/ahrr_ecel1p1_former_10yr_split.pdf")
for(p in probes){
  ahrr <- data.frame(beta2[,p])
  matching_rows <- grep(p, rownames(t(beta2)), value = TRUE)
  cpg_data <- t(beta2)[matching_rows, , drop = FALSE]
  cpg_long <- pivot_longer(as.data.frame(cpg_data),
                           cols = everything(),
                           names_to = "Sample_ID",
                           values_to = "Beta_Value")
  unique_samples <- unique(cpg_long$Sample_ID)
  ahrr<-left_join(cpg_long,pheno,by="Sample_ID") 

  ahrr = filter(ahrr,former_10yrs_split!="Missing")
  gene <- fm %>% dplyr::select(Probe_ID,genesUniq) %>% filter(Probe_ID == p) %>%
    dplyr::select(genesUniq) 
  print(ggboxplot(ahrr, x = "former_10yrs_split", y = "Beta_Value",
                  color = "former_10yrs_split", palette ="jama",order=c('Current Smoker','<=10yrs','11-20yrs','21-30yrs','>30yrs','Never Smoker'),
                  add = "jitter",xlab="Smoking Status", ylab = "Beta", title=paste0("Probe: ",p," Gene: ",gene,"\n by Years Since Quitting") )+
          theme(plot.title = element_text(hjust = 0.5)))
  
}
dev.off()


pdf("output/ahrr_ecel1p1_former_5yr_split.pdf")
for(p in probes){
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
  gene <- fm %>% dplyr::select(Probe_ID,genesUniq) %>% filter(Probe_ID == p) %>%
    dplyr::select(genesUniq) 
  print(ggboxplot(ahrr, x = "former_5yrs_split", y = "Beta_Value",
                  color = "former_5yrs_split", palette ="jama",order=c('Current Smoker','<=5yrs','6-10yrs','11-15yrs','>15yrs','Never Smoker'),
                  add = "jitter",xlab="Smoking Status", ylab = "Beta", title=paste0("Probe: ",p," Gene: ",gene,"\n by Years Since Quitting") )+
          theme(plot.title = element_text(hjust = 0.5)))
  
}
dev.off()




#####################################################################################################################################
## Grab T-statistics 
m<-c("1_un","3_","10_un","12_")
pdf("output/ahrr_ecel1p1_former_5yr_tstat.pdf",width=10)
for (i in m) {
    fm2 <- fread(paste0('five_year/model',i,'weighted.2.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = 'Current')
    fm3 <- fread(paste0('five_year/model',i,'weighted.3.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '<=5yr')
    fm4 <- fread(paste0('five_year/model',i,'weighted.4.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '6-10yr')
    fm5 <- fread(paste0('five_year/model',i,'weighted.5.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '11-15yr')
    fm6 <- fread(paste0('five_year/model',i,'weighted.6.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '>15yr')
    fm <- bind_rows(list(fm2,fm3,fm4,fm5,fm6))
    fm$smoke <- factor(fm$smoke,levels = c('Current','<=5yr','6-10yr','11-15yr','>15yr'))
    fm <- fm %>% left_join(.,probies,by=c('V1'='probeID'))
    i= substring(i,1,2)
    ahrr <- data.frame(fm)
    print(ggplot(data=ahrr, aes(x=smoke, y=tstat,group=V1)) + geom_line(aes(color=genesUniq,linetype=msa_only)) +
      ggtitle(paste0('Model ',i,'long-term')) + xlab('Years after quitting') + ylab('T-statistics') + theme_classic())
    
}    
dev.off()
pdf("output/ahrr_ecel1p1_former_10yr_tstat.pdf",width=10)
for (i in m) {
    fm2 <- fread(paste0('ten_year/model',i,'weighted.2.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = 'Current')
    fm3 <- fread(paste0('ten_year/model',i,'weighted.3.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '<=10yr')
    fm4 <- fread(paste0('ten_year/model',i,'weighted.4.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '11-20yr')
    fm5 <- fread(paste0('ten_year/model',i,'weighted.5.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '21-30yr')
    fm6 <- fread(paste0('ten_year/model',i,'weighted.6.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '>30yr')
    fm <- bind_rows(list(fm2,fm3,fm4,fm5,fm6))
    fm$smoke <- factor(fm$smoke,levels = c('Current','<=10yr','11-20yr','21-30yr','>30yr'))
    fm <- fm %>% left_join(.,probies,by=c('V1'='probeID'))
    i= substring(i,1,2)
    ahrr <- data.frame(fm)
    print(ggplot(data=ahrr, aes(x=smoke, y=tstat,group=V1)) + geom_line(aes(color=genesUniq,linetype=msa_only)) +
      ggtitle(paste0('Model ',i,'long-term')) + xlab('Years after quitting') + ylab('T-statistics') + theme_classic() )
    
}
dev.off()




######################################################################################################################
## Overlaps | Venn Diagrams

library(ggvenn)
m<-c(1,3,10,12)
m1 <- fread('output/model_output/svyglm_sginificant_cpgs/unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
m10 <- fread('output/model_output/svyglm_sginificant_cpgs/unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')
m3 <- fread('output/model_output/svyglm_sginificant_cpgs/weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
m12 <- fread('output/model_output/svyglm_sginificant_cpgs/weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')

s1 <- fread('output/significant_cpgs/svyglm/smokers/imputed/unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
s10 <- fread('output/significant_cpgs/svyglm/smokers/imputed/unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')
s3 <- fread('output/significant_cpgs/svyglm/smokers/imputed/weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
s12 <- fread('output/significant_cpgs/svyglm/smokers/imputed/weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')
mm=list(m1,'',m3,'','','','','','',m10,'',m12)
ss<-list(s1,'',s3,'','','','','','',s10,'',s12)
tt<-c('','','<= 5 years','6-10 years','11-15 years','> 15 years')
pdf('output/venn_five_year.pdf')
for (i in m) {
  fm2 <- fread(paste0('five_year/model',i,'.2.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Current')
  fm3 <- fread(paste0('five_year/model',i,'.3.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '<=5yr')
  fm4 <- fread(paste0('five_year/model',i,'.4.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '6-10yr')
  fm5 <- fread(paste0('five_year/model',i,'.5.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '11-15yr')
  fm6 <- fread(paste0('five_year/model',i,'.6.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '>15yr')
  fog <- mm[[i]]
  sog<- ss[[i]]
  x = list(Original_Former_Smoker_Model=fog$Probe_ID,
           Original_Smoker_Model=sog$Probe_ID,
           Less_Five_yr = fm3$V1,
           Six_Ten_yr=fm4$V1,
           Eleven_Fifteen_yr=fm5$V1,
           Greater_Fifteen_yr=fm6$V1)
  for (j in 3:6) {
      y = x[c(1,2,j)]
      print(ggvenn(
        y, 
        fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF","#0003B3DD"),
        stroke_size = 0.5, set_name_size = 3
      ) + ggtitle(paste0('Model ',i,' ',tt[j])))
  }
}
dev.off()



tt<-c('','','<= 10 years','11-20 years','21-30 years','> 30 years')
pdf('output/venn_ten_year.pdf')
for (i in m) {
  fm2 <- fread(paste0('ten_year/model',i,'.2.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Current')
  fm3 <- fread(paste0('ten_year/model',i,'.3.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '<=10yr')
  fm4 <- fread(paste0('ten_year/model',i,'.4.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '11-20yr')
  fm5 <- fread(paste0('ten_year/model',i,'.5.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '21-30yr')
  fm6 <- fread(paste0('ten_year/model',i,'.6.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '>30 yr')
  fog <- mm[[i]]
  sog<- ss[[i]]
  x = list(Original_Former_Smoker_Model=fog$Probe_ID,
           Original_Smoker_Model=sog$Probe_ID,
           Less_Ten_yr = fm3$V1,
           Eleven_Twenty_yr=fm4$V1,
           TwentyOne_Thirty_yr=fm5$V1,
           Greater_Thirty_yr=fm6$V1)
  for (j in 3:6) {
    y = x[c(1,2,j)]
    print(ggvenn(
      y, 
      fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF","#0003B3DD"),
      stroke_size = 0.5, set_name_size = 3
    ) + ggtitle(paste0('Model ',i,' ', tt[j])))
  }
}
dev.off()

tt<-c('','','Unweighted Euro','Weighted Euro','Unweighted Trans-ancestry','Weighted Trans-ancestry')
pdf('output/venn_continuous.pdf')
for (i in m) {
  #fm2 <- fread(paste0('ten_year/model',i,'.2.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Current')
  fm3 <- fread(paste0('continuous/former_cont_EU_unweighted.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Unweighted Euro')
  fm4 <- fread(paste0('continuous/former_cont_EU_weighted.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Weighted Euro')
  fm5 <- fread(paste0('continuous/former_cont_TA_unweighted.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Unweighted TA')
  fm6 <- fread(paste0('continuous/former_cont_TA_weighted.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Weighted TA')
  fog <- mm[[i]]
  sog<- ss[[i]]
  x = list(Original_Former_Smoker_Model=fog$Probe_ID,
           Original_Smoker_Model=sog$Probe_ID,
           Unweighted_Euro = fm3$V1,
           Weighted_Euro=fm4$V1,
           Unweighted_TA=fm5$V1,
           Weighted_TA=fm6$V1)
  for (j in 3:6) {
    y = x[c(1,2,j)]
    print(ggvenn(
      y, 
      fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF","#0003B3DD"),
      stroke_size = 0.5, set_name_size = 3
    ) + ggtitle(paste0('Model ',i,' ', tt[j])))
  }
}
dev.off()




############################################################################################
## Plot Betas vs continuous time since smoked


pdf("output/ahrr_ecel1p1_scatter_plot.pdf")
for(p in probes){
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
  gene <- manifest %>% dplyr::filter(probeID == p) %>% dplyr::select(genesUniq)
  print(ggplot(ahrr, aes(x = LH_TimeSinceSmked_v1.2, y = Beta_Value)) + geom_point(aes(color=smoke_status)) + 
          geom_smooth(method='loess', formula= y~x) +
          theme(plot.title = element_text(hjust = 0.5)) + theme_classic() + ggtitle(paste0('Probe: ',p,' Gene: ',gene)))
  
}
dev.off()
################
## Just former smokers
pdf("output/ahrr_ecel1p1_scatter_plot_former_only.pdf")
for(p in probes){
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
                                                         LH_smoke_status == 1 ~ 'Former Smoker')) %>% filter(smoke_status=='Former Smoker')
  gene <- manifest %>% dplyr::filter(probeID == p) %>% dplyr::select(genesUniq)
  print(ggplot(ahrr, aes(x = LH_TimeSinceSmked_v1.2, y = Beta_Value)) + geom_point(aes(color=smoke_status), show.legend = FALSE) + 
          geom_smooth(method='loess', formula= y~x, ) +
          theme(plot.title = element_text(hjust = 0.5)) + theme_classic() + ggtitle(paste0('Probe: ',p,' Gene: ',gene)))
  
}
dev.off()

####################
## With avg lines
pdf("output/ahrr_ecel1p1_scatter_plot_former_only.pdf")
for(p in probes){
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
  nev = ahrr %>% filter(smoke_status == 'Never Smoker')
  nevb = mean(nev$Beta_Value)
  ahrr = ahrr %>% filter(smoke_status=='Former Smoker')
  gene <- manifest %>% dplyr::filter(probeID == p) %>% dplyr::select(genesUniq)
  print(ggplot(ahrr, aes(x = LH_TimeSinceSmked_v1.2, y = Beta_Value)) + geom_point(aes(color=smoke_status), show.legend = FALSE) + 
          geom_smooth(method='loess', formula= y~x, ) + geom_hline(yintercept=currb,linetype=2,color = "#CD534CFF") + 
          geom_hline(yintercept=nevb,linetype=2,color = "#0073C2FF") + 
          theme(plot.title = element_text(hjust = 0.5)) + theme_classic() + ggtitle(paste0('Probe: ',p,' Gene: ',gene)))
  
}
dev.off()

pdf("output/ahrr_ecel1p1_scatter_plot_former_only_with4splines_wo_se.pdf")
for(p in probes){
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
  nev = ahrr %>% filter(smoke_status == 'Never Smoker')
  nevb = mean(nev$Beta_Value)
  ahrr = ahrr %>% filter(smoke_status=='Former Smoker')
  gene <- manifest %>% dplyr::filter(probeID == p) %>% dplyr::select(genesUniq)
  print(ggplot(ahrr, aes(x = LH_TimeSinceSmked_v1.2, y = Beta_Value)) + geom_point(aes(color=smoke_status), show.legend = FALSE) + 
          geom_smooth(method='lm', formula= y~bs(x,4),se=F ) + geom_hline(yintercept=currb,linetype=2,color = "#CD534CFF") + 
          geom_hline(yintercept=nevb,linetype=2,color = "#0073C2FF") + 
          theme(plot.title = element_text(hjust = 0.5)) + theme_classic() + ggtitle(paste0('Probe: ',p,' Gene: ',gene)))
  
}
dev.off()






ahrr_ecel1p1 <- data.frame(beta_probes) %>% tibble::rownames_to_column(.,'Sample_ID') %>% 
  left_join(.,pheno,by='Sample_ID') %>% dplyr::mutate(LH_TimeSinceSmked_v1.2 = case_when(LH_smoke_status == 0 ~ 67,
                                                                                         LH_smoke_status == 2 ~ 0,
                                                                                         .default = LH_TimeSinceSmked_v1),
                                                      smoke_status = case_when(LH_smoke_status == 0 ~ 'Never Smoker',
                                                                               LH_smoke_status == 2 ~ 'Current Smoker',
                                                                               LH_smoke_status == 1 ~ 'Former Smoker')) %>% 
  dplyr::select(Sample_ID,probes,LH_TimeSinceSmked_v1,LH_TimeSinceSmked_v1.2,LH_smoke_status,smoke_status,former_5yrs_split,former_10yrs_split)
save(ahrr_ecel1p1,file='ahrr_ecel1p1_betas.rdata')
beta_probes <- as.data.frame(beta2) %>% dplyr::select(probes)

for(meth in colnames(beta_probes)){
  d1 <- data.frame(y = beta_probes[,1], pheno)
  dat <- svydesign(id = ~Sample_ID, data = d1)
  mdl.formula <- formula(paste("y ~ LH_smoke_status  + Gender + LH_BMI + AgeAtHomeVisit + STATE + B + NK + CD4T + CD8T + Mono + Neutro  +
  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  mdl <- svyglm(formula = mdl.formula, design = dat)
  d1 <- d1 %>% dplyr::mutate(p = predict(mdl))
  d1 = d1 %>% dplyr::mutate(LH_TimeSinceSmked_v1.2 = case_when(LH_smoke_status == 0 ~ 75,
                                                                   LH_smoke_status == 2 ~ 0,
                                                                   .default = LH_TimeSinceSmked_v1),
                                smoke_status = case_when(LH_smoke_status == 0 ~ 'Never Smoker',
                                                         LH_smoke_status == 2 ~ 'Current Smoker',
                                                         LH_smoke_status == 1 ~ 'Former Smoker'))
  ggplot(d1, aes(x = LH_TimeSinceSmked_v1, y = p)) + geom_point(aes(color=smoke_status)) + geom_smooth(method='loess', formula= y~x)
  
}

########################################################################################################################
## Predicted Betas for former smokers with >0 yrs since smoking

preds <- fread('probe_predicted_beta.csv')
probes<-names(preds)[3:length(names(preds))]
pdf("output/ahrr_ecel1p1_predicted_scatter_plot.pdf",width=11)
for(p in probes){
  ##gsub to get correct names
  ##we need to grab current and never betas to get median and %
  cp <- substr(p, 1, nchar(p) - 5)
  ahrr_pred <- data.frame(preds[,c('Sample_ID',..p)])
  rownames(ahrr_pred) <- ahrr_pred[,1]
  ahrr_pred[,1] <- NULL
  matching_rows <- grep(p, rownames(t(preds)), value = TRUE)
  cpg_data <- t(ahrr_pred)[matching_rows, , drop = FALSE]
  cpg_long <- pivot_longer(as.data.frame(cpg_data),
                           cols = everything(),
                           names_to = "Sample_ID",
                           values_to = "Beta_Value")
  unique_samples <- unique(cpg_long$Sample_ID)
  ahrr_pred<-left_join(cpg_long,pheno,by="Sample_ID")
  
  ahrr_pred = filter(ahrr_pred,former_5yrs_split!="Missing")
  ahrr_pred = ahrr_pred %>% dplyr::mutate(LH_TimeSinceSmked_v1.2 = case_when(LH_smoke_status == 0 ~ 67,
                                                                   LH_smoke_status == 2 ~ 0,
                                                                   .default = LH_TimeSinceSmked_v1),
                                smoke_status = case_when(LH_smoke_status == 0 ~ 'Never Smoker',
                                                         LH_smoke_status == 2 ~ 'Current Smoker',
                                                         LH_smoke_status == 1 ~ 'Former Smoker')) 
  ahrr <- data.frame(beta2[,cp])
  matching_rows <- grep(cp, rownames(t(beta2)), value = TRUE)
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
  gene <- manifest %>% dplyr::filter(probeID == cp) %>% dplyr::select(genesUniq)
  print(ggplot(ahrr_pred, aes(x = LH_TimeSinceSmked_v1.2, y = Beta_Value)) + geom_point(aes(color=smoke_status))  + scale_color_manual(values = c('Former Smoker'='#DF8F44FF')) + labs(color='Smoke Status') + 
          geom_smooth(method='loess', formula= y~x,color = "#374E55FF") + geom_hline(aes(yintercept=currb,linetype="Current Smokers"),color = "#B24745FF") +  
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = curr5, ymax = curr95, fill = "#B24745FF", alpha = .1, color = NA) +
          geom_hline(aes(yintercept=nevb,linetype="Never Smokers"),color = "#00A1D5FF") + annotate("rect", xmin = -Inf, xmax = Inf, ymin = nev5, ymax = nev95, fill = "#00A1D5FF", alpha = .1, color = NA) +
          scale_linetype_manual(name = "Mean Beta", 
                                values = c(2, 2)) +
          theme(plot.title = element_text(hjust = 0.5)) + theme_classic() + ggtitle(paste0('Probe: ',cp,' Gene: ',gene)))
}
dev.off()

pdf("output/ahrr_ecel1p1_median_scatter_plot.pdf",width=11)
for(p in probes){
  ##gsub to get correct names
  ##we need to grab current and never betas to get median and %
  cp <- substr(p, 1, nchar(p) - 5)
  ahrr_pred <- data.frame(preds[,c('Sample_ID',..p)])
  rownames(ahrr_pred) <- ahrr_pred[,1]
  ahrr_pred[,1] <- NULL
  matching_rows <- grep(p, rownames(t(preds)), value = TRUE)
  cpg_data <- t(ahrr_pred)[matching_rows, , drop = FALSE]
  cpg_long <- pivot_longer(as.data.frame(cpg_data),
                           cols = everything(),
                           names_to = "Sample_ID",
                           values_to = "Beta_Value")
  unique_samples <- unique(cpg_long$Sample_ID)
  ahrr_pred<-left_join(cpg_long,pheno,by="Sample_ID")
  
  ahrr_pred = filter(ahrr_pred,former_5yrs_split!="Missing")
  ahrr_pred = ahrr_pred %>% dplyr::mutate(LH_TimeSinceSmked_v1.2 = case_when(LH_smoke_status == 0 ~ 67,
                                                                             LH_smoke_status == 2 ~ 0,
                                                                             .default = LH_TimeSinceSmked_v1),
                                          smoke_status = case_when(LH_smoke_status == 0 ~ 'Never Smoker',
                                                                   LH_smoke_status == 2 ~ 'Current Smoker',
                                                                   LH_smoke_status == 1 ~ 'Former Smoker')) 
  ahrr <- data.frame(beta2[,cp])
  matching_rows <- grep(cp, rownames(t(beta2)), value = TRUE)
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
  gene <- manifest %>% dplyr::filter(probeID == cp) %>% dplyr::select(genesUniq)
  print(ggplot(ahrrf, aes(x = LH_TimeSinceSmked_v1.2, y = Beta_Value)) + geom_point(aes(color=smoke_status))  + scale_color_manual(values = c('Former Smoker'='#DF8F44FF')) + labs(color='Smoke Status') + 
          geom_smooth(method='loess', formula= y~x,color = "#374E55FF") + geom_hline(aes(yintercept=currb,linetype="Current Smokers"),color = "#B24745FF") +  
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = curr5, ymax = curr95, fill = "#B24745FF", alpha = .1, color = NA) +
          geom_hline(aes(yintercept=nevb,linetype="Never Smokers"),color = "#00A1D5FF") + annotate("rect", xmin = -Inf, xmax = Inf, ymin = nev5, ymax = nev95, fill = "#00A1D5FF", alpha = .1, color = NA) +
          scale_linetype_manual(name = "Median Beta", 
                                values = c(2, 2)) +
          theme(plot.title = element_text(hjust = 0.5)) + theme_classic() + ggtitle(paste0('Probe: ',cp,' Gene: ',gene)))
}
dev.off()
