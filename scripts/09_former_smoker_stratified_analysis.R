library(dplyr)
load('751_clean_pc.rdata')
dir.create("output/former_stratified", recursive=TRUE, showWarnings=FALSE)
dir.create("output/five_year", recursive=TRUE, showWarnings=FALSE)
dir.create("output/ten_year", recursive=TRUE, showWarnings=FALSE)
dir.create("output/continuous", recursive=TRUE, showWarnings=FALSE)
#beta2 = t(clean_betas) # methylation data

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
##pc only
#pheno=pheno%>%mutate(asthma_weights=ifelse(LHSSCREEN_12142015=="Case",w_asthma,w_nonasthma))
#pheno = pheno %>% filter(Gender == inrdsx)
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
pheno$former_5yrs_split <- factor(pheno$former_5yrs_split, levels = c(0,1,2,3,4,5,99))
pheno <- pheno %>% mutate(former_10yrs_split = case_when(LH_smoke_status == 0 ~ factor(0),
                                                       LH_smoke_status == 2 ~ factor(1),
                                                       LH_smoke_status == 1 & LH_TimeSinceSmked_v1 <= 10 ~ factor(2),
                                                       LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 10 & LH_TimeSinceSmked_v1 <= 20 ~ factor(3),
                                                       LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 20 & LH_TimeSinceSmked_v1 <= 30 ~ factor(4),
                                                       LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 30 ~ factor(5),
                                                       LH_smoke_status == 99 ~ factor(99)))
pheno$former_10yrs_split <- factor(pheno$former_10yrs_split, levels = c(0,1,2,3,4,5,99))
##need to keep never smokers
betap <- t(betas_miss_masked)
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
betap_former <- betap[row.names(betap) %in% pheno_former$Sample_ID,]
betap_formerpc <- betap[row.names(betap) %in% pheno_formerpc$Sample_ID,]
#######################################################################################
##########
##Model 1 unweighted euro
## intercept 1
## smoke_2 2
## smoke_3 3
## smoke_4 4
## smoke_99 5

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

results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_stratified/model1_unweighted.3.csv")

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

results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_stratified/model1_unweighted.4.csv")
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

results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.unwgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_stratified/model10_unweighted.3.csv")
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

results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.unwgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_stratified/model10_unweighted.4.csv")


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

results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_stratified/model3_weighted.3.csv")

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

results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_stratified/model3_weighted.4.csv")

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

results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.wgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_stratified/model12_weighted.3.csv")

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

results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.wgt, pheno)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/former_stratified/model12_weighted.4.csv")


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
  
  results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc,i)
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
  
  results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.unwgt, pheno,i)
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
  
  results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc,i)
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
  
  results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.wgt, pheno,i)
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
  
  results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.unwgt, pheno_pc,i)
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
  
  results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.unwgt, pheno,i)
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
  
  results <- parApply(cl, betap_pc, 2, f.RLM.Adj.Robust.wgt, pheno_pc,i)
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
  
  results <- parApply(cl, betap, 2, f.RLM.Adj.Robust.wgt, pheno,i)
  results <- as.data.frame(t(results))
  Nbeta <- sum(!is.na(results$P))
  write.csv(results, paste0("output/ten_year/model12_weighted.",i,".csv"))
}


############################################################################################################
##Formers only with cont years since quit
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
results <- parApply(cl, betap_formerpc, 2, f.RLM.Adj.Robust.wgt, pheno_formerpc)
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

results <- parApply(cl, betap_formerpc, 2, f.RLM.Adj.Robust.wgt, pheno_formerpc)
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
results <- parApply(cl, betap_former, 2, f.RLM.Adj.Robust.wgt, pheno_former)
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

results <- parApply(cl, betap_former, 2, f.RLM.Adj.Robust.wgt, pheno_former)
results <- as.data.frame(t(results))
Nbeta <- sum(!is.na(results$P))
write.csv(results, "output/continuous/former_cont_TA_weighted.csv")







################################################################################################################
library(data.table)

bonf <- 0.05/230325
manifest <- fread("references/MSA.hg38.manifest.gencode.v41.tsv.gz")

analysis_dirs <- c(
  "output/former_stratified",
  "output/five_year",
  "output/ten_year",
  "output/continuous"
)

annotate_significant <- function(input_dir) {
  files <- list.files(input_dir, pattern="\\.csv$", full.names=TRUE)
  
  for (f in files) {
    results <- fread(f)
    
    if ("V1" %in% names(results)) setnames(results, "V1", "Probe_ID")
    if (!"Probe_ID" %in% names(results)) stop("Probe ID column not found in: ", f)
    if (!"P" %in% names(results)) stop("P-value column not found in: ", f)
    
    significant <- results[P <= bonf]
    significant <- merge(significant, manifest, by.x="Probe_ID", by.y="probeID", all.x=TRUE, sort=FALSE)
    
    output_file <- file.path(input_dir, paste0("significant_", tools::file_path_sans_ext(basename(f)), ".tsv"))
    fwrite(significant, output_file, sep="\t", na="")
    
    message(basename(f), ": ", nrow(significant), " significant probes")
  }
}

for (d in analysis_dirs) annotate_significant(d)


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
load('751_clean_pc.rdata')
beta2 = t(betas) # methylation data
pheno = sampsheetGF

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
fm <- bind_rows(list(fm1,fm3,fm10,fm12)) %>% distinct()
probes <- fm$Probe_ID


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

pheno <- pheno %>% mutate(former_5yrs_split = case_when(LH_smoke_status == 0 ~ factor('Never Smoker'),
                                                        LH_smoke_status == 2 ~ factor('Current Smoker'),
                                                        LH_smoke_status == 1 & LH_TimeSinceSmked_v1 <= 5 ~ factor('<=5yrs'),
                                                        LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 5 & LH_TimeSinceSmked_v1 <= 10 ~ factor('6-10yrs'),
                                                        LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 10 & LH_TimeSinceSmked_v1 <= 15 ~ factor('11-15yrs'),
                                                        LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 15 ~ factor('>15yrs'),
                                                        LH_smoke_status == 99 ~ factor('Missing')))
pheno <- pheno %>% mutate(former_10yrs_split = case_when(LH_smoke_status == 0 ~ factor('Never Smoker'),
                                                         LH_smoke_status == 2 ~ factor('Current Smoker'),
                                                         LH_smoke_status == 1 & LH_TimeSinceSmked_v1 <= 10 ~ factor('<=10yrs'),
                                                         LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 10 & LH_TimeSinceSmked_v1 <= 20 ~ factor('11-20yrs'),
                                                         LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 20 & LH_TimeSinceSmked_v1 <= 30 ~ factor('21-30yrs'),
                                                         LH_smoke_status == 1 & LH_TimeSinceSmked_v1 > 30 ~ factor('>30yrs'),
                                                         LH_smoke_status == 99 ~ factor('Missing')))
epic <- fread('references/EPIC.hg38.manifest.gencode.v36.tsv.gz') %>% mutate(epicT =T) %>% dplyr::select(probeID,epicT)
fourfity<- fread('references/HM450.hg38.manifest.tsv.gz')%>% mutate(ffiftyT =T) %>% dplyr::select(Probe_ID,ffiftyT)
##reduce probes to cpgs
probies <- manifest %>% filter(probeID %in% probes) %>% dplyr::select(probeID,genesUniq) %>%
  mutate(cpg = gsub(x = probeID, '\\_.*','')) %>% 
  left_join(epic,by=c('cpg'='probeID')) %>% left_join(fourfity,by=c('cpg'='Probe_ID')) %>% replace_na(list(epicT=FALSE,ffiftyT=FALSE)) %>%
  mutate(msa_only=case_when(epicT==TRUE | ffiftyT==TRUE ~ FALSE,.default=TRUE))
pdf("output/former_stratified/ahrr_ecel1p1_former_10yr_split.pdf")
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


pdf("output/former_stratified/ahrr_ecel1p1_former_5yr_split.pdf")
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
pdf("output/former_stratified/ahrr_ecel1p1_former_5yr_tstat.pdf",width=10)
for (i in m) {
    fm2 <- fread(paste0('output/five_year/model',i,'weighted.2.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = 'Current')
    fm3 <- fread(paste0('output/five_year/model',i,'weighted.3.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '<=5yr')
    fm4 <- fread(paste0('output/five_year/model',i,'weighted.4.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '6-10yr')
    fm5 <- fread(paste0('output/five_year/model',i,'weighted.5.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '11-15yr')
    fm6 <- fread(paste0('output/five_year/model',i,'weighted.6.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '>15yr')
    fm <- bind_rows(list(fm2,fm3,fm4,fm5,fm6))
    fm$smoke <- factor(fm$smoke,levels = c('Current','<=5yr','6-10yr','11-15yr','>15yr'))
    fm <- fm %>% left_join(.,probies,by=c('V1'='probeID'))
    i= substring(i,1,2)
    ahrr <- data.frame(fm)
    print(ggplot(data=ahrr, aes(x=smoke, y=tstat,group=V1)) + geom_line(aes(color=genesUniq,linetype=msa_only)) +
      ggtitle(paste0('Model ',i,'long-term')) + xlab('Years after quitting') + ylab('T-statistics') + theme_classic())
    
}    
dev.off()
pdf("output/former_stratified/ahrr_ecel1p1_former_10yr_tstat.pdf",width=10)
for (i in m) {
    fm2 <- fread(paste0('output/ten_year/model',i,'weighted.2.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = 'Current')
    fm3 <- fread(paste0('output/ten_year/model',i,'weighted.3.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '<=10yr')
    fm4 <- fread(paste0('output/ten_year/model',i,'weighted.4.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '11-20yr')
    fm5 <- fread(paste0('output/ten_year/model',i,'weighted.5.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '21-30yr')
    fm6 <- fread(paste0('output/ten_year/model',i,'weighted.6.csv')) %>% dplyr::filter(V1 %in% probes) %>% dplyr::mutate(smoke = '>30yr')
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
m1 <- fread('output/model_output/svyglm_significant_cpgs/unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
m10 <- fread('output/model_output/svyglm_significant_cpgs/unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')
m3 <- fread('output/model_output/svyglm_significant_cpgs/weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
m12 <- fread('output/model_output/svyglm_significant_cpgs/weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')

s1 <- fread('output/model_output/svyglm_significant_cpgs/unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
s10 <- fread('output/model_output/svyglm_significant_cpgs/unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')
s3 <- fread('output/model_output/svyglm_significant_cpgs/weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
s12 <- fread('output/model_output/svyglm_significant_cpgs/weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')
mm=list(m1,'',m3,'','','','','','',m10,'',m12)
ss<-list(s1,'',s3,'','','','','','',s10,'',s12)
tt<-c('','','<= 5 years','6-10 years','11-15 years','> 15 years')
pdf('output/five_year/venn_five_year.pdf')
for (i in m) {
  fm2 <- fread(paste0('output/five_year/model',i,'.2.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Current')
  fm3 <- fread(paste0('output/five_year/model',i,'.3.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '<=5yr')
  fm4 <- fread(paste0('output/five_year/model',i,'.4.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '6-10yr')
  fm5 <- fread(paste0('output/five_year/model',i,'.5.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '11-15yr')
  fm6 <- fread(paste0('output/five_year/model',i,'.6.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '>15yr')
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
pdf('output/ten_year/venn_ten_year.pdf')
for (i in m) {
  fm2 <- fread(paste0('output/ten_year/model',i,'.2.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Current')
  fm3 <- fread(paste0('output/ten_year/model',i,'.3.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '<=10yr')
  fm4 <- fread(paste0('output/ten_year/model',i,'.4.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '11-20yr')
  fm5 <- fread(paste0('output/ten_year/model',i,'.5.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '21-30yr')
  fm6 <- fread(paste0('output/ten_year/model',i,'.6.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = '>30 yr')
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
pdf('output/continuous/venn_continuous.pdf')
for (i in m) {
  #fm2 <- fread(paste0('output/ten_year/model',i,'.2.csv')) %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Current')
  fm3 <- fread('output/continuous/former_cont_EU_unweighted.csv') %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Unweighted Euro')
  fm4 <- fread('output/continuous/former_cont_EU_weighted.csv') %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Weighted Euro')
  fm5 <- fread('output/continuous/former_cont_TA_unweighted.csv') %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Unweighted TA')
  fm6 <- fread('output/continuous/former_cont_TA_weighted.csv') %>% dplyr::filter(P<=bonf) %>% dplyr::mutate(smoke = 'Weighted TA')
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


pdf("output/former_stratified/ahrr_ecel1p1_scatter_plot.pdf")
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
pdf("output/former_stratified/ahrr_ecel1p1_scatter_plot_former_only.pdf")
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
pdf("output/former_stratified/ahrr_ecel1p1_scatter_plot_former_only.pdf")
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

pdf("output/former_stratified/ahrr_ecel1p1_scatter_plot_former_only_with4splines_wo_se.pdf")
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
save(ahrr_ecel1p1,file='output/former_stratified/ahrr_ecel1p1_betas.rdata')
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

preds <- fread('output/former_stratified/probe_predicted_beta.csv')
probes<-names(preds)[3:length(names(preds))]
pdf("output/former_stratified/ahrr_ecel1p1_predicted_scatter_plot.pdf",width=11)
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

pdf("output/former_stratified/ahrr_ecel1p1_median_scatter_plot.pdf",width=11)
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
