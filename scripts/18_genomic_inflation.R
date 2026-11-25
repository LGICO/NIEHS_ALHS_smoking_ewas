library(bacon)
library(dplyr)
library(data.table)
bonf = 0.05/230325
dir.create("output/genomic_inflation", recursive=TRUE, showWarnings=FALSE)
dir.create("output/genomic_inflation/significant", recursive=TRUE, showWarnings=FALSE)
dir.create("output/genomic_inflation/full", recursive=TRUE, showWarnings=FALSE)
wes <- fread('output/svyglm_anno/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv') %>% 
              mutate(Z = Beta/SE,
                     chisq = Z^2)
wef <- fread('output/svyglm_anno/unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv') %>% 
               mutate(Z = Beta/SE,
                     chisq = Z^2)

wts <- fread('output/svyglm_anno/unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv') %>% 
              mutate(Z = Beta/SE,
                     chisq = Z^2)
wtsu <- fread('output/svyglm_anno/unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv') %>% 
  mutate(Z = Beta/SE,
         chisq = Z^2)
wtf <- fread('output/svyglm_anno/unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv') %>% 
              mutate(Z = Beta/SE,
                     chisq = Z^2)
lambda_gc <- function(beta, se) {
  z <- beta / se
  z <- z[is.finite(z)]
  
  median(z^2) / qchisq(0.5, df = 1)
}
lambda_gc_p <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  chisq <- qchisq(p, df = 1, lower.tail = FALSE)
  
  median(chisq) / qchisq(0.5, df = 1)
}
lambda_gc_z <- function(z) {
  z <- z[is.finite(z)]
  median(z^2) / qchisq(0.5, df = 1)
}

wtsu_bc <-bacon(teststatistics = wtsu$Z,effectsizes = wtsu$Beta, standarderrors = wtsu$SE)
inflation(wtsu_bc)

wes_bc <- bacon(teststatistics = wes$Z,effectsizes = wes$Beta, standarderrors = wes$SE)
inflation(wes_bc)
#sigma.0 
#1.173726 
wes$Beta_bacon <- es(wes_bc, corrected = TRUE)
wes$SE_bacon   <- se(wes_bc, corrected = TRUE)
wes$P_bacon    <- pval(wes_bc, corrected = TRUE)
c(
  uncorrected = lambda_gc_z(tstat(wes_bc, corrected = FALSE)),
  bacon_corrected = lambda_gc_z(tstat(wes_bc, corrected = TRUE))
)
wes_sig <- wes %>%
  filter(P_bacon <= bonf)
lambda_gc(wes$Beta,wes$SE)
estimates(wes_bc)
pdf(file = 'output/genomic_inflation/weighted_euro_smoker_GI.pdf',height = 18,width = 20)

print(traces(wes_bc,burnin = F))
print(posteriors(wes_bc))
print(plot(wes_bc,type='hist'))
print(plot(wes_bc,type='qq'))
dev.off()


wef_bc <- bacon(teststatistics = wef$Z,effectsizes = wef$Beta, standarderrors = wef$SE)
inflation(wef_bc)

#sigma.0 
#1.066241
wef$Beta_bacon <- es(wef_bc, corrected = TRUE)
wef$SE_bacon   <- se(wef_bc, corrected = TRUE)
wef$P_bacon    <- pval(wef_bc, corrected = TRUE)
c(
  uncorrected = lambda_gc_z(tstat(wef_bc, corrected = FALSE)),
  bacon_corrected = lambda_gc_z(tstat(wef_bc, corrected = TRUE))
)
wef_sig <- wef %>%
  filter(P_bacon <= bonf)
lambda_gc(wef$Beta,wef$SE)
pdf(file = 'output/genomic_inflation/weighted_euro_formersmoker_GI.pdf',height = 18,width = 20)
estimates(wef_bc)
traces(wef_bc,burnin = F)
posteriors(wef_bc)
plot(wef_bc,type='hist')
plot(wef_bc,type='qq')
dev.off()

wts_bc <- bacon(teststatistics = wts$Z,effectsizes = wts$Beta, standarderrors = wts$SE)
inflation(wts_bc)
# sigma.0 
#1.180236 
wts$Beta_bacon <- es(wts_bc, corrected = TRUE)
wts$SE_bacon   <- se(wts_bc, corrected = TRUE)
wts$P_bacon    <- pval(wts_bc, corrected = TRUE)

c(
  uncorrected = lambda_gc_z(tstat(wts_bc, corrected = FALSE)),
  bacon_corrected = lambda_gc_z(tstat(wts_bc, corrected = TRUE))
)

wts_sig <- wts %>%
  filter(P_bacon <= bonf)
lambda_gc(wts$Beta,wts$SE)
pdf(file = 'output/genomic_inflation/weighted_transanc_smoker_GI.pdf',height = 18,width = 20)
estimates(wts_bc)
traces(wts_bc,burnin = F)
posteriors(wts_bc)
plot(wts_bc,type='hist')
plot(wts_bc,type='qq')
dev.off()

wtf_bc <- bacon(teststatistics = wtf$Z,effectsizes = wtf$Beta, standarderrors = wtf$SE)
inflation(wtf_bc)
# sigma.0 
#1.086478 
wtf$Beta_bacon <- es(wtf_bc, corrected = TRUE)
wtf$SE_bacon   <- se(wtf_bc, corrected = TRUE)
wtf$P_bacon    <- pval(wtf_bc, corrected = TRUE)

c(
  uncorrected = lambda_gc_z(tstat(wtf_bc, corrected = FALSE)),
  bacon_corrected = lambda_gc_z(tstat(wtf_bc, corrected = TRUE))
)

wtf_sig <- wtf %>%
  filter(P_bacon <= bonf)
lambda_gc(wtf$Beta,wtf$SE)
estimates(wtf_bc)
pdf(file = 'output/genomic_inflation/weighted_transanc_formersmoker_GI.pdf',height = 18,width = 20)
traces(wtf_bc,burnin = F)
posteriors(wtf_bc)
plot(wtf_bc,type='hist')
plot(wtf_bc,type='qq')
dev.off()



library(data.table)
library(dplyr)

manifest_dir <- "references"

epic_manifest <- fread(
  file.path(manifest_dir, "EPIC.hg38.manifest.tsv.gz"),
  select = "Probe_ID"
)

hm450_manifest <- fread(
  file.path(manifest_dir, "HM450.hg38.manifest.tsv.gz"),
  select = "Probe_ID"
)

epic_probes <- unique(na.omit(epic_manifest$Probe_ID))
hm450_probes <- unique(na.omit(hm450_manifest$Probe_ID))

wes_sig <- wes_sig %>%
  mutate(
    CpG_ID  = sub("_.*$", "", V1),
    on_EPIC = CpG_ID %in% epic_probes,
    on_450K = CpG_ID %in% hm450_probes,
    array_overlap = case_when(
      on_EPIC & on_450K ~ "EPIC_AND_450K",
      on_EPIC           ~ "EPIC_ONLY",
      on_450K           ~ "450K_ONLY",
      TRUE              ~ "MSA_ONLY"
    )
  )

table(wes_sig$array_overlap)

wef_sig <- wef_sig %>%
  mutate(
    CpG_ID  = sub("_.*$", "", V1),
    on_EPIC = CpG_ID %in% epic_probes,
    on_450K = CpG_ID %in% hm450_probes,
    array_overlap = case_when(
      on_EPIC & on_450K ~ "EPIC_AND_450K",
      on_EPIC           ~ "EPIC_ONLY",
      on_450K           ~ "450K_ONLY",
      TRUE              ~ "MSA_ONLY"
    )
  )

table(wef_sig$array_overlap)


wts_sig <- wts_sig %>%
  mutate(
    CpG_ID  = sub("_.*$", "", V1),
    on_EPIC = CpG_ID %in% epic_probes,
    on_450K = CpG_ID %in% hm450_probes,
    array_overlap = case_when(
      on_EPIC & on_450K ~ "EPIC_AND_450K",
      on_EPIC           ~ "EPIC_ONLY",
      on_450K           ~ "450K_ONLY",
      TRUE              ~ "MSA_ONLY"
    )
  )

table(wts_sig$array_overlap)


wtf_sig <- wtf_sig %>%
  mutate(
    CpG_ID  = sub("_.*$", "", V1),
    on_EPIC = CpG_ID %in% epic_probes,
    on_450K = CpG_ID %in% hm450_probes,
    array_overlap = case_when(
      on_EPIC & on_450K ~ "EPIC_AND_450K",
      on_EPIC           ~ "EPIC_ONLY",
      on_450K           ~ "450K_ONLY",
      TRUE              ~ "MSA_ONLY"
    )
  )

table(wtf_sig$array_overlap)


all_sig <- bind_rows(
  "Smokers, Euro"  = wes_sig,
  "Former, Euro"   = wef_sig,
  "Smokers, Trans" = wts_sig,
  "Former, Trans"  = wtf_sig,
  .id = "Model"
)

msa_only_unique <- all_sig %>%
  filter(array_overlap == "MSA_ONLY") %>%
  distinct(CpG_ID, .keep_all = TRUE)

dim(msa_only_unique)
msa_only_unique


write.table(wes_sig,file='output/genomic_inflation/significant/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv',
            sep='\t',quote=F,na='',row.names = F)

write.table(wef_sig,file='output/genomic_inflation/significant/unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv',
            sep='\t',quote=F,na='',row.names = F)

write.table(wts_sig,file='output/genomic_inflation/significant/unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv',
            sep='\t',quote=F,na='',row.names = F)

write.table(wtf_sig,file='output/genomic_inflation/significant/unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv',
            sep='\t',quote=F,na='',row.names = F)


write.table(wes,file='output/genomic_inflation/full/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv',
            sep='\t',quote=F,na='',row.names = F)

write.table(wef,file='output/genomic_inflation/full/unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv',
            sep='\t',quote=F,na='',row.names = F)

write.table(wts,file='output/genomic_inflation/full/unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv',
            sep='\t',quote=F,na='',row.names = F)

write.table(wtf,file='output/genomic_inflation/full/unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv',
            sep='\t',quote=F,na='',row.names = F)




lambda_table <- data.frame(
  Model=c("Smokers, Euro", "Former smokers, Euro", "Smokers, Transancestry", "Former smokers, Transancestry"),
  Lambda_GC=c(1.480447, 1.208825, 1.491948, 1.255474)
)

write.table(lambda_table, "output/genomic_inflation/lambda_GC_summary.tsv", sep="\t", row.names=FALSE, quote=FALSE)
###################
## Run permutations

###################
## Run 100 permutations for Euro and Trans

library(survey)
library(lmtest)
library(parallel)
library(data.table)

perm_dir <- "output/perm"
dir.create(perm_dir, recursive = TRUE, showWarnings = FALSE)

n_perm <- 100

## Assign the data used in each observed model
beta_euro <- beta_pc
pheno_euro <- pheno_pc

beta_trans <- beta2
pheno_trans <- pheno

lambda_gc <- function(beta, se) {
  z <- beta / se
  median(z^2, na.rm = TRUE) / qchisq(0.5, df = 1)
}

f.RLM.Adj.Robust.wgt <- function(meth, dat1) {
  dat <- svydesign(id = ~Sample_ID, weights = ~asthma_weights, data = data.frame(y = meth, dat1))
  
  mdl <- svyglm(
    y ~ LH_smoke_status + ANYETS+ Gender + LH_BMI + AgeAtHomeVisit + STATE +
      B + NK + CD4T + CD8T + Mono + Neutro + Eosino +
      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
    design = dat
  )
  
  cf <- coeftest(mdl, vcov = vcov(mdl))
  coef_name <- "LH_smoke_status2"
  
  out <- c(
    Beta = cf[coef_name, "Estimate"],
    SE = cf[coef_name, "Std. Error"],
    P = cf[coef_name, "Pr(>|z|)"],
    N = sum(!is.na(meth)),
    N.exp = sum(!is.na(meth) & dat$variables$LH_smoke_status == 2),
    N.unexp = sum(!is.na(meth) & dat$variables$LH_smoke_status == 0)
  )
  
  return(out)
}

run_permutations <- function(beta_data, pheno_data, population, seed_offset, n_perm = 100) {
  output_dir <- file.path(perm_dir, population)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  stopifnot(nrow(beta_data) == nrow(pheno_data))
  
  lambda_perm <- data.frame(
    population = population,
    permutation = seq_len(n_perm),
    n_labels_changed = NA_integer_,
    lambda = NA_real_
  )
  
  for (b in seq_len(n_perm)) {
    set.seed(seed_offset + b)
    
    pheno_perm <- pheno_data
    smoke_perm <- pheno_data$LH_smoke_status
    
    valid <- !is.na(smoke_perm) & smoke_perm != "99"
    perm_groups <- interaction(pheno_data$LHSSCREEN_12142015, drop = TRUE)
    
    for (ii in split(which(valid), perm_groups[valid])) {
      smoke_perm[ii] <- sample(smoke_perm[ii])
    }
    
    pheno_perm$LH_smoke_status <- factor(smoke_perm, levels = levels(pheno_data$LH_smoke_status))
    pheno_perm$LH_smoke_status <- relevel(pheno_perm$LH_smoke_status, ref = "0")
    
    lambda_perm$n_labels_changed[b] <- sum(pheno_perm$LH_smoke_status != pheno_data$LH_smoke_status, na.rm = TRUE)
    
    perm_results <- parApply(cl, beta_data, 2, f.RLM.Adj.Robust.wgt, dat1 = pheno_perm)
    perm_results <- as.data.frame(t(perm_results))
    
    lambda_perm$lambda[b] <- lambda_gc(perm_results$Beta, perm_results$SE)
    
    fwrite(
      perm_results,
      file.path(output_dir, sprintf("weighted_%s_smokers_perm_%03d.tsv", population, b)),
      sep = "\t"
    )
    
    fwrite(lambda_perm[seq_len(b), ], file.path(output_dir, "permutation_lambdas.tsv"), sep = "\t")
    
    message(population, " permutation ", b, "/", n_perm, ": lambda = ", round(lambda_perm$lambda[b], 4))
  }
  
  return(lambda_perm)
}

## Run all 100 Euro permutations
lambda_euro <- run_permutations(
  beta_data = beta_euro,
  pheno_data = pheno_euro,
  population = "euro",
  seed_offset = 1000,
  n_perm = n_perm
)

## Run all 100 Trans permutations
lambda_trans <- run_permutations(
  beta_data = beta_trans,
  pheno_data = pheno_trans,
  population = "trans",
  seed_offset = 2000,
  n_perm = n_perm
)

## Combine both populations
lambda_all <- rbind(lambda_euro, lambda_trans)

fwrite(lambda_all, file.path(perm_dir, "all_permutation_lambdas.tsv"), sep = "\t")

## Summarize using the updated observed lambdas
observed_lambdas <- data.table(
  population = c("euro", "trans"),
  observed_lambda = c(1.480447, 1.491948)
)

lambda_summary <- merge(as.data.table(lambda_all), observed_lambdas, by = "population")

lambda_summary <- lambda_summary[, .(
  n_permutations = .N,
  observed_lambda = first(observed_lambda),
  permutation_min = min(lambda, na.rm = TRUE),
  permutation_median = median(lambda, na.rm = TRUE),
  permutation_mean = mean(lambda, na.rm = TRUE),
  permutation_q025 = quantile(lambda, 0.025, na.rm = TRUE),
  permutation_q975 = quantile(lambda, 0.975, na.rm = TRUE),
  permutation_max = max(lambda, na.rm = TRUE),
  n_at_or_above_observed = sum(lambda >= first(observed_lambda), na.rm = TRUE),
  empirical_P = (1 + sum(lambda >= first(observed_lambda), na.rm = TRUE)) / (.N + 1)
), by = population]

fwrite(lambda_summary, file.path(perm_dir, "permutation_summary.tsv"), sep = "\t")

print(lambda_summary)

lambda_trans <- fread("output/perm/trans/permutation_lambdas.tsv")
data.frame(
 Model="Smokers, Transancestry",
 Observed_lambda=1.491948,
 Permutation_median=median(lambda_trans$lambda),
 Permutation_mean=mean(lambda_trans$lambda),
 Permutation_max=max(lambda_trans$lambda),
 N_at_or_above_observed=sum(lambda_trans$lambda >= 1.491948),
 Empirical_P=(1 + sum(lambda_trans$lambda >= 1.491948)) / 101)

lambda_euro <- fread("output/perm/euro/permutation_lambdas.tsv")
data.frame(
   Model="Smokers, Euro",
   Observed_lambda=1.480447,
   Permutation_median=median(lambda_euro$lambda),
   Permutation_mean=mean(lambda_euro$lambda),
   Permutation_max=max(lambda_euro$lambda),
   N_at_or_above_observed=sum(lambda_euro$lambda >= 1.480447),
   Empirical_P=(1 + sum(lambda_euro$lambda >= 1.480447)) / 101)
 
library(data.table)

perm_dir <- "output/perm"
lambda_euro <- fread("output/perm/euro/permutation_lambdas.tsv")
lambda_trans <- fread("output/perm/trans/permutation_lambdas.tsv")

lambda_euro[, population := "euro"]
lambda_trans[, population := "trans"]

lambda_euro <- unique(lambda_euro, by="permutation")
lambda_trans <- unique(lambda_trans, by="permutation")
stopifnot(nrow(lambda_euro)==100, nrow(lambda_trans)==100)

lambda_all <- rbindlist(list(lambda_euro, lambda_trans), fill=TRUE)
