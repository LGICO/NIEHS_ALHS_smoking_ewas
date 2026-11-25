library(dplyr)
load('751_clean_pc.rdata')

### EWAS
library(survey)
library(MASS)
library(lmtest)
library(sandwich)
library(parallel)

dir.create("output/reversion", recursive=TRUE, showWarnings=FALSE)

## A former-versus-never difference within +/- 0.01 beta units is treated as equivalent.
## Change this only if a different biologically meaningful equivalence margin is chosen.
equivalence_margin <- 0.01
equivalence_alpha <- 0.05

cl <- makeCluster(getOption("cl.cores", 60))
clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(sandwich))
clusterEvalQ(cl,library(lmtest))
clusterEvalQ(cl,library(survey))

################################
## Phenotype preparation
new_pheno <- data.table::fread('samplesheets/lhs28901_Methylation_Phenotype_2025APR29.txt',sep=',')
np <- new_pheno %>% filter(LHID %in% pheno$LHID, COMMENT != 'Replicate')
x <- colnames(np)[!colnames(np) %in% colnames(pheno)]
np <- np %>% dplyr::select(LHID, partid, all_of(x))
pheno <- left_join(pheno, np, by=c('LHID','partid'))

pheno <- pheno %>% mutate(LEL19=ifelse(LEL19=='Yes',1,ifelse(LEL19=='No',0,99)))
pheno$LEL19 <- as.factor(pheno$LEL19)

pheno <- pheno %>% mutate(LH_smoke_status=case_when(LH_smoke_status=='Never'~0,
                                                    LH_smoke_status=='Former'~1,
                                                    LH_smoke_status=='Current'~2,
                                                    LH_smoke_status=='M'~99))
pheno$LH_smoke_status <- factor(pheno$LH_smoke_status, levels=c(0,1,2,99))

pheno <- pheno %>% mutate(European=case_when(European=='European'~1,
                                             European=='Non_European'~0,
                                             European=='GWAS sample failure'~99))
pheno$European <- factor(pheno$European, levels=c(0,1,99))
pheno$Gender <- as.factor(pheno$Gender)
pheno$inrdsx <- as.factor(pheno$inrdsx)

## Preserve the unknown ANYETS category instead of removing those samples.
pheno <- pheno %>% mutate(ANYETS=case_when(ANYETS=='No'~0, ANYETS=='Yes'~1, .default=99))
pheno <- pheno %>% mutate(LHSSCREEN_12142015=case_when(LHSSCREEN_12142015=='Control'~0,
                                                       LHSSCREEN_12142015=='Case'~1))
pheno <- pheno %>% mutate(Gender=case_when(Gender=='male'~0, Gender=='female'~1))
pheno <- pheno %>% mutate(STATE=case_when(STATE=='IA'~0, STATE=='NC'~1))
pheno$LHSSCREEN_12142015 <- as.factor(pheno$LHSSCREEN_12142015)
pheno$ANYETS <- factor(pheno$ANYETS, levels=c(0,1,99))
pheno$Gender <- as.factor(pheno$Gender)
pheno$STATE <- as.factor(pheno$STATE)

## One joint smoking variable with never smokers as the reference.
## Four cessation intervals are used for the primary analysis. Combining
## 11-15 and 16-20 years avoids the small European 11-15-year group (n=18).
pheno <- pheno %>% mutate(reversion_group=case_when(
  LH_smoke_status==0 ~ 'Never',
  LH_smoke_status==2 ~ 'Current',
  LH_smoke_status==1 & LH_TimeSinceSmked_v1<=10 ~ 'Former_le10',
  LH_smoke_status==1 & LH_TimeSinceSmked_v1>10 & LH_TimeSinceSmked_v1<=20 ~ 'Former_11to20',
  LH_smoke_status==1 & LH_TimeSinceSmked_v1>20 & LH_TimeSinceSmked_v1<=30 ~ 'Former_21to30',
  LH_smoke_status==1 & LH_TimeSinceSmked_v1>30 ~ 'Former_gt30',
  LH_smoke_status==99 ~ 'Unknown'))
pheno$reversion_group <- factor(pheno$reversion_group,
                                levels=c('Never','Current','Former_le10','Former_11to20',
                                         'Former_21to30','Former_gt30','Unknown'))


former <- pheno %>%
  filter(as.character(LH_smoke_status) %in% c("Former", "1")) %>%
  mutate(
    years_since_quit=as.numeric(as.character(LH_TimeSinceSmked_v1)),
    cessation_interval=case_when(
      years_since_quit <= 10 ~ "≤10 years",
      years_since_quit <= 20 ~ "11–20 years",
      years_since_quit <= 30 ~ "21–30 years",
      years_since_quit > 30 ~ ">30 years"
    )
  )

## Overall summary
former %>%
  summarise(
    N=dplyr::n(),
    N_with_cessation_time=sum(!is.na(years_since_quit)),
    Missing=sum(is.na(years_since_quit)),
    Median=median(years_since_quit, na.rm=TRUE),
    Q1=quantile(years_since_quit, 0.25, na.rm=TRUE),
    Q3=quantile(years_since_quit, 0.75, na.rm=TRUE),
    Minimum=min(years_since_quit, na.rm=TRUE),
    Maximum=max(years_since_quit, na.rm=TRUE)
  )

## Cessation-interval counts overall and among European participants
cessation_counts <- bind_rows(
  former %>% dplyr::mutate(Population="Overall"),
  former %>% dplyr::filter(as.character(European) %in% c("European", "1")) %>% dplyr::mutate(Population="European")
) %>%
  dplyr::filter(!is.na(cessation_interval)) %>%
  dplyr::count(Population, cessation_interval) %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(Percent=100*n/sum(n)) %>%
  dplyr::ungroup()

print(cessation_counts)

cessation_bin_counts <- pheno %>% dplyr::filter(LH_smoke_status==1) %>% dplyr::count(reversion_group,name='N')
print(cessation_bin_counts)
write.csv(cessation_bin_counts,'output/reversion/cessation_bin_counts.csv',row.names=FALSE)

cessation_bin_counts_by_population <- bind_rows(
  pheno %>% dplyr::filter(LH_smoke_status==1) %>% dplyr::count(reversion_group,name='N') %>% dplyr::mutate(Population='Transancestry'),
  pheno %>% dplyr::filter(LH_smoke_status==1,European==1) %>% dplyr::count(reversion_group,name='N') %>% dplyr::mutate(Population='European')
) %>% dplyr::select(Population,reversion_group,N)
write.csv(cessation_bin_counts_by_population,'output/reversion/cessation_bin_counts_by_population.csv',row.names=FALSE)

## Unimputed methylation values, aligned explicitly to phenotype order.
betap <- t(betas_miss_masked)
stopifnot(all(pheno$Sample_ID %in% rownames(betap)))
betap <- betap[match(pheno$Sample_ID, rownames(betap)), , drop=FALSE]

pheno_pc <- pheno %>% filter(European==1, !is.na(PC1))
betap_pc <- betap[match(pheno_pc$Sample_ID, rownames(betap)), , drop=FALSE]

stopifnot(identical(rownames(betap), pheno$Sample_ID))
stopifnot(identical(rownames(betap_pc), pheno_pc$Sample_ID))

## Probe annotations added to every filtered reversion output before it is written.
manifest_anno <- data.table::fread(
  'references/MSA.hg38.manifest.gencode.v41.tsv.gz'
) %>%
  transmute(Probe_ID=probeID,CpG_chrm,CpG_beg,CpG_end,probe_strand,genesUniq,
            geneNames,transcriptTypes,transcriptIDs,distToTSS) %>%
  distinct(Probe_ID,.keep_all=TRUE)

############################################################################################################
## Joint reversion model
fit_reversion <- function(meth, dat1, weighted=FALSE, trans=FALSE,
                          eq_margin=0.01, eq_alpha=0.05) {
  
  former_levels <- c('Former_le10','Former_11to20','Former_21to30','Former_gt30')
  previous_levels <- c('Current','Former_le10','Former_11to20','Former_21to30')
  out_names <- c('N','N_never','N_current','N_former_le10','N_former_11to20',
                 'N_former_21to30','N_former_gt30',
                 'Current_vs_never_Beta','Current_vs_never_SE',
                 'Current_vs_never_tstat','Current_vs_never_P')
  
  for (g in former_levels) {
    out_names <- c(out_names,
                   paste0(g,'_vs_never_Beta'), paste0(g,'_vs_never_SE'),
                   paste0(g,'_vs_never_tstat'), paste0(g,'_vs_never_P'),
                   paste0(g,'_vs_never_CI90_lower'), paste0(g,'_vs_never_CI90_upper'),
                   paste0(g,'_vs_never_TOST_P'), paste0(g,'_equivalent_to_never'),
                   paste0(g,'_vs_current_Beta'), paste0(g,'_vs_current_SE'),
                   paste0(g,'_vs_current_tstat'), paste0(g,'_vs_current_P'),
                   paste0(g,'_vs_current_CI90_lower'), paste0(g,'_vs_current_CI90_upper'),
                   paste0(g,'_vs_current_TOST_P'), paste0(g,'_equivalent_to_current'),
                   paste0(g,'_vs_previous_Beta'), paste0(g,'_vs_previous_SE'),
                   paste0(g,'_vs_previous_tstat'), paste0(g,'_vs_previous_P'),
                   paste0(g,'_percent_reversion'))
  }
  
  empty <- setNames(rep(NA_real_, length(out_names)), out_names)
  
  tryCatch({
    model_data <- data.frame(y=as.numeric(meth), dat1, check.names=FALSE)
    
    if (weighted) {
      dat <- svydesign(id=~Sample_ID, weights=~asthma_weights, data=model_data)
    } else {
      dat <- svydesign(id=~Sample_ID, data=model_data)
    }
    
    rhs <- c('reversion_group','ANYETS')
    if (trans) rhs <- c(rhs,'European')
    if (!weighted || !trans) rhs <- c(rhs,'LHSSCREEN_12142015')
    rhs <- c(rhs,'Gender','LH_BMI','AgeAtHomeVisit','STATE',
             'B','NK','CD4T','CD8T','Mono','Neutro','Eosino',
             paste0('PC',1:10))
    
    mdl.formula <- reformulate(rhs, response='y')
    mdl <- svyglm(formula=mdl.formula, design=dat)
    b <- coef(mdl)
    V <- vcov(mdl)
    df <- mdl$df.residual
    
    p_two_sided <- function(tvalue) {
      if (is.finite(df) && df>0) 2*pt(abs(tvalue),df=df,lower.tail=FALSE) else 2*pnorm(abs(tvalue),lower.tail=FALSE)
    }
    
    p_upper_tail <- function(tvalue) {
      if (is.finite(df) && df>0) pt(tvalue,df=df,lower.tail=FALSE) else pnorm(tvalue,lower.tail=FALSE)
    }
    
    p_lower_tail <- function(tvalue) {
      if (is.finite(df) && df>0) pt(tvalue,df=df,lower.tail=TRUE) else pnorm(tvalue,lower.tail=TRUE)
    }
    
    critical90 <- if (is.finite(df) && df>0) qt(0.95,df=df) else qnorm(0.95)
    
    coefficient_result <- function(coef_name) {
      if (!coef_name %in% names(b)) return(c(Beta=NA,SE=NA,tstat=NA,P=NA))
      est <- unname(b[coef_name])
      se <- sqrt(V[coef_name,coef_name])
      tstat <- est/se
      c(Beta=est,SE=se,tstat=tstat,P=p_two_sided(tstat))
    }
    
    contrast_result <- function(first_name, second_name) {
      if (!all(c(first_name,second_name) %in% names(b))) return(c(Beta=NA,SE=NA,tstat=NA,P=NA))
      L <- setNames(rep(0,length(b)),names(b))
      L[first_name] <- 1
      L[second_name] <- -1
      est <- sum(L*b)
      se <- sqrt(as.numeric(t(L)%*%V%*%L))
      tstat <- est/se
      c(Beta=est,SE=se,tstat=tstat,P=p_two_sided(tstat))
    }
    
    out <- empty
    out['N'] <- nobs(mdl)
    out['N_never'] <- sum(!is.na(model_data$y) & model_data$reversion_group=='Never')
    out['N_current'] <- sum(!is.na(model_data$y) & model_data$reversion_group=='Current')
    out['N_former_le10'] <- sum(!is.na(model_data$y) & model_data$reversion_group=='Former_le10')
    out['N_former_11to20'] <- sum(!is.na(model_data$y) & model_data$reversion_group=='Former_11to20')
    out['N_former_21to30'] <- sum(!is.na(model_data$y) & model_data$reversion_group=='Former_21to30')
    out['N_former_gt30'] <- sum(!is.na(model_data$y) & model_data$reversion_group=='Former_gt30')
    
    current_name <- 'reversion_groupCurrent'
    current <- coefficient_result(current_name)
    out[c('Current_vs_never_Beta','Current_vs_never_SE',
          'Current_vs_never_tstat','Current_vs_never_P')] <- current
    
    for (i in seq_along(former_levels)) {
      g <- former_levels[i]
      former_name <- paste0('reversion_group',g)
      previous_name <- if (previous_levels[i]=='Current') current_name else paste0('reversion_group',previous_levels[i])
      former <- coefficient_result(former_name)
      former_current <- contrast_result(former_name,current_name)
      former_previous <- contrast_result(former_name,previous_name)
      
      ci90_lower <- former['Beta']-critical90*former['SE']
      ci90_upper <- former['Beta']+critical90*former['SE']
      tost_lower <- (former['Beta']+eq_margin)/former['SE']
      tost_upper <- (former['Beta']-eq_margin)/former['SE']
      tost_p <- max(p_upper_tail(tost_lower),p_lower_tail(tost_upper))
      equivalent <- as.numeric(!is.na(tost_p) && tost_p<eq_alpha)
      current_ci90_lower <- former_current['Beta']-critical90*former_current['SE']
      current_ci90_upper <- former_current['Beta']+critical90*former_current['SE']
      current_tost_lower <- (former_current['Beta']+eq_margin)/former_current['SE']
      current_tost_upper <- (former_current['Beta']-eq_margin)/former_current['SE']
      current_tost_p <- max(p_upper_tail(current_tost_lower),p_lower_tail(current_tost_upper))
      equivalent_current <- as.numeric(!is.na(current_tost_p) && current_tost_p<eq_alpha)
      percent_reversion <- 100*(1-former['Beta']/current['Beta'])
      
      out[paste0(g,'_vs_never_Beta')] <- former['Beta']
      out[paste0(g,'_vs_never_SE')] <- former['SE']
      out[paste0(g,'_vs_never_tstat')] <- former['tstat']
      out[paste0(g,'_vs_never_P')] <- former['P']
      out[paste0(g,'_vs_never_CI90_lower')] <- ci90_lower
      out[paste0(g,'_vs_never_CI90_upper')] <- ci90_upper
      out[paste0(g,'_vs_never_TOST_P')] <- tost_p
      out[paste0(g,'_equivalent_to_never')] <- equivalent
      out[paste0(g,'_vs_current_Beta')] <- former_current['Beta']
      out[paste0(g,'_vs_current_SE')] <- former_current['SE']
      out[paste0(g,'_vs_current_tstat')] <- former_current['tstat']
      out[paste0(g,'_vs_current_P')] <- former_current['P']
      out[paste0(g,'_vs_current_CI90_lower')] <- current_ci90_lower
      out[paste0(g,'_vs_current_CI90_upper')] <- current_ci90_upper
      out[paste0(g,'_vs_current_TOST_P')] <- current_tost_p
      out[paste0(g,'_equivalent_to_current')] <- equivalent_current
      out[paste0(g,'_vs_previous_Beta')] <- former_previous['Beta']
      out[paste0(g,'_vs_previous_SE')] <- former_previous['SE']
      out[paste0(g,'_vs_previous_tstat')] <- former_previous['tstat']
      out[paste0(g,'_vs_previous_P')] <- former_previous['P']
      out[paste0(g,'_percent_reversion')] <- percent_reversion
    }
    
    out
  }, error=function(e) empty)
}

############################################################################################################
## Create candidate-DMP and reversion-classification files
interval_info <- data.frame(
  prefix=c('Former_le10','Former_11to20','Former_21to30','Former_gt30'),
  previous_prefix=c('Current','Former_le10','Former_11to20','Former_21to30'),
  interval=c('Former <=10 years','Former 11-20 years','Former 21-30 years','Former >30 years'),
  previous_interval=c('Current','Former <=10 years','Former 11-20 years','Former 21-30 years'),
  interval_order=1:4,
  stringsAsFactors=FALSE
)

create_reversion_tables <- function(results, model_name) {
  results$Probe_ID <- rownames(results)
  bonf <- 0.05/nrow(results)
  candidates <- results %>% filter(Current_vs_never_P<=bonf)
  
  long_list <- lapply(seq_len(nrow(interval_info)), function(i) {
    g <- interval_info$prefix[i]
    x <- data.frame(
      Probe_ID=candidates$Probe_ID,
      Model=model_name,
      Interval=interval_info$interval[i],
      Previous_interval=interval_info$previous_interval[i],
      Interval_order=interval_info$interval_order[i],
      Current_vs_never_Beta=candidates$Current_vs_never_Beta,
      Current_vs_never_SE=candidates$Current_vs_never_SE,
      Current_vs_never_P=candidates$Current_vs_never_P,
      Former_vs_never_Beta=candidates[[paste0(g,'_vs_never_Beta')]],
      Former_vs_never_SE=candidates[[paste0(g,'_vs_never_SE')]],
      Former_vs_never_P=candidates[[paste0(g,'_vs_never_P')]],
      Former_vs_never_CI90_lower=candidates[[paste0(g,'_vs_never_CI90_lower')]],
      Former_vs_never_CI90_upper=candidates[[paste0(g,'_vs_never_CI90_upper')]],
      Equivalence_TOST_P=candidates[[paste0(g,'_vs_never_TOST_P')]],
      Former_vs_current_Beta=candidates[[paste0(g,'_vs_current_Beta')]],
      Former_vs_current_SE=candidates[[paste0(g,'_vs_current_SE')]],
      Former_vs_current_P=candidates[[paste0(g,'_vs_current_P')]],
      Former_vs_current_CI90_lower=candidates[[paste0(g,'_vs_current_CI90_lower')]],
      Former_vs_current_CI90_upper=candidates[[paste0(g,'_vs_current_CI90_upper')]],
      Current_equivalence_TOST_P=candidates[[paste0(g,'_vs_current_TOST_P')]],
      Former_vs_previous_Beta=candidates[[paste0(g,'_vs_previous_Beta')]],
      Former_vs_previous_SE=candidates[[paste0(g,'_vs_previous_SE')]],
      Former_vs_previous_P=candidates[[paste0(g,'_vs_previous_P')]],
      Percent_reversion=candidates[[paste0(g,'_percent_reversion')]],
      stringsAsFactors=FALSE
    )
    
    x$Former_vs_never_FDR <- p.adjust(x$Former_vs_never_P,method='BH')
    x$Equivalence_TOST_FDR <- p.adjust(x$Equivalence_TOST_P,method='BH')
    x$Former_vs_current_FDR <- p.adjust(x$Former_vs_current_P,method='BH')
    x$Current_equivalence_TOST_FDR <- p.adjust(x$Current_equivalence_TOST_P,method='BH')
    x$Former_vs_previous_FDR <- p.adjust(x$Former_vs_previous_P,method='BH')
    x$Recovery_direction <- with(x, Former_vs_current_Beta*Current_vs_never_Beta<0 &
                                   abs(Former_vs_never_Beta)<abs(Current_vs_never_Beta))
    x$Equivalent_to_never <- x$Equivalence_TOST_FDR<0.05
    x$Equivalent_to_current <- x$Current_equivalence_TOST_FDR<0.05
    x$Different_from_current <- x$Former_vs_current_FDR<0.05
    
    x$Reversion_classification <- case_when(
      x$Equivalent_to_never & x$Different_from_current & x$Recovery_direction ~ 'Full/near-complete reversion',
      !x$Equivalent_to_never & x$Different_from_current & x$Recovery_direction & x$Former_vs_never_FDR<0.05 ~ 'Partial reversion with residual difference',
      !x$Equivalent_to_never & x$Different_from_current & x$Recovery_direction ~ 'Partial reversion; residual unresolved',
      x$Former_vs_never_FDR<0.05 & x$Equivalent_to_current ~ 'Persistent/equivalent to current',
      x$Former_vs_never_FDR<0.05 & !x$Different_from_current ~ 'Persistent/no detectable recovery',
      !x$Recovery_direction & abs(x$Former_vs_never_Beta)>=abs(x$Current_vs_never_Beta) ~ 'Non-monotonic/no attenuation',
      TRUE ~ 'Indeterminate'
    )
    x
  })
  
  reversion_long <- bind_rows(long_list) %>% arrange(Probe_ID,Interval_order)
  
  trajectory <- reversion_long %>% group_by(Probe_ID) %>% summarise(
    Monotonic_toward_never=all(diff(abs(Former_vs_never_Beta))<=0,na.rm=TRUE),
    First_equivalent_interval=ifelse(any(Reversion_classification=='Full/near-complete reversion'),
                                     Interval[min(which(Reversion_classification=='Full/near-complete reversion'))],NA_character_),
    .groups='drop')
  
  candidates <- candidates %>% left_join(trajectory,by='Probe_ID')
  summary_table <- reversion_long %>% count(Model,Interval,Reversion_classification,name='N')
  
  results_annotated <- results %>% left_join(manifest_anno,by='Probe_ID') %>%
    dplyr::select(Probe_ID,CpG_chrm,CpG_beg,CpG_end,probe_strand,genesUniq,dplyr::everything())
  candidates_annotated <- candidates %>% left_join(manifest_anno,by='Probe_ID') %>%
    dplyr::select(Probe_ID,CpG_chrm,CpG_beg,CpG_end,probe_strand,genesUniq,dplyr::everything())
  reversion_long_annotated <- reversion_long %>% left_join(manifest_anno,by='Probe_ID') %>%
    dplyr::select(Probe_ID,CpG_chrm,CpG_beg,CpG_end,probe_strand,genesUniq,dplyr::everything())
  persistent_annotated <- reversion_long_annotated %>%
    filter(Former_vs_never_FDR<0.05, !Equivalent_to_never) %>%
    arrange(Interval_order,genesUniq,Probe_ID)
  no_recovery_annotated <- reversion_long_annotated %>%
    filter(grepl('^Persistent',Reversion_classification)) %>%
    arrange(Interval_order,genesUniq,Probe_ID)
  focus_annotated <- reversion_long_annotated %>%
    filter(grepl('AHRR|ECEL1P1',genesUniq,ignore.case=TRUE)) %>%
    arrange(genesUniq,Probe_ID,Interval_order)
  
  write.csv(results_annotated,paste0('output/reversion/',model_name,'_all_cpgs_annotated.csv'),row.names=FALSE)
  write.csv(candidates_annotated,paste0('output/reversion/',model_name,'_current_smoking_dmps_annotated.csv'),row.names=FALSE)
  write.csv(reversion_long_annotated,paste0('output/reversion/',model_name,'_reversion_long_annotated.csv'),row.names=FALSE)
  write.csv(persistent_annotated,paste0('output/reversion/',model_name,'_persistent_residual_annotated.csv'),row.names=FALSE)
  write.csv(no_recovery_annotated,paste0('output/reversion/',model_name,'_no_detectable_recovery_annotated.csv'),row.names=FALSE)
  write.csv(focus_annotated,paste0('output/reversion/',model_name,'_AHRR_ECEL1P1_annotated.csv'),row.names=FALSE)
  write.csv(summary_table,paste0('output/reversion/',model_name,'_classification_summary.csv'),row.names=FALSE)
  
  summary_table
}

run_reversion_model <- function(meth, dat1, model_name, weighted=FALSE, trans=FALSE) {
  results <- parApply(cl,meth,2,fit_reversion,dat1=dat1,weighted=weighted,trans=trans,
                      eq_margin=equivalence_margin,eq_alpha=equivalence_alpha)
  results <- as.data.frame(t(results))
  create_reversion_tables(results,model_name)
}

############################################################################################################
## Model 1: unweighted European
summary_model1 <- run_reversion_model(betap_pc,pheno_pc,'model1_unweighted_euro',
                                      weighted=FALSE,trans=FALSE)

############################################################################################################
## Model 10: unweighted transancestry
summary_model10 <- run_reversion_model(betap,pheno,'model10_unweighted_trans',
                                       weighted=FALSE,trans=TRUE)

############################################################################################################
## Asthma weights used in the existing European weighted model
p_asthma_overall <- 0.06852481
n_asthma <- 186
n_nonasthma <- 4579
p_asthma_nested <- n_asthma/(n_asthma+n_nonasthma)
p_nonasthma_overall <- 1-p_asthma_overall
p_nonasthma_nested <- 1-p_asthma_nested
w_asthma <- 1/(p_asthma_nested/p_asthma_overall)
w_nonasthma <- 1/(p_nonasthma_nested/p_nonasthma_overall)
pheno_pc$asthma_weights[pheno_pc$LHSSCREEN_12142015==0] <- w_nonasthma
pheno_pc$asthma_weights[pheno_pc$LHSSCREEN_12142015==1] <- w_asthma

############################################################################################################
## Model 3: weighted European
summary_model3 <- run_reversion_model(betap_pc,pheno_pc,'model3_weighted_euro',
                                      weighted=TRUE,trans=FALSE)

############################################################################################################
## Asthma weights used in the existing transancestry weighted model
n_asthma <- sum(pheno$LHSSCREEN_12142015==1)
n_nonasthma <- sum(pheno$LHSSCREEN_12142015==0)
p_asthma_nested <- n_asthma/(n_asthma+n_nonasthma)
p_nonasthma_overall <- 1-p_asthma_overall
p_nonasthma_nested <- 1-p_asthma_nested
w_asthma <- 1/(p_asthma_nested/p_asthma_overall)
w_nonasthma <- 1/(p_nonasthma_nested/p_nonasthma_overall)
pheno$asthma_weights[pheno$LHSSCREEN_12142015==0] <- w_nonasthma
pheno$asthma_weights[pheno$LHSSCREEN_12142015==1] <- w_asthma

############################################################################################################
## Model 12: weighted transancestry
summary_model12 <- run_reversion_model(betap,pheno,'model12_weighted_trans',
                                       weighted=TRUE,trans=TRUE)

############################################################################################################
## Focused investigation of the weighted models
## Probe-design tables retain suffixes such as _TC11. Unique-CpG tables collapse
## those designs to the base cg identifier and keep the design with the smallest FDR.
read_weighted_reversion <- function(model_name,model_label) {
  read.csv(paste0('output/reversion/',model_name,'_reversion_long_annotated.csv'),
           check.names=FALSE) %>% mutate(Weighted_model=model_label)
}

weighted_long <- bind_rows(
  read_weighted_reversion('model3_weighted_euro','Weighted European'),
  read_weighted_reversion('model12_weighted_trans','Weighted transancestry')
) %>% mutate(
  CpG_ID=sub('_.*$','',Probe_ID),
  Residual_difference=!is.na(Former_vs_never_FDR) & Former_vs_never_FDR<0.05 & !Equivalent_to_never,
  Focus_gene=case_when(
    grepl('ECEL1P1',genesUniq,ignore.case=TRUE) ~ 'ECEL1P1',
    grepl('AHRR',genesUniq,ignore.case=TRUE) ~ 'AHRR',
    TRUE ~ NA_character_))

weighted_focus_all <- weighted_long %>% filter(!is.na(Focus_gene)) %>%
  arrange(Weighted_model,Focus_gene,CpG_ID,Probe_ID,Interval_order)

weighted_focus_residual_designs <- weighted_focus_all %>%
  filter(Residual_difference) %>%
  arrange(Weighted_model,Focus_gene,CpG_ID,Probe_ID,Interval_order)

weighted_focus_residual_unique <- weighted_focus_residual_designs %>%
  group_by(Weighted_model,Focus_gene,CpG_ID,Interval,Interval_order) %>%
  slice_min(Former_vs_never_FDR,n=1,with_ties=FALSE) %>% ungroup() %>%
  arrange(Weighted_model,Focus_gene,CpG_ID,Interval_order)

weighted_focus_summary <- weighted_focus_residual_designs %>%
  group_by(Weighted_model,Focus_gene,Interval,Interval_order) %>%
  summarise(Probe_designs=n_distinct(Probe_ID),Unique_CpGs=n_distinct(CpG_ID),.groups='drop') %>%
  arrange(Weighted_model,Focus_gene,Interval_order)

weighted_focus_gt30 <- weighted_focus_all %>%
  filter(Interval=='Former >30 years') %>%
  arrange(Weighted_model,Focus_gene,CpG_ID,Probe_ID) %>%
  dplyr::select(Weighted_model,Focus_gene,CpG_ID,Probe_ID,genesUniq,
                Current_vs_never_Beta,Former_vs_never_Beta,Former_vs_never_P,
                Former_vs_never_FDR,Equivalent_to_never,Former_vs_current_Beta,
                Former_vs_current_FDR,Percent_reversion,Reversion_classification,dplyr::everything())

write.csv(weighted_focus_all,'output/reversion/weighted_AHRR_ECEL1P1_all_probe_designs.csv',row.names=FALSE)
write.csv(weighted_focus_residual_designs,'output/reversion/weighted_AHRR_ECEL1P1_residual_probe_designs.csv',row.names=FALSE)
write.csv(weighted_focus_residual_unique,'output/reversion/weighted_AHRR_ECEL1P1_residual_unique_cpgs.csv',row.names=FALSE)
write.csv(weighted_focus_summary,'output/reversion/weighted_AHRR_ECEL1P1_interval_summary.csv',row.names=FALSE)
write.csv(weighted_focus_gt30,'output/reversion/weighted_AHRR_ECEL1P1_gt30_status.csv',row.names=FALSE)

## Compare residual unique CpGs between the two weighted populations at each interval.
weighted_residual_unique <- weighted_long %>% filter(Residual_difference) %>%
  distinct(Weighted_model,Interval,Interval_order,CpG_ID)

weighted_euro_residual <- weighted_residual_unique %>%
  filter(Weighted_model=='Weighted European') %>%
  transmute(Interval,Interval_order,CpG_ID,European=1L)

weighted_trans_residual <- weighted_residual_unique %>%
  filter(Weighted_model=='Weighted transancestry') %>%
  transmute(Interval,Interval_order,CpG_ID,Transancestry=1L)

weighted_residual_overlap <- full_join(weighted_euro_residual,weighted_trans_residual,
                                       by=c('Interval','Interval_order','CpG_ID')) %>%
  mutate(European=ifelse(is.na(European),0L,European),
         Transancestry=ifelse(is.na(Transancestry),0L,Transancestry),
         Overlap=European==1L & Transancestry==1L) %>%
  arrange(Interval_order,CpG_ID)

weighted_residual_overlap_summary <- weighted_residual_overlap %>%
  group_by(Interval,Interval_order) %>%
  summarise(European=sum(European),Transancestry=sum(Transancestry),
            Overlap=sum(Overlap),Union=n(),
            Percent_European_retained=100*Overlap/European,.groups='drop') %>%
  arrange(Interval_order)

write.csv(weighted_residual_overlap,'output/reversion/weighted_residual_unique_cpg_overlap.csv',row.names=FALSE)
write.csv(weighted_residual_overlap_summary,'output/reversion/weighted_residual_unique_cpg_overlap_summary.csv',row.names=FALSE)

############################################################################################################
## Combined classification summary
reversion_summary <- bind_rows(summary_model1,summary_model10,summary_model3,summary_model12)
write.csv(reversion_summary,'output/reversion/reversion_classification_summary.csv',row.names=FALSE)

stopCluster(cl)


library(dplyr)
library(ggplot2)
library(cowplot)

root <- '.'
reversion_dir <- file.path(root,'output','reversion')

interval_levels <- c('Former <=10 years','Former 11-20 years','Former 21-30 years','Former >30 years')

############################################################################################################
## Panel A: distribution of reversion classifications in the weighted models
reversion_summary <- read.csv(file.path(reversion_dir,'reversion_classification_summary.csv'),check.names=FALSE)
bin_counts <- read.csv(file.path(reversion_dir,'cessation_bin_counts_by_population.csv'),check.names=FALSE)

plot_a_data <- reversion_summary %>%
  filter(Model %in% c('model3_weighted_euro','model12_weighted_trans')) %>%
  mutate(Population=ifelse(Model=='model3_weighted_euro','European','Transancestry'),
         Interval=factor(Interval,levels=interval_levels),
         Classification=case_when(
           Reversion_classification=='Full/near-complete reversion' ~ 'Formal reversion',
           Reversion_classification=='Partial reversion with residual difference' ~ 'Partial reversion; residual remains',
           grepl('^Persistent',Reversion_classification) ~ 'Persistent/no detected recovery',
           TRUE ~ 'Unresolved/other')) %>%
  group_by(Population,Interval,Classification) %>% summarise(N=sum(N),.groups='drop') %>%
  group_by(Population,Interval) %>% mutate(Percent=100*N/sum(N)) %>% ungroup()

bin_counts <- bin_counts %>%
  mutate(Interval=case_when(
    reversion_group=='Former_le10' ~ 'Former <=10 years',
    reversion_group=='Former_11to20' ~ 'Former 11-20 years',
    reversion_group=='Former_21to30' ~ 'Former 21-30 years',
    reversion_group=='Former_gt30' ~ 'Former >30 years'),
    Interval=factor(Interval,levels=interval_levels))

plot_a <- ggplot(plot_a_data,aes(x=Interval,y=Percent,fill=Classification)) +
  geom_col(width=0.78,color='white',linewidth=0.3) +
  geom_text(data=plot_a_data %>% filter(Percent>=7),aes(label=sprintf('%.0f%%',Percent)),
            position=position_stack(vjust=0.5),size=3,color='white',fontface='bold') +
  geom_text(data=bin_counts,aes(x=Interval,y=104,label=paste0('n=',N)),inherit.aes=FALSE,size=3.2) +
  facet_wrap(~Population,nrow=1) +
  scale_fill_manual(values=c('Formal reversion'='#0072B2','Partial reversion; residual remains'='#E69F00',
                             'Persistent/no detected recovery'='#D55E00','Unresolved/other'='#8C8C8C')) +
  scale_y_continuous(breaks=seq(0,100,20),limits=c(0,108),expand=c(0,0)) +
  labs(x='Years since smoking cessation',y='Current-smoking-associated probe designs (%)',fill=NULL,
       title='A. Reversion classification across cessation intervals') +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=25,hjust=1),strip.background=element_blank(),
        strip.text=element_text(face='bold'),legend.position='bottom',plot.title=element_text(face='bold'))

ggsave(file.path(reversion_dir,'weighted_reversion_classification.pdf'),plot_a,width=11,height=5.8)
ggsave(file.path(reversion_dir,'weighted_reversion_classification.png'),plot_a,width=11,height=5.8,dpi=400)

############################################################################################################
## Panel B: AHRR and ECEL1P1 trajectories relative to never smokers
focus <- read.csv(file.path(reversion_dir,'weighted_AHRR_ECEL1P1_all_probe_designs.csv'),check.names=FALSE) %>%
  filter(CpG_ID %in% c('cg05575921','cg04284171','cg21566642'))

former_trajectory <- focus %>%
  transmute(Weighted_model,Focus_gene,CpG_ID,Probe_ID,Interval,Interval_order,
            Estimate=Former_vs_never_Beta,SE=Former_vs_never_SE)

current_trajectory <- focus %>%
  distinct(Weighted_model,Focus_gene,CpG_ID,Probe_ID,Current_vs_never_Beta,Current_vs_never_SE) %>%
  transmute(Weighted_model,Focus_gene,CpG_ID,Probe_ID,Interval='Current',Interval_order=0,
            Estimate=Current_vs_never_Beta,SE=Current_vs_never_SE)

trajectory <- bind_rows(current_trajectory,former_trajectory) %>%
  mutate(Interval=factor(Interval,levels=c('Current',interval_levels)),
         CpG_label=paste0(Focus_gene,'  ',CpG_ID),
         Lower=Estimate-1.96*SE,Upper=Estimate+1.96*SE)

plot_b <- ggplot(trajectory,aes(x=Interval,y=Estimate,group=Probe_ID,color=Probe_ID)) +
  geom_hline(yintercept=0,color='black',linetype='dashed',linewidth=0.4) +
  geom_errorbar(aes(ymin=Lower,ymax=Upper),width=0.08,linewidth=0.45,alpha=0.7) +
  geom_line(linewidth=0.7) + geom_point(size=2) +
  facet_grid(CpG_label~Weighted_model) +
  labs(x=NULL,y='Adjusted beta-value difference from never smokers',color='MSA probe design',
       title='B. AHRR and ECEL1P1 methylation trajectories after smoking cessation') +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=30,hjust=1),strip.background=element_blank(),
        strip.text=element_text(face='bold'),legend.position='bottom',plot.title=element_text(face='bold'))

ggsave(file.path(reversion_dir,'weighted_AHRR_ECEL1P1_trajectories.pdf'),plot_b,width=12,height=9)
ggsave(file.path(reversion_dir,'weighted_AHRR_ECEL1P1_trajectories.png'),plot_b,width=12,height=9,dpi=400)

############################################################################################################
## Combined manuscript figure
combined <- plot_grid(plot_a,plot_b,ncol=1,rel_heights=c(0.8,1.25),align='v',axis='lr')
ggsave(file.path(reversion_dir,'weighted_reversion_combined_figure.pdf'),combined,width=12,height=14)
ggsave(file.path(reversion_dir,'weighted_reversion_combined_figure.png'),combined,width=12,height=14,dpi=400)


