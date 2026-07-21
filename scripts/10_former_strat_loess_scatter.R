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
load('../output/cleaned_MSA_niehs_alhs.rdata')
beta2 = t(clean_betas) # methylation data
pheno = clean_sampsheet
tsv_path <- "../MSA.hg38.manifest.gencode.v41.tsv.gz"
manifest = fread(tsv_path)

fm1 <- fread('../output/significant_cpgs/svyglm/formersmokers/imputed/unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv') %>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1')) %>%
  dplyr::select(Probe_ID,genesUniq)
fm3 <- fread('../output/significant_cpgs/svyglm/formersmokers/imputed/unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1')) %>%
  dplyr::select(Probe_ID,genesUniq)
fm10 <- fread('../output/significant_cpgs/svyglm/formersmokers/imputed/weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1')) %>%
  dplyr::select(Probe_ID,genesUniq)
fm12 <- fread('../output/significant_cpgs/svyglm/formersmokers/imputed/weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')%>%
  filter(stringr::str_detect(genesUniq,'ECEL1P1')) %>%
  dplyr::select(Probe_ID,genesUniq)
fm <- bind_rows(list(fm1,fm3,fm10,fm12)) %>% distinct()
fm <- rbind(
  fm,
  data.table(
    Probe_ID = c(
      "cg13193840_TC11"
    ),
    genesUniq = c(
      "ECEL1P1"
    )
  )
)
probes <- fm$Probe_ID


new_pheno <- data.table::fread('../extra_vars/lhs28901_Methylation_Phenotype_2025APR29.txt',sep=',')
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
epic <- fread('../EPIC.hg38.manifest.gencode.v36.tsv.gz') %>% mutate(epicT =T) %>% dplyr::select(probeID,epicT)
fourfity<- fread('../HM450.hg38.manifest.tsv.gz')%>% mutate(ffiftyT =T) %>% dplyr::select(Probe_ID,ffiftyT)
##reduce probes to cpgs
probies <- manifest %>% filter(probeID %in% probes) %>% dplyr::select(probeID,genesUniq) %>%
  mutate(cpg = gsub(x = probeID, '\\_.*','')) %>% 
  left_join(epic,by=c('cpg'='probeID')) %>% left_join(fourfity,by=c('cpg'='Probe_ID')) %>% replace_na(list(epicT=FALSE,ffiftyT=FALSE)) %>%
  mutate(msa_only=case_when(epicT==TRUE | ffiftyT==TRUE ~ FALSE,.default=TRUE))


####################
## With avg lines
pdf("~/Desktop/PER_PROBE_small_ahrr_ecel1p1_mean_scatter_plot.pdf")
probes<-c("cg21566642_TC13" ,"cg21566642_TC11" , "cg21566642_TC12","cg13193840_TC11")
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
  print(ggplot(ahrrf, aes(x = LH_TimeSinceSmked_v1.2, y = Beta_Value)) + geom_point(aes(color=smoke_status))   + labs(color='Smoke Status') + 
          geom_smooth(method='loess', formula= y~x,color = "#374E55FF") + geom_hline(aes(yintercept=currb,linetype="Current Smokers"),color = "#B24745FF") + 
          #annotate("rect", xmin = -Inf, xmax = Inf, ymin = curr5, ymax = curr95, fill = "#B24745FF", alpha = .1, color = NA) +
          geom_hline(aes(yintercept=nevb,linetype="Never Smokers"),color = "#00A1D5FF") + #annotate("rect", xmin = -Inf, xmax = Inf, ymin = nev5, ymax = nev95, fill = "#00A1D5FF", alpha = .1, color = NA) +
          scale_linetype_manual(name = "Mean Beta", 
                                values = c(2, 2)) +
          theme(plot.title = element_text(hjust = 0.5)) + theme_classic() + ggtitle(paste0('Probe: ',p,' Gene: ',gene)))
  
}
dev.off()
beta2_subset <- beta2[, probes, drop = FALSE]

save(manifest,probes,beta2_subset,pheno,file='for_charles_time_since.rdata')
#################
## With avg lines
pdf("~/Desktop/Per_CPG_ahrr_ecel1p1_mean_scatter_plot.pdf", width = 11)

cpgs <- unique(sub("_.*", "", probes))

for (cpg in cpgs) {
  
  ## all probe columns that belong to this CpG
  matching_probes <- grep(paste0("^", cpg, "_"), colnames(beta2), value = TRUE)
  
  ## mean beta per sample across probes for this CpG
  cpg_df <- data.frame(
    Sample_ID = rownames(beta2),
    Beta_Value = rowMeans(beta2[, matching_probes, drop = FALSE], na.rm = TRUE)
  )
  
  ahrr <- left_join(cpg_df, pheno, by = "Sample_ID")
  
  ahrr <- ahrr %>%
    filter(former_5yrs_split != "Missing") %>%
    mutate(
      LH_TimeSinceSmked_v1.2 = case_when(
        LH_smoke_status == 0 ~ 67,
        LH_smoke_status == 2 ~ 0,
        .default = LH_TimeSinceSmked_v1
      ),
      smoke_status = case_when(
        LH_smoke_status == 0 ~ "Never Smoker",
        LH_smoke_status == 2 ~ "Current Smoker",
        LH_smoke_status == 1 ~ "Former Smoker"
      )
    )
  
  curr <- ahrr %>% filter(smoke_status == "Current Smoker")
  currb <- mean(curr$Beta_Value, na.rm = TRUE)
  curr95 <- quantile(curr$Beta_Value, 0.75, na.rm = TRUE)
  curr5  <- quantile(curr$Beta_Value, 0.25, na.rm = TRUE)
  
  nev <- ahrr %>% filter(smoke_status == "Never Smoker")
  nevb <- mean(nev$Beta_Value, na.rm = TRUE)
  nev95 <- quantile(nev$Beta_Value, 0.75, na.rm = TRUE)
  nev5  <- quantile(nev$Beta_Value, 0.25, na.rm = TRUE)
  
  ahrrf <- ahrr %>% filter(smoke_status == "Former Smoker")
  
  gene <- manifest %>%
    dplyr::filter(sub("_.*", "", probeID) == cpg) %>%
    dplyr::distinct(genesUniq) %>%
    dplyr::pull(genesUniq) %>%
    paste(collapse = ", ")
  
  print(
    ggplot(ahrrf, aes(x = LH_TimeSinceSmked_v1.2, y = Beta_Value)) +
      geom_point(aes(color = smoke_status)) +
      #scale_color_manual(values = c("Former Smoker" = "#DF8F44FF")) +
      labs(color = "Smoke Status") +
      geom_smooth(method = "loess", formula = y ~ x, color = "#374E55FF") +
      geom_hline(aes(yintercept = currb, linetype = "Current Smokers"), color = "#B24745FF") +
      geom_hline(aes(yintercept = nevb, linetype = "Never Smokers"), color = "#00A1D5FF") +
      #annotate("rect", xmin = -Inf, xmax = Inf, ymin = curr5, ymax = curr95,
      #         fill = "#B24745FF", alpha = 0.1, color = NA) +
      #annotate("rect", xmin = -Inf, xmax = Inf, ymin = nev5, ymax = nev95,
      #         fill = "#00A1D5FF", alpha = 0.1, color = NA) +
      scale_linetype_manual(name = "Mean Beta", values = c(2, 2)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0("CpG: ", cpg, " Gene: ", gene))
  )
}

dev.off()

















####################
## With loess degree 0 defaults to moving average
pdf("~/Desktop/ahrr_ecel1p1_moving_avg_scatter_plot.pdf",width=11)
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

pdf("~/Desktop/ahrr_ecel1p1_median_scatter_plot.pdf",width=11)
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
pdf("ahrr_ecel1p1_former_5yr_mean_beta.pdf",width=10)
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
