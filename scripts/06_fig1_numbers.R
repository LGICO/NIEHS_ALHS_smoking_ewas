library(data.table)
library(dplyr)

load('output/cleaned_MSA_niehs_alhs.rdata')
beta2 = t(clean_betas) # methylation data
betap <- t(clean_betas_pval_masked)
pheno = clean_sampsheet# phenotype data
setwd('output/model_output/')

table(pheno$European)
dim(beta2)

model1 = fread('svyglm_significant_cpgs/unimputed_unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
model3 = fread('svyglm_significant_cpgs/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
model10 = fread('svyglm_significant_cpgs/unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')
model12 = fread('svyglm_significant_cpgs/unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')


m1p = model1 %>% select(Probe_ID)
m3p = model3 %>% select(Probe_ID)
dplyr::full_join(m1p,m3p) %>% distinct

m10p = model10 %>% select(Probe_ID)
m12p =  model12 %>% select(Probe_ID)
full_join(m10p,m12p) %>% distinct

sp <- full_join(m1p,m3p) %>% full_join(.,m10p) %>% full_join(.,m12p) %>% distinct

###########
model1 = fread('svyglm_significant_cpgs/unimputed_unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
model3 = fread('svyglm_significant_cpgs/unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
model10 = fread('svyglm_significant_cpgs/unimputed_unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')
model12 = fread('svyglm_significant_cpgs/unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')


mf1p = model1 %>% select(Probe_ID)
mf3p = model3 %>% select(Probe_ID)
dplyr::full_join(mf1p,mf3p) %>% distinct

mf10p = model10 %>% select(Probe_ID)
mf12p =  model12 %>% select(Probe_ID)
full_join(mf10p,mf12p) %>% distinct

fp <- full_join(mf1p,mf3p) %>% full_join(.,mf10p) %>% full_join(.,mf12p) %>% distinct

msf <- full_join(m1p,m3p) %>% full_join(.,m10p) %>% full_join(.,m12p) %>% full_join(.,mf1p) %>% full_join(.,mf3p) %>% full_join(.,mf10p) %>% full_join(.,mf12p) %>% distinct

bonf<- 0.05/230325
bonf_epic = 0.05/817235
chrs<-c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
        "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
path = ""
mani_epic <- fread(paste0("references/EPIC.hg38.manifest.gencode.v36.tsv.gz")) %>%
  dplyr::select(probeID,genesUniq) %>% mutate(CpG=gsub("\\_.*","",probeID)) %>% dplyr::select(CpG,genesUniq) %>% distinct()

## Model 1 ###########################################
smokers=fread(file=paste0('svyglm_anno/unimputed_unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("ProbeID","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
formers=fread(file=paste0('svyglm_anno/unimputed_unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')) %>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("ProbeID","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
model1 <- bind_rows(smoker = smokers, former_smoker = formers,.id="status") %>% 
  mutate(CpG=gsub("\\_.*","",ProbeID))%>%
  filter(pvalue<=bonf) %>%
  dplyr::mutate(EPIC = case_when(CpG %in% mani_epic$CpG ~ TRUE, .default = FALSE)) %>%
  arrange(pvalue)
model1 %>% filter(EPIC== TRUE) %>% select(ProbeID) %>% distinct%>% dim
m1p <- model1 %>% filter(EPIC== TRUE) %>% select(ProbeID) %>% distinct
## Model 3 ###########################################
smokers=fread(file=paste0('svyglm_anno/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("ProbeID","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
formers=fread(file=paste0('svyglm_anno/unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')) %>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("ProbeID","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
model3 <- bind_rows(smoker = smokers, former_smoker = formers,.id="status") %>% 
  mutate(CpG=gsub("\\_.*","",ProbeID))%>%
  filter(pvalue<=bonf) %>%
  dplyr::mutate(EPIC = case_when(CpG %in% mani_epic$CpG ~ TRUE, .default = FALSE)) %>%
  arrange(pvalue)
model3 %>% filter(EPIC== TRUE)%>% select(ProbeID) %>% distinct %>% dim
m3p <- model3 %>% filter(EPIC== TRUE)%>% select(ProbeID) %>% distinct
## Model 10 ###########################################
smokers=fread(file=paste0('svyglm_anno/unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("ProbeID","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
formers=fread(file=paste0('svyglm_anno/unimputed_unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')) %>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("ProbeID","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
model10 <- bind_rows(smoker = smokers, former_smoker = formers,.id="status") %>% 
  mutate(CpG=gsub("\\_.*","",ProbeID))%>%
  filter(pvalue<=bonf) %>%
  dplyr::mutate(EPIC = case_when(CpG %in% mani_epic$CpG ~ TRUE, .default = FALSE)) %>%
  arrange(pvalue)
model10 %>% filter(EPIC== TRUE) %>% select(ProbeID) %>% distinct %>% dim
m10p <- model10 %>% filter(EPIC== TRUE) %>% select(ProbeID) %>% distinct
## Model 12 ###########################################
smokers=fread(file=paste0('svyglm_anno/unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("ProbeID","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
formers=fread(file=paste0('svyglm_anno/unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')) %>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("ProbeID","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
model12 <- bind_rows(smoker = smokers, former_smoker = formers,.id="status") %>% 
  mutate(CpG=gsub("\\_.*","",ProbeID))%>%
  filter(pvalue<=bonf) %>%
  dplyr::mutate(EPIC = case_when(CpG %in% mani_epic$CpG ~ TRUE, .default = FALSE)) %>%
  arrange(pvalue)
model12 %>% filter(EPIC== TRUE) %>% select(ProbeID) %>% distinct %>% dim
m12p <- model12 %>% filter(EPIC== TRUE) %>% select(ProbeID) %>% distinct

##81 -- 55
##77 -- 52
##110 --184
##91 -- 61
