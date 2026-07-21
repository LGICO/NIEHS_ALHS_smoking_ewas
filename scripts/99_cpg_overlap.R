tsv_path <- "references/MSA.hg38.manifest.tsv.gz"
path = ""
bonf<- 0.05/230325
bonf_epic = 0.05/817235
#sesameDataCache()
manifest = fread(tsv_path)
mani = fread(paste0(path,"MSA.hg38.manifest.gencode.v41.tsv.gz")) %>% 
  dplyr::select(probeID,genesUniq)
mani_epic <- fread(file=paste0(path,"EPIC.hg38.manifest.gencode.v36.tsv.gz")) %>%
  dplyr::select(probeID,genesUniq) %>% mutate(CpG=gsub("\\_.*","",probeID)) %>% dplyr::select(CpG,genesUniq) %>% distinct()

chrs<-c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
        "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
ma<-mani %>% filter(str_detect(genesUniq, regex("\\bAHRR", ignore_case = TRUE)))%>% mutate(probeID=gsub("\\_.*","",probeID))
mae<-mani_epic %>% filter(str_detect(genesUniq, regex("\\bAHRR", ignore_case = TRUE)))
api <- ma %>% filter(!probeID %in% mae$CpG) %>% select(probeID)

lond_sig <- fread(paste0(path,"/output/EWAS_ALHS_EU_ALL_smok1vs3_e.cov11_817235x1619.csv.s"),sep = ",") 
names(lond_sig)<-c("CpG","Beta","SE","P","N","N.exp","N.unexp")
lond_sig <- left_join(lond_sig,mani_epic,by='CpG')
lond_sig1 <- lond_sig %>% filter(pvalue<=bonf_epic)
x <- lond_sig %>% select("CpG", "Beta", "SE", "P", "N", "N.exp", "N.unexp", "CpG_chrm", 
                      "CpG_beg", "CpG_end", "probe_strand","genesUniq", 
                      "geneNames", "transcriptTypes", "transcriptIDs", "distToTSS")

lond_sig <- fread(paste0(path,"/output/EWAS_ALHS_EU_ALL_smok1vs2_e.cov11_2190_2025-04-23.csv"),sep = ",") 
names(lond_sig)<-c("CpG","Beta","SE","P","N","N.exp","N.unexp")
lond_sig <- left_join(lond_sig,mani_epic,by='CpG')

x <- lond_sig %>% select("CpG", "Beta", "SE", "P", "N", "N.exp", "N.unexp", "CpG_chrm", 
                         "CpG_beg", "CpG_end", "probe_strand","genesUniq", 
                         "geneNames", "transcriptTypes", "transcriptIDs", "distToTSS")
## Model 1 ###########################################
smokers=fread(file=paste0(path,'output/svyglm/unimputed_unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("ProbeID","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
formers=fread(file=paste0(path,'output/svyglm/unimputed_unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')) %>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("ProbeID","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
model1 <- bind_rows(smoker = smokers, former_smoker = formers,.id="status") %>% 
  mutate(CpG=gsub("\\_.*","",ProbeID))%>%
  filter(pvalue<=bonf) %>%
  dplyr::mutate(EPIC = case_when(CpG %in% mani_epic$CpG ~ TRUE, .default = FALSE),
                London_EWAS = case_when(CpG %in% lond_sig1$CpG ~ TRUE, .default = FALSE)) %>%
  arrange(pvalue)
  

## Model 3 ###########################################
smokers=fread(file=paste0(path,'output/svyglm/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("ProbeID","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
formers=fread(file=paste0(path,'output/svyglm/unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')) %>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("ProbeID","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
model3 <- bind_rows(smoker = smokers, former_smoker = formers,.id="status") %>% 
  mutate(CpG=gsub("\\_.*","",ProbeID))%>%
  filter(pvalue<=bonf) %>%
  dplyr::mutate(EPIC = case_when(CpG %in% mani_epic$CpG ~ TRUE, .default = FALSE),
                London_EWAS = case_when(CpG %in% lond_sig1$CpG ~ TRUE, .default = FALSE)) %>%
  arrange(pvalue)


## Model 10 ###########################################
smokers=fread(file=paste0(path,'output/svyglm/unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("ProbeID","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
formers=fread(file=paste0(path,'output/svyglm/unimputed_unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')) %>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("ProbeID","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
model10 <- bind_rows(smoker = smokers, former_smoker = formers,.id="status") %>% 
  mutate(CpG=gsub("\\_.*","",ProbeID))%>%
  filter(pvalue<=bonf) %>%
  dplyr::mutate(EPIC = case_when(CpG %in% mani_epic$CpG ~ TRUE, .default = FALSE),
                London_EWAS = case_when(CpG %in% lond_sig1$CpG ~ TRUE, .default = FALSE)) %>%
  arrange(pvalue)


## Model 12 ###########################################
smokers=fread(file=paste0(path,'output/svyglm/unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("ProbeID","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
formers=fread(file=paste0(path,'output/svyglm/unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')) %>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("ProbeID","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
model12 <- bind_rows(smoker = smokers, former_smoker = formers,.id="status") %>% 
  mutate(CpG=gsub("\\_.*","",ProbeID))%>%
  filter(pvalue<=bonf) %>%
  dplyr::mutate(EPIC = case_when(CpG %in% mani_epic$CpG ~ TRUE, .default = FALSE),
                London_EWAS = case_when(CpG %in% lond_sig1$CpG ~ TRUE, .default = FALSE)) %>%
  arrange(pvalue)

save(model1,model3,model10,model12,file = paste0(path,'output/sig_cpg_msa_epic_overlap.rdata'))




