###############################################################################################################
## Plot AHRR and ECEL1P1
## We need to get the AHRR and ECEL1P1 probes that were significant in former smokers in the original models
## Subset beta df to former smokers and add new groups
## First we can generate boxplots by groups to see if there is a trend
## Next we need to generate CB's plot with x-axis being years since 
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)

############################################################################################################
## Current analysis paths
analysis_root <- ""
significant_dir <- file.path(analysis_root,'output','significant_models')
figure_dir <- file.path(analysis_root,'output','former_stratified','figures')
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)

load('751_clean_pc.rdata')
manifest <- fread(file.path(analysis_root,'references','MSA.hg38.manifest.gencode.v41.tsv.gz'))

## Unimputed, QC-masked methylation values: samples x probe designs
beta2 <- t(betas_miss_masked)
stopifnot(all(pheno$Sample_ID %in% rownames(beta2)))
beta2 <- beta2[match(pheno$Sample_ID,rownames(beta2)),,drop=FALSE]
stopifnot(identical(rownames(beta2),pheno$Sample_ID))

############################################################################################################
## AHRR and ECEL1P1 probe designs significant in any unimputed former-smoker model
former_files <- c(
  'unimputed_unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv',
  'unimputed_unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv',
  'unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv',
  'unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'
)

probe_results <- bind_rows(lapply(former_files,function(x) fread(file.path(significant_dir,x)))) %>%
  filter(grepl('AHRR|ECEL1P1',genesUniq,ignore.case=TRUE)) %>%
  dplyr::select(Probe_ID,genesUniq) %>% distinct()

probes <- intersect(probe_results$Probe_ID,colnames(beta2))
if (length(probes)==0) stop('No significant AHRR or ECEL1P1 probe designs matched the unimputed methylation matrix.')
write.csv(probe_results %>% filter(Probe_ID %in% probes),
          file.path(figure_dir,'AHRR_ECEL1P1_unimputed_probes_plotted.csv'),row.names=FALSE)

############################################################################################################
## Add continuous time since cessation if it is not already present in pheno
if (!'LH_TimeSinceSmked_v1' %in% names(pheno)) {
  new_pheno <- fread(file.path(analysis_root,'samplesheets','lhs28901_Methylation_Phenotype_2025APR29.txt'),sep=',')
  new_pheno <- new_pheno %>% filter(LHID %in% pheno$LHID,COMMENT!='Replicate') %>%
    dplyr::select(LHID,partid,LH_TimeSinceSmked_v1)
  pheno <- left_join(pheno,new_pheno,by=c('LHID','partid'))
}

## Convert smoking status only if it is still stored as text.
if (any(as.character(pheno$LH_smoke_status) %in% c('Never','Former','Current','M'))) {
  pheno <- pheno %>% mutate(LH_smoke_status=case_when(
    LH_smoke_status=='Never' ~ 0,
    LH_smoke_status=='Former' ~ 1,
    LH_smoke_status=='Current' ~ 2,
    LH_smoke_status=='M' ~ 99))
}
pheno$LH_smoke_status <- as.character(pheno$LH_smoke_status)

pheno <- pheno %>% mutate(
  smoking_group=case_when(
    LH_smoke_status=='0' ~ 'Never smoker',
    LH_smoke_status=='1' ~ 'Former smoker',
    LH_smoke_status=='2' ~ 'Current smoker',
    TRUE ~ 'Missing/unknown'),
  cessation_group=case_when(
    smoking_group=='Current smoker' ~ 'Current',
    smoking_group=='Former smoker' & LH_TimeSinceSmked_v1<=10 ~ 'Former <=10 years',
    smoking_group=='Former smoker' & LH_TimeSinceSmked_v1>10 & LH_TimeSinceSmked_v1<=20 ~ 'Former 11-20 years',
    smoking_group=='Former smoker' & LH_TimeSinceSmked_v1>20 & LH_TimeSinceSmked_v1<=30 ~ 'Former 21-30 years',
    smoking_group=='Former smoker' & LH_TimeSinceSmked_v1>30 ~ 'Former >30 years',
    smoking_group=='Never smoker' ~ 'Never'))

group_levels <- c('Current','Former <=10 years','Former 11-20 years','Former 21-30 years','Former >30 years','Never')
pheno$cessation_group <- factor(pheno$cessation_group,levels=group_levels)

############################################################################################################
## Figure 1: one page per unimputed MSA probe design
pdf(file.path(figure_dir,'AHRR_ECEL1P1_unimputed_per_probe_continuous.pdf'),width=8,height=6,onefile=TRUE)
for (p in probes) {
  plot_data <- data.frame(Sample_ID=rownames(beta2),Beta_Value=beta2[,p]) %>%
    left_join(pheno,by='Sample_ID')
  
  current_mean <- plot_data %>% filter(smoking_group=='Current smoker') %>% summarise(x=mean(Beta_Value,na.rm=TRUE)) %>% pull(x)
  never_mean <- plot_data %>% filter(smoking_group=='Never smoker') %>% summarise(x=mean(Beta_Value,na.rm=TRUE)) %>% pull(x)
  former_data <- plot_data %>% filter(smoking_group=='Former smoker',!is.na(LH_TimeSinceSmked_v1),LH_TimeSinceSmked_v1>=0)
  gene <- manifest %>% filter(probeID==p) %>% distinct(genesUniq) %>% pull(genesUniq) %>% paste(collapse=', ')
  
  g <- ggplot(former_data,aes(x=LH_TimeSinceSmked_v1,y=Beta_Value)) +
    geom_point(color='#F8766D',alpha=0.55,size=1.5) +
    geom_smooth(method='loess',formula=y~x,color='#0072B2',fill='grey75',linewidth=0.8) +
    geom_hline(aes(yintercept=current_mean,linetype='Current-smoker mean'),color='#FF0000',linewidth=0.7) +
    geom_hline(aes(yintercept=never_mean,linetype='Never-smoker mean'),color='#0000FF',linewidth=0.7) +
    scale_linetype_manual(values=c('Current-smoker mean'='dashed','Never-smoker mean'='dotted'),name=NULL) +
    coord_cartesian(ylim=c(0,1)) +
    labs(title=paste0(gene,': ',p),x='Years since smoking cessation',y='Methylation beta value') +
    theme_classic(base_size=12) + theme(plot.title=element_text(hjust=0.5,face='bold'),legend.position='bottom')
  print(g)
}
dev.off()

############################################################################################################
## Figure 2: average methylation across all MSA designs for each base CpG
cpgs <- unique(sub('_.*$','',probes))
pdf(file.path(figure_dir,'AHRR_ECEL1P1_unimputed_per_cpg_continuous.pdf'),width=8,height=6,onefile=TRUE)
for (cpg in cpgs) {
  matching_probes <- grep(paste0('^',cpg,'_'),colnames(beta2),value=TRUE)
  cpg_mean <- rowMeans(beta2[,matching_probes,drop=FALSE],na.rm=TRUE)
  cpg_mean[is.nan(cpg_mean)] <- NA_real_
  plot_data <- data.frame(Sample_ID=rownames(beta2),Mean_Beta_Value=cpg_mean) %>% left_join(pheno,by='Sample_ID')
  
  current_mean <- plot_data %>% filter(smoking_group=='Current smoker') %>% summarise(x=mean(Mean_Beta_Value,na.rm=TRUE)) %>% pull(x)
  never_mean <- plot_data %>% filter(smoking_group=='Never smoker') %>% summarise(x=mean(Mean_Beta_Value,na.rm=TRUE)) %>% pull(x)
  former_data <- plot_data %>% filter(smoking_group=='Former smoker',!is.na(LH_TimeSinceSmked_v1),LH_TimeSinceSmked_v1>=0)
  gene <- manifest %>% filter(sub('_.*$','',probeID)==cpg) %>% distinct(genesUniq) %>% pull(genesUniq) %>% paste(collapse=', ')
  
  g <- ggplot(former_data,aes(x=LH_TimeSinceSmked_v1,y=Mean_Beta_Value)) +
    geom_point(color='#F8766D',alpha=0.55,size=1.5) +
    geom_smooth(method='loess',formula=y~x,color='#0072B2',fill='grey75',linewidth=0.8) +
    geom_hline(aes(yintercept=current_mean,linetype='Current-smoker mean'),color='#FF0000',linewidth=0.7) +
    geom_hline(aes(yintercept=never_mean,linetype='Never-smoker mean'),color='#0000FF',linewidth=0.7) +
    scale_linetype_manual(values=c('Current-smoker mean'='dashed','Never-smoker mean'='dotted'),name=NULL) +
    coord_cartesian(ylim=c(0,1)) +
    labs(title=paste0(gene,': ',cpg,' (mean of ',length(matching_probes),' MSA designs)'),
         x='Years since smoking cessation',y='Mean methylation beta value across MSA probe designs') +
    theme_classic(base_size=12) + theme(plot.title=element_text(hjust=0.5,face='bold'),legend.position='bottom')
  print(g)
}
dev.off()

############################################################################################################
## Figure 3: observed mean beta value by the four reversion-analysis intervals
group_plot_data <- bind_rows(lapply(cpgs,function(cpg) {
  matching_probes <- grep(paste0('^',cpg,'_'),colnames(beta2),value=TRUE)
  cpg_mean <- rowMeans(beta2[,matching_probes,drop=FALSE],na.rm=TRUE)
  cpg_mean[is.nan(cpg_mean)] <- NA_real_
  gene <- manifest %>% filter(sub('_.*$','',probeID)==cpg) %>% distinct(genesUniq) %>% pull(genesUniq) %>% paste(collapse=', ')
  data.frame(Sample_ID=rownames(beta2),CpG_ID=cpg,Gene=gene,Beta_Value=cpg_mean) %>% left_join(pheno,by='Sample_ID')
})) %>% filter(!is.na(cessation_group)) %>%
  group_by(Gene,CpG_ID,cessation_group) %>%
  summarise(N=sum(!is.na(Beta_Value)),Mean_Beta=mean(Beta_Value,na.rm=TRUE),
            SE=sd(Beta_Value,na.rm=TRUE)/sqrt(N),.groups='drop')

plot_groups <- ggplot(group_plot_data,aes(x=cessation_group,y=Mean_Beta,group=1)) +
  geom_line(color='black',linewidth=0.6) + geom_point(color='#0072B2',size=2.2) +
  geom_errorbar(aes(ymin=Mean_Beta-1.96*SE,ymax=Mean_Beta+1.96*SE),width=0.12,color='#0072B2') +
  facet_wrap(~paste0(Gene,': ',CpG_ID),scales='free_y') +
  labs(x='Smoking group and years since cessation',y='Observed mean methylation beta value (95% CI)') +
  theme_classic(base_size=11) + theme(axis.text.x=element_text(angle=30,hjust=1),strip.background=element_blank(),strip.text=element_text(face='bold'))

ggsave(file.path(figure_dir,'AHRR_ECEL1P1_unimputed_group_means.pdf'),plot_groups,width=12,height=8)
ggsave(file.path(figure_dir,'AHRR_ECEL1P1_unimputed_group_means.png'),plot_groups,width=12,height=8,dpi=400)

#save(manifest,probes,beta2,pheno,file=file.path(figure_dir,'AHRR_ECEL1P1_unimputed_plot_data.rdata'))
