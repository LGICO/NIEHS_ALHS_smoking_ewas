#########################################################################################
## Annotate output

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
dir.create("output/model_output/svyglm_anno", recursive=TRUE, showWarnings=FALSE)
dir.create("output/model_output/svyglm_significant_cpgs", recursive=TRUE, showWarnings=FALSE)
dir.create("output/boxplots", recursive=TRUE, showWarnings=FALSE)
fin <- list.files("output/model_output/svyglm", pattern="\\.csv$", full.names=TRUE)
bonf<- 0.05/230325 
tsv_path <- "references/MSA.hg38.manifest.gencode.v41.tsv.gz"
manifest = fread(tsv_path)

for(f in fin){
  uws <- fread(f)
  siguws <- uws %>% arrange(P)#%>% dplyr::filter(uws$P<=bonf)
  siguw_anno <- left_join(siguws,manifest,by=c("V1"="probeID"))
  f <- tools::file_path_sans_ext(basename(f))
  write.table(siguw_anno,file=file.path("output/model_output/svyglm_anno",paste0(f,".tsv")),sep='\t',quote = F,na = "",row.names = F)
}     
fin <- list.files("output/model_output/svyglm_anno", pattern="\\.tsv$", full.names=TRUE)

## filter by sig
for(f in fin){
  uws <- fread(f,header = T)
  print(f)
  colnames(uws)<- c('Probe_ID','Beta','SE','P','N','N.exp','N.unexp','CpG_chrm','CpG_beg','CpG_end','probe_strand','genesUniq','geneNames','transcriptTypes','transcriptIDs','distToTSS')
  siguws <- uws %>% dplyr::filter(uws$P<=bonf) %>% arrange(P) 
  f <- basename(f)
  write.table(siguws,file=file.path("output/model_output/svyglm_significant_cpgs",f),sep='\t',quote = F,na = "",row.names = F)
}    


##############################
## Weighted plots
###########
siguws <- fread('output/model_output/svyglm_significant_cpgs/weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')
siguws %>% arrange(P)
probes<-siguws %>% arrange(P) %>% dplyr::select(V1)
pdf("output/boxplots/weighted_significant_smoking_probes_europc_only.pdf")
for(p in as.list(probes$V1)){
  ahrr <- data.frame(beta_pc[,p])
  matching_rows <- grep(p, rownames(t(beta_pc)), value = TRUE)
  cpg_data <- t(beta_pc)[matching_rows, , drop = FALSE]
  cpg_long <- pivot_longer(as.data.frame(cpg_data),
                           cols = everything(),
                           names_to = "Sample_ID",
                           values_to = "Beta_Value")
  unique_samples <- unique(cpg_long$Sample_ID)
  ahrr<-left_join(cpg_long,pheno_pc,by="Sample_ID")
  
  ahrr=ahrr%>%mutate(LH_smoke_status=case_when(LH_smoke_status==0~"Never",
                                               LH_smoke_status==1~"Former",
                                               LH_smoke_status==2~"Current",
                                               LH_smoke_status==99~"Missing"))
  ahrr$LH_smoke_status=as.factor(ahrr$LH_smoke_status)
  ahrr = filter(ahrr,LH_smoke_status!="Missing")
  gene <- ungroup(siguws) %>% dplyr::select(V1,annot.symbol) %>% filter(V1 == p) %>%
    dplyr::select(annot.symbol) %>% separate(annot.symbol,into = 'x',extra='drop')
  #gene="AHRR;PDCD6-AHRR"
  print(ggboxplot(ahrr, x = "LH_smoke_status", y = "Beta_Value",
                  color = "LH_smoke_status", palette ="jco",order=c('Never','Former','Current','Missing'),
                  add = "jitter", title=paste0("Probe: ",p," Gene: ",gene," by Smoking Status") ))
  
}
dev.off()

siguws <- fread('output/model_output/svyglm_significant_cpgs/weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv')
siguws %>% arrange(P)
probes<-siguws %>% arrange(P) %>% dplyr::select(V1)
pdf("output/boxplots/weighted_significant_smoking_probes_transancestry.pdf")
for(p in as.list(probes$V1)){
  ahrr <- data.frame(beta2[,p])
  matching_rows <- grep(p, rownames(t(beta2)), value = TRUE)
  cpg_data <- t(beta2)[matching_rows, , drop = FALSE]
  cpg_long <- pivot_longer(as.data.frame(cpg_data),
                           cols = everything(),
                           names_to = "Sample_ID",
                           values_to = "Beta_Value")
  unique_samples <- unique(cpg_long$Sample_ID)
  ahrr<-left_join(cpg_long,pheno,by="Sample_ID")
  
  ahrr=ahrr%>%mutate(LH_smoke_status=case_when(LH_smoke_status==0~"Never",
                                               LH_smoke_status==1~"Former",
                                               LH_smoke_status==2~"Current",
                                               LH_smoke_status==99~"Missing"))
  ahrr$LH_smoke_status=as.factor(ahrr$LH_smoke_status)
  ahrr = filter(ahrr,LH_smoke_status!="Missing")
  gene <- ungroup(siguws) %>% dplyr::select(V1,annot.symbol) %>% filter(V1 == p) %>%
    dplyr::select(annot.symbol) %>% separate(annot.symbol,into = 'x',extra='drop')
  print(ggboxplot(ahrr, x = "LH_smoke_status", y = "Beta_Value",
                  color = "LH_smoke_status", palette ="jco",order=c('Never','Former','Current','Missing'),
                  add = "jitter", title=paste0("Probe: ",p," Gene: ",gene," by Smoking Status") ))
  
}
dev.off()
