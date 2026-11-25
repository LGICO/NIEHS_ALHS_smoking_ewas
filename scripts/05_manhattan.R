library(data.table)
library(sesame)
library(dplyr)
library(gtools)
library(ggplot2)
library(ggrepel)
library(miamiplot)
library(stringr)
library(gdata)

dir.create("output/miami", recursive=TRUE, showWarnings=FALSE)

tsv_path <- "references/MSA.hg38.manifest.tsv.gz"
path = ""
## Run from the repository parent directory.
bonf<- 0.05/230325
bonf_epic = 0.05/817235
#sesameDataCache()
manifest = fread(tsv_path)
mani = fread(paste0(path,"references/MSA.hg38.manifest.gencode.v41.tsv.gz")) %>% 
  dplyr::select(probeID,genesUniq)
mani_epic <- fread(paste0(path,"references/EPIC.hg38.manifest.gencode.v36.tsv.gz")) %>%#EPICv2.hg38.manifest.gencode.v41.tsv.gz
  dplyr::select(probeID,genesUniq) %>% mutate(CpG=gsub("\\_.*","",probeID)) %>% dplyr::select(CpG,genesUniq) %>% distinct()
mani_epic2 <- fread(paste0(path,"references/EPICv2.hg38.manifest.gencode.v41.tsv.gz")) %>%#EPICv2.hg38.manifest.gencode.v41.tsv.gz
  dplyr::select(probeID,genesUniq) %>% mutate(CpG=gsub("\\_.*","",probeID)) %>% dplyr::select(CpG,genesUniq) %>% distinct()
mani_450 <- fread(paste0(path,"references/HM450.hg38.manifest.tsv.gz")) %>%
   mutate(CpG=gsub("\\_.*","",Probe_ID)) %>% dplyr::select(CpG) %>% distinct()
chrs<-c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
        "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
ma<-mani %>% filter(str_detect(genesUniq, regex("\\bAHRR", ignore_case = TRUE)))%>% mutate(probeID=gsub("\\_.*","",probeID))
mae<-mani_epic %>% filter(str_detect(genesUniq, regex("\\bAHRR", ignore_case = TRUE))) %>% dplyr::select(CpG)
cpge4 <- rbind(mani_epic,mani_450,fill=T) %>%dplyr::select(CpG) %>% distinct()
only450<-anti_join(mani_450,mani_epic,by='CpG')
api <- ma %>% filter(!probeID %in% mae$CpG) %>% dplyr::select(probeID)

lond_sig <- fread("output/EWAS_ALHS_EU_ALL_smok1vs3_e.cov11_817235x1619.csv.s",sep = ",") %>%
  dplyr::select(V1,P)
names(lond_sig)<-c("CpG","pvalue")
lond_sig <- left_join(lond_sig,mani_epic,by='CpG')
lond_sig1 <- lond_sig %>% filter(pvalue<=bonf_epic)





wes <- fread('output/svyglm_anno/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv') %>% 
  mutate(Z = Beta/SE,
         chisq = Z^2)
wef <- fread('output/svyglm_anno/unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv') %>% 
  mutate(Z = Beta/SE,
         chisq = Z^2)

wts <- fread('output/svyglm_anno/unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv') %>% 
  mutate(Z = Beta/SE,
         chisq = Z^2)
wtf <- fread('output/svyglm_anno/unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv') %>% 
  mutate(Z = Beta/SE,
         chisq = Z^2)
known_cpgs <- unique(c(mani_epic$CpG, mani_epic2$CpG, mani_450$CpG))
msa_specific <- bind_rows(list("Smokers, Euro"=wes, "Former, Euro"=wef, "Smokers, Trans"=wts, "Former, Trans"=wtf), .id="Model") %>%
  filter(P <= bonf) %>% mutate(CpG_ID=sub("_.*$", "", V1)) %>%
  filter(!CpG_ID %in% known_cpgs) %>%  
  group_by(Model, CpG_ID) %>%
  slice_min(P, n=1, with_ties=FALSE) %>%
  ungroup()
bind_rows(msa_specific %>% count(Model, name="MSA_specific_CpGs"),
          msa_specific %>% distinct(CpG_ID) %>% summarise(Model="Unique across all models", MSA_specific_CpGs=n()))
write.table(msa_specific,file="output/miami/msa_specific_probes.tsv",sep='\t',na='',row.names = F,quote = F)
## 34 unique cpgs across the 4 models 75 with overlap across the 4 models

###################
#################
#####################
plot_miami<- function(smoke, form, path,output_pdf,mani_epic,bonf) {
    plot_df <- bind_rows(smoke,form) %>% 
      group_by(chr) %>% 
      # Compute chromosome size
      summarise(chrlength = max(pos)) %>%  
      # Calculate cumulative position of each chromosome
      mutate(cumulativechrlength = cumsum(as.numeric(chrlength))-chrlength) %>% 
      dplyr::select(-chrlength) %>%
      # Temporarily add the cumulative length of each chromosome to the initial 
      # dataset 
      left_join(bind_rows(top = smokers,bottom = formers,.id="status"), ., by=c("chr"="chr")) %>%
      # Sort by chr then position 
      arrange(chr, pos) %>%
      # Add the position to the cumulative chromosome length to get the position of 
      # this probe relative to all other probes
      mutate(rel_pos = pos + cumulativechrlength) %>%
      # Calculate the logged p-value too
      mutate(logged_p = -log10(pvalue)) %>%
      dplyr::select(-cumulativechrlength)
    
    axis_df <- plot_df %>% 
      group_by(chr) %>% 
      summarize(chr_center=(max(rel_pos) + min(rel_pos)) / 2)
    maxp <- ceiling(max(plot_df$logged_p, na.rm = TRUE))
    
    plot_df <- plot_df %>% mutate(CpGe=gsub("\\_.*","",CpG))
    
    upper_labels<- plot_df %>%
      filter(status == "top") %>%
      filter(!CpGe %in% lond_sig1$CpG) %>%
      filter(!CpGe %in% only450$CpG) %>%
      ##filter(!gene %in% lond_sig1$genesUniq) %>%
      filter(pvalue<=bonf) %>%
      group_by(chr) %>%
      arrange(desc(logged_p)) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      arrange(desc(logged_p)) %>%
      dplyr::slice(1:8) %>%
      dplyr::select(CpG)
    upper_labels2 <- plot_df %>%
      filter(status=='top') %>%
      filter(!CpGe %in% cpge4$CpG) %>%
      filter(pvalue<=1e-20)
    ahrr<- plot_df %>%
      filter(status=='top') %>%
      filter(CpGe %in% api$probeID) %>%
      filter(pvalue<=bonf)
    upper_labels2<- bind_rows(upper_labels2,ahrr)
    lower_labels <- plot_df %>%
      filter(status == "bottom") %>%
      filter(!CpGe %in% lond_sig1$CpG)%>%
      filter(!CpGe %in% only450$CpG) %>%
      #filter(!gene %in% lond_sig1$genesUniq) %>%
      filter(pvalue<=bonf) %>%
      group_by(chr) %>%
      arrange(desc(logged_p)) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      arrange(desc(logged_p)) %>%
      dplyr::slice(1:8) %>%
      dplyr::select(CpG)
    lower_labels2 <- plot_df %>%
      filter(status=='bottom') %>%
      filter(!CpGe %in% cpge4$CpG) %>%
      filter(pvalue<=1e-20)
    
    label_list <- rbind(upper_labels, lower_labels)
    
    plot_df <- plot_df %>%
      mutate(label = case_when(CpG %in% label_list$CpG ~ 
                                 paste0(CpG, "\n", gene), TRUE ~ ""))
    upper_labels<- bind_rows(upper_labels,upper_labels2)
    lower_labels<- bind_rows(lower_labels,lower_labels2)
    ahrr<- plot_df %>%
      filter(status=='bottom') %>%
      filter(CpGe %in% api$probeID) %>%
      filter(pvalue<=bonf)
    lower_labels2<- bind_rows(lower_labels2,ahrr)
    label_list <- rbind(upper_labels, lower_labels)
    plot_df <- plot_df %>%
      mutate(label2 = case_when(CpG %in% label_list$CpG ~ 
                                  paste0(CpG, "\n", gene), TRUE ~ "")) %>% 
      mutate(label2 = case_when(label != "" ~ label,.default = label2))
    upper_plot <- ggplot() + 
      geom_point(data = plot_df[which(plot_df$status == "top" ),], 
                 aes(x = rel_pos, y = logged_p, color = as.factor(chr)), 
                 size = 1) +
      scale_color_manual(values = rep(c("#3261cf", "#88a1db"), nrow(axis_df))) +
      scale_x_continuous(labels = axis_df$chr, 
                         breaks = axis_df$chr_center, 
                         expand = expansion(mult = 0.01)) +
      scale_y_continuous(limits = c(0, maxp+5), 
                         expand = expansion(mult = c(0.02, 0))) + 
      geom_hline(yintercept = -log10(1e-5), color = "grey", linetype = "dashed", 
                 size = 0.5) +
      geom_hline(yintercept = -log10(2.170808e-07), color = "#150c6b", linetype = "solid", 
                 size = 0.5) +
      labs(x = "", y = bquote(atop('Current smoker vs. Never smoker', '-log'[10]*'(p)'))) + 
      theme_classic() +
      theme(legend.position = "none", 
            axis.title.x = element_blank(),
            plot.margin = margin(b = 0, l = 10),
            axis.text=element_text(size=16),
            axis.title=element_text(size=24,face="bold")) + 
      geom_text_repel(data = plot_df%>% filter(plot_df$status == "top",CpG %in% upper_labels$CpG ),
                       aes(x = rel_pos, y = logged_p, label = label2),
                       size = 6, segment.size = 0.2, point.padding = 0.3,
                       ylim = c(maxp / 3, NA), box.padding = 2,force=4,force_pull=0.1) 
    
    #upper_plot
    
    lower_plot <- ggplot() + 
      geom_point(data = plot_df[which(plot_df$status == "bottom" ),], 
                 aes(x = rel_pos, y = logged_p, color = as.factor(chr)), 
                 size = 1) +
      scale_color_manual(values = rep(c("#a33636", "#d99898"), nrow(axis_df))) +
      scale_x_continuous(breaks = axis_df$chr_center, 
                         position = "top",
                         expand = expansion(mult = 0.01)) +
      scale_y_reverse(limits = c(maxp+5, 0), 
                      expand = expansion(mult = c(0, 0.02))) + 
      geom_hline(yintercept = -log10(1e-5), color = "grey", linetype = "dashed", 
                 size = 0.5) +
      geom_hline(yintercept = -log10(2.170808e-07), color = "#150c6b", linetype = "solid", 
                 size = 0.5) +
      labs(x = "", y = bquote(atop('Former smoker vs. Never smoker', '-log'[10]*'(p)'))) + 
      theme_classic() +
      theme(legend.position = "none", 
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = margin(t = 0, l = 10),
            axis.text=element_text(size=16),
            axis.title=element_text(size=24,face="bold")) + 
      geom_text_repel(data = 
                         plot_df%>% filter(plot_df$status == "bottom",CpG %in% lower_labels$CpG ),
                       aes(x = rel_pos, y = logged_p, label = label2),
                       size = 6, segment.size = 0.2, point.padding = 0.3,
                       ylim = c(NA, -(maxp / 3)), box.padding = 2,force=4,force_pull=0.1)
    
    #lower_plot
    pdf(file = paste0(path,output_pdf),height = 18,width = 20)
    p <- gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
    
    print(p)
    dev.off()
}


###########################################################################################################################
## Model 1
smokers=fread(file=paste0('output/svyglm_anno/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("CpG","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
CpGs <- smokers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

#smokers$chromosome <- factor(smokers$chromosome, chrs)

formers=fread(file=paste0('output/svyglm_anno/unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')) %>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("CpG","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
CpGf <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

op<-'msa_no_overlap_london_model1_unweighted_euro_miami.pdf'
plot_miami(smoke = smokers,form = formers,path = 'output/miami/',output_pdf = op,mani_epic = mani_epic,bonf = bonf)
smoke_sig <- smokers %>% filter(pvalue<=bonf) %>%dplyr::select(CpG)
form_sig <- formers %>% filter(pvalue<=bonf) %>%dplyr::select(CpG)
model1_form_biomarker <- anti_join(form_sig,smoke_sig)

###########################################################################################################################
## Model 3
smokers=fread(paste0('output/svyglm_anno/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("CpG","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
CpGs <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

#smokers$chromosome <- factor(smokers$chromosome, chrs)

formers=fread(paste0('output/svyglm_anno/unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("CpG","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
CpGf <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

op<-'msa_no_overlap_london_model3_weighted_euro_miami.pdf'
plot_miami(smoke = smokers,form = formers,path = 'output/miami/',output_pdf = op,mani_epic = mani_epic,bonf = bonf)

smoke_sig <- smokers %>% filter(pvalue<=bonf) %>%dplyr::select(CpG)
form_sig <- formers %>% filter(pvalue<=bonf) %>%dplyr::select(CpG)
model3_form_biomarker <- anti_join(form_sig,smoke_sig)

##############################################################################################################
## Transancestry
#############################################################################################################
## Model 10
smokers=fread(paste0('output/svyglm_anno/unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("CpG","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
CpGs <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

#smokers$chromosome <- factor(smokers$chromosome, chrs)

formers=fread(paste0('output/svyglm_anno/unimputed_unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("CpG","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
CpGf <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

op<-'msa_no_overlap_london_model10_unweighted_transanc_miami.pdf'
plot_miami(smoke = smokers,form = formers,path = 'output/miami/',output_pdf = op,mani_epic = mani_epic,bonf = bonf)

smoke_sig <- smokers %>% filter(pvalue<=bonf) %>%dplyr::select(CpG)
form_sig <- formers %>% filter(pvalue<=bonf) %>%dplyr::select(CpG)
model10_form_biomarker <- anti_join(form_sig,smoke_sig)


###########################################################################################################################

smokers=fread(paste0('output/svyglm_anno/unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("CpG","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
CpGs <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

#smokers$chromosome <- factor(smokers$chromosome, chrs)

formers=fread(paste0('output/svyglm_anno/unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("CpG","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
CpGf <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

op<-'msa_no_overlap_london_model12_weighted_transanc_miami.pdf'
plot_miami(smoke = smokers,form = formers,path = 'output/miami/',output_pdf = op,mani_epic = mani_epic,bonf = bonf)
smoke_sig <- smokers %>% filter(pvalue<=bonf) %>%dplyr::select(CpG)
form_sig <- formers %>% filter(pvalue<=bonf) %>%dplyr::select(CpG)
model12_form_biomarker <- anti_join(form_sig,smoke_sig)


form_biomarker <- vctrs::vec_c(model1_form_biomarker,model3_form_biomarker, model10_form_biomarker,model12_form_biomarker) %>% left_join(.,mani,by=c("CpG"="probeID"))
write.table(form_biomarker,file = "output/miami/former_smoker_biomarker.tsv",sep='\t',quote = F,row.names=F)
#################################################################################################################
################# Only on MSA
##################################################################################################################
plot_miami_new_msa<- function(smoke, form, path,output_pdf,mani_epic,bonf) {
  previous_array_cpgs <- unique(c(mani_epic$CpG, mani_epic2$CpG, mani_450$CpG))
  plot_df <- bind_rows(smoke,form) %>% 
    group_by(chr) %>% 
    # Compute chromosome size
    summarise(chrlength = max(pos)) %>%  
    # Calculate cumulative position of each chromosome
    mutate(cumulativechrlength = cumsum(as.numeric(chrlength))-chrlength) %>% 
    dplyr::select(-chrlength) %>%
    # Temporarily add the cumulative length of each chromosome to the initial 
    # dataset 
    left_join(bind_rows(top = smokers,bottom = formers,.id="status"), ., by=c("chr"="chr")) %>%
    # Sort by chr then position 
    arrange(chr, pos) %>%
    # Add the position to the cumulative chromosome length to get the position of 
    # this probe relative to all other probes
    mutate(rel_pos = pos + cumulativechrlength) %>%
    # Calculate the logged p-value too
    mutate(logged_p = -log10(pvalue)) %>%
    dplyr::select(-cumulativechrlength)
  
  axis_df <- plot_df %>% 
    group_by(chr) %>% 
    summarize(chr_center=(max(rel_pos) + min(rel_pos)) / 2)
  maxp <- ceiling(max(plot_df$logged_p, na.rm = TRUE))
  
  plot_df <- plot_df %>% mutate(CpGe=gsub("\\_.*","",CpG))
  
  upper_labels <- plot_df %>%
    filter(status == "top") %>%
    filter(!CpGe %in% previous_array_cpgs) %>%
    # filter(!CpGe %in% lond_sig1$CpG) %>%
    # filter(!gene %in% lond_sig1$genesUniq) %>%
    filter(pvalue <= bonf) %>%
    group_by(chr) %>%
    arrange(desc(logged_p)) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    arrange(desc(logged_p)) %>%
    dplyr::slice(1:12) %>%
    dplyr::select(CpG)
  
  lower_labels <- plot_df %>%
    filter(status == "bottom") %>%
    filter(!CpGe %in% previous_array_cpgs) %>%
    filter(pvalue <= bonf) %>%
    group_by(chr) %>%
    arrange(desc(logged_p)) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    arrange(desc(logged_p)) %>%
    dplyr::slice(1:12) %>%
    dplyr::select(CpG)
  
  label_list <- rbind(upper_labels, lower_labels)
  
  plot_df <- plot_df %>%
    mutate(label = case_when(CpG %in% label_list$CpG ~ 
                               paste0(CpG, "\n", gene), TRUE ~ ""))
  upper_plot <- ggplot() + 
    geom_point(data = plot_df[which(plot_df$status == "top" ),], 
               aes(x = rel_pos, y = logged_p, color = as.factor(chr)), 
               size = 1) +
    scale_color_manual(values = rep(c("#3261cf", "#88a1db"), nrow(axis_df))) +
    scale_x_continuous(labels = axis_df$chr, 
                       breaks = axis_df$chr_center, 
                       expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(0, maxp+5), 
                       expand = expansion(mult = c(0.02, 0))) + 
    geom_hline(yintercept = -log10(1e-5), color = "grey", linetype = "dashed", 
               size = 0.5) +
    geom_hline(yintercept = -log10(2.170808e-07), color = "#150c6b", linetype = "solid", 
               size = 0.5) +
    labs(x = "", y = bquote(atop('Current smoker vs. Never smoker', '-log'[10]*'(p)'))) + 
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(),
          plot.margin = margin(b = 0, l = 10),
          axis.text=element_text(size=16),
          axis.title=element_text(size=24,face="bold")) + 
    geom_text_repel(data = plot_df%>% filter(plot_df$status == "top",CpG %in% upper_labels$CpG ),
                     aes(x = rel_pos, y = logged_p, label = label),
                     size = 6, segment.size = 0.2, point.padding = 0.3,
                     ylim = c(maxp / 3, NA), box.padding = 2,force=4,force_pull=0.1) 
  lower_plot <- ggplot() + 
    geom_point(data = plot_df[which(plot_df$status == "bottom" ),], 
               aes(x = rel_pos, y = logged_p, color = as.factor(chr)), 
               size = 1) +
    scale_color_manual(values = rep(c("#a33636", "#d99898"), nrow(axis_df))) +
    scale_x_continuous(breaks = axis_df$chr_center, 
                       position = "top",
                       expand = expansion(mult = 0.01)) +
    scale_y_reverse(limits = c(maxp+5, 0), 
                    expand = expansion(mult = c(0, 0.02))) + 
    geom_hline(yintercept = -log10(1e-5), color = "grey", linetype = "dashed", 
               size = 0.5) +
    geom_hline(yintercept = -log10(2.170808e-07), color = "#150c6b", linetype = "solid", 
               size = 0.5) +
    labs(x = "", y = bquote(atop('Former smoker vs. Never smoker', '-log'[10]*'(p)'))) + 
    theme_classic() +
    theme(legend.position = "none", 
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = margin(t = 0, l = 10),
          axis.text=element_text(size=16),
          axis.title=element_text(size=24,face="bold")) + 
    geom_text_repel(data = 
                       plot_df%>% filter(plot_df$status == "bottom",CpG %in% lower_labels$CpG ),
                     aes(x = rel_pos, y = logged_p, label = label),
                     size = 6, segment.size = 0.2, point.padding = 0.3,
                     ylim = c(NA, -(maxp / 3)), box.padding = 2,force=4,force_pull=0.1)     +
    annotate(
      "text",
      x = Inf,
      y = maxp + 4,
      label = sprintf("Bonferroni cutoff: p = %.3e", bonf),
      hjust = 1.1,
      vjust = 0,
      size = 4,
      color = "#150c6b"
    )
  
  #lower_plot
  pdf(file = paste0(path,output_pdf),height = 18,width = 20)
  p <- gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
  
  print(p)
  dev.off()
}



###############################################################################################################################
### Model 1

smokers=fread(file=paste0('output/svyglm_anno/unimputed_unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("CpG","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
CpGs <- smokers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

#smokers$chromosome <- factor(smokers$chromosome, chrs)

formers=fread(file=paste0('output/svyglm_anno/unimputed_unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')) %>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("CpG","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
CpGf <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

op <-'msa_new_cpg_model1_unweighted_euro_miami.pdf'
plot_miami_new_msa(smoke = smokers,form = formers,path = 'output/miami/',output_pdf = op,mani_epic = mani_epic,bonf = bonf)
###########################################################################################################################
## Model 3
smokers=fread(paste0('output/svyglm_anno/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("CpG","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
CpGs <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

#smokers$chromosome <- factor(smokers$chromosome, chrs)

formers=fread(paste0('output/svyglm_anno/unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("CpG","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
CpGf <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)
op <-'msa_new_cpg_model3_weighted_euro_miami.pdf'
plot_miami_new_msa(smoke = smokers,form = formers,path = 'output/miami/',output_pdf = op,mani_epic = mani_epic,bonf = bonf)

##############################################################################################################
## Transancestry
#############################################################################################################
## Model 10

smokers=fread(paste0('output/svyglm_anno/unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("CpG","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
CpGs <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

#smokers$chromosome <- factor(smokers$chromosome, chrs)

formers=fread(paste0('output/svyglm_anno/unimputed_unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("CpG","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
CpGf <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)
op <-'msa_new_cpg_model10_unweighted_transanc_miami.pdf'
plot_miami_new_msa(smoke = smokers,form = formers,path = 'output/miami/',output_pdf = op,mani_epic = mani_epic,bonf = bonf)



###########################################################################################################################
## Model 12
smokers=fread(paste0('output/svyglm_anno/unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("CpG","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
CpGs <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

#smokers$chromosome <- factor(smokers$chromosome, chrs)

formers=fread(paste0('output/svyglm_anno/unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("CpG","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
CpGf <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)
op <-'msa_new_cpg_model12_weighted_transanc_miami.pdf'
plot_miami_new_msa(smoke = smokers,form = formers,path = 'output/miami/',output_pdf = op,mani_epic = mani_epic,bonf = bonf)


###########################################################################################################################
###########################################################################################################################
## MSA only stacked

plot_miami_new_msa<- function(smoke, form, path,output_pdf,mani_epic,bonf) {
  previous_array_cpgs <- unique(c(mani_epic$CpG, mani_epic2$CpG, mani_450$CpG))
  
  plot_df <- bind_rows(smoke,form) %>% 
    group_by(chr) %>% 
    # Compute chromosome size
    summarise(chrlength = max(pos)) %>%  
    # Calculate cumulative position of each chromosome
    mutate(cumulativechrlength = cumsum(as.numeric(chrlength))-chrlength) %>% 
    dplyr::select(-chrlength) %>%
    # Temporarily add the cumulative length of each chromosome to the initial 
    # dataset 
    left_join(bind_rows(top = smokers,bottom = formers,.id="status"), ., by=c("chr"="chr")) %>%
    # Sort by chr then position 
    arrange(chr, pos) %>%
    # Add the position to the cumulative chromosome length to get the position of 
    # this probe relative to all other probes
    mutate(rel_pos = pos + cumulativechrlength) %>%
    # Calculate the logged p-value too
    mutate(logged_p = -log10(pvalue)) %>%
    dplyr::select(-cumulativechrlength)
  
  axis_df <- plot_df %>% 
    group_by(chr) %>% 
    summarize(chr_center=(max(rel_pos) + min(rel_pos)) / 2)
  
  maxp <- ceiling(max(plot_df$logged_p[plot_df$status == "top"], na.rm = TRUE))
  maxp_lower <- ceiling(max(plot_df$logged_p[plot_df$status == "bottom"], na.rm = TRUE))
  
  plot_df <- plot_df %>% mutate(CpGe=gsub("\\_.*","",CpG))
  
  upper_labels<- plot_df %>%
    filter(status == "top") %>%
    filter(!CpGe %in% previous_array_cpgs) %>%
    #filter(!CpGe %in% lond_sig1$CpG) %>%
    #filter(!gene %in% lond_sig1$genesUniq) %>%
    filter(pvalue<=bonf) %>%
    group_by(chr) %>%
    arrange(desc(logged_p)) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    arrange(desc(logged_p)) %>%
    dplyr::slice(1:12) %>%
    dplyr::select(CpG)
  
  lower_labels <- plot_df %>%
    filter(status == "bottom") %>%
    filter(!CpGe %in% previous_array_cpgs) %>%
    filter(pvalue<=bonf) %>%
    group_by(chr) %>%
    arrange(desc(logged_p)) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    arrange(desc(logged_p)) %>%
    dplyr::slice(1:12) %>%
    dplyr::select(CpG)
  
  label_list <- rbind(upper_labels, lower_labels)
  
  plot_df <- plot_df %>%
    mutate(label = case_when(CpG %in% label_list$CpG ~ 
                               paste0(CpG, "\n", gene), TRUE ~ ""))
  
  upper_plot <- ggplot() + 
    geom_point(data = plot_df[which(plot_df$status == "top" ),], 
               aes(x = rel_pos, y = logged_p, color = as.factor(chr)), 
               size = 1) +
    scale_color_manual(values = rep(c("#3261cf", "#88a1db"), nrow(axis_df))) +
    scale_x_continuous(labels = axis_df$chr, 
                       breaks = axis_df$chr_center, 
                       expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(0, maxp+5), 
                       expand = expansion(mult = c(0.02, 0))) + 
    geom_hline(yintercept = -log10(1e-5), color = "grey", linetype = "dashed", 
               size = 0.5) +
    geom_hline(yintercept = -log10(bonf), color = "#150c6b", linetype = "solid", 
               size = 0.5) +
    labs(x = "", y = bquote(atop('Current smoker vs. Never smoker', '-log'[10]*'(p)'))) + 
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(),
          plot.margin = margin(b = 0, l = 10),
          axis.text=element_text(size=16),
          axis.title=element_text(size=24,face="bold")) + 
    geom_text_repel(data = plot_df%>% filter(plot_df$status == "top",CpG %in% upper_labels$CpG ),
                    aes(x = rel_pos, y = logged_p, label = label),
                    size = 6, segment.size = 0.2, point.padding = 0.3,
                    ylim = c(maxp / 3, NA), box.padding = 2,force=4,force_pull=0.1) +
    annotate("text", x = Inf, y = maxp + 4,
             label = sprintf("Bonferroni cutoff: p = %.3e", bonf),
             hjust = 1.1, vjust = 1, size = 4, color = "#150c6b")
  
  lower_plot <- ggplot() + 
    geom_point(data = plot_df[which(plot_df$status == "bottom" ),], 
               aes(x = rel_pos, y = logged_p, color = as.factor(chr)), 
               size = 1) +
    scale_color_manual(values = rep(c("#a33636", "#d99898"), nrow(axis_df))) +
    scale_x_continuous(labels = axis_df$chr,
                       breaks = axis_df$chr_center, 
                       expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(0, maxp_lower+5), 
                       expand = expansion(mult = c(0.02, 0))) + 
    geom_hline(yintercept = -log10(1e-5), color = "grey", linetype = "dashed", 
               size = 0.5) +
    geom_hline(yintercept = -log10(bonf), color = "#150c6b", linetype = "solid", 
               size = 0.5) +
    labs(x = "", y = bquote(atop('Former smoker vs. Never smoker', '-log'[10]*'(p)'))) + 
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(),
          plot.margin = margin(t = 0, l = 10),
          axis.text=element_text(size=16),
          axis.title=element_text(size=24,face="bold")) + 
    geom_text_repel(data = 
                      plot_df%>% filter(plot_df$status == "bottom",CpG %in% lower_labels$CpG ),
                    aes(x = rel_pos, y = logged_p, label = label),
                    size = 6, segment.size = 0.2, point.padding = 0.3,
                    ylim = c(maxp_lower / 3, NA), box.padding = 2,force=4,force_pull=0.1)
  
  #lower_plot
  pdf(file = paste0(path,output_pdf),height = 18,width = 20)
  p <- gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
  
  print(p)
  dev.off()
}


###############################################################################################################################
### Model 1

smokers=fread(file=paste0('output/svyglm_anno/unimputed_unweighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("CpG","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
CpGs <- smokers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

#smokers$chromosome <- factor(smokers$chromosome, chrs)

formers=fread(file=paste0('output/svyglm_anno/unimputed_unweighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv')) %>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("CpG","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
CpGf <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

op <-'msa_new_cpg_model1_unweighted_euro_stacked.pdf'
plot_miami_new_msa(smoke = smokers,form = formers,path = 'output/miami/',output_pdf = op,mani_epic = mani_epic,bonf = bonf)

###########################################################################################################################
## Model 3

smokers=fread(paste0('output/svyglm_anno/unimputed_weighted_smokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("CpG","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
CpGs <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

#smokers$chromosome <- factor(smokers$chromosome, chrs)

formers=fread(paste0('output/svyglm_anno/unimputed_weighted_formersmokers_EWAS_Results_europc_only_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("CpG","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
CpGf <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

op <-'msa_new_cpg_model3_weighted_euro_stacked.pdf'
plot_miami_new_msa(smoke = smokers,form = formers,path = 'output/miami/',output_pdf = op,mani_epic = mani_epic,bonf = bonf)

##############################################################################################################
## Transancestry
#############################################################################################################
## Model 10

smokers=fread(paste0('output/svyglm_anno/unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("CpG","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
CpGs <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

#smokers$chromosome <- factor(smokers$chromosome, chrs)

formers=fread(paste0('output/svyglm_anno/unimputed_unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("CpG","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
CpGf <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

op <-'msa_new_cpg_model10_unweighted_transanc_stacked.pdf'
plot_miami_new_msa(smoke = smokers,form = formers,path = 'output/miami/',output_pdf = op,mani_epic = mani_epic,bonf = bonf)

###########################################################################################################################
## Model 12

smokers=fread(paste0('output/svyglm_anno/unimputed_weighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(smokers)<- c("CpG","chr","pos","pvalue","gene")
smokers$chr<-smokers$chr %>% str_replace("chr","") %>% as.numeric()
CpGs <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

#smokers$chromosome <- factor(smokers$chromosome, chrs)

formers=fread(paste0('output/svyglm_anno/unimputed_weighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)%>% filter(CpG_chrm %in% chrs)
names(formers)<-c("CpG","chr","pos","pvalue","gene")
formers$chr<-formers$chr %>% str_replace("chr","") %>% as.numeric()
CpGf <- formers %>% filter(pvalue<=0.0000002170808) %>% dplyr::select(CpG)

op <-'msa_new_cpg_model12_weighted_transanc_stacked.pdf'
plot_miami_new_msa(smoke = smokers,form = formers,path = 'output/miami/',output_pdf = op,mani_epic = mani_epic,bonf = bonf)

