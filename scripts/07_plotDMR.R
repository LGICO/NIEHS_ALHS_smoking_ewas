#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("karyoploteR")
library(karyoploteR)
library(AnnotationHub)
library(data.table)
#library(sesame)
library(dplyr)
library(gtools)
library(ggplot2)
library(ggrepel)
library(stringr)
library(gdata)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(ChIPseeker)
library(plyranges)
#ahub <- AnnotationHub()
#ahub["AH5086"]

tsv_path <- "references/MSA.hg38.manifest.tsv.gz"
path = ""
bonf<- 0.05/230325
bonf_epic = 0.05/817235
#sesameDataCache()
manifest = fread(tsv_path)
load(paste0(path,'output/cleaned_MSA_niehs_alhs.rdata'))
mani = fread(paste0(path,"references/MSA.hg38.manifest.gencode.v41.tsv.gz")) %>% dplyr::select(CpG_chrm, CpG_beg, CpG_end, probe_strand,probeID) %>%
            dplyr::filter(probeID %in% rownames(clean_betas))
smokers=fread(paste0(path,'output/svyglm/unimputed_unweighted_smokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)
names(smokers)<- c("probeID","chr","pos","smokers_pvalue","gene")
smokers <- smokers %>% dplyr::select(probeID,smokers_pvalue)

formers=fread(paste0(path,'output/svyglm/unimputed_unweighted_formersmokers_EWAS_Results_transanc_smoke_gender_BMI_age_state_cells_pcs.tsv'))%>%
  dplyr::select(V1,CpG_chrm,CpG_beg,P,genesUniq)
names(formers)<- c("probeID","chr","pos","formers_pvalue","gene")
formers <- formers %>% dplyr::select(probeID,formers_pvalue)
pvals <- dplyr::left_join(smokers,formers)
clean_betas <- data.frame(probeID = row.names(clean_betas), clean_betas)
cb <- dplyr::left_join(pvals,clean_betas,by = 'probeID')



mani_cpg <- mani %>% mutate(probeID=gsub("\\_.*","",probeID)) %>% distinct()
## We need to make a granges object
mani_gr <- makeGRangesFromDataFrame(mani,keep.extra.columns = T,ignore.strand = F,
                                    seqnames.field = "CpG_chrm",start.field = "CpG_beg",
                                    end.field = "CpG_end",strand.field = "probe_strand",na.rm=T)
names(mani_gr)<-mani_gr$probeID
rownames(cb) <- cb$probeID
mcols(mani_gr)<-cb[names(mani_gr), ,on=.(probeID)]
mani_cpg <- makeGRangesFromDataFrame(mani,keep.extra.columns = T,ignore.strand = F,
                                    seqnames.field = "CpG_chrm",start.field = "CpG_beg",
                                    end.field = "CpG_end",strand.field = "probe_strand",na.rm=T)
names(mani_cpg)<-mani_cpg$probeID

mani_ug = fread(paste0(path,"MSA.hg38.manifest.gencode.v41.tsv.gz")) %>% dplyr::select(CpG_chrm, CpG_beg, CpG_end, probe_strand,probeID,genesUniq)
#chr2:232418299-232418301 --- chr2:232420896-232420898  ECEL1P1
library(GenomicState)
txdb <- gencode_txdb("47", "hg38", chrs = "chr2")

ecel.region<-toGRanges("2:232416025-232419910",genome = "hg38")
ecel1p1 <- subsetByOverlaps(mani_gr, ecel.region)
ecel_probe <- ecel1p1$probeID
ecel_probe <- c(ecel_probe,'cg21566642_TC11','cg21566642_TC12','cg21566642_TC13')
ecel_probe<-c(
  "cg04284215_BC11", "cg21566642_TC13", "cg21566642_TC11",
  "cg04284185_TC21", "cg21566642_TC12", "cg04284186_TC11",
  "cg04284209_BC11", "cg04284216_BC11", "cg04284207_BC11",
  "cg04284205_TC11", "cg06644428_BC13", "cg06644428_BC11",
  "cg06644428_BC12", "cg06644428_BC14", "cg04284171_TC11",
  "cg03329539_BC21", "cg04284170_TC11", "cg17087741_TC21",
  "cg21566642_TC11", "cg21566642_TC12", "cg21566642_TC13"
) 
## We are going to collapse everything into CpGs from probeIDs
cb_ecel <- cb %>% dplyr::filter(probeID %in% ecel_probe) %>% dplyr::mutate(probeID=gsub("\\_.*","",probeID)) %>% group_by(probeID) %>% summarise(across(everything(), list(mean)))
rownames(cb_ecel) <- cb_ecel$probeID
mani_ecel <-mani_cpg
mcols(mani_ecel)<-cb[names(mani_ecel), ,on=.(probeID)]
mani_ecel<- subsetByOverlaps(mani_ecel, ecel.region)
mani_cpg<- subsetByOverlaps(mani_ecel, ecel.region)
## Annotate
mani_cpg<- annotatePeak(mani_cpg, tssRegion=c(-3000, 3000),
                        TxDb=txdb, annoDb="org.Hs.eg.db",level = 'gene') %>% as_tibble()


######
ss<- clean_sampsheet %>% dplyr::select(Sample_ID,LH_smoke_status)
mcpg <- manifest %>% dplyr::mutate(cpg=gsub("\\_.*","",Probe_ID)) %>% dplyr::select(-Probe_ID)%>%
         dplyr::select(cpg,CpG_chrm,CpG_beg) %>% dplyr::distinct()
cb_ecel_pval <- cb_ecel %>% select(probeID,smokers_pvalue_1,formers_pvalue_1)
cb_ecel_betas <- cb_ecel %>% select(-smokers_pvalue_1,-formers_pvalue_1)
cb_ecel_betas_long <- cb_ecel_betas %>% tidyr::pivot_longer(cols=-probeID) %>%
        dplyr::mutate(name = substr(name,1,nchar(name)-2)) %>%
        dplyr::left_join(ss, by = c('name'='Sample_ID'))  %>% dplyr::filter(LH_smoke_status!="M") %>% 
        dplyr::left_join(mcpg,by=c('probeID'='cpg')) %>% dplyr::group_by(probeID) %>% arrange(CpG_beg)

names(cb_ecel_betas_long) <- c('cpg','sample','beta','smoke_status','chr','start')
cb_ecel_betas_long <- cb_ecel_betas_long %>% arrange(start)
library(forcats)

# order CpGs by genomic position
cpg_order <- cb_ecel_betas_long %>%
  distinct(cpg, start) %>%
  arrange(start) %>%
  pull(cpg)

cb_ecel_betas_long <- cb_ecel_betas_long %>%
  mutate(
    cpg = factor(cpg, levels = cpg_order)
  )
status_order <- cb_ecel_betas_long %>%
  count(smoke_status) %>%
  arrange(desc(n)) %>%
  pull(smoke_status)

cb_ecel_betas_long <- cb_ecel_betas_long %>%
  mutate(smoke_status = factor(smoke_status, levels = status_order))
cb_ecel_betas_long <- cb_ecel_betas_long %>%
  arrange(smoke_status, cpg)

custom.col <- c("","#7fcdbb","")
#custom.col <- c("#290101","#a33636","#db9a9a")
custom.col <- c(
  "Never" = "#db9a9a",
  "Former" = "#a33636",
  "Current" = "#290101"
)
#cb_ecel_betas_long$smoke_status <- factor(cb_ecel_betas_long$smoke_status,levels=c('Never','Former','Current'))
#group.colors = c(Current = "#3261cf", Former = "#a33636",Never)


p<- ggplot() +
  # largest group first (back)
  geom_dotplot(
    data = filter(cb_ecel_betas_long, smoke_status == "Never"),
    aes(x = cpg, y = beta, fill = smoke_status),
    binwidth = 0.006, stackratio = 0.5, stackdir = "center",
    binaxis = "y", stroke = NA, alpha = 0.9
  ) +
  geom_dotplot(
    data = filter(cb_ecel_betas_long, smoke_status == "Former"),
    aes(x = cpg, y = beta, fill = smoke_status),
    binwidth = 0.006, stackratio = 0.5, stackdir = "center",
    binaxis = "y", stroke = NA, alpha = 0.9
  ) +
  geom_dotplot(
    data = filter(cb_ecel_betas_long, smoke_status == "Current"),
    aes(x = cpg, y = beta, fill = smoke_status),
    binwidth = 0.006, stackratio = 0.5, stackdir = "center",
    binaxis = "y", stroke = NA, alpha = 0.9
  ) +
  geom_smooth(
    data = cb_ecel_betas_long,
    aes(x = cpg, y = beta, color = smoke_status, group = smoke_status,fill=smoke_status),
    method = "loess", size = 1.5, alpha=0.6
  ) +
  scale_fill_manual(values = custom.col) +
  scale_color_manual(values = custom.col) + theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) + scale_color_manual(
    values = custom.col,
    breaks = c("Never", "Former", "Current")
  ) +
  scale_fill_manual(
    values = custom.col,
    breaks = c("Never", "Former", "Current")
  )

p

pdf(paste0(path,"updated_dmr.pdf"), width = 10, height = 6)

p

dev.off()




library(ggpubr)
ggplot(cb_ecel_betas_long,aes(x=cpg,y=beta,dot.size=0.03,color=smoke_status,fill=smoke_status)) +
        geom_dotplot(binwidth=0.006,stackratio =0.5, stackdir='center',position='identity',binaxis = "y",stroke=NA, alpha=0.9) +
        geom_smooth(aes(group=smoke_status),method='loess') +
        labs(x = "CpG", y = bquote(atop('Beta Value'))) + 
        scale_fill_manual(values=custom.col,guide = "legend",name = "Smoking Status",aesthetics = c("colour", "fill"))+
        theme_classic() + ggtitle("Methylation beta values across ECEL1P1 on the MSA",)



ggplot(cb_ecel_betas_long,aes(x=cpg,y=beta,dot.size=0.03,color=smoke_status,fill=smoke_status)) +
  ggviolin(binwidth=0.006,stackratio =0.1, stackdir='centerwhole',position='dodge',binaxis = "y",stroke=NA, alpha=0.5) +
  geom_smooth(aes(group=smoke_status),method='loess') +
  labs(x = "CpG", y = bquote(atop('Beta Value'))) + 
  scale_color_discrete(guide = "legend",name = "Smoking Status",aesthetics = c("colour", "fill"))+
  theme_classic() + ggtitle("Methylation beta values across ECEL1P1 on the MSA")



save(cb_ecel_betas_long,file=paste0(path,'output/ecel1p1_beta_long.rdata'))













