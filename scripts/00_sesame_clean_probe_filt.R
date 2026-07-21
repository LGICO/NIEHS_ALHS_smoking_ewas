library(data.table)
library(tidyverse)
library(readxl)
library(knitr)
library(sesame)
library(wheatmap)
library(multtest)
library(limma)
library(RColorBrewer)
library(ggplot2)
library(minfi)
library(ChAMP)
library(impute)
#library(parallel)

setwd('')
## Run to initiate sesame
sesameDataCache()

## Read in manifest and snp files
tsv_path <- "references/MSA.hg38.manifest.tsv.gz"
msa_snps <- fread('references/MSA.hg38.snp.tsv.gz')
msa_masks <- fread('references/MSA.hg38.mask.tsv.gz')
manifest = sesameAnno_readManifestTSV(tsv_path)
manifest_450 = sesameAnno_readManifestTSV('references/HM450.hg38.manifest.tsv.gz')
manifest_epic = sesameAnno_readManifestTSV('references/EPIC.hg38.manifest.tsv.gz')
gencode = fread('references/MSA.hg38.manifest.gencode.v41.tsv.gz')

## read manifest to address file (for openSesame or readIDATpair)
addr = sesameAnno_buildAddressFile(tsv_path)

## Input directory and read in and modify sample sheet
IDAT_dir = "input_idats"
sampsheet <- read_csv("samplesheets/LHS_Sample_Sheet.csv") %>% mutate(Sample_ID = str_replace(Sample_ID,"-","_"))
sampsheet <- sampsheet %>% mutate(idat_prefix = paste0(Sentrix_ID,'_',Sentrix_Position))

##idat directory
idats <- paste0("input_idats/",sampsheet$idat_prefix)


## read in idats and perform GQCDPBM
betasG = openSesame(idats,prep="GQCDPBM", manifest = addr, BPPARAM = BiocParallel::MulticoreParam(50))
## Read in p values
pvals = openSesame(idats, func = pOOBAH,manifest = addr, return.pval=TRUE , BPPARAM = BiocParallel::MulticoreParam(50))

## load individuals with sex differences infersex()
sdfs <- openSesame(idats, prep = "GQCDPBM", func = NULL, manifest = addr,platform = "MSA")
inrdsx <- vapply(sdfs, function(sdf) inferSex(getBetas(sdf), platform = "MSA"), character(1))

sx<-rownames_to_column(data.frame(inrdsx),var="Sample_ID")
rm(inrdsx)
sampsheet <- left_join(sampsheet,sx,by="Sample_ID")
sampsheet <- sampsheet %>% mutate(sx = tolower(inrdsx))
sampsheet$Gender <- as.factor(sampsheet$Gender)
sampsheet$sx <-as.factor(sampsheet$sx)
ev<-fread('genotype/plink.eigenvec')

## link all IDs
ids = fread('linkage_IDs_msa.txt')

## attach PCs for all individuals except for 5 gwas failures
sampsheet<- sampsheet %>% left_join(ids,by='LHID')
ev<-ev %>% mutate(BestArrayID = gsub('\\.','-',paste0(V1,'_',V2))) %>%  select(BestArrayID,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20,V21,V22)
colnames(ev)<-c('BestArrayID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','PC16','PC17','PC18','PC19','PC20')
sampsheet <- sampsheet %>% left_join(ev,by='BestArrayID')

## Remove sex mismatch, gwas failures, and normalized density outliers
sampsheet = sampsheet %>% filter(Gender==sx , !is.na(PC1))
sampsheet = sampsheet %>% filter(idat_prefix != "207854510070_R04C01")
betasGF<- betasG[,colnames(betasG) %in% sampsheet$idat_prefix] #284317    843/2

## Determine missings
miss_sampleG <- apply(betasGF, 2, function(x) sum(is.na(x))/length(x) )
miss_probe_blankG <- apply(betasGF, 1, function(x) sum(is.na(x))/length(x) )

sum(miss_probe_blankG >=0.95) #16,661 | 16,845
masked_missing_probeG = miss_probe_blankG[miss_probe_blankG >=0.95]

## Filter missing 50%
betasGF <- betasGF[miss_probe_blankG < 0.50,] # 249885    842
##double check only ch and cg probes
as.data.frame(betasGF) %>% filter(grepl("cg", rownames(betasGF))) %>% dim(.)

## Now filter by 10% per sample
miss_sampleGF <- apply(betasGF, 2, function(x) sum(is.na(x))/length(x) )
sum(miss_sampleGF >= 0.10) ##843-78=765
betasGF <- betasGF[,miss_sampleGF <= 0.10] #249885    764
dim(betasGF)
#249892/284317 = 88%


## Mask XY probes and SNP probes and Zhou recommended probes
masksT = msa_masks %>% filter(M_general == TRUE) %>% select(Probe_ID)
betasGF <- betasGF[!rownames(betasGF) %in% masksT$Probe_ID,] # 236758    764
dim(betasGF)
xy <- manifest %>% filter(CpG_chrm == 'chrX' | CpG_chrm == 'chrY') %>% select(Probe_ID)
betasGF <- betasGF[!rownames(betasGF) %in% xy$Probe_ID,]
dim(betasGF) # 230325    764

##### ID density outlier
#	Coerce object to matrix
E <- as.matrix(betasGF)
narray <- ncol(E)
npoint <- 512
X <- Y <- matrix(0,npoint,narray)
for (a in 1:ncol(E)) {
	d <- density(E[,a],n=npoint,na.rm=TRUE)
	X[,a] <- d$x
	Y[,a] <- d$y
	}
for( i in 1:dim(X)[2]){
  for(j in 1:dim(X)[1]){
    a = X[j,i]
    b = Y[j,i]
    if((a<0.2 &  b > 2)){print(i)}
  }
}
colnames(E)[99] #"207854510070_R04C01"

#24,104,110,115,129,138,165,221,306,308,318,339,345,389,486,563,698,731,737,769,795,800


##cg Match rows and columns & generate mvals
############################
ssG <- data.frame(cbind(colnames(betasGF),sampsheet$idat_prefix[match(colnames(betasGF),sampsheet$idat_prefix)])) %>% filter(!is.na(X2))
betasGF <- betasGF %>% as.data.frame(.) %>% select(ssG$X2)
sampsheetGF <- sampsheet[match(colnames(betasGF),sampsheet$idat_prefix),]
colnames(betasGF) <- sampsheetGF$Sample_ID
MvalsGF <- BetaValueToMValue(betasGF)
colnames(MvalsGF) <- sampsheetGF$Sample_ID
row.names(sampsheetGF) <- sampsheetGF$Sample_ID

manifestGF <- manifest[match(rownames(betasGF), manifest$Probe_ID),]
sum(rownames(betasGF) != manifestGF$Probe_ID)
row.names(manifestGF) <- manifestGF$Probe_ID

## Create temporary ExpressionSet
eset_betasGF <- ExpressionSet(assayData = as.matrix(betasGF),phenoData = AnnotatedDataFrame(sampsheetGF),featureData = AnnotatedDataFrame(manifestGF))
eset_MvalsGF <- ExpressionSet(assayData = as.matrix(MvalsGF),phenoData = AnnotatedDataFrame(sampsheetGF),featureData = AnnotatedDataFrame(manifestGF))

## Print missing per sample
pdatGF <- pData(eset_betasGF)
eGF <- exprs(eset_betasGF)
pdatMF <- pData(eset_MvalsGF)
eMF <- exprs(eset_MvalsGF)
##Missings
miss_sampleGF <- apply(eGF, 2, function(x) sum(is.na(x))/length(x) )
miss_probeGF <- apply(eGF, 1, function(x) sum(is.na(x))/length(x))
miss_probe_fGF <- cut(miss_probeGF, breaks = seq(0, 1, 0.1))
missingPerc_sampleGF <- data.frame(Sample_ID = names(miss_sampleGF), smoke_status = pdatGF$LH_smoke_status, miss_sample = miss_sampleGF*100, stringsAsFactors = F)
missingPerc_sampleGF %>% kable()
write.table(missingPerc_sampleGF,'output/missing_perc_764.tsv',quote = F,row.names = F,col.names = T,na='',sep = '\t')

missing_sample <- miss_sampleGF
miss_probe <- miss_probeGF
miss_probe_f <- miss_probe_fGF
e <- eGF

## MDS Plot
plot_MDS <- function(mylabel, myLabel2,myeset,n,n2){
  group <- pData(myeset)[[mylabel]]
  col.group <- as.factor(group)
  levels(col.group) <- brewer.pal(nlevels(col.group), "Dark2")
  col.group.ch <- as.character(col.group)
  plotMDS(myeset, top = 2000,col = col.group.ch, gene.selection = "common", main = paste(mylabel,myLabel2 ), pch = 16,dim.plot = c(n,n2))
  legend("topleft", fill = levels(col.group), legend = levels(as.factor(group)))
}

## For SVD plot
svdin<- sampsheetGF %>% mutate_at(c("Sentrix_ID", "Sentrix_Position", "Sample_Plate", "Sample_Well", "Sample_Group",
                          "Sample_Name", "BOX", "LHID", "BSI_ID", "Species", "Tissue.Source", "STATE", "Gender",
                          "European", "AgeAtHomeVisit", "LHSSCREEN_12142015", "LH_smoke_status", "LEL19","ANYETS"),as.factor) %>%
              select("Sample_ID","Sentrix_ID", "Sentrix_Position", "Sample_Plate", "Sample_Well", "Sample_Group",
                     "Sample_Name", "BOX", "LHID", "BSI_ID", "Species", "Tissue.Source", "STATE", "Gender",
                     "European", "AgeAtHomeVisit", "LHSSCREEN_12142015", "LH_smoke_status", "LEL19","ANYETS")

system("mkdir -p output/")
##Plot densities
pdf("output/764_normalized_densities.pdf")
data.frame(missing_interval = miss_probe_f) %>%
  filter(!is.na(missing_interval)) %>%
  ggplot(aes(missing_interval)) +
  geom_bar()

#boxplot(e, las = 2, ylab = "Betas")
#abline(h = median(e), col = "blue")
plotDensities(eset_betasGF, main = "Betas n = 764 | 230,325", legend = FALSE)
plotDensities(eset_MvalsGF, main = "M values n = 764 | 230,325", legend = FALSE)
dev.off()

## KNN impute 10 nearest neighbors
e_knn_impGF <- impute.knn(eGF,k=10)
eset_betas_knnGF <- ExpressionSet(assayData = as.matrix(e_knn_impGF$data),phenoData = AnnotatedDataFrame(sampsheetGF),featureData = AnnotatedDataFrame(manifestGF))
e_knn_impMF <- impute.knn(eMF,k=10)
eset_mvals_knnGF <- ExpressionSet(assayData = as.matrix(e_knn_impMF$data),phenoData = AnnotatedDataFrame(sampsheetGF),featureData = AnnotatedDataFrame(manifestGF))


## Run SVD
pdf("output/764_svd_imp.pdf")
champ.SVD(beta = as.data.frame(e_knn_impGF$data),pd = as.data.frame(svdin))
plot_MDS(mylabel = "LH_smoke_status", myLabel2="knn_imp", myeset = eset_betas_knnGF,1,2)
plot_MDS(mylabel = "LH_smoke_status", myLabel2="knn_imp", myeset = eset_betas_knnGF,2,3)
plot_MDS(mylabel = "LHSSCREEN_12142015", myLabel2="knn_imp", myeset = eset_betas_knnGF,1,2)
plot_MDS(mylabel = "LHSSCREEN_12142015", myLabel2="knn_imp", myeset = eset_betas_knnGF,2,3)
dev.off()

## Run combat to remove batch effect
cinGF<- sampsheetGF %>% mutate_at(c("Sentrix_ID", "Sentrix_Position", "Sample_Plate", "Sample_Well",
                                     "Sample_Name", "BOX", "LHID", "BSI_ID", "Species", "Tissue.Source", "STATE", "Gender",
                                     "European", "AgeAtHomeVisit", "LHSSCREEN_12142015", "LH_smoke_status", "LEL19","ANYETS"),as.factor) %>%
                    select("Sample_ID","Sentrix_ID", "Sentrix_Position", "Sample_Plate", "Sample_Well", "BOX",
                           "Sample_Name","STATE", "Gender",
                           "European", "AgeAtHomeVisit", "LHSSCREEN_12142015", "LH_smoke_status", "LEL19")
myCombatMF <- champ.runCombat(beta = as.data.frame(e_knn_impMF$data),pd = as.data.frame(cinGF),batchname=c("Sentrix_ID","Sentrix_Position","BOX"),variablename = "LH_smoke_status",logitTrans=F)
eset_betas_knn_combatMF <- ExpressionSet(assayData = myCombatMF,phenoData = AnnotatedDataFrame(sampsheetGF),featureData = AnnotatedDataFrame(manifestGF))


myCombatGF <- champ.runCombat(beta = as.data.frame(e_knn_impGF$data),pd = as.data.frame(cinGF),batchname=c("Sentrix_ID","Sentrix_Position","BOX"),variablename = "LH_smoke_status")
eset_betas_knn_combatGF <- ExpressionSet(assayData = myCombatGF,phenoData = AnnotatedDataFrame(sampsheetGF),featureData = AnnotatedDataFrame(manifestGF))




## Rerun plots

pdf("output/764_svd_knn_combat.pdf")
champ.SVD(beta = as.data.frame(myCombatGF),pd = as.data.frame(svdin))
dev.off()
pdf("output/764_knn_imp_combat.pdf")
plot_MDS(mylabel = "LH_smoke_status", myLabel2="knn_imp_combat", myeset = eset_betas_knn_combatGF,1,2)
plot_MDS(mylabel = "LH_smoke_status", myLabel2="knn_imp_combat", myeset = eset_betas_knn_combatGF,2,3)
plot_MDS(mylabel = "LH_smoke_status", myLabel2="knn_imp_combat", myeset = eset_betas_knn_combatGF,3,4)

plot_MDS(mylabel = "Gender", myLabel2="knn_imp_combat", myeset = eset_betas_knn_combatGF,1,2)
plot_MDS(mylabel = "European", myLabel2="knn_imp_combat", myeset = eset_betas_knn_combatGF,1,2)
plot_MDS(mylabel = "STATE", myLabel2="knn_imp_combat", myeset = eset_betas_knn_combatGF,1,2)
plot_MDS(mylabel = "LHSSCREEN_12142015", myLabel2="knn_imp_combat", myeset = eset_betas_knn_combatGF,1,2)
plot_MDS(mylabel = "LEL19", myLabel2="knn_imp_combat", myeset = eset_betas_knn_combatGF,1,2)
plot_MDS(mylabel = "ANYETS", myLabel2="knn_imp_combat", myeset = eset_betas_knn_combatGF,1,2)
dev.off()

## Run epidish to get cell proportions
library(EpiDISH)
betas_knn_combat_epic = mLiftOver(exprs(eset_betas_knn_combatGF),target_platform='EPIC', source_platform='MSA',celltype='blood',BPPARAM=BiocParallel::MulticoreParam(2))
miss_probe_blank_epic <- apply(betas_knn_combat_epic, 1, function(x) sum(is.na(x))/length(x) )
betas_knn_combat_epic = betas_knn_combat_epic[miss_probe_blank_epic < 0.95,]
leuks <- estimateLeukocyte(betas_knn_combat_epic,platform='EPIC')
leuks = as.data.frame(leuks,row.names=names(leuks))
leuks <- rownames_to_column(leuks,"Sample_ID")
sampsheetGF = left_join(sampsheetGF, leuks,by = 'Sample_ID')
BloodFrac.m <- epidish(beta.m = betas_knn_combat_epic, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
sampsheetGF = left_join(sampsheetGF,rownames_to_column(as.data.frame(BloodFrac.m)),by=c("Sample_ID"="rowname"))
row.names(sampsheetGF)<-sampsheetGF$Sample_ID
eset_betas_knn_combatGF <- ExpressionSet(assayData = myCombatGF,phenoData = AnnotatedDataFrame(sampsheetGF),featureData = AnnotatedDataFrame(manifestGF))



## Final variables of interest

clean_betas = exprs(eset_betas_knn_combatGF)
mvals = exprs(eset_betas_knn_combatMF)
mvals2 = t(mvals)
#pheno = sampsheetGF
filt_pvals<- pvals[,colnames(pvals) %in% sampsheetGF$idat_prefix]
filt_pvals <- filt_pvals[rownames(filt_pvals) %in% rownames(betas),]
colnames(filt_pvals) <- sampsheetGF$Sample_Name[match(colnames(filt_pvals), sampsheetGF$idat_prefix)]
clean_sampsheet <- sampsheetGF
## P value mask matched pvalue and beta dataframes
betas_miss_masked <- clean_betas
mvals_miss_masked <- mvals
identical(colnames(betas_miss_masked),colnames(eGF))
identical(rownames(betas_miss_masked),rownames(eGF))
betas_miss_masked[is.na(eGF)] <-NA
mvals_miss_masked[is.na(eMF)] <-NA
clean_betas_pval_masked<- betas_miss_masked
system("mkdir -p output/rdata")
save(addr,clean_betas,clean_sampsheet,pvals,filt_pvals,clean_betas_pval_masked,BloodFrac.m,manifest,file='output/rdata/cleaned_MSA_niehs_alhs.rdata')
