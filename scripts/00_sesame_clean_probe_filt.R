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

dir.create("output", recursive=TRUE, showWarnings=FALSE)

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
bs = c("WG0342734_DNA_C08", "WG0342734_DNA_D10", "WG0342742_DNA_G11", 
       "WG0342736_DNA_D02", "WG0342736_DNA_G11", "WG0342737_DNA_C03", 
       "WG0342736_DNA_E07", "WG0342737_DNA_C01", "WG0342742_DNA_G02", 
       "WG0342737_DNA_B05", "WG0342737_DNA_B06", "WG0342734_DNA_B01", 
       "WG0342738_DNA_E03")
sampsheet <- sampsheet %>% filter(!Sample_ID %in% bs)
##idat directory
idats <- paste0("input_idats/",sampsheet$idat_prefix)


## read in idats and perform GQCDPBM
betasG = openSesame(idats,prep="GQCDPBM", manifest = addr, BPPARAM = BiocParallel::MulticoreParam(50))
## Read in p values
pvals = openSesame(idats, func = pOOBAH,manifest = addr, return.pval=TRUE , BPPARAM = BiocParallel::MulticoreParam(50))

## load individuals with sex differences infersex()
load('output/inferredsex.rdata')
sx<-rownames_to_column(data.frame(inrdsx),var="Sample_ID")
rm(inrdsx)
sampsheet <- left_join(sampsheet,sx,by="Sample_ID")
sampsheet <- sampsheet %>% mutate(sx = tolower(inrdsx))
sampsheet$Gender <- as.factor(sampsheet$Gender)
sampsheet$sx <-as.factor(sampsheet$sx)
#ev<-fread('genotype/for_Jason_11-2024/plink.eigenvec')

## link all IDs
ids = fread('samplesheets/linkage_IDs_msa.txt')

## attach PCs for all individuals except for 5 gwas failures
sampsheet<- sampsheet %>% left_join(ids,by='LHID')
#ev<-ev %>% mutate(BestArrayID = gsub('\\.','-',paste0(V1,'_',V2))) %>%  select(BestArrayID,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20,V21,V22)
#colnames(ev)<-c('BestArrayID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','PC16','PC17','PC18','PC19','PC20')
#sampsheet <- sampsheet %>% left_join(ev,by='BestArrayID')
ev <- fread("samplesheets/all_projected.sscore") %>%
  mutate(BestArrayID=gsub("\\.", "-", paste0(`#FID`, "_", IID))) %>%
  dplyr::select(BestArrayID, all_of(paste0("PC", 1:20, "_AVG")))

colnames(ev) <- c("BestArrayID", paste0("PC", 1:20))

sampsheet <- sampsheet %>%
  dplyr::select(-any_of(paste0("PC", 1:20))) %>%
  left_join(ev, by="BestArrayID")

colSums(is.na(sampsheet[paste0("PC", 1:20)]))
## Remove sex mismatch, gwas failures, and normalized density outliers
sampsheet = sampsheet %>% filter(Gender==sx , !is.na(PC1))
sampsheet = sampsheet %>% filter(idat_prefix != "207854510070_R04C01")
betasGF<- betasG[,colnames(betasG) %in% sampsheet$idat_prefix] #284317    829

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
betasGF <- betasGF[,miss_sampleGF <= 0.10] #249885    249826    751
dim(betasGF)
#249892/284317 = 88%


## Mask XY probes and SNP probes and Zhou recommended probes
masksT = msa_masks %>% filter(M_general == TRUE) %>% select(Probe_ID)
betasGF <- betasGF[!rownames(betasGF) %in% masksT$Probe_ID,] # 236758    765
dim(betasGF)
xy <- manifest %>% filter(CpG_chrm == 'chrX' | CpG_chrm == 'chrY') %>% select(Probe_ID)
betasGF <- betasGF[!rownames(betasGF) %in% xy$Probe_ID,]
dim(betasGF) # 230275    751

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
write.table(missingPerc_sampleGF,'output/missing_perc_751_pc_clean.tsv',quote = F,row.names = F,col.names = T,na='',sep = '\t')

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
##Plot densities
pdf("output/751_clean_pc_normalized_densities.pdf")
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
pdf("output/751_clean_pc_svd_imp.pdf")
champ.SVD(beta = as.data.frame(e_knn_impGF$data),pd = as.data.frame(svdin))
plot_MDS(mylabel = "LH_smoke_status", myLabel2="knn_imp", myeset = eset_betas_knnGF,1,2)
plot_MDS(mylabel = "LH_smoke_status", myLabel2="knn_imp", myeset = eset_betas_knnGF,2,3)
plot_MDS(mylabel = "LHSSCREEN_12142015", myLabel2="knn_imp", myeset = eset_betas_knnGF,1,2)
plot_MDS(mylabel = "LHSSCREEN_12142015", myLabel2="knn_imp", myeset = eset_betas_knnGF,2,3)
dev.off()
biological_vars <- c(
  "LH_smoke_status",
  "ANYETS",
  "LHSSCREEN_12142015",
  "Gender",
  "LH_BMI",
  "AgeAtHomeVisit",
  "STATE",
  paste0("PC", 1:20)
)
batch_vars <- c(
  "Sentrix_ID",
  "Sentrix_Position",
  "BOX"
)

setdiff(
  c(biological_vars, batch_vars),
  names(pdGF)
)
## Run combat to remove batch effect
cinGF<- sampsheetGF %>% mutate_at(c("Sentrix_ID", "Sentrix_Position", "Sample_Plate", "Sample_Well",
                                     "Sample_Name", "BOX", "LHID", "BSI_ID", "Species", "Tissue.Source", "STATE", "Gender",
                                     "European", "AgeAtHomeVisit", "LHSSCREEN_12142015", "LH_smoke_status", "LEL19","ANYETS"),as.factor) %>%
                    select("Sample_ID","Sentrix_ID", "Sentrix_Position", "Sample_Plate", "Sample_Well", "BOX",
                           "Sample_Name","STATE", "Gender",
                           "European", "AgeAtHomeVisit", "LHSSCREEN_12142015", "LH_smoke_status", "LEL19")
pdf("output/751_clean_pc_svd_imp.pdf")
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
# myCombatMF <- champ.runCombat(beta = as.data.frame(e_knn_impMF$data),pd = as.data.frame(cinGF),batchname=c("Sentrix_ID","Sentrix_Position","BOX"),variablename = "LH_smoke_status",logitTrans=F)
# eset_betas_knn_combatMF <- ExpressionSet(assayData = myCombatMF,phenoData = AnnotatedDataFrame(sampsheetGF),featureData = AnnotatedDataFrame(manifestGF))
# 
# 
# myCombatGF <- champ.runCombat(beta = as.data.frame(e_knn_impGF$data),pd = as.data.frame(cinGF),batchname=c("Sentrix_ID","Sentrix_Position","BOX"),variablename = "LH_smoke_status")
# eset_betas_knn_combatGF <- ExpressionSet(assayData = myCombatGF,phenoData = AnnotatedDataFrame(sampsheetGF),featureData = AnnotatedDataFrame(manifestGF))

library(Biobase)
library(sva)

## Fresh extraction: all 751 samples
betaGF <- Biobase::exprs(eset_betas_knnGF)
pdGF   <- Biobase::pData(eset_betas_knnGF)

stopifnot(
  ncol(betaGF) == 751,
  nrow(pdGF) == 751,
  all(colnames(betaGF) == rownames(pdGF))
)

## Preserve these variables during batch correction
## ANYETS is excluded because it has 36 missing values
biological_vars <- c(
  "LH_smoke_status",
  "LHSSCREEN_12142015",
  "Gender",
  "LH_BMI",
  "AgeAtHomeVisit",
  "STATE",
  paste0("PC", 1:20)
)

batch_vars <- c(
  "Sentrix_ID",
  "Sentrix_Position",
  "BOX"
)

required_vars <- c(
  biological_vars,
  batch_vars
)

stopifnot(
  all(required_vars %in% names(pdGF)),
  all(colSums(is.na(pdGF[, required_vars, drop = FALSE])) == 0)
)

## Convert beta values to M-values
eps <- 1e-6
betaGF <- pmin(pmax(betaGF, eps), 1 - eps)
mGF <- log2(betaGF / (1 - betaGF))

## Sequential batch correction
for (i in seq_along(batch_vars)) {

  later_batches <- if (i < length(batch_vars)) {
    batch_vars[(i + 1):length(batch_vars)]
  } else {
    character(0)
  }

  mod <- model.matrix(
    reformulate(
      c(biological_vars, later_batches)
    ),
    data = pdGF
  )

  stopifnot(
    nrow(mod) == 751,
    qr(mod)$rank == ncol(mod)
  )

  message("Correcting: ", batch_vars[i])

  mGF <- sva::ComBat(
    dat = mGF,
    batch = factor(pdGF[[batch_vars[i]]]),
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )
}

## Return to beta-value scale
myCombatGF <- plogis(mGF * log(2))

stopifnot(
  ncol(myCombatGF) == 751,
  all(colnames(myCombatGF) == rownames(pdGF))
)

dim(myCombatGF)
range(myCombatGF)

eset_betas_knn_combatGF <- ExpressionSet(assayData = myCombatGF,phenoData = AnnotatedDataFrame(sampsheetGF),featureData = AnnotatedDataFrame(manifestGF))
eset_mvals_knn_combatGF <- ExpressionSet(assayData = mGF,phenoData = AnnotatedDataFrame(sampsheetGF),featureData = AnnotatedDataFrame(manifestGF))




## Rerun plots

pdf("output/751_clean_pc_svd_knn_combat.pdf")
champ.SVD(beta = as.data.frame(myCombatGF),pd = as.data.frame(svdin))
dev.off()
pdf("output/751_clean_pc_knn_imp_combat.pdf")
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

betas = exprs(eset_betas_knn_combatGF)
beta2 = t(betas) # methylation data
mvals = exprs(eset_mvals_knn_combatGF)
mvals2 = t(mvals)
pheno = sampsheetGF
filt_pvals<- pvals[,colnames(pvals) %in% sampsheetGF$idat_prefix]
filt_pvals <- filt_pvals[rownames(filt_pvals) %in% rownames(betas),]
colnames(filt_pvals) <- sampsheetGF$Sample_Name[match(colnames(filt_pvals), sampsheetGF$idat_prefix)]

## P value mask matched pvalue and beta dataframes
betas_miss_masked <- betas
mvals_miss_masked <- mvals
identical(colnames(betas_miss_masked),colnames(eGF))
identical(rownames(betas_miss_masked),rownames(eGF))
betas_miss_masked[is.na(eGF)] <-NA
mvals_miss_masked[is.na(mGF)] <-NA


save(addr,betas,beta2,sampsheetGF,pheno,pvals,filt_pvals,betas_miss_masked,mvals_miss_masked,mvals,BloodFrac.m,manifest,file='751_clean_pc.rdata')
