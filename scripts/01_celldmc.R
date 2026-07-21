library(EpiDISH)
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
load('output/rdata/cleaned_MSA_niehs_alhs.rdata')
system("mkdir -p output/dmc")
pheno<- clean_sampsheet
#pheno=pheno%>%mutate(LEL19=ifelse(LEL19=="Yes",1,ifelse(LEL19=="No",0,99)))
#pheno$LEL19=as.factor(pheno$LEL19)


pheno$LH_smoke_status=factor(pheno$LH_smoke_status, levels = c("Never","Former","Current","M"))
pheno$European=factor(pheno$European,levels = c("European","Non_European"))
pheno <- pheno %>% mutate( ANYETS = case_when(ANYETS =="No"~"No",ANYETS=="Yes"~"Yes",.default = "M"))
pheno$ANYETS <- factor(pheno$ANYETS,levels=c("No","Yes","M"))
pheno$LHSSCREEN_12142015 = factor(pheno$LHSSCREEN_12142015,levels=c("Control","Case"))
pheno$Gender <- factor(pheno$Gender,levels = c("male","female"))
pheno$STATE <- factor(pheno$STATE, levels=c('IA','NC'))

## Split the data into never former and never current
nf = pheno %>% filter(LH_smoke_status == "Never" | LH_smoke_status == "Former")
nc = pheno %>% filter(LH_smoke_status == "Never" | LH_smoke_status == "Current")
fc = pheno %>% filter(LH_smoke_status == "Former" | LH_smoke_status == "Current")

betas_nf = clean_betas[,as.vector(nf$Sample_ID)]
betas_nc = clean_betas[,as.vector(nc$Sample_ID)]
betas_fc = clean_betas[,as.vector(fc$Sample_ID)]

ed_nf = nf %>% select(B, NK, CD4T, CD8T, Mono, Neutro)
ed_nc = nc %>% select(B, NK, CD4T, CD8T, Mono, Neutro)
ed_fc = fc %>% select(B, NK, CD4T, CD8T, Mono, Neutro)

cdmc_nf = CellDMC(beta.m=betas_nf, pheno.v=nf$LH_smoke_status,frac.m=as.matrix(ed_nf),
               cov.mod = model.matrix(~ Gender + AgeAtHomeVisit + STATE + ANYETS +
                                        LH_BMI + European + LHSSCREEN_12142015 +
                                        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
                                        PC8 + PC9 + PC10,data = nf),mc.cores=16)

cdmc_nc = CellDMC(beta.m=betas_nc, pheno.v=nc$LH_smoke_status,frac.m=as.matrix(ed_nc),
               cov.mod = model.matrix(~ Gender + AgeAtHomeVisit + STATE + ANYETS + 
                                        LH_BMI + European + LHSSCREEN_12142015 +
                                        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
                                        PC8 + PC9 + PC10,data = nc),mc.cores=16)

cdmc_fc = CellDMC(beta.m=betas_fc, pheno.v=fc$LH_smoke_status,frac.m=as.matrix(ed_fc),
                  cov.mod = model.matrix(~ Gender + AgeAtHomeVisit + STATE + ANYETS + 
                                           LH_BMI + European + LHSSCREEN_12142015 +
                                           PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
                                           PC8 + PC9 + PC10,data = fc),mc.cores=16)

cdmc_nf$dmct %>% as_tibble(.,rownames="CpG") %>% filter(DMC == 1) 
nccpg <- cdmc_nc$dmct %>% as_tibble(.,rownames="CpG") %>% filter(DMC == 1) %>% select(CpG)
fccpg <- cdmc_fc$dmct %>% as_tibble(.,rownames="CpG") %>% filter(DMC == 1) %>% select(CpG)
library(janitor)
ed_nc = nc %>% select(Sample_ID,B, NK, CD4T, CD8T, Mono, Neutro,LH_smoke_status)
ed_fc = fc %>% select(Sample_ID,B, NK, CD4T, CD8T, Mono, Neutro,LH_smoke_status)
ncbetas <- as_tibble(subset(betas_nc,rownames(betas) %in% nccpg$CpG),rownames="CpG") %>% t %>%
  as.data.frame %>% row_to_names(row_number = 1) 
ncbetas <- left_join(rownames_to_column(ncbetas),ed_nc,by =c("rowname"= 'Sample_ID'))
fcbetas <- as_tibble(subset(betas_fc,rownames(betas) %in% fccpg$CpG),rownames="CpG")%>% t %>%
  as.data.frame %>% row_to_names(row_number = 1) 
fcbetas <- left_join(rownames_to_column(fcbetas),ed_fc,by =c("rowname"= 'Sample_ID'))

library(cowplot)
library(egg)
pdf(file = "output/dmc/never_current_beta_dmc.pdf", width=12, height=4)
for(i in 2:(length(colnames(ncbetas))-7)){
    B <- ggplot(ncbetas, aes(x=as.numeric(B), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
      geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "B")
    NK <- ggplot(ncbetas, aes(x=as.numeric(NK), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
      geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "NK")
    CD4T <- ggplot(ncbetas, aes(x=as.numeric(CD4T), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
      geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "CD4T")
    CD8T<- ggplot(ncbetas, aes(x=as.numeric(CD8T), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
      geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "CD8T")
    Mono <- ggplot(ncbetas, aes(x=as.numeric(Mono), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
      geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "Mono")
    Neutro<- ggplot(ncbetas, aes(x=as.numeric(Neutro), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
      geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "Neutro")
    p <- plot_grid(B,NK,CD4T,CD8T,Mono,Neutro, labels = "AUTO",scale = 1)
    print(p)
}
dev.off()

pdf(file = "output/dmc/former_current_beta_dmc.pdf", width=12, height=4)
for(i in 2:(length(colnames(fcbetas))-7)){
  B <- ggplot(fcbetas, aes(x=as.numeric(B), y=as.numeric(!!sym(colnames(fcbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(fcbetas)[i],y = "B")
  NK <- ggplot(fcbetas, aes(x=as.numeric(NK), y=as.numeric(!!sym(colnames(fcbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(fcbetas)[i],y = "NK")
  CD4T <- ggplot(fcbetas, aes(x=as.numeric(CD4T), y=as.numeric(!!sym(colnames(fcbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(fcbetas)[i],y = "CD4T")
  CD8T<- ggplot(fcbetas, aes(x=as.numeric(CD8T), y=as.numeric(!!sym(colnames(fcbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(fcbetas)[i],y = "CD8T")
  Mono <- ggplot(fcbetas, aes(x=as.numeric(Mono), y=as.numeric(!!sym(colnames(fcbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(fcbetas)[i],y = "Mono")
  Neutro<- ggplot(fcbetas, aes(x=as.numeric(Neutro), y=as.numeric(!!sym(colnames(fcbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(fcbetas)[i],y = "Neutro")
  p <- plot_grid(B,NK,CD4T,CD8T,Mono,Neutro, labels = "AUTO",scale=1)
  print(p)
}
dev.off()


#####################
# Adjusting beta values for each CpG by cell proportions


betas_nc <- as_tibble(betas_nc,rownames="CpG")
# Load necessary libraries
library(dplyr)

# Create the beta values data frame
beta_values <- data.frame(
  CpG = c("cg00000109_BC21", "cg00000658_TC21", "cg00000714_TC21", 
          "cg00000721_TC21", "cg00001099_BC21"),
  A1 = c(0.8435673, 0.8267204, 0.2360347, 0.9349584, 0.6986932),
  A2 = c(0.8497762, 0.8516334, 0.2540831, 0.9464118, 0.8425044)
)

# Create the cell proportions data frame
cell_proportions <- data.frame(
  Sample_ID = c("A1", "A2"),
  B = c(0.0329, 0.0396),
  NK = c(0.0248, 0.0215),
  CD4T = c(0.164, 0.0723),
  CD8T = c(0.0779, 0.0380),
  Mono = c(0.0777, 0.0887),
  Neutro = c(0.622, 0.740)
)

# Function to adjust beta values for a specific sample
adjust_beta_values <- function(sample_name) {
  # Extract the beta values for the specified sample
  beta_values_sample <- betas_nc[[sample_name]]
  
  # Get the corresponding cell proportions for the specified sample
  cell_props_sample <- ed_nc %>%
    filter(Sample_ID == sample_name) %>%
    select(-Sample_ID,-LH_smoke_status)  # Remove Sample_ID for modeling
  
  # Create a data frame for modeling
  model_data <- data.frame(beta = beta_values_sample, cell_props_sample)
  
  # Fit a linear model for the current sample
  model <- lm(beta ~ ., data = model_data)
  
  # Calculate adjusted values
  adjusted_values <- predict(model)  # Use the model to predict adjusted values
  
  # Return adjusted values
  return(adjusted_values)
}

# Adjust beta values for sample A1
adjusted_A1 <- adjust_beta_values("A1")

# Adjust beta values for sample A2
adjusted_A2 <- adjust_beta_values("A2")

# Combine adjusted values into a data frame
adjusted_betas <- data.frame(
  CpG = beta_values$CpG,
  adjusted_A1 = adjusted_A1,
  adjusted_A2 = adjusted_A2
)

# View the adjusted beta values
print(adjusted_betas)

pdf(file = "output/dmc/never_current_estimate_dmc.pdf", width=12, height=4)
for(i in 2:(length(colnames(ncbetas))-7)){
  B <- ggplot(ncbetas, aes(x=as.numeric(B), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "B")
  NK <- ggplot(ncbetas, aes(x=as.numeric(NK), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "NK")
  CD4T <- ggplot(ncbetas, aes(x=as.numeric(CD4T), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "CD4T")
  CD8T<- ggplot(ncbetas, aes(x=as.numeric(CD8T), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "CD8T")
  Mono <- ggplot(ncbetas, aes(x=as.numeric(Mono), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "Mono")
  Neutro<- ggplot(ncbetas, aes(x=as.numeric(Neutro), y=as.numeric(!!sym(colnames(ncbetas)[i])), shape=LH_smoke_status, color=LH_smoke_status)) +
    geom_point() + geom_smooth(aes(group=LH_smoke_status),method='loess') + labs(x = colnames(ncbetas)[i],y = "Neutro")
  p <- plot_grid(B,NK,CD4T,CD8T,Mono,Neutro, labels = "AUTO",scale = 1)
  print(p)
}
dev.off()






write.table(cdmc_nc$dmct,sep='\t',file= "output/dmc/cell_dmc_never_current.tsv",quote = F)
write.table(cdmc_fc$dmct,sep='\t',file= "output/dmc/cell_dmc_former_current.tsv",quote = F)

wb <- wb_workbook()
wb <- wb$add_worksheet("DMCT")
wb$add_data("DMCT", cdmc_fc$dmct,row_names = T)
wb$add_worksheet("B")
wb$add_data("B", cdmc_fc$coe$B,row_names = T)
wb$add_worksheet("NK")
wb$add_data("NK", cdmc_fc$coe$NK, row_names = T)
wb$add_worksheet("CD4T")
wb$add_data("CD4T", cdmc_fc$coe$CD4T,row_names = T)
wb$add_worksheet("CD8T")
wb$add_data("CD8T", cdmc_fc$coe$CD8T,row_names = T)
wb$add_worksheet("Mono")
wb$add_data("Mono", cdmc_fc$coe$Mono,row_names = T)
wb$add_worksheet("Neutro")
wb$add_data("Neutro", cdmc_fc$coe$Neutro, row_names = T)

wb_save(wb,file = 'output/dmc/former_current_dmc.xlsx')

