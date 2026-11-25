ss <- fread('samplesheets/LHS_Sample_Sheet.csv')
load('751_clean_pc.rdata')

table(ss$LH_smoke_status)
table(ss$LH_smoke_status)/887*100
ss %>% filter(SM3Mo != "" & LH_smoke_status == "Never") %>% dim
ss %>% filter(SM3Mo != "" & LH_smoke_status == "Former") %>% dim
ss %>% filter(SM3Mo != "" & LH_smoke_status == "Current") %>% dim
ss %>% filter(SM3Mo != "" & LH_smoke_status == "M") %>% dim

ss %>% filter(SM3Yr != "" & LH_smoke_status == "Never") %>% dim
ss %>% filter(SM3Yr != "" & LH_smoke_status == "Former") %>% dim
ss %>% filter(SM3Yr != "" & LH_smoke_status == "Current") %>% dim
ss %>% filter(SM3Yr != "" & LH_smoke_status == "M") %>% dim

table(pheno$LH_smoke_status)
table(pheno$LH_smoke_status)/764*100
pheno %>% filter(SM3Mo != "" & LH_smoke_status == 0) %>% dim
pheno %>% filter(SM3Mo != "" & LH_smoke_status == 1) %>% dim
pheno %>% filter(SM3Mo != "" & LH_smoke_status == 2) %>% dim
pheno %>% filter(SM3Mo != "" & LH_smoke_status == 99) %>% dim

pheno %>% filter(SM3Yr != "" & LH_smoke_status == 0) %>% dim
pheno %>% filter(SM3Yr != "" & LH_smoke_status == 1) %>% dim
pheno %>% filter(SM3Yr != "" & LH_smoke_status == 2) %>% dim
pheno %>% filter(SM3Yr != "" & LH_smoke_status == 99) %>% dim


pheno %>% dplyr::filter(LH_smoke_status ==0) %>% pull(LH_BMI) %>% mean()
pheno %>% dplyr::filter(LH_smoke_status ==0) %>% pull(LH_BMI) %>% sd()
pheno %>% dplyr::filter(LH_smoke_status ==1) %>% pull(LH_BMI) %>% mean()
pheno %>% dplyr::filter(LH_smoke_status ==1) %>% pull(LH_BMI) %>% sd()
pheno %>% dplyr::filter(LH_smoke_status ==2) %>% pull(LH_BMI) %>% mean()
pheno %>% dplyr::filter(LH_smoke_status ==2) %>% pull(LH_BMI) %>% sd()
pheno %>% dplyr::filter(LH_smoke_status ==99) %>% pull(LH_BMI) %>% mean()
pheno %>% dplyr::filter(LH_smoke_status ==99) %>% pull(LH_BMI) %>% sd()

pheno %>% dplyr::filter(LH_smoke_status ==0) %>% pull(LHSSCREEN_12142015) %>% table()
pheno %>% dplyr::filter(LH_smoke_status ==1) %>% pull(LHSSCREEN_12142015) %>% table()
pheno %>% dplyr::filter(LH_smoke_status ==2) %>% pull(LHSSCREEN_12142015) %>% table()
pheno %>% dplyr::filter(LH_smoke_status ==99) %>% pull(LHSSCREEN_12142015) %>% table()

pheno %>% dplyr::filter(LH_smoke_status ==0) %>% pull(STATE) %>% table()
pheno %>% dplyr::filter(LH_smoke_status ==1) %>% pull(STATE) %>% table()
pheno %>% dplyr::filter(LH_smoke_status ==2) %>% pull(STATE) %>% table()
pheno %>% dplyr::filter(LH_smoke_status ==99) %>% pull(STATE) %>% table()
