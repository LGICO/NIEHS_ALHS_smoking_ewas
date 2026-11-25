library(dplyr)
library(data.table)

load('751_clean_pc.rdata')

output_dir <- "output/table1"
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

## Standardize labels
table1_data <- pheno %>% mutate(
  Smoking_group=case_when(
    as.character(LH_smoke_status) %in% c("2", "Current") ~ "Current",
    as.character(LH_smoke_status) %in% c("1", "Former") ~ "Former",
    as.character(LH_smoke_status) %in% c("0", "Never") ~ "Never",
    TRUE ~ "Missing/unknown"
  ),
  Gender_table=case_when(
    tolower(as.character(Gender)) %in% c("1", "female", "f") ~ "Female",
    tolower(as.character(Gender)) %in% c("0", "male", "m") ~ "Male",
    TRUE ~ "Missing/unknown"
  ),
  State_table=case_when(
    as.character(STATE) %in% c("0", "IA", "Iowa") ~ "Iowa",
    as.character(STATE) %in% c("1", "NC", "North Carolina") ~ "North Carolina",
    TRUE ~ "Missing/unknown"
  ),
  Screening_table=case_when(
    as.character(LHSSCREEN_12142015) %in% c("1", "Case") ~ "Case",
    as.character(LHSSCREEN_12142015) %in% c("0", "Control") ~ "Control",
    TRUE ~ "Missing/unknown"
  ),
  Ancestry_table=case_when(
    as.character(European) %in% c("1", "European") ~ "European",
    as.character(European) %in% c("0", "Non_European", "Non-European") ~ "Non-European",
    TRUE ~ "Missing/unknown"
  ),
  ANYETS_table=case_when(
    as.character(ANYETS) %in% c("0", "No") ~ "No",
    as.character(ANYETS) %in% c("1", "Yes") ~ "Yes",
    TRUE ~ "Missing/unknown"
  )
)

group_order <- c("Current", "Former", "Never", "Missing/unknown")
test_data <- table1_data %>% filter(Smoking_group %in% c("Current", "Former", "Never"))
test_data$Smoking_group <- factor(test_data$Smoking_group, levels=c("Never", "Former", "Current"))

format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

continuous_info <- data.frame(
  Variable=c("LH_BMI", "AgeAtHomeVisit", "B", "NK", "CD4T", "CD8T", "Mono", "Neutro"),
  Characteristic=c("BMI, kg/m²", "Age at home visit, years", "% B cells", "% NK cells",
                   "% CD4 T cells", "% CD8 T cells", "% monocytes", "% neutrophils"),
  Digits=c(2, 1, 2, 2, 2, 2, 2, 2),
  Multiply=c(1, 1, 100, 100, 100, 100, 100, 100),
  stringsAsFactors=FALSE
)

categorical_info <- list(
  "Sex"=list(variable="Gender_table", levels=c("Female", "Male")),
  "State"=list(variable="State_table", levels=c("Iowa", "North Carolina")),
  "Asthma screening status"=list(variable="Screening_table", levels=c("Case", "Control")),
  "Ancestry group"=list(variable="Ancestry_table", levels=c("European", "Non-European")),
  "Secondhand smoke exposure"=list(variable="ANYETS_table", levels=c("No", "Yes", "Missing/unknown"))
)

## Continuous-variable rows
continuous_rows <- bind_rows(lapply(seq_len(nrow(continuous_info)), function(i) {
  v <- continuous_info$Variable[i]
  label <- continuous_info$Characteristic[i]
  digits <- continuous_info$Digits[i]
  multiplier <- continuous_info$Multiply[i]
  
  values <- sapply(group_order, function(g) {
    x <- table1_data[table1_data$Smoking_group == g, ][[v]] * multiplier
    sprintf(
      paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"),
      median(x, na.rm=TRUE),
      quantile(x, 0.25, na.rm=TRUE),
      quantile(x, 0.75, na.rm=TRUE)
    )
  })
  
  x_all <- table1_data[[v]] * multiplier
  overall <- sprintf(
    paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"),
    median(x_all, na.rm=TRUE),
    quantile(x_all, 0.25, na.rm=TRUE),
    quantile(x_all, 0.75, na.rm=TRUE)
  )
  
  p <- kruskal.test(test_data[[v]] ~ test_data$Smoking_group)$p.value
  
  data.frame(
    Characteristic=label,
    Level="",
    Current=values["Current"],
    Former=values["Former"],
    Never=values["Never"],
    `Missing/unknown`=values["Missing/unknown"],
    Overall=overall,
    Overall_P=format_p(p),
    check.names=FALSE
  )
}))

## Categorical-variable rows
categorical_rows <- bind_rows(lapply(names(categorical_info), function(label) {
  v <- categorical_info[[label]]$variable
  levels_to_show <- categorical_info[[label]]$levels
  
  test_subset <- test_data
  if(v == "ANYETS_table") test_subset <- test_subset %>% filter(ANYETS_table %in% c("No", "Yes"))
  tab_test <- table(test_subset$Smoking_group, test_subset[[v]], useNA="no")
  expected <- suppressWarnings(chisq.test(tab_test)$expected)
  overall_test <- if(any(expected < 5)) fisher.test(tab_test) else chisq.test(tab_test, correct=FALSE)
  overall_p <- format_p(overall_test$p.value)
  
  bind_rows(lapply(seq_along(levels_to_show), function(j) {
    level <- levels_to_show[j]
    
    values <- sapply(group_order, function(g) {
      group_values <- table1_data[table1_data$Smoking_group == g, ][[v]]
      n <- sum(group_values == level, na.rm=TRUE)
      denominator <- length(group_values)
      sprintf("%d (%.1f%%)", n, 100 * n / denominator)
    })
    
    all_values <- table1_data[[v]]
    overall_n <- sum(all_values == level, na.rm=TRUE)
    
    data.frame(
      Characteristic=ifelse(j == 1, label, ""),
      Level=level,
      Current=values["Current"],
      Former=values["Former"],
      Never=values["Never"],
      `Missing/unknown`=values["Missing/unknown"],
      Overall=sprintf("%d (%.1f%%)", overall_n, 100 * overall_n / nrow(table1_data)),
      Overall_P=ifelse(j == 1, overall_p, ""),
      check.names=FALSE
    )
  }))
}))

## Arrange the final Table 1
table1_final <- bind_rows(
  continuous_rows %>% filter(Characteristic %in% c("BMI, kg/m²", "Age at home visit, years")),
  categorical_rows,
  continuous_rows %>% filter(!Characteristic %in% c("BMI, kg/m²", "Age at home visit, years"))
)

## Add sample sizes to column headings
group_n <- table(table1_data$Smoking_group)

names(table1_final)[names(table1_final) == "Current"] <- paste0("Current (n=", group_n["Current"], ")")
names(table1_final)[names(table1_final) == "Former"] <- paste0("Former (n=", group_n["Former"], ")")
names(table1_final)[names(table1_final) == "Never"] <- paste0("Never (n=", group_n["Never"], ")")
names(table1_final)[names(table1_final) == "Missing/unknown"] <- paste0("Missing/unknown (n=", group_n["Missing/unknown"], ")")
names(table1_final)[names(table1_final) == "Overall"] <- paste0("Overall (N=", nrow(table1_data), ")")

## Pairwise comparisons
comparisons <- list(c("Never", "Former"), c("Never", "Current"), c("Former", "Current"))

continuous_pairwise <- bind_rows(lapply(seq_len(nrow(continuous_info)), function(i) {
  v <- continuous_info$Variable[i]
  label <- continuous_info$Characteristic[i]
  
  bind_rows(lapply(comparisons, function(groups) {
    dat <- test_data %>% filter(Smoking_group %in% groups)
    if(v == "ANYETS_table") dat <- dat %>% filter(ANYETS_table %in% c("No", "Yes"))
    
    data.frame(
      Characteristic=label,
      Comparison=paste(groups, collapse=" vs "),
      Test="Wilcoxon rank-sum",
      P=wilcox.test(dat[[v]] ~ droplevels(dat$Smoking_group), exact=FALSE)$p.value
    )
  }))
}))

categorical_pairwise <- bind_rows(lapply(names(categorical_info), function(label) {
  v <- categorical_info[[label]]$variable
  
  bind_rows(lapply(comparisons, function(groups) {
    dat <- test_data %>% filter(Smoking_group %in% groups)
    tab <- table(droplevels(dat$Smoking_group), dat[[v]], useNA="no")
    expected <- suppressWarnings(chisq.test(tab)$expected)
    test <- if(any(expected < 5)) fisher.test(tab) else chisq.test(tab, correct=FALSE)
    
    data.frame(
      Characteristic=label,
      Comparison=paste(groups, collapse=" vs "),
      Test=ifelse(any(expected < 5), "Fisher exact", "Chi-squared"),
      P=test$p.value
    )
  }))
}))

pairwise_tests <- bind_rows(continuous_pairwise, categorical_pairwise) %>%
  group_by(Characteristic) %>%
  mutate(P_Holm=p.adjust(P, method="holm"), P_Holm_formatted=format_p(P_Holm)) %>%
  ungroup()

## Save results
fwrite(table1_final, file.path(output_dir, "Table1_final.tsv"), sep="\t")
fwrite(pairwise_tests, file.path(output_dir, "Table1_pairwise_comparisons.tsv"), sep="\t")

if(requireNamespace("openxlsx", quietly=TRUE)) {
  openxlsx::write.xlsx(
    list(
      `Table 1`=table1_final,
      `Pairwise comparisons`=pairwise_tests
    ),
    file.path(output_dir, "Table1_final.xlsx"),
    overwrite=TRUE
  )
}

print(table1_final, row.names=FALSE)
print(pairwise_tests, n=Inf)
