#!/usr/bin/env Rscript

# Extract MSA SNP genotypes and inferred sex from paired IDAT files.
#
# Usage:
#   Rscript extract_msa_idat_calls.R SAMPLE_SHEET.csv IDAT_DIR OUT_DIR \
#       [SAMPLE_ID_COLUMN] [RECORDED_SEX_COLUMN] [GENOME]
#
# GENOME defaults to hg38. SAMPLE_ID_COLUMN defaults to Sample_Name, and
# RECORDED_SEX_COLUMN defaults to Sex. The sample sheet must contain either a
# Basename column or both Sentrix_ID and Sentrix_Position.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 6L) {
  stop(paste(
    "Usage: Rscript extract_msa_idat_calls.R SAMPLE_SHEET.csv IDAT_DIR OUT_DIR",
    "[SAMPLE_ID_COLUMN] [RECORDED_SEX_COLUMN] [GENOME]"
  ))
}

sample_sheet_path <- normalizePath(args[[1]], mustWork = TRUE)
idat_dir <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- args[[3]]
sample_id_col <- if (length(args) >= 4L) args[[4]] else "Sample_Name"
recorded_sex_col <- if (length(args) >= 5L) args[[5]] else "Sex"
genome <- if (length(args) >= 6L) args[[6]] else "hg38"

if (!requireNamespace("sesame", quietly = TRUE)) {
  stop("The Bioconductor package 'sesame' is required.")
}
if (!requireNamespace("sesameData", quietly = TRUE)) {
  stop("The Bioconductor package 'sesameData' is required.")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
vcf_dir <- file.path(out_dir, "msa_vcf")
dir.create(vcf_dir, recursive = TRUE, showWarnings = FALSE)

read_illumina_sheet <- function(path) {
  lines <- readLines(path, warn = FALSE)
  data_marker <- which(trimws(lines) == "[Data]")
  if (length(data_marker) > 0L) {
    text <- paste(lines[(data_marker[[1]] + 1L):length(lines)], collapse = "\n")
    read.csv(text = text, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

normalise_sex <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[y %in% c("M", "MALE", "1", "XY")] <- "M"
  y[y %in% c("F", "FEMALE", "2", "XX")] <- "F"
  y[!y %in% c("M", "F")] <- "UNKNOWN"
  y
}

safe_name <- function(x) {
  y <- gsub("[^A-Za-z0-9_.-]", "_", as.character(x))
  ifelse(nchar(y) == 0L, "unnamed", y)
}

sheet <- read_illumina_sheet(sample_sheet_path)
required <- c(sample_id_col, recorded_sex_col)
missing_cols <- setdiff(required, names(sheet))
if (length(missing_cols) > 0L) {
  stop("Missing sample-sheet column(s): ", paste(missing_cols, collapse = ", "))
}
if (!("Basename" %in% names(sheet)) &&
    !all(c("Sentrix_ID", "Sentrix_Position") %in% names(sheet))) {
  stop("Sample sheet needs Basename, or both Sentrix_ID and Sentrix_Position.")
}

prefix_for_row <- function(i) {
  if ("Basename" %in% names(sheet) && nzchar(trimws(sheet[["Basename"]][[i]]))) {
    b <- sheet[["Basename"]][[i]]
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", b)) b <- file.path(idat_dir, b)
    return(sub("_(Grn|Red)\\.idat(\\.gz)?$", "", b, ignore.case = TRUE))
  }
  file.path(idat_dir, paste0(
    sheet[["Sentrix_ID"]][[i]], "_", sheet[["Sentrix_Position"]][[i]]
  ))
}

annotation_name <- paste0("MSA.", genome, ".snp.tsv.gz")
message("Loading SeSAMe SNP annotation: ", annotation_name)
annotation_path <- 'references/MSA.hg38.snp.tsv.gz'
anno <- utils::read.delim(
  gzfile(annotation_path), check.names = FALSE, stringsAsFactors = FALSE
)

qc <- vector("list", nrow(sheet))
for (i in seq_len(nrow(sheet))) {
  sample_id <- as.character(sheet[[sample_id_col]][[i]])
  prefix <- prefix_for_row(i)
  out_vcf <- file.path(vcf_dir, paste0(safe_name(sample_id), ".vcf"))
  message(sprintf("[%d/%d] %s", i, nrow(sheet), sample_id))

  result <- tryCatch({
    sdf <- sesame::readIDATpair(prefix, platform = "MSA")
    # Sex inference consumes methylation beta values. Genotype calling below
    # deliberately uses the raw SigDF because channel inference can erase
    # informative Type-I SNP signal.
    betas <- sesame::openSesame(sdf)
    inferred <- as.character(sesame::inferSex(betas, platform = "MSA"))
    sesame::formatVCF(sdf, anno = anno, vcf = out_vcf, genome = genome)
    data.frame(
      msa_id = sample_id,
      recorded_sex = normalise_sex(sheet[[recorded_sex_col]][[i]]),
      msa_inferred_sex = normalise_sex(inferred),
      sex_concordant = normalise_sex(sheet[[recorded_sex_col]][[i]]) ==
        normalise_sex(inferred),
      vcf_path = out_vcf,
      status = "OK",
      error = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      msa_id = sample_id,
      recorded_sex = normalise_sex(sheet[[recorded_sex_col]][[i]]),
      msa_inferred_sex = "UNKNOWN",
      sex_concordant = FALSE,
      vcf_path = out_vcf,
      status = "ERROR",
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  qc[[i]] <- result
}

qc_df <- do.call(rbind, qc)
utils::write.table(
  qc_df, file.path(out_dir, "msa_sex_qc.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA"
)
message("Wrote ", file.path(out_dir, "msa_sex_qc.tsv"))
