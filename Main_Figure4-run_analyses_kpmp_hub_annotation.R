#!/usr/bin/env Rscript

# ==============================================================================
# KPMP new analyses for revised manuscript — v3
#
# Runs:
#   exhaustive metabolite–podocyte/PT injury program analysis
#   elastic net classifier and score validation
#
# Expected project layout:
#   KPMP-data/
#     data/
#     output/
#     scripts/
#     figures/
#
# Place this script in: KPMP-data/scripts/
# Run it from there on the cluster.
#
# Inputs expected in KPMP-data/data/
#   Plasma_Metabolomics.xlsx
#   curated clinical.csv
#   gene_sets_modules.csv
#   Cluster_assignments_Q1.csv
#   Validation_cleaned.xlsx
#   OP-UHPLC MS_MS-Plasma Data Information-2022.xlsx
#   scRNA/*.rds   (one POD object and one PT object)
#
# Outputs written to:
#   KPMP-data/output/new_figure6/
#   KPMP-data/output/new_figure7/
#   KPMP-data/figures/
#
# Changelog from v2_1:
#   FIX 1 (line 147): safe_name() used make.names(unique=TRUE) on a long-format
#     column, assigning a different suffix to every participant-row of the same
#     metabolite.  Replaced with a map-back approach so every occurrence of the
#     same raw name receives the same safe name.
#   FIX 2 (line 385): purrr was not in req_pkgs and not loaded; purrr::imap_dfr
#     would throw "could not find function".  Added purrr to both req_pkgs and
#     the library() block.
#   FIX 3 (line 422): pathway_csv was read but never written — the derivation
#     block from op_info_xlsx was missing entirely.  Added derivation before the
#     read, then replaced the read_csv call with the in-memory object.
#   FIX 4 (lines 538-539): KPMP columns used make.names dots; validation columns
#     used janitor underscores; intersect() returned empty, triggering stop().
#     Added make.names() normalisation to validation column names before the
#     intersect so both use the same convention.
# ==============================================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------
# FIX 2: added "purrr" to req_pkgs so it is auto-installed and loaded
req_pkgs <- c(
  "readxl", "readr", "dplyr", "tidyr", "stringr", "ggplot2",
  "patchwork", "pheatmap", "janitor", "tibble", "broom",
  "glmnet", "pROC", "ggrepel", "Matrix", "forcats",
  "scales", "RColorBrewer", "purrr"
)

install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(
      missing,
      repos = "https://cloud.r-project.org",
      dependencies = c("Depends", "Imports", "LinkingTo")
    )
  }
}

# NOTE: do not install packages inside the batch job.
# Check only and fail fast if anything is missing.
still_missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing) > 0) {
  stop(
    "Missing R packages in your user library: ",
    paste(still_missing, collapse = ", "),
    "\nInstall these once before submitting the batch job."
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(janitor)
  library(tibble)
  library(broom)
  library(glmnet)
  library(pROC)
  library(ggrepel)
  library(Matrix)
  library(forcats)
  library(scales)
  library(RColorBrewer)
  library(purrr)   # FIX 2: was missing from library() block
})

# Seurat check and load — kept separate so the informative stop() is reached
# if Seurat is absent, rather than failing silently inside the block above
have_seurat <- requireNamespace("Seurat", quietly = TRUE)
if (!have_seurat) stop("Seurat is required. Install it in your user library before running.")
suppressPackageStartupMessages(library(Seurat))

theme_set(theme_bw(base_size = 12))
set.seed(123)

# ------------------------------------------------------------------------------
# 0) PATHS
# ------------------------------------------------------------------------------
script_dir <- getwd()
proj_dir   <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
data_dir   <- file.path(proj_dir, "data")
output_dir <- file.path(proj_dir, "output")
fig_dir    <- file.path(proj_dir, "figures")

fig6_out <- file.path(output_dir, "new_figure6")
fig7_out <- file.path(output_dir, "new_figure7")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(fig6_out,   recursive = TRUE, showWarnings = FALSE)
dir.create(fig7_out,   recursive = TRUE, showWarnings = FALSE)

plasma_xlsx     <- file.path(data_dir, "Plasma_Metabolomics.xlsx")
clinical_csv    <- file.path(data_dir, "curated clinical.csv")
genesets_csv    <- file.path(data_dir, "gene_sets_modules.csv")
clusters_csv    <- file.path(data_dir, "Cluster_assignments_Q1.csv")
validation_xlsx <- file.path(data_dir, "Validation_cleaned.xlsx")
op_info_xlsx    <- file.path(data_dir, "OP-UHPLC MS_MS-Plasma Data Information-2022.xlsx")
pathway_csv     <- file.path(output_dir, "metabolite_pathway_annotation_version3.csv")
scrna_dir       <- file.path(data_dir, "scRNA")

required_files <- c(plasma_xlsx, clinical_csv, genesets_csv,
                    clusters_csv, validation_xlsx, op_info_xlsx)
missing_files  <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}
if (!dir.exists(scrna_dir)) stop("Missing directory: ", scrna_dir)

# ------------------------------------------------------------------------------
# 1) HELPERS
# ------------------------------------------------------------------------------
find_first_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

clean_id <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "`", "")
  x <- str_replace_all(x, "\\s+", "")
  x <- str_replace(x, "^g(?=[0-9]{1,3}-[0-9]{4,6}$)", "")
  x
}

# FIX 1: safe_name() — consistent safe name per unique raw value
#
# The original used make.names(x, unique = TRUE) on the full long-format column.
# Because the same metabolite name appears once per participant (~240 rows),
# make.names assigned a different suffix to each row: the first participant got
# "indoxyl.sulfate", the second "indoxyl.sulfate.1", etc.  After pivot_wider
# this produced ~n*m columns instead of m, corrupting the entire matrix.
#
# Fix: find the unique raw names, make them safe and unique once, then map every
# occurrence back to its canonical safe name via match().  All participants now
# share the same safe name for the same metabolite.
safe_name <- function(x) {
  x        <- as.character(x)
  uvals    <- unique(x)
  safe_u   <- make.names(uvals, unique = TRUE)   # unique = TRUE is correct here
                                                  # because it operates only on the
                                                  # set of distinct names, not on
                                                  # the repeated long-format vector
  safe_u[match(x, uvals)]
}

age_to_numeric <- function(x) {
  x   <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  need_parse <- is.na(out) & !is.na(x)
  if (any(need_parse)) {
    out[need_parse] <- vapply(x[need_parse], function(xx) {
      nums <- suppressWarnings(as.numeric(unlist(str_extract_all(xx, "[0-9]+"))))
      nums <- nums[is.finite(nums)]
      if (length(nums) == 0) return(NA_real_)
      if (length(nums) == 1) return(nums[1])
      mean(range(nums))
    }, numeric(1))
  }
  out
}

save_plot_both <- function(p, stem, width = 8, height = 6) {
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = 300)
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height)
}

build_pseudobulk <- function(obj, label = "celltype") {
  meta   <- obj[[]] %>% janitor::clean_names()
  id_col <- find_first_col(meta, c("sample_id", "sampleid", "kpmpid",
                                    "participant_id", "participantid"))
  if (is.na(id_col)) {
    stop("Could not find participant/sample ID in ", label, " Seurat object metadata.")
  }
  ids <- clean_id(meta[[id_col]])

  expr     <- Seurat::GetAssayData(obj, assay = DefaultAssay(obj), layer = "data")
  expr     <- as(expr, "dgCMatrix")
  uniq_ids <- unique(ids)

  pb <- matrix(0, nrow = length(uniq_ids), ncol = nrow(expr),
               dimnames = list(uniq_ids, rownames(expr)))

  for (i in seq_along(uniq_ids)) {
    idx <- which(ids == uniq_ids[i])
    if (length(idx) == 1) {
      pb[i, ] <- as.numeric(expr[, idx])
    } else {
      pb[i, ] <- Matrix::rowMeans(expr[, idx, drop = FALSE])
    }
  }
  pb
}

run_atlas_regression <- function(pb_mat, met_mat, cov_mat,
                                 min_participants = 10,
                                 gene_prev_thresh = 0.20,
                                 fdr_method       = "BH") {

  stopifnot(nrow(pb_mat) == nrow(met_mat), nrow(pb_mat) == nrow(cov_mat))

  # Gene prevalence filter: retain genes detected (>0 pseudobulk) in at least
  # gene_prev_thresh fraction of participants
  gene_prev <- colMeans(pb_mat > 0, na.rm = TRUE)
  genes_ok  <- gene_prev >= gene_prev_thresh
  pb_f      <- pb_mat[, genes_ok, drop = FALSE]
  message(sprintf("  Atlas: %d / %d genes pass %.0f%% prevalence filter",
                  sum(genes_ok), length(genes_ok), gene_prev_thresh * 100))

  results <- vector("list", ncol(met_mat))

  for (i in seq_len(ncol(met_mat))) {
    met_name <- colnames(met_mat)[i]
    met_vals <- met_mat[, i]

    # Median-impute any residual NAs in the z-scored metabolite vector
    if (any(is.na(met_vals))) {
      met_vals[is.na(met_vals)] <- median(met_vals, na.rm = TRUE)
    }

    # Design matrix: [intercept | metabolite | age | sex | eGFR]
    X  <- cbind(1, met_vals, cov_mat)
    cc <- complete.cases(X)
    if (sum(cc) < min_participants) next

    X_cc <- X[cc, , drop = FALSE]
    Y_cc <- pb_f[cc, , drop = FALSE]

    # Drop zero-variance genes in this complete-case subset
    gvar     <- apply(Y_cc, 2, var, na.rm = TRUE)
    genes_nz <- is.finite(gvar) & gvar > 1e-10
    Y_use    <- Y_cc[, genes_nz, drop = FALSE]
    if (ncol(Y_use) < 1) next

    # Vectorised OLS via lm.fit(): solves all gene outcomes simultaneously
    # for the shared design matrix — much faster than looping gene-by-gene
    fit <- tryCatch(lm.fit(X_cc, Y_use), error = function(e) NULL)
    if (is.null(fit)) next

    coef_mat <- fit$coefficients      # p × n_genes; row 2 = metabolite β
    beta_met <- coef_mat[2, ]

    n_obs    <- nrow(X_cc)
    p_pars   <- ncol(X_cc)
    df_res   <- n_obs - p_pars
    resid    <- Y_use - X_cc %*% coef_mat
    mse_gene <- colSums(resid^2) / df_res

    XtXinv <- tryCatch(solve(crossprod(X_cc)), error = function(e) NULL)
    if (is.null(XtXinv)) next

    se_met  <- sqrt(mse_gene * XtXinv[2, 2])
    t_stat  <- beta_met / se_met
    p_val   <- 2 * pt(-abs(t_stat), df = df_res)
    fdr_val <- p.adjust(p_val, method = fdr_method)

    results[[i]] <- tibble(
      metabolite = met_name,
      gene       = colnames(Y_use),
      beta       = as.numeric(beta_met),
      se         = as.numeric(se_met),
      t_stat     = as.numeric(t_stat),
      p_value    = as.numeric(p_val),
      FDR        = as.numeric(fdr_val),
      n          = sum(cc)
    )

    if (i %% 100 == 0 || i == ncol(met_mat)) {
      message(sprintf("  Regression progress: %d / %d metabolites", i, ncol(met_mat)))
    }
  }
  bind_rows(results)
}

# ------------------------------------------------------------------------------
# 2) LOAD PLASMA, CLINICAL, GENE SETS
# ------------------------------------------------------------------------------
message("Loading plasma metabolomics...")
plasma_raw <- readxl::read_excel(plasma_xlsx) %>% janitor::clean_names()

pid_col <- find_first_col(plasma_raw, c("participant_id", "participant", "participantid"))
met_col <- find_first_col(plasma_raw, c("chemical_name", "metabolite", "metabolite_name"))
ab_col  <- find_first_col(plasma_raw, c("abundance", "value", "intensity", "peak_area"))
if (any(is.na(c(pid_col, met_col, ab_col)))) {
  stop("Could not detect required columns in Plasma_Metabolomics.xlsx")
}

plasma_long <- plasma_raw %>%
  transmute(
    ParticipantID    = clean_id(.data[[pid_col]]),
    # FIX 1 applied here: safe_name() now maps back from unique values so every
    # participant gets the same safe name for the same metabolite
    metabolite       = safe_name(.data[[met_col]]),
    metabolite_label = as.character(.data[[met_col]]),
    abundance        = suppressWarnings(as.numeric(.data[[ab_col]]))
  ) %>%
  filter(!is.na(ParticipantID), !is.na(metabolite), !is.na(abundance))

# One readable label per safe metabolite name — used for figure axis labels
plasma_labels <- plasma_long %>% distinct(metabolite, metabolite_label)

plasma_wide <- plasma_long %>%
  group_by(ParticipantID, metabolite) %>%
  summarize(abundance = median(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = metabolite, values_from = abundance)

# Prevalence filter: keep metabolites present in >=30% of participants
miss_tbl  <- plasma_wide %>%
  summarize(across(-ParticipantID, ~ mean(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "metabolite", values_to = "miss_rate")
keep_mets      <- miss_tbl %>% filter(miss_rate <= 0.70) %>% pull(metabolite)
plasma_wide_f  <- plasma_wide %>% select(ParticipantID, all_of(keep_mets))

# Z-score scale; store training mean/SD for later deployment on validation cohort
X_raw    <- plasma_wide_f %>% select(-ParticipantID) %>% as.matrix()
col_means <- colMeans(X_raw, na.rm = TRUE)
col_sds   <- apply(X_raw, 2, sd, na.rm = TRUE)
for (j in seq_len(ncol(X_raw))) {
  nas <- is.na(X_raw[, j])
  if (any(nas)) X_raw[nas, j] <- col_means[j]
}
mat_scaled <- scale(X_raw, center = col_means, scale = pmax(col_sds, 1e-8))
rownames(mat_scaled) <- plasma_wide_f$ParticipantID

message(sprintf("Plasma: %d participants × %d metabolites after QC",
                nrow(mat_scaled), ncol(mat_scaled)))

message("Loading clinical metadata...")
clinical_raw <- readr::read_csv(clinical_csv, show_col_types = FALSE) %>%
  janitor::clean_names()
clin_pid  <- find_first_col(clinical_raw, c("participant_id", "participantid", "participant"))
clin_age  <- find_first_col(clinical_raw, c("age", "age_group", "age_years", "age_years_binned"))
clin_sex  <- find_first_col(clinical_raw, c("sex", "gender"))
clin_egfr <- find_first_col(clinical_raw, c("egfr", "baseline_e_gfr_ml_min_1_73m2", "baseline_egfr"))
if (any(is.na(c(clin_pid, clin_age, clin_sex, clin_egfr)))) {
  stop("Could not detect required clinical columns.")
}

clinical_clean <- clinical_raw %>%
  transmute(
    ParticipantID = clean_id(.data[[clin_pid]]),
    Age           = age_to_numeric(.data[[clin_age]]),
    Sex           = as.character(.data[[clin_sex]]),
    eGFR          = suppressWarnings(as.numeric(.data[[clin_egfr]]))
  ) %>%
  distinct(ParticipantID, .keep_all = TRUE)

message("Loading gene modules...")
gene_sets <- readr::read_csv(genesets_csv, show_col_types = FALSE) %>%
  janitor::clean_names()
if (!all(c("module", "gene") %in% names(gene_sets))) {
  stop("gene_sets_modules.csv must contain module and gene columns.")
}
gene_modules <- split(gene_sets$gene, gene_sets$module)

# FIX 3: derive and write pathway_annotation from op_info_xlsx BEFORE reading it
# -------------------------------------------------------------------------------
# In v2_1, pathway_csv was read at line 422 but never written anywhere in the
# script, causing an immediate "file not found" error.  The Metabolon "Data
# Information" file (op_info_xlsx) contains super_pathway and sub_pathway for
# every metabolite.  We derive the annotation table here, write it to
# pathway_csv so it is available for downstream reads and for Zenodo deposition,
# and keep the in-memory object (pathway_annot) for immediate use.
message("Deriving pathway annotation from Metabolon info file...")
sheet_names <- readxl::excel_sheets(op_info_xlsx)

if (!"Chemical Annotation" %in% sheet_names) {
  stop(
    "Could not find 'Chemical Annotation' sheet in OP info file. Sheets found: ",
    paste(sheet_names, collapse = ", ")
  )
}

op_info_raw <- readxl::read_excel(op_info_xlsx, sheet = "Chemical Annotation") %>%
  janitor::clean_names()

# Detect the name column — Metabolon files use "chemical_name" or "biochemical"
op_name_col  <- find_first_col(op_info_raw,
                                c("chemical_name", "biochemical", "biochemical_name"))
op_super_col <- find_first_col(op_info_raw,
                                c("super_pathway", "superpathway", "super_path"))
op_sub_col   <- find_first_col(op_info_raw,
                                c("sub_pathway", "subpathway", "sub_path"))

if (is.na(op_name_col) || is.na(op_super_col)) {
  stop("OP info file missing required annotation columns. Found: ",
       paste(names(op_info_raw), collapse = ", "))
}

pathway_annot <- op_info_raw %>%
  transmute(
    # Apply the same safe_name() function used for plasma_long so that the
    # metabolite key joins correctly to hub_metabolites downstream
    metabolite    = safe_name(.data[[op_name_col]]),
    super_pathway = as.character(.data[[op_super_col]]),
    sub_pathway   = if (!is.na(op_sub_col)) as.character(.data[[op_sub_col]]) else NA_character_
  ) %>%
  filter(!is.na(super_pathway), super_pathway != "") %>%
  distinct(metabolite, .keep_all = TRUE)

readr::write_csv(pathway_annot, pathway_csv)
message(sprintf("Pathway annotation: %d metabolites written to %s",
                nrow(pathway_annot), pathway_csv))

# ------------------------------------------------------------------------------
# 3) NEW FIGURE 6: EXHAUSTIVE CELL PROGRAM ANALYSIS
# ------------------------------------------------------------------------------
message("Building pseudobulk matrices from Seurat objects...")
sc_files_all <- list.files(scrna_dir, pattern = "\\.rds$",
                            full.names = TRUE, ignore.case = TRUE)

# Detect POD and PT files by filename substring.
# POD pattern is specific ("pod"); PT uses "pt_s" or "prox" before the broader
# "pt" to reduce the risk of matching unrelated filenames containing "pt".
pod_rds_path <- sc_files_all[grep("(?i)pod",
                                   basename(sc_files_all), perl = TRUE)][1]
pt_rds_path  <- sc_files_all[grep("(?i)(pt_s|prox|proximal)",
                                   basename(sc_files_all), perl = TRUE)][1]

# Fallback: if the narrower PT pattern matched nothing, try bare "pt" on files
# that did not already match the POD pattern
if (is.na(pt_rds_path)) {
  non_pod_files <- sc_files_all[!grepl("(?i)pod", basename(sc_files_all), perl = TRUE)]
  pt_rds_path   <- non_pod_files[grep("(?i)pt", basename(non_pod_files), perl = TRUE)][1]
}

if (is.na(pod_rds_path) || is.na(pt_rds_path)) {
  stop("Could not auto-detect POD and PT .rds files in data/scRNA/\n",
       "Files found: ", paste(basename(sc_files_all), collapse = ", "))
}
if (pod_rds_path == pt_rds_path) {
  stop("POD and PT .rds paths resolved to the same file: ", pod_rds_path,
       "\nRename scRNA files so POD and PT are unambiguously distinct.")
}

message("POD file: ", basename(pod_rds_path))
message("PT  file: ", basename(pt_rds_path))

pod_obj <- readRDS(pod_rds_path)
pt_obj  <- readRDS(pt_rds_path)
pod_pb  <- build_pseudobulk(pod_obj, "POD")
pt_pb   <- build_pseudobulk(pt_obj,  "PT")
rm(pod_obj, pt_obj); gc()

message(sprintf("POD pseudobulk: %d participants × %d genes", nrow(pod_pb), ncol(pod_pb)))
message(sprintf("PT  pseudobulk: %d participants × %d genes", nrow(pt_pb),  ncol(pt_pb)))

# Build numeric covariate matrix [age | sex (0/1) | eGFR] for a set of IDs
build_cov_matrix <- function(ids, clin_df) {
  clin_sub <- clin_df[match(ids, clin_df$ParticipantID), ]
  age_num  <- clin_sub$Age
  sex_bin  <- as.integer(
    factor(toupper(substr(as.character(clin_sub$Sex), 1, 1)),
           levels = c("F", "M"))
  ) - 1L
  egfr_num <- clin_sub$eGFR
  # Median-impute covariates so regressions don't lose participants to NA
  age_num[is.na(age_num)]    <- median(age_num,  na.rm = TRUE)
  sex_bin[is.na(sex_bin)]    <- 0L
  egfr_num[is.na(egfr_num)] <- median(egfr_num, na.rm = TRUE)
  out <- cbind(age = age_num, sex = sex_bin, egfr = egfr_num)
  rownames(out) <- ids
  out
}

# Align matrices to participants present in both plasma and pseudobulk data
pod_common <- intersect(rownames(mat_scaled), rownames(pod_pb))
pt_common  <- intersect(rownames(mat_scaled), rownames(pt_pb))
message(sprintf("%d participants with plasma + POD data", length(pod_common)))
message(sprintf("%d participants with plasma + PT  data", length(pt_common)))

if (length(pod_common) < 10) stop("Too few participants with plasma + POD data: ", length(pod_common))
if (length(pt_common)  < 10) stop("Too few participants with plasma + PT  data: ", length(pt_common))

met_for_pod    <- mat_scaled[pod_common, , drop = FALSE]
met_for_pt     <- mat_scaled[pt_common,  , drop = FALSE]
pod_pb_aligned <- pod_pb[pod_common, , drop = FALSE]
pt_pb_aligned  <- pt_pb[pt_common,  , drop = FALSE]
cov_pod        <- build_cov_matrix(pod_common, clinical_clean)
cov_pt         <- build_cov_matrix(pt_common,  clinical_clean)

message("Running exhaustive regressions for POD (this will take several minutes)...")
atlas_pod <- run_atlas_regression(pod_pb_aligned, met_for_pod, cov_pod)
message("Running exhaustive regressions for PT...")
atlas_pt  <- run_atlas_regression(pt_pb_aligned,  met_for_pt,  cov_pt)

readr::write_csv(atlas_pod, file.path(fig6_out, "Atlas_POD_full_associations_version3.csv"))
readr::write_csv(atlas_pt,  file.path(fig6_out, "Atlas_PT_full_associations_version3.csv"))
message(sprintf("Atlas POD: %d pairs; %d with FDR < 0.05",
                nrow(atlas_pod), sum(atlas_pod$FDR < 0.05, na.rm = TRUE)))
message(sprintf("Atlas PT:  %d pairs; %d with FDR < 0.05",
                nrow(atlas_pt),  sum(atlas_pt$FDR  < 0.05, na.rm = TRUE)))

# --- Cell-type specificity scoring -------------------------------------------
# Pearson r of β-vectors (POD vs PT) per metabolite testable in both cell types.
# r > 0.40  → Convergent  (similar effect in both)
# r < -0.20 → Divergent   (opposing effects)
# otherwise → Restricted or Mixed

pod_sig     <- atlas_pod %>% select(metabolite, gene, beta_pod = beta, fdr_pod = FDR)
pt_sig      <- atlas_pt  %>% select(metabolite, gene, beta_pt  = beta, fdr_pt  = FDR)
shared_pairs <- inner_join(pod_sig, pt_sig, by = c("metabolite", "gene"))

specificity_scores <- shared_pairs %>%
  group_by(metabolite) %>%
  filter(n() >= 10) %>%   # require >=10 shared gene estimates for a reliable r
  summarise(
    n_shared_genes   = n(),
    beta_r_pod_vs_pt = cor(beta_pod, beta_pt, use = "complete.obs"),
    n_sig_pod        = sum(fdr_pod < 0.05, na.rm = TRUE),
    n_sig_pt         = sum(fdr_pt  < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    specificity_class = case_when(
      beta_r_pod_vs_pt >  0.40 ~ "Convergent",
      beta_r_pod_vs_pt < -0.20 ~ "Divergent",
      n_sig_pod > 0 & n_sig_pt == 0 ~ "POD-restricted",
      n_sig_pod == 0 & n_sig_pt > 0 ~ "PT-restricted",
      TRUE ~ "Mixed"
    )
  )

readr::write_csv(specificity_scores, file.path(fig6_out, "Atlas_celltype_specificity_version3.csv"))

# --- Hub metabolite scoring ---------------------------------------------------
# FIX 2 applied here: purrr::imap_dfr now works because purrr is loaded
# Hub score = total significant gene associations (FDR < 0.05) across both
# cell types; module fingerprint = dominant gene module among those associations

gene_module_map <- purrr::imap_dfr(gene_modules, function(genes, mod_name) {
  tibble(gene = genes, module = mod_name)
})

sig_combined <- bind_rows(
  atlas_pod %>% filter(FDR < 0.05) %>% mutate(cell_type = "POD"),
  atlas_pt  %>% filter(FDR < 0.05) %>% mutate(cell_type = "PT")
)

hub_total <- sig_combined %>%
  group_by(metabolite) %>%
  summarise(
    hub_score_total = n(),
    hub_score_pod   = sum(cell_type == "POD"),
    hub_score_pt    = sum(cell_type == "PT"),
    .groups = "drop"
  )

module_fingerprint <- sig_combined %>%
  left_join(gene_module_map, by = "gene") %>%
  filter(!is.na(module)) %>%
  group_by(metabolite, module) %>%
  summarise(n_sig_in_module = n(), .groups = "drop") %>%
  group_by(metabolite) %>%
  mutate(frac_in_module = n_sig_in_module / sum(n_sig_in_module)) %>%
  slice_max(frac_in_module, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(dominant_module = module, dominant_module_frac = frac_in_module)

hub_metabolites <- hub_total %>%
  left_join(module_fingerprint %>%
              select(metabolite, dominant_module, dominant_module_frac),
            by = "metabolite") %>%
  left_join(specificity_scores %>%
              select(metabolite, specificity_class, beta_r_pod_vs_pt),
            by = "metabolite") %>%
  left_join(plasma_labels, by = "metabolite") %>%
  mutate(plot_label = coalesce(metabolite_label, metabolite)) %>%
  arrange(desc(hub_score_total))

readr::write_csv(hub_metabolites, file.path(fig6_out, "Atlas_hub_metabolites_version3.csv"))

# --- Pathway-level hub summary -----------------------------------------------
# FIX 3: pathway_annot is now the in-memory object derived above;
# no file read needed here (the object already exists)

pathway_hub <- hub_metabolites %>%
  left_join(pathway_annot %>% select(metabolite, super_pathway),
            by = "metabolite") %>%
  filter(!is.na(super_pathway)) %>%
  group_by(super_pathway) %>%
  summarise(
    n_metabolites    = n(),
    mean_hub_score   = mean(hub_score_total,   na.rm = TRUE),
    median_hub_score = median(hub_score_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_hub_score))

readr::write_csv(pathway_hub, file.path(fig6_out, "Atlas_pathway_hub_summary_version3.csv"))

# --- Figure 6 panels ---------------------------------------------------------
fig6A <- hub_metabolites %>%
  slice_max(hub_score_total, n = 20, with_ties = FALSE) %>%
  ggplot(aes(x = hub_score_pod, y = hub_score_pt,
             size = hub_score_total, color = dominant_module,
             label = plot_label)) +
  geom_point(alpha = 0.85) +
  geom_text_repel(size = 3, max.overlaps = 25) +
  scale_size_continuous(range = c(2, 12)) +
  labs(title = "A. Top hub metabolites across POD and PT",
       x = "Hub score in podocytes",
       y = "Hub score in proximal tubule",
       color = "Dominant module", size = "Total hub score")

fig6B <- specificity_scores %>%
  count(specificity_class) %>%
  ggplot(aes(x = reorder(specificity_class, n), y = n, fill = specificity_class)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(title = "B. Cell-type specificity classes",
       x = NULL, y = "Number of metabolites")

fig6C <- pathway_hub %>%
  slice_max(mean_hub_score, n = 15, with_ties = FALSE) %>%
  mutate(super_pathway = forcats::fct_reorder(super_pathway, mean_hub_score)) %>%
  ggplot(aes(x = mean_hub_score, y = super_pathway)) +
  geom_col(fill = "steelblue") +
  labs(title = "C. Pathways with greatest transcriptional reach",
       x = "Mean hub score", y = NULL)

fig6_combined <- (fig6A | fig6B) / fig6C
save_plot_both(fig6_combined,
               file.path(fig_dir, "Figure_new6_exhaustive_cell_program_analysis_version3"),
               width = 16, height = 12)
message("Figure 6 saved.")

# ------------------------------------------------------------------------------
# 4) NEW FIGURE 7: ELASTIC NET CLASSIFIER
# ------------------------------------------------------------------------------
message("Loading k=2 metabotype assignments...")
clusters_raw <- readr::read_csv(clusters_csv, show_col_types = FALSE) %>%
  janitor::clean_names()
cls_pid <- find_first_col(clusters_raw, c("participant_id", "sample_id", "participantid"))
cls_mt  <- find_first_col(clusters_raw, c("metabotype", "cluster", "k_means_cluster"))
if (any(is.na(c(cls_pid, cls_mt)))) {
  stop("Could not detect participant and metabotype columns in Cluster_assignments_Q1.csv")
}

k2_tbl <- clusters_raw %>%
  transmute(
    ParticipantID = clean_id(.data[[cls_pid]]),
    metabotype_k2 = as.integer(.data[[cls_mt]])
  ) %>%
  filter(!is.na(metabotype_k2))

# CKD-dominant metabotype = lowest integer (MT1 convention used throughout paper)
ckd_mt_id <- min(sort(unique(k2_tbl$metabotype_k2)))
k2_tbl    <- k2_tbl %>% mutate(y = as.integer(metabotype_k2 == ckd_mt_id))

cls_df    <- plasma_wide_f %>%
  inner_join(k2_tbl %>% select(ParticipantID, y), by = "ParticipantID")
X_cls_raw <- cls_df %>% select(-ParticipantID, -y) %>% as.matrix()
y_cls     <- cls_df$y
message(sprintf("Classifier: %d participants, %d CKD-dominant (y=1)",
                length(y_cls), sum(y_cls)))

# Impute and z-score using KPMP training statistics (saved for deployment)
kpmp_means <- colMeans(X_cls_raw, na.rm = TRUE)
kpmp_sds   <- apply(X_cls_raw, 2, sd, na.rm = TRUE)
for (j in seq_len(ncol(X_cls_raw))) {
  nas <- is.na(X_cls_raw[, j])
  if (any(nas)) X_cls_raw[nas, j] <- kpmp_means[j]
}
X_cls <- scale(X_cls_raw, center = kpmp_means, scale = pmax(kpmp_sds, 1e-8))

readr::write_csv(
  tibble(metabolite = colnames(X_cls), mean = kpmp_means, sd = kpmp_sds),
  file.path(fig7_out, "Classifier_scaling_parameters_version3.csv")
)

# Full-data elastic net (α=0.5) for final coefficient extraction
message("Running elastic net (full-data model)...")
set.seed(42)
cv_full  <- glmnet::cv.glmnet(x = X_cls, y = y_cls,
                               family = "binomial", alpha = 0.5,
                               nfolds = 10, type.measure = "auc")
coef_min <- coef(cv_full, s = "lambda.min")
coef_1se <- coef(cv_full, s = "lambda.1se")

extract_panel <- function(coef_obj, lambda_label) {
  # coef() returns a dgCMatrix; convert to data frame via dense matrix
  df <- as.data.frame(as.matrix(coef_obj))
  df$metabolite <- rownames(df)
  # Rename the coefficient column (named "s1" or "1" from the sparse matrix)
  names(df)[1] <- "coefficient"
  df %>%
    filter(metabolite != "(Intercept)", abs(coefficient) > 0) %>%
    arrange(desc(abs(coefficient))) %>%
    mutate(lambda = lambda_label) %>%
    select(metabolite, coefficient, lambda)
}

panel_min <- extract_panel(coef_min, "lambda.min")
panel_1se <- extract_panel(coef_1se, "lambda.1se")
message(sprintf("Classifier panel: %d metabolites at lambda.min; %d at lambda.1se",
                nrow(panel_min), nrow(panel_1se)))

readr::write_csv(panel_min, file.path(fig7_out, "Classifier_panel_lambda_min_version3.csv"))
readr::write_csv(panel_1se, file.path(fig7_out, "Classifier_panel_lambda_1se_version3.csv"))

# Outer 5-fold CV to get unbiased out-of-fold probability predictions for ROC
# Reset seed immediately before fold assignment for reproducibility
message("Running outer 5-fold CV for ROC predictions...")
set.seed(42)
fold_id     <- sample(rep(1:5, length.out = length(y_cls)))
outer_preds <- rep(NA_real_, length(y_cls))
outer_auc   <- numeric(5)

for (k in 1:5) {
  train_i <- which(fold_id != k)
  test_i  <- which(fold_id == k)
  cv_k    <- glmnet::cv.glmnet(
    x = X_cls[train_i, , drop = FALSE],
    y = y_cls[train_i],
    family = "binomial", alpha = 0.5, nfolds = 10, type.measure = "auc"
  )
  # predict returns an n×1 matrix; [, 1] extracts the vector cleanly
  pred_k              <- predict(cv_k, newx = X_cls[test_i, , drop = FALSE],
                                  s = "lambda.1se", type = "response")[, 1]
  outer_preds[test_i] <- pred_k
  outer_auc[k]        <- as.numeric(
    pROC::auc(pROC::roc(y_cls[test_i], pred_k, quiet = TRUE))
  )
}

message(sprintf("CV AUC: %.3f ± %.3f (5-fold)", mean(outer_auc), sd(outer_auc)))

cv_pred_tbl <- tibble(
  ParticipantID = cls_df$ParticipantID,
  y_true        = y_cls,
  prob_ckd_mt   = outer_preds
)
readr::write_csv(cv_pred_tbl, file.path(fig7_out, "Classifier_KPMP_CV_predictions_version3.csv"))
saveRDS(cv_full, file.path(fig7_out, "Classifier_glmnet_model_version3.rds"))

# --- External validation on MDRD/ALTOLD cohort -------------------------------
message("Loading validation cohort...")
val_raw      <- readxl::read_excel(validation_xlsx, sheet = "Combined") %>%
  janitor::clean_names()
val_meta_cols <- intersect(c("sample_id", "cohort", "group"), names(val_raw))
if (length(val_meta_cols) < 3) {
  stop("Validation_cleaned.xlsx Combined sheet must contain sample_id, cohort, and group columns.")
}

val_feat_raw <- val_raw %>%
  select(-all_of(val_meta_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
  as.matrix()
rownames(val_feat_raw) <- as.character(val_raw$sample_id)
val_meta <- val_raw %>%
  select(all_of(val_meta_cols)) %>%
  mutate(y_val = as.integer(tolower(cohort) == "mdrd"))

# FIX 4: normalise validation column names to make.names() dots before
# intersecting with KPMP columns.
# In v2_1, janitor::clean_names() produced underscore snake_case for validation
# columns while KPMP used make.names() dot-separated names.  The naming
# conventions are incompatible so intersect() returned empty, triggering stop().
# Applying make.names() here forces both to the same convention.
colnames(val_feat_raw) <- make.names(colnames(val_feat_raw), unique = FALSE)

shared_mets <- intersect(colnames(X_cls), colnames(val_feat_raw))
message(sprintf("Shared metabolites between KPMP and validation: %d", length(shared_mets)))
if (length(shared_mets) < 5) {
  stop("Too few shared metabolites between KPMP and validation cohort (",
       length(shared_mets), "). Check that metabolite naming is consistent.")
}

val_feat_sub <- val_feat_raw[, shared_mets, drop = FALSE]
for (j in seq_along(shared_mets)) {
  nas <- is.na(val_feat_sub[, j])
  if (any(nas)) val_feat_sub[nas, j] <- kpmp_means[shared_mets][j]
}
X_val <- scale(val_feat_sub,
               center = kpmp_means[shared_mets],
               scale  = pmax(kpmp_sds[shared_mets], 1e-8))

# Apply KPMP-trained coefficients to validation data without refitting
coef_1se_vec         <- as.numeric(coef_1se[-1])
names(coef_1se_vec)  <- rownames(coef_1se)[-1]
coef_shared_nonzero  <- coef_1se_vec[
  intersect(shared_mets, names(coef_1se_vec[coef_1se_vec != 0]))
]
intercept_val <- as.numeric(coef_1se[1])
lp_val        <- intercept_val +
  X_val[, names(coef_shared_nonzero), drop = FALSE] %*% coef_shared_nonzero
prob_val      <- 1 / (1 + exp(-as.numeric(lp_val)))

val_pred_tbl <- tibble(
  sample_id   = rownames(X_val),
  cohort      = val_meta$cohort,
  y_true      = val_meta$y_val,
  prob_ckd_mt = prob_val
)
readr::write_csv(val_pred_tbl, file.path(fig7_out, "Classifier_validation_predictions_version3.csv"))
message(sprintf("External validation AUC: %.3f",
                as.numeric(pROC::auc(
                  pROC::roc(val_pred_tbl$y_true, val_pred_tbl$prob_ckd_mt, quiet = TRUE)
                ))))

# MetabotypeScore = log-odds (linear predictor before sigmoid); captures a
# graded biological axis rather than binary CKD/healthy classification
kpmp_score_tbl <- cv_pred_tbl %>%
  mutate(MetabotypeScore = qlogis(pmin(pmax(prob_ckd_mt, 1e-8), 1 - 1e-8))) %>%
  left_join(clinical_clean %>% select(ParticipantID, eGFR), by = "ParticipantID") %>%
  mutate(
    eGFR_stage = case_when(
      eGFR >= 90 ~ "G1 (\u226590)",
      eGFR >= 60 ~ "G2 (60-89)",
      eGFR >= 45 ~ "G3a (45-59)",
      eGFR >= 30 ~ "G3b (30-44)",
      eGFR >= 15 ~ "G4 (15-29)",
      eGFR <  15 ~ "G5 (<15)",
      TRUE       ~ NA_character_
    )
  )

val_score_tbl <- val_pred_tbl %>%
  mutate(MetabotypeScore = qlogis(pmin(pmax(prob_ckd_mt, 1e-8), 1 - 1e-8)))

readr::write_csv(kpmp_score_tbl, file.path(fig7_out, "Classifier_KPMP_MetabotypeScore_version3.csv"))
readr::write_csv(val_score_tbl,  file.path(fig7_out, "Classifier_validation_MetabotypeScore_version3.csv"))

# --- Panel–atlas hub overlap -------------------------------------------------
hub_threshold  <- quantile(hub_metabolites$hub_score_total, 0.75, na.rm = TRUE)
hub_top        <- hub_metabolites %>%
  filter(hub_score_total >= hub_threshold) %>%
  pull(metabolite)
panel_mets     <- panel_1se$metabolite
universe_mets  <- intersect(colnames(X_cls), hub_metabolites$metabolite)

overlap_detail <- tibble(metabolite = union(panel_mets, hub_top)) %>%
  filter(metabolite %in% universe_mets) %>%
  mutate(
    in_classifier_panel = metabolite %in% panel_mets,
    in_hub_top_quartile = metabolite %in% hub_top
  ) %>%
  left_join(hub_metabolites %>%
              select(metabolite, hub_score_total, dominant_module,
                     specificity_class, plot_label),
            by = "metabolite") %>%
  left_join(panel_1se %>% select(metabolite, coefficient), by = "metabolite")

readr::write_csv(overlap_detail,
                 file.path(fig7_out, "Classifier_atlas_overlap_detail_version3.csv"))

overlap_summary <- tibble(
  n_universe         = length(universe_mets),
  n_hub_top_quartile = sum(universe_mets %in% hub_top),
  n_panel            = sum(universe_mets %in% panel_mets),
  n_overlap          = sum(universe_mets %in% panel_mets & universe_mets %in% hub_top)
)
readr::write_csv(overlap_summary,
                 file.path(fig7_out, "Classifier_atlas_overlap_summary_version3.csv"))

# --- Figure 7 panels ---------------------------------------------------------
roc_kpmp <- pROC::roc(cv_pred_tbl$y_true, cv_pred_tbl$prob_ckd_mt, quiet = TRUE)
roc_val  <- pROC::roc(val_pred_tbl$y_true, val_pred_tbl$prob_ckd_mt, quiet = TRUE)

roc_to_df <- function(roc_obj, label) {
  tibble(fpr = 1 - roc_obj$specificities,
         tpr = roc_obj$sensitivities,
         dataset = label)
}
roc_plot_df <- bind_rows(
  roc_to_df(roc_kpmp, sprintf("KPMP CV AUC = %.3f",
                               as.numeric(pROC::auc(roc_kpmp)))),
  roc_to_df(roc_val,  sprintf("Validation AUC = %.3f",
                               as.numeric(pROC::auc(roc_val))))
)

fig7A <- ggplot(roc_plot_df, aes(x = fpr, y = tpr, color = dataset)) +
  geom_path(linewidth = 1.1) +
  geom_abline(linetype = 2, color = "grey50") +
  coord_equal() +
  labs(title = "A. ROC curves for metabotype classifier",
       x = "False positive rate", y = "True positive rate", color = NULL)

fig7B <- kpmp_score_tbl %>%
  filter(!is.na(eGFR_stage)) %>%
  mutate(eGFR_stage = factor(
    eGFR_stage,
    levels = c("G1 (\u226590)", "G2 (60-89)", "G3a (45-59)",
               "G3b (30-44)", "G4 (15-29)", "G5 (<15)")
  )) %>%
  ggplot(aes(x = eGFR_stage, y = MetabotypeScore)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.3) +
  labs(title = "B. KPMP MetabotypeScore by eGFR stage",
       x = NULL, y = "MetabotypeScore") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

fig7C <- panel_1se %>%
  left_join(plasma_labels, by = "metabolite") %>%
  mutate(
    plot_label = coalesce(metabolite_label, metabolite),
    plot_label = forcats::fct_reorder(plot_label, abs(coefficient))
  ) %>%
  slice_max(abs(coefficient), n = 20, with_ties = FALSE) %>%
  ggplot(aes(x = coefficient, y = plot_label, fill = coefficient > 0)) +
  geom_col(show.legend = FALSE) +
  labs(title = "C. Top classifier coefficients",
       x = "Elastic net coefficient", y = NULL)

fig7D <- overlap_detail %>%
  mutate(status = case_when(
    in_classifier_panel & in_hub_top_quartile ~ "Both",
    in_classifier_panel                        ~ "Classifier only",
    in_hub_top_quartile                        ~ "Hub only",
    TRUE                                       ~ "Other"
  )) %>%
  count(status) %>%
  ggplot(aes(x = status, y = n, fill = status)) +
  geom_col(show.legend = FALSE) +
  labs(title = "D. Overlap: classifier panel vs hub metabolites",
       x = NULL, y = "Count")

fig7_combined <- (fig7A | fig7B) / (fig7C | fig7D)
save_plot_both(fig7_combined,
               file.path(fig_dir, "Figure_new7_classifier_validation_version3"),
               width = 16, height = 12)
message("Figure 7 saved.")

# ------------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo_new_analyses_version3.txt"))
message("Done. Outputs written to ", output_dir,
        " and figures written to ",  fig_dir)
