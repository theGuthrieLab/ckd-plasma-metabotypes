#!/usr/bin/env Rscript

################################################################
# Figure 5: Metabolite-associated podocyte and proximal tubule
#           transcriptional programs in KPMP
#
# 5A: Top hub metabolites across podocyte and proximal tubule
#     transcriptional association networks
# 5B: Cell-type specificity classes for metabolite-associated
#     gene programs
# 5C: Metabolite pathways with greatest transcriptional reach
#
################################################################

options(stringsAsFactors = FALSE)

# ---- Packages -------------------------------------------------
req_pkgs <- c(
  "readxl", "readr", "dplyr", "tidyr", "stringr", "ggplot2",
  "patchwork", "janitor", "tibble", "broom", "ggrepel",
  "Matrix", "forcats", "scales", "RColorBrewer", "purrr"
)

missing_pkgs <- req_pkgs[
  !vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall these before running."
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
  library(janitor)
  library(tibble)
  library(broom)
  library(ggrepel)
  library(Matrix)
  library(forcats)
  library(scales)
  library(RColorBrewer)
  library(purrr)
})

have_seurat <- requireNamespace("Seurat", quietly = TRUE)

if (!have_seurat) {
  stop("Seurat is required. Install it in your user library before running.")
}

suppressPackageStartupMessages({
  library(Seurat)
})

theme_set(theme_bw(base_size = 12))
set.seed(123)

# ================================================================
# EDIT THIS ONE LINE to point to the project folder
# ================================================================
base_dir <- "."   # e.g. "D:/UC Berkeley/KPMP paper/repo/Figure5"

data_dir   <- file.path(base_dir, "data")
output_dir <- file.path(base_dir, "Figure5_output")
figure_dir <- file.path(base_dir, "figuresFinal")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Input files ----------------------------------------------
plasma_xlsx  <- file.path(data_dir, "Plasma_Metabolomics.xlsx")
clinical_csv <- file.path(data_dir, "curated clinical.csv")
genesets_csv <- file.path(data_dir, "gene_sets_modules.csv")
op_info_xlsx <- file.path(data_dir, "OP-UHPLC MS_MS-Plasma Data Information-2022.xlsx")
scrna_dir    <- file.path(data_dir, "scRNA")

pathway_csv <- file.path(output_dir, "metabolite_pathway_annotation_Figure5.csv")

required_files <- c(
  plasma_xlsx,
  clinical_csv,
  genesets_csv,
  op_info_xlsx
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

if (!dir.exists(scrna_dir)) {
  stop("Missing scRNA directory: ", scrna_dir)
}

# ---- Helper functions -----------------------------------------
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

safe_name <- function(x) {
  x <- as.character(x)
  uvals <- unique(x)
  safe_u <- make.names(uvals, unique = TRUE)
  safe_u[match(x, uvals)]
}

age_to_numeric <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))

  need_parse <- is.na(out) & !is.na(x)

  if (any(need_parse)) {
    out[need_parse] <- vapply(
      x[need_parse],
      function(xx) {
        nums <- suppressWarnings(
          as.numeric(unlist(str_extract_all(xx, "[0-9]+")))
        )

        nums <- nums[is.finite(nums)]

        if (length(nums) == 0) return(NA_real_)
        if (length(nums) == 1) return(nums[1])

        mean(range(nums))
      },
      numeric(1)
    )
  }

  out
}

save_plot_both <- function(plot, stem, width = 8, height = 6) {
  ggsave(
    paste0(stem, ".png"),
    plot,
    width = width,
    height = height,
    dpi = 300
  )

  ggsave(
    paste0(stem, ".pdf"),
    plot,
    width = width,
    height = height
  )
}

save_tiff <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(figure_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    compression = "lzw"
  )
}

build_pseudobulk <- function(obj, label = "celltype") {
  meta <- obj[[]] %>%
    janitor::clean_names()

  id_col <- find_first_col(
    meta,
    c(
      "sample_id",
      "sampleid",
      "kpmpid",
      "participant_id",
      "participantid"
    )
  )

  if (is.na(id_col)) {
    stop(
      "Could not find participant/sample ID in ",
      label,
      " Seurat object metadata."
    )
  }

  ids <- clean_id(meta[[id_col]])

  expr <- Seurat::GetAssayData(
    obj,
    assay = DefaultAssay(obj),
    layer = "data"
  )

  expr <- as(expr, "dgCMatrix")
  uniq_ids <- unique(ids)

  pb <- matrix(
    0,
    nrow = length(uniq_ids),
    ncol = nrow(expr),
    dimnames = list(uniq_ids, rownames(expr))
  )

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

build_cov_matrix <- function(ids, clin_df) {
  clin_sub <- clin_df[match(ids, clin_df$ParticipantID), ]

  age_num <- clin_sub$Age

  sex_bin <- as.integer(
    factor(
      toupper(substr(as.character(clin_sub$Sex), 1, 1)),
      levels = c("F", "M")
    )
  ) - 1L

  egfr_num <- clin_sub$eGFR

  age_num[is.na(age_num)] <- median(age_num, na.rm = TRUE)
  sex_bin[is.na(sex_bin)] <- 0L
  egfr_num[is.na(egfr_num)] <- median(egfr_num, na.rm = TRUE)

  out <- cbind(
    age = age_num,
    sex = sex_bin,
    egfr = egfr_num
  )

  rownames(out) <- ids

  out
}

run_atlas_regression <- function(
  pb_mat,
  met_mat,
  cov_mat,
  min_participants = 10,
  gene_prev_thresh = 0.20,
  fdr_method = "BH"
) {
  stopifnot(
    nrow(pb_mat) == nrow(met_mat),
    nrow(pb_mat) == nrow(cov_mat)
  )

  gene_prev <- colMeans(pb_mat > 0, na.rm = TRUE)
  genes_ok <- gene_prev >= gene_prev_thresh
  pb_f <- pb_mat[, genes_ok, drop = FALSE]

  message(
    sprintf(
      "  Atlas: %d / %d genes pass %.0f%% prevalence filter",
      sum(genes_ok),
      length(genes_ok),
      gene_prev_thresh * 100
    )
  )

  results <- vector("list", ncol(met_mat))

  for (i in seq_len(ncol(met_mat))) {
    met_name <- colnames(met_mat)[i]
    met_vals <- met_mat[, i]

    if (any(is.na(met_vals))) {
      met_vals[is.na(met_vals)] <- median(met_vals, na.rm = TRUE)
    }

    X <- cbind(1, met_vals, cov_mat)
    cc <- complete.cases(X)

    if (sum(cc) < min_participants) next

    X_cc <- X[cc, , drop = FALSE]
    Y_cc <- pb_f[cc, , drop = FALSE]

    gvar <- apply(Y_cc, 2, var, na.rm = TRUE)
    genes_nz <- is.finite(gvar) & gvar > 1e-10
    Y_use <- Y_cc[, genes_nz, drop = FALSE]

    if (ncol(Y_use) < 1) next

    fit <- tryCatch(
      lm.fit(X_cc, Y_use),
      error = function(e) NULL
    )

    if (is.null(fit)) next

    coef_mat <- fit$coefficients
    beta_met <- coef_mat[2, ]

    n_obs <- nrow(X_cc)
    p_pars <- ncol(X_cc)
    df_res <- n_obs - p_pars

    resid <- Y_use - X_cc %*% coef_mat
    mse_gene <- colSums(resid^2) / df_res

    XtXinv <- tryCatch(
      solve(crossprod(X_cc)),
      error = function(e) NULL
    )

    if (is.null(XtXinv)) next

    se_met <- sqrt(mse_gene * XtXinv[2, 2])
    t_stat <- beta_met / se_met
    p_val <- 2 * pt(-abs(t_stat), df = df_res)
    fdr_val <- p.adjust(p_val, method = fdr_method)

    results[[i]] <- tibble(
      metabolite = met_name,
      gene = colnames(Y_use),
      beta = as.numeric(beta_met),
      se = as.numeric(se_met),
      t_stat = as.numeric(t_stat),
      p_value = as.numeric(p_val),
      FDR = as.numeric(fdr_val),
      n = sum(cc)
    )

    if (i %% 100 == 0 || i == ncol(met_mat)) {
      message(
        sprintf(
          "  Regression progress: %d / %d metabolites",
          i,
          ncol(met_mat)
        )
      )
    }
  }

  bind_rows(results)
}

# ================================================================
# LOAD PLASMA METABOLOMICS
# ================================================================

message("Loading plasma metabolomics...")

plasma_raw <- readxl::read_excel(plasma_xlsx) %>%
  janitor::clean_names()

pid_col <- find_first_col(
  plasma_raw,
  c("participant_id", "participant", "participantid")
)

met_col <- find_first_col(
  plasma_raw,
  c("chemical_name", "metabolite", "metabolite_name")
)

ab_col <- find_first_col(
  plasma_raw,
  c("abundance", "value", "intensity", "peak_area")
)

if (any(is.na(c(pid_col, met_col, ab_col)))) {
  stop("Could not detect required columns in Plasma_Metabolomics.xlsx")
}

plasma_long <- plasma_raw %>%
  transmute(
    ParticipantID = clean_id(.data[[pid_col]]),
    metabolite = safe_name(.data[[met_col]]),
    metabolite_label = as.character(.data[[met_col]]),
    abundance = suppressWarnings(as.numeric(.data[[ab_col]]))
  ) %>%
  filter(
    !is.na(ParticipantID),
    !is.na(metabolite),
    !is.na(abundance)
  )

plasma_labels <- plasma_long %>%
  distinct(metabolite, metabolite_label)

plasma_wide <- plasma_long %>%
  group_by(ParticipantID, metabolite) %>%
  summarize(
    abundance = median(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = metabolite,
    values_from = abundance
  )

miss_tbl <- plasma_wide %>%
  summarize(across(-ParticipantID, ~ mean(is.na(.x)))) %>%
  pivot_longer(
    everything(),
    names_to = "metabolite",
    values_to = "miss_rate"
  )

keep_mets <- miss_tbl %>%
  filter(miss_rate <= 0.70) %>%
  pull(metabolite)

plasma_wide_f <- plasma_wide %>%
  select(ParticipantID, all_of(keep_mets))

X_raw <- plasma_wide_f %>%
  select(-ParticipantID) %>%
  as.matrix()

col_means <- colMeans(X_raw, na.rm = TRUE)
col_sds <- apply(X_raw, 2, sd, na.rm = TRUE)

for (j in seq_len(ncol(X_raw))) {
  nas <- is.na(X_raw[, j])

  if (any(nas)) {
    X_raw[nas, j] <- col_means[j]
  }
}

mat_scaled <- scale(
  X_raw,
  center = col_means,
  scale = pmax(col_sds, 1e-8)
)

rownames(mat_scaled) <- plasma_wide_f$ParticipantID

message(
  sprintf(
    "Plasma: %d participants x %d metabolites after QC",
    nrow(mat_scaled),
    ncol(mat_scaled)
  )
)

# ================================================================
# LOAD CLINICAL METADATA
# ================================================================

message("Loading clinical metadata...")

clinical_raw <- readr::read_csv(
  clinical_csv,
  show_col_types = FALSE
) %>%
  janitor::clean_names()

clin_pid <- find_first_col(
  clinical_raw,
  c("participant_id", "participantid", "participant")
)

clin_age <- find_first_col(
  clinical_raw,
  c("age", "age_group", "age_years", "age_years_binned")
)

clin_sex <- find_first_col(
  clinical_raw,
  c("sex", "gender")
)

clin_egfr <- find_first_col(
  clinical_raw,
  c("egfr", "baseline_e_gfr_ml_min_1_73m2", "baseline_egfr")
)

if (any(is.na(c(clin_pid, clin_age, clin_sex, clin_egfr)))) {
  stop("Could not detect required clinical columns.")
}

clinical_clean <- clinical_raw %>%
  transmute(
    ParticipantID = clean_id(.data[[clin_pid]]),
    Age = age_to_numeric(.data[[clin_age]]),
    Sex = as.character(.data[[clin_sex]]),
    eGFR = suppressWarnings(as.numeric(.data[[clin_egfr]]))
  ) %>%
  distinct(ParticipantID, .keep_all = TRUE)

# ================================================================
# LOAD GENE MODULES
# ================================================================

message("Loading gene modules...")

gene_sets <- readr::read_csv(
  genesets_csv,
  show_col_types = FALSE
) %>%
  janitor::clean_names()

if (!all(c("module", "gene") %in% names(gene_sets))) {
  stop("gene_sets_modules.csv must contain module and gene columns.")
}

gene_modules <- split(gene_sets$gene, gene_sets$module)

# ================================================================
# DERIVE METABOLITE PATHWAY ANNOTATION
# ================================================================

message("Deriving pathway annotation from Metabolon info file...")

sheet_names <- readxl::excel_sheets(op_info_xlsx)

if (!"Chemical Annotation" %in% sheet_names) {
  stop(
    "Could not find 'Chemical Annotation' sheet in OP info file. Sheets found: ",
    paste(sheet_names, collapse = ", ")
  )
}

op_info_raw <- readxl::read_excel(
  op_info_xlsx,
  sheet = "Chemical Annotation"
) %>%
  janitor::clean_names()

op_name_col <- find_first_col(
  op_info_raw,
  c("chemical_name", "biochemical", "biochemical_name")
)

op_super_col <- find_first_col(
  op_info_raw,
  c("super_pathway", "superpathway", "super_path")
)

op_sub_col <- find_first_col(
  op_info_raw,
  c("sub_pathway", "subpathway", "sub_path")
)

if (is.na(op_name_col) || is.na(op_super_col)) {
  stop(
    "OP info file missing required annotation columns. Found: ",
    paste(names(op_info_raw), collapse = ", ")
  )
}

pathway_annot <- op_info_raw %>%
  transmute(
    metabolite = safe_name(.data[[op_name_col]]),
    super_pathway = as.character(.data[[op_super_col]]),
    sub_pathway = if (!is.na(op_sub_col)) {
      as.character(.data[[op_sub_col]])
    } else {
      NA_character_
    }
  ) %>%
  filter(
    !is.na(super_pathway),
    super_pathway != ""
  ) %>%
  distinct(metabolite, .keep_all = TRUE)

readr::write_csv(pathway_annot, pathway_csv)

message(
  sprintf(
    "Pathway annotation: %d metabolites written to %s",
    nrow(pathway_annot),
    pathway_csv
  )
)

# ================================================================
# BUILD PSEUDOBULK MATRICES FROM SEURAT OBJECTS
# ================================================================

message("Building pseudobulk matrices from Seurat objects...")

sc_files_all <- list.files(
  scrna_dir,
  pattern = "\\.rds$",
  full.names = TRUE,
  ignore.case = TRUE
)

pod_rds_path <- sc_files_all[
  grep("(?i)pod", basename(sc_files_all), perl = TRUE)
][1]

pt_rds_path <- sc_files_all[
  grep("(?i)(pt_s|prox|proximal)", basename(sc_files_all), perl = TRUE)
][1]

if (is.na(pt_rds_path)) {
  non_pod_files <- sc_files_all[
    !grepl("(?i)pod", basename(sc_files_all), perl = TRUE)
  ]

  pt_rds_path <- non_pod_files[
    grep("(?i)pt", basename(non_pod_files), perl = TRUE)
  ][1]
}

if (is.na(pod_rds_path) || is.na(pt_rds_path)) {
  stop(
    "Could not auto-detect POD and PT .rds files in data/scRNA/\n",
    "Files found: ",
    paste(basename(sc_files_all), collapse = ", ")
  )
}

if (pod_rds_path == pt_rds_path) {
  stop("POD and PT .rds paths resolved to the same file: ", pod_rds_path)
}

message("POD file: ", basename(pod_rds_path))
message("PT file: ", basename(pt_rds_path))

pod_obj <- readRDS(pod_rds_path)
pt_obj <- readRDS(pt_rds_path)

pod_pb <- build_pseudobulk(pod_obj, "POD")
pt_pb <- build_pseudobulk(pt_obj, "PT")

rm(pod_obj, pt_obj)
gc()

message(
  sprintf(
    "POD pseudobulk: %d participants x %d genes",
    nrow(pod_pb),
    ncol(pod_pb)
  )
)

message(
  sprintf(
    "PT pseudobulk: %d participants x %d genes",
    nrow(pt_pb),
    ncol(pt_pb)
  )
)

# ================================================================
# ALIGN PLASMA, PSEUDOBULK, AND COVARIATES
# ================================================================

pod_common <- intersect(rownames(mat_scaled), rownames(pod_pb))
pt_common <- intersect(rownames(mat_scaled), rownames(pt_pb))

message(sprintf("%d participants with plasma + POD data", length(pod_common)))
message(sprintf("%d participants with plasma + PT data", length(pt_common)))

if (length(pod_common) < 10) {
  stop("Too few participants with plasma + POD data: ", length(pod_common))
}

if (length(pt_common) < 10) {
  stop("Too few participants with plasma + PT data: ", length(pt_common))
}

met_for_pod <- mat_scaled[pod_common, , drop = FALSE]
met_for_pt <- mat_scaled[pt_common, , drop = FALSE]

pod_pb_aligned <- pod_pb[pod_common, , drop = FALSE]
pt_pb_aligned <- pt_pb[pt_common, , drop = FALSE]

cov_pod <- build_cov_matrix(pod_common, clinical_clean)
cov_pt <- build_cov_matrix(pt_common, clinical_clean)

# ================================================================
# RUN EXHAUSTIVE METABOLITE-GENE REGRESSIONS
# ================================================================

message("Running exhaustive regressions for POD...")

atlas_pod <- run_atlas_regression(
  pod_pb_aligned,
  met_for_pod,
  cov_pod
)

message("Running exhaustive regressions for PT...")

atlas_pt <- run_atlas_regression(
  pt_pb_aligned,
  met_for_pt,
  cov_pt
)

readr::write_csv(
  atlas_pod,
  file.path(output_dir, "Figure5_POD_full_associations.csv")
)

readr::write_csv(
  atlas_pt,
  file.path(output_dir, "Figure5_PT_full_associations.csv")
)

message(
  sprintf(
    "Atlas POD: %d pairs; %d with FDR < 0.05",
    nrow(atlas_pod),
    sum(atlas_pod$FDR < 0.05, na.rm = TRUE)
  )
)

message(
  sprintf(
    "Atlas PT: %d pairs; %d with FDR < 0.05",
    nrow(atlas_pt),
    sum(atlas_pt$FDR < 0.05, na.rm = TRUE)
  )
)

# ================================================================
# CELL-TYPE SPECIFICITY SUMMARY
# ================================================================

pod_sig <- atlas_pod %>%
  select(
    metabolite,
    gene,
    beta_pod = beta,
    fdr_pod = FDR
  )

pt_sig <- atlas_pt %>%
  select(
    metabolite,
    gene,
    beta_pt = beta,
    fdr_pt = FDR
  )

shared_pairs <- inner_join(
  pod_sig,
  pt_sig,
  by = c("metabolite", "gene")
)

specificity_scores <- shared_pairs %>%
  group_by(metabolite) %>%
  filter(n() >= 10) %>%
  summarise(
    n_shared_genes = n(),
    beta_r_pod_vs_pt = cor(beta_pod, beta_pt, use = "complete.obs"),
    n_sig_pod = sum(fdr_pod < 0.05, na.rm = TRUE),
    n_sig_pt = sum(fdr_pt < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    specificity_class = case_when(
      beta_r_pod_vs_pt > 0.40 ~ "Convergent",
      beta_r_pod_vs_pt < -0.20 ~ "Divergent",
      n_sig_pod > 0 & n_sig_pt == 0 ~ "POD-restricted",
      n_sig_pod == 0 & n_sig_pt > 0 ~ "PT-restricted",
      TRUE ~ "Mixed"
    )
  )

readr::write_csv(
  specificity_scores,
  file.path(output_dir, "Figure5_celltype_specificity.csv")
)

# ================================================================
# HUB METABOLITE SUMMARY
# ================================================================

gene_module_map <- purrr::imap_dfr(
  gene_modules,
  function(genes, mod_name) {
    tibble(
      gene = genes,
      module = mod_name
    )
  }
)

sig_combined <- bind_rows(
  atlas_pod %>%
    filter(FDR < 0.05) %>%
    mutate(cell_type = "POD"),
  atlas_pt %>%
    filter(FDR < 0.05) %>%
    mutate(cell_type = "PT")
)

hub_total <- sig_combined %>%
  group_by(metabolite) %>%
  summarise(
    hub_score_total = n(),
    hub_score_pod = sum(cell_type == "POD"),
    hub_score_pt = sum(cell_type == "PT"),
    .groups = "drop"
  )

module_fingerprint <- sig_combined %>%
  left_join(gene_module_map, by = "gene") %>%
  filter(!is.na(module)) %>%
  group_by(metabolite, module) %>%
  summarise(
    n_sig_in_module = n(),
    .groups = "drop"
  ) %>%
  group_by(metabolite) %>%
  mutate(
    frac_in_module = n_sig_in_module / sum(n_sig_in_module)
  ) %>%
  slice_max(
    frac_in_module,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  rename(
    dominant_module = module,
    dominant_module_frac = frac_in_module
  )

hub_metabolites <- hub_total %>%
  left_join(
    module_fingerprint %>%
      select(
        metabolite,
        dominant_module,
        dominant_module_frac
      ),
    by = "metabolite"
  ) %>%
  left_join(
    specificity_scores %>%
      select(
        metabolite,
        specificity_class,
        beta_r_pod_vs_pt
      ),
    by = "metabolite"
  ) %>%
  left_join(plasma_labels, by = "metabolite") %>%
  mutate(
    plot_label = coalesce(metabolite_label, metabolite)
  ) %>%
  arrange(desc(hub_score_total))

readr::write_csv(
  hub_metabolites,
  file.path(output_dir, "Figure5_hub_metabolites.csv")
)

# ================================================================
# PATHWAY HUB SUMMARY
# ================================================================

pathway_hub <- hub_metabolites %>%
  left_join(
    pathway_annot %>%
      select(metabolite, super_pathway),
    by = "metabolite"
  ) %>%
  filter(!is.na(super_pathway)) %>%
  group_by(super_pathway) %>%
  summarise(
    n_metabolites = n(),
    mean_hub_score = mean(hub_score_total, na.rm = TRUE),
    median_hub_score = median(hub_score_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_hub_score))

readr::write_csv(
  pathway_hub,
  file.path(output_dir, "Figure5_pathway_hub_summary.csv")
)

# ================================================================
# PANEL 5A: Top hub metabolites across POD and PT
# ================================================================

fig5A <- hub_metabolites %>%
  slice_max(
    hub_score_total,
    n = 20,
    with_ties = FALSE
  ) %>%
  ggplot(
    aes(
      x = hub_score_pod,
      y = hub_score_pt,
      size = hub_score_total,
      color = dominant_module,
      label = plot_label
    )
  ) +
  geom_point(alpha = 0.85) +
  ggrepel::geom_text_repel(
    size = 3,
    max.overlaps = 25
  ) +
  scale_size_continuous(range = c(2, 12)) +
  labs(
    title = "A",
    x = "Hub score in podocytes",
    y = "Hub score in proximal tubule",
    color = "Dominant module",
    size = "Total hub score"
  ) +
  theme_cowplot() +
  panel_border(color = "black") +
  theme(legend.position = "right")

save_tiff(
  fig5A,
  "Fig5A_top_hub_metabolites_POD_PT.tiff",
  width = 7,
  height = 5.5
)

message("Saved Fig 5A")

# ================================================================
# PANEL 5B: Cell-type specificity classes
# ================================================================

fig5B <- specificity_scores %>%
  count(specificity_class) %>%
  ggplot(
    aes(
      x = reorder(specificity_class, n),
      y = n,
      fill = specificity_class
    )
  ) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(
    title = "B",
    x = NULL,
    y = "Number of metabolites"
  ) +
  theme_cowplot() +
  panel_border(color = "black")

save_tiff(
  fig5B,
  "Fig5B_celltype_specificity_classes.tiff",
  width = 5,
  height = 4
)

message("Saved Fig 5B")

# ================================================================
# PANEL 5C: Pathways with greatest transcriptional reach
# ================================================================

fig5C <- pathway_hub %>%
  slice_max(
    mean_hub_score,
    n = 15,
    with_ties = FALSE
  ) %>%
  mutate(
    super_pathway = forcats::fct_reorder(
      super_pathway,
      mean_hub_score
    )
  ) %>%
  ggplot(
    aes(
      x = mean_hub_score,
      y = super_pathway
    )
  ) +
  geom_col(fill = "steelblue") +
  labs(
    title = "C",
    x = "Mean hub score",
    y = NULL
  ) +
  theme_cowplot() +
  panel_border(color = "black")

save_tiff(
  fig5C,
  "Fig5C_pathways_transcriptional_reach.tiff",
  width = 6,
  height = 5
)

message("Saved Fig 5C")

# ================================================================
# OPTIONAL COMBINED FIGURE 5
# ================================================================

fig5_combined <- (fig5A | fig5B) / fig5C

save_tiff(
  fig5_combined,
  "Figure5_combined_panels_A_C.tiff",
  width = 12,
  height = 10
)

save_plot_both(
  fig5_combined,
  file.path(figure_dir, "Figure5_combined_panels_A_C"),
  width = 12,
  height = 10
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo_Figure5.txt")
)

message("\n========== Figure 5 complete ==========")
message("Outputs in: ", output_dir)
message("Figures in: ", figure_dir)
