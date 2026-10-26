#!/usr/bin/env Rscript

################################################################
# Figure 1: KPMP metabotype discovery, classifier performance,
#           and validation-cohort alignment
#
# 1B: KPMP k-means Hartigan-Wong metabotype discovery and PCA
#     colored by eGFR-relabelled metabotype
# 1C: Five-fold cross-validated elastic-net ROC curves for KPMP
#     and validation classifiers
# 1D: Cross-validated metabotype score by metabotype and cohort
# 1E: Top KPMP classifier metabolites from the lambda.1se model
# 1F: Top validation classifier metabolites from the lambda.1se model
# 1G: Concordance of top classifier coefficients between KPMP
#     and validation cohorts
#
################################################################

options(stringsAsFactors = FALSE)

# ---- Packages -------------------------------------------------
req_pkgs <- c(
  "readxl", "readr", "dplyr", "tidyr", "stringr", "janitor",
  "ggplot2", "patchwork", "glmnet", "pROC", "ggrepel",
  "tibble", "forcats", "cluster", "vegan", "cowplot"
)

missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(janitor)
  library(ggplot2)
  library(patchwork)
  library(glmnet)
  library(pROC)
  library(ggrepel)
  library(tibble)
  library(forcats)
  library(cluster)
  library(vegan)
  library(cowplot)
})

theme_set(theme_bw(base_size = 12))
set.seed(123)

# ================================================================
# EDIT THIS ONE LINE to point to the project folder
# ================================================================
base_dir <- "."   # e.g. "D:/UC Berkeley/KPMP paper/repo/Figure1"

data_dir   <- file.path(base_dir, "data")
output_dir <- file.path(base_dir, "Figure1_output")
figure_dir <- file.path(base_dir, "figuresFinal")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Input files ----------------------------------------------
kpmp_plasma_path      <- file.path(data_dir, "Plasma_Metabolomics.xlsx")
kpmp_clinical_path    <- file.path(data_dir, "Clinical_CKD_HealthyReference_only.xlsx")
validation_path       <- file.path(data_dir, "Validation_cleaned.xlsx")
validation_cluster_path <- file.path(data_dir, "Table_cluster_assignments.csv")

required_files <- c(
  kpmp_plasma_path,
  kpmp_clinical_path,
  validation_path,
  validation_cluster_path
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

# ---- Palettes -------------------------------------------------
group_colors <- c(
  "CKD" = "#fc8d62",
  "Healthy Reference" = "#66c2a5"
)

metabotype_colors <- c(
  "1"  = "#1b9e77",
  "2"  = "#7570b3",
  "3"  = "#d95f02",
  "4"  = "#e7298a",
  "5"  = "#66a61e",
  "6"  = "#e6ab02",
  "7"  = "#a6761d",
  "8"  = "#666666",
  "9"  = "#1f78b4",
  "10" = "#b2df8a"
)

origin_colors1 <- c(
  "microbiome-dependent" = "#1b9e77",
  "host-derived"         = "#7570b3",
  "diet-derived"         = "#d95f02",
  "xenobiotic"           = "#666666",
  "unclassified"         = "#999999"
)

cohort_colors <- c(
  "KPMP" = "#2C6E7F",
  "Validation" = "#8C6D5C"
)

# ---- Helper functions -----------------------------------------
clean_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "`", "")
  x <- stringr::str_replace_all(x, "\\s+", "")
  x <- stringr::str_replace(x, "^g(?=[0-9]{1,3}-[0-9]{4,6}$)", "")
  x
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

save_png <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(figure_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
}

find_first_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

extract_bracket_id <- function(x) {
  x <- as.character(x)
  stringr::str_match(x, "\\[(\\d+)\\]")[, 2]
}

extract_validation_name <- function(x) {
  x <- as.character(x)
  stringr::str_trim(stringr::str_remove(x, "\\s*\\[[0-9]+\\]$"))
}

impute_scale_numeric_matrix <- function(X_raw) {
  X_raw <- as.matrix(X_raw)

  for (j in seq_len(ncol(X_raw))) {
    X_raw[, j] <- suppressWarnings(as.numeric(X_raw[, j]))
  }

  # Log-transform only non-negative matrices.
  if (all(X_raw >= 0 | is.na(X_raw))) {
    X_raw <- log2(X_raw + 1)
  }

  means <- colMeans(X_raw, na.rm = TRUE)
  sds   <- apply(X_raw, 2, sd, na.rm = TRUE)

  means[!is.finite(means)] <- 0
  sds[!is.finite(sds) | sds == 0] <- 1

  for (j in seq_len(ncol(X_raw))) {
    nas <- is.na(X_raw[, j])
    if (any(nas)) X_raw[nas, j] <- means[j]
  }

  scale(X_raw, center = means, scale = sds)
}

impute_scale_matrix <- function(df, id_col, y_col) {
  X_raw <- df %>%
    select(-all_of(c(id_col, y_col))) %>%
    as.matrix()

  for (j in seq_len(ncol(X_raw))) {
    X_raw[, j] <- suppressWarnings(as.numeric(X_raw[, j]))
  }

  means <- colMeans(X_raw, na.rm = TRUE)
  sds   <- apply(X_raw, 2, sd, na.rm = TRUE)

  means[!is.finite(means)] <- 0
  sds[!is.finite(sds) | sds == 0] <- 1

  for (j in seq_len(ncol(X_raw))) {
    nas <- is.na(X_raw[, j])
    if (any(nas)) X_raw[nas, j] <- means[j]
  }

  X <- scale(X_raw, center = means, scale = sds)

  list(
    X = X,
    mean = means,
    sd = sds,
    features = colnames(X_raw)
  )
}

# ================================================================
# COHORT-SPECIFIC PREPROCESSING
# ================================================================

prep_kpmp_matrix <- function(path_xlsx) {
  message("Loading KPMP plasma metabolomics...")

  plasma_raw <- readxl::read_excel(path_xlsx) %>%
    janitor::clean_names()

  pid_col  <- find_first_col(plasma_raw, c("participant_id", "participant", "participantid"))
  met_col  <- find_first_col(plasma_raw, c("metabolite"))
  ab_col   <- find_first_col(plasma_raw, c("abundance", "value", "intensity", "peak_area"))
  name_col <- find_first_col(plasma_raw, c("chemical_name", "chemical_name_", "chem_name", "biochemical"))

  if (any(is.na(c(pid_col, met_col, ab_col)))) {
    stop("Could not detect required KPMP columns in Plasma_Metabolomics.xlsx")
  }

  plasma_long <- plasma_raw %>%
    transmute(
      ParticipantID    = clean_id(.data[[pid_col]]),
      metabolite       = as.character(.data[[met_col]]),
      metabolite_label = if (!is.na(name_col)) as.character(.data[[name_col]]) else as.character(.data[[met_col]]),
      abundance        = suppressWarnings(as.numeric(.data[[ab_col]]))
    ) %>%
    filter(!is.na(ParticipantID), !is.na(metabolite), !is.na(abundance))

  labels <- plasma_long %>%
    distinct(metabolite, metabolite_label) %>%
    mutate(
      metabolite_label = ifelse(
        is.na(metabolite_label) | metabolite_label == "",
        metabolite,
        metabolite_label
      )
    )

  wide <- plasma_long %>%
    group_by(ParticipantID, metabolite) %>%
    summarize(abundance = median(abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = metabolite, values_from = abundance)

  miss_tbl <- wide %>%
    summarize(across(-ParticipantID, ~ mean(is.na(.x)))) %>%
    pivot_longer(everything(), names_to = "metabolite", values_to = "miss_rate")

  keep_mets <- miss_tbl %>%
    filter(miss_rate <= 0.70) %>%
    pull(metabolite)

  wide <- wide %>%
    select(ParticipantID, all_of(keep_mets))

  list(wide = wide, labels = labels)
}

prep_kpmp_clinical <- function(path_xlsx) {
  message("Loading KPMP clinical metadata...")

  clinical_raw <- readxl::read_excel(path_xlsx) %>%
    janitor::clean_names()

  pid_col <- find_first_col(
    clinical_raw,
    c("participant_id", "participant", "participantid")
  )

  sample_col <- find_first_col(
    clinical_raw,
    c("sample_id", "sampleid", "specimen_id")
  )

  group_col <- find_first_col(
    clinical_raw,
    c("group", "enrollment_category", "disease_group")
  )

  egfr_col <- find_first_col(
    clinical_raw,
    c(
      "baseline_e_gfr_ml_min_1_73m2",
      "baseline_egfr_ml_min_1_73m2",
      "baseline_egfr_ml_min_1_73_m2",
      "egfr"
    )
  )

  if (any(is.na(c(pid_col, group_col, egfr_col)))) {
    stop("Could not detect participant ID, group, or eGFR columns in clinical file")
  }

  clinical_raw %>%
    transmute(
      ParticipantID = clean_id(.data[[pid_col]]),
      sample_id = if (!is.na(sample_col)) clean_id(.data[[sample_col]]) else ParticipantID,
      group = as.character(.data[[group_col]]),
      egfr = suppressWarnings(as.numeric(.data[[egfr_col]]))
    ) %>%
    distinct(ParticipantID, .keep_all = TRUE)
}

prep_validation_matrix <- function(path_xlsx) {
  message("Loading validation matrix...")

  combined_raw <- readxl::read_excel(path_xlsx, sheet = "Combined")

  req_cols <- c("sample_id", "cohort", "group")
  if (!all(req_cols %in% names(combined_raw))) {
    stop(
      "Validation_cleaned.xlsx Combined sheet must contain: ",
      paste(req_cols, collapse = ", ")
    )
  }

  meta <- combined_raw %>%
    transmute(
      sample_id = clean_id(sample_id),
      cohort = as.character(cohort),
      group = as.character(group)
    )

  feature_df <- combined_raw %>%
    select(-sample_id, -cohort, -group)

  original_feature_names <- colnames(feature_df)
  metabolite_ids   <- extract_bracket_id(original_feature_names)
  metabolite_names <- extract_validation_name(original_feature_names)

  keep_cols <- !is.na(metabolite_ids)
  feature_df <- feature_df[, keep_cols, drop = FALSE]
  metabolite_ids <- metabolite_ids[keep_cols]
  metabolite_names <- metabolite_names[keep_cols]

  feature_df <- feature_df %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.))))

  X0 <- as.matrix(feature_df)

  base_ids <- metabolite_ids
  uniq_ids <- unique(base_ids)

  X_list <- lapply(uniq_ids, function(id) {
    idx <- which(base_ids == id)
    vals <- X0[, idx, drop = FALSE]

    if (ncol(vals) == 1) {
      as.numeric(vals[, 1])
    } else {
      apply(vals, 1, function(v) median(v, na.rm = TRUE))
    }
  })

  X <- do.call(cbind, X_list)

  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
    colnames(X) <- uniq_ids[1]
  } else {
    colnames(X) <- uniq_ids
  }

  rownames(X) <- meta$sample_id

  sd_val   <- apply(X, 2, sd, na.rm = TRUE)
  prev_val <- 100 * colSums(!is.na(X) & X != 0) / nrow(X)
  keep     <- sd_val > 0 & !is.na(sd_val) & prev_val >= 30

  X <- X[, keep, drop = FALSE]
  X[is.na(X)] <- 0

  label_tbl <- tibble(
    metabolite = uniq_ids,
    validation_label = vapply(
      uniq_ids,
      function(id) metabolite_names[match(id, base_ids)][1],
      character(1)
    )
  ) %>%
    filter(metabolite %in% colnames(X)) %>%
    distinct(metabolite, .keep_all = TRUE)

  out <- tibble(sample_id = meta$sample_id) %>%
    bind_cols(as.data.frame(X))

  list(meta = meta, wide = out, labels = label_tbl)
}

prep_cluster_labels <- function(path_csv, id_candidates, cluster_candidates, out_id = "id") {
  cl <- readr::read_csv(path_csv, show_col_types = FALSE) %>%
    janitor::clean_names()

  id_col <- find_first_col(cl, id_candidates)
  mt_col <- find_first_col(cl, cluster_candidates)

  if (any(is.na(c(id_col, mt_col)))) {
    stop("Could not detect ID/metabotype columns in ", basename(path_csv))
  }

  out <- cl %>%
    transmute(
      !!out_id := clean_id(.data[[id_col]]),
      metabotype = suppressWarnings(as.integer(.data[[mt_col]]))
    ) %>%
    filter(!is.na(metabotype)) %>%
    distinct()

  ckd_id <- min(sort(unique(out$metabotype)))

  out %>%
    mutate(y = as.integer(metabotype == ckd_id))
}

extract_panel <- function(coef_obj, lambda_label) {
  df <- as.data.frame(as.matrix(coef_obj))
  df$metabolite <- rownames(df)
  names(df)[1] <- "coefficient"

  df %>%
    filter(metabolite != "(Intercept)", abs(coefficient) > 0) %>%
    arrange(desc(abs(coefficient))) %>%
    mutate(lambda = lambda_label) %>%
    select(metabolite, coefficient, lambda)
}

run_classifier <- function(X, y, nfolds_outer = 5, alpha = 0.5, seed = 42) {
  set.seed(seed)

  cv_full <- glmnet::cv.glmnet(
    x = X,
    y = y,
    family = "binomial",
    alpha = alpha,
    nfolds = 10,
    type.measure = "auc"
  )

  coef_min <- coef(cv_full, s = "lambda.min")
  coef_1se <- coef(cv_full, s = "lambda.1se")

  panel_min <- extract_panel(coef_min, "lambda.min")
  panel_1se <- extract_panel(coef_1se, "lambda.1se")

  fold_id <- sample(rep(seq_len(nfolds_outer), length.out = length(y)))
  outer_preds <- rep(NA_real_, length(y))
  outer_auc   <- numeric(nfolds_outer)

  for (k in seq_len(nfolds_outer)) {
    train_i <- which(fold_id != k)
    test_i  <- which(fold_id == k)

    cv_k <- glmnet::cv.glmnet(
      x = X[train_i, , drop = FALSE],
      y = y[train_i],
      family = "binomial",
      alpha = alpha,
      nfolds = 10,
      type.measure = "auc"
    )

    pred_k <- predict(
      cv_k,
      newx = X[test_i, , drop = FALSE],
      s = "lambda.1se",
      type = "response"
    )[, 1]

    outer_preds[test_i] <- pred_k
    outer_auc[k] <- as.numeric(
      pROC::auc(
        pROC::roc(y[test_i], pred_k, quiet = TRUE)
      )
    )
  }

  score_tbl <- tibble(
    y_true = y,
    prob_ckd_mt = outer_preds,
    MetabotypeScore = qlogis(pmin(pmax(outer_preds, 1e-8), 1 - 1e-8))
  )

  list(
    cv_full = cv_full,
    panel_min = panel_min,
    panel_1se = panel_1se,
    cv_predictions = score_tbl,
    auc_mean = mean(outer_auc),
    auc_sd = sd(outer_auc),
    roc = pROC::roc(y, outer_preds, quiet = TRUE)
  )
}

roc_to_df <- function(roc_obj, label) {
  tibble(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    dataset = label
  )
}

add_metabolite_labels <- function(panel_tbl, kpmp_labels, val_labels = NULL) {
  out <- panel_tbl %>%
    left_join(kpmp_labels, by = "metabolite")

  if (!is.null(val_labels)) {
    out <- out %>% left_join(val_labels, by = "metabolite")
  }

  out %>%
    mutate(
      plot_label = dplyr::coalesce(
        metabolite_label,
        validation_label,
        metabolite
      )
    )
}

load_metabolite_origin <- function(data_dir) {
  possible_paths <- c(
    file.path(data_dir, "metabolite_origin.csv"),
    file.path(data_dir, "metabolite_origins.csv"),
    file.path(data_dir, "metabolite_origin_lookup.csv")
  )

  origin_path <- possible_paths[file.exists(possible_paths)][1]

  if (is.na(origin_path)) {
    message(
      "No metabolite-origin lookup file found. ",
      "Using 'unclassified' for Figure 1G colors."
    )
    return(NULL)
  }

  origin_raw <- readr::read_csv(origin_path, show_col_types = FALSE) %>%
    janitor::clean_names()

  met_col <- find_first_col(
    origin_raw,
    c("metabolite", "metabolite_id", "chemical_id", "chem_id", "plot_label", "metabolite_name")
  )

  origin_col <- find_first_col(
    origin_raw,
    c("metabolite_origin", "origin", "source", "classification")
  )

  if (any(is.na(c(met_col, origin_col)))) {
    message("Metabolite-origin file found but required columns were not detected. Using unclassified.")
    return(NULL)
  }

  origin_raw %>%
    transmute(
      origin_key = as.character(.data[[met_col]]),
      metabolite_origin = as.character(.data[[origin_col]])
    ) %>%
    distinct(origin_key, .keep_all = TRUE)
}

# ================================================================
# LOAD KPMP INPUTS
# ================================================================

kpmp <- prep_kpmp_matrix(kpmp_plasma_path)
kpmp_clinical <- prep_kpmp_clinical(kpmp_clinical_path)

kpmp_meta <- kpmp$wide %>%
  select(ParticipantID) %>%
  inner_join(kpmp_clinical, by = "ParticipantID") %>%
  transmute(
    participant_id = ParticipantID,
    sample_id = sample_id,
    group = group,
    egfr = egfr
  )

kpmp_wide_for_clustering <- kpmp$wide %>%
  inner_join(
    kpmp_meta %>%
      transmute(ParticipantID = participant_id),
    by = "ParticipantID"
  )

metabolite_cols_kpmp <- setdiff(names(kpmp_wide_for_clustering), "ParticipantID")

X_scaled <- impute_scale_numeric_matrix(
  kpmp_wide_for_clustering %>%
    select(all_of(metabolite_cols_kpmp))
)

rownames(X_scaled) <- kpmp_wide_for_clustering$ParticipantID
d <- dist(X_scaled, method = "euclidean")

meta <- kpmp_meta %>%
  slice(match(kpmp_wide_for_clustering$ParticipantID, participant_id))

# ================================================================
# PANEL 1B: K-means Hartigan-Wong metabotype discovery and PCA
# ================================================================

k_grid <- 2:10

sil_scores <- sapply(k_grid, function(k) {
  tryCatch({
    set.seed(123)

    km_tmp <- kmeans(
      X_scaled,
      centers = k,
      nstart = 50,
      algorithm = "Hartigan-Wong"
    )

    mean(cluster::silhouette(km_tmp$cluster, d)[, "sil_width"])
  }, error = function(e) NA_real_)
})

valid_idx <- which(!is.na(sil_scores))

if (length(valid_idx) == 0) {
  stop("No valid silhouette scores were calculated.")
}

k_best <- k_grid[valid_idx[which.max(sil_scores[valid_idx])]]

set.seed(123)
km_final <- kmeans(
  X_scaled,
  centers = k_best,
  nstart = 50,
  algorithm = "Hartigan-Wong"
)

km_clusters <- km_final$cluster

# Relabel metabotypes by ascending mean eGFR:
# Metabotype 1 = lowest mean eGFR.
egfr_for_relabel <- meta$egfr
cluster_egfr <- tapply(egfr_for_relabel, km_clusters, mean, na.rm = TRUE)
rank_order <- order(cluster_egfr)
label_map <- setNames(seq_along(rank_order), rank_order)
km_clusters <- as.integer(label_map[as.character(km_clusters)])

cat(
  "K-means Hartigan-Wong: best k =",
  k_best,
  ", silhouette =",
  round(max(sil_scores[valid_idx]), 4),
  "\n"
)

cat(
  "Cluster sizes relabeled by eGFR:",
  paste(table(km_clusters), collapse = ", "),
  "\n"
)

sil_df <- data.frame(
  k = k_grid,
  silhouette = sil_scores
)

fig1B_silhouette <- ggplot(sil_df, aes(x = k, y = silhouette)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_point(
    data = sil_df %>% filter(k == k_best),
    color = "red",
    size = 5
  ) +
  labs(
    title = "Figure 1B: Silhouette scores",
    x = "k",
    y = "Mean silhouette width"
  ) +
  theme_cowplot() +
  panel_border(color = "black")

save_tiff(
  fig1B_silhouette,
  "Fig1B_silhouette_k_selection.tiff",
  width = 5,
  height = 4
)

results <- meta %>%
  mutate(metabotype = factor(km_clusters))

write.csv(
  results,
  file.path(output_dir, "Cluster_assignments_Q1.csv"),
  row.names = FALSE
)

patient_list_q1 <- results %>%
  select(participant_id, sample_id, group, metabotype)

write.csv(
  patient_list_q1,
  file.path(output_dir, "Patient_list_by_metabotype_Q1.csv"),
  row.names = FALSE
)

perm_cluster <- vegan::adonis2(d ~ km_clusters, permutations = 999)
perm_R2 <- round(perm_cluster$R2[1], 4)
perm_p  <- perm_cluster$`Pr(>F)`[1]

perm_subtitle <- paste0(
  "PERMANOVA: R\u00B2 = ",
  perm_R2,
  ", p = ",
  ifelse(perm_p < 0.001, "<0.001", round(perm_p, 3))
)

cat("=== PERMANOVA: K-means HW clusters ===\n")
cat(
  "R2 =",
  perm_R2,
  "  F =",
  round(perm_cluster$F[1], 2),
  "  p =",
  perm_p,
  "\n"
)

write.csv(
  data.frame(
    algorithm = "K-means",
    variant = "Hartigan-Wong",
    k = k_best,
    R2 = perm_R2,
    F_stat = round(perm_cluster$F[1], 2),
    p_value = perm_p
  ),
  file.path(output_dir, "Table_PERMANOVA_cluster_significance.csv"),
  row.names = FALSE
)

valid_group <- meta$group %in% c("CKD", "Healthy Reference")
has_two_groups <- length(unique(meta$group[valid_group])) >= 2

perm_group_sub <- ""

if (sum(valid_group) >= 10 && has_two_groups) {
  d_group <- dist(X_scaled[valid_group, ], method = "euclidean")

  perm_group <- vegan::adonis2(
    d_group ~ meta$group[valid_group],
    permutations = 999
  )

  perm_group_R2 <- round(perm_group$R2[1], 4)
  perm_group_p  <- perm_group$`Pr(>F)`[1]

  perm_group_sub <- paste0(
    "Group PERMANOVA: R\u00B2 = ",
    perm_group_R2,
    ", p = ",
    ifelse(perm_group_p < 0.001, "<0.001", round(perm_group_p, 3))
  )

  cat("=== PERMANOVA: CKD vs Healthy Reference ===\n")
  cat(
    "R2 =",
    perm_group_R2,
    "  F =",
    round(perm_group$F[1], 2),
    "  p =",
    perm_group_p,
    "\n"
  )

  write.csv(
    data.frame(
      biospecimen = "plasma",
      test = "CKD_vs_HealthyReference",
      R2 = perm_group_R2,
      F_stat = round(perm_group$F[1], 2),
      p_value = perm_group_p
    ),
    file.path(output_dir, "Table_PERMANOVA_group_significance.csv"),
    row.names = FALSE
  )
}

pca <- prcomp(
  X_scaled,
  center = TRUE,
  scale. = FALSE
)

var_explained <- round(
  100 * pca$sdev^2 / sum(pca$sdev^2),
  1
)

pca_base <- data.frame(
  participant_id = meta$participant_id,
  sample_id  = meta$sample_id,
  group      = meta$group,
  metabotype = factor(km_clusters),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

fig1B <- ggplot(
  pca_base,
  aes(PC1, PC2, color = metabotype, shape = group)
) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = metabotype_colors) +
  labs(
    title = "B",
    subtitle = perm_subtitle,
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)"),
    color = "Metabotype",
    shape = "Group"
  ) +
  theme_cowplot() +
  panel_border(color = "black") +
  theme(legend.position = "right")

save_tiff(
  fig1B,
  "Fig1B_KPMP_PCA_by_metabotype.tiff",
  width = 6.5,
  height = 5
)

message("Saved Fig 1B")

# ================================================================
# PREPARE CLUSTER LABELS FOR CLASSIFIER ANALYSES
# ================================================================

# Use the Figure 1B-derived KPMP cluster assignments for the KPMP classifier.
kpmp_cluster_path <- file.path(output_dir, "Cluster_assignments_Q1.csv")

kpmp_clusters <- prep_cluster_labels(
  kpmp_cluster_path,
  id_candidates = c("participant_id", "sample_id", "participantid"),
  cluster_candidates = c("metabotype", "cluster", "k_means_cluster"),
  out_id = "ParticipantID"
)

val <- prep_validation_matrix(validation_path)

val_clusters <- prep_cluster_labels(
  validation_cluster_path,
  id_candidates = c("sample_id", "participant_id", "participantid"),
  cluster_candidates = c("metabotype", "cluster", "k_means_cluster"),
  out_id = "sample_id"
)

# ================================================================
# MERGE LABELS AND RESTRICT TO SHARED METABOLITE IDs
# ================================================================

kpmp_df <- kpmp$wide %>%
  inner_join(
    kpmp_clusters %>%
      select(ParticipantID, metabotype, y),
    by = "ParticipantID"
  )

val_df <- val$wide %>%
  inner_join(
    val_clusters %>%
      select(sample_id, metabotype, y),
    by = "sample_id"
  )

shared_mets <- intersect(
  setdiff(colnames(kpmp_df), c("ParticipantID", "metabotype", "y")),
  setdiff(colnames(val_df), c("sample_id", "metabotype", "y"))
)

message("Shared metabolites by ID: ", length(shared_mets))

if (length(shared_mets) < 10) {
  stop(
    "Too few shared metabolites between KPMP and validation after ID harmonization: ",
    length(shared_mets)
  )
}

kpmp_df <- kpmp_df %>%
  select(ParticipantID, metabotype, y, all_of(shared_mets))

val_df <- val_df %>%
  select(sample_id, metabotype, y, all_of(shared_mets))

kpmp_scaled <- impute_scale_matrix(
  kpmp_df %>% select(-metabotype),
  "ParticipantID",
  "y"
)

val_scaled <- impute_scale_matrix(
  val_df %>% select(-metabotype),
  "sample_id",
  "y"
)

readr::write_csv(
  tibble(
    metabolite = kpmp_scaled$features,
    mean = kpmp_scaled$mean,
    sd = kpmp_scaled$sd
  ),
  file.path(output_dir, "KPMP_scaling_parameters.csv")
)

readr::write_csv(
  tibble(
    metabolite = val_scaled$features,
    mean = val_scaled$mean,
    sd = val_scaled$sd
  ),
  file.path(output_dir, "Validation_scaling_parameters.csv")
)

# ================================================================
# FIT SEPARATE KPMP AND VALIDATION CLASSIFIERS
# ================================================================

message("Running KPMP classifier...")
kpmp_fit <- run_classifier(kpmp_scaled$X, kpmp_df$y)

message("Running validation classifier...")
val_fit <- run_classifier(val_scaled$X, val_df$y)

saveRDS(
  kpmp_fit$cv_full,
  file.path(output_dir, "KPMP_glmnet_model.rds")
)

saveRDS(
  val_fit$cv_full,
  file.path(output_dir, "Validation_glmnet_model.rds")
)

readr::write_csv(
  kpmp_fit$panel_min,
  file.path(output_dir, "KPMP_panel_lambda_min.csv")
)

readr::write_csv(
  kpmp_fit$panel_1se,
  file.path(output_dir, "KPMP_panel_lambda_1se.csv")
)

readr::write_csv(
  val_fit$panel_min,
  file.path(output_dir, "Validation_panel_lambda_min.csv")
)

readr::write_csv(
  val_fit$panel_1se,
  file.path(output_dir, "Validation_panel_lambda_1se.csv")
)

kpmp_score_tbl <- bind_cols(
  tibble(
    ParticipantID = kpmp_df$ParticipantID,
    metabotype = kpmp_df$metabotype,
    cohort = "KPMP"
  ),
  kpmp_fit$cv_predictions
)

val_score_tbl <- bind_cols(
  tibble(
    sample_id = val_df$sample_id,
    metabotype = val_df$metabotype,
    cohort = "Validation"
  ),
  val_fit$cv_predictions
)

readr::write_csv(
  kpmp_score_tbl,
  file.path(output_dir, "KPMP_CV_predictions_and_score.csv")
)

readr::write_csv(
  val_score_tbl,
  file.path(output_dir, "Validation_CV_predictions_and_score.csv")
)

summary_tbl <- tibble(
  cohort = c("KPMP", "Validation"),
  n_samples = c(nrow(kpmp_df), nrow(val_df)),
  n_ckd_dominant = c(sum(kpmp_df$y), sum(val_df$y)),
  auc_mean = c(kpmp_fit$auc_mean, val_fit$auc_mean),
  auc_sd = c(kpmp_fit$auc_sd, val_fit$auc_sd),
  n_panel_lambda_min = c(nrow(kpmp_fit$panel_min), nrow(val_fit$panel_min)),
  n_panel_lambda_1se = c(nrow(kpmp_fit$panel_1se), nrow(val_fit$panel_1se))
)

readr::write_csv(
  summary_tbl,
  file.path(output_dir, "classifier_summary_table.csv")
)

# ================================================================
# PANEL 1C: Cross-validated ROC curves
# ================================================================

roc_df <- bind_rows(
  roc_to_df(kpmp_fit$roc, "KPMP"),
  roc_to_df(val_fit$roc, "Validation")
)

auc_labels <- summary_tbl %>%
  mutate(
    label = paste0(
      cohort,
      " AUC = ",
      round(auc_mean, 2),
      " \u00b1 ",
      round(auc_sd, 2)
    )
  )

fig1C <- ggplot(
  roc_df,
  aes(x = fpr, y = tpr, color = dataset)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey70"
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = cohort_colors) +
  annotate(
    "text",
    x = 0.65,
    y = 0.25,
    label = paste(auc_labels$label, collapse = "\n"),
    hjust = 0,
    size = 3.5
  ) +
  coord_equal() +
  labs(
    title = "C",
    x = "False positive rate",
    y = "True positive rate",
    color = "Cohort"
  ) +
  theme_cowplot() +
  panel_border(color = "black") +
  theme(legend.position = "bottom")

save_tiff(
  fig1C,
  "Fig1C_classifier_ROC_curves.tiff",
  width = 5,
  height = 5
)

message("Saved Fig 1C")

# ================================================================
# PANEL 1D: Cross-validated metabotype score by cohort/metabotype
# ================================================================

score_plot_df <- bind_rows(
  kpmp_score_tbl %>%
    transmute(
      sample_id = ParticipantID,
      cohort = cohort,
      metabotype = factor(metabotype),
      MetabotypeScore = MetabotypeScore
    ),
  val_score_tbl %>%
    transmute(
      sample_id = sample_id,
      cohort = cohort,
      metabotype = factor(metabotype),
      MetabotypeScore = MetabotypeScore
    )
)

fig1D <- ggplot(
  score_plot_df,
  aes(x = metabotype, y = MetabotypeScore, fill = metabotype)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.55) +
  facet_wrap(~ cohort, scales = "free_x") +
  scale_fill_manual(values = metabotype_colors) +
  labs(
    title = "D",
    x = "Metabotype",
    y = "Cross-validated metabotype score"
  ) +
  theme_cowplot() +
  panel_border(color = "black") +
  theme(legend.position = "none")

save_tiff(
  fig1D,
  "Fig1D_cross_validated_metabotype_score.tiff",
  width = 6,
  height = 4
)

message("Saved Fig 1D")

# ================================================================
# PANEL 1E: Top KPMP lambda.1se classifier coefficients
# ================================================================

kpmp_panel_plot <- kpmp_fit$panel_1se %>%
  add_metabolite_labels(kpmp$labels) %>%
  slice_max(abs(coefficient), n = 20, with_ties = FALSE) %>%
  mutate(
    direction = ifelse(coefficient > 0, "CKD-dominant metabotype", "Other metabotype"),
    plot_label = forcats::fct_reorder(plot_label, coefficient)
  )

fig1E <- ggplot(
  kpmp_panel_plot,
  aes(x = coefficient, y = plot_label, fill = direction)
) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  labs(
    title = "E",
    x = "KPMP classifier coefficient",
    y = NULL,
    fill = "Direction"
  ) +
  theme_cowplot() +
  panel_border(color = "black") +
  theme(legend.position = "bottom")

save_tiff(
  fig1E,
  "Fig1E_KPMP_top_classifier_coefficients.tiff",
  width = 6,
  height = 6
)

message("Saved Fig 1E")

# ================================================================
# PANEL 1F: Top validation lambda.1se classifier coefficients
# ================================================================

val_panel_plot <- val_fit$panel_1se %>%
  add_metabolite_labels(kpmp$labels, val$labels) %>%
  slice_max(abs(coefficient), n = 20, with_ties = FALSE) %>%
  mutate(
    direction = ifelse(coefficient > 0, "CKD-dominant metabotype", "Other metabotype"),
    plot_label = forcats::fct_reorder(plot_label, coefficient)
  )

fig1F <- ggplot(
  val_panel_plot,
  aes(x = coefficient, y = plot_label, fill = direction)
) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  labs(
    title = "F",
    x = "Validation classifier coefficient",
    y = NULL,
    fill = "Direction"
  ) +
  theme_cowplot() +
  panel_border(color = "black") +
  theme(legend.position = "bottom")

save_tiff(
  fig1F,
  "Fig1F_validation_top_classifier_coefficients.tiff",
  width = 6,
  height = 6
)

message("Saved Fig 1F")

# ================================================================
# PANEL 1G: Cross-cohort coefficient concordance
# ================================================================

align_tbl <- full_join(
  kpmp_fit$panel_1se %>%
    select(metabolite, coefficient_kpmp = coefficient),
  val_fit$panel_1se %>%
    select(metabolite, coefficient_validation = coefficient),
  by = "metabolite"
) %>%
  left_join(kpmp$labels, by = "metabolite") %>%
  left_join(val$labels, by = "metabolite") %>%
  mutate(
    plot_label = dplyr::coalesce(
      metabolite_label,
      validation_label,
      metabolite
    ),
    in_kpmp = !is.na(coefficient_kpmp),
    in_validation = !is.na(coefficient_validation),
    overlap_class = case_when(
      in_kpmp & in_validation ~ "Both",
      in_kpmp ~ "KPMP only",
      in_validation ~ "Validation only",
      TRUE ~ "Neither"
    )
  )

readr::write_csv(
  align_tbl,
  file.path(output_dir, "KPMP_validation_metabolite_alignment.csv")
)

shared_align1 <- align_tbl %>%
  filter(in_kpmp, in_validation)

readr::write_csv(
  shared_align1,
  file.path(output_dir, "KPMP_validation_shared_predictive_metabolites.csv")
)

origin_lookup <- load_metabolite_origin(data_dir)

if (!is.null(origin_lookup)) {
  shared_align1 <- shared_align1 %>%
    left_join(
      origin_lookup,
      by = c("metabolite" = "origin_key")
    )

  # If the origin lookup keyed by numeric metabolite ID did not match,
  # try matching by plot label.
  if (all(is.na(shared_align1$metabolite_origin))) {
    shared_align1 <- shared_align1 %>%
      select(-metabolite_origin) %>%
      left_join(
        origin_lookup,
        by = c("plot_label" = "origin_key")
      )
  }
} else {
  shared_align1 <- shared_align1 %>%
    mutate(metabolite_origin = "unclassified")
}

shared_align1 <- shared_align1 %>%
  mutate(
    metabolite_origin = ifelse(
      is.na(metabolite_origin) | metabolite_origin == "",
      "unclassified",
      metabolite_origin
    )
  )

top_metabolites <- shared_align1 %>%
  slice_max(
    abs(coefficient_kpmp) + abs(coefficient_validation),
    n = 25,
    with_ties = FALSE
  )

cor_res <- top_metabolites %>%
  summarise(
    rho = cor(
      coefficient_kpmp,
      coefficient_validation,
      method = "spearman",
      use = "complete.obs"
    ),
    p_value = suppressWarnings(
      cor.test(
        coefficient_kpmp,
        coefficient_validation,
        method = "spearman",
        exact = FALSE
      )$p.value
    )
  )

rho_lab <- paste0(
  "Spearman \u03c1 = ",
  round(cor_res$rho, 2),
  "\np = ",
  signif(cor_res$p_value, 2)
)

fig1G <- top_metabolites %>%
  ggplot(
    aes(
      x = coefficient_kpmp,
      y = coefficient_validation,
      label = plot_label,
      color = metabolite_origin
    )
  ) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
  geom_point(alpha = 0.9, size = 3) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "steelblue",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  ggrepel::geom_text_repel(
    size = 3,
    max.overlaps = 30
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = rho_lab,
    hjust = 1.1,
    vjust = -0.2,
    size = 4
  ) +
  scale_color_manual(values = origin_colors1) +
  labs(
    title = "G",
    x = "KPMP coefficient\nDiscovery cohort",
    y = "Validation coefficient\nValidation cohort",
    color = "Metabolite origin"
  ) +
  theme_cowplot() +
  panel_border(color = "black") +
  theme(legend.position = "bottom")

save_tiff(
  fig1G,
  "Fig1G_KPMP_validation_coefficient_concordance.tiff",
  width = 6,
  height = 5.5
)

message("Saved Fig 1G")

# ================================================================
# OPTIONAL COMBINED FIGURE 1 PANELS B-G
# ================================================================

combined_fig1_BG <- (
  fig1B | fig1C
) / (
  fig1D | fig1G
)

save_tiff(
  combined_fig1_BG,
  "Fig1_combined_panels_B_C_D_G.tiff",
  width = 12,
  height = 10
)

message("\n========== Figure 1 complete ==========")
message("Outputs in: ", output_dir)
message("Figures in: ", figure_dir)
