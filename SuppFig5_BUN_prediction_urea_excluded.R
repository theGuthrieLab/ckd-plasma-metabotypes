#!/usr/bin/env Rscript

################################################################
# Supplemental Figure 5 (Reviewer 3, comment 4)
# Figure 3 sensitivity analysis: is the improvement in BUN
# prediction driven by plasma urea?
#
# Data assembly, cross-validation, and the R2 definition are
# identical to Figure3_plasma_prediction_serum_filtration_markers.R (mean fold-level squared correlation,
# averaged over N_REPEATS random 5-fold partitions). With the same
# settings, row (i) reproduces the Figure 3A BUN values exactly.
#
# It also reports plasma urea's rank by permutation importance and
# plots the top-25 BUN predictors with urea removed (Supplemental Figure 5).
#
# INPUTS (data/): same as the Figure 3 script
#   Plasma_Metabolomics.xlsx
#   Serum-biomarkers-KPMP-CKD_Healthy_with_Baseline_eGFR.csv
#   metabolite_origin_main_final.csv
################################################################

options(stringsAsFactors = FALSE)

# ---- Packages -------------------------------------------------
req_pkgs <- c("tidyverse", "readxl", "janitor", "ranger")

missing_pkgs <- req_pkgs[
  !vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall these before running.")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(ranger)
})

# ================================================================
# EDIT THIS ONE LINE to point to the folder containing the inputs
# ================================================================
base_dir <- "."   # folder that contains data/, e.g. "D:/UC Berkeley/KPMP paper/repo"

data_dir   <- file.path(base_dir, "data")
output_dir <- file.path(base_dir, "SuppFig5_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Analysis settings (keep identical to the Figure 3 script) -
N_FOLDS   <- 5
N_REPEATS <- 10
NUM_TREES <- 1000
SEED      <- 1

# ---- Input files ----------------------------------------------
plasma_path <- file.path(data_dir, "Plasma_Metabolomics.xlsx")
serum_path  <- file.path(data_dir, "Serum-biomarkers-KPMP-CKD_Healthy_with_Baseline_eGFR.csv")
origin_path <- file.path(data_dir, "metabolite_origin_main_final.csv")

for (p in c(plasma_path, serum_path)) {
  if (!file.exists(p)) stop("Required input not found: ", p)
}

origin_colors <- c(
  "diet-derived"         = "#E66101",  # orange
  "host-derived"         = "#9590C8",  # purple
  "microbiome-dependent" = "#1B9E77",  # teal
  "xenobiotic"           = "#999999",  # gray
  "unclassified"         = "#D9D9D9"   # light gray (no origin assigned)
)

# ================================================================
# 1. LOAD AND ASSEMBLE DATA
# (identical in this script and the Figure 3 script)
# ================================================================
covars   <- c("age_years_binned", "sex", "baseline_e_gfr_ml_min_1_73m2")
outcomes <- c(BUN          = "blood_urea_nitrogen_mg_d_l",
              Creatinine   = "creatinine_mg_d_l",
              `Cystatin C` = "cystatin_c_mg_l")

serum <- read_csv(serum_path, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(participant_id = as.character(participant_id))

missing_cols <- setdiff(c("participant_id", covars, outcomes), names(serum))
if (length(missing_cols) > 0) {
  stop("Columns not found in the serum file after clean_names(): ",
       paste(missing_cols, collapse = ", "))
}

serum <- serum %>%
  mutate(
    sex              = factor(sex),
    age_years_binned = factor(age_years_binned),
    baseline_e_gfr_ml_min_1_73m2 =
      suppressWarnings(as.numeric(baseline_e_gfr_ml_min_1_73m2)),
    across(all_of(unname(outcomes)), ~ suppressWarnings(as.numeric(.x)))
  )

plasma_raw <- read_excel(plasma_path) %>% clean_names()

need_plasma <- c("participant_id", "metabolite", "abundance", "chemical_name")
missing_plasma <- setdiff(need_plasma, names(plasma_raw))
if (length(missing_plasma) > 0) {
  stop("Columns not found in the plasma file after clean_names(): ",
       paste(missing_plasma, collapse = ", "))
}

plasma_long <- plasma_raw %>%
  transmute(
    participant_id = as.character(participant_id),
    metabolite     = as.character(metabolite),
    abundance      = suppressWarnings(as.numeric(abundance))
  ) %>%
  filter(!is.na(participant_id), !is.na(metabolite))

met_lookup <- plasma_raw %>%
  transmute(metabolite    = as.character(metabolite),
            chemical_name = as.character(chemical_name)) %>%
  distinct(metabolite, .keep_all = TRUE)

plasma_wide <- plasma_long %>%
  group_by(participant_id, metabolite) %>%
  summarise(abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = metabolite, values_from = abundance)

met_ids  <- setdiff(names(plasma_wide), "participant_id")
features <- make.names(met_ids, unique = TRUE)
names(plasma_wide) <- c("participant_id", features)

feature_map <- tibble(feature = features, metabolite = met_ids) %>%
  left_join(met_lookup, by = "metabolite") %>%
  mutate(label = coalesce(chemical_name, metabolite))

dat <- plasma_wide %>%
  inner_join(serum, by = "participant_id") %>%
  as.data.frame()

metabolite_cols <- features

for (m in metabolite_cols) {
  v <- dat[[m]]
  med <- median(v, na.rm = TRUE)
  if (!is.finite(med)) med <- 0
  v[is.na(v)] <- med
  dat[[m]] <- v
}

message("Merged data: ", nrow(dat), " participants, ",
        length(metabolite_cols), " metabolites")

# ================================================================
# 2. REPEATED CROSS-VALIDATION
# (identical in this script and the Figure 3 script)
# ================================================================
repeated_cv_r2 <- function(dat, outcome, feats,
                           n_folds = N_FOLDS, n_repeats = N_REPEATS,
                           num_trees = NUM_TREES, seed = SEED) {
  required <- unique(c(outcome, intersect(feats, covars)))
  d <- dat[complete.cases(dat[, required, drop = FALSE]),
           c(outcome, feats), drop = FALSE]
  y <- d[[outcome]]
  x <- d[, feats, drop = FALSE]
  n <- nrow(d)
  if (n < 20) stop("Too few complete rows for cross-validation: ", n)

  per_repeat <- vapply(seq_len(n_repeats), function(r) {
    set.seed(seed + r)
    fold <- sample(rep_len(seq_len(n_folds), n))
    fold_r2 <- vapply(seq_len(n_folds), function(k) {
      tr <- fold != k
      te <- fold == k
      fit <- ranger(x = x[tr, , drop = FALSE], y = y[tr],
                    num.trees = num_trees, seed = seed + r)
      pred <- predict(fit, data = x[te, , drop = FALSE])$predictions
      suppressWarnings(cor(y[te], pred))^2
    }, numeric(1))
    mean(fold_r2, na.rm = TRUE)
  }, numeric(1))

  list(per_repeat = per_repeat, n = n)
}

# ================================================================
# 3. EXCLUSION SETS (exact chemical-name matches)
# Exact matching avoids catching related analytes such as
# 3-ureidopropionate or phosphocreatine.
# ================================================================
find_feature <- function(name) feature_map$feature[feature_map$label == name]

urea_f       <- find_feature("urea")
creatinine_f <- find_feature("creatinine")
creatine_f   <- find_feature("creatine")

if (length(urea_f) == 0) stop("No metabolite named exactly 'urea' was found.")

excl_primary      <- urea_f
excl_conservative <- c(urea_f, creatinine_f, creatine_f)

message("\nExcluded (primary):      ",
        paste(feature_map$label[feature_map$feature %in% excl_primary], collapse = ", "))
message("Excluded (conservative): ",
        paste(feature_map$label[feature_map$feature %in% excl_conservative], collapse = ", "))

met_sets <- list(
  "All metabolites"                         = metabolite_cols,
  "Plasma urea excluded"                    = setdiff(metabolite_cols, excl_primary),
  "Urea, creatinine, and creatine excluded" = setdiff(metabolite_cols, excl_conservative)
)

# ================================================================
# 4. CROSS-VALIDATED R2 FOR BUN
# ================================================================
outcome_bun <- outcomes[["BUN"]]

message("\nRepeated cross-validation for BUN (", N_REPEATS, " x ", N_FOLDS, "-fold) ...")

base_res <- repeated_cv_r2(dat, outcome_bun, covars)
message(sprintf("  Baseline clinical: CV R2 = %.3f (SD %.3f)",
                mean(base_res$per_repeat), sd(base_res$per_repeat)))

by_repeat <- purrr::imap_dfr(met_sets, function(mets, set_name) {
  both <- repeated_cv_r2(dat, outcome_bun, c(covars, mets))
  only <- repeated_cv_r2(dat, outcome_bun, mets)
  message(sprintf("  %-40s baseline+met = %.3f  met-only = %.3f",
                  set_name, mean(both$per_repeat), mean(only$per_repeat)))
  tibble(
    feature_set             = set_name,
    n_metabolites           = length(mets),
    repeat_id               = seq_len(N_REPEATS),
    cv_r2_baseline          = base_res$per_repeat,
    cv_r2_baseline_plus_met = both$per_repeat,
    cv_r2_metabolites_only  = only$per_repeat,
    delta_r2                = both$per_repeat - base_res$per_repeat
  )
})

sens_summary <- by_repeat %>%
  group_by(feature_set, n_metabolites) %>%
  summarise(
    baseline_mean          = mean(cv_r2_baseline),
    baseline_sd            = sd(cv_r2_baseline),
    baseline_plus_met_mean = mean(cv_r2_baseline_plus_met),
    baseline_plus_met_sd   = sd(cv_r2_baseline_plus_met),
    metabolites_only_mean  = mean(cv_r2_metabolites_only),
    metabolites_only_sd    = sd(cv_r2_metabolites_only),
    delta_r2_mean          = mean(delta_r2),
    delta_r2_sd            = sd(delta_r2),
    n_repeats_positive     = sum(delta_r2 > 0),
    .groups = "drop"
  ) %>%
  mutate(feature_set = factor(feature_set, levels = names(met_sets))) %>%
  arrange(feature_set)

write.csv(by_repeat, file.path(output_dir, "R3C4_BUN_urea_sensitivity_by_repeat.csv"),
          row.names = FALSE)
write.csv(sens_summary, file.path(output_dir, "R3C4_BUN_urea_sensitivity_summary.csv"),
          row.names = FALSE)

message("\nSummary (mean and SD across repeats):")
print(as.data.frame(sens_summary %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))

# ================================================================
# 5. PERMUTATION IMPORTANCE, WITH AND WITHOUT UREA
# ================================================================
fit_importance <- function(mets) {
  d <- dat[!is.na(dat[[outcome_bun]]), c(outcome_bun, mets)]
  fit <- ranger(x = d[, mets, drop = FALSE], y = d[[outcome_bun]],
                num.trees = NUM_TREES, importance = "permutation",
                sample.fraction = 0.632, seed = SEED)
  tbl <- tibble(feature = names(fit$variable.importance),
                importance = as.numeric(fit$variable.importance)) %>%
    left_join(feature_map, by = "feature") %>%
    arrange(desc(importance)) %>%
    mutate(rank = row_number())
  list(tbl = tbl, oob_r2 = cor(d[[outcome_bun]], fit$predictions)^2)
}

imp_all <- fit_importance(met_sets[["All metabolites"]])
imp_noU <- fit_importance(met_sets[["Plasma urea excluded"]])

urea_rank <- imp_all$tbl$rank[imp_all$tbl$feature %in% urea_f]
message(sprintf("\nOut-of-bag R2, all metabolites:      %.3f", imp_all$oob_r2))
message(sprintf("Out-of-bag R2, plasma urea excluded: %.3f", imp_noU$oob_r2))
message(sprintf("Plasma urea importance rank: %d of %d",
                urea_rank[1], nrow(imp_all$tbl)))

write.csv(
  bind_rows(
    imp_all$tbl %>% mutate(feature_set = "All metabolites"),
    imp_noU$tbl %>% mutate(feature_set = "Plasma urea excluded")
  ) %>% select(feature_set, rank, label, metabolite, importance),
  file.path(output_dir, "R3C4_BUN_importance_with_without_urea.csv"),
  row.names = FALSE
)

write.csv(
  tibble(oob_r2_all_metabolites = imp_all$oob_r2,
         oob_r2_urea_excluded   = imp_noU$oob_r2,
         urea_rank = urea_rank[1],
         n_metabolites = nrow(imp_all$tbl)),
  file.path(output_dir, "R3C4_BUN_oob_and_urea_rank.csv"),
  row.names = FALSE
)

# ================================================================
# 6. TOP-METABOLITE PANEL, UREA EXCLUDED
# ================================================================
origin_lookup <- NULL
if (file.exists(origin_path)) {
  origin_lookup <- read_csv(origin_path, show_col_types = FALSE) %>%
    transmute(label = metabolite, origin = origin) %>%
    distinct(label, .keep_all = TRUE)
}


# ---- Vertical importance plot (matches the submitted Figure 3C style) ----
TOP_N <- 25

# Shorter display names for long metabolite labels
short_labels <- c(
  "N,N,N-trimethyl-alanylproline betaine (TMAP)"     = "TMAP",
  "2,3-dihydroxy-5-methylthio-4-pentenoate (DMTPA)*" = "DMTPA*"
)

plot_importance_vertical <- function(df, title = NULL) {
  df <- df %>%
    mutate(display = ifelse(label %in% names(short_labels),
                            unname(short_labels[label]), label))
  df$display <- factor(df$display, levels = df$display[order(df$importance)])

  # Show all four origins in the legend; add "unclassified" only if present
  legend_levels <- names(origin_colors)
  if (!any(df$origin == "unclassified")) {
    legend_levels <- setdiff(legend_levels, "unclassified")
  }

  legend_theme <- if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside",
          legend.position.inside = c(0.02, 0.98),
          legend.justification = c(0, 1))
  } else {
    theme(legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1))
  }

  ggplot(df, aes(display, importance, fill = origin)) +
    geom_col(width = 0.9, key_glyph = "point", show.legend = TRUE) +
    scale_fill_manual(values = origin_colors, limits = legend_levels,
                      name = "Metabolite origin", drop = FALSE) +
    guides(fill = guide_legend(
      override.aes = list(shape = 21, size = 4, colour = NA))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, x = NULL, y = "Permutation importance") +
    theme_classic(base_size = 13) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
          plot.title = element_text(face = "bold", size = 14),
          legend.background = element_blank(),
          legend.key = element_blank()) +
    legend_theme
}

top_imp <- imp_noU$tbl %>%
  filter(!grepl("^X\\s*-\\s*[0-9]", label)) %>%   # drop unnamed Metabolon compounds
  slice_head(n = TOP_N)

if (!is.null(origin_lookup)) {
  top_imp <- top_imp %>% left_join(origin_lookup, by = "label")
} else {
  top_imp$origin <- NA_character_
}
top_imp <- top_imp %>%
  mutate(origin = ifelse(is.na(origin), "unclassified", origin))

n_uncl <- sum(top_imp$origin == "unclassified")
if (n_uncl > 0) {
  message("NOTE: ", n_uncl, " of the top-", TOP_N, " metabolites have no origin assignment:")
  message("  ", paste(top_imp$label[top_imp$origin == "unclassified"], collapse = "\n  "))
}

write.csv(top_imp %>% select(rank, label, importance, origin),
          file.path(output_dir, "R3C4_top_BUN_urea_excluded_with_origin.csv"),
          row.names = FALSE)

p_top <- plot_importance_vertical(
  top_imp, title = "Top plasma metabolites predicting BUN (plasma urea excluded)")

ggsave(file.path(output_dir, "R3C4_top_BUN_urea_excluded.png"),
       p_top, width = 10, height = 7, dpi = 300)
ggsave(file.path(output_dir, "R3C4_top_BUN_urea_excluded.pdf"),
       p_top, width = 10, height = 7)
ggsave(file.path(output_dir, "Supplemental_Figure_5.tiff"),
       p_top, width = 10, height = 7, dpi = 600, compression = "lzw", bg = "white")

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo_Figure3_urea_sensitivity.txt"))

message("\n========== Sensitivity analysis complete ==========")
message("Outputs in: ", output_dir)
