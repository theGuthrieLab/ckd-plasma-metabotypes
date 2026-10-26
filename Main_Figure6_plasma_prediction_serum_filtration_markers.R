#!/usr/bin/env Rscript

################################################################
# Figure 6: Plasma metabolite prediction of serum filtration markers
#
# 6A: Cross-validated R2 for BUN, creatinine, and cystatin C under
#     three random forest configurations:
#       - Baseline clinical (age, sex, baseline eGFR)
#       - Baseline + metabolites
#       - Metabolites only
#     Bars show the mean across repeated 5-fold cross-validation;
#     error bars show the SD across repeats.
# 6B: Random forest out-of-bag observed vs predicted BUN
#     (metabolites only); R2 is the squared correlation of the
#     out-of-bag predictions with observed BUN.
# 6C: Top 25 permutation-importance metabolites for BUN, colored by
#     metabolite origin.
#
#
# INPUTS (data/)
#   Plasma_Metabolomics.xlsx
#       long format: Participant ID, Metabolite, Abundance,
#       CHEMICAL_NAME, ...
#   Serum-biomarkers-KPMP-CKD_Healthy_with_Baseline_eGFR.csv
#       participant-level serum biomarkers and covariates
#   metabolite_origin_main_final.csv
#       columns: metabolite (chemical name), origin, source
#
################################################################

options(stringsAsFactors = FALSE)

# ---- Packages -------------------------------------------------
req_pkgs <- c("tidyverse", "readxl", "janitor", "ranger", "patchwork")

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
  library(patchwork)
})

# ================================================================
# EDIT THIS ONE LINE to point to the folder that contains data/
# ================================================================
base_dir <- "."   # e.g. "D:/UC Berkeley/KPMP paper/repo"

data_dir   <- file.path(base_dir, "data")
output_dir <- file.path(base_dir, "Figure6_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Analysis settings ----------------------------------------
N_FOLDS   <- 5
N_REPEATS <- 10
NUM_TREES <- 1000
SEED      <- 1
TOP_N     <- 25

# ---- Figure settings ------------------------------------------
TIFF_DPI    <- 600
BASE_SIZE   <- 18   # panels A and B
C_BASE_SIZE <- 18   # panel C
C_LABEL     <- 13   # panel C metabolite names
TAG_SIZE    <- 26   # panel letters

# ---- Input files ----------------------------------------------
plasma_path <- file.path(data_dir, "Plasma_Metabolomics.xlsx")
serum_path  <- file.path(data_dir, "Serum-biomarkers-KPMP-CKD_Healthy_with_Baseline_eGFR.csv")
origin_path <- file.path(data_dir, "metabolite_origin_main_final.csv")

for (p in c(plasma_path, serum_path, origin_path)) {
  if (!file.exists(p)) stop("Required input not found: ", p)
}

# ---- Palettes -------------------------------------------------
bar_colors <- c(
  "Baseline clinical"      = "#BDBDBD",  # light gray
  "Baseline + metabolites" = "#2C6E7F",  # teal
  "Metabolites only"       = "#8C6D5C"   # brown
)

origin_colors <- c(
  "diet-derived"         = "#E66101",  # orange
  "host-derived"         = "#9590C8",  # purple
  "microbiome-dependent" = "#1B9E77",  # teal
  "xenobiotic"           = "#999999",  # gray
  "unclassified"         = "#D9D9D9"   # light gray (no origin assigned)
)

box_theme <- theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))

# ================================================================
# 1. LOAD AND ASSEMBLE DATA
# ================================================================
covars   <- c("age_years_binned", "sex", "baseline_e_gfr_ml_min_1_73m2")
outcomes <- c(BUN          = "blood_urea_nitrogen_mg_d_l",
              Creatinine   = "creatinine_mg_d_l",
              `Cystatin C` = "cystatin_c_mg_l")

# ---- Serum biomarkers and covariates --------------------------
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

# ---- Plasma metabolomics (long -> wide) -----------------------
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

# Metabolite ID -> chemical name
met_lookup <- plasma_raw %>%
  transmute(metabolite    = as.character(metabolite),
            chemical_name = as.character(chemical_name)) %>%
  distinct(metabolite, .keep_all = TRUE)

plasma_wide <- plasma_long %>%
  group_by(participant_id, metabolite) %>%
  summarise(abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = metabolite, values_from = abundance)

# Syntactically valid feature names (numeric IDs become e.g. X100000007)
met_ids  <- setdiff(names(plasma_wide), "participant_id")
features <- make.names(met_ids, unique = TRUE)
names(plasma_wide) <- c("participant_id", features)

feature_map <- tibble(feature = features, metabolite = met_ids) %>%
  left_join(met_lookup, by = "metabolite") %>%
  mutate(label = coalesce(chemical_name, metabolite))

# ---- Merge ----------------------------------------------------
dat <- plasma_wide %>%
  inner_join(serum, by = "participant_id") %>%
  as.data.frame()

metabolite_cols <- features

# Median-impute metabolite values
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
# 3. CROSS-VALIDATED R2 (PANEL 6A VALUES)
# ================================================================
model_sets <- list(
  "Baseline clinical"      = covars,
  "Baseline + metabolites" = c(covars, metabolite_cols),
  "Metabolites only"       = metabolite_cols
)

message("\nRepeated cross-validation: ", N_REPEATS, " x ", N_FOLDS,
        "-fold, ", length(outcomes) * length(model_sets), " models ...")

cv_long <- purrr::imap_dfr(outcomes, function(ocol, oname) {
  purrr::imap_dfr(model_sets, function(feats, mname) {
    res <- repeated_cv_r2(dat, ocol, feats)
    message(sprintf("  %-11s | %-23s n = %3d  CV R2 = %.3f (SD %.3f)",
                    oname, mname, res$n,
                    mean(res$per_repeat), sd(res$per_repeat)))
    tibble(outcome = oname, model = mname, n = res$n,
           repeat_id = seq_along(res$per_repeat),
           cv_r2 = res$per_repeat)
  })
})

cv_summary <- cv_long %>%
  group_by(outcome, model, n) %>%
  summarise(mean_cv_r2 = mean(cv_r2), sd_cv_r2 = sd(cv_r2), .groups = "drop")

# Paired delta R2 (baseline + metabolites minus baseline, same folds)
delta_summary <- cv_long %>%
  filter(model %in% c("Baseline clinical", "Baseline + metabolites")) %>%
  select(outcome, model, repeat_id, cv_r2) %>%
  pivot_wider(names_from = model, values_from = cv_r2) %>%
  mutate(delta_r2 = `Baseline + metabolites` - `Baseline clinical`) %>%
  group_by(outcome) %>%
  summarise(mean_delta_r2 = mean(delta_r2), sd_delta_r2 = sd(delta_r2),
            n_repeats_positive = sum(delta_r2 > 0), .groups = "drop")

write.csv(cv_long, file.path(output_dir, "Fig6A_CV_R2_by_repeat.csv"),
          row.names = FALSE)
write.csv(cv_summary, file.path(output_dir, "Fig6A_CV_R2_summary.csv"),
          row.names = FALSE)
write.csv(delta_summary, file.path(output_dir, "Fig6A_delta_R2_summary.csv"),
          row.names = FALSE)

message("\nCross-validated R2 (mean, SD across repeats):")
print(as.data.frame(cv_summary %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
message("\nDelta R2, baseline + metabolites vs baseline:")
print(as.data.frame(delta_summary %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))


# ================================================================
# 4. PANEL 6A
# ================================================================
plot_a <- cv_summary %>%
  mutate(model   = factor(model, levels = names(bar_colors)),
         outcome = factor(outcome, levels = names(outcomes)))

p6a <- ggplot(plot_a, aes(outcome, mean_cv_r2, fill = model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(aes(ymin = pmax(mean_cv_r2 - sd_cv_r2, 0),
                    ymax = mean_cv_r2 + sd_cv_r2),
                position = position_dodge(width = 0.75), width = 0.2, linewidth = 0.6) +
  scale_fill_manual(values = bar_colors) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = expression(paste("Cross-validated ", R^2)), fill = "Model") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_classic(base_size = BASE_SIZE) +
  box_theme +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(size = BASE_SIZE))

# ================================================================
# 5. PANELS 6B & 6C: BUN random forest (metabolites only)
# ================================================================
outcome_bun <- outcomes[["BUN"]]

d_bun <- dat[!is.na(dat[[outcome_bun]]), c("participant_id", outcome_bun, metabolite_cols)]

rf_bun <- ranger(
  x = d_bun[, metabolite_cols, drop = FALSE],
  y = d_bun[[outcome_bun]],
  num.trees = NUM_TREES,
  importance = "permutation",
  sample.fraction = 0.632,
  seed = SEED
)

obs  <- d_bun[[outcome_bun]]
pred <- rf_bun$predictions
oob_r2 <- cor(obs, pred)^2
message(sprintf("\nBUN out-of-bag R2 (metabolites only): %.3f", oob_r2))

write.csv(tibble(participant_id = d_bun$participant_id,
                 observed_bun = obs, oob_predicted_bun = pred),
          file.path(output_dir, "Fig6B_BUN_observed_vs_oob_predicted.csv"), row.names = FALSE)

p6b <- ggplot(tibble(obs = obs, pred = pred), aes(obs, pred)) +
  geom_point(color = "black", alpha = 0.8, size = 2.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +
  labs(subtitle = sprintf("R\u00b2 = %.2f", oob_r2),
       x = "Observed BUN (mg/dL)", y = "Predicted BUN (mg/dL)") +
  theme_classic(base_size = BASE_SIZE) +
  box_theme +
  theme(axis.text = element_text(colour = "black"),
        plot.subtitle = element_text(size = BASE_SIZE, colour = "black"))

# ---- Panel 6C -------------------------------------------------
origin_lookup <- read_csv(origin_path, show_col_types = FALSE) %>%
  transmute(label  = trimws(metabolite),
            origin = tolower(trimws(origin))) %>%
  filter(!is.na(label), label != "")

bad_origin <- setdiff(unique(origin_lookup$origin), names(origin_colors))
if (length(bad_origin) > 0) {
  stop("Unexpected origin values in ", basename(origin_path), ": ", paste(bad_origin, collapse = ", "))
}
origin_lookup <- distinct(origin_lookup, label, .keep_all = TRUE)
message("Loaded metabolite origin lookup (", nrow(origin_lookup), " entries)")

imp_all <- tibble(feature = names(rf_bun$variable.importance),
                  importance = as.numeric(rf_bun$variable.importance)) %>%
  left_join(feature_map, by = "feature") %>%
  arrange(desc(importance)) %>%
  mutate(rank = row_number())

write.csv(imp_all %>% select(rank, label, metabolite, importance),
          file.path(output_dir, "Fig6C_BUN_importance_all.csv"),
          row.names = FALSE)

short_labels <- c(
  "N,N,N-trimethyl-alanylproline betaine (TMAP)"     = "TMAP",
  "2,3-dihydroxy-5-methylthio-4-pentenoate (DMTPA)*" = "DMTPA*"
)

top_imp <- imp_all %>%
  mutate(label = trimws(label)) %>%
  filter(!grepl("^X\\s*-\\s*[0-9]", label)) %>%   # drop unnamed Metabolon compounds
  slice_head(n = TOP_N) %>%
  left_join(origin_lookup, by = "label") %>%
  mutate(origin = ifelse(is.na(origin), "unclassified", origin))

n_uncl <- sum(top_imp$origin == "unclassified")
if (n_uncl > 0) {
  message("NOTE: ", n_uncl, " of the top-", TOP_N, " metabolites have no origin assignment:")
  message("  ", paste(top_imp$label[top_imp$origin == "unclassified"], collapse = "\n  "))
}

write.csv(top_imp %>% select(rank, label, importance, origin),
          file.path(output_dir, "Fig3C_top_BUN_with_origin.csv"),
          row.names = FALSE)

plot_c <- top_imp %>%
  mutate(display = ifelse(label %in% names(short_labels), unname(short_labels[label]), label))
plot_c$display <- factor(plot_c$display, levels = plot_c$display[order(plot_c$importance)])

legend_levels <- names(origin_colors)
if (!any(plot_c$origin == "unclassified")) legend_levels <- setdiff(legend_levels, "unclassified")

legend_inside <- if (utils::packageVersion("ggplot2") >= "3.5.0") {
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1))
} else {
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
}

p6c <- ggplot(plot_c, aes(display, importance, fill = origin)) +
  geom_col(width = 0.9, key_glyph = "point", show.legend = TRUE) +
  scale_fill_manual(values = origin_colors, limits = legend_levels,
                    name = "Metabolite origin", drop = FALSE) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5.5, colour = NA))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Permutation importance") +
  theme_classic(base_size = C_BASE_SIZE) +
  box_theme +
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = C_LABEL),
        legend.text = element_text(size = C_BASE_SIZE - 2),
        legend.background = element_blank(),
        legend.key = element_blank()) +
  legend_inside

# ================================================================
# 6. SAVE PANELS AND COMBINED FIGURE (PNG 300 dpi, PDF, TIFF LZW)
# ================================================================
save_all <- function(plot, stem, width, height) {
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300)
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(output_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         dpi = TIFF_DPI, compression = "lzw", bg = "white")
}

tag <- function(p, letter) p + labs(tag = letter) +
  theme(plot.tag = element_text(face = "bold", size = TAG_SIZE))

save_all(tag(p6a, "A"), "Figure6A_CV_R2",             width = 9,    height = 6.5)
save_all(tag(p6b, "B"), "Figure6B_BUN_obs_pred",      width = 6.5,  height = 6.5)
save_all(tag(p6c, "C"), "Figure6C_top_BUN_by_origin", width = 15.5, height = 9)
