################################################################
# Figure 3: Plasma metabolite prediction of serum filtration markers
#
# 3A: Five-fold cross-validated R2 (baseline / metabolites / combined)
#     for BUN, creatinine, cystatin C                        
# 3B: Random forest out-of-bag observed vs. predicted, BUN   
# 3C: Top 20 permutation-importance metabolites for BUN,
#     colored by metabolite origin                          
#
################################################################

# ---- Packages -------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
  library(janitor)
})

set.seed(1)

# ================================================================
# EDIT THIS ONE LINE to point to the folder containing the inputs
# ================================================================
base_dir <- "."   # e.g. "D:/UC Berkeley/KPMP paper/repo/Figure3"

input_dir  <- base_dir
output_dir <- file.path(base_dir, "Figure3_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Palettes -------------------------------------------------
# Panel 3A bar chart (three model configurations)
bar_colors <- c(
  "Baseline clinical"      = "#BDBDBD",  # light gray
  "Baseline + metabolites" = "#2C6E7F",  # teal
  "Metabolites only"       = "#8C6D5C"   # brown
)

# Panel 3C metabolite origin (consistent with Figures 4 & 5)
origin_colors <- c(
  "diet-derived"         = "#E66101",  # orange
  "host-derived"         = "#9590C8",  # purple
  "microbiome-dependent" = "#1B9E77",  # teal
  "xenobiotic"           = "#999999"   # gray
)

# ---- Load data ------------------------------------------------
merged_path <- file.path(input_dir, "plasma_merged_wide.csv")
if (!file.exists(merged_path)) stop("File not found: ", merged_path)

dat_raw <- read_csv(merged_path, show_col_types = FALSE)
message("Loaded: ", nrow(dat_raw), " participants, ", ncol(dat_raw), " columns")

# clean-name -> original-name map (for importance labels)
name_map <- tibble(
  clean    = janitor::make_clean_names(names(dat_raw)),
  original = names(dat_raw)
)

dat <- janitor::clean_names(dat_raw)

# ---- Column bookkeeping --------------------------------------
serum_cols <- c(
  "participant_id", "enrollment_category", "sex",
  "age_years_binned", "race", "sample_id",
  "baseline_egfr_ml_min_1_73m2",
  "egfrcr", "egfrcys", "egfrcrcys",
  "albumin_g_d_l", "alkaline_phosphatase_u_l",
  "alanine_aminotransferase_u_l", "aspartate_aminotransferase_u_l",
  "carbon_dioxide_m_eq_l", "total_bilirubin_mg_d_l",
  "direct_bilirubin_mg_d_l", "blood_urea_nitrogen_mg_d_l",
  "creatinine_mg_d_l", "glucose_mg_d_l", "phosphate_mg_d_l",
  "c_reactive_protein_high_sensitivity_mg_l",
  "magnesium_mg_d_l", "total_protein_g_d_l", "uric_acid_mg_d_l",
  "sodium_m_eq_l", "potassium_m_eq_l", "chloride_m_eq_l",
  "calcium_mg_d_l", "cystatin_c_mg_l"
)

outcome_creat <- "creatinine_mg_d_l"
outcome_cys   <- "cystatin_c_mg_l"
outcome_bun   <- "blood_urea_nitrogen_mg_d_l"

metabolite_cols <- setdiff(names(dat), serum_cols)
message("Plasma metabolites: ", length(metabolite_cols))

baseline_cols <- c("age_years_binned", "sex", "baseline_egfr_ml_min_1_73m2")
dat$sex                         <- as.factor(dat$sex)
dat$age_years_binned            <- as.factor(dat$age_years_binned)
dat$baseline_egfr_ml_min_1_73m2 <- suppressWarnings(as.numeric(dat$baseline_egfr_ml_min_1_73m2))

for (o in c(outcome_creat, outcome_cys, outcome_bun))
  dat[[o]] <- suppressWarnings(as.numeric(dat[[o]]))

for (m in metabolite_cols) {
  v <- suppressWarnings(as.numeric(dat[[m]]))
  med <- median(v, na.rm = TRUE); if (!is.finite(med)) med <- 0
  v[is.na(v)] <- med
  dat[[m]] <- v
}

# ================================================================
# PANEL 3A: FIVE-FOLD CROSS-VALIDATED R2  (computed)
# Pooled held-out predictions: R2 = 1 - SS_res / SS_tot
# ================================================================
cv_r2 <- function(dat, outcome, features, k = 5, num.trees = 1000, seed = 1) {
  set.seed(seed)
  d <- dat[, c(outcome, features), drop = FALSE]
  d <- d[complete.cases(d[[outcome]]), , drop = FALSE]
  n <- nrow(d)
  if (n < 20) stop("Too few rows for CV: ", n)
  fold <- sample(rep_len(1:k, n))
  obs  <- d[[outcome]]
  pred <- rep(NA_real_, n)
  for (i in seq_len(k)) {
    tr <- d[fold != i, , drop = FALSE]
    te <- d[fold == i, , drop = FALSE]
    fit <- ranger(x = tr[, features, drop = FALSE], y = tr[[outcome]], num.trees = num.trees)
    pred[fold == i] <- predict(fit, te[, features, drop = FALSE])$predictions
  }
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
}

message("\nComputing 5-fold CV R2 (9 models) ...")
cfg <- list(
  "Baseline clinical"      = baseline_cols,
  "Baseline + metabolites" = c(baseline_cols, metabolite_cols),
  "Metabolites only"       = metabolite_cols
)
outcomes <- c(BUN = outcome_bun, Creatinine = outcome_creat, `Cystatin C` = outcome_cys)

cv_results <- purrr::imap_dfr(outcomes, function(ocol, oname) {
  purrr::imap_dfr(cfg, function(feats, mname) {
    r2 <- cv_r2(dat, ocol, feats)
    message(sprintf("  %-12s | %-24s CV R2 = %.3f", oname, mname, r2))
    tibble(outcome = oname, model = mname, cv_r2 = r2)
  })
})
write.csv(cv_results, file.path(output_dir, "Fig3A_CV_R2_computed.csv"), row.names = FALSE)

cv_results$model   <- factor(cv_results$model, levels = names(bar_colors))
cv_results$outcome <- factor(cv_results$outcome, levels = c("BUN", "Creatinine", "Cystatin C"))

p3a <- ggplot(cv_results, aes(outcome, cv_r2, fill = model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = bar_colors) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "A", x = NULL, y = expression(CV~R^2), fill = "Model") +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 13))

ggsave(file.path(output_dir, "Fig3A_CV_R2_barplot.png"), p3a, width = 8, height = 6, dpi = 300)
message("Saved Fig 3A")

# ================================================================
# PANEL 3B & 3C: BUN out-of-bag RF (scatter + permutation importance)
# ================================================================
rf_bun <- local({
  d <- dat[, c(outcome_bun, metabolite_cols), drop = FALSE]
  d <- d[complete.cases(d[[outcome_bun]]), , drop = FALSE]
  fit <- ranger(x = d[, metabolite_cols, drop = FALSE], y = d[[outcome_bun]],
                num.trees = 1000, importance = "permutation", sample.fraction = 0.632)
  obs <- d[[outcome_bun]]; pred <- fit$predictions
  list(fit = fit, obs = obs, pred = pred,
       r2 = cor(obs, pred, use = "complete.obs")^2,
       rmse = sqrt(mean((obs - pred)^2, na.rm = TRUE)))
})
message(sprintf("BUN OOB: R2=%.3f RMSE=%.3f", rf_bun$r2, rf_bun$rmse))

# ---- Panel 3B: observed vs predicted -------------------------
p3b <- ggplot(tibble(obs = rf_bun$obs, pred = rf_bun$pred), aes(obs, pred)) +
  geom_point(color = "black", alpha = 0.8, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "B",
       subtitle = sprintf("R\u00b2 = %.2f", rf_bun$r2),
       x = "Observed BUN (mg/dL)", y = "Predicted BUN (mg/dL)") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12, color = "gray30"))

ggsave(file.path(output_dir, "Fig3B_RF_BUN_obs_pred.png"), p3b, width = 6, height = 5, dpi = 300)
message("Saved Fig 3B")

# ---- Panel 3C: top 20 importance, colored by origin ----------
# Optional origin lookup: metabolite_origin.csv with columns metabolite, origin
origin_path <- file.path(input_dir, "metabolite_origin.csv")
origin_lookup <- NULL
if (file.exists(origin_path)) {
  origin_lookup <- read_csv(origin_path, show_col_types = FALSE)
  oc <- intersect(c("metabolite", "Metabolite", "chemical_name"), names(origin_lookup))[1]
  gc <- intersect(c("origin", "metabolite_origin", "Origin"), names(origin_lookup))[1]
  origin_lookup <- origin_lookup %>% transmute(label = .data[[oc]], origin = .data[[gc]])
  message("Loaded metabolite origin lookup (", nrow(origin_lookup), " entries)")
} else {
  message("NOTE: metabolite_origin.csv not found - panel 3C bars shown in single color. ",
          "Add it (columns: metabolite, origin) to color by origin.")
}

imp <- rf_bun$fit$variable.importance
imp_df <- tibble(clean = names(imp), importance = as.numeric(imp)) %>%
  left_join(name_map, by = "clean") %>%
  mutate(label = coalesce(original, clean)) %>%
  filter(!grepl("^X-", label)) %>%
  arrange(desc(importance)) %>%
  slice_head(n = 20)

if (!is.null(origin_lookup)) {
  imp_df <- imp_df %>% left_join(origin_lookup, by = "label") %>%
    mutate(origin = ifelse(is.na(origin), "host-derived", origin))
} else {
  imp_df$origin <- "host-derived"
}

imp_df <- imp_df %>% mutate(label = factor(label, levels = rev(label)))

p3c <- ggplot(imp_df, aes(importance, label, fill = origin)) +
  geom_col() +
  scale_fill_manual(values = origin_colors, name = "Metabolite origin") +
  labs(title = "C", x = "Permutation importance\n(% increase in MSE)", y = NULL) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 10),
        legend.position = "top")

ggsave(file.path(output_dir, "Fig3C_top20_BUN_by_origin.png"), p3c, width = 9, height = 7, dpi = 300)
message("Saved Fig 3C")

message("\n========== Figure 3 complete ==========")
message("Outputs in: ", output_dir)
