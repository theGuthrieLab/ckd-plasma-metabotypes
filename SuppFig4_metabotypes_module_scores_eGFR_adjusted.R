#!/usr/bin/env Rscript

################################################################
# Supplemental Figure 4 (Reviewer 2, comment 3)
# Transcriptional module scores by discovery-cohort plasma metabotype before and after adjustment for eGFR
#
# Uses the discovery-cohort (Figure 1, k = 2) metabotype labels and the
# participant-level module scores from Figure 5D, fitted under the
# same three specifications:
#   (1) Unadjusted:       score ~ metabotype
#   (2) Age + sex:        score ~ metabotype + age + sex
#   (3) Age + sex + eGFR: score ~ metabotype + age + sex + eGFR
#
# For each module and model it reports the overall metabotype effect
# (nested F-test; with two metabotypes this equals the M2 vs M1 test)
# and the M2 vs M1 difference with a 95% CI. BH-adjusted P values
# across the six modules are added within each model.
#
# Age is supplied in 10-year bins and converted to the bin midpoint.
#
# INPUTS (paths set below)
#   Cluster_assignments_Q1.csv                      (Figure 1 labels)
#   Clinical_CKD_HealthyReference_only.xlsx
#   Figure5D_molecular_programs_plotting_table.csv  (module scores)
################################################################

options(stringsAsFactors = FALSE)

req_pkgs <- c("tidyverse", "readxl", "broom")
missing_pkgs <- req_pkgs[
  !vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(broom)
})

# ================================================================
# EDIT THESE PATHS
# ================================================================
base_dir <- "."   # folder that contains data/, e.g. "D:/UC Berkeley/KPMP paper/repo"

data_dir     <- file.path(base_dir, "data")
cluster_path <- file.path(data_dir, "Cluster_assignments_Q1.csv")
clin_path    <- file.path(data_dir, "Clinical_CKD_HealthyReference_only.xlsx")
module_path  <- file.path(data_dir, "Figure5D_molecular_programs_plotting_table.csv")


out_dir <- file.path(base_dir, "SuppFig8_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (p in c(cluster_path, clin_path, module_path)) {
  if (!file.exists(p)) stop("Input not found: ", p)
}

module_levels <- c("AhR_response", "Fibrosis_ECM", "Inflammation",
                   "Oxidative_stress", "Podocyte_identity", "Tubular_injury")

# ================================================================
# HELPERS 
# ================================================================
clean_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "\\s+", "")
  stringr::str_replace(x, "^g(?=[0-9]{1,3}-[0-9]{4,6}$)", "")
}

bin_midpoint <- function(x) {
  nums <- stringr::str_match(as.character(x), "(\\d+)\\s*-\\s*(\\d+)")
  rowMeans(cbind(as.numeric(nums[, 2]), as.numeric(nums[, 3])))
}

pick_col <- function(df, candidates, label) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    stop("Could not find the ", label, " column. Columns present:\n  ",
         paste(names(df), collapse = "\n  "))
  }
  hit[1]
}

# ================================================================
# 1. FIGURE 1 METABOTYPE ASSIGNMENTS
# ================================================================
clus_raw <- read_csv(cluster_path, show_col_types = FALSE)

id_col  <- pick_col(clus_raw, c("participant_id", "ParticipantID",
                                "Participant ID"), "participant ID")
mt_col  <- pick_col(clus_raw, c("metabotype", "Metabotype", "cluster"),
                    "metabotype")
grp_col <- pick_col(clus_raw, c("group", "Group", "Enrollment Category"),
                    "group")

clusters <- clus_raw %>%
  transmute(participant_id = clean_id(.data[[id_col]]),
            metabotype     = as.character(.data[[mt_col]]),
            group          = as.character(.data[[grp_col]])) %>%
  filter(!is.na(metabotype), metabotype != "") %>%
  distinct(participant_id, .keep_all = TRUE)

lev <- sort(unique(clusters$metabotype))
if (all(grepl("^[0-9]+$", lev))) {
  lev <- as.character(sort(as.numeric(lev)))
  clusters$metabotype <- factor(paste0("M", clusters$metabotype),
                                levels = paste0("M", lev))
} else {
  clusters$metabotype <- factor(clusters$metabotype, levels = lev)
}

# ================================================================
# 2. CLINICAL COVARIATES
# ================================================================
clin_raw <- read_excel(clin_path)
egfr_cols <- grep("^Baseline eGFR", names(clin_raw), value = TRUE)
egfr_col  <- egfr_cols[!grepl("Binned", egfr_cols)][1]
if (is.na(egfr_col)) stop("Could not find an unbinned baseline eGFR column.")

clin <- clin_raw %>%
  transmute(
    participant_id = clean_id(`Participant ID`),
    age            = bin_midpoint(`Age (Years) (Binned)`),
    sex            = factor(Sex),
    egfr           = suppressWarnings(as.numeric(.data[[egfr_col]]))
  ) %>%
  distinct(participant_id, .keep_all = TRUE)

# ================================================================
# 3. MODULE SCORES 
# ================================================================
mod_raw <- read_csv(module_path, show_col_types = FALSE)
mod_id  <- pick_col(mod_raw, c("SampleID", "ParticipantID", "participant_id"),
                    "module-table participant ID")

modules <- mod_raw %>%
  transmute(participant_id = clean_id(.data[[mod_id]]),
            outcome        = as.character(module),
            value          = as.numeric(module_score))

# ================================================================
# 4. MERGE
# ================================================================
n_mod <- n_distinct(modules$participant_id)
dat <- modules %>%
  inner_join(clusters, by = "participant_id") %>%
  left_join(clin, by = "participant_id") %>%
  mutate(outcome = factor(outcome, levels = module_levels))

n_matched <- n_distinct(dat$participant_id)
message("Participants with module scores: ", n_mod,
        " | matched to Figure 1 metabotypes: ", n_matched)
if (n_matched < n_mod) {
  warning(n_mod - n_matched, " module-score participants have no Figure 1 ",
          "label; check participant IDs.")
}

participants <- dat %>% distinct(participant_id, .keep_all = TRUE)
message("\nGroup x Figure 1 metabotype (participants with module scores):")
print(table(participants$group, participants$metabotype, useNA = "ifany"))
message("Missing eGFR: ",
        paste(names(table(participants$group[is.na(participants$egfr)])),
              table(participants$group[is.na(participants$egfr)]),
              collapse = ", "))

write.csv(participants %>% select(participant_id, group, metabotype, age, sex, egfr),
          file.path(out_dir, "modules_analysis_sample.csv"),
          row.names = FALSE)


complete_cc <- function(d) d %>% filter(!is.na(age), !is.na(sex), !is.na(egfr))
samples <- list(
  "All participants" = complete_cc(dat)
)
for (s in names(samples)) {
  pp <- samples[[s]] %>% distinct(participant_id, .keep_all = TRUE)
  message(s, ": n = ", nrow(pp), " (",
          paste(names(table(pp$metabotype)), table(pp$metabotype),
                sep = " = ", collapse = ", "), ")")
}

# ================================================================
# 5. MODELS
# ================================================================
model_specs <- list(
  "Unadjusted"       = character(0),
  "Age + sex"        = c("age", "sex"),
  "Age + sex + eGFR" = c("age", "sex", "egfr")
)

fit_one <- function(d, outcome_name, spec_name, covs) {
  dd <- d %>% filter(outcome == outcome_name) %>%
    select(all_of(c("value", "metabotype", covs)))
  dd$metabotype <- droplevels(dd$metabotype)
  if (nrow(dd) < 10 || nlevels(dd$metabotype) < 2) {
    return(list(coefs = NULL, test = NULL))
  }

  full <- lm(as.formula(paste("value ~",
                              paste(c("metabotype", covs), collapse = " + "))),
             data = dd)
  reduced <- if (length(covs) == 0) lm(value ~ 1, data = dd) else
    lm(as.formula(paste("value ~", paste(covs, collapse = " + "))), data = dd)
  av <- anova(reduced, full)

  test <- tibble(outcome = outcome_name, model = spec_name, n = nrow(dd),
                 df_num = av$Df[2], df_den = av$Res.Df[2],
                 F_statistic = av$F[2], p_metabotype = av$`Pr(>F)`[2])

  coefs <- broom::tidy(full, conf.int = TRUE) %>%
    filter(grepl("^metabotype", term)) %>%
    transmute(outcome = outcome_name, model = spec_name, n = nrow(dd),
              contrast = paste0(sub("^metabotype", "", term), " vs ",
                                levels(dd$metabotype)[1]),
              estimate, conf.low, conf.high, std.error, p.value)

  list(coefs = coefs, test = test)
}

run_sample <- function(d, label) {
  res <- purrr::map(module_levels, function(o) {
    purrr::imap(model_specs, function(covs, spec) fit_one(d, o, spec, covs))
  })
  flat <- unlist(res, recursive = FALSE)
  list(coefs = bind_rows(lapply(flat, `[[`, "coefs")) %>% mutate(sample = label),
       tests = bind_rows(lapply(flat, `[[`, "test"))  %>% mutate(sample = label))
}

message("\nFitting models ...")
all_res <- purrr::imap(samples, run_sample)

coefs_all <- bind_rows(lapply(all_res, `[[`, "coefs")) %>%
  group_by(sample, model) %>%
  mutate(p_BH_across_modules = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  mutate(sample = factor(sample, levels = names(samples))) %>%
  select(sample, everything())

tests_all <- bind_rows(lapply(all_res, `[[`, "tests")) %>%
  group_by(sample, model) %>%
  mutate(p_BH_across_modules = p.adjust(p_metabotype, method = "BH")) %>%
  ungroup() %>%
  mutate(sample = factor(sample, levels = names(samples)),
         model  = factor(model, levels = names(model_specs)),
         outcome = factor(outcome, levels = module_levels)) %>%
  select(sample, everything()) %>%
  arrange(sample, outcome, model)

write.csv(coefs_all, file.path(out_dir, "modules_adjusted.csv"),
          row.names = FALSE)
write.csv(tests_all, file.path(out_dir, "modules_overall_tests.csv"),
          row.names = FALSE)

message("\nOverall metabotype effect on module scores (nested F-test):")
print(as.data.frame(tests_all %>%
  mutate(F_statistic = round(F_statistic, 2),
         p_metabotype = signif(p_metabotype, 3),
         p_BH_across_modules = signif(p_BH_across_modules, 3))))

# ================================================================
# 6. FOREST PLOTS 
# ================================================================
pal <- c("M2 vs M1" = "#B2182B")

plot_forest <- function(label) {
  pd <- coefs_all %>%
    filter(sample == label) %>%
    mutate(model   = factor(model, levels = rev(names(model_specs))),
           outcome = factor(outcome, levels = module_levels))
  n_lab <- unique(pd$n)

  p <- ggplot(pd, aes(estimate, model, color = contrast)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
    facet_wrap(~ outcome, scales = "free_x", ncol = 3) +
    scale_color_manual(values = pal, name = "Contrast") +
    labs(title = "Transcriptional module scores by plasma metabotype, before and after adjustment",
         subtitle = paste0("n = ", paste(n_lab, collapse = "/")),
         x = "Difference vs M1 (95% CI)", y = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  stem <- paste0("modules_forest_", gsub(" ", "_", label))
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = 11, height = 6.5, dpi = 300)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = 11, height = 6.5)
}

for (s in names(samples)) plot_forest(s)

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_modules.txt"))

message("\n========== Done ==========")
message("Outputs in: ", out_dir)
