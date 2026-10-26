############################################################################################################
# Figure 2: Metabotype score, eGFR associations, and pathway enrichment
#
# 2A: Association between metabotype score and baseline eGFR
# 2B: Top plasma metabolites associated with eGFR, colored by origin
# 2C: Pathway enrichment of metabolites ranked by correlation
#     with metabotype score
# 2D: Structural kidney injury by discovery-cohort plasma metabotype before and after adjustment for eGFR
#
############################################################################################################

# ---- Packages -------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(readr)
  library(janitor)
  library(stringr)
  library(fgsea)
  library(forcats)
  library(cowplot)
})

set.seed(1)

# ================================================================
# EDIT THIS ONE LINE to point to the folder containing the inputs
# ================================================================
base_dir <- "."   # e.g. "D:/UC Berkeley/KPMP paper/repo/Figure2"

input_dir  <- file.path(base_dir, "data")
output_dir <- file.path(base_dir, "Figure2_output")
figure_dir <- file.path(base_dir, "figuresFinal")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helper functions -----------------------------------------

clean_id <- function(x) {
  str_replace_all(str_trim(as.character(x)), "\\s+", "")
}

theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
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

# ---- Palettes -------------------------------------------------

group_colors <- c(
  "CKD" = "#fc8d62",
  "Healthy Reference" = "#66c2a5"
)

origin_colors <- c(
  "microbiome-dependent" = "#1b9e77",
  "host-derived"         = "#7570b3",
  "diet-derived"         = "#d95f02",
  "xenobiotic"           = "#666666",
  "unclassified"         = "#999999"
)

fgsea_colors <- c(
  "FALSE" = "#D3D3D3",
  "TRUE"  = "#D3D3D3"
)

# ---- Input files ----------------------------------------------

clinical_path <- file.path(input_dir, "Clinical_CKD_HealthyReference_only.xlsx")
clusters_path <- file.path(input_dir, "Cluster_assignments_Q1.csv")
scores_path   <- file.path(input_dir, "KPMP_CV_predictions_and_score.csv")
panel_path    <- file.path(input_dir, "KPMP_panel_lambda_1se.csv")
plasma_path   <- file.path(input_dir, "Plasma_Metabolomics.xlsx")
annot_path    <- file.path(input_dir, "Annotation-chem-plasma-metabolites.xlsx")

required_files <- c(
  clinical_path,
  clusters_path,
  scores_path,
  panel_path,
  plasma_path,
  annot_path
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

# ---- Load data ------------------------------------------------

clinical <- read_excel(clinical_path) %>%
  clean_names()

clusters <- read_csv(clusters_path, show_col_types = FALSE) %>%
  clean_names()

scores <- read_csv(scores_path, show_col_types = FALSE) %>%
  clean_names()

kpmp_panel <- read_csv(panel_path, show_col_types = FALSE)

plasma <- read_excel(plasma_path) %>%
  clean_names()

annot_raw <- read_excel(annot_path) %>%
  clean_names() %>%
  rename(metabolite = chem_id)

message("Loaded clinical rows: ", nrow(clinical))
message("Loaded cluster rows: ", nrow(clusters))
message("Loaded score rows: ", nrow(scores))
message("Loaded plasma rows: ", nrow(plasma))
message("Loaded annotation rows: ", nrow(annot_raw))

# ================================================================
# PREPARE SHARED DATA OBJECTS
# ================================================================

plasma_df <- plasma %>%
  rename(
    metabolite = metabolite,
    value = abundance,
    metabolite_name = chemical_name
  ) %>%
  mutate(
    participant_id = clean_id(participant_id),
    metabolite = as.character(metabolite)
  )

met_annotation <- annot_raw %>%
  transmute(
    metabolite = as.character(metabolite),
    metabolite_name = plot_name,
    super_pathway = super_pathway,
    superclass = superclass,
    class = class
  ) %>%
  distinct(metabolite, .keep_all = TRUE)

df <- clinical %>%
  mutate(participant_id = clean_id(participant_id)) %>%
  inner_join(
    clusters %>%
      mutate(participant_id = clean_id(participant_id)),
    by = "participant_id"
  ) %>%
  left_join(
    scores %>%
      mutate(participant_id = clean_id(participant_id)),
    by = "participant_id"
  ) %>%
  rename(egfr = baseline_e_gfr_ml_min_1_73m2)

message("Merged clinical/cluster/score rows: ", nrow(df))

# ================================================================
# PANEL 2A: Metabotype score vs. baseline eGFR
# ================================================================

cor_res <- cor.test(
  df$metabotype_score,
  df$egfr,
  method = "spearman"
)

fig2A <- ggplot(
  df,
  aes(
    x = metabotype_score,
    y = egfr
  )
) +
  geom_point(
    aes(color = group),
    size = 2.8,
    alpha = 0.9
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0(
      "rho = ",
      round(cor_res$estimate, 2)
    ),
    hjust = 1.1,
    vjust = 1.1
  ) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "A",
    x = "Metabotype score",
    y = "eGFR",
    color = "Group"
  ) +
  theme_cowplot() +
  panel_border(color = "black")

save_tiff(
  fig2A,
  "Fig2A.tiff",
  width = 6,
  height = 4
)

message("Saved Fig 2A")

# ================================================================
# PANEL 2B: Top eGFR-associated metabolites
# ================================================================
# Uses the top 50 classifier metabolites by absolute coefficient,
# then ranks those metabolites by association with eGFR.

top_mets <- kpmp_panel %>%
  slice_max(abs(coefficient), n = 50) %>%
  mutate(metabolite = as.character(metabolite)) %>%
  pull(metabolite)

plasma_subset <- plasma_df %>%
  filter(metabolite %in% top_mets)

res_met <- plasma_subset %>%
  left_join(df, by = "participant_id") %>%
  group_by(metabolite) %>%
  group_modify(~{
    tmp <- .x %>%
      filter(!is.na(value), !is.na(egfr))

    if (nrow(tmp) < 3) {
      return(tibble(beta = NA_real_, p = NA_real_))
    }

    fit <- lm(egfr ~ value, data = tmp)
    coef_tbl <- summary(fit)$coefficients

    tibble(
      beta = coef_tbl[2, 1],
      p = coef_tbl[2, 4]
    )
  }) %>%
  ungroup()

top_met <- res_met %>%
  filter(!is.na(p)) %>%
  arrange(p) %>%
  slice_head(n = 20) %>%
  left_join(met_annotation, by = "metabolite")

write.csv(
  top_met,
  file.path(output_dir, "top_metabolite_associations_egfr.csv"),
  row.names = FALSE
)

# Optional manually classified metabolite-origin table.
# If present, it should contain the columns:
# metabolite_name, beta, metabolite_origin
classified_path <- file.path(output_dir, "top_metabolite_associations_egfr_classified.csv")

if (file.exists(classified_path)) {
  top_met_plot <- read_csv(classified_path, show_col_types = FALSE) %>%
    clean_names()
  message("Loaded manually classified metabolite-origin table")
} else {
  top_met_plot <- top_met %>%
    mutate(metabolite_origin = "unclassified")
  message(
    "No classified metabolite-origin table found. ",
    "Using generated top_met table with origin set to 'unclassified'."
  )
}

fig2B <- ggplot(
  top_met_plot,
  aes(
    x = beta,
    y = reorder(metabolite_name, beta),
    fill = metabolite_origin
  )
) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_fill_manual(values = origin_colors) +
  labs(
    title = "B",
    x = "Association with eGFR",
    y = NULL,
    fill = "Metabolite origin"
  ) +
  theme_cowplot() +
  panel_border(color = "black") +
  theme(legend.position = "bottom")

save_tiff(
  fig2B,
  "Fig2B.tiff",
  width = 4,
  height = 5
)

message("Saved Fig 2B")

# ================================================================
# PANEL 2C: Pathway enrichment by metabotype-score correlation
# ================================================================
# Metabolites are z-scaled across participants, ranked by Pearson
# correlation with the metabotype score, and used for fgsea.

fgsea_df <- plasma %>%
  rename(value = abundance) %>%
  mutate(
    participant_id = clean_id(participant_id),
    metabolite = as.character(metabolite)
  ) %>%
  left_join(
    scores %>%
      mutate(participant_id = clean_id(participant_id)),
    by = "participant_id"
  ) %>%
  left_join(
    clusters %>%
      mutate(participant_id = clean_id(participant_id)),
    by = "participant_id"
  ) %>%
  group_by(metabolite) %>%
  mutate(value = as.numeric(scale(value))) %>%
  ungroup()

res_ranked <- fgsea_df %>%
  group_by(metabolite) %>%
  group_modify(~{
    tmp <- .x %>%
      filter(!is.na(value), !is.na(metabotype_score))

    if (nrow(tmp) < 3) {
      return(tibble(stat = NA_real_))
    }

    tibble(
      stat = cor(
        tmp$value,
        tmp$metabotype_score,
        use = "complete.obs"
      )
    )
  }) %>%
  ungroup() %>%
  filter(!is.na(stat))

ranked_stats <- res_ranked$stat
names(ranked_stats) <- res_ranked$metabolite
ranked_stats <- sort(ranked_stats, decreasing = TRUE)

pathways_super <- met_annotation %>%
  filter(!is.na(super_pathway)) %>%
  split(.$super_pathway) %>%
  purrr::map("metabolite")

pathways_class <- met_annotation %>%
  filter(!is.na(class)) %>%
  split(.$class) %>%
  purrr::map("metabolite")

pathways_superclass <- met_annotation %>%
  filter(!is.na(superclass)) %>%
  split(.$superclass) %>%
  purrr::map("metabolite")

fgsea_super <- fgsea(
  pathways = pathways_super,
  stats = ranked_stats
)

fgsea_class <- fgsea(
  pathways = pathways_class,
  stats = ranked_stats
)

fgsea_superclass <- fgsea(
  pathways = pathways_superclass,
  stats = ranked_stats
)

write.csv(
  as_tibble(fgsea_super),
  file.path(output_dir, "Fig2C_fgsea_super_pathway.csv"),
  row.names = FALSE
)

write.csv(
  as_tibble(fgsea_class),
  file.path(output_dir, "Fig2C_fgsea_class.csv"),
  row.names = FALSE
)

write.csv(
  as_tibble(fgsea_superclass),
  file.path(output_dir, "Fig2C_fgsea_superclass.csv"),
  row.names = FALSE
)

plot_fgsea <- function(res, title, n = 12) {
  res %>%
    as_tibble() %>%
    filter(!is.na(padj), !is.na(NES)) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    arrange(padj) %>%
    slice_head(n = n) %>%
    mutate(
      pathway = as.character(pathway),
      pathway = forcats::fct_reorder(pathway, NES)
    ) %>%
    ggplot(aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_col(width = 0.7) +
    geom_vline(
      xintercept = 0,
      linetype = 2,
      color = "grey60"
    ) +
    scale_fill_manual(values = fgsea_colors) +
    labs(
      title = title,
      x = "Enrichment Score",
      y = NULL
    ) +
    theme_cowplot() +
    panel_border(color = "black") +
    theme(legend.position = "none")
}

# Main Figure 2C panel: class-level enrichment
fig2C <- plot_fgsea(
  fgsea_class,
  title = "C"
)

save_tiff(
  fig2C,
  "Fig2C.tiff",
  width = 5,
  height = 5
)

message("Saved Fig 2C")

# Optional supplementary pathway-enrichment panels
fig2C_super_pathway <- plot_fgsea(
  fgsea_super,
  title = "Super pathway enrichment"
)

fig2C_superclass <- plot_fgsea(
  fgsea_superclass,
  title = "Superclass enrichment"
)

save_tiff(
  fig2C_super_pathway,
  "Fig2C_super_pathway.tiff",
  width = 5,
  height = 5
)

save_tiff(
  fig2C_superclass,
  "Fig2C_superclass.tiff",
  width = 5,
  height = 5
)

message("\n========== Figure 2 complete ==========")
message("Outputs in: ", output_dir)
message("Figures in: ", figure_dir)

################################################################
# 2D: Structural kidney injury by discovery-cohort plasma metabotype before and after adjustment for eGFR (Reviewer 2, comment 3)
# Figure 1 metabotypes: do they carry information beyond eGFR?
#
# Refits the association between discovery-cohort metabotype and
# structural kidney injury under three specifications:
#   (1) Unadjusted:       outcome ~ metabotype
#   (2) Age + sex:        outcome ~ metabotype + age + sex
#   (3) Age + sex + eGFR: outcome ~ metabotype + age + sex + eGFR
#
# Outcomes: cortical interstitial fibrosis (%), tubular atrophy (%),
# and their sum (IFTA, %).
#
# INPUTS (data/)
#   Cluster_assignments_Q1.csv
#   Clinical_CKD_HealthyReference_only.xlsx
#   KPMP TIV Descriptor Scores 03.30.2025.xlsx   (headers on row 4)
#
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
tiv_path     <- file.path(data_dir, "KPMP TIV Descriptor Scores 03.30.2025.xlsx")

out_dir <- file.path(base_dir, "Fig2D_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TIV_HEADER_ROW <- 4   # headers are on row 4 of the TIV workbook

for (p in c(cluster_path, clin_path, tiv_path)) {
  if (!file.exists(p)) stop("Input not found: ", p)
}

# ================================================================
# HELPERS
# ================================================================
# Strip whitespace and a leading "g" from IDs such as g27-10088
clean_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "\\s+", "")
  stringr::str_replace(x, "^g(?=[0-9]{1,3}-[0-9]{4,6}$)", "")
}

# "60-69 Years" -> 64.5
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
# 1. METABOTYPE ASSIGNMENTS
# ================================================================
clus_raw <- read_csv(cluster_path, show_col_types = FALSE)

id_col  <- pick_col(clus_raw, c("participant_id", "ParticipantID",
                                "Participant ID", "participantid"),
                    "participant ID")
mt_col  <- pick_col(clus_raw, c("metabotype", "Metabotype", "cluster"),
                    "metabotype")
grp_col <- intersect(c("group", "Group", "enrollment_category",
                       "Enrollment Category"), names(clus_raw))

clusters <- clus_raw %>%
  transmute(participant_id = clean_id(.data[[id_col]]),
            metabotype     = as.character(.data[[mt_col]]),
            group          = if (length(grp_col) > 0) {
                               as.character(.data[[grp_col[1]]])
                             } else NA_character_) %>%
  filter(!is.na(metabotype), metabotype != "") %>%
  distinct(participant_id, .keep_all = TRUE)

# Label metabotypes consistently (M1, M2, ... in sorted order)
lev <- sort(unique(clusters$metabotype))
if (all(grepl("^[0-9]+$", lev))) {
  lev <- as.character(sort(as.numeric(lev)))
  clusters$metabotype <- factor(paste0("M", clusters$metabotype),
                                levels = paste0("M", lev))
} else {
  clusters$metabotype <- factor(clusters$metabotype, levels = lev)
}

message("Metabotype assignments: ", nrow(clusters), " participants")
print(table(clusters$metabotype))

# ================================================================
# 2. CLINICAL COVARIATES
# ================================================================
clin_raw <- read_excel(clin_path)

# The workbook contains both "Baseline eGFR (ml/min/1.73m2)" and a
# "(Binned)" version; take the unbinned numeric column.
egfr_cols <- grep("^Baseline eGFR", names(clin_raw), value = TRUE)
egfr_col  <- egfr_cols[!grepl("Binned", egfr_cols)][1]
if (is.na(egfr_col)) stop("Could not find an unbinned baseline eGFR column.")
message("Using eGFR column: ", egfr_col)

clin <- clin_raw %>%
  transmute(
    participant_id = clean_id(`Participant ID`),
    age            = bin_midpoint(`Age (Years) (Binned)`),
    sex            = factor(Sex),
    egfr           = suppressWarnings(as.numeric(.data[[egfr_col]]))
  ) %>%
  distinct(participant_id, .keep_all = TRUE)

message("Clinical rows: ", nrow(clin),
        " | missing eGFR: ", sum(is.na(clin$egfr)),
        " | missing age: ", sum(is.na(clin$age)))

# ================================================================
# 3. HISTOPATHOLOGY (TIV descriptor scores)
# ================================================================
tiv_raw <- read_excel(tiv_path, skip = TIV_HEADER_ROW - 1)
names(tiv_raw) <- names(tiv_raw) %>%
  stringr::str_replace_all("\\s+", " ") %>% stringr::str_trim()

tiv_id   <- pick_col(tiv_raw, c("Participant ID", "ParticipantID"),
                     "TIV participant ID")
fib_col  <- pick_col(tiv_raw, c("Cortex: Interstitial Fibrosis: Percent"),
                     "interstitial fibrosis")
atr_col  <- pick_col(tiv_raw,
                     c("Cortex: Tubular Atrophy: Common Type: Percent"),
                     "tubular atrophy")

tiv <- tiv_raw %>%
  transmute(
    participant_id        = clean_id(.data[[tiv_id]]),
    `Interstitial Fibrosis` = suppressWarnings(as.numeric(.data[[fib_col]])),
    `Tubular Atrophy`       = suppressWarnings(as.numeric(.data[[atr_col]]))
  ) %>%
  mutate(IFTA = `Interstitial Fibrosis` + `Tubular Atrophy`) %>%
  filter(!is.na(`Interstitial Fibrosis`) | !is.na(`Tubular Atrophy`)) %>%
  distinct(participant_id, .keep_all = TRUE)

message("Participants with histopathology: ", nrow(tiv))

# ================================================================
# 4. MERGE
# ================================================================
dat <- clusters %>%
  inner_join(clin, by = "participant_id") %>%
  inner_join(tiv,  by = "participant_id")

message("\nAnalysis sample (metabotype + clinical + pathology): ", nrow(dat))
print(table(dat$metabotype))
message("Missing eGFR within the analysis sample: ", sum(is.na(dat$egfr)))

if (nrow(dat) < 20) {
  stop("Too few participants after merging. Check that participant IDs ",
       "match across the three files.")
}

dat$metabotype <- droplevels(dat$metabotype)
ref_level <- levels(dat$metabotype)[1]
message("Reference metabotype: ", ref_level)

write.csv(dat, file.path(out_dir, "R2C3_Figure1_analysis_sample.csv"),
          row.names = FALSE)

# ---- Analysis samples -----------------------------------------
# Primary: CKD participants only, so the contrast is unambiguously
# within CKD. Sensitivity: all participants with histopathology,
# which may include a small number of healthy reference donors
# (living-donor pre-perfusion biopsies).
if (all(is.na(dat$group))) {
  message("\nNOTE: no group column found; running the full sample only.")
  samples <- list("All participants" = dat)
} else {
  ckd <- dat %>% filter(group == "CKD")
  message("\nGroup composition of the analysis sample:")
  print(table(dat$group, dat$metabotype))
  message("Primary (CKD only): n = ", nrow(ckd),
          " | Sensitivity (all): n = ", nrow(dat))
  samples <- list("CKD only" = ckd, "All participants" = dat)
}

make_long <- function(d) {
  d %>% pivot_longer(cols = c("IFTA", "Interstitial Fibrosis", "Tubular Atrophy"),
                     names_to = "outcome", values_to = "value")
}

# ================================================================
# 5. MODELS
# ================================================================
model_specs <- list(
  "Unadjusted"       = character(0),
  "Age + sex"        = c("age", "sex"),
  "Age + sex + eGFR" = c("age", "sex", "egfr")
)

fit_one <- function(long, outcome_name, spec_name, covs) {
  dd <- long %>%
    filter(outcome == outcome_name) %>%
    select(all_of(c("value", "metabotype", covs))) %>%
    tidyr::drop_na()
  dd$metabotype <- droplevels(dd$metabotype)

  if (nrow(dd) < 10 || nlevels(dd$metabotype) < 2) {
    return(list(coefs = NULL, test = NULL))
  }

  full <- lm(as.formula(paste("value ~",
                              paste(c("metabotype", covs), collapse = " + "))),
             data = dd)
  reduced <- if (length(covs) == 0) {
    lm(value ~ 1, data = dd)
  } else {
    lm(as.formula(paste("value ~", paste(covs, collapse = " + "))), data = dd)
  }

  av <- anova(reduced, full)

  test <- tibble(
    outcome = outcome_name, model = spec_name, n = nrow(dd),
    df_num = av$Df[2], df_den = av$Res.Df[2],
    F_statistic = av$F[2], p_metabotype = av$`Pr(>F)`[2]
  )

  coefs <- broom::tidy(full, conf.int = TRUE) %>%
    filter(grepl("^metabotype", term)) %>%
    transmute(
      outcome = outcome_name, model = spec_name, n = nrow(dd),
      contrast = paste0(sub("^metabotype", "", term), " vs ",
                        levels(dd$metabotype)[1]),
      estimate, conf.low, conf.high, std.error, p.value
    )

  list(coefs = coefs, test = test)
}

message("\nFitting models ...")
run_sample <- function(d, sample_label) {
  long <- make_long(d)
  res <- purrr::map(unique(long$outcome), function(o) {
    purrr::imap(model_specs, function(covs, spec) fit_one(long, o, spec, covs))
  })
  flat <- unlist(res, recursive = FALSE)
  list(
    coefs = bind_rows(lapply(flat, `[[`, "coefs")) %>% mutate(sample = sample_label),
    tests = bind_rows(lapply(flat, `[[`, "test"))  %>% mutate(sample = sample_label)
  )
}

all_res <- purrr::imap(samples, function(d, lab) run_sample(d, lab))

coefs_all <- bind_rows(lapply(all_res, `[[`, "coefs")) %>%
  mutate(sample = factor(sample, levels = names(samples))) %>%
  select(sample, everything())
tests_all <- bind_rows(lapply(all_res, `[[`, "tests")) %>%
  mutate(sample = factor(sample, levels = names(samples)),
         model  = factor(model, levels = names(model_specs))) %>%
  select(sample, everything()) %>%
  arrange(sample, outcome, model)

write.csv(coefs_all,
          file.path(out_dir, "R2C3_Figure1_structural_adjusted.csv"),
          row.names = FALSE)
write.csv(tests_all,
          file.path(out_dir, "R2C3_Figure1_overall_metabotype_tests.csv"),
          row.names = FALSE)

message("\nOverall metabotype effect (nested F-test):")
print(as.data.frame(tests_all %>%
  mutate(F_statistic = round(F_statistic, 2),
         p_metabotype = signif(p_metabotype, 3))))

message("\nContrast estimates (difference vs ", ref_level, ", 95% CI):")
print(as.data.frame(coefs_all %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ round(.x, 2)),
         p.value = signif(p.value, 3)) %>%
  select(outcome, model, contrast, n, estimate, conf.low, conf.high, p.value)))

# ================================================================
# 6. FOREST PLOT
# ================================================================
contrast_levels <- sort(unique(coefs_all$contrast))
pal <- setNames(
  c("#B2182B", "#999999", "#2166AC", "#1B9E77")[seq_along(contrast_levels)],
  contrast_levels)

primary_label <- names(samples)[1]
plot_df <- coefs_all %>%
  filter(sample == primary_label) %>%
  mutate(model = factor(model, levels = rev(names(model_specs))),
         outcome = factor(outcome,
                          levels = c("IFTA", "Interstitial Fibrosis",
                                     "Tubular Atrophy")))

p_forest <- ggplot(plot_df, aes(estimate, model, color = contrast)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~ outcome, scales = "free_x") +
  scale_color_manual(values = pal, name = "Contrast") +
  labs(title = paste0("Structural kidney injury by plasma metabotype, ",
                      "before and after adjustment (", primary_label, ")"),
       x = paste0("Difference vs ", ref_level, " (95% CI)"), y = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(out_dir, "R2C3_Figure1_forest.png"), p_forest,
       width = 11, height = 5, dpi = 300)
ggsave(file.path(out_dir, "R2C3_Figure1_forest.pdf"), p_forest,
       width = 11, height = 5)

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_R2C3_Figure1.txt"))

message("\n========== Done ==========")
message("Outputs in: ", out_dir)

