################################################################
# Figure 2: Metabotype score, eGFR associations, and pathway enrichment
#
# 2A: Association between metabotype score and baseline eGFR
# 2B: Top plasma metabolites associated with eGFR, colored by origin
# 2C: Pathway enrichment of metabolites ranked by correlation
#     with metabotype score
#
################################################################

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
