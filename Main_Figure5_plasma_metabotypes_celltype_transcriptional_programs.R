################################################################
# Figure 5 (revised): panels A-D only, Subgroup A/B/C labels
#
# 5A: PCA of plasma subgroups (k = 3) in the scRNA-matched subcohort
# 5B: Clinical kidney function markers by subgroup
# 5C: Structural kidney injury by subgroup
# 5D: Transcriptional module scores by subgroup
#
#
# Labels: the stored data files keep the original M1-M3 labels.
# SUBGROUP_LABELS below is the single point of translation:
#   M1 = Subgroup A = blue (#2166AC)
#   M2 = Subgroup B = red  (#B2182B)
#   M3 = Subgroup C = grey (#999999)
#
# INPUTS (data/ or scRNA_01_output/)
#   Figure5_PanelA_PCA_scores.csv, Figure5B_clinical_markers_plotting_table.csv,
#   Figure5C_structural_injury_plotting_table.csv,
#   Figure5D_molecular_programs_plotting_table.csv
################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

set.seed(1)

# ================================================================
# EDIT THIS ONE LINE to point to the folder that contains data/
# ================================================================
base_dir <- "."   # e.g. "D:/UC Berkeley/KPMP paper/repo"

data_dir <- file.path(base_dir, "data")
out_dir  <- file.path(base_dir, "Figure6_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TIFF_DPI <- 600   # combined figure TIFF; LZW-compressed

# ---- Locate inputs ---------------------------------------------
# Plotting tables are read from data/ (Dryad) 
find_input <- function(fname) {
  cands <- c(file.path(data_dir, fname), file.path(base_dir, "scRNA_01_output", fname))
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0) stop("Could not find ", fname, ". Looked in:\n  ", paste(cands, collapse = "\n  "))
  message("  ", fname, " -> ", hit[1])
  hit[1]
}

message("Locating input files...")
f_pca    <- find_input("Figure5_PanelA_PCA_scores.csv")
f_clin   <- find_input("Figure5B_clinical_markers_plotting_table.csv")
f_struct <- find_input("Figure5C_structural_injury_plotting_table.csv")
f_mol    <- find_input("Figure5D_molecular_programs_plotting_table.csv")

# ---- Label mapping (stored label -> displayed label) ----------
SUBGROUP_LABELS <- c("M1" = "Subgroup A",
                     "M2" = "Subgroup B",
                     "M3" = "Subgroup C")
subgroup_levels <- unname(SUBGROUP_LABELS)

# Short tick labels for the x-axis (axis title says "Plasma subgroup")
TICK_LABELS <- setNames(sub("^Subgroup ", "", subgroup_levels), subgroup_levels)

relabel <- function(x) {
  out <- unname(SUBGROUP_LABELS[as.character(x)])
  bad <- unique(as.character(x)[is.na(out) & !is.na(x)])
  if (length(bad) > 0) {
    stop("Unexpected metabotype labels in input: ", paste(bad, collapse = ", "),
         "\n  Expected: ", paste(names(SUBGROUP_LABELS), collapse = ", "))
  }
  factor(out, levels = subgroup_levels)
}

# Display names for facets; values not listed are kept as they are
# (with a message), so a renamed marker never silently disappears.
pretty_factor <- function(x, map) {
  x <- as.character(x)
  extra <- setdiff(unique(x), names(map))
  if (length(extra) > 0) message("  NOTE: no display name for: ", paste(extra, collapse = ", "))
  out <- ifelse(x %in% names(map), unname(map[x]), x)
  factor(out, levels = c(unname(map), extra))
}

clin_names <- c("eGFR" = "eGFR", "creatinine" = "Creatinine",
                "cystatinC" = "Cystatin C", "BUN" = "BUN")
struct_names <- c("IFTA" = "IFTA",
                  "Interstitial Fibrosis" = "Interstitial fibrosis",
                  "Tubular Atrophy" = "Tubular atrophy")
module_names <- c("AhR_response" = "AhR response", "Fibrosis_ECM" = "Fibrosis/ECM",
                  "Inflammation" = "Inflammation", "Oxidative_stress" = "Oxidative stress",
                  "Podocyte_identity" = "Podocyte identity", "Tubular_injury" = "Tubular injury")

# ---- Subgroup palette (identical to manuscript version) -------
met_colors <- setNames(c("#2166AC", "#B2182B", "#999999"), subgroup_levels)

# Shared box-and-points layer for panels B-D
box_layers <- function(box_width, point_size) {
  list(
    geom_boxplot(outlier.shape = NA, width = box_width, alpha = 0.55),
    geom_jitter(aes(color = subgroup), width = 0.15, height = 0,
                size = point_size, alpha = 0.85),
    scale_fill_manual(values = met_colors),
    scale_color_manual(values = met_colors),
    scale_x_discrete(labels = TICK_LABELS)
  )
}

# ================================================================
# FIGURE 5A: PCA
# ================================================================
pca_scores <- read_csv(f_pca, show_col_types = FALSE) %>%
  mutate(subgroup = relabel(metabotype))

# variance explained: from scRNA_01 output if available
pc1_var <- 13.37; pc2_var <- 6.98
f_var <- file.path(dirname(f_pca), "Figure5A_PCA_variance_explained.csv")
if (file.exists(f_var)) {
  pv <- read_csv(f_var, show_col_types = FALSE)
  pc1_var <- pv$percent_variance[pv$PC == "PC1"]; pc2_var <- pv$percent_variance[pv$PC == "PC2"]
}

p_pca <- ggplot(pca_scores, aes(PC1, PC2, color = subgroup)) +
  stat_ellipse(aes(fill = subgroup), geom = "polygon",
               alpha = 0.12, color = NA, level = 0.95) +
  stat_ellipse(aes(color = subgroup), linewidth = 0.8, level = 0.95) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = met_colors) +
  scale_fill_manual(values = met_colors) +
  labs(title = "Plasma subgroups (k = 3)",
       x = paste0("PC1 (", pc1_var, "% variance)"),
       y = paste0("PC2 (", pc2_var, "% variance)"),
       color = "Subgroup", fill = "Subgroup") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
        legend.position = "right", legend.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "Figure5A_PCA.png"), p_pca, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "Figure5A_PCA.pdf"), p_pca, width = 7, height = 5.5)
message("Saved Figure 5A")

# ================================================================
# FIGURE 5B: Clinical markers by subgroup
# ================================================================
clin_long <- read_csv(f_clin, show_col_types = FALSE) %>%
  mutate(subgroup = relabel(metabotype),
         marker   = pretty_factor(marker, clin_names))

p_clin <- ggplot(clin_long, aes(subgroup, value, fill = subgroup)) +
  box_layers(box_width = 0.6, point_size = 2) +
  facet_wrap(~ marker, scales = "free_y", ncol = 2, drop = FALSE) +
  labs(title = "Clinical kidney function markers by plasma subgroup",
       x = "Plasma subgroup", y = NULL) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "none")

ggsave(file.path(out_dir, "Figure5B_clinical_markers.png"), p_clin, width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Figure5B_clinical_markers.pdf"), p_clin, width = 10, height = 7)
message("Saved Figure 5B")

# ================================================================
# FIGURE 5C: Structural injury by subgroup
# ================================================================
struct_long <- read_csv(f_struct, show_col_types = FALSE) %>%
  filter(!is.na(value)) %>%   # one participant has no pathology scores
  mutate(subgroup = relabel(metabotype),
         marker   = pretty_factor(marker, struct_names))

p_struct <- ggplot(struct_long, aes(subgroup, value, fill = subgroup)) +
  box_layers(box_width = 0.65, point_size = 2) +
  facet_wrap(~ marker, scales = "free_y", ncol = 3) +
  labs(title = "Structural kidney injury by plasma subgroup",
       x = "Plasma subgroup", y = "Percent (%)") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "none")

ggsave(file.path(out_dir, "Figure5C_structural_injury.png"), p_struct, width = 12, height = 6, dpi = 300)
ggsave(file.path(out_dir, "Figure5C_structural_injury.pdf"), p_struct, width = 12, height = 6)
message("Saved Figure 5C")

# ================================================================
# FIGURE 5D: Transcriptional module scores by subgroup
# ================================================================
mol_long <- read_csv(f_mol, show_col_types = FALSE) %>%
  mutate(subgroup = relabel(metabotype),
         module   = pretty_factor(module, module_names))

p_mol <- ggplot(mol_long, aes(subgroup, module_score, fill = subgroup)) +
  box_layers(box_width = 0.65, point_size = 1.6) +
  facet_wrap(~ module, scales = "free_y", ncol = 3) +
  labs(title = "Transcriptional module scores by plasma subgroup",
       x = "Plasma subgroup", y = "Mean module score") +
  theme_bw(base_size = 16) +
  theme(strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        legend.title = element_blank(), legend.position = "right")

ggsave(file.path(out_dir, "Figure5D_module_scores.png"), p_mol, width = 12, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Figure5D_module_scores.pdf"), p_mol, width = 12, height = 7)
message("Saved Figure 5D")

