################################################################
# Figure 6: CKD plasma metabotypes & cell-type transcriptional programs
#
# 6A: PCA of plasma metabotypes (k = 3) in the scRNA-matched subcohort
# 6B: Clinical kidney function markers by metabotype
# 6C: Structural kidney injury by metabotype
# 6D: Transcriptional module scores by metabotype
# 6E: POD metabotype-defining metabolite-gene heatmap
# 6F: PT  metabotype-defining metabolite-gene heatmap
################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(pheatmap)
  library(RColorBrewer)
})

set.seed(1)

# ================================================================
# EDIT THIS ONE LINE to point to the folder containing the inputs
# ================================================================
base_dir <- "."   # e.g. "D:/UC Berkeley/KPMP paper/repo/Figure6"

out_dir <- file.path(base_dir, "Figure6_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Metabotype palette (consistent across panels) ------------
met_colors <- c("M1" = "#2166AC", "M2" = "#B2182B", "M3" = "#999999")

# ================================================================
# FIGURE 6A: PCA
# ================================================================
pca_scores <- read_csv(file.path(base_dir, "Figure6_PanelA_PCA_scores.csv"),
                       show_col_types = FALSE) %>%
  mutate(metabotype = factor(metabotype, levels = c("M1", "M2", "M3")))

pc1_var <- 13.37; pc2_var <- 6.98   # variance explained (from PCA object)

p_pca <- ggplot(pca_scores, aes(PC1, PC2, color = metabotype)) +
  stat_ellipse(aes(fill = metabotype), geom = "polygon",
               alpha = 0.12, color = NA, level = 0.95) +
  stat_ellipse(aes(color = metabotype), linewidth = 0.8, level = 0.95) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = met_colors) +
  scale_fill_manual(values = met_colors) +
  labs(title = "Plasma metabotypes (k = 3)",
       x = paste0("PC1 (", pc1_var, "% variance)"),
       y = paste0("PC2 (", pc2_var, "% variance)"),
       color = "metabotype", fill = "metabotype") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
        legend.position = "right", legend.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "Figure6A_PCA.png"), p_pca, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "Figure6A_PCA.pdf"), p_pca, width = 7, height = 5.5)
message("Saved Figure 6A")

# ================================================================
# FIGURE 6B: Clinical markers by metabotype
# ================================================================
clin_long <- read_csv(file.path(base_dir, "Figure6B_clinical_markers_plotting_table.csv"),
                      show_col_types = FALSE) %>%
  mutate(metabotype = factor(metabotype, levels = c("M1", "M2", "M3")),
         marker     = factor(marker, levels = c("eGFR", "creatinine", "cystatinC", "BUN")))

p_clin <- ggplot(clin_long, aes(metabotype, value, fill = metabotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.55) +
  geom_jitter(aes(color = metabotype), width = 0.15, height = 0, size = 2, alpha = 0.85) +
  facet_wrap(~ marker, scales = "free_y", ncol = 2, drop = FALSE) +
  scale_fill_manual(values = met_colors) +
  scale_color_manual(values = met_colors) +
  labs(title = "Clinical kidney function markers by plasma metabotype",
       x = "Plasma metabotype", y = NULL) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "none")

ggsave(file.path(out_dir, "Figure6B_clinical_markers.png"), p_clin, width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Figure6B_clinical_markers.pdf"), p_clin, width = 10, height = 7)
message("Saved Figure 6B")

# ================================================================
# FIGURE 6C: Structural injury by metabotype
# ================================================================
struct_long <- read_csv(file.path(base_dir, "Figure6C_structural_injury_plotting_table.csv"),
                        show_col_types = FALSE) %>%
  mutate(metabotype = factor(metabotype, levels = c("M1", "M2", "M3")),
         marker     = factor(marker, levels = c("IFTA", "Interstitial Fibrosis", "Tubular Atrophy")))

p_struct <- ggplot(struct_long, aes(metabotype, value, fill = metabotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.55) +
  geom_jitter(aes(color = metabotype), width = 0.15, height = 0, size = 2, alpha = 0.85) +
  facet_wrap(~ marker, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = met_colors) +
  scale_color_manual(values = met_colors) +
  labs(title = "Structural kidney injury by plasma metabotype",
       x = "Plasma metabotype", y = "Percent (%)") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "none")

ggsave(file.path(out_dir, "Figure6C_structural_injury.png"), p_struct, width = 12, height = 6, dpi = 300)
ggsave(file.path(out_dir, "Figure6C_structural_injury.pdf"), p_struct, width = 12, height = 6)
message("Saved Figure 6C")

# ================================================================
# FIGURE 6D: Transcriptional module scores by metabotype
# ================================================================
mol_long <- read_csv(file.path(base_dir, "Figure6D_molecular_programs_plotting_table.csv"),
                     show_col_types = FALSE) %>%
  mutate(metabotype = factor(metabotype, levels = c("M1", "M2", "M3")),
         module     = factor(module, levels = c("AhR_response", "Fibrosis_ECM",
                                                "Inflammation", "Oxidative_stress",
                                                "Podocyte_identity", "Tubular_injury")))

p_mol <- ggplot(mol_long, aes(metabotype, module_score, fill = metabotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.65) +
  geom_jitter(aes(color = metabotype), width = 0.15, height = 0, size = 1.6, alpha = 0.85) +
  facet_wrap(~ module, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = met_colors) +
  scale_color_manual(values = met_colors) +
  labs(title = "Transcriptional module scores by plasma metabotype",
       x = "Plasma metabotype", y = "Mean module score") +
  theme_bw(base_size = 16) +
  theme(strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        legend.title = element_blank(), legend.position = "right")

ggsave(file.path(out_dir, "Figure6D_module_scores.png"), p_mol, width = 12, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Figure6D_module_scores.pdf"), p_mol, width = 12, height = 7)
message("Saved Figure 6D")

# ================================================================
# FIGURE 6E & 6F: Metabotype-defining metabolite-gene heatmaps
# ================================================================
pod_annot   <- read.csv(file.path(base_dir, "POD_gene_met_assoc_nonDrug_FULL_annotated.csv"), stringsAsFactors = FALSE)
pt_annot    <- read.csv(file.path(base_dir, "PT_gene_met_assoc_nonDrug_FULL_annotated.csv"),  stringsAsFactors = FALSE)
nondrug_tbl <- read.csv(file.path(base_dir, "FINAL_top_metabolites_nonDrug_nonX_withClass.csv"), stringsAsFactors = FALSE)
gene_sets   <- read.csv(file.path(base_dir, "gene_sets_modules.csv"), stringsAsFactors = FALSE)

# Standardize module table -> module, gene
if (all(c("module", "gene") %in% names(gene_sets))) {
  gene_sets2 <- gene_sets %>% transmute(module = as.character(module), gene = as.character(gene))
} else if (all(c("Module", "Gene") %in% names(gene_sets))) {
  gene_sets2 <- gene_sets %>% transmute(module = as.character(Module), gene = as.character(Gene))
} else stop("gene_sets_modules.csv must contain columns module,gene (or Module,Gene)")

gene_modules <- split(gene_sets2$gene, gene_sets2$module)

# One row per metabotype x metabolite (best FDR)
nondrug_tbl <- nondrug_tbl %>%
  mutate(metabolite = as.character(metabolite),
         metabotype_high = factor(metabotype_high, levels = c("M1", "M2", "M3"))) %>%
  group_by(metabotype_high, metabolite) %>%
  slice_min(FDR, n = 1, with_ties = FALSE) %>%
  ungroup()

make_combined_heatmap <- function(assoc_annot, celltype, nondrug_tbl, gene_modules,
                                  top_genes_per_module = 5,
                                  cell_w = 45, cell_h = 18, font_row = 9, font_col = 10) {
  met_ids <- unique(as.character(nondrug_tbl$metabolite))
  modules <- names(gene_modules)

  # Top genes per module (lowest p)
  all_genes_df <- data.frame()
  for (m in modules) {
    dd <- assoc_annot %>%
      filter(gene %in% gene_modules[[m]], metabolite %in% met_ids) %>%
      group_by(gene, metabolite) %>% slice_min(p, n = 1, with_ties = FALSE) %>% ungroup()
    if (nrow(dd) == 0) next
    gr <- dd %>% filter(is.finite(p)) %>% group_by(gene) %>%
      summarize(best_p = min(p), .groups = "drop") %>%
      arrange(best_p) %>% slice_head(n = top_genes_per_module) %>% mutate(module = m)
    all_genes_df <- bind_rows(all_genes_df, gr)
  }
  if (nrow(all_genes_df) == 0) stop("No genes found for ", celltype)

  all_genes_df <- all_genes_df %>% group_by(gene) %>%
    slice_min(best_p, n = 1, with_ties = FALSE) %>% ungroup()
  all_genes_df$module <- factor(all_genes_df$module, levels = modules)
  all_genes_df <- all_genes_df %>% arrange(module, best_p)
  gene_keep <- all_genes_df$gene

  # Beta matrix (genes x metabolites)
  mat <- assoc_annot %>%
    filter(gene %in% gene_keep, metabolite %in% met_ids) %>%
    group_by(gene, metabolite) %>% slice_min(p, n = 1, with_ties = FALSE) %>% ungroup() %>%
    select(gene, metabolite, beta) %>%
    pivot_wider(names_from = metabolite, values_from = beta) %>%
    as.data.frame()
  rownames(mat) <- mat$gene; mat$gene <- NULL
  mat <- as.matrix(mat); mat[is.na(mat)] <- 0

  # Column order: metabotype then FDR
  col_order_tbl <- nondrug_tbl %>% arrange(metabotype_high, FDR) %>%
    distinct(metabolite, .keep_all = TRUE) %>% filter(metabolite %in% colnames(mat))
  mat <- mat[, as.character(col_order_tbl$metabolite), drop = FALSE]

  # Clean metabolite labels
  col_labels <- col_order_tbl$metabolite_name_clean
  names(col_labels) <- col_order_tbl$metabolite
  colnames(mat) <- unname(col_labels[colnames(mat)])

  # Keep row order; drop constant rows
  mat <- mat[gene_keep[gene_keep %in% rownames(mat)], , drop = FALSE]
  mat <- mat[apply(mat, 1, sd) != 0, , drop = FALSE]
  all_genes_df <- all_genes_df %>% filter(gene %in% rownames(mat))

  # Manual row-scaling (matches Figure 4C-D)
  mat_scaled <- t(scale(t(mat))); mat_scaled[is.na(mat_scaled)] <- 0
  lim    <- max(abs(mat_scaled), na.rm = TRUE)
  breaks <- seq(-lim, lim, length.out = 101)
  pal    <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

  # Annotations
  ann_row <- data.frame(
    Module = as.character(all_genes_df$module[match(rownames(mat_scaled), all_genes_df$gene)]),
    row.names = rownames(mat_scaled), check.names = FALSE)
  ann_col <- data.frame(
    Metabotype = as.character(col_order_tbl$metabotype_high),
    `Metabolite Class` = as.character(col_order_tbl$met_class),
    row.names = col_labels, check.names = FALSE)

  metabotype_colors <- met_colors
  class_levels <- sort(unique(ann_col$`Metabolite Class`))
  class_pal <- if (length(class_levels) <= 8)
    brewer.pal(max(3, length(class_levels)), "Set2")[seq_along(class_levels)]
  else colorRampPalette(brewer.pal(8, "Set2"))(length(class_levels))
  names(class_pal) <- class_levels

  module_levels <- levels(all_genes_df$module)
  module_levels <- module_levels[module_levels %in% unique(all_genes_df$module)]
  module_pal <- if (length(module_levels) <= 8)
    brewer.pal(max(3, length(module_levels)), "Dark2")[seq_along(module_levels)]
  else colorRampPalette(brewer.pal(8, "Dark2"))(length(module_levels))
  names(module_pal) <- module_levels

  annotation_colors <- list(Metabotype = metabotype_colors,
                            `Metabolite Class` = class_pal, Module = module_pal)

  # Gaps
  row_gaps <- cumsum(table(factor(ann_row$Module, levels = module_levels)))
  row_gaps <- row_gaps[-length(row_gaps)]
  col_gaps <- cumsum(table(factor(ann_col$Metabotype, levels = c("M1", "M2", "M3"))))[1:2]

  panel  <- ifelse(celltype == "POD", "E", "F")
  ttl    <- paste0(celltype, " \u2014 metabotype-defining metabolite\u2013gene associations")
  out_png <- file.path(out_dir, paste0("Figure6", panel, "_", celltype, "_heatmap.png"))
  out_pdf <- file.path(out_dir, paste0("Figure6", panel, "_", celltype, "_heatmap.pdf"))

  for (f in c(out_png, out_pdf)) {
    pheatmap(mat_scaled, scale = "none", color = pal, breaks = breaks,
             cluster_rows = FALSE, cluster_cols = FALSE,
             annotation_col = ann_col, annotation_row = ann_row,
             annotation_colors = annotation_colors,
             fontsize_row = font_row, fontsize_col = font_col, angle_col = 45,
             border_color = "grey85", main = ttl,
             gaps_row = row_gaps, gaps_col = col_gaps,
             cellwidth = cell_w, cellheight = cell_h, filename = f)
  }
  message("Saved Figure 6", panel, " (", celltype, "): ",
          nrow(mat_scaled), " genes x ", ncol(mat_scaled), " metabolites")
}

message("\n========== POD heatmap (6E) ==========")
make_combined_heatmap(pod_annot, "POD", nondrug_tbl, gene_modules, top_genes_per_module = 5)

message("\n========== PT heatmap (6F) ==========")
make_combined_heatmap(pt_annot, "PT", nondrug_tbl, gene_modules, top_genes_per_module = 5)

message("\n========== Figure 6 complete ==========")
message("Outputs in: ", out_dir)
