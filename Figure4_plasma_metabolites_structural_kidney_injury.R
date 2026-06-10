################################################################
# Figure 4: Plasma metabolites & structural kidney injury
#
# 4A: Volcano  - interstitial fibrosis (univariate LM, colored by origin)
# 4B: Volcano  - tubular atrophy (univariate LM, colored by origin)
# 4C: Heatmap  - PT pseudobulk gene expression vs top-20 fibrosis RF metabolites
# 4D: Heatmap  - PT pseudobulk gene expression vs top-20 atrophy RF metabolites
#
# Volcano inputs are univariate linear-regression results (beta, FDR).
# Heatmap inputs are participant-level metabolite-gene beta matrices.
################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

# ================================================================
# EDIT THIS ONE LINE to point to the folder containing the inputs
# ================================================================
base_dir <- "."   # e.g. "D:/UC Berkeley/KPMP paper/repo/Figure4"

out_dir <- file.path(base_dir, "Figure4_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Input files (place all in base_dir) ----------------------
fib_lm_file  <- file.path(base_dir, "Q1_Plasma_LM_Fibrosis_full_results.csv")
ta_lm_file   <- file.path(base_dir, "Plasma_TubularAtrophy_linear_regression_results.csv")
fib_mat_file <- file.path(base_dir, "PT_fibrosis_heatmap_matrix.csv")
ta_mat_file  <- file.path(base_dir, "PT_tubular_atrophy_heatmap_matrix.csv")

# ---- Metabolite origin palette --------------------------------
origin_colors <- c(
  "Host-derived"         = "#7B3294",  # purple
  "Diet-derived"         = "#E66101",  # orange
  "Microbiome-dependent" = "#1B9E77",  # teal
  "Xenobiotic"           = "#999999",  # gray
  "Not significant"      = "grey80"
)

# ================================================================
# FIGURE 4A & 4B: VOLCANO PLOTS (colored by metabolite origin)
# ================================================================
make_volcano <- function(path, title_text, outfile, n_labels = 10) {
  res <- read_csv(path, show_col_types = FALSE)

  # Expect: estimate, FDR, chemical_name; metabolite_origin optional
  origin_col <- intersect(c("metabolite_origin", "origin"), names(res))
  res <- res %>%
    mutate(
      estimate      = as.numeric(estimate),
      FDR           = as.numeric(FDR),
      neg_log10_FDR = -log10(FDR)
    )

  if (length(origin_col) == 1) {
    res$plot_group <- ifelse(res$FDR < 0.05, res[[origin_col[1]]], "Not significant")
  } else {
    # Fallback: significance only (origin annotation can be merged in later)
    res$plot_group <- ifelse(res$FDR < 0.05, "Host-derived", "Not significant")
    message("NOTE: no metabolite_origin column in ", basename(path),
            " - points shown as significant/not. Merge origin labels for final figure.")
  }

  labels_df <- res %>%
    filter(FDR < 0.05, !grepl("^X-", chemical_name)) %>%
    arrange(desc(neg_log10_FDR)) %>%
    slice_head(n = n_labels)

  p <- ggplot(res, aes(estimate, neg_log10_FDR, color = plot_group)) +
    geom_point(alpha = 0.75, size = 1.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = labels_df, aes(label = chemical_name),
                    color = "black", size = 3.2, fontface = "bold",
                    box.padding = 0.4, max.overlaps = Inf,
                    segment.color = "grey50", segment.size = 0.3) +
    scale_color_manual(values = origin_colors, name = "Metabolite origin") +
    labs(title = title_text, x = "Effect size (\u03B2)",
         y = expression(-log[10](FDR))) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")

  ggsave(outfile, p, width = 9, height = 6, dpi = 300)
  n_sig <- sum(res$FDR < 0.05, na.rm = TRUE)
  message("Saved: ", basename(outfile), "  (", n_sig, " significant at FDR<0.05)")
}

make_volcano(fib_lm_file,
             "Plasma metabolites associated with interstitial fibrosis",
             file.path(out_dir, "Figure4A_volcano_fibrosis.png"))

make_volcano(ta_lm_file,
             "Plasma metabolites associated with tubular atrophy",
             file.path(out_dir, "Figure4B_volcano_atrophy.png"))

# ================================================================
# FIGURE 4C & 4D: PROXIMAL TUBULE HEATMAPS
# Row-scaled beta matrices; blue-white-red palette
# ================================================================
read_heatmap_matrix <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  rn <- df[[1]]
  mat <- as.matrix(df[, -1, drop = FALSE])
  mat <- apply(mat, 2, as.numeric)
  rownames(mat) <- rn
  message("Loaded matrix: ", nrow(mat), " genes x ", ncol(mat),
          " metabolites from ", basename(path))
  mat
}

save_heatmap <- function(mat, title_text, out_png, out_pdf) {
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[is.na(mat_scaled)] <- 0
  lim    <- max(abs(mat_scaled), na.rm = TRUE)
  breaks <- seq(-lim, lim, length.out = 101)
  pal    <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

  for (f in c(out_png, out_pdf)) {
    pheatmap(mat_scaled, color = pal, breaks = breaks,
             cluster_rows = TRUE, cluster_cols = TRUE,
             border_color = "grey85", main = title_text,
             fontsize_row = 8, fontsize_col = 9, angle_col = 45,
             filename = f, width = 8, height = 8)
  }
  message("Saved: ", basename(out_png), " + .pdf")
}

save_heatmap(read_heatmap_matrix(fib_mat_file),
             "PT \u2014 fibrosis RF metabolites",
             file.path(out_dir, "Figure4C_PT_fibrosis_heatmap.png"),
             file.path(out_dir, "Figure4C_PT_fibrosis_heatmap.pdf"))

save_heatmap(read_heatmap_matrix(ta_mat_file),
             "PT \u2014 tubular atrophy RF metabolites",
             file.path(out_dir, "Figure4D_PT_tubular_atrophy_heatmap.png"),
             file.path(out_dir, "Figure4D_PT_tubular_atrophy_heatmap.pdf"))

message("\n========== Figure 4 complete ==========")
message("Outputs in: ", out_dir)
