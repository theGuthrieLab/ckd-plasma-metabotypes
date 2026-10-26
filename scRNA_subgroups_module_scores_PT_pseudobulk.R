#!/usr/bin/env Rscript

################################################################
# scRNA-seq step 1: plasma subgroups in the scRNA-seq subcohort and
#                   participant-level transcriptional module scores
#
#
# LABELS
# Clusters are stored as M1-M3. In the revised manuscript they are
# displayed as Subgroups A-C (M1 = A, M2 = B, M3 = C); the figure
# scripts apply that mapping at plotting time.
#
# INPUTS (data/)
#   scRNA/*.rds                        KPMP cell-type subset Seurat objects
#                                      (not redistributed; see README)
#   Plasma_Metabolomics.xlsx           long format: Participant ID,
#                                      Metabolite, Abundance, ...
#   curated_clinical.csv               (or "curated clinical.csv")
#   Serum-biomarkers-KPMP-CKD_Healthy_with_Baseline_eGFR.csv
#                                      creatinine, cystatin C, BUN
#   KPMP TIV Descriptor Scores 03.30.2025.xlsx
#                                      KPMP histopathology descriptors
#   gene_sets_modules.csv              columns: module, gene
#
################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(broom)
})

set.seed(1)

#############################################
# 0) Paths + parameters
#############################################

# ================================================================
# EDIT THIS ONE LINE to point to the folder that contains data/
# ================================================================
base_dir <- "."   # e.g. "D:/UC Berkeley/KPMP paper/repo"

data_dir <- file.path(base_dir, "data")
out_dir  <- file.path(base_dir, "scRNA_01_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) stop("Input not found. Looked for:\n  ", paste(paths, collapse = "\n  "))
  hit[1]
}

sc_folder    <- file.path(data_dir, "scRNA")
plasma_xlsx  <- first_existing(file.path(data_dir, "Plasma_Metabolomics.xlsx"))
serum_csv    <- first_existing(file.path(data_dir, "Serum-biomarkers-KPMP-CKD_Healthy_with_Baseline_eGFR.csv"))
clinical_csv <- first_existing(file.path(data_dir, c("curated_clinical.csv", "curated clinical.csv")))
tiv_xlsx     <- first_existing(file.path(data_dir, "KPMP TIV Descriptor Scores 03.30.2025.xlsx"))
genesets_csv <- first_existing(file.path(data_dir, "gene_sets_modules.csv"))

# --- Clustering settings ---
k <- 3                                  # number of clusters
nstart <- 50                            # k-means restarts
min_nonmissing_per_metabolite <- 0.70   # keep metabolites measured in >= 70% of participants

#############################################
# 1) Load scRNA Seurat subset objects (.rds) + inventory metadata
#############################################

sc_files <- list.files(sc_folder, pattern = "\\.rds$", full.names = TRUE)
if (length(sc_files) == 0) stop("No .rds files found in ", sc_folder)
message("scRNA-seq objects found: ", length(sc_files))

is_seurat <- function(x) inherits(x, "Seurat")

get_meta_cols <- function(obj){
  if(!is_seurat(obj)) return(NA_character_)
  paste(colnames(obj@meta.data), collapse = ";")
}

get_ids_from_file <- function(f){
  obj <- readRDS(f)
  if(!is_seurat(obj)) return(NULL)

  md <- obj@meta.data

  # SampleID in metadata; fallbacks in order of preference
  if("SampleID" %in% colnames(md)) sample_id <- md$SampleID
  else if("ParticipantID" %in% colnames(md)) sample_id <- md$ParticipantID
  else if("Participant.ID" %in% colnames(md)) sample_id <- md$Participant.ID
  else if("KPMPID" %in% colnames(md)) sample_id <- md$KPMPID
  else if("SpecimenID" %in% colnames(md)) sample_id <- md$SpecimenID
  else if("InternalID" %in% colnames(md)) sample_id <- md$InternalID
  else {
    # last resort: prefix before the first "_" of the cell name
    cn <- colnames(obj)
    sample_id <- str_extract(cn, "^[^_]+")
  }

  # cell population label: metadata if present, otherwise the file name
  cell_pop <- if("cell_pop" %in% colnames(md)) as.character(md$cell_pop) else basename(f) |> str_replace("\\.rds$", "")

  tibble(
    SampleID = as.character(sample_id),
    cell_pop = cell_pop
  )
}

meta_inventory <- tibble(
  file = basename(sc_files),
  is_seurat = NA,
  n_cells = NA_integer_,
  n_genes = NA_integer_,
  meta_cols = NA_character_
)

for(i in seq_along(sc_files)){
  obj <- readRDS(sc_files[i])
  meta_inventory$is_seurat[i] <- is_seurat(obj)
  if(is_seurat(obj)){
    meta_inventory$n_cells[i] <- ncol(obj)
    meta_inventory$n_genes[i] <- nrow(obj)
    meta_inventory$meta_cols[i] <- get_meta_cols(obj)
  }
  rm(obj); invisible(gc())
}

write.csv(meta_inventory, file.path(out_dir, "scRNA_meta_inventory.csv"), row.names = FALSE)

# All participant IDs present in any scRNA-seq object
all_sc_sample_ids <- bind_rows(lapply(sc_files, get_ids_from_file)) %>%
  filter(!is.na(SampleID)) %>%
  distinct(SampleID) %>%
  pull(SampleID)

message("Participants with scRNA-seq data: ", length(all_sc_sample_ids))

#############################################
# 2) Plasma metabolomics (long) -> subgroups (k-means)
#############################################

plasma <- read_excel(plasma_xlsx)

stopifnot(all(c("Participant ID", "Metabolite", "Abundance") %in% colnames(plasma)))

plasma2 <- plasma %>%
  mutate(
    ParticipantID = as.character(`Participant ID`),
    Metabolite = as.character(Metabolite),
    Abundance_num = suppressWarnings(as.numeric(Abundance))
  ) %>%
  select(ParticipantID, Metabolite, Abundance_num)

# restrict to participants with scRNA-seq data
plasma2 <- plasma2 %>% filter(ParticipantID %in% all_sc_sample_ids)

# participant x metabolite matrix (median of duplicate measurements)
plasma_wide <- plasma2 %>%
  group_by(ParticipantID, Metabolite) %>%
  summarize(Abundance = median(Abundance_num, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Metabolite, values_from = Abundance)

# filter metabolites by missingness
miss_rate <- plasma_wide %>%
  summarize(across(-ParticipantID, ~ mean(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "Metabolite", values_to = "missing_rate")

keep_mets <- miss_rate %>%
  filter(missing_rate <= (1 - min_nonmissing_per_metabolite)) %>%
  pull(Metabolite)

plasma_wide_f <- plasma_wide %>%
  select(ParticipantID, all_of(keep_mets))

mat <- plasma_wide_f %>% select(-ParticipantID) %>% as.data.frame()
mat <- as.matrix(mat)

# drop all-NA and zero-variance metabolites
mat <- mat[, colSums(!is.na(mat)) > 0, drop = FALSE]
v <- apply(mat, 2, function(x) var(x, na.rm = TRUE))
mat <- mat[, is.finite(v) & v > 0, drop = FALSE]

# impute missing values with the metabolite median
col_meds <- apply(mat, 2, function(x) median(x, na.rm = TRUE))
for(j in seq_len(ncol(mat))){
  idx <- which(is.na(mat[, j]))
  if(length(idx) > 0) mat[idx, j] <- col_meds[j]
}

mat_scaled <- scale(mat)
stopifnot(sum(!is.finite(mat_scaled)) == 0)

message("Clustering ", nrow(mat_scaled), " participants on ", ncol(mat_scaled), " metabolites")

km <- kmeans(mat_scaled, centers = k, nstart = nstart)

metabotype_tbl <- tibble(
  ParticipantID = plasma_wide_f$ParticipantID,
  metabotype = paste0("M", km$cluster)
)

write.csv(metabotype_tbl, file.path(out_dir, "plasma_metabotypes_kmeans.csv"), row.names = FALSE)
print(table(metabotype_tbl$metabotype))

# PCA of the same scaled matrix (Figure 5A)
pca <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_var <- round(100 * var_explained[1], 2)
pc2_var <- round(100 * var_explained[2], 2)

pca_df <- tibble(
  ParticipantID = plasma_wide_f$ParticipantID,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
) %>% left_join(metabotype_tbl, by = "ParticipantID")

write.csv(pca_df, file.path(out_dir, "Figure5_PanelA_PCA_scores.csv"), row.names = FALSE)
write.csv(tibble(PC = c("PC1", "PC2"), percent_variance = c(pc1_var, pc2_var)),
          file.path(out_dir, "Figure5A_PCA_variance_explained.csv"), row.names = FALSE)
message("PC1 = ", pc1_var, "% variance; PC2 = ", pc2_var, "% variance")

#############################################
# 3) Kidney function markers by subgroup (Figure 5B)
#############################################

clinical <- read.csv(clinical_csv, stringsAsFactors = FALSE)

if("Participant.ID" %in% colnames(clinical)) clinical <- clinical %>% rename(ParticipantID = Participant.ID)
if("Participant ID" %in% colnames(clinical)) clinical <- clinical %>% rename(ParticipantID = `Participant ID`)
if(!"ParticipantID" %in% colnames(clinical)) stop("clinical file missing ParticipantID column")

clinical2 <- clinical %>% select(ParticipantID, Enrollment.Category, Sex, Age..Years...Binned., Baseline.eGFR..ml.min.1.73m2.) %>%
  rename(
    Sex = Sex,
    Age = Age..Years...Binned.,
    eGFR = Baseline.eGFR..ml.min.1.73m2.
  )

metabotype_clin <- metabotype_tbl %>%
  left_join(clinical2, by = "ParticipantID")

serum <- read.csv(serum_csv, stringsAsFactors = FALSE)

serum_sub <- serum %>%
  rename(ParticipantID = Participant.ID) %>%
  transmute(
    ParticipantID = as.character(ParticipantID),
    creatinine = suppressWarnings(as.numeric(`Creatinine..mg.dL.`)),
    cystatinC  = suppressWarnings(as.numeric(`Cystatin.C..mg.L.`)),
    BUN        = suppressWarnings(as.numeric(`Blood.urea.nitrogen..mg.dL.`))
  )

metabotype_clin2 <- metabotype_clin %>%
  left_join(serum_sub, by = "ParticipantID")

# Unadjusted and age- and sex-adjusted models; M1 is the reference
fit_marker_models <- function(df, marker){
  df2 <- df %>% filter(!is.na(.data[[marker]])) %>%
    mutate(
      metabotype = factor(metabotype, levels = c("M1","M2","M3")),
      Sex = factor(Sex),
      Age = factor(Age)
    )

  m_unadj <- lm(reformulate("metabotype", response = marker), data = df2)
  m_adj   <- lm(reformulate(c("metabotype","Age","Sex"), response = marker), data = df2)

  tibble(
    marker = marker,
    model = c("unadjusted", "adjusted_AgeSex"),
    n = nrow(df2),
    r2 = c(summary(m_unadj)$r.squared, summary(m_adj)$r.squared),
    adj_r2 = c(summary(m_unadj)$adj.r.squared, summary(m_adj)$adj.r.squared),
    f_p = c(anova(m_unadj)[["Pr(>F)"]][1], anova(m_adj)[["Pr(>F)"]][1])
  ) %>%
    left_join(
      bind_rows(
        broom::tidy(m_unadj) %>% mutate(model = "unadjusted"),
        broom::tidy(m_adj) %>% mutate(model = "adjusted_AgeSex")
      ) %>% select(model, term, estimate, std.error, statistic, p.value),
      by = "model"
    )
}

markers_clin <- c("eGFR","creatinine","cystatinC","BUN")
clin_tbl <- bind_rows(lapply(markers_clin, \(m) fit_marker_models(metabotype_clin2, m)))
write.csv(clin_tbl, file.path(out_dir, "Q5_clinical_marker_regression_tables.csv"), row.names = FALSE)

clin_long <- metabotype_clin2 %>%
  select(ParticipantID, metabotype, all_of(markers_clin)) %>%
  pivot_longer(cols = all_of(markers_clin), names_to = "marker", values_to = "value") %>%
  mutate(metabotype = factor(metabotype, levels = c("M1","M2","M3")))

write.csv(clin_long, file.path(out_dir, "Figure5B_clinical_markers_plotting_table.csv"), row.names = FALSE)

#############################################
# 4) Histopathology by subgroup (Figure 5C)
#############################################

# The KPMP descriptor workbook has title rows above the header; find the
# row that holds "Participant ID" and use it as the header.
tiv_raw <- suppressMessages(read_excel(tiv_xlsx, col_names = FALSE))
hdr_row <- which(apply(tiv_raw, 1, function(r) any(trimws(as.character(r)) == "Participant ID")))[1]
if (is.na(hdr_row)) stop("Could not find a 'Participant ID' header row in ", basename(tiv_xlsx))
tiv <- read_excel(tiv_xlsx, skip = hdr_row - 1)

names(tiv) <- names(tiv) %>%
  str_replace_all("\\s+", " ") %>%
  str_trim()

tiv <- tiv %>%
  rename(
    ParticipantID = `Participant ID`,
    cortex_interstitial_fibrosis_percent = `Cortex: Interstitial Fibrosis: Percent`,
    cortex_tubular_atrophy_common_type_percent = `Cortex: Tubular Atrophy: Common Type: Percent`
  ) %>%
  mutate(
    ParticipantID = as.character(ParticipantID),
    cortex_interstitial_fibrosis_percent = suppressWarnings(as.numeric(cortex_interstitial_fibrosis_percent)),
    cortex_tubular_atrophy_common_type_percent = suppressWarnings(as.numeric(cortex_tubular_atrophy_common_type_percent))
  ) %>%
  mutate(
    IFTA_percent = cortex_interstitial_fibrosis_percent + cortex_tubular_atrophy_common_type_percent
  )

metabotype_struct <- metabotype_clin %>%
  left_join(tiv %>% select(ParticipantID, IFTA_percent, cortex_interstitial_fibrosis_percent, cortex_tubular_atrophy_common_type_percent),
            by = "ParticipantID") %>%
  mutate(metabotype = factor(metabotype, levels = c("M1","M2","M3")))

struct_markers <- c("IFTA_percent","cortex_interstitial_fibrosis_percent","cortex_tubular_atrophy_common_type_percent")
struct_tbl <- bind_rows(lapply(struct_markers, \(m) fit_marker_models(metabotype_struct, m)))
write.csv(struct_tbl, file.path(out_dir, "Q5_structural_regression_tables.csv"), row.names = FALSE)

struct_long <- metabotype_struct %>%
  select(ParticipantID, metabotype, all_of(struct_markers)) %>%
  pivot_longer(cols = all_of(struct_markers), names_to = "marker", values_to = "value") %>%
  mutate(
    marker = recode(marker,
                    IFTA_percent = "IFTA",
                    cortex_interstitial_fibrosis_percent = "Interstitial Fibrosis",
                    cortex_tubular_atrophy_common_type_percent = "Tubular Atrophy")
  )

write.csv(struct_long, file.path(out_dir, "Figure5C_structural_injury_plotting_table.csv"), row.names = FALSE)

#############################################
# 5) Participant-level module scores (Figure 5D)
#############################################

gene_sets <- read.csv(genesets_csv, stringsAsFactors = FALSE)
stopifnot(all(c("module","gene") %in% colnames(gene_sets)))

modules <- unique(gene_sets$module)
gene_list <- lapply(modules, function(m) unique(gene_sets$gene[gene_sets$module == m]))
names(gene_list) <- modules

# Score one cell-type object: AddModuleScore per cell, then the mean
# per participant. Modules with < 5 genes present in the object are skipped.
score_one_obj <- function(f){
  obj <- readRDS(f)
  if(!is_seurat(obj)) return(NULL)

  md <- obj@meta.data
  cell_pop <- basename(f) |> str_replace("\\.rds$", "")

  if("SampleID" %in% colnames(md)) md$SampleID <- as.character(md$SampleID)
  else if("ParticipantID" %in% colnames(md)) md$SampleID <- as.character(md$ParticipantID)
  else md$SampleID <- str_extract(colnames(obj), "^[^_]+")

  feats <- rownames(obj)
  gl <- lapply(gene_list, \(g) intersect(g, feats))

  keep <- sapply(gl, length) >= 5
  gl <- gl[keep]
  if(length(gl) == 0) return(NULL)

  obj <- AddModuleScore(obj, features = gl, name = "module_", assay = DefaultAssay(obj))

  # AddModuleScore adds module_1, module_2, ... in the order of gl
  score_cols <- grep("^module_\\d+$", colnames(obj@meta.data), value = TRUE)

  module_map <- tibble(
    module = names(gl),
    score_col = score_cols
  )

  out <- obj@meta.data %>%
    as_tibble(rownames = "cell") %>%
    select(cell, SampleID, any_of(c("cell_type","celltype","cell_type_ontology_term_id","cluster"))) %>%
    mutate(cell_pop = cell_pop) %>%
    bind_cols(obj@meta.data %>% as_tibble() %>% select(all_of(score_cols)))

  out2 <- out %>%
    group_by(SampleID, cell_pop) %>%
    summarize(across(all_of(score_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(cols = all_of(score_cols), names_to = "score_col", values_to = "module_score") %>%
    left_join(module_map, by = "score_col") %>%
    select(SampleID, cell_pop, module, module_score)

  rm(obj); invisible(gc())
  out2
}

sc_scores_all <- bind_rows(lapply(sc_files, score_one_obj)) %>%
  filter(!is.na(module_score))

saveRDS(sc_scores_all, file.path(out_dir, "scRNA_donor_level_module_scores_all.rds"))

clinical_cov <- clinical2 %>% select(ParticipantID, Age, Sex, eGFR)

sc_met <- sc_scores_all %>%
  left_join(metabotype_tbl, by = c("SampleID" = "ParticipantID")) %>%
  left_join(clinical_cov, by = c("SampleID" = "ParticipantID")) %>%
  mutate(
    metabotype = factor(metabotype, levels = c("M1","M2","M3")),
    Sex = factor(Sex),
    Age = factor(Age)
  )

# One value per participant and module: mean across cell types
sc_plot_df <- sc_met %>%
  group_by(SampleID, metabotype, module) %>%
  summarize(module_score = mean(module_score, na.rm = TRUE), .groups = "drop")

write.csv(sc_plot_df, file.path(out_dir, "Figure5D_molecular_programs_plotting_table.csv"), row.names = FALSE)

fit_sc_module_models <- function(df, mod){
  df2 <- df %>% filter(module == mod, !is.na(module_score)) %>%
    mutate(metabotype = factor(metabotype, levels = c("M1","M2","M3")))

  m_unadj <- lm(module_score ~ metabotype, data = df2)
  m_adj   <- lm(module_score ~ metabotype + Age + Sex + eGFR, data = df2)

  bind_rows(
    broom::tidy(m_unadj) %>% mutate(module = mod, model = "unadjusted"),
    broom::tidy(m_adj) %>% mutate(module = mod, model = "adjusted_AgeSexeGFR")
  ) %>%
    mutate(n = nrow(df2))
}

# Covariates joined for the regression models
sc_model_df <- sc_plot_df %>%
  left_join(clinical_cov, by = c("SampleID" = "ParticipantID")) %>%
  mutate(Sex = factor(Sex), Age = factor(Age))

molecular_tbl <- bind_rows(lapply(unique(sc_model_df$module), \(m) fit_sc_module_models(sc_model_df, m)))
write.csv(molecular_tbl, file.path(out_dir, "Q5_molecular_program_regression_tables.csv"), row.names = FALSE)

#############################################
# 6) Combined regression tables
#############################################

combined_tbl <- bind_rows(
  clin_tbl %>% mutate(section = "clinical"),
  struct_tbl %>% mutate(section = "structural"),
  molecular_tbl %>% rename(marker = module) %>% mutate(section = "molecular")
)

write.csv(combined_tbl, file.path(out_dir, "Q5_all_regression_tables_combined.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_scRNA_01.txt"))

message("DONE. Outputs written to: ", out_dir)


################################################################
# scRNA-seq step 2: proximal tubule pseudobulk associations with the
#                   random forest-prioritized fibrosis and tubular
#                   atrophy metabolites (Figure 3C-D inputs)
#
#   1. Participant-level pseudobulk for proximal tubule S1/S2 cells:
#      mean of log-normalized expression across all cells from each
#      participant.
#   2. Plasma metabolomics as a participant x metabolite (chemical
#      name) matrix, median of duplicate measurements.
#   3. Top 20 named metabolites (unnamed X- compounds removed) by
#      random forest permutation importance for interstitial
#      fibrosis and for tubular atrophy (Supplemental Table 5).
#   4. For each metabolite and gene: expression ~ metabolite (lm),
#      BH-FDR across genes within each metabolite.
#   5. Heatmap matrix of beta coefficients: genes ranked by their best
#      FDR across the 20 metabolites (then largest |beta|), with
#      canonical PT markers always included; 40 genes in total.
#      MT-, S100 and ACTB genes are excluded.
#
#
# INPUTS (data/)
#   scRNA/PT_S1_S2_scRNAseq.rds        KPMP Seurat object (not redistributed)
#   Plasma_Metabolomics.xlsx           long format with CHEMICAL_NAME
#   Fibrosis_permutation_importance_annotated.csv
#   TA_permutation_importance_annotated.csv
#                                      columns: importance, chemical_name
################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(readxl)
  library(Matrix)
})

########################
# 0) PATHS AND PARAMETERS
########################

# ================================================================
# EDIT THIS ONE LINE to point to the folder that contains data/
# ================================================================
base_dir <- "."   # e.g. "D:/UC Berkeley/KPMP paper/repo"

data_dir <- file.path(base_dir, "data")
out_dir  <- file.path(base_dir, "scRNA_03_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) stop("Input not found. Looked for:\n  ", paste(paths, collapse = "\n  "))
  hit[1]
}

pt_rds      <- first_existing(file.path(data_dir, "scRNA", "PT_S1_S2_scRNAseq.rds"))
plasma_file <- first_existing(file.path(data_dir, "Plasma_Metabolomics.xlsx"))
plasma_sheet <- 1
fib_rf_file <- first_existing(file.path(data_dir, "Fibrosis_permutation_importance_annotated.csv"))
ta_rf_file  <- first_existing(file.path(data_dir, "TA_permutation_importance_annotated.csv"))

top_n_rf_metabolites <- 20
top_n_genes_heatmap  <- 40
exclude_gene_regex   <- "^(MT-|S100|ACTB$)"
pt_gene_must <- c("SLC34A1","CUBN","LRP2","SLC22A6","SLC22A8","SLC5A2","SLC12A1")

########################
# 1) HELPERS
########################

clean_met <- function(x){
  x <- as.character(x)
  x <- str_replace_all(x, "`", "")
  x <- str_trim(x)
  x <- str_replace(x, "\\s*\\[[0-9]+\\]\\s*$", "")   # drop trailing "[12345]"
  x <- str_replace_all(x, "\\s+", " ")
  x
}

# g27-10088 -> 27-10088
clean_id <- function(x){
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "`", "")
  x <- str_replace_all(x, "\\s+", "")
  x <- str_replace(x, "^g(?=[0-9]{1,3}-[0-9]{4,6}$)", "")
  x
}

# keep only IDs of the KPMP form "29-10006"
id_to_kpmp_like <- function(x){
  x <- clean_id(x)
  ok <- str_detect(x, "^\\d{1,3}-\\d{4,6}$")
  ifelse(ok, x, NA_character_)
}

########################
# 2) PLASMA: PARTICIPANT x METABOLITE MATRIX
########################

read_plasma <- function(path){
  ext <- tools::file_ext(path)
  if (tolower(ext) %in% c("xlsx","xls")) {
    df <- readxl::read_excel(path, sheet = plasma_sheet)
  } else if (tolower(ext) %in% c("csv")) {
    df <- readr::read_csv(path, show_col_types = FALSE)
  } else {
    stop("Unsupported plasma file extension: ", ext)
  }
  
  names(df) <- names(df) %>% str_replace_all("\\s+", " ") %>% str_trim()
  
  get_col <- function(candidates){
    hit <- intersect(candidates, names(df))
    if (length(hit) == 0) return(NA_character_)
    hit[[1]]
  }
  
  col_pid <- get_col(c("Participant ID","ParticipantID","participant_id","Participant"))
  col_ab  <- get_col(c("Abundance","abundance","Value","value"))
  col_met <- get_col(c("CHEMICAL_NAME","Chemical Name","chemical_name","ChemicalName"))
  
  if (is.na(col_pid) || is.na(col_ab) || is.na(col_met)) {
    stop("Plasma data is missing required columns. Found columns: ", paste(names(df), collapse=", "))
  }
  
  df %>%
    transmute(
      SampleID_clean = clean_id(.data[[col_pid]]),
      met_clean      = clean_met(.data[[col_met]]),
      Abundance      = as.numeric(.data[[col_ab]])
    )
}

plasma_long <- read_plasma(plasma_file)

plasma_matrix <- plasma_long %>%
  group_by(SampleID_clean, met_clean) %>%
  summarize(Abundance = median(Abundance, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = met_clean, values_from = Abundance)

message("Plasma matrix (participants x metabolites): ", nrow(plasma_matrix), " x ", (ncol(plasma_matrix)-1))

########################
# 3) PT PSEUDOBULK
########################

# participants x genes: mean of the log-normalized "data" layer per participant
build_pseudobulk <- function(obj, id_priority = c("SampleID","KPMPID","DonorID","SpecimenID","orig.ident")){
  md <- obj@meta.data
  if (is.null(md)) stop("Could not find meta.data in scRNA object.")
  
  id_col <- intersect(id_priority, colnames(md))
  if (length(id_col) == 0) stop("No usable ID column found in meta.data. Available columns: ", paste(colnames(md), collapse=", "))
  id_col <- id_col[[1]]
  ids_kpmp <- id_to_kpmp_like(md[[id_col]])
  
  # if most IDs do not parse, try the other ID columns
  if (mean(is.na(ids_kpmp)) > 0.5 && length(id_priority) > 1){
    for (cand in id_priority[-1]){
      if (cand %in% colnames(md)){
        tmp <- id_to_kpmp_like(md[[cand]])
        if (mean(is.na(tmp)) < mean(is.na(ids_kpmp))){
          ids_kpmp <- tmp
          id_col <- cand
        }
      }
    }
  }
  if (all(is.na(ids_kpmp))) stop("Failed to derive KPMP participant IDs from meta.data.")
  
  assay <- Seurat::DefaultAssay(obj)
  expr <- tryCatch(Seurat::GetAssayData(obj, assay = assay, layer = "data"), error = function(e) NULL)
  if (is.null(expr) || ncol(expr) == 0) expr <- Seurat::GetAssayData(obj, assay = assay, layer = "counts")
  
  if (length(ids_kpmp) != ncol(expr)) {
    stop("meta.data rows (", length(ids_kpmp), ") do not match expression ncol (", ncol(expr), ").")
  }
  
  donor_levels <- sort(unique(na.omit(ids_kpmp)))
  donor_index <- match(ids_kpmp, donor_levels)
  keep <- !is.na(donor_index)
  
  expr_keep <- expr[, keep, drop = FALSE]
  donor_index <- donor_index[keep]
  
  M <- Matrix::sparseMatrix(
    i = donor_index,
    j = seq_along(donor_index),
    x = 1,
    dims = c(length(donor_levels), length(donor_index))
  )
  # expr_keep is genes x cells; M is participants x cells (1 = cell belongs
  # to participant). expr_keep %*% t(M) sums each gene over a participant's cells.
  gene_by_donor_sum <- Matrix::t(expr_keep %*% Matrix::t(M))  # participants x genes (sum)
  donor_counts <- Matrix::rowSums(M)
  
  donor_by_gene_mean <- as.matrix(gene_by_donor_sum / donor_counts)
  rownames(donor_by_gene_mean) <- donor_levels
  
  message("Pseudobulk built using ID column: ", id_col, " | participants x genes: ",
          nrow(donor_by_gene_mean), " x ", ncol(donor_by_gene_mean))
  donor_by_gene_mean
}

if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required to read the scRNA-seq objects.")
pt_obj  <- readRDS(pt_rds)
pt_expr <- build_pseudobulk(pt_obj)
rm(pt_obj); invisible(gc())

pt_ids <- intersect(plasma_matrix$SampleID_clean, rownames(pt_expr))
message("Participants with PT pseudobulk and plasma: ", length(pt_ids))

pt_expr2 <- pt_expr[pt_ids, , drop = FALSE]

plasma_pt_wide <- plasma_matrix %>%
  tibble::column_to_rownames("SampleID_clean") %>%
  as.data.frame()
plasma_pt_wide <- plasma_pt_wide[pt_ids, , drop = FALSE]

########################
# 4) RANDOM FOREST METABOLITE LISTS
########################

get_rf_mets <- function(path, top_n){
  df <- readr::read_csv(path, show_col_types = FALSE)
  names(df) <- str_to_lower(names(df))
  col_name <- if ("chemical_name" %in% names(df)) "chemical_name" else if ("metabolite" %in% names(df)) "metabolite" else NA_character_
  col_imp  <- if ("importance" %in% names(df)) "importance" else if ("perm_importance" %in% names(df)) "perm_importance" else NA_character_
  if (is.na(col_name) || is.na(col_imp)) stop("RF file missing required columns: ", path)
  
  df %>%
    transmute(chemical_name = clean_met(.data[[col_name]]),
              importance = as.numeric(.data[[col_imp]])) %>%
    filter(!is.na(importance)) %>%
    arrange(desc(importance)) %>%
    filter(!str_detect(chemical_name, "^X-")) %>%
    distinct(chemical_name, .keep_all = TRUE) %>%
    slice_head(n = top_n) %>%
    pull(chemical_name)
}

fib_mets <- get_rf_mets(fib_rf_file, top_n_rf_metabolites)
ta_mets  <- get_rf_mets(ta_rf_file,  top_n_rf_metabolites)
message("RF metabolites: fibrosis = ", length(fib_mets), " | tubular atrophy = ", length(ta_mets))

########################
# 5) GENE ~ METABOLITE REGRESSION
########################

run_gene_met_regression <- function(expr_donor_by_gene, plasma_sub, met_vec, label){
  stopifnot(all(rownames(plasma_sub) == rownames(expr_donor_by_gene)))
  
  met_present <- intersect(met_vec, colnames(plasma_sub))
  if (length(met_present) == 0) stop("No requested metabolites found in plasma matrix for ", label)
  if (length(met_present) < length(met_vec)) {
    message("NOTE (", label, "): not in plasma matrix: ", paste(setdiff(met_vec, met_present), collapse = "; "))
  }
  
  genes <- colnames(expr_donor_by_gene)
  out_list <- vector("list", length(met_present))
  names(out_list) <- met_present
  
  for (m in met_present) {
    x <- plasma_sub[[m]]
    ok <- !is.na(x)
    if (sum(ok) < 10) {
      message("Skipping ", label, " metabolite = ", m, " (fewer than 10 participants)")
      next
    }
    
    res_m <- lapply(seq_along(genes), function(j){
      y <- expr_donor_by_gene[, j]
      ok2 <- ok & !is.na(y)
      if (sum(ok2) < 10) return(NULL)
      fit <- stats::lm(y[ok2] ~ x[ok2])
      sm <- summary(fit)
      if (nrow(sm$coefficients) < 2) return(NULL)   # gene constant across participants
      data.frame(gene = genes[j],
                 beta = unname(coef(fit)[2]),
                 p = unname(sm$coefficients[2, 4]),
                 n = sum(ok2),
                 metabolite = m,
                 stringsAsFactors = FALSE)
    })
    res_m <- bind_rows(res_m)
    if (nrow(res_m) == 0) next
    res_m$FDR <- p.adjust(res_m$p, method = "fdr")
    out_list[[m]] <- res_m
    message(label, " regression done: ", m, " | genes = ", nrow(res_m))
  }
  
  bind_rows(out_list)
}

# Genes for the heatmap: best FDR across metabolites, then largest |beta|;
# must-include markers first, 40 genes in total
pick_top_genes <- function(reg_df, met_vec, top_n, must_include = character(0), exclude_regex = NULL){
  dd <- reg_df %>%
    filter(metabolite %in% met_vec) %>%
    mutate(gene = as.character(gene), beta = as.numeric(beta), FDR = as.numeric(FDR))
  if (!is.null(exclude_regex) && nzchar(exclude_regex)) dd <- dd %>% filter(!str_detect(gene, exclude_regex))
  
  best <- dd %>%
    group_by(gene) %>%
    summarize(best_FDR = min(FDR, na.rm = TRUE),
              best_abs_beta = max(abs(beta), na.rm = TRUE), .groups = "drop") %>%
    arrange(best_FDR, desc(best_abs_beta))
  
  genes <- unique(c(must_include, best$gene))
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) > top_n) {
    keep_n <- max(top_n - length(must_include), 0)
    genes <- unique(c(must_include, best$gene[seq_len(min(keep_n, nrow(best)))]))
  }
  genes
}

reg_to_beta_mat <- function(reg_df, genes, mets){
  dd <- reg_df %>%
    filter(gene %in% genes, metabolite %in% mets) %>%
    select(gene, metabolite, beta) %>%
    distinct()
  mat <- tidyr::pivot_wider(dd, names_from = metabolite, values_from = beta) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  mat[, intersect(mets, colnames(mat)), drop = FALSE]
}

run_block <- function(met_vec, label, matrix_file){
  res <- run_gene_met_regression(pt_expr2, plasma_pt_wide, met_vec, paste0("PT_", label))
  write_csv(res, file.path(out_dir, paste0("PT_gene_by_gene_assoc_RF_", label, ".csv")))
  
  met_present <- intersect(met_vec, colnames(plasma_pt_wide))
  genes <- pick_top_genes(res, met_present, top_n_genes_heatmap,
                          must_include = pt_gene_must, exclude_regex = exclude_gene_regex)
  mat <- reg_to_beta_mat(res, genes, met_present)
  
  write.csv(data.frame(Gene = rownames(mat), mat, check.names = FALSE),
            file.path(out_dir, matrix_file), row.names = FALSE)
  message("Saved ", matrix_file, ": ", nrow(mat), " genes x ", ncol(mat), " metabolites")
  invisible(mat)
}

########################
# 6) RUN
########################

run_block(fib_mets, "Fibrosis",       "PT_fibrosis_heatmap_matrix.csv")
run_block(ta_mets,  "TubularAtrophy", "PT_tubular_atrophy_heatmap_matrix.csv")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_scRNA_02.txt"))
message("DONE. Outputs saved to: ", out_dir)

