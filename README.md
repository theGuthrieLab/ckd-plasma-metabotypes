# ckd-plasma-metabotypes

Analysis code for the study *"Host, diet, microbiome, and xenobiotic plasma metabolites associate with distinct kidney cell-type transcriptional programs in chronic kidney disease"* (Nowshad,F., Su, M., and Guthrie, L.).

This repository contains R scripts that reproduce the figures from the paper by integrating plasma metabolomics with kidney single-cell transcriptomics in the Kidney Precision Medicine Project (KPMP) cohort and two validation cohorts (MDRD and ALTOLD). Each script (`Main_Figure1.R`–`Main_Figure6.R`) regenerates the corresponding manuscript figure from the input files deposited in Dryad (see Data availability below).

## Data availability

All input data files required to run these scripts are deposited in the Dryad repository: https://doi.org/10.5061/dryad.18931zdcc

To reproduce the figures, download the files from Dryad and place them in a `data/` subfolder within each script's working directory (set via the `base_dir` variable at the top of each script). The Figure 3 and 5 script additionally requires author-constructed single-cell Seurat objects in `data/scRNA/`; these are not deposited owing to file size (see the Dryad README for instructions).

Source data originate from the Kidney Precision Medicine Project (KPMP; https://www.kpmp.org) and the MDRD and ALTOLD validation cohorts.
