library(argparse)
library(Seurat)

### PREREQUISITE: PrepSCTFindMarkers must have been run if using the SCT assay.
### For the SCENICplus AUC assays used here this is not needed.

# ---- Parameters (passed in by Snakemake from config.yaml) ------------------
# These replace the hard-coded Ident_to_use / Ident_pair_N / directories of the
# original script. identity_pairs arrives as "A,B;C,D" (pairs split by ';').
parser <- ArgumentParser()
parser$add_argument("--rds",              required = TRUE)
parser$add_argument("--assay",            default  = "SCENICplus_Gene_AUC")
parser$add_argument("--ident_col",        required = TRUE)
parser$add_argument("--identity_pairs",   required = TRUE)
parser$add_argument("--out_dir",          required = TRUE)
parser$add_argument("--min_cells",        type = "integer", default = 10L)
parser$add_argument("--protein_only",     action = "store_true", default = FALSE)
parser$add_argument("--protein_gene_dir", default = "")
parser$add_argument("--species",          default = "mouse")
args <- parser$parse_args()

main_dir            <- args$out_dir
Ident_to_use        <- args$ident_col
min_cells_per_ident <- args$min_cells

# Identity pairs, e.g. list(c("CThPN_mutant","CThPN_control"), c(...), ...)
pairs_raw      <- strsplit(args$identity_pairs, ";")[[1]]
identity_pairs <- lapply(trimws(pairs_raw), function(p) trimws(strsplit(p, ",")[[1]]))

################################################################################################
### Load object
################################################################################################
seurat_obj <- readRDS(args$rds)

DefaultAssay(seurat_obj) <- args$assay
Idents(seurat_obj)       <- Ident_to_use

################################################################################################
### Read protein coding genes (optional)
################################################################################################
prot_genes <- NULL
if (args$protein_only && nzchar(args$protein_gene_dir)) {
  file_name  <- paste0(args$species, "_protein_coding_symbols_uniprot.txt")
  prot_genes <- readLines(file.path(args$protein_gene_dir, file_name))
}

seurat_features <- rownames(seurat_obj)
features_use    <- if (args$protein_only && !is.null(prot_genes))
                     intersect(prot_genes, seurat_features) else NULL

################################################################################################
### Folder preparation
################################################################################################
dir_save_total_markers <- file.path(main_dir, "All_genes/")
dir_save_rnk_files     <- file.path(main_dir, "rnks/")
dir_save_DEG_markers   <- file.path(main_dir, "Significant_DEGs/")

for (d in c(dir_save_total_markers, dir_save_rnk_files, dir_save_DEG_markers)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

################################################################################################
### FindMarkers and save csv
################################################################################################
markers_list <- list()

for (pair in identity_pairs) {
  ident_1st <- pair[1]
  ident_2nd <- pair[2]
  markers_name <- paste0("markers_", ident_1st, "_vs_", ident_2nd)

  # how many cells per identity?
  cells1 <- tryCatch(WhichCells(seurat_obj, idents = ident_1st), error = function(e) character(0))
  cells2 <- tryCatch(WhichCells(seurat_obj, idents = ident_2nd), error = function(e) character(0))

  if (length(cells1) < min_cells_per_ident || length(cells2) < min_cells_per_ident) {
    message(sprintf("Skipping %s vs %s: %d vs %d cells (< %d).",
                    ident_1st, ident_2nd, length(cells1), length(cells2), min_cells_per_ident))
    next
  }

  message("Finding markers: ", ident_1st, " vs ", ident_2nd)
  markers <- FindMarkers(
    object = seurat_obj,
    ident.1 = ident_1st,
    ident.2 = ident_2nd,
    features = features_use,
    min.cells.group = min_cells_per_ident,
    logfc.threshold = 0,
    min.pct = 1e-4
  )

  # drop non-activating eRegulons (-/-, +/-, -/+); keep only (+/+)
  exclude_pattern <- "-/-|\\+/-|-/\\+"
  markers <- markers[!grepl(exclude_pattern, rownames(markers)), ]

  markers_list[[markers_name]] <- markers
}

for (markers_name in names(markers_list)) {
  write.csv(markers_list[[markers_name]],
            file = paste0(dir_save_total_markers, markers_name, ".csv"))
}

################################################################################################
### Make and save rnk files
################################################################################################
for (markers_name in names(markers_list)) {

  df <- markers_list[[markers_name]]

  GSEA_table <- df[, c("avg_log2FC", "p_val")]

  GSEA_table$sign <- sign(GSEA_table$avg_log2FC)

  GSEA_table <- GSEA_table[, c("p_val", "sign")]

  GSEA_table$preranked <- -10 * log10(GSEA_table[, "p_val"]) * GSEA_table[, "sign"]

  smallest_value <- min(GSEA_table$preranked[is.finite(GSEA_table$preranked)], na.rm = TRUE)
  highest_value  <- max(GSEA_table$preranked[is.finite(GSEA_table$preranked)], na.rm = TRUE)

  GSEA_table$preranked[GSEA_table$preranked == Inf]  <- highest_value + 0.000000000000000001
  GSEA_table$preranked[GSEA_table$preranked == -Inf] <- smallest_value - 0.000000000000000001

  GSEA_table <- GSEA_table[order(GSEA_table$preranked, decreasing = TRUE), ]

  GSEA_table <- GSEA_table[, "preranked", drop = FALSE]

  rnk_path <- paste0(dir_save_rnk_files, markers_name, ".rnk")
  write.table(GSEA_table, file = rnk_path, sep = "\t", quote = FALSE, col.names = FALSE)

  cat("Wrote", markers_name, "to", rnk_path, "\n")
}

################################################################################################
### Make and save significant DEG csv (up / down)
################################################################################################
for (markers_name in names(markers_list)) {

  df <- markers_list[[markers_name]]

  up_name   <- sub("_vs_", "_up_vs_",   markers_name, fixed = TRUE)
  down_name <- sub("_vs_", "_down_vs_", markers_name, fixed = TRUE)

  write.csv(subset(df, avg_log2FC > 0 & p_val_adj <= 0.05),
            file = paste0(dir_save_DEG_markers, up_name, ".csv"))
  write.csv(subset(df, avg_log2FC < 0 & p_val_adj <= 0.05),
            file = paste0(dir_save_DEG_markers, down_name, ".csv"))
}
