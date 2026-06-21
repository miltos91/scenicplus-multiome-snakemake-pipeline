library(argparse)
library(Seurat)
library(readr)

# ---- Parameters (passed in by Snakemake from config.yaml) ------------------
parser <- ArgumentParser()
parser$add_argument("--rds",            required = TRUE)
parser$add_argument("--gene_auc_csv",   required = TRUE)
parser$add_argument("--region_auc_csv", required = TRUE)
parser$add_argument("--out_rds",        required = TRUE)
args <- parser$parse_args()

# Load the Seurat object (single-object .rds file)
seurat_obj <- readRDS(args$rds)

gene_based   <- as.data.frame(read_csv(args$gene_auc_csv,   show_col_types = FALSE))
region_based <- as.data.frame(read_csv(args$region_auc_csv, show_col_types = FALSE))

# The first column holds the cell barcodes; name it "Cell"
colnames(gene_based)[1]   <- "Cell"
colnames(region_based)[1] <- "Cell"

gene_based$Cell   <- gsub("___cisTopic", "", gene_based$Cell)
region_based$Cell <- gsub("___cisTopic", "", region_based$Cell)

rownames(gene_based)   <- gene_based$Cell
rownames(region_based) <- region_based$Cell

gene_based$Cell <- NULL
region_based$Cell <- NULL

gene_based   <- t(gene_based)
region_based <- t(region_based)

# Restore TF(+/+) naming that was flattened to TF_+_+ during GMT export
rownames(gene_based)   <- sub("_(?=[^_]*$)", ")", rownames(gene_based),   perl = TRUE)
rownames(region_based) <- sub("_(?=[^_]*$)", ")", rownames(region_based), perl = TRUE)
rownames(gene_based)   <- sub("_(?=[^_]*$)", "(", rownames(gene_based),   perl = TRUE)
rownames(region_based) <- sub("_(?=[^_]*$)", "(", rownames(region_based), perl = TRUE)

gene_based   <- as.matrix(gene_based)
region_based <- as.matrix(region_based)

gene_based   <- gene_based[,   colnames(seurat_obj), drop = FALSE]
region_based <- region_based[, colnames(seurat_obj), drop = FALSE]

stopifnot(identical(colnames(seurat_obj), colnames(gene_based)))
stopifnot(identical(colnames(seurat_obj), colnames(region_based)))

seurat_obj[["SCENICplus_Gene_AUC"]]    <- CreateAssayObject(data = gene_based)
seurat_obj[["SCENICplus_Regions_AUC"]] <- CreateAssayObject(data = region_based)

rm(gene_based, region_based)

# Save the updated object to the output path (does not touch the input file)
dir.create(dirname(args$out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(seurat_obj, args$out_rds)
