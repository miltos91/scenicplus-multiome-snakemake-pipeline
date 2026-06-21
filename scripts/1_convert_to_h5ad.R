library(argparse)
library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)
library(Signac)
library(Matrix)

# ---- Parameters (passed in by Snakemake from config.yaml) ------------------
# These replace the hard-coded paths of the original script.
parser <- ArgumentParser()
parser$add_argument("--rds",                 required = TRUE)
parser$add_argument("--rna_outdir",          required = TRUE)
parser$add_argument("--atac_outdir",         required = TRUE)
args <- parser$parse_args()

rna_dir  <- args$rna_outdir
atac_dir <- args$atac_outdir
dir.create(rna_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(atac_dir, recursive = TRUE, showWarnings = FALSE)

# Load the Seurat object. The input is a single-object .rds file, so readRDS
# returns it directly (no name to look up).
seurat_obj <- readRDS(args$rds)

######################################################################################
############################## Convert RNAseq to h5ad ################################
######################################################################################

DefaultAssay(seurat_obj) <- "SCT"

sce <- as.SingleCellExperiment(seurat_obj, assay = "SCT")

zellkonverter::writeH5AD(sce, file = file.path(rna_dir, "adata_with_SCT.h5ad"), X_name = "logcounts")

anndata <- import("anndata", convert = FALSE)

in_file  <- file.path(rna_dir, "adata_with_SCT.h5ad")
out_file <- file.path(rna_dir, "adata_with_SCT_withraw.h5ad")

adata <- anndata$read_h5ad(in_file)
adata$raw <- adata$copy()
adata$write_h5ad(out_file)

######################################################################################
############################# Convert ATACseq to MTX ################################
######################################################################################

DefaultAssay(seurat_obj) <- "ATAC"
mat <- GetAssayData(seurat_obj, layer = "counts")

all_zero <- Matrix::rowSums(mat) == 0
mat_keep <- mat[!all_zero, , drop = FALSE]
rownames(mat_keep) <- sub("-", ":", rownames(mat_keep))
## Important: to avoid problems later, replace the style chr1-1324-34546 with chr1:1324-34546

Matrix::writeMM(mat_keep, file.path(atac_dir, "counts.mtx"))
write.table(rownames(mat_keep), file.path(atac_dir, "regions.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(colnames(mat_keep), file.path(atac_dir, "cells.tsv"),   row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

message("Step 1 done.")
