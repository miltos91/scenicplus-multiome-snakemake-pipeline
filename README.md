# SCENIC+ Multiome Pipeline — eRegulon Networks & Differential eRegulon Analysis

*End-to-end Snakemake pipeline that runs SCENIC+ on a Seurat v5 multiome object to infer TF–enhancer–gene eRegulons and differential eRegulon activity.*

## Main goal of the pipeline
This pipeline takes a processed **Seurat v5 multiome object** (paired scRNA-seq + scATAC-seq) as input and runs, as a single Snakemake workflow, every step needed to:
1. Convert the Seurat object into the RNA (`.h5ad`) and ATAC (sparse matrix) inputs SCENIC+ expects.
2. Learn regulatory topics from the ATAC data with cisTopic (MALLET LDA) and turn differentially accessible regions (DARs) into region sets.
3. Build a **dataset-specific cisTarget motif database** from those regions.
4. Run **SCENIC+** to infer eRegulons (TF → enhancer → gene) and per-cell eRegulon activity (AUC).
5. Add the eRegulon AUC back into the Seurat object and run **differential eRegulon analysis** (FindMarkers + preranked `.rnk` files for GSEA).

Everything is configured in one file (`config.yaml`) and launched with one command (`bash submit.sh`), either locally or on an HPC via SLURM.

## Essential packages and versions

### Python environment (tested on 3.11.0)
SCENIC+ has complex dependencies so install it by following the official instructions
(https://github.com/aertslab/scenicplus). This pipeline was developed against a conda env named `scenicplus1a2`.

Required
- scenicplus 1.0a2
- pycisTopic 2.0a0
- pycistarget 1.1
- snakemake 8.5.5

Recommended (as tested)
- numpy 1.26.4
- pandas 1.5.3
- scipy 1.12.0
- anndata 0.10.5.post1
- mudata 0.2.3
- pyarrow 10.0.1
- scikit-learn 1.3.2
- PyYAML 6.0.1

### R environment (tested on 4.5.3)
Required
- Seurat >= 5.0.0 (tested 5.3.1)
- Signac 1.15.0

Recommended
- SingleCellExperiment 1.31.1
- zellkonverter 1.19.2
- reticulate 1.43.0
- Matrix 1.7.4
- argparse 2.3.1
- readr 2.1.5

### Command-line tools
- MALLET 2.0.8, topic modelling (Java); **included in this repo** under `mallet/`
- bedtools, consensus region merging
- Cluster-Buster (`cbust`) + the aertslab motif collection, downloaded automatically by step 4A

## Installation

### 1) Clone the repository
```bash
git clone git@github.com:miltos91/scenicplus-multiome-pipeline.git
cd scenicplus-multiome-pipeline
```

### 2) Build the Python environment (SCENIC+)
SCENIC+ is best installed from its official repository:
```bash
conda create -n scenicplus1a2 python=3.11 -y
conda activate scenicplus1a2

git clone https://github.com/aertslab/scenicplus
cd scenicplus
pip install .
```
To pin the exact environment used here:
```bash
conda env export -n scenicplus1a2 > environment.yml
```

### 3) Install R libraries
```bash
Rscript -e 'pkgs <- c("Seurat","Signac","reticulate","Matrix","argparse","readr"); \
miss <- setdiff(pkgs, rownames(installed.packages())); \
if (length(miss)) install.packages(miss, repos="https://cloud.r-project.org")'

# Bioconductor packages:
Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); \
BiocManager::install(c("SingleCellExperiment","zellkonverter"))'
```

### 4) MALLET (included)
MALLET 2.0.8 is bundled with this repository under `mallet/`, it is the specific
build this pipeline relies on, so there is nothing to download. Just point
`mallet_bin` in `config.yaml` at it, for example:
```yaml
mallet_bin: "/your/clone/path/mallet/mallet-2.0.8/bin/mallet"
```

### 5) cisTarget database assets (one-time)
On the first run, step 4A automatically downloads Cluster-Buster and the motif collection into `ctx_db_base`. You only need to provide the genome FASTA + chrom sizes for your species (set in `config.yaml`).

## Prerequisites
Before running, the input Seurat object must already have:
- an **SCT assay** (normalised RNA), and
- **DAR analysis** completed (Signac `FindDAR`), one file per population.

Save the object as a single `.rds`:
```r
saveRDS(seurat_obj, "/path/to/seurat_obj.rds")
```

## How to use

### 1) Configure config.yaml
`config.yaml` has all paths and parameters to set there:
- **Data**: the `.rds` object path, `output_base`, cistarget DB folder, genome files.
- **Populations**: the main population column, the DAR population names, and the FindMarkers identity pairs.
- **Compute**: `n_cpu`, `use_slurm` (local vs SLURM), MALLET memory.
- **SCENIC+**: all motif-enrichment and network-inference parameters.

The file has explanatory notes inline.

### 2) Run the full pipeline
```bash
bash submit.sh
```
- `use_slurm: true` → submits one SLURM job (HPC).
- `use_slurm: false` → runs Snakemake locally in your terminal.

To resume after an interruption, simply run `bash submit.sh`. Snakemake continues from wherever the existing output files leave off.

### Outputs (by stage), under `output_base/`
- `1_object_files/` RNA `.h5ad` + ATAC sparse matrix
- `2_generated_files/` topic models, region sets, and SCENIC+ `outs/` (eRegulons, AUC `.h5mu`)
- `3_scplus_pipeline/` the SCENIC+ Snakemake workflow
- `4_auc_output/` eRegulon AUC CSVs, GMT, and the Seurat object with AUC assays
- `5_findmarkers/` differential eRegulon tables + `.rnk` files for GSEA

### Reset / re-run from scratch
There is no hidden state. Snakemake decides what to skip from the output files that exist. To start over, delete the outputs:
```bash
rm -rf <output_base>
rm -rf <ctx_db_base>/<run_name>_db     # only if you want to rebuild the motif database
rm -rf .snakemake
```
Keep the rest of `ctx_db_base`. those are the one-time downloads.

## Complete steps of the pipeline

### 1) 1_convert_to_h5ad.R Seurat → SCENIC+ inputs
- Load the single-object `.rds` Seurat file
- Export the SCT RNA assay as AnnData (`.h5ad`), adding a raw copy
- Export the ATAC counts as a sparse matrix (`counts.mtx` + `regions.tsv` + `cells.tsv`)
- Rename ATAC regions from `chr1-100-200` to `chr1:100-200`

### 2) 2_topic_modelling.py cisTopic LDA (MALLET)
- Build a cisTopic object from the ATAC matrix
- Run LDA topic modelling with MALLET (one or several topic numbers)
- Select the best model (or use it directly when only one was trained)
- Binarize topics into region sets (otsu + top-3k) as BED files

### 3) 3_dars_to_bed.R DARs → BED
- For each population, convert Signac `FindDAR` CSVs into BED files of differentially accessible regions (or copy ready-made BEDs)

### 4A) 4A_one_time_setup.sh one-time database assets
- Download Cluster-Buster (`cbust`), the aertslab motif collection, and the cisTarget helper repo (skipped on later runs via a sentinel file)

### 4B) 4B_build_cistarget_db.sh dataset-specific cisTarget database
- Merge topic + DAR BEDs into a consensus region set
- Pad regions and build a background FASTA
- Score motifs with Cluster-Buster → cisTarget database (`.feather` rankings + scores)

### 5) generate_scenicplus_config.py SCENIC+ setup
- Initialize the SCENIC+ Snakemake workflow (`scenicplus init_snakemake`)
- Generate its `config.yaml` automatically from this pipeline's config

### 6) SCENIC+ Snakemake network inference
- Run SCENIC+: motif enrichment (cisTarget + DEM), TF→gene and region→gene links, eRegulon construction, and per-cell AUC scoring → `scplusmdata.h5mu`

### 7) 7_export_auc.py export AUC + GMT
- Extract gene-based and region-based eRegulon AUC matrices to CSV
- Export the activating (`+/+`) eRegulons as a GMT file

### 8) 8_add_auc_to_seurat.R AUC back into Seurat
- Add `SCENICplus_Gene_AUC` and `SCENICplus_Regions_AUC` assays to the Seurat object, saved as `.rds`

### 9) 9_findmarkers.R differential eRegulon analysis
- Run `FindMarkers` on the eRegulon AUC assay for each identity pair
- Export full marker tables, significant up/down eRegulons, and preranked `.rnk` files for GSEA
