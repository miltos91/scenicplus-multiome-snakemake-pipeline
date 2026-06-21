configfile: "config.yaml"

import os

MASTER_CONFIG = os.path.abspath("config.yaml")

RUN      = config["run_name"]
OUT      = config["output_base"]
CONDA_PY = config["conda_env_python"]
CONDA_R  = config.get("conda_env_r", "")

# Prefix for R steps: use conda run if an env is specified, else call Rscript directly
R_RUN = f"conda run --no-capture-output -n {CONDA_R} " if CONDA_R else ""

# ── Directory layout ─────────────────────────────────────────
RNA_DIR     = os.path.join(OUT, "1_object_files", "scRNAseq")
ATAC_DIR    = os.path.join(OUT, "1_object_files", "scATACseq")
GEN_DIR     = os.path.join(OUT, "2_generated_files")
MALLET_DIR  = os.path.join(GEN_DIR, "mallet_files")
RS_DIR      = os.path.join(GEN_DIR, "region_sets")
DAR_DIR     = os.path.join(RS_DIR, "DARs")
OUTS_DIR    = os.path.join(GEN_DIR, "outs")
TMP_DIR     = os.path.join(OUT, "tmp")
DB_PREFIX   = f"{RUN}_db"
DB_DIR      = os.path.join(config["ctx_db_base"], DB_PREFIX)
SCPLUS_DIR  = os.path.join(OUT, "3_scplus_pipeline")
AUC_DIR     = os.path.join(OUT, "4_auc_output")
MARKERS_DIR = os.path.join(OUT, "5_findmarkers")

SCPLUS_CONFIG = os.path.join(SCPLUS_DIR, "Snakemake", "config", "config.yaml")
SCPLUS_H5MU   = os.path.join(OUTS_DIR, f"scplusmdata_{RUN}.h5mu")

# ── Helper to format identity pairs for shell arg ────────────
def fmt_pairs(pairs):
    return ";".join(f"{p[0]},{p[1]}" for p in pairs)

# ── Final target ─────────────────────────────────────────────
FINAL = (
    os.path.join(MARKERS_DIR, ".done")
    if config.get("run_findmarkers", False)
    else os.path.join(AUC_DIR, "seurat_with_auc.rds")
)

rule all:
    input: FINAL

# ─────────────────────────────────────────────────────────────
# Step 1: Seurat → h5ad + ATACseq MTX
# ─────────────────────────────────────────────────────────────
rule convert_to_h5ad:
    input:
        rds = config["seurat_rds"]
    output:
        rna_h5ad     = os.path.join(RNA_DIR,  "adata_with_SCT_withraw.h5ad"),
        atac_counts  = os.path.join(ATAC_DIR, "counts.mtx"),
        atac_regions = os.path.join(ATAC_DIR, "regions.tsv"),
        atac_cells   = os.path.join(ATAC_DIR, "cells.tsv"),
    shell:
        """
        {R_RUN}Rscript scripts/1_convert_to_h5ad.R \
            --rds                {input.rds} \
            --rna_outdir         {RNA_DIR} \
            --atac_outdir        {ATAC_DIR}
        """

# ─────────────────────────────────────────────────────────────
# Step 2: Topic modelling (MALLET LDA)
# ─────────────────────────────────────────────────────────────
rule topic_modelling:
    input:
        counts  = os.path.join(ATAC_DIR, "counts.mtx"),
        regions = os.path.join(ATAC_DIR, "regions.tsv"),
        cells   = os.path.join(ATAC_DIR, "cells.tsv"),
    output:
        best_model  = os.path.join(MALLET_DIR, "best_model.pkl"),
        topics_top  = directory(os.path.join(RS_DIR, "Topics_top3k")),
        topics_otsu = directory(os.path.join(RS_DIR, "Topics_otsu")),
    params:
        n_topics = " ".join(str(t) for t in config["n_topics"]),
    shell:
        """
        python scripts/2_topic_modelling.py \
            --atac_dir      {ATAC_DIR} \
            --out_dir       {GEN_DIR} \
            --mallet_bin    {config[mallet_bin]} \
            --mallet_memory {config[mallet_memory]} \
            --n_topics      {params.n_topics} \
            --n_cpu         {config[n_cpu]} \
            --n_iter        {config[n_iter]} \
            --seed          {config[seed]}
        """

# ─────────────────────────────────────────────────────────────
# Step 3: DARs → BED files
# ─────────────────────────────────────────────────────────────
DAR_IS_BED = bool(config.get("dar_input_is_bed", False))
DAR_FORMAT = "bed" if DAR_IS_BED else "csv"

def dar_input_files():
    pattern = "{pop}.bed" if DAR_IS_BED else "DARs_{pop}.csv"
    return expand(os.path.join(config["dar_input_dir"], pattern),
                  pop = config["dar_populations"])

rule dars_to_bed:
    input:
        dar_files = dar_input_files()
    output:
        beds = expand(os.path.join(DAR_DIR, "{pop}.bed"),
                      pop = config["dar_populations"])
    params:
        populations = ",".join(config["dar_populations"]),
        fmt         = DAR_FORMAT,
    shell:
        """
        {R_RUN}Rscript scripts/3_dars_to_bed.R \
            --input_dir    {config[dar_input_dir]} \
            --output_dir   {DAR_DIR} \
            --populations  "{params.populations}" \
            --input_format {params.fmt}
        """

# ─────────────────────────────────────────────────────────────
# Step 4A: One-time cistarget setup (skipped if already done)
# ─────────────────────────────────────────────────────────────
rule one_time_setup:
    output:
        sentinel = os.path.join(config["ctx_db_base"], ".setup_complete")
    shell:
        "bash scripts/4A_one_time_setup.sh {config[ctx_db_base]}"

# ─────────────────────────────────────────────────────────────
# Step 4B: Build cistarget database
# ─────────────────────────────────────────────────────────────
rule build_cistarget_db:
    input:
        sentinel    = os.path.join(config["ctx_db_base"], ".setup_complete"),
        topics_top  = os.path.join(RS_DIR, "Topics_top3k"),
        topics_otsu = os.path.join(RS_DIR, "Topics_otsu"),
        dar_beds    = expand(os.path.join(DAR_DIR, "{pop}.bed"),
                             pop = config["dar_populations"]),
    output:
        rankings = os.path.join(DB_DIR, f"{DB_PREFIX}.regions_vs_motifs.rankings.feather"),
        scores   = os.path.join(DB_DIR, f"{DB_PREFIX}.regions_vs_motifs.scores.feather"),
    shell:
        """
        # positional args — order must match the block at the top of 4B_build_cistarget_db.sh:
        #   ctx_db_base db_dir db_prefix topics_top topics_otsu dar_dir genome_fa chromsizes n_threads
        bash scripts/4B_build_cistarget_db.sh \
            {config[ctx_db_base]} \
            {DB_DIR} \
            {DB_PREFIX} \
            {RS_DIR}/Topics_top3k \
            {RS_DIR}/Topics_otsu \
            {DAR_DIR} \
            {config[genome_fa]} \
            {config[chromsizes]} \
            {config[n_cpu]}
        """

# ─────────────────────────────────────────────────────────────
# Step 5: Init SCENIC+ Snakemake + generate its config
# ─────────────────────────────────────────────────────────────
rule init_scenicplus:
    input:
        best_model = os.path.join(MALLET_DIR, "best_model.pkl"),
        rna_h5ad   = os.path.join(RNA_DIR, "adata_with_SCT_withraw.h5ad"),
        rankings   = os.path.join(DB_DIR, f"{DB_PREFIX}.regions_vs_motifs.rankings.feather"),
        scores     = os.path.join(DB_DIR, f"{DB_PREFIX}.regions_vs_motifs.scores.feather"),
    output:
        scplus_cfg = SCPLUS_CONFIG
    shell:
        """
        mkdir -p {SCPLUS_DIR}
        cd {SCPLUS_DIR}
        scenicplus init_snakemake --out_dir Snakemake

        python scripts/generate_scenicplus_config.py \
            --master_config {MASTER_CONFIG} \
            --out           {output.scplus_cfg} \
            --best_model    {input.best_model} \
            --rna_h5ad      {input.rna_h5ad} \
            --region_sets   {RS_DIR} \
            --rankings      {input.rankings} \
            --scores        {input.scores} \
            --motif_annot   {config[motif_annotations]} \
            --outs_dir      {OUTS_DIR} \
            --tmp_dir       {TMP_DIR}
        """

# ─────────────────────────────────────────────────────────────
# Step 6: Run SCENIC+ Snakemake
# ─────────────────────────────────────────────────────────────
rule run_scenicplus:
    input:
        scplus_cfg = SCPLUS_CONFIG
    output:
        h5mu = SCPLUS_H5MU
    shell:
        """
        cd {SCPLUS_DIR}/Snakemake
        snakemake \
            --snakefile     workflow/Snakefile \
            --configfile    config/config.yaml \
            --cores         {config[n_cpu]} \
            --use-conda \
            --rerun-incomplete \
            --nolock \
            --printshellcmds
        """

# ─────────────────────────────────────────────────────────────
# Step 7: Export AUC matrices + GMT
# ─────────────────────────────────────────────────────────────
rule export_auc:
    input:
        h5mu = SCPLUS_H5MU
    output:
        gene_auc   = os.path.join(AUC_DIR, "eRegulon_gene_AUC.csv"),
        region_auc = os.path.join(AUC_DIR, "eRegulon_regions_AUC.csv"),
        gmt        = os.path.join(AUC_DIR, "Gene_based_eRegulons.gmt"),
    shell:
        """
        python scripts/7_export_auc.py \
            --work_dir {OUTS_DIR} \
            --out_dir  {AUC_DIR} \
            --run_name {RUN}
        """

# ─────────────────────────────────────────────────────────────
# Step 8: Add AUC assays to Seurat object
# ─────────────────────────────────────────────────────────────
rule add_auc_to_seurat:
    input:
        rds        = config["seurat_rds"],
        gene_auc   = os.path.join(AUC_DIR, "eRegulon_gene_AUC.csv"),
        region_auc = os.path.join(AUC_DIR, "eRegulon_regions_AUC.csv"),
    output:
        out_rds = os.path.join(AUC_DIR, "seurat_with_auc.rds")
    shell:
        """
        {R_RUN}Rscript scripts/8_add_auc_to_seurat.R \
            --rds            {input.rds} \
            --gene_auc_csv   {input.gene_auc} \
            --region_auc_csv {input.region_auc} \
            --out_rds        {output.out_rds}
        """

# ─────────────────────────────────────────────────────────────
# Step 9: FindMarkers + rnk files (conditional)
# ─────────────────────────────────────────────────────────────
if config.get("run_findmarkers", False):
    rule findmarkers:
        input:
            rds = os.path.join(AUC_DIR, "seurat_with_auc.rds")
        output:
            done = os.path.join(MARKERS_DIR, ".done")
        params:
            pairs        = fmt_pairs(config["identity_pairs"]),
            protein_flag = "--protein_only" if config.get("check_only_protein_coding") else "",
            protein_dir  = config.get("protein_coding_genes_dir", ""),
            species_r    = config.get("species_r", "mouse"),
        shell:
            """
            {R_RUN}Rscript scripts/9_findmarkers.R \
                --rds              {input.rds} \
                --assay            SCENICplus_Gene_AUC \
                --ident_col        {config[col_type_and_group]} \
                --identity_pairs   "{params.pairs}" \
                --out_dir          {MARKERS_DIR} \
                --min_cells        {config[min_cells_per_ident]} \
                {params.protein_flag} \
                --protein_gene_dir "{params.protein_dir}" \
                --species          {params.species_r}
            touch {output.done}
            """
