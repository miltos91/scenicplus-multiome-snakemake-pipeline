# ============================================================================
# Writes the SCENIC+ inner config.yaml (the original step 5A_config.yaml).
# Original 5A was a static YAML with hard-coded paths; here the same YAML is
# built from the master config.yaml + the file paths produced by earlier steps,
# so nothing has to be edited by hand. The dict below mirrors 5A section by
# section: input_data / output_data / params_general / params_data_preparation
# / params_motif_enrichment / params_inference.
# ============================================================================
import argparse, os, yaml

# SCENIC+ uses one species name for biomart (e.g. mmusculus) and another for the
# motif annotations (e.g. mus_musculus); map between them.
SPECIES_MOTIF_MAP = {
    "mmusculus":   "mus_musculus",
    "hsapiens":    "homo_sapiens",
    "rnorvegicus": "rattus_norvegicus",
}

# ---- Parameters (paths from Snakemake, everything else from config.yaml) ----
parser = argparse.ArgumentParser()
parser.add_argument("--master_config", required=True,  help="Path to pipeline/config.yaml")
parser.add_argument("--out",           required=True)
parser.add_argument("--best_model",    required=True)
parser.add_argument("--rna_h5ad",      required=True)
parser.add_argument("--region_sets",   required=True)
parser.add_argument("--rankings",      required=True)
parser.add_argument("--scores",        required=True)
parser.add_argument("--motif_annot",   required=True)
parser.add_argument("--outs_dir",      required=True)
parser.add_argument("--tmp_dir",       required=True)
args = parser.parse_args()

with open(args.master_config) as f:
    mc = yaml.safe_load(f)

def get(key, default=None):
    return mc.get(key, default)

RUN = get("run_name")
OUT = args.outs_dir
os.makedirs(OUT, exist_ok=True)
os.makedirs(args.tmp_dir, exist_ok=True)
os.makedirs(os.path.dirname(args.out), exist_ok=True)

is_multiome    = str(get("is_multiome", True)).lower() not in ("false", "0", "no")
species_motif  = SPECIES_MOTIF_MAP.get(get("species", "mmusculus"), get("species", "mmusculus"))

cfg = {
    "input_data": {
        "cisTopic_obj_fname":        args.best_model,
        "GEX_anndata_fname":         args.rna_h5ad,
        "region_set_folder":         args.region_sets,
        "ctx_db_fname":              args.rankings,
        "dem_db_fname":              args.scores,
        "path_to_motif_annotations": args.motif_annot,
    },
    "output_data": {
        "combined_GEX_ACC_mudata":      os.path.join(OUT, f"ACC_GEX_{RUN}.h5mu"),
        "dem_result_fname":             os.path.join(OUT, f"dem_results_{RUN}.hdf5"),
        "ctx_result_fname":             os.path.join(OUT, f"ctx_results_{RUN}.hdf5"),
        "output_fname_dem_html":        os.path.join(OUT, "dem_results.html"),
        "output_fname_ctx_html":        os.path.join(OUT, "ctx_results.html"),
        "cistromes_direct":             os.path.join(OUT, f"cistromes_direct_{RUN}.h5ad"),
        "cistromes_extended":           os.path.join(OUT, f"cistromes_extended_{RUN}.h5ad"),
        "tf_names":                     os.path.join(OUT, f"tf_names_{RUN}.txt"),
        "genome_annotation":            os.path.join(OUT, "genome_annotation.tsv"),
        "chromsizes":                   os.path.join(OUT, "chromsizes.tsv"),
        "search_space":                 os.path.join(OUT, f"search_space_{RUN}.tsv"),
        "tf_to_gene_adjacencies":       os.path.join(OUT, f"tf_to_gene_adj_{RUN}.tsv"),
        "region_to_gene_adjacencies":   os.path.join(OUT, f"region_to_gene_adj_{RUN}.tsv"),
        "eRegulons_direct":             os.path.join(OUT, f"eRegulons_direct_{RUN}.tsv"),
        "eRegulons_extended":           os.path.join(OUT, f"eRegulons_extended_{RUN}.tsv"),
        "AUCell_direct":                os.path.join(OUT, f"AUCell_direct_{RUN}.h5mu"),
        "AUCell_extended":              os.path.join(OUT, f"AUCell_extended_{RUN}.h5mu"),
        "scplus_mdata":                 os.path.join(OUT, f"scplusmdata_{RUN}.h5mu"),
    },
    "params_general": {
        "temp_dir": args.tmp_dir,
        "n_cpu":    get("n_cpu", 32),
        "seed":     get("seed", 666),
    },
    "params_data_preparation": {
        "bc_transform_func":          "\"lambda x: f'{x}___cisTopic'\"",
        "is_multiome":                is_multiome,
        "key_to_group_by":            "",
        "nr_cells_per_metacells":     get("nr_cells_per_metacells", 10),
        "direct_annotation":          "Direct_annot",
        "extended_annotation":        "Orthology_annot",
        "species":                    get("species", "mmusculus"),
        "biomart_host":               get("biomart_host", "nov2020.archive.ensembl.org"),
        "search_space_upstream":      get("search_space_upstream", "1000 150000"),
        "search_space_downstream":    get("search_space_downstream", "1000 150000"),
        "search_space_extend_tss":    get("search_space_extend_tss", "10 10"),
    },
    "params_motif_enrichment": {
        "species":                             species_motif,
        "annotation_version":                  get("annotation_version", "v10nr_clust"),
        "motif_similarity_fdr":                get("motif_similarity_fdr", 0.001),
        "orthologous_identity_threshold":      get("orthologous_identity_threshold", 0.0),
        "annotations_to_use":                  get("annotations_to_use", "Direct_annot Orthology_annot"),
        "fraction_overlap_w_dem_database":     get("fraction_overlap_w_dem_database", 0.4),
        "dem_max_bg_regions":                  get("dem_max_bg_regions", 500),
        "dem_balance_number_of_promoters":     get("dem_balance_number_of_promoters", True),
        "dem_promoter_space":                  get("dem_promoter_space", 1000),
        "dem_adj_pval_thr":                    get("dem_adj_pval_thr", 0.05),
        "dem_log2fc_thr":                      get("dem_log2fc_thr", 1.0),
        "dem_mean_fg_thr":                     get("dem_mean_fg_thr", 0.0),
        "dem_motif_hit_thr":                   get("dem_motif_hit_thr", 3.0),
        "fraction_overlap_w_ctx_database":     get("fraction_overlap_w_ctx_database", 0.4),
        "ctx_auc_threshold":                   get("ctx_auc_threshold", 0.005),
        "ctx_nes_threshold":                   get("ctx_nes_threshold", 3.0),
        "ctx_rank_threshold":                  get("ctx_rank_threshold", 0.05),
    },
    "params_inference": {
        "tf_to_gene_importance_method":       get("tf_to_gene_importance_method", "GBM"),
        "region_to_gene_importance_method":   get("region_to_gene_importance_method", "GBM"),
        "region_to_gene_correlation_method":  get("region_to_gene_correlation_method", "SR"),
        "order_regions_to_genes_by":          get("order_regions_to_genes_by", "importance"),
        "order_TFs_to_genes_by":              get("order_TFs_to_genes_by", "importance"),
        "gsea_n_perm":                        get("gsea_n_perm", 1000),
        "quantile_thresholds_region_to_gene": get("quantile_thresholds_region_to_gene", "0.85 0.90 0.95"),
        "top_n_regionTogenes_per_gene":       get("top_n_regionTogenes_per_gene", "5 10 15"),
        "top_n_regionTogenes_per_region":     get("top_n_regionTogenes_per_region", ""),
        "min_regions_per_gene":               get("min_regions_per_gene", 0),
        "rho_threshold":                      get("rho_threshold", 0.05),
        "min_target_genes":                   get("min_target_genes", 10),
    },
}

with open(args.out, "w") as f:
    yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)

print(f"SCENIC+ config written to: {args.out}")
