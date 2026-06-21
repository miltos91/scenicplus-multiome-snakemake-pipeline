import argparse, os
import anndata as ad
import mudata
import pandas as pd

# ---- Parameters (passed in by Snakemake from config.yaml) ------------------
parser = argparse.ArgumentParser()
parser.add_argument("--work_dir",  required=True)
parser.add_argument("--out_dir",   required=True)
parser.add_argument("--run_name",  required=True)
args = parser.parse_args()

work_dir = args.work_dir
out_dir  = args.out_dir
run_name = args.run_name
os.makedirs(out_dir, exist_ok=True)

################################################################

scplus_mdata = mudata.read(os.path.join(work_dir, f"scplusmdata_{run_name}.h5mu"))

eRegulon_gene_AUC = ad.concat(
    [scplus_mdata["direct_gene_based_AUC"],
     scplus_mdata["extended_gene_based_AUC"]],
    axis=1,
)

eRegulon_regions_AUC = ad.concat(
    [scplus_mdata["direct_region_based_AUC"],
     scplus_mdata["extended_region_based_AUC"]],
    axis=1,
)

gene_auc_df = pd.DataFrame(
    eRegulon_gene_AUC.X,
    index=eRegulon_gene_AUC.obs_names,
    columns=eRegulon_gene_AUC.var_names
)

region_auc_df = pd.DataFrame(
    eRegulon_regions_AUC.X,
    index=eRegulon_regions_AUC.obs_names,
    columns=eRegulon_regions_AUC.var_names
)

gene_auc_df.to_csv(os.path.join(out_dir, "eRegulon_gene_AUC.csv"))
region_auc_df.to_csv(os.path.join(out_dir, "eRegulon_regions_AUC.csv"))


##########################################################
### Save metadata eRegulon list
#########################

# 1. Retrieve the metadata DataFrames
direct_meta   = scplus_mdata.uns['direct_e_regulon_metadata']
extended_meta = scplus_mdata.uns['extended_e_regulon_metadata']

# 2. Combine them into one DataFrame, keep only activating (+/+) eRegulons
combined_meta = pd.concat([direct_meta, extended_meta])
filtered_meta = combined_meta[combined_meta['Gene_signature_name'].str.contains("+/+", regex=False)].copy()
filtered_meta['Gene_signature_name'] = filtered_meta['Gene_signature_name'].str.replace('/', '_', regex=False)
filtered_meta['Gene_signature_name'] = filtered_meta['Gene_signature_name'].str.replace(r'[()]', '', regex=True)

print(f"Original rows: {len(combined_meta)}")
print(f"Rows after keeping only (+/+): {len(filtered_meta)}")
print(filtered_meta['Gene_signature_name'].head())

# 3. Define the output filename
gmt_filename = os.path.join(out_dir, "Gene_based_eRegulons.gmt")

# 4. Process and write to GMT
print(f"Exporting to: {gmt_filename}")

with open(gmt_filename, 'w') as f:
    for sig_name, group in filtered_meta.groupby('Gene_signature_name'):
        genes = group['Gene'].unique().tolist()
        # Build the tab-separated string first to avoid an f-string backslash error
        genes_str = "\t".join(genes)
        line = f"{sig_name}\tNA\t{genes_str}\n"
        f.write(line)

print("Export complete.")

# Optional: print the first few lines of the generated file to verify
with open(gmt_filename, 'r') as f:
    print("\nFirst 3 lines of the GMT file:")
    for i in range(3):
        print(f.readline().strip()[:100] + "...")  # truncated for display
