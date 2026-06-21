#!/bin/bash
########## ONE-TIME SETUP: cistarget database tools ##########
# Downloads cbust + the motif collection once. Safe to call repeatedly:
# it skips everything if the sentinel file already exists.
set -euo pipefail

CTX_DB_BASE="${1:?Usage: 4A_one_time_setup.sh <ctx_db_base_dir>}"
SENTINEL="$CTX_DB_BASE/.setup_complete"

if [ -f "$SENTINEL" ]; then
    echo "cistarget setup already done (sentinel: $SENTINEL). Skipping."
    exit 0
fi

# 1) Base folder for all cistarget stuff
mkdir -p "$CTX_DB_BASE"
cd "$CTX_DB_BASE"
echo "Using ctx DB folder: $CTX_DB_BASE"

# 2) Clone helper repo (once)
if [ ! -d "$CTX_DB_BASE/create_cisTarget_databases" ]; then
    git clone https://github.com/aertslab/create_cisTarget_databases
fi

# 3) Cluster-Buster binary
if [ ! -x "$CTX_DB_BASE/cbust" ]; then
    wget -O "$CTX_DB_BASE/cbust" \
        https://resources.aertslab.org/cistarget/programs/cbust
    chmod a+x "$CTX_DB_BASE/cbust"
fi

# 4) Motif collection v10nr_clust_public
if [ ! -d "$CTX_DB_BASE/aertslab_motif_collection" ]; then
    mkdir -p "$CTX_DB_BASE/aertslab_motif_collection"
    wget -O "$CTX_DB_BASE/aertslab_motif_collection/v10nr_clust_public.zip" \
        https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip
    (cd "$CTX_DB_BASE/aertslab_motif_collection" && unzip -q v10nr_clust_public.zip)
fi

touch "$SENTINEL"
echo "One-time setup complete in: $CTX_DB_BASE"
