import argparse, os, pickle, random
import numpy as np
from scipy.io import mmread
from pycisTopic.cistopic_class import create_cistopic_object
from pycisTopic.lda_models import run_cgs_models_mallet, evaluate_models
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.utils import region_names_to_coordinates

# ---- Parameters (passed in by Snakemake from config.yaml) ------------------
parser = argparse.ArgumentParser()
parser.add_argument("--atac_dir",      required=True)
parser.add_argument("--out_dir",       required=True)
parser.add_argument("--mallet_bin",    required=True)
parser.add_argument("--mallet_memory", default="180g")
parser.add_argument("--n_topics",      nargs="+", type=int, default=[60])
parser.add_argument("--n_cpu",         type=int, default=32)
parser.add_argument("--n_iter",        type=int, default=500)
parser.add_argument("--seed",          type=int, default=666)
args = parser.parse_args()

random.seed(args.seed)
np.random.seed(args.seed)

# --- 0) Project paths --------------------------------------------------------
atac_dir   = args.atac_dir
out_dir    = args.out_dir
mallet_dir = os.path.join(out_dir, "mallet_files")
rs_dir     = os.path.join(out_dir, "region_sets")
os.makedirs(mallet_dir, exist_ok=True)
os.makedirs(os.path.join(rs_dir, "Topics_otsu"),  exist_ok=True)
os.makedirs(os.path.join(rs_dir, "Topics_top3k"), exist_ok=True)

# --- 1) cisTopic object ------------------------------------------------------
# Build a sparse matrix of cells x regions from the Signac ATAC counts exported
# in step 1, so we can run topic modelling on it.
X       = mmread(os.path.join(atac_dir, "counts.mtx")).tocsr().astype(np.int32)   # ensure integer counts
regions = [l.strip() for l in open(os.path.join(atac_dir, "regions.tsv"))]
cells   = [l.strip() for l in open(os.path.join(atac_dir, "cells.tsv"))]

assert X.shape == (len(regions), len(cells))

cobj = create_cistopic_object(
    fragment_matrix = X,
    region_names    = regions,
    cell_names      = cells
)

# --- 2) Topic modeling -------------------------------------------------------
# Identifies peaks that go together (topics) and, for each cell, which topics
# appear in it. Gives a weight of peaks per topic and a weight of topic per cell.
os.environ["MALLET_MEMORY"] = args.mallet_memory
os.environ["MALLET_OPTS"]   = "-XX:+UseG1GC"
os.makedirs("/tmp/mallet_tmp", exist_ok=True)

models = run_cgs_models_mallet(
    cobj,
    n_topics     = args.n_topics,            # try several numbers of topics, then pick the best
    n_cpu        = args.n_cpu,
    n_iter       = args.n_iter,
    random_state = args.seed,
    tmp_path     = "/tmp/mallet_tmp",
    save_path    = None,
    mallet_path  = args.mallet_bin
)

with open(os.path.join(mallet_dir, "Mallet_models.pkl"), "wb") as f:
    pickle.dump(models, f, protocol=pickle.HIGHEST_PROTOCOL)

# evaluate_models rescales each metric across the candidate models. With a single
# model there is nothing to compare (the rescaling divides by zero), so when only
# one number of topics was requested we just use that model directly.
if len(models) == 1:
    best_model = models[0]
else:
    best_model = evaluate_models(
        models,
        select_model = None,
        return_model = True,
        metrics      = ["Arun_2010", "Cao_Juan_2009", "Minmo_2011", "loglikelihood"],
        plot_metrics = False,
        save         = os.path.join(mallet_dir, "model_selection.pdf")
    )

cobj.add_LDA_model(best_model)

with open(os.path.join(mallet_dir, "best_model.pkl"), "wb") as f:
    pickle.dump(cobj, f, protocol=pickle.HIGHEST_PROTOCOL)

# --- 3) Binarize topics -> region sets --------------------------------------
bin_otsu = binarize_topics(cobj, method="otsu", plot=False)
bin_topk = binarize_topics(cobj, method="ntop", ntop=3000, plot=False)

for topic, df in bin_otsu.items():
    region_names_to_coordinates(df.index).sort_values(["Chromosome", "Start", "End"]).to_csv(
        os.path.join(rs_dir, "Topics_otsu", f"{topic}.bed"), sep="\t", header=False, index=False)

for topic, df in bin_topk.items():
    region_names_to_coordinates(df.index).sort_values(["Chromosome", "Start", "End"]).to_csv(
        os.path.join(rs_dir, "Topics_top3k", f"{topic}.bed"), sep="\t", header=False, index=False)

print("Region sets written to:", rs_dir)
