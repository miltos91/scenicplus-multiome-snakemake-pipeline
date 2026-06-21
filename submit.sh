#!/bin/bash
# Launches the full SCENICplus pipeline.
# Edit config.yaml for all parameters first, then: bash submit.sh

# Read parameters from config.yaml so they stay in sync
PARTITION=$(python3 -c "import yaml; c=yaml.safe_load(open('config.yaml')); print(c['slurm_partition'])")
MEM=$(python3       -c "import yaml; c=yaml.safe_load(open('config.yaml')); print(c['slurm_mem'])")
TIME=$(python3      -c "import yaml; c=yaml.safe_load(open('config.yaml')); print(c['slurm_time'])")
NCPU=$(python3      -c "import yaml; c=yaml.safe_load(open('config.yaml')); print(c['n_cpu'])")
MAIL=$(python3      -c "import yaml; c=yaml.safe_load(open('config.yaml')); print(c['slurm_mail'])")
RUN=$(python3       -c "import yaml; c=yaml.safe_load(open('config.yaml')); print(c['run_name'])")
ENVPY=$(python3     -c "import yaml; c=yaml.safe_load(open('config.yaml')); print(c['conda_env_python'])")
CONDA_SH=$(python3  -c "import yaml, os; c=yaml.safe_load(open('config.yaml')); print(os.path.expanduser(c.get('conda_sh', '~/miniforge3/etc/profile.d/conda.sh')))")
USE_SLURM=$(python3 -c "import yaml; c=yaml.safe_load(open('config.yaml')); print(str(c.get('use_slurm', True)).lower())")

if [ "$USE_SLURM" != "true" ]; then
    echo "use_slurm: false — running locally on $(hostname)"
    source "${CONDA_SH}"
    conda activate "${ENVPY}"
    export PYTHONUNBUFFERED=1 
    set -euo pipefail
    snakemake \
        --snakefile     Snakefile \
        --configfile    config.yaml \
        --cores         ${NCPU} \
        --rerun-incomplete \
        --printshellcmds
    echo "Pipeline complete: $(date)"
    exit 0
fi

# ── HPC run
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=SCENICplus_${RUN}
#SBATCH --partition=${PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${NCPU}
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME}
#SBATCH --output=SCENICplus_${RUN}_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${MAIL}

echo "Job started: \$(date)"
echo "Node: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"

source "${CONDA_SH}"
conda activate ${ENVPY}
export PYTHONUNBUFFERED=1   # stream python prints live into the SLURM log

set -euo pipefail

cd "\$(dirname "\$0")"

snakemake \
    --snakefile     Snakefile \
    --configfile    config.yaml \
    --cores         \$SLURM_CPUS_PER_TASK \
    --rerun-incomplete \
    --printshellcmds

echo "Pipeline complete: \$(date)"
EOF
