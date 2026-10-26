#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# HGPS cross-dataset meta-analysis + IQM QAOA benchmark — run every analysis script in pipeline order.
#
#   bash run_all.sh                 # all stages: 01 01test 02 03 04
#   bash run_all.sh 02 03           # selected stages only
#   QUICK=1 bash run_all.sh         # shorter QAOA optimisations (verification only)
#
# Stages
#   01      raw-RNA-seq downstream analyses (DESeq2, rMATS tables, enrichment, signatures,
#           pathway/TF, isoform switch, WGCNA) from the included count/rMATS/salmon tables
#   01test  raw-data steps (download → QC → trim → STAR → RSeQC → featureCounts → salmon →
#           rMATS → dupRadar) on the tiny test subsets; writes only to $WORK_ROOT
#   02      per-dataset DE of GEO-processed data (GSE231761, GSE306264)
#   03      k=3 human random-effects meta-analysis, mouse comparison, GSEA,
#           GSE231761 donor-aware sensitivity analysis (03_meta_analysis)
#   04      QUBO construction, classical benchmark, QAOA simulation, IQM hardware runs/scoring
#
# Outputs are written next to the scripts (as in the original project), i.e. INTO this
# directory tree. To keep the deposit pristine, run on a copy:
#   cp -a HGPS_meta_QAOA /work/hgps_run && bash /work/hgps_run/run_all.sh
#
# Hardware: credential variables are UNSET by default, so every hardware script uses its
# local-simulator fallback (outputs carry _SIMULATOR_FALLBACK). Set USE_HARDWARE=1 (and the
# IQM variables in .env.example) to submit to the real IQM devices (costs money).
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
DEP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DEP"
export PROGERIA_DEPOSIT="$DEP"

CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
export CONDA_ENVS="${CONDA_ENVS:-$CONDA_BASE/envs}"
RS="${RS:-$CONDA_ENVS/lmna_r/bin/Rscript}"                  # R 4.5.3   (envs/lmna_r.environment.yml)
PY_IQM="${PY_IQM:-$CONDA_ENVS/iqm_env/bin/python}"           # py 3.11   (envs/iqm_env.environment.yml)
WORK_ROOT="${WORK_ROOT:-$(mktemp -d "${TMPDIR:-/tmp}/progeria_test.XXXXXX")}"

export PYTHONHASHSEED=0          # networkx subgraph node order (= qubit↔gene order) depends on it
export MPLBACKEND=Agg
export STRING_CACHE="${STRING_CACHE:-1}"   # reuse the shipped STRING tables (live STRING has drifted)
if [[ "${USE_HARDWARE:-0}" != "1" ]]; then
    unset IQM_TOKEN IQM_SERVER_URL 2>/dev/null || true
fi

log() { echo "[$(date '+%F %T')] $*"; }
run() { log "▶ $*"; local t0=$SECONDS; "$@"; log "✔ ($((SECONDS - t0)) s) $1 ${2:-}"; }

G1=01_rnaseq_processing/GSE165409_mouse_aorta
G2=01_rnaseq_processing/GSE198983_human_AoSMC
M=04_quantum_qaoa/meta_QUBO_IQM/scripts

stage_01() {   # downstream analyses of the two re-processed raw datasets (full, from included tables)
    run "$RS" $G1/scripts/08_deseq2.R
    run "$RS" $G1/scripts/09_rmats_figures.R        # writes significant_events_*.tsv used by 12
    run "$RS" $G1/scripts/12_enrichment.R
    run "$RS" $G1/scripts/11_isoformswitch.R        # ~10–60 min
    ( set -a; source $G2/scripts/00_config.sh; set +a
      run "$RS" $G2/scripts/09_deseq2.R
      run "$RS" $G2/scripts/10_enrichment.R         # network: KEGG REST
      run "$RS" $G2/scripts/11_signatures.R
      run "$RS" $G2/scripts/12_pathway_tf.R
      run "$RS" $G2/scripts/16_isoformswitch.R      # ~10–60 min
      run "$RS" $G2/scripts/17_wgcna.R )
    DRY_RUN=1 run "$RS" $G2/scripts/install_r_pkgs.R  # report-only (never installs in DRY_RUN mode)
}

stage_01test() {   # raw-data scripts on the tiny test subsets (outputs only in $WORK_ROOT)
    log "test work dirs under $WORK_ROOT"
    THREADS="${THREADS:-4}" RUN_DOWNSTREAM=1 run bash $G1/test/run_test.sh "$WORK_ROOT/GSE165409"
    PROJ="$WORK_ROOT/GSE165409" TINY_TEST=1 THREADS="${THREADS:-4}" run bash $G1/scripts/run_remaining.sh
    PROJ="$WORK_ROOT/GSE165409" TINY_TEST=1 THREADS="${THREADS:-4}" run bash $G1/scripts/run_remaining2.sh
    WORK="$WORK_ROOT/GSE198983" THREADS="${THREADS:-4}" run bash $G2/test/run_test.sh
}

stage_02() {   # per-dataset DE (GEO processed data)
    run "$RS" 02_per_dataset_DE/GSE231761/scripts/01_deseq2_GSE231761.R
    run "$RS" 02_per_dataset_DE/GSE231761/scripts/02_enrichment_activity.R
    run "$RS" 02_per_dataset_DE/GSE231761/scripts/03_passage_positional_cnv.R
    run "$RS" 02_per_dataset_DE/GSE306264/scripts/01_deseq2_GSE306264.R
    run "$RS" 02_per_dataset_DE/GSE306264/scripts/02_enrichment.R
    run "$RS" 02_per_dataset_DE/GSE306264/scripts/03_activity_heterogeneity.R
    run "$RS" 02_per_dataset_DE/GSE306264/scripts/04_positional_LAD.R
    run "$RS" 02_per_dataset_DE/GSE306264/scripts/05_wgcna.R
}

stage_03() {   # meta-analysis (06 must precede 10: it writes GSE198983 19_positional table)
    for s in 01_ortholog_and_preliminary_meta 02_meta_unshrunken 05_meta_k3 06_k3_figures \
             07_candidate_biology 08_common_pathways 10_fourway_common 21_core287_GO_ORA \
             22_meta_GSEA 27_gene_annotation 28_sensitivity_GSE231761_donor; do
        run "$RS" 03_meta_analysis/scripts/$s.R
    done
}

stage_04() {   # meta-analysis genes → QUBO → classical / QAOA simulation → IQM hardware (fallback)
    run "$RS" $M/00_export_vst.R
    run "$PY_IQM" $M/01_build_network_qubo.py
    run "$PY_IQM" $M/02_build_coexpression_qubo.py
    run "$PY_IQM" $M/03_classical_benchmark.py
    run "$PY_IQM" $M/04_qaoa_simulator.py            # ~30–90 min (QUICK=1: shorter)
    run "$PY_IQM" $M/05_verify_p3_and_figures.py
    run "$PY_IQM" $M/06_iqm_device_probe.py
    run "$PY_IQM" $M/07_iqm_plan_and_predict.py      # long (QUICK=1: ~40 min)
    run "$PY_IQM" $M/08_iqm_submit.py --dry-run
    run "$PY_IQM" $M/08_iqm_submit.py --yes
    run "$PY_IQM" $M/10_iqm_connect_test.py
    run "$PY_IQM" $M/11_qrisp_qaoa_iqm.py --mode mock --instance coexp_ising_16 --depth 1
    run "$PY_IQM" $M/11_qrisp_qaoa_iqm.py --mode hw --instance coexp_ising_8    --depth 1 --repeat 5 --yes
    run "$PY_IQM" $M/11_qrisp_qaoa_iqm.py --mode hw --instance coexp_ising_16   --depth 1 --repeat 5 --yes
    run "$PY_IQM" $M/11_qrisp_qaoa_iqm.py --mode hw --instance coexp_ising_16   --depth 2 --repeat 5 --yes
    run "$PY_IQM" $M/11_qrisp_qaoa_iqm.py --mode hw --instance string_maxcut_19 --depth 1 --repeat 5 --yes
    run "$PY_IQM" $M/11_qrisp_qaoa_iqm.py --mode hw --instance string_maxcut_19 --depth 2 --repeat 5 --yes
    run "$PY_IQM" $M/11_qrisp_qaoa_iqm.py --mode hw --instance coexp_ising_34   --depth 1 --device emerald --theta 1.7171,1.0954 --yes
    run "$PY_IQM" $M/18_gate_counts.py
    run "$PY_IQM" $M/13_benchmark_agreement.py
    run "$PY_IQM" $M/14_collect_iqm_results.py
    run "$PY_IQM" $M/16_fig4_benchmark_and_partition_stability.py   # Fig 4E-H values, per-gene stability
    run "$PY_IQM" $M/17_energy_landscape.py
    # identical copies archived with the hardware results:
    run "$PY_IQM" 04_quantum_qaoa/IQM_hardware_results/06_scripts/10_iqm_connect_test.py
    run "$PY_IQM" 04_quantum_qaoa/IQM_hardware_results/06_scripts/13_benchmark_agreement.py
    # NOT run: 11_qrisp_qaoa_iqm.py --mode local (re-trains angles and overwrites the shipped
    #          qrisp_angles_*.json that were used on the real hardware).
}

STAGES=("$@"); [[ ${#STAGES[@]} -gt 0 ]] || STAGES=(01 01test 02 03 04)
for s in "${STAGES[@]}"; do
    log "════════ stage $s ════════"
    "stage_$s"
done
log "════════ all requested stages finished ════════"
