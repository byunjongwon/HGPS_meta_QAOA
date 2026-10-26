# VERIFICATION_LOG — HGPS meta-analysis and IQM QAOA benchmark code

**Date**: 2026-10-19 – 2026-10-22. **Platform**: Linux 6.11 aarch64, 20 cores, 119 GB RAM (shared, heavily loaded by
other jobs during parts of the runs, so wall times are inflated 2–5×).
**Credentials**: none. Every command ran with `IQM_TOKEN`, `IQM_SERVER_URL` (and all other vendor tokens) unset,
`PYTHONHASHSEED=0`, `MPLBACKEND=Agg`, `STRING_CACHE=1`.

Interpreters: `R` = envs/lmna_r Rscript (R 4.5.3);
`PY_IQM` = envs/iqm_env python 3.11.15. Bioinformatics tools: see `envs/TOOL_VERSIONS.tsv`.

Verification levels: **full** = the script ran end-to-end on the real (included) inputs;
**full (QUICK=1)** = full run with the optional reduced-iteration switch (the default-parameter run was also done, see column);
**test-subset** = raw-data step run end-to-end on the included tiny test subset (real input is 100+ GB);
**fallback** = IQM hardware path replaced by a local IQM fake-backend / Aer simulation because no token was set;
**dry-run** = report-only mode.

## A. Integrated pipeline run (final deposit, pipeline order)

Scratch copy of the deposit, driven by `run_all.sh` (stages run as `01`, `01test`, `02 03 04` with `QUICK=1`;
the chains ran concurrently). This run predates the public release; rows for code that is not part of the
release (in-house MEF data set, GSE137083) have been removed from the table below — see section E.
Start 2026-10-22 14:30, last stage finished 16:12. Every step returned exit code 0.
After the run, all shipped IQM raw-hardware files (`IQM_hardware_results/02_raw_hardware/*`, `meta_QUBO_IQM/results/qrisp_hw_*`,
`qrisp_angles_*`) were compared with `cmp`: **unchanged** (fallback outputs went to `*_SIMULATOR_FALLBACK*`).

| # | script (path from deposit root) | command / args | env | runtime | status | level |
|---|---|---|---|---|---|---|
| 01-1 | 01_rnaseq_processing/GSE165409_mouse_aorta/scripts/08_deseq2.R | `$R …` | lmna_r | 25 s | PASS | full |
| 01-2 | …/GSE165409_mouse_aorta/scripts/09_rmats_figures.R | `$R …` | lmna_r | 9 s | PASS | full |
| 01-3 | …/GSE165409_mouse_aorta/scripts/12_enrichment.R | `$R …` | lmna_r | 95 s | PASS | full |
| 01-4 | …/GSE165409_mouse_aorta/scripts/11_isoformswitch.R | `$R …` | lmna_r | 721 s | PASS | full |
| 01-5 | 01_rnaseq_processing/GSE198983_human_AoSMC/scripts/09_deseq2.R | `set -a; source scripts/00_config.sh; $R …` | lmna_r | 151 s | PASS | full |
| 01-6 | …/GSE198983_human_AoSMC/scripts/10_enrichment.R | same (network: KEGG REST) | lmna_r | 967 s | PASS | full |
| 01-7 | …/GSE198983_human_AoSMC/scripts/11_signatures.R | same | lmna_r | 26 s | PASS | full |
| 01-8 | …/GSE198983_human_AoSMC/scripts/12_pathway_tf.R | same | lmna_r | 55 s | PASS | full |
| 01-9 | …/GSE198983_human_AoSMC/scripts/16_isoformswitch.R | same | lmna_r | 1592 s | PASS | full |
| 01-10 | …/GSE198983_human_AoSMC/scripts/17_wgcna.R | same | lmna_r | 50 s | PASS | full |
| 01-11 | …/GSE198983_human_AoSMC/scripts/install_r_pkgs.R | `DRY_RUN=1 $R …` | lmna_r | 5 s | PASS | dry-run (never installs; would modify the env) |
| T-1 | 01_rnaseq_processing/GSE165409_mouse_aorta/test/run_test.sh | `THREADS=4 RUN_DOWNSTREAM=1 bash … $WORK/GSE165409` → build_test_reference.sh (STAR index, 03b_make_bed12.R), 01_fastqc_raw, 02_trim, 03_star_align, 04_post_align_qc, 05_featurecounts (+author-count comparison), 06_salmon, 07_rmats, 13_dupradar, 08_deseq2 (TINY_TEST), 09 | lmna_rmats / rnaseq / lmna_tools / lmna_r | 168 s | PASS | test-subset |
| T-2 | …/GSE165409_mouse_aorta/scripts/run_remaining.sh | `PROJ=$WORK/GSE165409 TINY_TEST=1 bash …` (04→05→06→07→08) | bash | 14 s | PASS | test-subset |
| T-3 | …/GSE165409_mouse_aorta/scripts/run_remaining2.sh | same (08→09→12→11→13) | bash | 61 s | PASS | test-subset |
| T-4 | 01_rnaseq_processing/GSE198983_human_AoSMC/test/run_test.sh | `WORK=$WORK/GSE198983 THREADS=4 bash …` → 01_download_fastq (tiny public ENA run SRR19598819: GET+md5, then SKIP), 02_reference (mini genome, STAR + salmon decoy index), 02b_make_bed12, 03_fastqc_raw, 04_trim, 05_star_align, 06_post_align_qc (+06b_dupradar.R, Picard, Qualimap), 07_featurecounts, 08_salmon, 13_rmats, run_pipeline.sh ONLY=07 | lmna_rmats / rnaseq / lmna_tools / lmna_r | 184 s | PASS | test-subset |
| 02-4 | …/GSE231761/scripts/01_deseq2_GSE231761.R | `$R …` | lmna_r | 20 s | PASS | full |
| 02-5 | …/GSE231761/scripts/02_enrichment_activity.R | `$R …` | lmna_r | 48 s | PASS | full |
| 02-6 | …/GSE231761/scripts/03_passage_positional_cnv.R | `$R …` | lmna_r | 23 s | PASS | full |
| 02-7 | …/GSE306264/scripts/01_deseq2_GSE306264.R | `$R …` | lmna_r | 37 s | PASS | full |
| 02-8 | …/GSE306264/scripts/02_enrichment.R | `$R …` | lmna_r | 82 s | PASS | full |
| 02-9 | …/GSE306264/scripts/03_activity_heterogeneity.R | `$R …` | lmna_r | 23 s | PASS | full |
| 02-10 | …/GSE306264/scripts/04_positional_LAD.R | `$R …` | lmna_r | 32 s | PASS | full |
| 02-11 | …/GSE306264/scripts/05_wgcna.R | `$R …` | lmna_r | 58 s | PASS | full |
| 03-1 | 03_meta_analysis/scripts/01_ortholog_and_preliminary_meta.R | `$R …` | lmna_r | 49 s | PASS | full |
| 03-2 | …/02_meta_unshrunken.R | `$R …` | lmna_r | 6 s | PASS | full |
| 03-3 | …/05_meta_k3.R | `$R …` | lmna_r | 41 s | PASS | full |
| 03-4 | …/06_k3_figures.R (computes GSE198983 positional/LAD table) | `$R …` | lmna_r | 18 s | PASS | full |
| 03-5 | …/07_candidate_biology.R | `$R …` | lmna_r | 380 s | PASS | full |
| 03-6 | …/08_common_pathways.R | `$R …` | lmna_r | 5 s | PASS | full |
| 03-7 | …/10_fourway_common.R | `$R …` | lmna_r | 49 s | PASS | full |
| 03-8 | …/21_core287_GO_ORA.R | `$R …` | lmna_r | 110 s | PASS | full |
| 03-9 | …/22_meta_GSEA.R | `$R …` | lmna_r | 96 s | PASS | full |
| 03-10 | …/27_gene_annotation.R | `$R …` | lmna_r | 6 s | PASS | full |
| 04-0 | 04_quantum_qaoa/meta_QUBO_IQM/scripts/00_export_vst.R | `$R …` | lmna_r | 12 s | PASS | full |
| 04-1 | …/01_build_network_qubo.py | `STRING_CACHE=1 $PY_IQM …` | iqm_env | 4 s | PASS | full |
| 04-2 | …/02_build_coexpression_qubo.py | `$PY_IQM …` | iqm_env | 1 s | PASS | full |
| 04-3 | …/03_classical_benchmark.py | `$PY_IQM …` | iqm_env | 58 s | PASS | full |
| 04-4 | …/04_qaoa_simulator.py | `QUICK=1 $PY_IQM …` | iqm_env | 2269 s | PASS | full (QUICK=1); default run: see B |
| 04-5 | …/05_verify_p3_and_figures.py | `$PY_IQM …` | iqm_env | 397 s | PASS | full |
| 04-6 | …/06_iqm_device_probe.py | `$PY_IQM …` | iqm_env | 4 s | PASS | fallback (IQMFakeAphrodite) |
| 04-7 | …/07_iqm_plan_and_predict.py | `QUICK=1 $PY_IQM …` | iqm_env | 263 s | PASS | full (QUICK=1); default run: see B |
| 04-8 | …/08_iqm_submit.py | `--dry-run` | iqm_env | 312 s | PASS | dry-run |
| 04-9 | …/08_iqm_submit.py | `--yes` | iqm_env | 321 s | PASS | fallback (IQM fake backends) |
| 04-10 | …/10_iqm_connect_test.py | `$PY_IQM …` | iqm_env | 3 s | PASS | fallback (IQMFakeApollo) |
| 04-11 | …/11_qrisp_qaoa_iqm.py | `--mode mock --instance coexp_ising_16 --depth 1` | iqm_env | 66 s | PASS | fallback |
| 04-12 | …/11_qrisp_qaoa_iqm.py | `--mode hw --instance coexp_ising_8 --depth 1 --repeat 5 --yes` | iqm_env | 7 s | PASS | fallback |
| 04-13 | …/11_qrisp_qaoa_iqm.py | `--mode hw --instance coexp_ising_16 --depth 1 --repeat 5 --yes` | iqm_env | 171 s | PASS | fallback |
| 04-14 | …/11_qrisp_qaoa_iqm.py | `--mode hw --instance coexp_ising_16 --depth 2 --repeat 5 --yes` | iqm_env | 367 s | PASS | fallback |
| 04-15 | …/11_qrisp_qaoa_iqm.py | `--mode hw --instance string_maxcut_19 --depth 1 --repeat 5 --yes` | iqm_env | 314 s | PASS | fallback |
| 04-16 | …/11_qrisp_qaoa_iqm.py | `--mode hw --instance string_maxcut_19 --depth 2 --repeat 5 --yes` | iqm_env | 404 s | PASS | fallback |
| 04-17 | …/11_qrisp_qaoa_iqm.py | `--mode hw --instance coexp_ising_34 --depth 1 --device emerald --theta 1.7171,1.0954 --yes` | iqm_env | 13 s | PASS | fallback (approximate MPS, bond dim 64; exact 34-q simulation needs ~4.5 TB) |
| 04-18 | …/18_gate_counts.py | `$PY_IQM …` | iqm_env | 4 s | PASS | full |
| 04-19 | …/13_benchmark_agreement.py (scores the real IQM raw counts) | `$PY_IQM …` | iqm_env | 3 s | PASS | full |
| 04-20 | …/14_collect_iqm_results.py | `$PY_IQM …` | iqm_env | 0 s | PASS | full |
| 04-21 | …/17_energy_landscape.py | `$PY_IQM …` (uses shipped landscape cache) | iqm_env | 1 s | PASS | full (no-cache recomputation: see B) |
| 04-22 | 04_quantum_qaoa/IQM_hardware_results/06_scripts/10_iqm_connect_test.py | `$PY_IQM …` | iqm_env | 4 s | PASS | fallback |
| 04-23 | 04_quantum_qaoa/IQM_hardware_results/06_scripts/13_benchmark_agreement.py | `$PY_IQM …` | iqm_env | 2 s | PASS | full |

## B. Section verification runs (done by section, before/after the scope changes)

Each section was additionally verified in its own scratch copy; results relevant beyond section A:

| script | what was additionally verified | runtime | result |
|---|---|---|---|
| 04 `04_qaoa_simulator.py` | default parameters (MAXITER 300) | 87 min (loaded machine) | PASS; `qaoa_noiseless.tsv`, `qaoa_noise_scan.tsv`, `resource_estimate.tsv` identical to the shipped files |
| 04 `07_iqm_plan_and_predict.py` | default parameters (8 restarts × 400) | 3.6 h (loaded) | PASS; noiseless ratios and `iqm_run_plan.json` identical; fake-backend noisy predictions differ in the 3rd decimal (unseeded) |
| 04 `17_energy_landscape.py` | without the shipped cache | 17 min | PASS; `.npz` identical to shipped |
| 04 `01_build_network_qubo.py` | default live-STRING mode | – | PASS (runs), but STRING has drifted (panelC 18 nodes instead of 19) → run_all uses `STRING_CACHE=1` |
| 04 `11_qrisp_qaoa_iqm.py` | `--mode local` (re-trains angles) | – | PASS; excluded from run_all because it overwrites the angle files used on the real hardware |
| 04 06_scripts/11 copy | `--mode hw --instance coexp_ising_8` from the IQM_hardware_results location | – | PASS (fallback); 06_scripts/{10,11,13} are byte-identical to scripts/ |
| 02/03 all scripts | regenerated tables/RDS compared with originals | ~25 min | all byte-identical (except `vst_for_quantum.tsv`, written by 04/00) |
| 01 GSE198983 09/10/11/12/16/17 | compared with originals | – | DE/vst/WGCNA/DTU/signature tables byte-identical; KEGG tables differ (live KEGG REST) |
| 01 GSE165409 08/09/11/12 | compared with originals | 11: 59 min | DE top genes/padj identical; switching gene counts identical |
| 01 GSE165409 test steps 01–07, 13 individually | each script alone with `PROJ=<work>` | 3–202 s each | PASS (test-subset); rMATS finds the Lmna exon-11 A5SS (progerin) ΔPSI −0.825 |
| 01 GSE165409 `run_remaining3/4/5.sh` | resume subsets of run_remaining2 on the test work dir | 2–8 min | PASS (test-subset) |
| 01 GSE198983 `make_test_subset.sh`, GSE165409 `make_test_subset.sh` | rebuilt the test subsets from the original BAMs/reference (needs `SRC_PROJ` = original full-size run, not shipped) | 19 s / 75 s | PASS |
| all `.sh` | `bash -n` | – | OK |

## C. Not fully verified (and why)

- **Raw-data steps on the real data** (download of ~100 GB FASTQ, full GRCh38/GRCm38 index builds, STAR/salmon/rMATS on all samples):
  verified on the tiny test subsets only (level test-subset). The full pipeline was run originally on the same scripts.
- **Real IQM hardware submission** (06, 08 `--yes`, 10, 11 hw/mock): not executed (no credentials, by design);
  the fallback paths were executed. The real-hardware code paths are unchanged apart from minimal credential handling.
- **`install_r_pkgs.R`** install mode not executed (would modify the R environment); `DRY_RUN=1` mode verified.
- **QUICK=1** was used for 04 and 07 in section A; the default-parameter runs were done in section B.
- Stochastic parts (IQM fake-backend sampling, KEGG live API) reproduce only approximately.

## D. Not included

Figure-assembly / presentation code, code for quantum platforms other than IQM, the in-house MEF data set,
GSE137083 (screened and rejected) and the GENCODE reference files (`fetch_references.sh`) are not part of this
repository.

## E. Release verification (public repository, 2026-10-23)

Run on a scratch copy of the release tree with `bash run_all.sh 02 03` (R 4.5.3, `lmna_r`), after removal of the
MEF data set and GSE137083, the `run_all.sh` / `06_k3_figures.R` edits and the addition of
`28_sensitivity_GSE231761_donor.R`. Every step exited 0. All result files of stages 02 and 03 that are shipped
in the repository were regenerated **byte-identical** (`diff -rq`); the only new outputs were the three
`sensitivity_GSE231761_*.tsv` files, which were then added to the release.

| # | script | runtime | status | level |
|---|---|---|---|---|
| R-02 | stage 02 (GSE231761 01–03, GSE306264 01–05) | 5 min | PASS, outputs identical | full |
| R-03 | stage 03 (01, 02, 05, 06, 07, 08, 10, 21, 22, 27) | 5 min | PASS, outputs identical | full |
| R-28 | `03_meta_analysis/scripts/28_sensitivity_GSE231761_donor.R` | 4 min | PASS; internal checks: primary model reproduces `meta_k3_final.tsv` (324 / 287 / 97 / 34) and `GSEA_meta_Hallmark.tsv` (10 Hallmark sets FDR<0.05; 106 nonredundant GO:BP terms) | full |

| R-16 | `04_quantum_qaoa/meta_QUBO_IQM/scripts/16_fig4_benchmark_and_partition_stability.py` (added 2026-10-26; deposit-relative paths and a per-gene stability table were the only changes to the original analysis script) | 1 min | PASS; reproduces the reported per-gene stability: co-expression 10/16 genes stable in all ten hardware runs, mean 93.8%; STRING 4/19, mean 81.6%, *WNT5A* 50% | full |

| R-12 | `01_rnaseq_processing/GSE165409_mouse_aorta/scripts/09_rmats_figures.R` → `12_enrichment.R` (re-run 2026-10-26) | 2 min | PASS; ORA tables byte-identical; GSEA tables (unseeded) within ±0.02 NES of the shipped originals, identical significance for the six recurrent pathways in all three contrasts. The original `results/12_enrichment/*.tsv` used for Fig. 3C were added to the release. | full |

Apart from script 16, stages 01 and 04 contain no code changed for the release and were verified as described in sections A–B.
