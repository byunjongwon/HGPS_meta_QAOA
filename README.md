# HGPS cross-dataset meta-analysis and IQM QAOA benchmark

Analysis code, intermediate results and raw quantum-hardware measurement counts for

> Jeong Y\*, Byun J\*, Kim Y, Moon S-H, Kwon T. **Cross-Dataset Transcriptomic Analysis of Shared HGPS
> Features and Exploratory Quantum Network Benchmarking.** (manuscript under review; citation will be
> updated on publication). \*Equal contribution.

The repository reproduces the complete analytical workflow of the manuscript:

1. dataset-specific differential expression of four public bulk RNA-seq datasets (three human, one mouse);
2. gene-wise random-effects (DerSimonian–Laird) meta-analysis of the three human datasets, including the
   donor-aware sensitivity analysis for GSE231761;
3. Hallmark / GO enrichment per dataset and on the integrated ranking, and cross-species comparison with
   *Lmna*<sup>G609G</sup> mouse aorta;
4. construction of the STRING Max-Cut and signed co-expression Ising QUBO instances, exact enumeration and
   classical heuristics, QAOA parameter optimisation, execution on IQM Garnet, and scoring of the raw
   hardware counts against exact optima and equal-budget random sampling.

Archived version: Zenodo DOI `10.5281/zenodo.23660505` (added after the first release).

---

## 1. Manuscript ↔ code map

Paths are relative to the repository root. `results/` folders contain the outputs as produced for the
manuscript; re-running a script regenerates them in place.

| Manuscript item | Script(s) | Main output(s) |
|---|---|---|
| Sample structure, PCA / correlation / distance heatmaps (Fig. 1B–G) | `02_per_dataset_DE/GSE306264/scripts/01_deseq2_GSE306264.R`, `02_per_dataset_DE/GSE231761/scripts/01_deseq2_GSE231761.R`, `01_rnaseq_processing/GSE198983_human_AoSMC/scripts/09_deseq2.R`, `01_rnaseq_processing/GSE165409_mouse_aorta/scripts/08_deseq2.R` | `results/` of each dataset (VST objects, PCA tables / diagnostic plots) |
| Raw-read reprocessing of GSE165409 and GSE198983 (FASTQ → fastp → STAR → featureCounts) | `01_rnaseq_processing/*/scripts/01…07` (+ `test/` mini data set) | `01_rnaseq_processing/*/results/05_counts`, `07_counts` |
| Per-dataset DESeq2 (Fig. 2A–D) | same four DE scripts | `DE_*.tsv` |
| GSE231761 model (`~ passage + condition`, technical duplicates collapsed) and exclusion of *BMPR1AP1*, *CA8*, *LINC02899*, *NPIPB2* | `02_per_dataset_DE/GSE231761/scripts/01_deseq2_GSE231761.R`, `03_passage_positional_cnv.R` | `GSE231761_objects.rds`, `CNV_polymorphism_candidates.tsv` |
| GSE198983 inferred two-level factor (factor A/B) as blocking covariate | `01_rnaseq_processing/GSE198983_human_AoSMC/scripts/09_deseq2.R` | `results/09_deseq2/dds.rds`, `DE_Progerin_vs_Control_blocked.tsv` |
| Meta-analysis: 12,290 genes → 324 (FDR<0.05) → 287 concordant (same direction, I²<75%) → 97 (\|log2FC\|≥0.5) → 34 (≥1); 16 / 7 cross-species genes; pairwise Spearman (Fig. 2E, F, H) | `03_meta_analysis/scripts/05_meta_k3.R` | `meta_k3_final.tsv`, `candidates_k3_final.tsv`, `candidates_k3_strict.tsv` |
| **Donor-aware sensitivity analysis for GSE231761** (revision) | `03_meta_analysis/scripts/28_sensitivity_GSE231761_donor.R` | `sensitivity_GSE231761_summary.tsv`, `_hallmark.tsv`, `_genes.tsv` |
| Integrated Hallmark and GO:BP GSEA (Fig. 2G) | `03_meta_analysis/scripts/22_meta_GSEA.R` | `GSEA_meta_Hallmark.tsv`, `GSEA_meta_GOBP.tsv` |
| Per-dataset Hallmark GSEA and recurrent pathways across the four contexts (Fig. 3A, B) | `02_per_dataset_DE/*/scripts/02_enrichment*.R`, `03_meta_analysis/scripts/08_common_pathways.R`, `10_fourway_common.R` | `TABLE2_hallmark_across_datasets.tsv`, `TABLE3_hallmark_4way.tsv` |
| Mouse contrasts *Lmna*<sup>G609G</sup> vs young WT, old vs young WT, *Lmna*<sup>G609G</sup> vs old WT (Fig. 3C) | `01_rnaseq_processing/GSE165409_mouse_aorta/scripts/08_deseq2.R`, `12_enrichment.R` | `results/08_deseq2/DE_*.tsv`, `results/12_enrichment/GSEA_hallmark_<contrast>.tsv` |
| Mouse aorta column of the four-context Hallmark comparison (Fig. 3A, B) | `03_meta_analysis/scripts/10_fourway_common.R` | `GSEA_GSE165409_Hallmark_full.tsv`, `TABLE3_hallmark_4way.tsv` |
| Contractile markers in AOSMC (Fig. 3D) | `01_rnaseq_processing/GSE198983_human_AoSMC/scripts/09_deseq2.R`, `11_signatures.R` | `results/09_deseq2`, `results/11_signatures` |
| Gene annotation for Table I | `03_meta_analysis/scripts/27_gene_annotation.R` | `qubo_gene_annotation.tsv` |
| STRING network (19 genes, 26 edges) → weighted Max-Cut QUBO (Fig. 4A, B) | `04_quantum_qaoa/meta_QUBO_IQM/scripts/01_build_network_qubo.py` (uses the shipped STRING query, see §5) | `string_maxcut_19_*`, `panel*_string_edges.tsv` |
| Signed co-expression network (16 genes, 47 edges) → Ising / QUBO (Fig. 4A, B) | `00_export_vst.R`, `02_build_coexpression_qubo.py` | `coexp_ising_16_*` |
| Exact enumeration, greedy, simulated annealing, genetic algorithm (Fig. 4D, E) | `03_classical_benchmark.py` | `classical_benchmark.tsv`, `*_classical.json` |
| Noiseless p=1 QAOA landscapes (Fig. 4C) | `17_energy_landscape.py` (+ `04_qaoa_simulator.py`) | `*_landscape_*.npz`, `qaoa_noiseless.tsv` |
| QAOA angle optimisation (statevector) and fixed-angle execution on IQM Garnet, 5 runs × 4,096 shots, p=1/2 | `11_qrisp_qaoa_iqm.py` (`--mode local` for angles, `--mode hw` for hardware) | `qrisp_angles_*.json`, `qrisp_hw_*.json`; raw counts in `04_quantum_qaoa/IQM_hardware_results/02_raw_hardware/` |
| Normalised solution quality, equal-budget random baseline, one-sample t-tests, mean partition agreement with the exact solution (Fig. 4D, F) | `13_benchmark_agreement.py`, `14_collect_iqm_results.py`, `18_gate_counts.py` | `benchmark_table.txt`, `benchmark_agreement.json`, `hardware_summary_table.json`, `gate_counts.json` |
| Runtime per solve (Fig. 4E), pairwise partition agreement (Fig. 4F), per-gene partition stability across the ten hardware runs, e.g. 10/16 and 4/19 genes stable, *WNT5A* 50% (Fig. 4G, H) | `16_fig4_benchmark_and_partition_stability.py` | `partition_stability_coexp_ising_16.tsv`, `partition_stability_string_maxcut_19.tsv`; panels in `figures/fig4/` |

Scripts not listed (e.g. isoform switching, rMATS, WGCNA, positional/LAD and transcription-factor
activity analyses, `07_candidate_biology.R`, `21_core287_GO_ORA.R`) are exploratory analyses run on the same
data; they are kept for completeness and are not required for any number reported in the manuscript.
Figure-assembly code is not included, with one exception: `16_fig4_benchmark_and_partition_stability.py`
draws Figure 4 and also computes the runtime, pairwise-agreement and per-gene stability values reported in
the text. Its runtimes are wall-clock measurements and differ slightly between executions.

## 2. Repository layout

```
HGPS_meta_QAOA/
├── README.md  LICENSE  CITATION.cff  run_all.sh  fetch_references.sh  VERIFICATION_LOG.md  .env.example
├── envs/                     conda environment exports (versions pinned), R sessionInfo, tool versions
├── 01_rnaseq_processing/
│   ├── GSE165409_mouse_aorta/       Lmna G609G vs WT 2 mo vs WT 24 mo aorta (FASTQ reprocessing + DE)
│   └── GSE198983_human_AoSMC/       FLAG-progerin vs FLAG aortic SMC (FASTQ reprocessing + DE)
│       each: scripts/  meta/ (sample sheets, ENA/SRA URLs, md5)  results/  test/ (tiny test data set)  ref/ (see fetch_references.sh)
├── 02_per_dataset_DE/
│   ├── GSE306264/                   HGPS (20) vs control (10) dermal fibroblast lines, GEO count table
│   └── GSE231761/                   HGPS (2) vs control (2) iPSC-VSMC lines × P7/P14, GEO count files
│       each: scripts/  meta/  raw_data/ (GEO supplementary count files)  results/
├── 03_meta_analysis/                scripts/ 01 … 28, results/
└── 04_quantum_qaoa/
    ├── meta_QUBO_IQM/               scripts/ 00 … 18, results/ (QUBO matrices, classical, simulation, hardware JSON)
    └── IQM_hardware_results/        02_raw_hardware (IQM raw counts), 04_settings (angles, run plan),
                                     05_reference (exact optima), 01_summary, 06_scripts (copies of 10/11/13)
```

## 3. Data sources

| Dataset | Content | Input used here |
|---|---|---|
| GSE306264 | primary dermal fibroblasts, 20 HGPS and 10 control donor lines, technical triplicates | GEO `GSE306264_featureCounts_all_gene_name.txt.gz` (`02_per_dataset_DE/GSE306264/raw_data`) |
| GSE231761 | iPSC-derived VSMCs, 2 HGPS and 2 control donor lines, passages 7 and 14, technical duplicates | GEO `GSE231761_RAW.tar` per-sample featureCounts files (`raw_data`) |
| GSE198983 (PRJNA817614) | human aortic SMCs, FLAG vs FLAG-progerin, 4 libraries per group | FASTQ from ENA (`meta/ena_filereport.tsv`), reprocessed; counts included |
| GSE165409 (SRP303066) | aorta of 2-month *Lmna*<sup>G609G</sup>, 2-month WT and 24-month WT mice, n = 6 per group | FASTQ from SRA/ENA (`meta/fastq_urls.txt`, `md5.txt`), reprocessed; counts included |
| STRING | protein association network (string-db.org API, version current at the time of the 2026 query), combined score ≥ 400 | query results cached in `04_quantum_qaoa/meta_QUBO_IQM/results/panel*_string_edges.tsv` |
| IQM Resonance | IQM Garnet (20 qubits) hardware runs | raw measurement counts in `04_quantum_qaoa/IQM_hardware_results/02_raw_hardware` |
| MSigDB Hallmark, GO | via `msigdbr` 26.1.0 and `org.Hs.eg.db` / `clusterProfiler` | versions in `envs/R_sessionInfo_lmna_r.txt` |

GSE137083 was screened and not retained (two donor lines with cell-line identity confounding); its analysis
is not part of the manuscript and is not included.

## 4. Software environment

Verified on Linux (aarch64, 20 cores, 119 GB RAM) with miniforge3. Exported environments are in `envs/`.

| conda env | used for | files |
|---|---|---|
| `lmna_r` (R 4.5.3; DESeq2 1.50.2, metafor, limma, edgeR, fgsea, msigdbr 26.1.0, clusterProfiler 4.18.4, WGCNA 1.74 …) | all R scripts | `lmna_r.environment.yml`, `R_sessionInfo_lmna_r.txt`, `R_installed_packages_lmna_r.tsv` |
| `iqm_env` (Python 3.11; qrisp 0.8.2, iqm-client 34.0.2, qiskit 2.4.1, qiskit-aer 0.17.2, networkx) | stage 04 | `iqm_env.environment.yml`, `iqm_env.pip_freeze.txt` |
| `lmna_rmats`, `rnaseq`, `lmna_tools` | STAR 2.7.11b, samtools, rMATS, FastQC, MultiQC, featureCounts 2.1.1, salmon, RSeQC, fastp 1.1.0 … | `*.environment.yml`, `TOOL_VERSIONS.tsv` |

```bash
conda env create -f envs/lmna_r.environment.yml      # likewise for the other environments
export CONDA_BASE=/path/to/miniforge3                # if not ~/miniforge3 (or set RS / PY_IQM / CONDA_ENVS)
```

The `*.environment.yml` files were exported on aarch64; on x86-64 a few build-specific packages may need to be
resolved by conda.

## 5. How to run

```bash
cp -a HGPS_meta_QAOA /work/hgps_run && cd /work/hgps_run   # outputs are written next to the scripts
bash run_all.sh                    # all stages: 01 01test 02 03 04
bash run_all.sh 02 03              # per-dataset DE + meta-analysis only (~15 min)
QUICK=1 bash run_all.sh 04         # stage 04 with shorter QAOA optimisations (verification mode)
bash fetch_references.sh           # GENCODE files, only needed for the isoform-switch scripts of stage 01
```

| Stage | Content |
|---|---|
| `01` | downstream analyses of the two reprocessed datasets from the included count / rMATS / salmon tables |
| `01test` | raw-read steps (download, QC, trimming, STAR, RSeQC, featureCounts, salmon, rMATS) on the tiny test subsets |
| `02` | DE of GSE231761 and GSE306264 from the GEO count files |
| `03` | k = 3 meta-analysis, cross-species comparison, pathway analyses, sensitivity analysis |
| `04` | QUBO construction, classical benchmark, QAOA simulation, IQM execution (simulator fallback), scoring of the shipped hardware counts |

Full raw-read reprocessing (~100 GB FASTQ) is not run by `run_all.sh`; place the FASTQ files in
`<dataset>/raw_data/fastq/`, run `fetch_references.sh`, and execute the numbered scripts in `scripts/`
(`run_pipeline.sh` / `run_remaining*.sh`).

### Quantum hardware

Without credentials every hardware script falls back to a local simulator (IQM fake backends
`IQMFakeApollo` ≈ Garnet, `IQMFakeAphrodite` ≈ Emerald, or Aer) and writes to files or folders suffixed
`_SIMULATOR_FALLBACK`, so the shipped hardware results are never overwritten. Scoring scripts (13, 14, 18)
always read the shipped hardware counts. To submit to IQM hardware (billed), copy `.env.example` to `.env`,
fill in `IQM_TOKEN` / `IQM_SERVER_URL`, and run with `USE_HARDWARE=1`.

### Reproducibility notes

- `PYTHONHASHSEED=0` is set by `run_all.sh`; the networkx node order (= qubit ↔ gene order) depends on it.
- `STRING_CACHE=1` (default) uses the shipped STRING query. The live STRING service has changed since the
  analysis (the 19-gene component becomes 18 genes), which would no longer match the hardware data.
- `11_qrisp_qaoa_iqm.py --mode local` re-optimises the QAOA angles and overwrites the angle files that were
  used on hardware, so `run_all.sh` does not call it.
- IQM fake-backend sampling and the live KEGG REST queries are not seeded; such values reproduce to about the
  third decimal. Wall-clock timings (Fig. 4E) depend on the machine.

## 6. Known limitations of the shipped analysis

- Only the eight-gene control instance gave a distribution-level signal above equal-budget random sampling
  (p = 5.5 × 10⁻⁴); the 16- and 19-gene instances are at random level at both depths (`benchmark_table.txt`),
  as reported in the manuscript.
- The genetic algorithm reached the exact STRING optimum in 17 of 20 restarts (best of 20 = exact optimum);
  simulated annealing reached it in 20 of 20 (`classical_benchmark.tsv`).
- The eight-gene control instance (`coexp_ising_8_*`: the seven mouse-FDR-supported genes plus *NDRG2*,
  edges from the same GSE306264 co-expression matrix) is shipped as an input; the call that generated it is
  not part of `02_build_coexpression_qubo.py`. The same applies to the unused `coexp_ising_10/12` files.
- `03_meta_analysis/scripts/06_k3_figures.R` also computes the GSE198983 positional table read by
  `10_fourway_common.R`, so it must run before 10 (as in `run_all.sh`).
- KEGG enrichment in `GSE198983/10_enrichment.R` depends on the live KEGG REST API.
- `GSE165409/12_enrichment.R` runs fgsea/gseGO without a fixed seed. The shipped
  `results/12_enrichment/` tables are the ones used for Fig. 3C; a re-run reproduces them to within about
  ±0.02 NES with the same significance calls for the six recurrent pathways.
- `11_isoformswitch.R` (GSE165409) `switchPlot()` calls fail as in the original run, producing empty PDFs.

## 7. Verification

`VERIFICATION_LOG.md` lists every script with command, environment, runtime, result and verification level.

## 8. License and citation

Code: MIT License (`LICENSE`). Result tables and the IQM measurement counts are released under CC BY 4.0.
The GEO/SRA datasets remain subject to the terms of their original submitters. Please cite the manuscript
and the Zenodo archive (`CITATION.cff`).

Contact: Taeho Kwon (kwon@kribb.re.kr), Sung-Hwan Moon (moonsh@cau.ac.kr).
