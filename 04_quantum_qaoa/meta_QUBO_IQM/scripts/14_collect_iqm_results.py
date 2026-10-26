#!/usr/bin/env python3
"""
실제 IQM 양자컴퓨터에서 나온 결과 일체를 바탕화면 한 폴더로 모은다.

  01_summary/       요약표 · README · 실행 이력
  02_raw_hardware/  실기 측정 원본 (counts 포함)
  03_figures/       실기 결과로 만든 그림
  04_settings/      학습된 각도 · 백엔드 정보 · 실행 계획
  05_reference/     비교용 고전 기준선과 문제 정의
"""
from __future__ import annotations

import glob
import json
import os
import shutil
from datetime import datetime

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
R = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")
FIG = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "figures")
DESK = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "figures")  # [deposit] was ~/Desktop
DEST = str(DEPOSIT / "04_quantum_qaoa" / "IQM_hardware_results")  # [deposit] was ~/Desktop/IQM_hardware_results


def _copy_keep(src, dst):
    """[deposit] the deposit's IQM_hardware_results/ already holds the archived real-hardware bundle
    (02_raw_hardware, 04_settings, 06_scripts, 00_README.txt are inputs). Never overwrite those files;
    only copy when the destination does not exist yet."""
    if os.path.exists(dst):
        print(f"  [deposit] keep existing {os.path.relpath(dst, DEST)}")
        return
    shutil.copy2(src, dst)

HW_RUNS = [
    ("qrisp_hw_coexp_ising_8_p1_garnet_x5.json", "coexp_ising_8", 8, 1, "Garnet", 5,
     "8 유전자 · 26 게이트 · 5회 반복  ← 신호가 검출된 유일한 크기"),
    ("qrisp_hw_coexp_ising_16_p1_garnet_x5.json", "coexp_ising_16", 16, 1, "Garnet", 5,
     "16 유전자 · 156 게이트 · 5회 반복  (4-way 핵심 7개 유전자 포함)"),
    ("qrisp_hw_coexp_ising_16_p2_garnet_x5.json", "coexp_ising_16", 16, 2, "Garnet", 5,
     "16 유전자 · p=2 · 5회 반복  (깊이 대조군)"),
    ("qrisp_hw_coexp_ising_16_p1_garnet.json", "coexp_ising_16", 16, 1, "Garnet", 1,
     "16 유전자 단일 실행 (최초 시도)"),
    ("qrisp_hw_string_maxcut_19_p1_garnet_x5.json", "string_maxcut_19", 19, 1, "Garnet", 5,
     "19 유전자 · p=1 · 5회 반복"),
    ("qrisp_hw_string_maxcut_19_p2_garnet_x5.json", "string_maxcut_19", 19, 2, "Garnet", 5,
     "19 유전자 · p=2 · 5회 반복"),
    ("qrisp_hw_string_maxcut_19_p2_garnet.json", "string_maxcut_19", 19, 2, "Garnet", 1,
     "19 유전자 p=2 단일 실행 — 이 회차만 최적해 도달(이후 5회 반복에서 재현 안 됨)"),
    ("qrisp_hw_string_maxcut_19_p1_garnet.json", "string_maxcut_19", 19, 1, "Garnet", 1,
     "19 유전자 p=1 단일 실행 (최초 실기 실행)"),
    ("qrisp_hw_coexp_ising_34_p1_emerald.json", "coexp_ising_34", 34, 1, "Emerald", 1,
     "34 유전자 · 504 게이트 · Emerald 54큐빗 — 한계 실증"),
]


def main() -> int:
    # [deposit] original: shutil.rmtree(DEST) then rebuild. In the deposit DEST contains the real IQM
    # raw-count inputs, so it is updated in place instead (no deletion).
    for d in ("01_summary", "02_raw_hardware", "03_figures", "04_settings", "05_reference"):
        os.makedirs(f"{DEST}/{d}", exist_ok=True)

    # ---------------------------------------------------------- 02 원본
    idx = []
    for fn, tag, n, p, dev, nrep, desc in HW_RUNS:
        src = f"{R}/{fn}"
        if not os.path.exists(src):
            print(f"  [없음] {fn}")
            continue
        d = json.load(open(src))
        out = f"{n:02d}genes_p{p}_{dev}_x{nrep}.json"
        _copy_keep(src, f"{DEST}/02_raw_hardware/{out}")  # [deposit] never overwrite raw HW counts
        s = d.get("summary")
        idx.append(dict(file=out, instance=tag, n_genes=n, depth=p, device=dev,
                        repeats=nrep, shots=d.get("n_shots", d.get("shots")),
                        ratio_mean=(s["ratio_mean_avg"] if s else d["ratio_mean"]),
                        ratio_mean_sd=(s["ratio_mean_sd"] if s else 0.0),
                        ratio_best=(s["ratio_best_avg"] if s else d["ratio_best"]),
                        ratio_best_sd=(s["ratio_best_sd"] if s else 0.0),
                        optimum_hits=(s["optimum_hits"] if s else int(d["reached_optimum"])),
                        genes=d.get("genes"), module_1=d.get("module_1"),
                        module_2=d.get("module_2"), theta=d.get("theta"),
                        wall_seconds=d.get("wall_s"), note=desc))
    json.dump(idx, open(f"{DEST}/01_summary/run_index.json", "w"), indent=2, ensure_ascii=False)

    # ---------------------------------------------------------- 03 그림
    for src, dst in [(f"{DESK}/Fig_quantum_vs_classical.png", "01_quantum_vs_classical.png"),
                     (f"{DESK}/Fig_hardware_results.png", "02_hardware_results_detail.png"),
                     (f"{DESK}/Fig_genes_to_quantum.png", "03_genes_to_quantum_problem.png"),
                     (f"{FIG}/F3_routing_overhead.png", "04_routing_overhead_sim.png"),
                     (f"{FIG}/F5_projected_fidelity.png", "05_projected_fidelity_sim.png")]:
        if os.path.exists(src):
            shutil.copy2(src, f"{DEST}/03_figures/{dst}")

    # ---------------------------------------------------------- 04 설정
    for f in [g for g in glob.glob(f"{R}/qrisp_angles_*.json") if "SIMULATOR_FALLBACK" not in g] + \
             [f"{R}/iqm_live_backend.json", f"{R}/iqm_run_plan.json",
              f"{R}/iqm_predicted.json"]:
        if os.path.exists(f):
            _copy_keep(f, f"{DEST}/04_settings/{os.path.basename(f)}")  # [deposit] keep HW settings

    # ---------------------------------------------------------- 05 참조
    for tag in ("coexp_ising_8", "coexp_ising_16", "string_maxcut_19", "coexp_ising_34"):
        for suf in ("_meta.json", "_classical.json"):
            src = f"{R}/{tag}{suf}"
            if os.path.exists(src):
                shutil.copy2(src, f"{DEST}/05_reference/{tag}{suf}")
    for f in ("benchmark_agreement.json", "benchmark_table.txt",
              "hardware_summary_table.json", "classical_benchmark.tsv",
              "qubo_instances.tsv", "gate_counts.json"):
        if os.path.exists(f"{R}/{f}"):
            shutil.copy2(f"{R}/{f}", f"{DEST}/01_summary/{f}")

    # 사용한 스크립트도 함께
    os.makedirs(f"{DEST}/06_scripts", exist_ok=True)
    for f in ("11_qrisp_qaoa_iqm.py", "10_iqm_connect_test.py",
              "13_benchmark_agreement.py", "12_hardware_results_figure.py"):
        src = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "scripts" / f)  # [deposit]
        if os.path.exists(src):
            _copy_keep(src, f"{DEST}/06_scripts/{f}")  # [deposit] both copies are maintained in the deposit

    # ---------------------------------------------------------- README
    bt = open(f"{R}/benchmark_table.txt").read() if os.path.exists(f"{R}/benchmark_table.txt") else ""
    _readme = f"{DEST}/00_README.txt"
    if os.path.exists(_readme):  # [deposit] keep the archived README; write the regenerated one alongside
        _readme = f"{DEST}/01_summary/00_README_regenerated.txt"
    with open(_readme, "w", encoding="utf-8") as fh:
        fh.write(f"""IQM 양자컴퓨터 실기 결과 일체
{'=' * 78}
생성 {datetime.now():%Y-%m-%d %H:%M}
장비 IQM Garnet (20 큐빗) / IQM Emerald (54 큐빗),  IQM Resonance 클라우드
소프트웨어 Qrisp QAOA + qiskit-iqm

무엇을 한 것인가
{'-' * 78}
조로증 환자 데이터 4종(사람 3 + 마우스 1)에서 공통으로 변한 유전자를 찾고,
그 유전자들의 공발현 관계를 Ising 모형으로 바꿔 실제 양자컴퓨터에서 풀었다.
유전자 1개 = 큐빗 1개. 답은 "유전자를 두 모듈로 나누는 최적 분할"이다.

핵심 설계: 각도 최적화는 로컬 시뮬레이터에서 끝내고,
          양자컴퓨터에서는 확정된 회로를 1회 측정만 했다.
          (변분 루프를 실기에 걸면 회로 수백 개가 필요해 예산 초과)

결과 요약
{'-' * 78}
{bt}
읽는 법
  quantum%  : 양자컴이 낸 분할이 고전 정답과 같은 비율
  random%   : 같은 4,096 shots 를 무작위로 뽑았을 때의 같은 지표
  greedy%   : 고전 그리디 알고리즘
  p-value   : 분포 기대값이 무작위(0)와 다른지에 대한 검정 (n=5)

네 줄 결론
{'-' * 78}
1. 양자컴 답은 고전 정답과 8유전자 100%, 16유전자 92.5%(p=1) / 95.0%(p=2),
   19유전자 83.2%(p=1) / 80.0%(p=2) 일치했다.

2. 그러나 같은 예산의 무작위 추출도 100% / 92.4% / 92.0% / 81.8% / 81.9% 로
   사실상 동일하다. 탐색공간이 작아 4,096번 뽑으면 무작위로도 좋은 답이 나오므로,
   일치율만으로는 양자가 기여했다고 말할 수 없다.

3. 무작위가 흉내낼 수 없는 지표는 '분포 전체의 기대값' 하나뿐이며,
   8유전자(트랜스파일 후 31 게이트)에서만 유의했고(0.071 +- 0.016, p=5.5e-4),
   16유전자(156 게이트) 이상에서는 모두 소실되었다
   (16 p=1 p=0.29 / 16 p=2 p=0.39 / 19 p=1 p=0.64 / 19 p=2 p=0.17).
   -> 게이트 수 순으로 보면 31개(유의) 다음이 76개(p=0.64)이므로,
      현행 초전도 하드웨어의 실용 한계는 2큐빗 게이트 31 ~ 76개 사이에 있다.
      그 사이를 좁히려면 10~14 유전자 인스턴스를 추가 측정하면 된다.

4. 깊이를 p=1 -> p=2 로 늘려도 두 인스턴스 모두에서 개선되지 않았다.
   무잡음 시뮬레이터에서는 깊이가 늘면 좋아져야 하므로, 이 역전은
   게이트 수 증가에 따른 잡음 누적이 알고리즘 이득을 상쇄했음을 뜻한다.
   16/19 유전자를 p=1, p=2 각각 5회씩 대칭 측정해 이 결론을 확인했다.

반드시 함께 보고할 것
{'-' * 78}
  시뮬레이티드 어닐링은 모든 문제를 1.3초에 20/20 최적해로 해결한다.
  이 규모에서 양자 우위는 없다. 본 실험의 목적은 우위 입증이 아니라
  '현행 하드웨어가 고전 답을 어디까지 재현하는가'의 정량화다.

폴더 구성
{'-' * 78}
  01_summary/       요약표, 실행 이력(run_index.json), 고전 기준선
  02_raw_hardware/  실기 측정 원본 7건 (모든 비트열과 확률 포함)
  03_figures/       그림 5장
  04_settings/      학습된 각도, 백엔드 정보, 실행 계획, 예측값
  05_reference/     문제 정의(유전자 목록·결합상수)와 고전 정답
  06_scripts/       실행에 사용한 스크립트

주의
{'-' * 78}
  garnet:mock 은 잡음 시뮬레이터가 아니라 스텁이다.
  8개 상태가 균일 확률로 반환되므로 수치를 신뢰하면 안 된다.
  배관 검증용으로만 사용했고, 이 폴더의 결과는 전부 실제 QPU 측정이다.
""")

    tot = sum(os.path.getsize(os.path.join(r, f))
              for r, _, fs in os.walk(DEST) for f in fs)
    n = sum(len(fs) for _, _, fs in os.walk(DEST))
    print(f"정리 완료: {DEST}")
    print(f"  파일 {n}개, {tot/1024/1024:.1f} MB")
    for d in sorted(os.listdir(DEST)):
        pp = f"{DEST}/{d}"
        if os.path.isdir(pp):
            print(f"  {d:20s} {len(os.listdir(pp)):2d}개")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
