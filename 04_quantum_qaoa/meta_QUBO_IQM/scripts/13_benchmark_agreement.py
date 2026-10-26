#!/usr/bin/env python3
"""
양자 vs 고전 벤치마킹 — 최종 정리

목표: "양자컴 답이 고전컴 답과 얼마나 같은가" 를 정량화한다.
      (양자 우위 주장이 아니라 재현 범위의 측정)

산출
  1) 일치율 표/그림  : 양자 · 무작위 · 고전 각각이 정답 분할을 얼마나 재현하는가
  2) 종합 비교표     : 완전탐색 / SA / GA / 그리디 / 무작위 / 양자
"""
from __future__ import annotations

import glob
import json
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
R = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")
OUT = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "figures")  # [deposit] was ~/Desktop
os.makedirs(OUT, exist_ok=True)  # [deposit] output dir may not exist yet
RNG = np.random.default_rng(11)
SHOTS = 4096

C_Q, C_RAND, C_SA, C_GA, C_EX = "#C44E52", "#9E9E9E", "#55A868", "#8172B3", "#4C72B0"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": .25})

# 2큐빗 게이트 수는 18번이 산출한 gate_counts.json 하나만 신뢰한다.
# (이전에는 스크립트마다 논리 게이트 수와 트랜스파일 후 게이트 수가 뒤섞여 있었다)
GC = json.load(open(f"{R}/gate_counts.json"))


def twoq_of(tag, p, fallback):
    k = f"{tag}_p{p}"
    return GC[k]["transpiled_2q"] if k in GC else fallback


# (tag, N, QAOA 깊이, 실기 결과 파일)
# 16/19 유전자는 p=1 과 p=2 를 모두 5회씩 측정했으므로 두 줄로 나란히 비교한다.
INST = [("coexp_ising_8", 8, 1, "qrisp_hw_coexp_ising_8_p1_garnet_x5.json"),
        ("coexp_ising_16", 16, 1, "qrisp_hw_coexp_ising_16_p1_garnet_x5.json"),
        ("coexp_ising_16", 16, 2, "qrisp_hw_coexp_ising_16_p2_garnet_x5.json"),
        ("string_maxcut_19", 19, 1, "qrisp_hw_string_maxcut_19_p1_garnet_x5.json"),
        ("string_maxcut_19", 19, 2, "qrisp_hw_string_maxcut_19_p2_garnet_x5.json"),
        # 34유전자는 Emerald 이고 statevector 트랜스파일 산출 대상이 아니어서 계획값을 쓴다
        ("coexp_ising_34", 34, 1, "qrisp_hw_coexp_ising_34_p1_emerald.json")]


def Eb(Q, X):
    return np.einsum("ij,jk,ik->i", X, Q, X)


def agree(a, b):
    """두 분할의 일치율. 전역 반전(0<->1)은 같은 분할이므로 큰 쪽을 취한다."""
    a, b = np.asarray(a), np.asarray(b)
    return float(max((a == b).mean(), (a == 1 - b).mean()))


def greedy(Q):
    n = Q.shape[0]
    x = np.zeros(n)
    改 = True
    while 改:
        改 = False
        base = float(x @ Q @ x)
        for i in range(n):
            x[i] = 1 - x[i]
            e = float(x @ Q @ x)
            if e < base - 1e-12:
                base, 改 = e, True
            else:
                x[i] = 1 - x[i]
    return float(x @ Q @ x), x.astype(int)


def main() -> int:
    rows = []
    print("=" * 96)
    print(" 양자 vs 고전 벤치마킹")
    print("=" * 96)

    for tag, n, depth, hwf in INST:
        twoq = twoq_of(tag, depth, 509)     # 34유전자 Emerald 계획값
        Q = np.load(f"{R}/{tag}_Q.npy")
        c = json.load(open(f"{R}/{tag}_classical.json"))
        xopt = np.array(c["x_opt"]) if c.get("x_opt") else None
        eopt, erand = c["E_ref"], c["E_random_mean"]
        den = erand - eopt

        hw = json.load(open(f"{R}/{hwf}")) if os.path.exists(f"{R}/{hwf}") else None
        runs = (hw.get("runs") or [hw]) if hw else []
        q_ratio = np.array([r["ratio_mean"] for r in runs]) if runs else np.array([np.nan])
        q_best = np.array([r["ratio_best"] for r in runs]) if runs else np.array([np.nan])
        q_hit = sum(r["reached_optimum"] for r in runs)
        q_agree = (np.array([agree(r["best_bits"], xopt) for r in runs])
                   if xopt is not None and runs else np.array([np.nan]))
        pval = (stats.ttest_1samp(q_ratio, 0).pvalue if len(q_ratio) > 1 else np.nan)

        # 무작위 기준선 (같은 shot 예산)
        rb, ra, rh = [], [], []
        for _ in range(200):
            X = RNG.integers(0, 2, (SHOTS, n)).astype(float)
            E = Eb(Q, X)
            k = int(np.argmin(E))
            rb.append((erand - E[k]) / den)
            rh.append(abs(E[k] - eopt) < 1e-9)
            if xopt is not None:
                ra.append(agree(X[k].astype(int), xopt))

        # 고전 기준선
        t0 = time.time()
        eg, xg = greedy(Q)
        t_gr = time.time() - t0
        gr_agree = agree(xg, xopt) if xopt is not None else np.nan

        rows.append(dict(
            instance=tag, n=n, twoq=twoq, depth=depth, space=2 ** n,
            E_opt=eopt, ref=c.get("ref_source", "exact"),
            q_mean=float(np.nanmean(q_ratio)), q_mean_sd=float(np.nanstd(q_ratio, ddof=1)) if len(q_ratio) > 1 else 0.0,
            q_best=float(np.nanmean(q_best)), q_best_sd=float(np.nanstd(q_best, ddof=1)) if len(q_best) > 1 else 0.0,
            q_hit=q_hit, q_runs=len(runs), q_pvalue=float(pval),
            q_agree=float(np.nanmean(q_agree)) * 100,
            q_agree_sd=(float(np.nanstd(q_agree, ddof=1)) * 100 if len(q_agree) > 1 else 0.0),
            rand_best=float(np.mean(rb)), rand_best_sd=float(np.std(rb, ddof=1)),
            rand_hit=float(np.mean(rh)) * 100,
            rand_agree=(float(np.mean(ra)) * 100 if ra else np.nan),
            greedy_ratio=(erand - eg) / den, greedy_agree=gr_agree * 100, t_greedy=t_gr,
            sa_hit=c.get("SA_hit"), ga_hit=c.get("GA_hit"),
            t_sa=c.get("t_SA"), t_ga=c.get("t_GA")))
        # 초기 인스턴스(16/19/34)의 SA/GA 결과는 03번이 별도 TSV 에 저장했다
        if rows[-1]["sa_hit"] is None:
            import csv
            tsv = f"{R}/classical_benchmark.tsv"
            if os.path.exists(tsv):
                for row in csv.DictReader(open(tsv), delimiter="\t"):
                    if row["instance"] == tag:
                        rows[-1]["sa_hit"] = int(row["SA_hit"])
                        rows[-1]["ga_hit"] = int(row["GA_hit"])
                        rows[-1]["t_sa"] = float(row["t_SA"])
                        rows[-1]["t_ga"] = float(row["t_GA"])

        print(f"\n{tag}  (N={n}, 2q {twoq}, 탐색공간 {2**n:,})")
        print(f"  양자   평균 {np.nanmean(q_ratio):+.4f}±{np.nanstd(q_ratio, ddof=1) if len(q_ratio)>1 else 0:.4f}"
              f"  최선 {np.nanmean(q_best):.4f}  최적해 {q_hit}/{len(runs)}"
              f"  일치율 {np.nanmean(q_agree)*100:.1f}%  p={pval:.3g}")
        print(f"  무작위 최선 {np.mean(rb):.4f}  최적해 {np.mean(rh)*100:.0f}%"
              f"  일치율 {(np.mean(ra)*100 if ra else float('nan')):.1f}%")
        print(f"  고전   SA {c.get('SA_hit')}/20  GA {c.get('GA_hit')}/20  "
              f"그리디 ratio {(erand-eg)/den:.4f} 일치율 {gr_agree*100:.1f}%")

    with open(f"{R}/benchmark_agreement.json", "w") as fh:
        json.dump(rows, fh, indent=2)

    # ==================================================================== 그림
    ok = [r for r in rows if not np.isnan(r["q_agree"])]
    x = np.arange(len(ok))
    w = .26
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4))
    fig.suptitle("Quantum versus classical on the progeria gene panels",
                 fontsize=15.5, fontweight="bold", y=0.99)
    fig.text(0.5, 0.925, "how closely does the quantum answer reproduce the exact classical answer?",
             ha="center", fontsize=10, color="#5A6672", style="italic")
    lab = [f"{r['n']} genes\n({r['twoq']} 2q gates)" for r in ok]

    # --- 1) 일치율
    ax = axes[0]
    ax.bar(x - w, [r["q_agree"] for r in ok], w, yerr=[r["q_agree_sd"] for r in ok],
           capsize=4, color=C_Q, alpha=.88, label="quantum (IQM Garnet)")
    ax.bar(x, [r["rand_agree"] for r in ok], w, color=C_RAND, alpha=.88,
           label="random sampling, same budget")
    ax.bar(x + w, [r["greedy_agree"] for r in ok], w, color=C_SA, alpha=.88,
           label="greedy (classical)")
    ax.axhline(100, color="k", ls="--", lw=1)
    ax.text(len(ok) - .5, 100.8, "identical to the exact answer", ha="right", fontsize=7.5)
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=8.5)
    ax.set_ylabel("genes assigned to the same module as the exact answer (%)")
    ax.set_ylim(40, 108); ax.legend(fontsize=8, loc="lower left")
    ax.set_title("A   Agreement with the exact classical answer",
                 fontsize=11.5, fontweight="bold", loc="left")
    for i, r in enumerate(ok):
        ax.text(i - w, r["q_agree"] + 1.4, f"{r['q_agree']:.0f}%", ha="center",
                fontsize=8, color=C_Q, fontweight="bold")

    # --- 2) 최적해 재현률
    ax = axes[1]
    ax.bar(x - w/2, [r["q_hit"] / r["q_runs"] * 100 for r in ok], w,
           color=C_Q, alpha=.88, label="quantum")
    ax.bar(x + w/2, [r["rand_hit"] for r in ok], w, color=C_RAND, alpha=.88,
           label="random sampling")
    for i, r in enumerate(ok):
        ax.text(i - w/2, r["q_hit"] / r["q_runs"] * 100 + 2,
                f"{r['q_hit']}/{r['q_runs']}", ha="center", fontsize=8, color=C_Q)
    ax.axhline(100, color=C_SA, ls="-", lw=2.2)
    ax.text(len(ok) - .5, 102, "simulated annealing:  20/20 in 1.3 s", ha="right",
            fontsize=8, color=C_SA, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=8.5)
    ax.set_ylabel("runs reaching the exact optimum (%)")
    ax.set_ylim(0, 118); ax.legend(fontsize=8, loc="upper right")
    ax.set_title("B   Reaching the exact optimum",
                 fontsize=11.5, fontweight="bold", loc="left")

    # --- 3) 분포 기대값 (유일하게 무작위와 구분되는 지표)
    ax = axes[2]
    m = np.array([r["q_mean"] for r in ok])
    sd = np.array([r["q_mean_sd"] for r in ok])
    cols = [C_Q if r["q_pvalue"] < .05 else "#C9A0A2" for r in ok]
    ax.bar(x, m, .5, yerr=sd, capsize=5, color=cols, alpha=.9)
    ax.axhline(0, color="k", lw=1.4)
    ax.text(-0.45, 0.004, "0 = same as random", fontsize=8, color="#5A6672")
    for i, r in enumerate(ok):
        tag_ = ("p = %.0e\nsignal present" % r["q_pvalue"]) if r["q_pvalue"] < .05 \
               else ("p = %.2f\nno signal" % r["q_pvalue"])
        ax.text(i, r["q_mean"] + r["q_mean_sd"] + 0.006, tag_, ha="center", fontsize=7.8,
                color="#3E7A4E" if r["q_pvalue"] < .05 else "#8A929A")
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=8.5)
    ax.set_ylabel("mean approximation ratio of the whole distribution")
    ax.set_ylim(-0.012, max(m + sd) * 1.55)
    ax.set_title("C   The only metric random cannot fake",
                 fontsize=11.5, fontweight="bold", loc="left")

    fig.tight_layout(rect=[0, 0.01, 1, 0.915])
    p = os.path.join(OUT, "Fig_quantum_vs_classical.png")
    fig.savefig(p, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"\n그림 저장: {p}")

    # ==================================================================== 표
    hdr = ("  N    p   2q gates   search space      quantum%   random%   greedy%   "
           "quantum opt   SA      GA      p-value")
    print("\n" + "=" * 111)
    print(" 종합 비교표")
    print("=" * 111)
    print(hdr)
    print("-" * 111)
    lines = []
    def fmt(v, w, d=1):
        return f"{v:>{w}.{d}f}" if v == v else f"{'n/a':>{w}}"

    for r in rows:
        pv = f"{r['q_pvalue']:.2g}" if r["q_pvalue"] == r["q_pvalue"] else "n/a"
        sa = f"{r['sa_hit']}/20" if r["sa_hit"] is not None else "n/a"
        ga = f"{r['ga_hit']}/20" if r["ga_hit"] is not None else "n/a"
        ln = (f"{r['n']:>3}{r['depth']:>5}{r['twoq']:>10}{r['space']:>15,}"
              f"{fmt(r['q_agree'], 14)}{fmt(r['rand_agree'], 10)}{fmt(r['greedy_agree'], 10)}"
              f"{(str(r['q_hit']) + '/' + str(r['q_runs'])):>14}"
              f"{sa:>8}{ga:>8}{pv:>11}")
        print(ln); lines.append(ln)
    with open(f"{R}/benchmark_table.txt", "w") as fh:
        fh.write(hdr + "\n" + "-" * 111 + "\n" + "\n".join(lines) + "\n")
    print(f"\n표 저장: {R}/benchmark_table.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
