#!/usr/bin/env python3
"""
(1) p=3 이상 현상 검증
    coexp_ising_16 에서 p=3 근사비(0.353)가 p=2(0.549)보다 낮게 나왔다.
    QAOA 는 원리상 p 가 커지면 단조 개선되어야 하므로, 이는 물리가 아니라
    COBYLA 가 6개 파라미터 지형에서 국소최소에 갇힌 '최적화 실패'일 가능성이 높다.
    재시작 수를 3 -> 12 로 늘려 확인한다.

(2) 벤치마크 그림 (PNG 300dpi)
    F1 근사비 vs QAOA 층수 p        (고전 SA/GA/그리디 기준선 함께)
    F2 최적해 확률 vs p
    F3 SWAP 오버헤드 / 회로 깊이 vs 문제 크기  <- 핵심 결과
    F4 잡음 스캔
    F5 예상 회로 충실도 vs 문제 크기 (실기 실행 가능 한계)
"""
from __future__ import annotations

import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize

from qiskit.circuit.library import QAOAAnsatz
from qiskit.quantum_info import SparsePauliOp, Statevector

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
OUT = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")
FIG = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "figures")
os.makedirs(FIG, exist_ok=True)
RNG = np.random.default_rng(23)

C_SA, C_GA, C_GR = "#55A868", "#8172B3", "#937860"
C_Q = "#C44E52"


def qubo_to_ising(Q):
    n = Q.shape[0]
    h = np.zeros(n); J = np.zeros((n, n))
    for i in range(n):
        h[i] += -Q[i, i] / 2.0
        for j in range(i + 1, n):
            w = Q[i, j] + Q[j, i]
            if w == 0:
                continue
            h[i] += -w / 4.0; h[j] += -w / 4.0; J[i, j] += w / 4.0
    terms = []
    for i in range(n):
        if abs(h[i]) > 1e-12:
            p = ["I"] * n; p[n - 1 - i] = "Z"; terms.append(("".join(p), h[i]))
    for i in range(n):
        for j in range(i + 1, n):
            if abs(J[i, j]) > 1e-12:
                p = ["I"] * n; p[n - 1 - i] = "Z"; p[n - 1 - j] = "Z"
                terms.append(("".join(p), J[i, j]))
    return SparsePauliOp.from_list(terms)


def all_energies(Q):
    n = Q.shape[0]
    idx = np.arange(1 << n, dtype=np.uint64)
    bits = np.arange(n, dtype=np.uint64)
    X = ((idx[:, None] >> bits) & 1).astype(np.float64)
    return np.einsum("ij,jk,ik->i", X, Q, X)


def qaoa_best(Q, p, E_all, restarts):
    op = qubo_to_ising(Q)
    ans = QAOAAnsatz(op, reps=p).decompose(reps=3)
    def f(th):
        return float(np.abs(Statevector.from_instruction(
            ans.assign_parameters(th)).data) ** 2 @ E_all)
    best = None
    for _ in range(restarts):
        x0 = RNG.uniform(0, np.pi, ans.num_parameters)
        r = minimize(f, x0, method="COBYLA", options={"maxiter": 400, "rhobeg": 0.4})
        if best is None or r.fun < best.fun:
            best = r
    probs = np.abs(Statevector.from_instruction(ans.assign_parameters(best.x)).data) ** 2
    return float(best.fun), float(probs[np.isclose(E_all, E_all.min(), atol=1e-9)].sum())


# --------------------------------------------------------------------------- #
def main() -> int:
    qa = pd.read_csv(os.path.join(OUT, "qaoa_noiseless.tsv"), sep="\t")
    cl = pd.read_csv(os.path.join(OUT, "classical_benchmark.tsv"), sep="\t")
    rs = pd.read_csv(os.path.join(OUT, "resource_estimate.tsv"), sep="\t")
    ns = pd.read_csv(os.path.join(OUT, "qaoa_noise_scan.tsv"), sep="\t")

    # ------------------------------------------------------------- (1) p=3 검증
    print("=" * 80)
    print(" (1) p=3 이상 현상 검증 : coexp_ising_16")
    print("=" * 80)
    tag = "coexp_ising_16"
    Q = np.load(os.path.join(OUT, f"{tag}_Q.npy"))
    c = json.load(open(os.path.join(OUT, f"{tag}_classical.json")))
    E_all = all_energies(Q)
    e_ref, e_rand = c["E_ref"], c["E_random_mean"]
    denom = e_rand - e_ref
    prev = qa[(qa.instance == tag) & (qa.p == 3)].iloc[0]
    print(f"  기존(재시작 3회): ratio = {prev.ratio_expect:.4f}, P_opt = {prev.P_optimal*100:.2f}%")
    t0 = time.time()
    e3, p3 = qaoa_best(Q, 3, E_all, restarts=12)
    r3 = (e_rand - e3) / denom
    print(f"  재검증(재시작 12회): ratio = {r3:.4f}, P_opt = {p3*100:.2f}%  ({time.time()-t0:.0f}초)")
    if r3 > prev.ratio_expect + 1e-3:
        print("  -> 개선됨. 원인은 물리가 아니라 COBYLA 국소최소 (최적화 실패) 로 확정.")
        qa.loc[(qa.instance == tag) & (qa.p == 3), ["E_expect", "ratio_expect", "P_optimal"]] = \
            [e3, r3, p3]
        qa.loc[(qa.instance == tag) & (qa.p == 3), "note"] = "re-optimized with 12 restarts"
        qa.to_csv(os.path.join(OUT, "qaoa_noiseless.tsv"), sep="\t", index=False)
    else:
        print("  -> 재시작을 늘려도 개선 없음. 파라미터 지형 자체가 험한 인스턴스로 기록.")

    # ------------------------------------------------------------- (2) 그림
    print("\n" + "=" * 80)
    print(" (2) 벤치마크 그림 생성")
    print("=" * 80)
    plt.rcParams.update({"figure.dpi": 300, "savefig.dpi": 300, "font.size": 9,
                         "axes.grid": True, "grid.alpha": .3})

    insts = sorted(qa.instance.unique())
    lab = {"string_maxcut_19": "STRING Max-Cut (N=19)",
           "coexp_ising_16": "co-expression Ising (N=16)"}

    # F1 근사비 vs p
    fig, axes = plt.subplots(1, len(insts), figsize=(5.2 * len(insts), 4.2))
    axes = np.atleast_1d(axes)
    for ax, ins in zip(axes, insts):
        d = qa[qa.instance == ins].sort_values("p")
        cc = cl[cl.instance == ins].iloc[0]
        ax.axhline(1.0, color="k", ls="-", lw=1.2, label=f"optimum (exact)" if cc.ref_source == "exact" else "best known")
        ax.axhline(cc.ratio_SA, color=C_SA, ls="--", lw=1.4, label=f"SA  ({cc.ratio_SA:.3f})")
        ax.axhline(cc.ratio_GA, color=C_GA, ls=":", lw=1.6, label=f"GA  ({cc.ratio_GA:.3f})")
        ax.axhline(cc.ratio_greedy, color=C_GR, ls="-.", lw=1.2, label=f"greedy ({cc.ratio_greedy:.3f})")
        ax.plot(d.p, d.ratio_expect, "o-", color=C_Q, lw=2, ms=7, label="QAOA (expectation)")
        ax.plot(d.p, d.ratio_best, "s--", color=C_Q, alpha=.5, ms=6, label="QAOA (best of 1024 samples)")
        ax.set_xticks(d.p); ax.set_xlabel("QAOA layers p"); ax.set_ylim(-0.05, 1.1)
        ax.set_ylabel("approximation ratio\n(random = 0, optimum = 1)")
        ax.set_title(lab.get(ins, ins)); ax.legend(fontsize=7, loc="lower right")
    fig.suptitle("QAOA vs classical solvers on progeria meta-signature QUBOs (noiseless simulator)",
                 fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "F1_ratio_vs_p.png")); plt.close(fig)

    # F2 최적해 확률
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for ins, mk in zip(insts, ["o-", "s-"]):
        d = qa[qa.instance == ins].sort_values("p")
        ax.plot(d.p, d.P_optimal * 100, mk, lw=2, ms=7, label=lab.get(ins, ins))
        n = int(d.n_qubits.iloc[0])
        ax.axhline(100 / 2 ** n, ls=":", lw=1, alpha=.6)
        ax.text(d.p.max(), 100 / 2 ** n, f"  random 1/2^{n}", va="bottom", fontsize=7, alpha=.7)
    ax.set_yscale("log"); ax.set_xticks(sorted(qa.p.unique()))
    ax.set_xlabel("QAOA layers p"); ax.set_ylabel("probability of sampling the optimum (%)")
    ax.set_title("Optimal-state probability grows with circuit depth")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "F2_optimal_probability.png")); plt.close(fig)

    # F3 SWAP 오버헤드 / 깊이  <- 핵심
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    for dev, mk, col in [("IQM Garnet", "o-", "#4C72B0"), ("IQM Emerald", "s-", "#C44E52")]:
        d = rs[(rs.device == dev) & (rs.p == 1)].sort_values("n_qubits")
        if not d.empty:
            axes[0].plot(d.n_qubits, d.swap_overhead, mk, color=col, lw=2, ms=8, label=f"{dev}  p=1")
            axes[1].plot(d.n_qubits, d.depth_grid, mk, color=col, lw=2, ms=8, label=f"{dev}  p=1")
        d2 = rs[(rs.device == dev) & (rs.p == 2)].sort_values("n_qubits")
        if not d2.empty:
            axes[0].plot(d2.n_qubits, d2.swap_overhead, mk.replace("-", "--"), color=col,
                         lw=1.5, ms=6, alpha=.6, label=f"{dev}  p=2")
            axes[1].plot(d2.n_qubits, d2.depth_grid, mk.replace("-", "--"), color=col,
                         lw=1.5, ms=6, alpha=.6, label=f"{dev}  p=2")
    axes[0].axhline(1.0, color="k", ls=":", lw=1)
    axes[0].set_xlabel("problem size N (qubits)")
    axes[0].set_ylabel("2-qubit gates on grid / all-to-all")
    axes[0].set_title("SWAP overhead from square-lattice routing")
    axes[0].legend(fontsize=7)
    axes[1].set_yscale("log"); axes[1].set_xlabel("problem size N (qubits)")
    axes[1].set_ylabel("transpiled circuit depth")
    axes[1].set_title("Circuit depth after routing")
    axes[1].legend(fontsize=7)
    fig.suptitle("Hardware cost of mapping the QUBO onto IQM-like square lattices\n"
                 "(grid approximation with the same qubit count, not the vendor coupling map)",
                 fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "F3_routing_overhead.png")); plt.close(fig)

    # F4 잡음 스캔
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    for ins, mk in zip(insts, ["o-", "s-"]):
        d = ns[ns.instance == ins].sort_values("err2q")
        ax.plot(d.err2q * 100, d.ratio_mean, mk, lw=2, ms=7, label=lab.get(ins, ins))
    ax.axvspan(0.3, 1.0, alpha=.12, color="grey")
    ax.text(0.62, ax.get_ylim()[1] * .96, "typical current\n2q error", ha="center", va="top", fontsize=7)
    ax.set_xlabel("two-qubit depolarizing error (%)")
    ax.set_ylabel("approximation ratio (mean over shots)")
    ax.set_title("Noise degradation at p = 1")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "F4_noise_scan.png")); plt.close(fig)

    # F5 예상 충실도 = (1-err)^(2q gates)
    fig, ax = plt.subplots(figsize=(7, 4.4))
    errs = [0.003, 0.005, 0.01]
    for err, ls in zip(errs, ["-", "--", ":"]):
        for dev, col in [("IQM Garnet", "#4C72B0"), ("IQM Emerald", "#C44E52")]:
            d = rs[(rs.device == dev) & (rs.p == 1)].sort_values("n_qubits")
            if d.empty:
                continue
            fid = (1 - err) ** d.twoq_grid
            ax.plot(d.n_qubits, fid, ls, marker="o", color=col, lw=1.6, ms=6,
                    label=f"{dev}, 2q err {err*100:.1f}%")
    ax.axhline(0.1, color="k", ls=":", lw=1)
    ax.text(ax.get_xlim()[1], 0.105, "10% — signal essentially lost", ha="right", fontsize=7)
    ax.set_yscale("log"); ax.set_xlabel("problem size N (qubits)")
    ax.set_ylabel(r"estimated circuit fidelity  $(1-\epsilon)^{N_{2q}}$")
    ax.set_title("Projected fidelity on hardware at p = 1\n"
                 "practical ceiling for these sparse biological QUBOs", fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "F5_projected_fidelity.png")); plt.close(fig)

    # 충실도 표 저장
    rows = []
    for _, r in rs[rs.p == 1].iterrows():
        for err in errs:
            rows.append(dict(device=r.device, n_qubits=r.n_qubits, p=1,
                             twoq_grid=r.twoq_grid, depth_grid=r.depth_grid,
                             err2q=err, fidelity=round(float((1 - err) ** r.twoq_grid), 5)))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "projected_fidelity.tsv"), sep="\t", index=False)

    print(f"  F1~F5 저장 -> {FIG}")
    for f in sorted(os.listdir(FIG)):
        print(f"    {f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
