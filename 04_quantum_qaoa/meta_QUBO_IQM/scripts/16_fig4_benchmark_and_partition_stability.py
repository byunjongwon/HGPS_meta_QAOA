#!/usr/bin/env python3
"""
Figure 4 계열 그림 + 분할 안정성 수치 (조로증 / IQM)
  [release] 원래 이름 16_asfv_style_figures.py (ASFV 논문과 같은 그림 형식이라는 뜻).
            Fig 4E 실행시간, Fig 4F-H 분할 일치도·유전자별 안정성 수치가 이 스크립트에서 나오므로
            공개본에 포함했다. 유전자별 안정성은 표로도 저장한다 (results/partition_stability_*.tsv).

  Fig A  네트워크 위상 (STRING PPI vs 공발현), 노드색 = 커뮤니티
  Fig B  QUBO 계수 행렬 히트맵
  Fig C  QAOA 에너지 지형 (gamma, beta)
  Fig D  해 품질 비교 : 완전탐색 / SA / GA / 그리디 / IQM QAOA
  Fig E  실행시간 비교
  Fig F  유전자별 커뮤니티 배정 일관성 + 정답 대비 일치율
  Fig G  분할 간 쌍별 일치도 히트맵
  Fig H  회로 깊이 대 해 품질
"""
from __future__ import annotations

import json
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.colors import TwoSlopeNorm

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
R = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")
OUT = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "figures" / "fig4")  # [deposit] was ~/Desktop/IQM_hardware_results/03_figures/asfv_style
RNG = np.random.default_rng(4)

NETS = [("string_maxcut_19", "STRING PPI", 19), ("coexp_ising_16", "Co-expression", 16)]
C_A, C_B = "#4C72B0", "#C44E52"
C_POS, C_NEG = "#C44E52", "#4C72B0"
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": .22,
                     "axes.titlesize": 13, "axes.labelsize": 11.5,
                     "xtick.labelsize": 10.5, "ytick.labelsize": 10.5,
                     "legend.fontsize": 10, "figure.dpi": 110})

SHORT = {"Brute force": "Exact\n(brute force)", "Simulated annealing": "SA",
         "Genetic algorithm": "GA", "Greedy": "Greedy"}


def short(k):
    return SHORT.get(k, k)


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    p = f"{OUT}/{name}"
    fig.savefig(p, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  {name}")


def load(tag):
    Q = np.load(f"{R}/{tag}_Q.npy")
    meta = json.load(open(f"{R}/{tag}_meta.json"))
    cls = json.load(open(f"{R}/{tag}_classical.json"))
    G = nx.Graph()
    G.add_nodes_from(meta["nodes"])
    for u, v, w in meta["edges"]:
        G.add_edge(u, v, weight=float(w))
    return Q, meta, cls, G


def Eb(Q, X):
    return np.einsum("ij,jk,ik->i", X, Q, X)


def agree(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(max((a == b).mean(), (a == 1 - b).mean()))


def sa(Q, steps=20000, seed=0):
    r = np.random.default_rng(seed)
    n = Q.shape[0]
    x = r.integers(0, 2, n).astype(float)
    e = float(x @ Q @ x)
    b, bx = e, x.copy()
    T0, T1 = max(1.0, np.abs(Q).max()), 1e-3
    for t in range(steps):
        T = T0 * (T1 / T0) ** (t / steps)
        i = r.integers(n)
        x[i] = 1 - x[i]
        en = float(x @ Q @ x)
        if en - e <= 0 or r.random() < np.exp(-(en - e) / max(T, 1e-12)):
            e = en
            if e < b:
                b, bx = e, x.copy()
        else:
            x[i] = 1 - x[i]
    return b, bx.astype(int)


def ga(Q, gens=300, pop=60, seed=0):
    r = np.random.default_rng(seed)
    n = Q.shape[0]
    P = r.integers(0, 2, (pop, n)).astype(float)
    for _ in range(gens):
        P = P[np.argsort(Eb(Q, P))]
        el = P[: pop // 4]
        ch = []
        while len(ch) < pop - len(el):
            a, b = el[r.integers(len(el))], el[r.integers(len(el))]
            c = r.integers(1, n)
            k = np.concatenate([a[:c], b[c:]])
            m = r.random(n) < .05
            k[m] = 1 - k[m]
            ch.append(k)
        P = np.vstack([el, np.array(ch)])
    E = Eb(Q, P)
    i = int(np.argmin(E))
    return float(E[i]), P[i].astype(int)


def greedy(Q):
    n = Q.shape[0]
    x = np.zeros(n)
    ch = True
    while ch:
        ch = False
        base = float(x @ Q @ x)
        for i in range(n):
            x[i] = 1 - x[i]
            e = float(x @ Q @ x)
            if e < base - 1e-12:
                base, ch = e, True
            else:
                x[i] = 1 - x[i]
    return float(x @ Q @ x), x.astype(int)


def hw_partitions(tag):
    """실기에서 얻은 최선해 (p 별)."""
    out = {}
    for fn, p in [(f"qrisp_hw_{tag}_p1_garnet_x5.json", 1),
                  (f"qrisp_hw_{tag}_p2_garnet_x5.json", 2),
                  (f"qrisp_hw_{tag}_p1_garnet.json", 1),
                  (f"qrisp_hw_{tag}_p2_garnet.json", 2)]:
        f = f"{R}/{fn}"
        if os.path.exists(f) and p not in out:
            d = json.load(open(f))
            s = d.get("summary")
            out[p] = dict(bits=np.array(d["best_bits"]),
                          ratio_best=(s["ratio_best_avg"] if s else d["ratio_best"]),
                          ratio_best_sd=(s["ratio_best_sd"] if s else 0.0),
                          ratio_mean=(s["ratio_mean_avg"] if s else d["ratio_mean"]),
                          ratio_mean_sd=(s["ratio_mean_sd"] if s else 0.0),
                          hits=(s["optimum_hits"] if s else int(d["reached_optimum"])),
                          runs=(s["n_runs"] if s else 1),
                          wall=d.get("wall_s"), E_best=d["E_best"],
                          # 회차별 원본 (Fig F 가 최고 1회만 그리지 않도록)
                          per_run=d.get("runs") or [dict(rep=1, best_bits=d["best_bits"],
                                                         ratio_best=d["ratio_best"],
                                                         reached_optimum=d["reached_optimum"])])
    return out


def qubo_to_ising_op(Q):
    from qiskit.quantum_info import SparsePauliOp
    n = Q.shape[0]
    h = np.zeros(n); J = np.zeros((n, n))
    for i in range(n):
        h[i] += -Q[i, i] / 2
        for j in range(i + 1, n):
            w = Q[i, j] + Q[j, i]
            if w:
                h[i] += -w / 4; h[j] += -w / 4; J[i, j] += w / 4
    t = []
    for i in range(n):
        if abs(h[i]) > 1e-12:
            s = ["I"] * n; s[n - 1 - i] = "Z"; t.append(("".join(s), h[i]))
    for i in range(n):
        for j in range(i + 1, n):
            if abs(J[i, j]) > 1e-12:
                s = ["I"] * n; s[n - 1 - i] = "Z"; s[n - 1 - j] = "Z"
                t.append(("".join(s), J[i, j]))
    return SparsePauliOp.from_list(t)


def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    D = {}
    print("고전 기준선 계산 중 …")
    for tag, name, n in NETS:
        Q, meta, cls, G = load(tag)
        nodes = meta["nodes"]
        ii = np.arange(1 << n, dtype=np.uint64); bb = np.arange(n, dtype=np.uint64)
        X = ((ii[:, None] >> bb) & 1).astype(float)
        t0 = time.time(); E = Eb(Q, X); t_bf = time.time() - t0
        eopt = float(E.min()); xopt = X[int(np.argmin(E))].astype(int)
        erand = float(E.mean()); den = erand - eopt
        t0 = time.time(); e_sa, x_sa = min((sa(Q, seed=100 + i) for i in range(20)), key=lambda v: v[0]); t_sa = time.time() - t0
        t0 = time.time(); e_ga, x_ga = min((ga(Q, seed=200 + i) for i in range(20)), key=lambda v: v[0]); t_ga = time.time() - t0
        t0 = time.time(); e_gr, x_gr = greedy(Q); t_gr = time.time() - t0
        D[tag] = dict(Q=Q, meta=meta, G=G, nodes=nodes, n=n, name=name,
                      eopt=eopt, xopt=xopt, erand=erand, den=den, E=E,
                      sol={"Brute force": (eopt, xopt, t_bf),
                           "Simulated annealing": (e_sa, x_sa, t_sa),
                           "Genetic algorithm": (e_ga, x_ga, t_ga),
                           "Greedy": (e_gr, x_gr, t_gr)},
                      hw=hw_partitions(tag))
        print(f"  {name} (N={n}) 완료")

    # ================================================================== Fig A
    fig, axes = plt.subplots(1, 2, figsize=(19, 9.2))
    for ax, (tag, name, n) in zip(axes, NETS):
        d = D[tag]; G = d["G"]; nodes = d["nodes"]
        part = d["sol"]["Simulated annealing"][1]
        # spring layout 은 음수 가중치를 처리하지 못해 노드가 뭉친다.
        # 배치에는 |가중치| 를 쓴 보조 그래프를 사용한다 (색·굵기는 원래 부호 유지).
        L = nx.Graph()
        L.add_nodes_from(G.nodes())
        for u_, v_, dd_ in G.edges(data=True):
            L.add_edge(u_, v_, weight=abs(dd_["weight"]))
        pos = nx.kamada_kawai_layout(L, weight="weight")
        pos = nx.spring_layout(L, pos=pos, weight="weight", seed=3, k=1.3, iterations=500)
        for u, v, dd in G.edges(data=True):
            w = dd["weight"]
            ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                    color=C_POS if w > 0 else C_NEG, lw=.6 + 2.2 * abs(w) /
                    max(abs(x["weight"]) for *_, x in G.edges(data=True)),
                    alpha=.5, zorder=1)
        for i, g in enumerate(nodes):
            ax.scatter(*pos[g], s=420, c=C_A if part[i] == 0 else C_B,
                       edgecolors="white", lw=1.8, zorder=2)
        # 라벨을 그래프 중심 반대 방향으로 밀어 노드/엣지와 겹치지 않게 한다
        cx = np.mean([p_[0] for p_ in pos.values()])
        cy = np.mean([p_[1] for p_ in pos.values()])
        for g in nodes:
            dx, dy = pos[g][0] - cx, pos[g][1] - cy
            nrm = np.hypot(dx, dy) or 1.0
            ox, oy = dx / nrm * .135, dy / nrm * .135
            ax.text(pos[g][0] + ox, pos[g][1] + oy, g,
                    ha="center" if abs(ox) < .09 else ("left" if ox > 0 else "right"),
                    va="center", fontsize=10, fontweight="bold", zorder=5, color="#22303C",
                    bbox=dict(boxstyle="round,pad=.20", fc="white", ec="#D5DBE1",
                              lw=.6, alpha=.94))
        ax.set_axis_off(); ax.margins(.30)
        ax.set_title(f"{name}\nn = {G.number_of_nodes()} genes,  m = {G.number_of_edges()} edges,  "
                     f"density = {nx.density(G):.3f}", fontsize=14, fontweight="bold", pad=14)
    axes[0].scatter([], [], s=180, c=C_A, label="community 1")
    axes[0].scatter([], [], s=180, c=C_B, label="community 2")
    axes[0].legend(fontsize=11, loc="lower left", framealpha=.95)
    fig.suptitle("A   Interaction network topology of the progeria QUBO gene subsets",
                 fontsize=17, fontweight="bold", y=1.02)
    fig.text(.5, .975, "node colour = two-community partition from simulated annealing;  "
             "edge colour = positive (red) / negative (blue) coupling",
             ha="center", fontsize=11.5, color="#5A6672", style="italic")
    save(fig, "A_network_topology.png")

    # ================================================================== Fig B
    # 좌측 컬러바 눈금이 우측 패널 유전자명과 겹쳐 간격을 넓힌다
    fig, axes = plt.subplots(1, 2, figsize=(18, 7.6))
    fig.subplots_adjust(wspace=.34)
    for ax, (tag, name, n) in zip(axes, NETS):
        d = D[tag]; Q = d["Q"]
        v = np.abs(Q).max()
        im = ax.imshow(Q, cmap="RdBu_r", norm=TwoSlopeNorm(0, -v, v))
        ax.set_xticks(range(n)); ax.set_yticks(range(n))
        ax.set_xticklabels(d["nodes"], rotation=90, fontsize=9.5)
        ax.set_yticklabels(d["nodes"], fontsize=9.5)
        ax.grid(False)
        nz = int((np.abs(Q) > 1e-12).sum())
        ax.set_title(f"QUBO matrix Q — {name}\n{n} qubits,  {nz} non-zero coefficients,  "
                     f"|Q|max = {v:.2f}", fontsize=13.5, fontweight="bold", pad=12)
        cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.03)
        cb.set_label("coefficient", fontsize=11)
    fig.suptitle("B   QUBO formulation:   minimise  H = x$^T$Q x", fontsize=17,
                 fontweight="bold", y=1.03)
    save(fig, "B_qubo_matrices.png")

    # ================================================================== Fig C
    from qiskit.circuit.library import QAOAAnsatz
    from qiskit.quantum_info import Statevector
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2))
    for ax, (tag, name, n) in zip(axes, NETS):
        d = D[tag]
        op = qubo_to_ising_op(d["Q"])
        ans = QAOAAnsatz(op, reps=1).decompose(reps=3)
        m = 26 if n <= 16 else 18
        gs = np.linspace(0, np.pi, m); bs = np.linspace(0, np.pi / 2, m)
        Z = np.zeros((m, m))
        for i, g in enumerate(gs):
            for j, b in enumerate(bs):
                pr = np.abs(Statevector.from_instruction(
                    ans.assign_parameters([g, b])).data) ** 2
                Z[j, i] = (d["erand"] - float(pr @ d["E"])) / d["den"]
        im = ax.imshow(Z, origin="lower", cmap="viridis",
                       extent=[gs[0], gs[-1], bs[0], bs[-1]], aspect="auto")
        k = np.unravel_index(np.argmax(Z), Z.shape)
        ax.plot(gs[k[1]], bs[k[0]], "r*", ms=18, label=f"best  ratio {Z.max():.3f}")
        ax.set_xlabel(r"QAOA angle  $\gamma$", fontsize=12)
        ax.set_ylabel(r"QAOA angle  $\beta$", fontsize=12)
        ax.set_title(f"{name}  ({n} genes,  p = 1)", fontsize=13.5, fontweight="bold", pad=12)
        ax.legend(fontsize=10.5, loc="upper right", framealpha=.92); ax.grid(False)
        cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.03)
        cb.set_label("approximation ratio", fontsize=11)
    fig.suptitle("C   QAOA energy landscape  (noiseless simulator, p = 1)", fontsize=17,
                 fontweight="bold", y=1.03)
    fig.text(.5, .955, "brighter = better;  the red star marks the angles used on hardware",
             ha="center", fontsize=11.5, color="#5A6672", style="italic")
    save(fig, "C_energy_landscape.png")

    # ================================================================== Fig D
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.8))
    for ax, (tag, name, n) in zip(axes, NETS):
        d = D[tag]
        labs, vals, errs, cols = [], [], [], []
        for k, (e, x, t) in d["sol"].items():
            r_ = (d["erand"] - e) / d["den"]
            labs.append(short(k)); vals.append(r_); errs.append(0)
            cols.append("#55A868" if k == "Brute force" else "#7F8C9A")
        lo_ = min(vals + [0])
        for p, h in sorted(d["hw"].items()):
            labs.append(f"IQM Garnet\nQAOA  p={p}")
            vals.append(h["ratio_best"]); errs.append(h["ratio_best_sd"]); cols.append(C_B)
        # p=2 는 19유전자에서만 실행했다. 16유전자 패널에서 그 칸이 그냥 빠져 있으면
        # '실패' 로 오해되므로, 점선 빈 슬롯을 두어 미실행임을 명시한다.
        ghost = 2 not in d["hw"]
        if ghost:
            labs.append("IQM Garnet\nQAOA  p=2"); vals.append(0.0)
            errs.append(0); cols.append("none")
        xx = np.arange(len(labs))
        ax.bar(xx, vals, .60, yerr=errs, capsize=5, color=cols, alpha=.9)
        if ghost:
            ax.bar(xx[-1], 1.0, .60, fc="none", ec="#B9C3CC", ls="--", lw=1.6, zorder=1)
            # 슬롯 폭이 좁으므로 안에는 짧게만 쓰고, 이유는 부제에 적는다
            ax.text(xx[-1], .50, "not\nrun", ha="center", va="center",
                    fontsize=13, color="#8B959E", style="italic", linespacing=1.5)
        ax.axhline(1.0, color="#55A868", ls="--", lw=2.0)
        ax.set_xticks(xx); ax.set_xticklabels(labs, fontsize=11)
        ax.set_ylabel("approximation ratio  (best solution)", fontsize=12)
        ax.set_ylim(min(0, lo_ * 1.25), 1.16)
        ax.set_title(f"{name}  ({n} genes)", fontsize=14, fontweight="bold", pad=12)
        for i, v in enumerate(vals):
            if cols[i] == "none":       # 미실행 슬롯에는 수치를 적지 않는다
                continue
            ax.text(i, v + .035, f"{v:.3f}", ha="center", fontsize=11, fontweight="bold",
                    color=cols[i] if cols[i] != "#7F8C9A" else "#42505C")
    fig.suptitle("D   Solution quality:  classical solvers versus IQM QAOA",
                 fontsize=17, fontweight="bold", y=1.03)
    fig.text(.5, .975, "green dashed line = exact optimum;   SA = simulated annealing,  "
             "GA = genetic algorithm;   IQM bars are mean ± SD over 5 hardware repeats",
             ha="center", fontsize=11.5, color="#5A6672", style="italic")
    fig.text(.5, .932, "both instances were measured at p = 1 and p = 2 with five independent "
             "hardware repeats each (4,096 shots per repeat)",
             ha="center", fontsize=11.5, color="#5A6672", style="italic")
    save(fig, "D_solution_quality.png")

    # ================================================================== Fig E
    # 우측 패널의 긴 y축 라벨이 좌측 패널 안으로 넘어가 겹쳤다.
    # 라벨을 짧게 줄이고("wall clock" 은 캡션에 이미 있다) 패널 간격을 넓힌다.
    fig, axes = plt.subplots(1, 2, figsize=(16.5, 5.8))
    fig.subplots_adjust(wspace=.34)
    for ax, (tag, name, n) in zip(axes, NETS):
        d = D[tag]
        keys = list(d["sol"].keys())
        labs = [short(k).replace("\n", " ") for k in keys]
        ts = [d["sol"][k][2] for k in keys]
        cols = ["#7F8C9A"] * len(labs)
        for p, h in sorted(d["hw"].items()):
            if h["wall"]:
                labs.append(f"IQM QAOA  p={p}")
                ts.append(h["wall"] / max(1, h["runs"]))
                cols.append(C_B)
        y = np.arange(len(labs))
        ax.barh(y, ts, .6, color=cols, alpha=.88)
        ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=11)
        ax.invert_yaxis(); ax.set_xscale("log")
        # 막대 끝 숫자가 축 밖으로 나가지 않도록 오른쪽 여유를 넉넉히 둔다
        ax.set_xlim(min(ts) * .30, max(ts) * 12)
        ax.set_xlabel("time per solve (seconds, log scale)", fontsize=12)
        ax.set_title(f"{name}  ({n} genes)", fontsize=14, fontweight="bold", pad=12)
        for i, t in enumerate(ts):
            ax.text(t * 1.30, i, f"{t:.3g} s", va="center", fontsize=10.5, fontweight="bold")
    fig.suptitle("E   Runtime comparison", fontsize=17, fontweight="bold", y=1.04)
    fig.text(.5, .955, "IQM times are wall-clock per hardware execution and exclude "
             "the classical angle optimisation performed on a simulator",
             ha="center", fontsize=11.5, color="#5A6672", style="italic")
    save(fig, "E_runtime.png")

    # ================================================================== Fig F
    for tag, name, n in NETS:
        d = D[tag]
        methods, mats = [], []
        for k, (e, x, t) in d["sol"].items():
            methods.append(short(k).replace("\n", " ")); mats.append(np.asarray(x))
        # 이전 판은 5회 중 가장 좋았던 1회(best_bits)만 "IQM" 이라고 표시해
        # 100% 일치처럼 보였다. 회차를 모두 개별 행으로 보여준다.
        n_cls = len(methods)
        for p, h in sorted(d["hw"].items()):
            for r_ in h["per_run"]:
                methods.append(f"IQM p={p}  run {r_['rep']}"
                               + ("  ★" if r_["reached_optimum"] else ""))
                mats.append(np.asarray(r_["best_bits"]))
        M = np.array(mats)
        ref = np.asarray(d["xopt"])
        for i in range(len(M)):
            if (M[i] == ref).mean() < .5:
                M[i] = 1 - M[i]
        # 안정성은 하드웨어 회차만으로 계산한다 (고전 4종은 전부 정답이라 희석된다)
        ag = (M[n_cls:] == ref).mean(0) * 100
        row_ag = (M == ref).mean(1) * 100
        # [release] save per-gene stability (previously only drawn in the figure)
        with open(f"{R}/partition_stability_{tag}.tsv", "w") as fh:
            fh.write("gene\tpct_hw_runs_matching_exact\n")
            for g_, v_ in zip(d["nodes"], ag):
                fh.write(f"{g_}\t{v_:.1f}\n")
        print(f"  [{tag}] hardware runs={len(M) - n_cls}: genes stable in every run = "
              f"{int((ag == 100).sum())}/{n}, mean per-gene agreement = {ag.mean():.1f}%")
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(18, max(6.6, .46 * max(n, len(methods)) + 3.0)),
                                     gridspec_kw={"width_ratios": [1.5, 1], "wspace": .30})
        a1.imshow(M, cmap="coolwarm", aspect="auto", vmin=0, vmax=1)
        a1.axhline(n_cls - .5, color="black", lw=2.2)      # 고전 / 양자 경계
        for r_, v in enumerate(row_ag):                     # 행별 정답 일치율
            a1.text(n - .35, r_, f"{v:.0f}%", ha="left", va="center", fontsize=10.5,
                    fontweight="bold", color="#42505C")
        a1.set_yticks(range(len(methods))); a1.set_yticklabels(methods, fontsize=11.5)
        a1.set_xticks(range(n)); a1.set_xticklabels(d["nodes"], rotation=90, fontsize=10)
        a1.grid(False)
        for r_ in range(M.shape[0]):
            for c_ in range(M.shape[1]):
                a1.text(c_, r_, str(int(M[r_, c_])), ha="center", va="center",
                        fontsize=8, color="white", fontweight="bold")
        a1.set_title(f"community assignment per gene — {name}\n"
                     "(aligned to the exact answer;  0 / 1 = the two modules)",
                     fontsize=13.5, fontweight="bold", pad=12)
        o = np.argsort(ag)
        a2.barh(np.arange(n), ag[o], .68,
                color=["#55A868" if v == 100 else ("#DD8452" if v >= 60 else C_B) for v in ag[o]],
                alpha=.88)
        a2.set_yticks(np.arange(n)); a2.set_yticklabels(np.array(d["nodes"])[o], fontsize=11)
        a2.set_xlabel("% of IQM hardware runs agreeing with the exact answer", fontsize=12)
        a2.set_xlim(0, 118); a2.invert_yaxis()
        for i, v in enumerate(ag[o]):
            a2.text(v + 2.5, i, f"{v:.0f}%", va="center", fontsize=10.5, fontweight="bold")
        a2.set_title("per-gene stability across IQM hardware runs", fontsize=13.5,
                     fontweight="bold", pad=12)
        fig.suptitle(f"F   Gene-level community assignment consistency — {name}",
                     fontsize=17, fontweight="bold", y=1.035)
        fig.text(.5, .975, "every hardware repeat is shown as its own row (★ = that run reached "
                 "the exact optimum);  % at the right of each row = agreement with the exact answer",
                 ha="center", fontsize=11.5, color="#5A6672", style="italic")
        save(fig, f"F_gene_consistency_{tag}.png")

    # ================================================================== Fig G
    fig, axes = plt.subplots(1, 2, figsize=(17.5, 7.6))
    for ax, (tag, name, n) in zip(axes, NETS):
        d = D[tag]
        methods, mats = [], []
        for k, (e, x, t) in d["sol"].items():
            methods.append(short(k).replace("\n", " ").replace("(brute force)", "")
                           .strip()); mats.append(np.asarray(x))
        # 각 방법을 '분할 집합' 으로 다룬다. 고전은 원소 1개, IQM 은 회차 5개.
        # best-of-5 를 쓰면 16유전자는 두 조건 모두 최적해를 찍어 전 칸이 1.00 이 되어
        # 아무 정보도 주지 못했다 (Fig F 와 같은 체리피킹).
        groups = [[m] for m in mats]
        for p, h in sorted(d["hw"].items()):
            methods.append(f"IQM p={p}\n({h['runs']} runs)")
            groups.append([np.asarray(r["best_bits"]) for r in h["per_run"]])
        methods.append("Random"); groups.append([RNG.integers(0, 2, n) for _ in range(5)])
        k_ = len(groups)
        A = np.zeros((k_, k_)); S = np.zeros((k_, k_))
        for i in range(k_):
            for j in range(k_):
                if i == j:
                    # 같은 조건 안의 서로 다른 회차끼리 — 재현성
                    v = [agree(groups[i][a], groups[i][b])
                         for a in range(len(groups[i])) for b in range(a + 1, len(groups[i]))]
                    v = v or [1.0]
                else:
                    v = [agree(a, b) for a in groups[i] for b in groups[j]]
                A[i, j] = float(np.mean(v))
                S[i, j] = float(np.std(v, ddof=1)) if len(v) > 1 else 0.0
        im = ax.imshow(A, cmap="YlGnBu", vmin=.5, vmax=1)
        ax.set_xticks(range(k_)); ax.set_yticks(range(k_))
        ax.set_xticklabels(methods, fontsize=10.5)
        ax.set_yticklabels(methods, fontsize=10.5)
        ax.grid(False)
        for i in range(k_):
            for j in range(k_):
                col = "white" if A[i, j] > .82 else "#22303C"
                if S[i, j] > 1e-9:
                    ax.text(j, i - .13, f"{A[i,j]:.2f}", ha="center", va="center",
                            fontsize=10.5, fontweight="bold", color=col)
                    ax.text(j, i + .20, f"±{S[i,j]:.2f}", ha="center", va="center",
                            fontsize=8.5, color=col)
                else:
                    ax.text(j, i, f"{A[i,j]:.2f}", ha="center", va="center", fontsize=11,
                            fontweight="bold", color=col)
        ax.set_title(f"{name}  ({n} genes)", fontsize=14, fontweight="bold", pad=12)
        cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.03)
        cb.set_label("partition agreement", fontsize=11)
    fig.suptitle("G   Pairwise agreement between classical and IQM partitions",
                 fontsize=17, fontweight="bold", y=1.04)
    fig.text(.5, .965, "each IQM entry averages all five hardware repeats (mean ± SD); the IQM "
             "diagonal is run-to-run reproducibility, not a self-comparison.   "
             "1.00 = identical partition,  ~0.5 = random-level agreement",
             ha="center", fontsize=11, color="#5A6672", style="italic")
    save(fig, "G_partition_agreement.png")

    # ================================================================== Fig H
    fig, ax = plt.subplots(figsize=(10, 6.8))
    # 게이트 수는 18번이 산출한 gate_counts.json 만 신뢰한다 (하드코딩 금지)
    _gc = json.load(open(f"{R}/gate_counts.json"))
    depth_map = {(v["instance"], v["depth"]): v["transpiled_2q"] for v in _gc.values()}
    # 두 시리즈가 156 게이트에서 겹치므로 라벨을 위/아래로 갈라 놓는다
    for (tag, name, n), col, mk, dy in zip(NETS, [C_A, C_B], ["o", "s"], [-40, 16]):
        d = D[tag]
        xs, ys, es = [], [], []
        for p, h in sorted(d["hw"].items()):
            xs.append(depth_map.get((tag, p), np.nan))
            ys.append(h["ratio_best"]); es.append(h["ratio_best_sd"])
        ax.errorbar(xs, ys, yerr=es, fmt=mk + "-", color=col, ms=11, capsize=5, lw=2,
                    label=f"{name} ({n} genes)")
        for x_, y_, p in zip(xs, ys, sorted(d["hw"])):
            ax.annotate(f"p = {p}\n{int(x_)} gates", (x_, y_), textcoords="offset points",
                        xytext=(14, dy), fontsize=11, color=col, fontweight="bold",
                        ha="left", linespacing=1.4)
    ax.axhline(1.0, color="#55A868", ls="--", lw=1.8, label="exact optimum")
    ax.set_xlabel("two-qubit gates after transpilation to the IQM lattice", fontsize=12.5)
    ax.set_ylabel("approximation ratio (best solution)", fontsize=12.5)
    ax.set_ylim(0.6, 1.10); ax.set_xlim(40, 400)
    ax.legend(fontsize=11, loc="lower left")
    ax.set_title("H   Circuit size versus solution quality on IQM hardware",
                 fontsize=16, fontweight="bold", pad=14)
    save(fig, "H_depth_vs_quality.png")

    # 요약 저장
    summ = {}
    for tag, name, n in NETS:
        d = D[tag]
        summ[tag] = dict(
            name=name, n=n, edges=d["G"].number_of_edges(),
            density=round(nx.density(d["G"]), 4), E_opt=d["eopt"],
            classical={k: dict(ratio=(d["erand"] - v[0]) / d["den"], seconds=v[2])
                       for k, v in d["sol"].items()},
            iqm={f"p{p}": dict(ratio_best=h["ratio_best"], sd=h["ratio_best_sd"],
                               ratio_mean=h["ratio_mean"], hits=h["hits"], runs=h["runs"])
                 for p, h in d["hw"].items()})
    json.dump(summ, open(f"{OUT}/summary.json", "w"), indent=2)
    print(f"\n저장 -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
