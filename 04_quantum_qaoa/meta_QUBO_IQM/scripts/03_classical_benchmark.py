#!/usr/bin/env python3
"""
고전 기준선 — QAOA 를 정직하게 평가하기 위한 참조값을 먼저 확정한다.

  1) 완전탐색(exact)      N <= 24 : 2^N 전수 열거로 진짜 최적해를 구한다.
                                    -> 근사비를 '추정'이 아니라 '사실'로 보고할 수 있다.
  2) Simulated Annealing  다중 재시작
  3) Genetic Algorithm    다중 재시작
  4) 무작위 / 그리디      바닥 기준선

에너지 정의는 모든 방법에서 동일하다:  E(x) = x^T Q x  (최소화)
근사비는 무작위 기대값을 0, 최적해를 1 로 정규화한 값으로 보고한다.
    ratio = (E_random_mean - E) / (E_random_mean - E_opt)
단순 E/E_opt 는 offset 에 따라 부풀려지므로 쓰지 않는다.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
OUT = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")
EXACT_MAX_N = 24
SA_RESTARTS, SA_STEPS = 20, 20_000
GA_RESTARTS, GA_GEN, GA_POP = 20, 300, 60
RNG = np.random.default_rng(7)


def energy(Q: np.ndarray, x: np.ndarray) -> float:
    return float(x @ Q @ x)


def energies_batch(Q: np.ndarray, X: np.ndarray) -> np.ndarray:
    """X: (m, n) 0/1 행렬 -> (m,) 에너지"""
    return np.einsum("ij,jk,ik->i", X, Q, X)


# ------------------------------------------------------------------ 완전탐색
def exact_solve(Q: np.ndarray, chunk: int = 1 << 20) -> tuple[float, np.ndarray, float]:
    n = Q.shape[0]
    total = 1 << n
    best_e, best_x = np.inf, None
    bits = np.arange(n, dtype=np.uint64)
    t0 = time.time()
    for s in range(0, total, chunk):
        idx = np.arange(s, min(s + chunk, total), dtype=np.uint64)
        X = ((idx[:, None] >> bits) & 1).astype(np.float64)
        E = energies_batch(Q, X)
        k = int(np.argmin(E))
        if E[k] < best_e:
            best_e, best_x = float(E[k]), X[k].copy()
    return best_e, best_x.astype(int), time.time() - t0


# --------------------------------------------------------------------- SA
def simulated_annealing(Q: np.ndarray, steps: int, seed: int) -> tuple[float, np.ndarray]:
    rng = np.random.default_rng(seed)
    n = Q.shape[0]
    x = rng.integers(0, 2, n).astype(np.float64)
    e = energy(Q, x)
    best_e, best_x = e, x.copy()
    # 이웃 뒤집기의 에너지 변화: dE = (1-2x_i)(2*sum_j Q_ij x_j - 2 Q_ii x_i + Q_ii)
    T0, T1 = max(1.0, np.abs(Q).max()), 1e-3
    for t in range(steps):
        T = T0 * (T1 / T0) ** (t / steps)
        i = rng.integers(n)
        xi = x[i]
        s = 1.0 - 2.0 * xi
        dE = s * (2.0 * (Q[i] @ x) - 2.0 * Q[i, i] * xi + Q[i, i] * s * 0 + Q[i, i])
        # 정확한 델타를 위해 직접 계산 (n 이 작아 비용 무시 가능)
        x[i] = 1 - xi
        e_new = energy(Q, x)
        dE = e_new - e
        if dE <= 0 or rng.random() < np.exp(-dE / max(T, 1e-12)):
            e = e_new
            if e < best_e:
                best_e, best_x = e, x.copy()
        else:
            x[i] = xi
    return best_e, best_x.astype(int)


# --------------------------------------------------------------------- GA
def genetic_algorithm(Q: np.ndarray, gens: int, pop_size: int, seed: int
                      ) -> tuple[float, np.ndarray]:
    rng = np.random.default_rng(seed)
    n = Q.shape[0]
    pop = rng.integers(0, 2, (pop_size, n)).astype(np.float64)
    for _ in range(gens):
        E = energies_batch(Q, pop)
        order = np.argsort(E)
        pop = pop[order]
        elite = pop[: pop_size // 4]
        children = []
        while len(children) < pop_size - len(elite):
            a, b = elite[rng.integers(len(elite))], elite[rng.integers(len(elite))]
            cut = rng.integers(1, n)
            c = np.concatenate([a[:cut], b[cut:]])
            mask = rng.random(n) < 0.05
            c[mask] = 1 - c[mask]
            children.append(c)
        pop = np.vstack([elite, np.array(children)])
    E = energies_batch(Q, pop)
    k = int(np.argmin(E))
    return float(E[k]), pop[k].astype(int)


# ------------------------------------------------------------------ 기준선
def random_baseline(Q: np.ndarray, m: int = 20000) -> tuple[float, float]:
    X = RNG.integers(0, 2, (m, Q.shape[0])).astype(np.float64)
    E = energies_batch(Q, X)
    return float(E.mean()), float(E.min())


def greedy(Q: np.ndarray) -> float:
    n = Q.shape[0]
    x = np.zeros(n)
    improved = True
    while improved:
        improved = False
        base = energy(Q, x)
        for i in range(n):
            x[i] = 1 - x[i]
            if energy(Q, x) < base - 1e-12:
                base = energy(Q, x)
                improved = True
            else:
                x[i] = 1 - x[i]
    return energy(Q, x)


# --------------------------------------------------------------------------- #
def main() -> int:
    inst = pd.read_csv(os.path.join(OUT, "qubo_instances.tsv"), sep="\t")
    print("=" * 88)
    print(" 고전 기준선 벤치마크")
    print("=" * 88)

    rows = []
    for _, r in inst.iterrows():
        tag = r.instance
        Q = np.load(os.path.join(OUT, f"{tag}_Q.npy"))
        n = Q.shape[0]
        print(f"\n---- {tag}  (N={n}, {r['kind']}) ----")

        e_rand_mean, e_rand_best = random_baseline(Q)
        e_greedy = greedy(Q)

        # 완전탐색
        if n <= EXACT_MAX_N:
            e_opt, x_opt, t_ex = exact_solve(Q)
            ref, ref_src = e_opt, "exact"
            print(f"  완전탐색   E = {e_opt:11.4f}   ({2**n:,} 상태, {t_ex:.1f}초)")
        else:
            e_opt, t_ex, ref_src = None, None, "best-known"
            print(f"  완전탐색   불가 (2^{n} = {2**n:.3g} 상태)")

        # SA
        t0 = time.time()
        sa = [simulated_annealing(Q, SA_STEPS, 100 + i) for i in range(SA_RESTARTS)]
        t_sa = time.time() - t0
        e_sa = min(s[0] for s in sa)
        sa_all = np.array([s[0] for s in sa])

        # GA
        t0 = time.time()
        ga = [genetic_algorithm(Q, GA_GEN, GA_POP, 200 + i) for i in range(GA_RESTARTS)]
        t_ga = time.time() - t0
        e_ga = min(g[0] for g in ga)
        ga_all = np.array([g[0] for g in ga])

        if ref_src == "best-known":
            ref = min(e_sa, e_ga, e_greedy)
            print(f"  참조해     E = {ref:11.4f}   (SA/GA 최고값, 최적 보장 없음)")

        def ratio(e: float) -> float:
            d = e_rand_mean - ref
            return float((e_rand_mean - e) / d) if abs(d) > 1e-12 else np.nan

        print(f"  무작위평균 E = {e_rand_mean:11.4f}   ratio {ratio(e_rand_mean):.4f}")
        print(f"  그리디     E = {e_greedy:11.4f}   ratio {ratio(e_greedy):.4f}")
        print(f"  SA  (x{SA_RESTARTS})  E = {e_sa:11.4f}   ratio {ratio(e_sa):.4f}   "
              f"최적도달 {int((sa_all <= ref + 1e-9).sum())}/{SA_RESTARTS}   {t_sa:.1f}초")
        print(f"  GA  (x{GA_RESTARTS})  E = {e_ga:11.4f}   ratio {ratio(e_ga):.4f}   "
              f"최적도달 {int((ga_all <= ref + 1e-9).sum())}/{GA_RESTARTS}   {t_ga:.1f}초")

        rec = dict(instance=tag, kind=r["kind"], n_qubits=n, n_edges=int(r.n_edges),
                   ref_source=ref_src, E_ref=ref, E_random_mean=e_rand_mean,
                   E_greedy=e_greedy, E_SA=e_sa, E_GA=e_ga,
                   ratio_random=ratio(e_rand_mean), ratio_greedy=ratio(e_greedy),
                   ratio_SA=ratio(e_sa), ratio_GA=ratio(e_ga),
                   SA_hit=int((sa_all <= ref + 1e-9).sum()), SA_restarts=SA_RESTARTS,
                   GA_hit=int((ga_all <= ref + 1e-9).sum()), GA_restarts=GA_RESTARTS,
                   t_exact=t_ex, t_SA=t_sa, t_GA=t_ga)
        rows.append(rec)

        payload = {"instance": tag, "E_ref": ref, "ref_source": ref_src,
                   "E_random_mean": e_rand_mean, "E_random_best": e_rand_best,
                   "E_greedy": e_greedy, "E_SA": e_sa, "E_GA": e_ga}
        if ref_src == "exact":
            payload["x_opt"] = x_opt.tolist()
        with open(os.path.join(OUT, f"{tag}_classical.json"), "w") as fh:
            json.dump(payload, fh, indent=2)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "classical_benchmark.tsv"), sep="\t", index=False)
    print("\n" + "=" * 88)
    print(" 요약 (ratio: 무작위=0, 참조해=1)")
    print("=" * 88)
    print(df[["instance", "n_qubits", "ref_source", "E_ref", "ratio_greedy",
              "ratio_SA", "ratio_GA", "SA_hit", "GA_hit"]].to_string(index=False))
    print(f"\n저장 -> {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
