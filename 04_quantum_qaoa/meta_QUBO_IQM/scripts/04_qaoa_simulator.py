#!/usr/bin/env python3
"""
QAOA 시뮬레이터 벤치마크 + IQM 토폴로지 자원 추정

(A) 무잡음 시뮬레이션 (Aer statevector)
      p = 1,2,3 에 대해 COBYLA 로 각도 최적화 -> 근사비, 최적해 확률
      완전탐색 최적해를 알고 있으므로 근사비를 정확히 계산한다.

(B) IQM 토폴로지 자원 추정
      Garnet 20q / Emerald 54q 는 정사각 격자형 초전도 칩이다.
      임의 연결 회로를 격자에 올리면 SWAP 이 삽입되어 깊이가 늘어난다.
      전결합(all-to-all) 가정 대비 실제 격자에서의 2큐빗 게이트 수·깊이를 비교해
      "이 문제를 실기에 올릴 수 있는가"를 정량화한다.
      * 실제 칩 커플링맵이 아니라 동일 큐빗 수의 격자 근사임을 명시한다.

(C) 잡음 모델 시뮬레이션
      2큐빗 게이트 오류율을 바꿔가며 근사비가 무너지는 지점을 찾는다.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import QAOAAnsatz
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit.transpiler import CouplingMap
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
OUT = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")
SIM_MAX_N = 22            # statevector 로 전수 기대값을 계산할 상한
P_LAYERS = [1, 2, 3]
MAXITER = 300
SEED = 11
RNG = np.random.default_rng(SEED)
# [deposit] optional QUICK=1 for a fast smoke test (60 instead of 300 COBYLA iterations; p=3 is kept
# because 05 needs it); default unchanged. Full run: ~30 min originally, 87 min on a loaded 20-core
# machine (string_maxcut_19 p=3 dominates).
if os.environ.get("QUICK") == "1":
    MAXITER = 60


# --------------------------------------------------------------- QUBO -> Ising
def qubo_to_ising(Q: np.ndarray) -> tuple[SparsePauliOp, float]:
    """E(x)=x^T Q x  (x in {0,1})  ->  Ising SparsePauliOp + 상수항."""
    n = Q.shape[0]
    h = np.zeros(n)
    J = np.zeros((n, n))
    const = 0.0
    for i in range(n):
        const += Q[i, i] / 2.0
        h[i] += -Q[i, i] / 2.0
        for j in range(i + 1, n):
            w = Q[i, j] + Q[j, i]
            if w == 0:
                continue
            const += w / 4.0
            h[i] += -w / 4.0
            h[j] += -w / 4.0
            J[i, j] += w / 4.0
    terms = []
    for i in range(n):
        if abs(h[i]) > 1e-12:
            p = ["I"] * n
            p[n - 1 - i] = "Z"
            terms.append(("".join(p), h[i]))
    for i in range(n):
        for j in range(i + 1, n):
            if abs(J[i, j]) > 1e-12:
                p = ["I"] * n
                p[n - 1 - i] = "Z"
                p[n - 1 - j] = "Z"
                terms.append(("".join(p), J[i, j]))
    return SparsePauliOp.from_list(terms), const


def all_energies(Q: np.ndarray) -> np.ndarray:
    n = Q.shape[0]
    idx = np.arange(1 << n, dtype=np.uint64)
    bits = np.arange(n, dtype=np.uint64)
    X = ((idx[:, None] >> bits) & 1).astype(np.float64)
    return np.einsum("ij,jk,ik->i", X, Q, X)


# ----------------------------------------------------------------------- QAOA
def run_qaoa(Q: np.ndarray, p: int, E_all: np.ndarray, e_ref: float,
             e_rand: float) -> dict:
    op, const = qubo_to_ising(Q)
    n = Q.shape[0]
    ansatz = QAOAAnsatz(op, reps=p)
    ansatz = ansatz.decompose(reps=3)

    nfev = {"n": 0}

    def expectation(theta: np.ndarray) -> float:
        qc = ansatz.assign_parameters(theta)
        sv = Statevector.from_instruction(qc)
        probs = np.abs(sv.data) ** 2
        nfev["n"] += 1
        return float(probs @ E_all)

    best = None
    for trial in range(3):                       # 다중 시작점
        x0 = RNG.uniform(0, np.pi, ansatz.num_parameters)
        res = minimize(expectation, x0, method="COBYLA",
                       options={"maxiter": MAXITER, "rhobeg": 0.5})
        if best is None or res.fun < best.fun:
            best = res

    qc = ansatz.assign_parameters(best.x)
    probs = np.abs(Statevector.from_instruction(qc).data) ** 2
    e_qaoa = float(probs @ E_all)
    p_opt = float(probs[np.isclose(E_all, E_all.min(), atol=1e-9)].sum())
    # 샘플링해서 최선을 고르는 실전 방식 (1024 shots 상당)
    k = min(1024, len(probs))
    top = np.argsort(-probs)[:k]
    e_best_sampled = float(E_all[top].min())

    denom = e_rand - e_ref
    return dict(p=p, E_expect=e_qaoa, E_best_sampled=e_best_sampled,
                ratio_expect=float((e_rand - e_qaoa) / denom),
                ratio_best=float((e_rand - e_best_sampled) / denom),
                P_optimal=p_opt, n_params=ansatz.num_parameters,
                n_fev=nfev["n"], depth_ideal=ansatz.depth(),
                cx_ideal=ansatz.count_ops().get("cx", 0)
                + ansatz.count_ops().get("rzz", 0))


# --------------------------------------------------- IQM 격자 토폴로지 자원추정
def grid_coupling(n_qubits: int) -> tuple[CouplingMap, tuple[int, int]]:
    """n 큐빗을 담는 가장 정사각에 가까운 격자."""
    r = int(np.floor(np.sqrt(n_qubits)))
    while r > 1 and n_qubits % r:
        r -= 1
    c = n_qubits // r
    return CouplingMap.from_grid(r, c), (r, c)


def resource_estimate(Q: np.ndarray, p: int, device_qubits: int, label: str) -> dict:
    op, _ = qubo_to_ising(Q)
    ansatz = QAOAAnsatz(op, reps=p).decompose(reps=3)
    cm, shape = grid_coupling(device_qubits)
    t0 = time.time()
    tq = transpile(ansatz, coupling_map=cm, basis_gates=["rz", "sx", "x", "cz"],
                   optimization_level=3, seed_transpiler=SEED)
    dt = time.time() - t0
    ops = tq.count_ops()
    ideal = ansatz.count_ops()
    return dict(device=label, grid=f"{shape[0]}x{shape[1]}", p=p,
                depth_alltoall=ansatz.depth(),
                twoq_alltoall=ideal.get("rzz", 0) + ideal.get("cx", 0),
                depth_grid=tq.depth(), twoq_grid=ops.get("cz", 0),
                swap_overhead=round(ops.get("cz", 0) /
                                    max(1, ideal.get("rzz", 0) + ideal.get("cx", 0)), 2),
                transpile_sec=round(dt, 1))


# ------------------------------------------------------------------- 잡음
def noisy_ratio(Q: np.ndarray, theta: np.ndarray, p: int, err2q: float,
                E_all: np.ndarray, e_ref: float, e_rand: float,
                shots: int = 8192) -> dict:
    op, _ = qubo_to_ising(Q)
    ansatz = QAOAAnsatz(op, reps=p).decompose(reps=3)
    qc = ansatz.assign_parameters(theta)
    qc.measure_all()
    nm = NoiseModel()
    nm.add_all_qubit_quantum_error(depolarizing_error(err2q, 2), ["cx", "cz"])
    nm.add_all_qubit_quantum_error(depolarizing_error(err2q / 10, 1), ["sx", "x", "rz"])
    sim = AerSimulator(noise_model=nm, basis_gates=["rz", "sx", "x", "cz"])
    tq = transpile(qc, sim, optimization_level=1, seed_transpiler=SEED)
    counts = sim.run(tq, shots=shots, seed_simulator=SEED).result().get_counts()
    n = Q.shape[0]
    tot, esum, ebest = 0, 0.0, np.inf
    for bits, c in counts.items():
        b = bits.replace(" ", "")[::-1][:n]
        i = int(b[::-1], 2)
        e = E_all[i]
        esum += e * c
        tot += c
        ebest = min(ebest, e)
    e_mean = esum / tot
    denom = e_rand - e_ref
    return dict(err2q=err2q, E_mean=e_mean,
                ratio_mean=float((e_rand - e_mean) / denom),
                ratio_best=float((e_rand - ebest) / denom))


# --------------------------------------------------------------------------- #
def main() -> int:
    inst = pd.read_csv(os.path.join(OUT, "qubo_instances.tsv"), sep="\t")
    print("=" * 90)
    print(" QAOA 시뮬레이터 벤치마크 + IQM 토폴로지 자원 추정")
    print("=" * 90)

    qaoa_rows, res_rows, noise_rows = [], [], []

    for _, r in inst.iterrows():
        tag = r.instance
        Q = np.load(os.path.join(OUT, f"{tag}_Q.npy"))
        n = Q.shape[0]
        cls = json.load(open(os.path.join(OUT, f"{tag}_classical.json")))
        e_ref, e_rand = cls["E_ref"], cls["E_random_mean"]

        print(f"\n{'='*90}\n {tag}  (N={n}, {r['kind']})   참조해 E={e_ref:.4f} ({cls['ref_source']})")

        # ---- (B) 자원 추정은 모든 크기에서 가능
        print(f"\n  [자원 추정] IQM 격자 토폴로지 (동일 큐빗 수 격자 근사)")
        for dev, dq in [("IQM Garnet", 20), ("IQM Emerald", 54)]:
            if n > dq:
                print(f"    {dev:12s} ({dq}q) : N={n} > {dq} — 올릴 수 없음")
                continue
            for p in [1, 2]:
                est = resource_estimate(Q, p, dq, dev)
                est["instance"] = tag
                est["n_qubits"] = n
                res_rows.append(est)
                print(f"    {dev:12s} p={p}  격자 {est['grid']:5s}  "
                      f"2q: 전결합 {est['twoq_alltoall']:4d} -> 격자 {est['twoq_grid']:5d} "
                      f"(x{est['swap_overhead']:.2f})   깊이 {est['depth_alltoall']:4d} -> {est['depth_grid']:5d}")

        # ---- (A) 무잡음 시뮬레이션
        if n > SIM_MAX_N:
            print(f"\n  [시뮬레이션] N={n} > {SIM_MAX_N} — statevector 불가 "
                  f"(2^{n} 진폭 = {2**n * 16 / 1e9:.0f} GB). 자원추정만 수행.")
            continue

        E_all = all_energies(Q)
        print(f"\n  [무잡음 QAOA]  (ratio: 무작위=0, 최적해=1)")
        best_theta = None
        for p in P_LAYERS:
            t0 = time.time()
            out = run_qaoa(Q, p, E_all, e_ref, e_rand)
            out["instance"] = tag
            out["n_qubits"] = n
            out["sec"] = round(time.time() - t0, 1)
            qaoa_rows.append(out)
            print(f"    p={p}  기대값 E={out['E_expect']:9.4f}  ratio={out['ratio_expect']:.4f}   "
                  f"최적해확률 {out['P_optimal']*100:6.3f}%   "
                  f"샘플최선 ratio={out['ratio_best']:.4f}   "
                  f"2q게이트 {out['cx_ideal']}  {out['sec']}초")

        # ---- (C) 잡음 스캔 (p=1 최적각도 재사용)
        print(f"\n  [잡음 스캔] p=1, 2큐빗 탈분극 오류율별")
        op, _ = qubo_to_ising(Q)
        a1 = QAOAAnsatz(op, reps=1).decompose(reps=3)
        best = None
        for _ in range(3):
            x0 = RNG.uniform(0, np.pi, a1.num_parameters)
            res = minimize(lambda th: float(np.abs(Statevector.from_instruction(
                a1.assign_parameters(th)).data) ** 2 @ E_all),
                x0, method="COBYLA", options={"maxiter": MAXITER})
            if best is None or res.fun < best.fun:
                best = res
        for err in [0.0, 0.001, 0.003, 0.005, 0.01, 0.02]:
            nr = noisy_ratio(Q, best.x, 1, err, E_all, e_ref, e_rand)
            nr["instance"] = tag
            nr["n_qubits"] = n
            noise_rows.append(nr)
            print(f"    2q err={err*100:5.2f}%   평균 ratio={nr['ratio_mean']:.4f}   "
                  f"샘플최선 ratio={nr['ratio_best']:.4f}")

    if qaoa_rows:
        pd.DataFrame(qaoa_rows).to_csv(os.path.join(OUT, "qaoa_noiseless.tsv"),
                                       sep="\t", index=False)
    pd.DataFrame(res_rows).to_csv(os.path.join(OUT, "resource_estimate.tsv"),
                                  sep="\t", index=False)
    if noise_rows:
        pd.DataFrame(noise_rows).to_csv(os.path.join(OUT, "qaoa_noise_scan.tsv"),
                                        sep="\t", index=False)
    print(f"\n저장 -> {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
