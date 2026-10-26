#!/usr/bin/env python3
"""
IQM 실기 실행 계획 수립 + 결과 예측  (토큰 불필요, QPU 사용량 0)

IQM 이 배포하는 fake backend 는 실제 소자의 커플링맵과 오류 프로파일을 담고 있다.
  IQMFakeApollo     20q  (Garnet 급)   엣지 30, 평균차수 3.00, 2q 오류 1.38%
  IQMFakeAphrodite  54q  (Emerald 급)  엣지 90, 평균차수 3.33, 2q 오류 1.50%
04번의 '정사각 격자 근사' 를 이것으로 교체해 실측 기반 자원추정을 한다.

핵심 설계 원칙 — 1분 예산
  하드웨어에서 변분 최적화를 절대 하지 않는다.
  COBYLA 루프는 수백 번의 작업 제출이 필요해 예산을 즉시 초과한다.
  각도는 시뮬레이터에서 미리 최적화하고, 하드웨어에서는 '확정된 회로 1개'만 실행한다.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from qiskit import transpile
from qiskit.circuit.library import QAOAAnsatz
from qiskit.quantum_info import SparsePauliOp, Statevector
from iqm.qiskit_iqm import IQMFakeApollo, IQMFakeAphrodite

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
OUT = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")
RNG = np.random.default_rng(31)

# shot 당 시간 (IQM 공개 사양 기준 보수적 가정)
T_1Q, T_2Q, T_READOUT, T_RESET = 40e-9, 100e-9, 3e-6, 250e-6
BUDGET_S = 60.0                      # 사용 가능한 QPU 시간
# [deposit] optional QUICK=1 smoke test: 2 restarts x 100 COBYLA iterations instead of 8 x 400
# (angle optimisation of the 19-qubit statevector dominates; multi-hour on a loaded machine). Default unchanged.
_QUICK = os.environ.get("QUICK") == "1"


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
    t = []
    for i in range(n):
        if abs(h[i]) > 1e-12:
            p = ["I"] * n; p[n - 1 - i] = "Z"; t.append(("".join(p), h[i]))
    for i in range(n):
        for j in range(i + 1, n):
            if abs(J[i, j]) > 1e-12:
                p = ["I"] * n; p[n - 1 - i] = "Z"; p[n - 1 - j] = "Z"
                t.append(("".join(p), J[i, j]))
    return SparsePauliOp.from_list(t)


def all_energies(Q):
    n = Q.shape[0]
    idx = np.arange(1 << n, dtype=np.uint64)
    bits = np.arange(n, dtype=np.uint64)
    X = ((idx[:, None] >> bits) & 1).astype(np.float64)
    return np.einsum("ij,jk,ik->i", X, Q, X)


def optimise_angles(Q, p, E_all, restarts=2 if _QUICK else 8):  # [deposit] QUICK
    op = qubo_to_ising(Q)
    ans = QAOAAnsatz(op, reps=p).decompose(reps=3)
    f = lambda th: float(np.abs(Statevector.from_instruction(
        ans.assign_parameters(th)).data) ** 2 @ E_all)
    best = None
    for _ in range(restarts):
        r = minimize(f, RNG.uniform(0, np.pi, ans.num_parameters),
                     method="COBYLA", options={"maxiter": 100 if _QUICK else 400, "rhobeg": 0.4})  # [deposit] QUICK
        if best is None or r.fun < best.fun:
            best = r
    return best.x, float(best.fun)


def energies_from_counts(counts, E_all, n):
    tot, esum, ebest = 0, 0.0, np.inf
    for bits, c in counts.items():
        b = bits.replace(" ", "")
        i = int(b[-n:], 2)
        e = E_all[i]
        esum += e * c
        tot += c
        ebest = min(ebest, e)
    return esum / tot, ebest


def main() -> int:
    inst = pd.read_csv(os.path.join(OUT, "qubo_instances.tsv"), sep="\t")
    devices = [("IQMFakeApollo (20q, Garnet급)", IQMFakeApollo(), 20),
               ("IQMFakeAphrodite (54q, Emerald급)", IQMFakeAphrodite(), 54)]

    print("=" * 100)
    print(" IQM 실측 토폴로지 기준 자원추정 + 잡음 예측 + 1분 예산 계획")
    print("=" * 100)

    rows, plan = [], []
    for _, r in inst.iterrows():
        tag, n = r.instance, int(r.n_qubits)
        Q = np.load(os.path.join(OUT, f"{tag}_Q.npy"))
        cls = json.load(open(os.path.join(OUT, f"{tag}_classical.json")))
        e_ref, e_rand = cls["E_ref"], cls["E_random_mean"]
        denom = e_rand - e_ref

        print(f"\n{'='*100}\n {tag}  N={n}  ({r['kind']})   참조해 {e_ref:.4f} [{cls['ref_source']}]")

        E_all = all_energies(Q) if n <= 22 else None
        for dev_name, dev, dq in devices:
            if n > dq:
                print(f"   {dev_name:34s} : N={n} > {dq}q — 실행 불가")
                continue
            for p in (1, 2):
                ans = QAOAAnsatz(qubo_to_ising(Q), reps=p).decompose(reps=3)
                ideal2q = sum(v for k, v in ans.count_ops().items() if k in ("rzz", "cx"))
                ans.measure_all()
                t0 = time.time()
                tq = transpile(ans, backend=dev, optimization_level=3, seed_transpiler=11)
                ops = tq.count_ops()
                n2q = ops.get("cz", 0)
                n1q = ops.get("r", 0)
                depth = tq.depth()
                # 오류율은 fake backend 프로파일에서
                try:
                    ep = dev.error_profile
                    vals = []
                    for d in ep.two_qubit_gate_depolarizing_error_parameters.values():
                        vals += [float(x) for x in (d.values() if hasattr(d, "values") else [d])]
                    err2q = float(np.median(vals))
                except Exception:
                    err2q = 0.014
                fid = (1 - err2q) ** n2q
                t_circ = depth * T_2Q                       # 보수적: 깊이 전부 2q 시간
                t_shot = t_circ + T_READOUT + T_RESET
                rows.append(dict(instance=tag, n_qubits=n, device=dev_name, p=p,
                                 twoq_ideal=ideal2q, twoq_hw=n2q, oneq_hw=n1q,
                                 swap_x=round(n2q / max(1, ideal2q), 2), depth=depth,
                                 err2q=err2q, fidelity=fid,
                                 t_shot_s=t_shot, transpile_s=round(time.time() - t0, 1)))
                print(f"   {dev_name:34s} p={p}  2q {ideal2q:4d}->{n2q:5d} (x{n2q/max(1,ideal2q):.2f})  "
                      f"깊이 {depth:5d}  충실도 {fid:.4f}  shot {t_shot*1e6:6.0f}us")

        # ---- 잡음 시뮬레이션 예측 (fake backend = 실제 오류 프로파일 시뮬레이터)
        if E_all is not None and n <= 20:
            dev = IQMFakeApollo()
            print(f"\n   [예측] IQMFakeApollo 잡음 시뮬레이션 (실기와 같은 오류 프로파일)")
            for p in (1, 2):
                theta, e_ideal = optimise_angles(Q, p, E_all)
                ans = QAOAAnsatz(qubo_to_ising(Q), reps=p).decompose(reps=3)
                qc = ans.assign_parameters(theta)
                qc.measure_all()
                tq = transpile(qc, backend=dev, optimization_level=3, seed_transpiler=11)
                res = dev.run(tq, shots=4096).result()
                e_mean, e_best = energies_from_counts(res.get_counts(), E_all, n)
                print(f"     p={p}  무잡음 ratio {(e_rand-e_ideal)/denom:.4f}  ->  "
                      f"잡음 평균 ratio {(e_rand-e_mean)/denom:.4f}   "
                      f"4096샷 최선 ratio {(e_rand-e_best)/denom:.4f}"
                      f"   {'최적해 도달' if abs(e_best-e_ref)<1e-9 else '최적해 미도달'}")
                plan.append(dict(instance=tag, n_qubits=n, p=p, device="Apollo/Garnet",
                                 theta=list(map(float, theta)),
                                 ratio_noiseless=float((e_rand - e_ideal) / denom),
                                 ratio_noisy_mean=float((e_rand - e_mean) / denom),
                                 ratio_noisy_best=float((e_rand - e_best) / denom),
                                 reached_optimum=bool(abs(e_best - e_ref) < 1e-9)))

    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(OUT, "iqm_resource_real.tsv"), sep="\t", index=False)
    if plan:
        with open(os.path.join(OUT, "iqm_predicted.json"), "w") as fh:
            json.dump(plan, fh, indent=2)

    # ---------------------------------------------------------------- 예산 계획
    print("\n" + "=" * 100)
    print(f" 1분({BUDGET_S:.0f}초) QPU 예산 계획")
    print("=" * 100)
    print("\n  원칙: 하드웨어에서 변분 최적화 금지. 시뮬레이터에서 각도를 확정하고")
    print("        하드웨어에서는 '고정 각도 회로 1개'만 실행한다.")
    print(f"        (COBYLA 루프를 하드웨어로 돌리면 회로 수백 개 -> 예산 즉시 초과)\n")

    cand = R[(R.p <= 2)].copy()
    cand["priority"] = np.where(cand.fidelity > 0.3, 1, np.where(cand.fidelity > 0.05, 2, 3))
    cand = cand.sort_values(["priority", "fidelity"], ascending=[True, False])

    SHOTS = 4096
    total = 0.0
    print(f"  {'순위':<4}{'인스턴스':<20}{'소자':<12}{'p':<3}{'충실도':<9}{'shot시간':<11}{'shots':<8}{'소요':<9}누적")
    print("  " + "-" * 96)
    chosen = []
    for _, c in cand.iterrows():
        dev_short = "Garnet(20q)" if "Apollo" in c.device else "Emerald(54q)"
        t = c.t_shot_s * SHOTS
        if total + t > BUDGET_S:
            continue
        total += t
        chosen.append(dict(instance=c.instance, device=dev_short, p=int(c.p),
                           shots=SHOTS, fidelity=round(float(c.fidelity), 4),
                           twoq=int(c.twoq_hw), depth=int(c.depth),
                           est_seconds=round(float(t), 2)))
        print(f"  {len(chosen):<4}{c.instance:<20}{dev_short:<12}{c.p:<3}"
              f"{c.fidelity:<9.4f}{c.t_shot_s*1e6:<11.0f}{SHOTS:<8}{t:<9.2f}{total:.2f}s")
    print(f"\n  총 예상 QPU 시간: {total:.1f}초 / {BUDGET_S:.0f}초  "
          f"(여유 {BUDGET_S-total:.1f}초)")
    print(f"  실행 회로 수: {len(chosen)}개")

    with open(os.path.join(OUT, "iqm_run_plan.json"), "w") as fh:
        json.dump({"budget_s": BUDGET_S, "shots": SHOTS,
                   "estimated_total_s": total, "jobs": chosen,
                   "timing_assumptions": {"t_1q_s": T_1Q, "t_2q_s": T_2Q,
                                          "t_readout_s": T_READOUT, "t_reset_s": T_RESET}},
                  fh, indent=2)
    print(f"\n저장 -> {OUT}/iqm_resource_real.tsv, iqm_run_plan.json, iqm_predicted.json")
    print("\n※ QPU 실행 없음. fake backend(로컬 시뮬레이터)만 사용했습니다.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
