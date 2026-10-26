#!/usr/bin/env python3
"""
IQM 실기 제출 스크립트

사용법
    # 1) 반드시 먼저: fake backend 로 전 과정 검증 (QPU 사용량 0)
    python 08_iqm_submit.py --dry-run

    # 2) 실제 제출 (--yes 없으면 제출하지 않음)
    source ~/.iqm_emerald.sh
    python 08_iqm_submit.py --yes                    # 전체 계획
    python 08_iqm_submit.py --yes --only 1,2,3       # 일부만

안전장치
  * --yes 를 주지 않으면 계획만 출력하고 종료한다.
  * 제출 전 예상 QPU 시간을 출력하고 예산(BUDGET_S)을 넘으면 중단한다.
  * 하드웨어에서 변분 최적화를 하지 않는다. 각도는 07번이 시뮬레이터로 확정한 값.
  * 각 작업 결과를 즉시 저장하므로 중간에 끊겨도 앞선 결과는 남는다.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

from qiskit import transpile
from qiskit.circuit.library import QAOAAnsatz
from qiskit.quantum_info import SparsePauliOp

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
OUT = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")
BUDGET_S = 60.0
SHOTS = 4096


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


def energy_of_bitstring(Q, bits: str) -> float:
    n = Q.shape[0]
    x = np.array([int(c) for c in bits[-n:][::-1]], dtype=float)
    return float(x @ Q @ x)


def analyse(counts: dict, Q: np.ndarray, e_ref: float, e_rand: float) -> dict:
    tot = sum(counts.values())
    esum, ebest, bbest = 0.0, np.inf, None
    for b, c in counts.items():
        e = energy_of_bitstring(Q, b.replace(" ", ""))
        esum += e * c
        if e < ebest:
            ebest, bbest = e, b
    denom = e_rand - e_ref
    return dict(shots=tot, E_mean=esum / tot, E_best=ebest, best_bitstring=bbest,
                ratio_mean=(e_rand - esum / tot) / denom,
                ratio_best=(e_rand - ebest) / denom,
                reached_optimum=bool(abs(ebest - e_ref) < 1e-9))


def get_angles(instance: str, p: int, pred: list) -> tuple[list[float], str]:
    for r in pred:
        if r["instance"] == instance and r["p"] == p:
            return r["theta"], "simulator-optimised"
    # N>22 는 statevector 최적화가 불가능하므로 같은 p 의 다른 인스턴스 각도를 전이한다
    for r in pred:
        if r["p"] == p:
            return r["theta"], f"transferred from {r['instance']}"
    raise SystemExit(f"각도를 찾을 수 없습니다: {instance} p={p}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--yes", action="store_true", help="실제 QPU 제출 (없으면 계획만 출력)")
    ap.add_argument("--dry-run", action="store_true", help="fake backend 로 검증 (QPU 사용 0)")
    ap.add_argument("--only", type=str, default="", help="실행할 작업 번호 (예: 1,2,4)")
    ap.add_argument("--shots", type=int, default=SHOTS)
    args = ap.parse_args()

    plan = json.load(open(os.path.join(OUT, "iqm_run_plan.json")))
    pred = json.load(open(os.path.join(OUT, "iqm_predicted.json")))
    jobs = plan["jobs"]
    if args.only:
        pick = {int(x) for x in args.only.split(",")}
        jobs = [j for i, j in enumerate(jobs, 1) if i in pick]

    est = sum(j["est_seconds"] for j in jobs) * (args.shots / plan["shots"])
    print("=" * 96)
    print(f" IQM 제출 계획 — 작업 {len(jobs)}개, shots {args.shots}, 예상 QPU {est:.1f}초 / {BUDGET_S:.0f}초")
    print("=" * 96)
    print(f"  {'#':<3}{'인스턴스':<20}{'소자':<14}{'p':<3}{'예상충실도':<11}{'예상초':<8}")
    for i, j in enumerate(jobs, 1):
        print(f"  {i:<3}{j['instance']:<20}{j['device']:<14}{j['p']:<3}{j['fidelity']:<11.4f}"
              f"{j['est_seconds']*args.shots/plan['shots']:<8.2f}")
    if est > BUDGET_S:
        sys.exit(f"\n중단: 예상 QPU 시간 {est:.1f}초 > 예산 {BUDGET_S:.0f}초. --only 로 줄이세요.")
    print(f"\n  예산 여유: {BUDGET_S - est:.1f}초")

    if not args.yes and not args.dry_run:
        print("\n계획만 출력했습니다. 실제 제출하려면 --yes, 검증하려면 --dry-run 을 주세요.")
        return 0

    # ------------------------------------------------------------- 백엔드 준비
    # [deposit] credential-less fallback: --yes without IQM_SERVER_URL/IQM_TOKEN runs the same jobs on
    # IQM fake backends (local Aer noise simulation) and writes *_SIMULATOR_FALLBACK files instead of
    # iqm_hw_* / iqm_hardware_summary.tsv, so real-hardware outputs are never overwritten.
    fallback = (not args.dry_run) and not (os.environ.get("IQM_SERVER_URL") and os.environ.get("IQM_TOKEN"))
    if fallback:
        print("\n[deposit] no IQM_TOKEN → running on local simulator (IQMFakeApollo / IQMFakeAphrodite); "
              "outputs *_SIMULATOR_FALLBACK")
    if args.dry_run or fallback:
        from iqm.qiskit_iqm import IQMFakeApollo, IQMFakeAphrodite
        backends = {"Garnet(20q)": IQMFakeApollo(), "Emerald(54q)": IQMFakeAphrodite()}
        mode = ("SIMULATOR FALLBACK (fake backend, QPU 사용량 0)" if fallback  # [deposit]
                else "DRY-RUN (fake backend, QPU 사용량 0)")
    else:
        from iqm.qiskit_iqm import IQMProvider
        url = os.environ.get("IQM_SERVER_URL")
        if not url or not os.environ.get("IQM_TOKEN"):
            sys.exit("IQM_SERVER_URL / IQM_TOKEN 이 없습니다. source ~/.iqm_emerald.sh")
        try:
            bk = IQMProvider(url).get_backend()
        except Exception as exc:                                    # noqa: BLE001
            sys.exit(f"백엔드 연결 실패: {type(exc).__name__}: {exc}\n"
                     f"  -> 먼저 python scripts/10_iqm_connect_test.py 로 진단하세요.")
        print(f"\n  연결됨: {bk.name} ({bk.num_qubits} 큐빗)")
        backends = {"Garnet(20q)": bk, "Emerald(54q)": bk}
        mode = f"실기 제출 -> {bk.name}"
        # 서버가 소자 하나만 제공하므로, 계획상 Garnet/Emerald 로 나뉜 같은 회로는
        # 실제로는 동일한 작업이 된다. (instance, p) 기준으로 중복을 제거한다.
        seen, uniq = set(), []
        for j in jobs:
            key = (j["instance"], j["p"])
            if key in seen:
                print(f"  중복 제거: {j['instance']} p={j['p']} ({j['device']}) — "
                      f"서버 소자가 하나뿐")
                continue
            seen.add(key)
            j = dict(j); j["device"] = "Garnet(20q)"
            uniq.append(j)
        jobs = uniq
        est = sum(x["est_seconds"] for x in jobs) * (args.shots / plan["shots"])
        print(f"  중복 제거 후 작업 {len(jobs)}개, 예상 QPU {est:.1f}초")
    print(f"\n모드: {mode}\n" + "=" * 96)

    results = []
    t_start = time.time()
    for i, j in enumerate(jobs, 1):
        tag, p, dev_key = j["instance"], j["p"], j["device"]
        Q = np.load(os.path.join(OUT, f"{tag}_Q.npy"))
        cls = json.load(open(os.path.join(OUT, f"{tag}_classical.json")))
        e_ref, e_rand = cls["E_ref"], cls["E_random_mean"]
        theta, src = get_angles(tag, p, pred)
        backend = backends[dev_key]
        if Q.shape[0] > backend.num_qubits:
            print(f"[{i}/{len(jobs)}] {tag} p={p} — N={Q.shape[0]} > "
                  f"{backend.num_qubits} 큐빗, 이 소자에서는 실행 불가. 건너뜀")
            continue

        ans = QAOAAnsatz(qubo_to_ising(Q), reps=p).decompose(reps=3)
        qc = ans.assign_parameters(theta)
        qc.measure_all()
        tq = transpile(qc, backend=backend, optimization_level=3, seed_transpiler=11)
        n2q = tq.count_ops().get("cz", 0)

        print(f"\n[{i}/{len(jobs)}] {tag}  p={p}  {dev_key}  각도출처={src}")
        print(f"      2q게이트 {n2q}, 깊이 {tq.depth()}, shots {args.shots}")
        t0 = time.time()
        try:
            counts = backend.run(tq, shots=args.shots).result().get_counts()
        except Exception as exc:                                    # noqa: BLE001
            print(f"      실패: {type(exc).__name__}: {str(exc)[:160]}")
            print("      -> 다음 작업으로 넘어갑니다 (앞선 결과는 저장되어 있습니다)")
            continue
        dt = time.time() - t0

        a = analyse(counts, Q, e_ref, e_rand)
        pm = next((r for r in pred if r["instance"] == tag and r["p"] == p), None)
        rec = dict(job=i, instance=tag, n_qubits=int(Q.shape[0]), p=p, device=dev_key,
                   twoq=n2q, depth=tq.depth(), wall_s=round(dt, 2),
                   angle_source=src, mode="dry-run" if args.dry_run else ("simulator-fallback" if fallback else "hardware"), **a)  # [deposit]
        if pm:
            rec["pred_ratio_mean"] = pm["ratio_noisy_mean"]
            rec["pred_ratio_best"] = pm["ratio_noisy_best"]
        results.append(rec)

        print(f"      평균 ratio {a['ratio_mean']:.4f}"
              + (f"  (예측 {pm['ratio_noisy_mean']:.4f})" if pm else ""))
        print(f"      최선 ratio {a['ratio_best']:.4f}"
              + (f"  (예측 {pm['ratio_noisy_best']:.4f})" if pm else "")
              + f"   최적해 {'도달 ✅' if a['reached_optimum'] else '미도달'}   {dt:.1f}초")

        suffix = "dryrun" if args.dry_run else "hw"
        fbs = "_SIMULATOR_FALLBACK" if fallback else ""  # [deposit]
        with open(os.path.join(OUT, f"iqm_{suffix}_{tag}_p{p}_{dev_key.split('(')[0]}{fbs}.json"), "w") as fh:
            json.dump({**rec, "counts": counts}, fh, indent=2)

    df = pd.DataFrame(results)
    fn = "iqm_dryrun_summary.tsv" if args.dry_run else "iqm_hardware_summary.tsv"
    if fallback:
        fn = "iqm_hardware_summary_SIMULATOR_FALLBACK.tsv"  # [deposit]
    df.to_csv(os.path.join(OUT, fn), sep="\t", index=False)

    print("\n" + "=" * 96)
    print(f" 완료 — 총 소요 {time.time()-t_start:.1f}초 (벽시계)")
    print("=" * 96)
    cols = ["job", "instance", "p", "device", "twoq", "ratio_mean", "ratio_best", "reached_optimum"]
    if "pred_ratio_mean" in df:
        cols.insert(6, "pred_ratio_mean")
    print(df[cols].to_string(index=False))
    print(f"\n저장 -> {OUT}/{fn}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
