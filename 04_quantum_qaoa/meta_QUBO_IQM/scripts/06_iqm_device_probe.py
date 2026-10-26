#!/usr/bin/env python3
"""
IQM 소자 정보 조회 (읽기 전용 — QPU 실행 없음, 큐 사용량 0)

목적
  1) 실제 커플링맵을 받아 04번의 '격자 근사' 자원추정을 실측으로 교체
  2) 게이트 시간 / 오류율로 shot 당 실행시간을 계산 -> 1분 예산 배분
  3) 각 인스턴스가 실기에 올라갈 수 있는지 최종 판정
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
OUT = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")


def main() -> int:
    from iqm.qiskit_iqm import IQMProvider

    url = os.environ.get("IQM_SERVER_URL")
    token = os.environ.get("IQM_TOKEN")
    # [deposit] credential-less fallback: probe IQM's fake Emerald-class backend (54q) instead of the
    # live device; outputs get the suffix _SIMULATOR_FALLBACK so real-device files are never overwritten.
    fb = "" if (url and token) else "_SIMULATOR_FALLBACK"
    if fb:
        from iqm.qiskit_iqm import IQMFakeAphrodite
        print("[deposit] no IQM_TOKEN/IQM_SERVER_URL → running on local simulator "
              "(IQMFakeAphrodite, 54q Emerald-class fake backend); outputs *_SIMULATOR_FALLBACK")
        backend = IQMFakeAphrodite()
    else:
        print(f"서버: {url}")
        provider = IQMProvider(url)
        backend = provider.get_backend()
    print(f"백엔드: {backend.name}")
    print(f"큐빗 수: {backend.num_qubits}")

    cm = backend.coupling_map
    edges = sorted({tuple(sorted(e)) for e in cm.get_edges()})
    deg = np.zeros(backend.num_qubits, int)
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    print(f"커플링: 엣지 {len(edges)}개, 평균차수 {deg.mean():.2f}, 최대차수 {deg.max()}")
    print(f"기본 게이트: {backend.operation_names}")

    info = dict(backend=backend.name, num_qubits=int(backend.num_qubits),
                n_edges=len(edges), mean_degree=float(deg.mean()),
                max_degree=int(deg.max()), edges=[list(e) for e in edges],
                basis_gates=list(backend.operation_names))

    # ---- 보정 데이터 (있으면)
    for attr in ("target", "_calibration_set_id"):
        if hasattr(backend, attr):
            print(f"  {attr}: 사용 가능")
    try:
        tgt = backend.target
        durs, errs = [], []
        for inst, props in tgt.items():
            for qargs, p in (props or {}).items():
                if p is None:
                    continue
                if len(qargs or ()) == 2 and getattr(p, "duration", None):
                    durs.append(p.duration)
                    if getattr(p, "error", None) is not None:
                        errs.append(p.error)
        if durs:
            info["twoq_duration_s"] = float(np.median(durs))
            print(f"  2큐빗 게이트 시간(중앙값): {np.median(durs)*1e9:.0f} ns")
        if errs:
            info["twoq_error"] = float(np.median(errs))
            print(f"  2큐빗 오류율(중앙값): {np.median(errs)*100:.3f} %")
    except Exception as exc:                      # noqa: BLE001
        print(f"  보정 데이터 조회 실패(무시): {exc}")

    with open(os.path.join(OUT, f"iqm_device_info{fb}.json"), "w") as fh:  # [deposit] fb suffix
        json.dump(info, fh, indent=2)

    # ---- 실제 커플링맵으로 재-transpile
    from qiskit import transpile
    from qiskit.circuit.library import QAOAAnsatz
    from qiskit.quantum_info import SparsePauliOp

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

    inst = pd.read_csv(os.path.join(OUT, "qubo_instances.tsv"), sep="\t")
    twoq_dur = info.get("twoq_duration_s", 100e-9)
    err2q = info.get("twoq_error", 0.005)
    READOUT_S, RESET_S = 2e-6, 200e-6           # 보수적 가정 (IQM 문서 범위)

    rows = []
    print("\n" + "=" * 92)
    print(" 실제 커플링맵 기준 자원 추정 + shot 시간")
    print("=" * 92)
    for _, r in inst.iterrows():
        n = int(r.n_qubits)
        if n > backend.num_qubits:
            print(f"  {r.instance:20s} N={n:3d}  > {backend.num_qubits}q — 올릴 수 없음")
            continue
        Q = np.load(os.path.join(OUT, f"{r.instance}_Q.npy"))
        for p in (1, 2):
            ans = QAOAAnsatz(qubo_to_ising(Q), reps=p).decompose(reps=3)
            ans.measure_all()
            tq = transpile(ans, backend=backend, optimization_level=3, seed_transpiler=11)
            ops = tq.count_ops()
            n2q = sum(v for k, v in ops.items() if k in ("cz", "move", "czprime"))
            depth = tq.depth()
            circ_s = depth * twoq_dur
            shot_s = circ_s + READOUT_S + RESET_S
            fid = (1 - err2q) ** n2q
            rows.append(dict(instance=r.instance, n_qubits=n, p=p, depth=depth,
                             twoq=n2q, ops=dict(ops), circuit_s=circ_s,
                             shot_s=shot_s, fidelity=fid,
                             shots_per_min=int(60 / shot_s)))
            print(f"  {r.instance:20s} p={p}  깊이 {depth:5d}  2q {n2q:5d}  "
                  f"충실도 {fid:.4f}  shot {shot_s*1e6:7.1f} us  "
                  f"-> 1분에 {int(60/shot_s):,} shots")

    pd.DataFrame(rows).to_json(os.path.join(OUT, f"iqm_resource_real{fb}.json"),  # [deposit] fb suffix
                               orient="records", indent=2)
    print(f"\n저장 -> {OUT}/iqm_device_info{fb}.json, iqm_resource_real{fb}.json")  # [deposit]
    print("\n※ 이 스크립트는 조회만 했습니다. QPU 실행 시간 사용량 0.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
