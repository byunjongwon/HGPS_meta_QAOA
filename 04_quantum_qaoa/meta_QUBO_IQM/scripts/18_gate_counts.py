#!/usr/bin/env python3
"""
모든 인스턴스 x 깊이 조합의 트랜스파일 후 2큐빗 게이트 수를 한 곳에서 산출한다.

왜 필요한가
  12 / 13 / 16 번 스크립트가 서로 다른 게이트 수를 쓰고 있었다
  (19유전자 p=2 를 156 이라 한 곳과 312 라 한 곳, 16유전자 p=1 을 156 / 164 로
  적은 곳이 섞여 있었다). "실용 한계 약 30개 2큐빗 게이트" 라는 결론이
  이 숫자에 직접 걸려 있으므로 하나의 산출물로 고정한다.

방법
  07번이 계획을 세울 때 쓴 것과 동일하게 IQMFakeApollo(Garnet급) 위에서
  optimization_level=3, seed_transpiler=11 로 트랜스파일한 결과를 센다.
  결과는 results/gate_counts.json 에 저장하고 다른 스크립트가 이것만 읽는다.
"""
from __future__ import annotations

import json
import os  # [deposit] needed for PROGERIA_DEPOSIT

import numpy as np
from iqm.qiskit_iqm import IQMFakeApollo
from qiskit import QuantumCircuit, transpile

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
R = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")
INST = ["coexp_ising_8", "coexp_ising_16", "string_maxcut_19"]
DEPTHS = [1, 2]


def qaoa_circuit(Q, p):
    """QAOA ansatz — 게이트 수만 셀 것이므로 각도는 임의값으로 둔다."""
    n = Q.shape[0]
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for layer in range(p):
        g = 0.3 * (layer + 1)
        for i in range(n):
            for j in range(i + 1, n):
                if abs(Q[i, j]) > 1e-12:
                    qc.cx(i, j)
                    qc.rz(2 * g * Q[i, j], j)
                    qc.cx(i, j)
            if abs(Q[i, i]) > 1e-12:
                qc.rz(g * Q[i, i], i)
        qc.rx(2 * 0.4 * (layer + 1), range(n))
    qc.measure_all()
    return qc


def main() -> int:
    dev = IQMFakeApollo()
    out = {}
    print(f"{'instance':20s} {'p':>2s} {'edges':>6s} {'logical 2q':>11s} "
          f"{'transpiled 2q':>14s} {'depth':>7s}")
    for tag in INST:
        Q = np.load(f"{R}/{tag}_Q.npy")
        n_edge = int((np.abs(np.triu(Q, 1)) > 1e-12).sum())
        for p in DEPTHS:
            qc = qaoa_circuit(Q, p)
            tq = transpile(qc, backend=dev, optimization_level=3, seed_transpiler=11)
            ops = tq.count_ops()
            n2q = sum(v for k, v in ops.items()
                      if k in ("cz", "cx", "ecr", "iswap", "cp"))
            out[f"{tag}_p{p}"] = dict(instance=tag, n_qubits=int(Q.shape[0]), depth=p,
                                      edges=n_edge, logical_2q=2 * n_edge * p,
                                      transpiled_2q=int(n2q),
                                      circuit_depth=int(tq.depth()))
            print(f"{tag:20s} {p:2d} {n_edge:6d} {2 * n_edge * p:11d} "
                  f"{n2q:14d} {tq.depth():7d}")
    json.dump(out, open(f"{R}/gate_counts.json", "w"), indent=2)
    print(f"\n저장 -> {R}/gate_counts.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
