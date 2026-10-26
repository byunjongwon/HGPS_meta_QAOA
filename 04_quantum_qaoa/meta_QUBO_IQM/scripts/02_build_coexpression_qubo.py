#!/usr/bin/env python3
"""
01번의 결과: 엄격 필터를 통과한 유전자들은 STRING PPI 에서 거의 연결되지 않는다
(panelA 16개 중 14개 고립).  이는 버그가 아니라 실질적 발견이다 —
"여러 데이터셋·종에서 재현되는 조로증 시그니처"는 하나의 단백질 복합체나
경로 모듈이 아니라 기능적으로 분산된 유전자 집합이다.

따라서 QUBO 인스턴스를 두 계열로 만든다.

  [계열 1] STRING PPI  : panelC 의 최대 연결성분 (19노드/26엣지) — 실제 상호작용 기반
  [계열 2] 공발현 Ising: 유전자 간 발현 상관을 결합상수 J_ij 로 쓰는 Ising 바닥상태 탐색

계열 2 의 물리적 정의
---------------------
    H(s) = - sum_{i<j} J_ij s_i s_j          s_i in {-1, +1},  J_ij = corr(g_i, g_j)

  양의 상관 유전자는 같은 부호(같은 모듈), 음의 상관은 반대 부호가 되려 한다.
  경쟁하는 제약이 섞이면 frustration 이 생기고 바닥상태 탐색은 NP-hard 다.
  = 부호 있는 네트워크의 최적 2-모듈 분할 문제이며, Ising/QUBO 로 정확히 매핑된다.

  희소화: 각 노드마다 |J| 상위 k개 엣지만 남긴다 (k-NN 그래프).
  -> 밀집 QUBO 의 SWAP 폭증 없이 격자형 칩에 올릴 수 있다.

QUBO 변환:  s = 1 - 2x  (x in {0,1})
    H = -sum J_ij (1-2x_i)(1-2x_j)
      = -sum J_ij + 2 sum J_ij (x_i + x_j) - 4 sum J_ij x_i x_j
"""
from __future__ import annotations

import json
import os
import sys

import networkx as nx
import numpy as np
import pandas as pd

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
META = str(DEPOSIT / "03_meta_analysis" / "results")
OUT = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")
os.makedirs(OUT, exist_ok=True)

VST_306264 = None  # 아래에서 R 산출물로부터 읽는다
KNN_K = 4          # 노드당 유지할 엣지 수 (희소성 제어)


def load_expression() -> pd.DataFrame:
    """GSE306264 VST (30 세포주) — 가장 표본이 큰 사람 데이터셋을 공발현 기준으로 쓴다."""
    p = str(DEPOSIT / "02_per_dataset_DE" / "GSE306264" / "results" / "vst_for_quantum.tsv")  # [deposit]
    if not os.path.exists(p):
        sys.exit(f"VST 파일이 없습니다: {p}\n  -> 03_export_vst.R 을 먼저 실행하세요")
    return pd.read_csv(p, sep="\t", index_col=0)


def ising_from_corr(C: np.ndarray, nodes: list[str], k: int) -> nx.Graph:
    """상관행렬 -> k-NN 희소 부호 그래프."""
    n = len(nodes)
    G = nx.Graph()
    G.add_nodes_from(nodes)
    keep: set[tuple[int, int]] = set()
    for i in range(n):
        order = np.argsort(-np.abs(C[i]))
        cnt = 0
        for j in order:
            if j == i:
                continue
            keep.add((min(i, j), max(i, j)))
            cnt += 1
            if cnt >= k:
                break
    for i, j in sorted(keep):
        G.add_edge(nodes[i], nodes[j], weight=float(C[i, j]))
    return G


def ising_to_qubo(G: nx.Graph, nodes: list[str]) -> tuple[np.ndarray, float]:
    """H = -sum J_ij s_i s_j,  s = 1-2x  ->  x^T Q x + offset."""
    idx = {g: i for i, g in enumerate(nodes)}
    n = len(nodes)
    Q = np.zeros((n, n))
    offset = 0.0
    for u, v, d in G.edges(data=True):
        J = d["weight"]
        i, j = idx[u], idx[v]
        offset += -J
        Q[i, i] += 2 * J
        Q[j, j] += 2 * J
        Q[i, j] += -2 * J
        Q[j, i] += -2 * J
    return Q, offset


def maxcut_to_qubo(G: nx.Graph, nodes: list[str]) -> tuple[np.ndarray, float]:
    idx = {g: i for i, g in enumerate(nodes)}
    Q = np.zeros((len(nodes), len(nodes)))
    for u, v, d in G.edges(data=True):
        w = d["weight"]
        i, j = idx[u], idx[v]
        Q[i, i] -= w
        Q[j, j] -= w
        Q[i, j] += w
        Q[j, i] += w
    return Q, 0.0


def describe(G: nx.Graph, tag: str) -> dict:
    deg = [d for _, d in G.degree()]
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    info = dict(instance=tag, n_qubits=G.number_of_nodes(), n_edges=G.number_of_edges(),
                density=round(nx.density(G), 4), mean_degree=round(float(np.mean(deg)), 2),
                max_degree=int(max(deg)), n_components=len(comps),
                largest_component=len(comps[0]))
    print(f"    {tag:26s} 노드 {info['n_qubits']:3d}  엣지 {info['n_edges']:4d}  "
          f"밀도 {info['density']:.3f}  평균차수 {info['mean_degree']:.2f}  "
          f"성분 {info['n_components']}")
    return info


def save_instance(tag: str, G: nx.Graph, nodes: list[str], Q: np.ndarray,
                  offset: float, kind: str, lfc: dict[str, float]) -> None:
    np.save(os.path.join(OUT, f"{tag}_Q.npy"), Q)
    nx.write_graphml(G, os.path.join(OUT, f"{tag}_graph.graphml"))
    with open(os.path.join(OUT, f"{tag}_meta.json"), "w") as fh:
        json.dump({"instance": tag, "kind": kind, "nodes": nodes,
                   "n_qubits": len(nodes), "offset": offset,
                   "lfc": {g: float(lfc.get(g, np.nan)) for g in nodes},
                   "edges": [[u, v, float(d["weight"])] for u, v, d in G.edges(data=True)]},
                  fh, indent=2)


def main() -> int:
    print("=" * 78)
    print(" QUBO 인스턴스 구축 : STRING PPI + 공발현 Ising")
    print("=" * 78)

    M = pd.read_csv(os.path.join(META, "meta_k3_final.tsv"), sep="\t")
    base = (M.meta_padj < 0.05) & (M.all_same_dir) & (M.I2 < 75)
    panels = {
        "panelA_16": M[base & (M.meta_lfc.abs() >= 1) & M.mouse_dir_only],
        "panelB_34": M[base & (M.meta_lfc.abs() >= 1)],
        "panelC_97": M[base & (M.meta_lfc.abs() >= 0.5)],
    }
    lfc_all = dict(zip(M.gene, M.meta_lfc))

    expr = load_expression()
    print(f"\n발현행렬: {expr.shape[0]} 유전자 x {expr.shape[1]} 세포주 (GSE306264 VST)")

    rows = []

    # ---------------------------------------------------------------- 계열 1
    print("\n[계열 1] STRING PPI 기반 Max-Cut")
    gml = os.path.join(OUT, "panelC_54_graph.graphml")
    if os.path.exists(gml):
        Gs = nx.read_graphml(gml)
        nodes = sorted(Gs.nodes())
        Q, off = maxcut_to_qubo(Gs, nodes)
        rows.append({**describe(Gs, "string_maxcut_19"), "kind": "STRING Max-Cut",
                     "source": "panelC_54 largest component"})
        save_instance("string_maxcut_19", Gs, nodes, Q, off, "maxcut_string", lfc_all)
    else:
        print("    panelC graphml 없음 — 01번을 먼저 실행하세요")

    # ---------------------------------------------------------------- 계열 2
    print(f"\n[계열 2] 공발현 Ising (노드당 |J| 상위 {KNN_K}개 엣지)")
    for name, df in panels.items():
        genes = [g for g in df.gene.tolist() if g in expr.index]
        miss = len(df) - len(genes)
        X = expr.loc[genes].to_numpy()
        C = np.corrcoef(X)
        np.fill_diagonal(C, 0.0)
        G = ising_from_corr(C, genes, KNN_K)
        tag = f"coexp_ising_{len(genes)}"
        info = describe(G, tag)
        if miss:
            print(f"      (발현행렬에 없어 제외된 유전자 {miss}개)")
        Q, off = ising_to_qubo(G, genes)
        save_instance(tag, G, genes, Q, off, "ising_coexpression", lfc_all)
        rows.append({**info, "kind": "co-expression Ising", "source": name})

        # 상관 분포 요약
        w = np.array([d["weight"] for _, _, d in G.edges(data=True)])
        print(f"      결합상수 J: 양수 {int((w > 0).sum())} / 음수 {int((w < 0).sum())}, "
              f"|J| 범위 {np.abs(w).min():.2f}~{np.abs(w).max():.2f}  "
              f"-> frustration 존재 {'예' if (w > 0).any() and (w < 0).any() else '아니오'}")

    S = pd.DataFrame(rows)
    S["brute_force"] = np.where(S.n_qubits <= 28, "가능 (2^N)", "불가")
    S["hardware"] = np.select(
        [S.n_qubits <= 20, S.n_qubits <= 54],
        ["Garnet(20q) / Emerald(54q)", "Emerald(54q)"], "축소 필요")
    S.to_csv(os.path.join(OUT, "qubo_instances.tsv"), sep="\t", index=False)
    print("\n" + "=" * 78)
    print(" QUBO 인스턴스 요약")
    print("=" * 78)
    print(S[["instance", "kind", "n_qubits", "n_edges", "mean_degree", "max_degree",
             "n_components", "brute_force", "hardware"]].to_string(index=False))
    print(f"\n저장 -> {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
