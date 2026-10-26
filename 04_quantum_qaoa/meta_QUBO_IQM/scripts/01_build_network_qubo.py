#!/usr/bin/env python3
"""
조로증 메타분석 후보 유전자 -> STRING 네트워크 -> QUBO 구성

문제 정의
---------
메타분석에서 살아남은 유전자들이 이루는 상호작용 네트워크에서, 서로 다른
기능 모듈로 나뉘는 최적 분할(Max-Cut)을 찾는다.  Max-Cut 은 QUBO/Ising 으로
정확히 매핑되고, 유전자 1개 = 큐빗 1개다.

  maximize  sum_{(i,j) in E} w_ij * (x_i + x_j - 2 x_i x_j)      x in {0,1}
  <=> minimize  sum_{(i,j) in E} w_ij * (2 x_i x_j - x_i - x_j)

Max-Cut 을 고른 이유
  - STRING 네트워크는 희소(sparse)하다. 밀집 QUBO 와 달리 격자형 초전도 칩
    (IQM Garnet 20q / Emerald 54q) 에 SWAP 폭증 없이 올릴 수 있다.
  - Max-Cut 은 QUBO 벤치마크의 표준 문제라 고전 알고리즘과 공정하게 비교된다.
  - 엄격 필터를 통과한 유전자 수(16~34)를 그대로 쓸 수 있다.
    (유전자 '선택' 문제로 잡으면 후보가 수백 개 필요해 하드웨어에 안 들어간다)

세 가지 문제 크기를 만든다.
  panelA  16 유전자  (사람 3개 일치 + |LFC|>=1 + 마우스 방향일치)  -> Garnet 20q
  panelB  34 유전자  (사람 3개 일치 + |LFC|>=1)                    -> Emerald 54q
  panelC  <=54 유전자 (|LFC|>=0.5 상위)                            -> Emerald 54q 상한
"""
from __future__ import annotations

import json
import os
import sys
import time

import networkx as nx
import numpy as np
import pandas as pd
import requests

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
META = str(DEPOSIT / "03_meta_analysis" / "results")
OUT = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")
os.makedirs(OUT, exist_ok=True)

STRING_API = "https://string-db.org/api"
SPECIES = 9606          # Homo sapiens
MIN_SCORE = 400         # medium confidence (STRING 기본 권장)
CALLER = "progeria_meta_qubo"


# --------------------------------------------------------------------------- #
def load_panels() -> dict[str, pd.DataFrame]:
    M = pd.read_csv(os.path.join(META, "meta_k3_final.tsv"), sep="\t")
    base = (M.meta_padj < 0.05) & (M.all_same_dir) & (M.I2 < 75)

    A = M[base & (M.meta_lfc.abs() >= 1) & (M.mouse_dir_only)].copy()
    B = M[base & (M.meta_lfc.abs() >= 1)].copy()
    C = M[base & (M.meta_lfc.abs() >= 0.5)].copy()
    C = C.reindex(C.meta_lfc.abs().sort_values(ascending=False).index).head(54)

    panels = {"panelA_16": A, "panelB_34": B, "panelC_54": C}
    for k, v in panels.items():
        print(f"  {k:12s} {len(v):3d} 유전자   "
              f"UP {int((v.meta_lfc > 0).sum())} / DOWN {int((v.meta_lfc < 0).sum())}")
    return panels


def string_network(genes: list[str]) -> pd.DataFrame:
    """STRING 에서 유전자 간 상호작용을 받아온다."""
    r = requests.post(
        f"{STRING_API}/tsv/network",
        data={"identifiers": "%0d".join(genes), "species": SPECIES,
              "required_score": MIN_SCORE, "caller_identity": CALLER},
        timeout=120)
    r.raise_for_status()
    from io import StringIO
    df = pd.read_csv(StringIO(r.text), sep="\t")
    if df.empty:
        return df
    # 입력 유전자쌍만 남기고 자기루프/중복 제거
    keep = set(genes)
    df = df[df.preferredName_A.isin(keep) & df.preferredName_B.isin(keep)]
    df = df[df.preferredName_A != df.preferredName_B]
    pair = df.apply(lambda r_: tuple(sorted((r_.preferredName_A, r_.preferredName_B))), axis=1)
    df = df.assign(_pair=pair).drop_duplicates("_pair").drop(columns="_pair")
    return df


def build_graph(genes: list[str], edges: pd.DataFrame, lfc: dict[str, float]) -> nx.Graph:
    G = nx.Graph()
    for g in genes:
        G.add_node(g, lfc=float(lfc[g]))
    if not edges.empty:
        for _, e in edges.iterrows():
            G.add_edge(e.preferredName_A, e.preferredName_B, weight=float(e.score))
    return G


def maxcut_qubo(G: nx.Graph, nodes: list[str]) -> tuple[np.ndarray, float]:
    """Max-Cut 을 최소화형 QUBO 로.  cut = -(x^T Q x) + offset."""
    n = len(nodes)
    idx = {g: i for i, g in enumerate(nodes)}
    Q = np.zeros((n, n))
    offset = 0.0
    for u, v, d in G.edges(data=True):
        w = d["weight"]
        i, j = idx[u], idx[v]
        # cut 기여 = w(x_i + x_j - 2 x_i x_j) -> 최소화형: -w x_i - w x_j + 2w x_i x_j
        Q[i, i] -= w
        Q[j, j] -= w
        Q[i, j] += w
        Q[j, i] += w
        offset += 0.0
    return Q, offset


def qubo_energy(Q: np.ndarray, x: np.ndarray) -> float:
    return float(x @ Q @ x)


def cut_value(G: nx.Graph, nodes: list[str], x: np.ndarray) -> float:
    idx = {g: i for i, g in enumerate(nodes)}
    return float(sum(d["weight"] for u, v, d in G.edges(data=True)
                     if x[idx[u]] != x[idx[v]]))


# --------------------------------------------------------------------------- #
def main() -> int:
    print("=" * 74)
    print(" 조로증 메타분석 후보 -> STRING 네트워크 -> Max-Cut QUBO")
    print("=" * 74)

    print("\n[1] 패널 로드")
    panels = load_panels()

    summary = []
    for name, df in panels.items():
        print(f"\n[2] STRING 조회 : {name} ({len(df)} 유전자)")
        genes = df.gene.tolist()
        lfc = dict(zip(df.gene, df.meta_lfc))
        try:
            # [deposit] STRING_CACHE=1 reuses the archived STRING response saved by the original run
            # (results/<panel>_string_edges.tsv). The live STRING API has drifted since (panelC LCC is
            # 18 nodes/22 edges instead of 19/26 on 2026-10-19), which would no longer match the
            # 19-gene IQM hardware runs. Default (unset) = original behaviour: query STRING online.
            cache_tsv = os.path.join(OUT, f"{name}_string_edges.tsv")
            if os.environ.get("STRING_CACHE") == "1" and os.path.exists(cache_tsv):
                edges = pd.read_csv(cache_tsv, sep="\t")
                print(f"    [deposit] STRING_CACHE=1 -> {os.path.basename(cache_tsv)} 사용 (온라인 조회 생략)")
            else:
                edges = string_network(genes)
        except Exception as exc:                       # noqa: BLE001
            print(f"    STRING 실패: {exc}")
            return 1
        time.sleep(1.0)                                # STRING 예의상 간격

        G = build_graph(genes, edges, lfc)
        iso = [g for g in G.nodes if G.degree(g) == 0]
        comps = sorted(nx.connected_components(G), key=len, reverse=True)
        deg = [d for _, d in G.degree()]
        dens = nx.density(G)
        print(f"    노드 {G.number_of_nodes()}, 엣지 {G.number_of_edges()}, "
              f"밀도 {dens:.3f}, 평균차수 {np.mean(deg):.2f}")
        print(f"    고립 노드 {len(iso)}개, 최대 연결성분 {len(comps[0])}개 노드")
        if iso:
            print(f"    고립: {', '.join(iso)}")

        # 고립 노드는 Max-Cut 에 기여하지 않으므로 큐빗 낭비다 -> 최대 연결성분만 사용
        core = sorted(comps[0])
        Gc = G.subgraph(core).copy()
        print(f"    -> QUBO 대상(최대 연결성분): {len(core)} 노드 / "
              f"{Gc.number_of_edges()} 엣지, 평균차수 {2*Gc.number_of_edges()/len(core):.2f}")

        Q, off = maxcut_qubo(Gc, core)
        np.save(os.path.join(OUT, f"{name}_Q.npy"), Q)
        nx.write_graphml(Gc, os.path.join(OUT, f"{name}_graph.graphml"))
        edges.to_csv(os.path.join(OUT, f"{name}_string_edges.tsv"), sep="\t", index=False)
        with open(os.path.join(OUT, f"{name}_nodes.json"), "w") as fh:
            json.dump({"nodes": core,
                       "lfc": {g: float(lfc[g]) for g in core},
                       "isolated_removed": iso,
                       "n_qubits": len(core)}, fh, indent=2)

        summary.append(dict(
            panel=name, n_genes=len(genes), n_qubits=len(core),
            n_edges_full=G.number_of_edges(), n_edges_core=Gc.number_of_edges(),
            density=round(nx.density(Gc), 4),
            mean_degree=round(2 * Gc.number_of_edges() / len(core), 2),
            max_degree=max(d for _, d in Gc.degree()),
            n_isolated=len(iso),
            brute_force=("가능" if len(core) <= 28 else "불가"),
            hardware=("Garnet(20q)" if len(core) <= 20 else
                      "Emerald(54q)" if len(core) <= 54 else "축소 필요"),
        ))

    S = pd.DataFrame(summary)
    S.to_csv(os.path.join(OUT, "panel_summary.tsv"), sep="\t", index=False)
    print("\n" + "=" * 74)
    print(" 패널 요약")
    print("=" * 74)
    print(S.to_string(index=False))
    print(f"\n저장 -> {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
