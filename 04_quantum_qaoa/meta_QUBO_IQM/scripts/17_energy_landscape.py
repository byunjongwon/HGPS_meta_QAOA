#!/usr/bin/env python3
"""
Fig C 재작성 — 제대로 된 QAOA p=1 에너지 지형

이전 판의 문제
  1) 격자가 18x18 로 너무 거칠어 계단처럼 보였다
  2) 좌측 컬러바 라벨과 우측 y축 라벨이 겹쳤다
  3) 별표를 "하드웨어에서 쓴 각도"라고 적었지만 실제로는 스캔 최댓값이었다

해결
  * 코스트 유니터리 exp(-i g H_C) 는 계산기저에서 대각이므로 위상만 곱하면 된다.
    믹서 exp(-i b sum X) 는 큐빗별 2x2 회전이라 reshape 로 적용한다.
    이 방식이면 statevector 를 매번 새로 만드는 것보다 수십 배 빠르므로
    72 x 54 격자를 현실적인 시간에 계산할 수 있다.
  * 별표는 '스캔으로 찾은 최적 각도' 로 정확히 표기하고,
    하드웨어에 실제로 쓴 qrisp 학습 각도는 별도 표기 없이 캡션에만 적는다
    (파라미터 규약이 달라 같은 좌표계로 겹쳐 그릴 수 없다).
"""
from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
R = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")
OUT = str(DEPOSIT / "04_quantum_qaoa" / "IQM_hardware_results" / "03_figures" / "asfv_style")
NETS = [("string_maxcut_19", "STRING PPI", 19, 60, 46),
        ("coexp_ising_16", "Co-expression", 16, 84, 64)]

plt.rcParams.update({"font.size": 11, "axes.titlesize": 14, "axes.labelsize": 12.5,
                     "xtick.labelsize": 11, "ytick.labelsize": 11})


def spectrum(Q):
    n = Q.shape[0]
    ii = np.arange(1 << n, dtype=np.uint64)
    bb = np.arange(n, dtype=np.uint64)
    X = ((ii[:, None] >> bb) & 1).astype(np.float64)
    return np.einsum("ij,jk,ik->i", X, Q, X)


def apply_mixer(psi, n, beta):
    """exp(-i beta sum_q X_q) 를 큐빗마다 2x2 회전으로 적용."""
    c, s = np.cos(beta), np.sin(beta)
    for q in range(n):
        psi = psi.reshape(-1, 2, 1 << q)
        a = psi[:, 0, :].copy()
        b = psi[:, 1, :].copy()
        psi[:, 0, :] = c * a - 1j * s * b
        psi[:, 1, :] = -1j * s * a + c * b
        psi = psi.reshape(-1)
    return psi


def landscape(Q, E, ng, nb):
    n = Q.shape[0]
    amp0 = np.full(1 << n, 1.0 / np.sqrt(1 << n), dtype=np.complex128)
    gs = np.linspace(0, 2 * np.pi, ng)
    bs = np.linspace(0, np.pi, nb)
    Z = np.empty((nb, ng))
    for i, g in enumerate(gs):
        psi_c = amp0 * np.exp(-1j * g * E)          # 대각 코스트 유니터리
        for j, b in enumerate(bs):
            psi = apply_mixer(psi_c.copy(), n, b)
            Z[j, i] = float((np.abs(psi) ** 2) @ E)
    return gs, bs, Z


def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    # constrained_layout + bbox_inches="tight" 조합은 제목/캡션 자리를 못 지켜줬다.
    # 여백을 직접 고정하고 tight 크롭 없이 저장한다.
    fig, axes = plt.subplots(1, 2, figsize=(17.5, 8.0))
    fig.subplots_adjust(left=.05, right=.955, top=.855, bottom=.195, wspace=.26)

    for ax, (tag, name, n, ng, nb) in zip(axes, NETS):
        Q = np.load(f"{R}/{tag}_Q.npy")
        c = json.load(open(f"{R}/{tag}_classical.json"))
        eopt, erand = c["E_ref"], c["E_random_mean"]
        den = erand - eopt
        # 지형 계산은 수 분 걸리므로 캐시해 두고 재실행 때는 그림만 다시 그린다
        cache = f"{R}/{tag}_landscape_{ng}x{nb}.npz"
        if os.path.exists(cache):
            z = np.load(cache)
            gs, bs, Ze = z["gs"], z["bs"], z["Z"]
            print(f"  {name} (N={n}) 캐시 사용 {os.path.basename(cache)}", flush=True)
        else:
            print(f"  {name} (N={n}) 지형 계산 중 … 격자 {ng} x {nb}", flush=True)
            gs, bs, Ze = landscape(Q, spectrum(Q), ng, nb)
            np.savez_compressed(cache, gs=gs, bs=bs, Z=Ze)
        Z = (erand - Ze) / den                       # 근사비로 환산

        im = ax.pcolormesh(gs, bs, Z, cmap="magma", shading="gouraud")
        cs = ax.contour(gs, bs, Z, levels=8, colors="white", linewidths=.6, alpha=.45)
        ax.clabel(cs, inline=True, fontsize=8, fmt="%.2f")

        k = np.unravel_index(np.argmax(Z), Z.shape)
        gb, bb_ = gs[k[1]], bs[k[0]]
        ax.plot(gb, bb_, marker="*", ms=26, color="#39FF6A", mec="black", mew=1.4,
                zorder=5, linestyle="none",
                label=f"scan optimum   ratio {Z.max():.3f}\n$\\gamma$={gb:.3f}, $\\beta$={bb_:.3f}")
        # 최적 봉우리가 좌상단에 있으므로 범례는 비어 있는 우상단으로 뺀다
        ax.legend(fontsize=10.5, loc="upper right", framealpha=.92, borderpad=.6)
        ax.set_xlabel(r"cost angle  $\gamma$")
        ax.set_ylabel(r"mixer angle  $\beta$")
        ax.set_title(f"{name}   ({n} genes,  {1 << n:,} states)", fontweight="bold", pad=12)
        ax.set_xticks(np.linspace(0, 2 * np.pi, 5))
        ax.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
        ax.set_yticks(np.linspace(0, np.pi, 5))
        ax.set_yticklabels(["0", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])
        cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.02)
        cb.set_label("approximation ratio", fontsize=11.5)
        print(f"    최적 ratio {Z.max():.4f} at gamma={gb:.4f}, beta={bb_:.4f}", flush=True)

    fig.text(.5, .945, "C   QAOA energy landscape at p = 1  (exact statevector, noiseless)",
             ha="center", va="center", fontsize=17, fontweight="bold")
    fig.text(.5, 0.062,
             "Every point is the exact expectation value of the full 2$^N$-state wavefunction. "
             "Bright = better.  White contours mark equal approximation ratio.\n"
             "Hardware runs used angles trained independently with Qrisp "
             "(γ = 1.717, β = 1.095 in Qrisp's convention), which is not the same "
             "parameterisation as this scan and so is not overlaid here.",
             ha="center", fontsize=10.5, color="#5A6672", style="italic", linespacing=1.5)
    p = f"{OUT}/C_energy_landscape.png"
    fig.savefig(p, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"\n저장 -> {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
