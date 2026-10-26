#!/usr/bin/env python3
"""
Qrisp QAOA — 조로증 공통 유전자 QUBO 를 IQM 실기에서 푼다

실행 순서 (반드시 이 순서로)
    1) --mode local   qrisp 로컬 시뮬레이터에서 각도 최적화 + 검증   (QPU 0)
    2) --mode mock    garnet:mock 에서 동일 회로 실행                (QPU 0, 무료)
    3) --mode hw      실제 garnet / emerald 제출                      (QPU 사용)

설계 원칙
    하드웨어에서 변분 최적화를 하지 않는다.
    각도는 1)에서 확정하고, 하드웨어에는 '고정 각도 회로 1개'만 보낸다.
    (QAOAProblem.run 을 하드웨어에 직접 걸면 COBYLA 가 회로를 수십~수백 번
     제출해 1분 예산을 즉시 초과한다)

사용법
    source ~/.iqm_emerald.sh
    python 11_qrisp_qaoa_iqm.py --mode local --instance coexp_ising_16 --depth 1
    python 11_qrisp_qaoa_iqm.py --mode mock  --instance coexp_ising_16 --depth 1
    python 11_qrisp_qaoa_iqm.py --mode hw    --instance coexp_ising_16 --depth 1 --yes
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import numpy as np

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
R = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")
BUDGET_S = 60.0


# --------------------------------------------------------------------------- #
def load_instance(tag: str):
    Q = np.load(os.path.join(R, f"{tag}_Q.npy"))
    meta = json.load(open(os.path.join(R, f"{tag}_meta.json")))
    cls = json.load(open(os.path.join(R, f"{tag}_classical.json")))
    return Q, meta, cls


def energy(Q: np.ndarray, x: np.ndarray) -> float:
    return float(x @ Q @ x)


def key_to_bits(key, n: int) -> np.ndarray:
    """qrisp QuantumArray 결과 키(OutcomeArray) -> 0/1 벡터.

    OutcomeArray(['0','1',...]) 형태이며 배열 순서가 곧 큐빗 순서이므로
    비트 뒤집기(endian) 문제가 없다.
    """
    try:
        flat = np.array(key).flatten()
        return np.array([int(str(v)) for v in flat], dtype=float)[:n]
    except Exception:                                               # noqa: BLE001
        b = str(key).replace(" ", "")[-n:]
        return np.array([int(c) for c in b], dtype=float)


def analyse_counts(counts: dict, Q: np.ndarray, e_ref: float, e_rand: float) -> dict:
    """qrisp 결과(dict: OutcomeArray -> 확률/카운트) 를 에너지로 환산."""
    n = Q.shape[0]
    tot = sum(counts.values())
    esum, ebest, bbest = 0.0, np.inf, None
    for key, w in counts.items():
        arr = key_to_bits(key, n)
        if arr.size != n:
            continue
        e = energy(Q, arr)
        esum += e * w
        if e < ebest:
            ebest, bbest = e, arr.astype(int).tolist()
    denom = e_rand - e_ref
    return dict(n_states=len(counts), shots=tot,
                E_mean=esum / tot, E_best=ebest, best_bits=bbest,
                ratio_mean=(e_rand - esum / tot) / denom,
                ratio_best=(e_rand - ebest) / denom,
                reached_optimum=bool(abs(ebest - e_ref) < 1e-9))


def best_orientation(counts, Q, e_ref, e_rand):
    return analyse_counts(counts, Q, e_ref, e_rand)


# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["local", "mock", "hw"], required=True)
    ap.add_argument("--instance", default="coexp_ising_16")
    ap.add_argument("--depth", type=int, default=1)
    ap.add_argument("--shots", type=int, default=4096)
    ap.add_argument("--max-iter", type=int, default=60, help="로컬 최적화 반복 수")
    ap.add_argument("--device", default="garnet", help="hw 모드에서 garnet 또는 emerald")
    ap.add_argument("--yes", action="store_true", help="hw 모드 실제 제출 확인")
    ap.add_argument("--repeat", type=int, default=1,
                    help="같은 회로를 N회 반복 측정 (오차막대 확보). 학습은 1회만 수행")
    ap.add_argument("--theta", type=str, default="",
                    help="로컬 학습을 건너뛰고 각도를 직접 지정 (예: 1.717,1.095). "
                         "N>22 라 statevector 시뮬레이션이 불가능한 인스턴스에 사용")
    args = ap.parse_args()

    Q, meta, cls = load_instance(args.instance)
    # [deposit] credential-less fallback: mock/hw mode without IQM_TOKEN runs on a local simulator and
    # every file it writes gets the suffix _SIMULATOR_FALLBACK (real-hardware files are never overwritten).
    fallback = args.mode in ("mock", "hw") and not os.environ.get("IQM_TOKEN", "").strip()
    fbs = "_SIMULATOR_FALLBACK" if fallback else ""
    n = Q.shape[0]
    e_ref, e_rand = cls["E_ref"], cls["E_random_mean"]
    genes = meta["nodes"]

    print("=" * 88)
    print(f" Qrisp QAOA   인스턴스 {args.instance}   N={n} 큐빗   depth={args.depth}   mode={args.mode}")
    print("=" * 88)
    print(f"  유전자: {', '.join(genes)}")
    print(f"  참조해 E = {e_ref:.4f} ({cls['ref_source']}),  무작위 평균 E = {e_rand:.4f}")

    from qrisp import QuantumArray, QuantumVariable
    from qrisp.qaoa import QUBO_problem, QUBO_obj, create_QUBO_cl_cost_function

    problem = QUBO_problem(Q)
    qarg = QuantumArray(qtype=QuantumVariable(1), shape=(n,))

    # ------------------------------------------------------- 1') 각도 직접 지정
    if args.theta:
        from qrisp import h
        th = [float(x) for x in args.theta.split(",")]
        if len(th) != 2 * args.depth:
            sys.exit(f"--theta 는 {2*args.depth}개 값이 필요합니다 (depth={args.depth})")
        print(f"\n[1] 로컬 학습 건너뜀 — 각도 직접 지정 {['%.4f' % t for t in th]}")
        print(f"    (N={n} 은 statevector 시뮬레이션이 불가능하므로 다른 인스턴스에서 전이)")

        def trained(qarg_gen):
            """train_function 의 circuit_generator 와 동일한 구조."""
            h(qarg_gen)
            for i in range(args.depth):
                problem.cost_operator(qarg_gen, th[i])
                problem.mixer(qarg_gen, th[i + args.depth])

        opt_theta = th
        a_loc = None
        if args.mode == "local":
            print("    local 모드에서는 할 일이 없습니다 (시뮬레이션 불가).")
            return 0
    else:
        _train_locally = True

    if not args.theta:
        pass
    # ------------------------------------------------------- 1) 로컬 학습
    # qrisp 의 train_function 은 '학습된 각도가 박힌 회로 생성기'를 돌려준다.
    # 학습은 로컬 시뮬레이터에서만 하고(무료), 하드웨어에서는 그 회로를 한 번 측정만 한다.
    if not args.theta:
      print(f"\n[1] 로컬 시뮬레이터에서 각도 학습 (max_iter={args.max_iter}) …")
      t0 = time.time()
      trained = problem.train_function(qarg, depth=args.depth,
                                       mes_kwargs={"shots": args.shots},
                                       max_iter=args.max_iter)
      dt_train = time.time() - t0

      opt_theta = None
      for cell in (trained.__closure__ or []):
          try:
              v = np.atleast_1d(cell.cell_contents)
              if v.size == 2 * args.depth and np.issubdtype(np.asarray(v).dtype, np.number):
                  opt_theta = [float(x) for x in v]
          except Exception:                                            # noqa: BLE001
              continue
      print(f"    학습 완료 {dt_train:.1f}초"
            + (f",  최적 각도 {['%.4f' % t for t in opt_theta]}" if opt_theta else ""))

      # 학습된 회로를 로컬에서 한 번 평가해 기준값 확보
      qloc = QuantumArray(qtype=QuantumVariable(1), shape=(n,))
      trained(qloc)
      res_local = qloc.get_measurement(shots=args.shots)
      a_loc = analyse_counts(res_local, Q, e_ref, e_rand)
      print(f"    [로컬] 기대 E {a_loc['E_mean']:.4f}  ratio {a_loc['ratio_mean']:.4f}   "
            f"최선 ratio {a_loc['ratio_best']:.4f}  "
            f"{'최적해 도달' if a_loc['reached_optimum'] else '최적해 미도달'}")

      with open(os.path.join(R, f"qrisp_angles_{args.instance}_p{args.depth}{fbs}.json"), "w") as fh:  # [deposit] fbs
          json.dump({"instance": args.instance, "depth": args.depth,
                     "theta": opt_theta, "train_seconds": round(dt_train, 1),
                     "local": {k: v for k, v in a_loc.items() if k != "best_bits"},
                     "best_bits": a_loc["best_bits"]}, fh, indent=2)

      arr = a_loc["best_bits"]
      g0 = [g for g, v in zip(genes, arr) if v == 0]
      g1 = [g for g, v in zip(genes, arr) if v == 1]
      print(f"    모듈 1 ({len(g0)}): {', '.join(g0)}")
      print(f"    모듈 2 ({len(g1)}): {', '.join(g1)}")
      if args.mode == "local":
          return 0

    # ------------------------------------------------------- 2) 백엔드 준비
    tok = os.environ.get("IQM_TOKEN", "").strip()
    device = "garnet:mock" if args.mode == "mock" else args.device
    if fallback:
        # [deposit] no token: IQM's fake Garnet-class IQMFakeApollo / Emerald-class IQMFakeAphrodite
        # (real coupling map + IQM noise model, executed locally by Aer) wrapped as a qrisp
        # VirtualBackend that calls fake.run() directly (qrisp's QiskitBackend/SamplerV2 path returns
        # all-zero counts for IQM fake backends). N>24 cannot be simulated with noise (or exactly: an
        # exact MPS needs ~4.5 TB for coexp_ising_34), so an APPROXIMATE noiseless matrix-product-state
        # AerSimulator with bond dimension capped at 64 is used instead.
        from qiskit import QuantumCircuit as _QC, transpile as _tp
        from qrisp.interface import VirtualBackend
        if n <= 24:
            from iqm.qiskit_iqm import IQMFakeAphrodite, IQMFakeApollo
            fake = IQMFakeAphrodite() if "emerald" in device else IQMFakeApollo()
        else:
            from qiskit_aer import AerSimulator
            fake = AerSimulator(method="matrix_product_state",
                                matrix_product_state_max_bond_dimension=64,
                                max_memory_mb=10**9)  # bypass the exact-MPS memory pre-check

        def _fake_run(qasm_str, shots=None, token=""):
            qc_ = _QC.from_qasm_str(qasm_str)
            flat = _QC(len(qc_.qubits), len(qc_.clbits))
            for ins in qc_:
                flat.append(ins.operation, [qc_.qubits.index(q) for q in ins.qubits],
                            [qc_.clbits.index(c) for c in ins.clbits])
            cnt = fake.run(_tp(flat, backend=fake), shots=shots or 1000).result().get_counts()
            return {k.replace(" ", ""): v for k, v in cnt.items()}

        print(f"\n[deposit] no IQM_TOKEN → running on local simulator ({fake.name}) instead of "
              f"device_instance={device}; outputs *_SIMULATOR_FALLBACK")
        backend = VirtualBackend(_fake_run)
    else:  # [deposit] original real-device path (unchanged, re-indented)
     if not tok:
        sys.exit("IQM_TOKEN 이 없습니다.  source ~/.iqm_emerald.sh")
     print(f"\n[2] Qrisp IQMBackend 연결  device_instance={device}")
     from qrisp.interface import IQMBackend
     # qrisp 는 api_token 인자를 IQM 클라이언트에 넘기는데, 동시에 IQM_TOKEN 환경변수가
     # 설정되어 있으면 "인증 파라미터 출처를 섞지 말라"는 오류가 난다.
     # 생성하는 동안만 환경변수를 잠시 치운다.
     saved = os.environ.pop("IQM_TOKEN", None)
     saved_url = os.environ.pop("IQM_SERVER_URL", None)
     try:
        backend = IQMBackend(api_token=tok, device_instance=device)
     except Exception as exc:                                        # noqa: BLE001
        sys.exit(f"    연결 실패: {type(exc).__name__}: {str(exc)[:250]}")
     finally:
        if saved is not None:
            os.environ["IQM_TOKEN"] = saved
        if saved_url is not None:
            os.environ["IQM_SERVER_URL"] = saved_url
     print("    연결 성공")

    if args.mode == "hw":
        est = args.shots * 300e-6
        print(f"\n    예상 QPU 시간 약 {est:.1f}초 (shots {args.shots}) / 예산 {BUDGET_S:.0f}초")
        if not args.yes:
            print("    --yes 가 없어 제출하지 않고 종료합니다.")
            return 0

    # ------------------------------------------------------- 3) 제출
    print(f"\n[3] {'실기' if args.mode=='hw' else 'mock'} 실행  shots={args.shots}  x{args.repeat}회")
    print("    학습된 고정 각도 회로를 반복 측정 (하드웨어 변분 최적화 없음)")
    runs, t_all = [], time.time()
    for rep in range(1, args.repeat + 1):
        qhw = QuantumArray(qtype=QuantumVariable(1), shape=(n,))
        trained(qhw)
        t0 = time.time()
        try:
            res = qhw.get_measurement(backend=backend, shots=args.shots)
        except Exception as exc:                                    # noqa: BLE001
            print(f"    [{rep}/{args.repeat}] 실패: {type(exc).__name__}: {str(exc)[:150]}")
            continue
        dt_r = time.time() - t0
        a_r = analyse_counts(res, Q, e_ref, e_rand)
        a_r["rep"], a_r["wall_s"] = rep, round(dt_r, 2)
        a_r["counts"] = {"".join(map(str, key_to_bits(k, n).astype(int))): float(v)
                         for k, v in res.items()}
        runs.append(a_r)
        print(f"    [{rep}/{args.repeat}] 평균 ratio {a_r['ratio_mean']:+.4f}   "
              f"최선 ratio {a_r['ratio_best']:.4f}   "
              f"{'최적해 ✅' if a_r['reached_optimum'] else '미도달  '}   "
              f"해 {a_r['n_states']}개   {dt_r:.1f}초")
    if not runs:
        sys.exit("    모든 반복이 실패했습니다.")
    dt = time.time() - t_all

    mr = np.array([r["ratio_mean"] for r in runs])
    br = np.array([r["ratio_best"] for r in runs])
    hit = sum(r["reached_optimum"] for r in runs)
    sd = lambda v: float(v.std(ddof=1)) if len(v) > 1 else 0.0
    print(f"\n    === {len(runs)}회 요약 ===")
    print(f"    평균 ratio : {mr.mean():+.4f} +- {sd(mr):.4f}   (범위 {mr.min():+.4f} ~ {mr.max():+.4f})")
    print(f"    최선 ratio : {br.mean():.4f} +- {sd(br):.4f}   (범위 {br.min():.4f} ~ {br.max():.4f})")
    print(f"    최적해 도달: {hit}/{len(runs)}회")
    print(f"    총 벽시계 {dt:.1f}초,  QPU 추정 {len(runs)*args.shots*300e-6:.1f}초")

    a = runs[int(np.argmin([r['E_best'] for r in runs]))]
    arr = a["best_bits"]
    g0 = [g for g, v in zip(genes, arr) if v == 0]
    g1 = [g for g, v in zip(genes, arr) if v == 1]
    print(f"\n    최고 성적 회차의 분할")
    print(f"    모듈 1 ({len(g0)}): {', '.join(g0)}")
    print(f"    모듈 2 ({len(g1)}): {', '.join(g1)}")
    if a_loc:
        print(f"    로컬 대비:  로컬 ratio {a_loc['ratio_mean']:.4f} -> {args.mode} {mr.mean():.4f}")

    pred = json.load(open(os.path.join(R, "iqm_predicted.json")))
    pm = next((r for r in pred if r["instance"] == args.instance and r["p"] == args.depth), None)
    if pm:
        print(f"    예측 대조: 평균 ratio 예측 {pm['ratio_noisy_mean']:.4f} / 실측 {mr.mean():.4f}")
        print(f"               최선 ratio 예측 {pm['ratio_noisy_best']:.4f} / 실측 {br.mean():.4f}")

    out = dict(engine="qrisp", mode=args.mode, device=device,
               simulator_fallback=fallback,  # [deposit]
               instance=args.instance, n_qubits=n, depth=args.depth,
               n_shots=args.shots, n_repeat=len(runs), wall_s=round(dt, 2),
               genes=genes, theta=opt_theta, module_1=g0, module_2=g1,
               local_ratio_mean=(a_loc["ratio_mean"] if a_loc else None),
               local_ratio_best=(a_loc["ratio_best"] if a_loc else None),
               summary=dict(ratio_mean_avg=float(mr.mean()), ratio_mean_sd=sd(mr),
                            ratio_best_avg=float(br.mean()), ratio_best_sd=sd(br),
                            optimum_hits=int(hit), n_runs=len(runs)),
               runs=[{k: v for k, v in r.items() if k != "counts"} for r in runs],
               counts_per_run=[r["counts"] for r in runs],
               **{k: v for k, v in a.items() if k not in ("counts", "wall_s", "rep")})
    suffix = f"_x{len(runs)}" if len(runs) > 1 else ""
    fn = os.path.join(R, f"qrisp_{args.mode}_{args.instance}_p{args.depth}_"
                         f"{device.replace(':','_')}{suffix}{fbs}.json")  # [deposit] fbs
    with open(fn, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\n저장 -> {fn}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
