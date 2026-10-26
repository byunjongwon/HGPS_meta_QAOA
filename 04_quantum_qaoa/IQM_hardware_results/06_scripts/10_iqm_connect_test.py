#!/usr/bin/env python3
"""
IQM 연결 진단 — 토큰을 넣은 뒤 가장 먼저 실행할 것

QPU 를 사용하지 않는다. 어떤 URL / 어떤 디바이스에 접근 가능한지만 확인한다.
토큰 값은 절대 출력하지 않는다.

사용법
    source ~/.iqm_emerald.sh
    python 10_iqm_connect_test.py
"""
from __future__ import annotations

import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")

CANDIDATE_URLS = [
    "https://resonance.meetiqm.com",
    "https://cocos.resonance.meetiqm.com/garnet",
    "https://cocos.resonance.meetiqm.com/emerald",
    "https://cocos.resonance.meetiqm.com/deneb",
]
# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
from pathlib import Path as _Path
DEPOSIT = _Path(os.environ.get("PROGERIA_DEPOSIT", _Path(__file__).resolve().parents[3]))
OUT = str(DEPOSIT / "04_quantum_qaoa" / "meta_QUBO_IQM" / "results")


def mask(tok: str) -> str:
    return f"<{len(tok)}자 토큰>" if tok else "<없음>"


def main() -> int:
    tok = os.environ.get("IQM_TOKEN", "").strip()
    url_env = os.environ.get("IQM_SERVER_URL", "").strip()
    print("=" * 78)
    print(" IQM 연결 진단  (QPU 실행 없음)")
    print("=" * 78)
    print(f"  IQM_TOKEN      : {mask(tok)}")
    print(f"  IQM_SERVER_URL : {url_env or '<없음>'}")
    if not tok:
        print("\n토큰이 비어 있습니다.  ~/.iqm_emerald.sh 의 IQM_TOKEN= 뒤에 값을 넣고")
        print("  source ~/.iqm_emerald.sh   후 다시 실행하세요.")
        # [deposit] credential-less fallback: report the properties of IQM's fake Garnet-class backend
        # (IQMFakeApollo, 20q) instead of a live device, write iqm_live_backend_SIMULATOR_FALLBACK.json
        # (never the real iqm_live_backend.json) and exit 0.
        print("\n[deposit] no IQM_TOKEN → running on local simulator (IQMFakeApollo, Garnet-class fake backend)")
        from iqm.qiskit_iqm import IQMFakeApollo
        bk = IQMFakeApollo()
        edges = sorted({tuple(sorted(e)) for e in bk.coupling_map.get_edges()})
        info = dict(url="local-fake-backend", backend=bk.name, num_qubits=int(bk.num_qubits),
                    n_edges=len(edges), basis_gates=list(bk.operation_names),
                    edges=[list(e) for e in edges], note="SIMULATOR_FALLBACK (no IQM_TOKEN)")
        os.makedirs(OUT, exist_ok=True)
        with open(os.path.join(OUT, "iqm_live_backend_SIMULATOR_FALLBACK.json"), "w") as fh:
            json.dump(info, fh, indent=2)
        print(f"  백엔드   : {bk.name}\n  큐빗     : {bk.num_qubits}\n  커플링   : 엣지 {len(edges)}개"
              f"\n  기본게이트: {bk.operation_names}")
        print(f"\n저장 -> {OUT}/iqm_live_backend_SIMULATOR_FALLBACK.json")
        return 0

    # ---------------------------------------------------- 1) 원시 HTTP 확인
    import requests
    print("\n[1] REST 엔드포인트 접근 확인")
    reachable = []
    for base in ["https://api.resonance.meetiqm.com", "https://resonance.meetiqm.com"]:
        for path in ["/quantum-computers", "/api/v1/quantum-computers"]:
            u = base + path
            try:
                r = requests.get(u, headers={"Authorization": f"Bearer {tok}"},
                                 timeout=25, allow_redirects=True)
                mark = "OK " if r.status_code == 200 else "   "
                print(f"  {mark} {r.status_code}  {u}")
                if r.status_code == 200:
                    reachable.append((u, r))
            except Exception as exc:                      # noqa: BLE001
                print(f"       ERR  {u}  ({type(exc).__name__})")

    for u, r in reachable:
        try:
            data = r.json()
            print(f"\n  {u} 응답:")
            print("   ", json.dumps(data, ensure_ascii=False)[:800])
        except Exception:
            print(f"\n  {u} 응답(텍스트): {r.text[:400]}")

    # ---------------------------------------------------- 2) Qiskit 어댑터
    print("\n[2] Qiskit IQM 어댑터로 백엔드 연결 시도")
    from iqm.qiskit_iqm import IQMProvider

    urls = ([url_env] if url_env else []) + [u for u in CANDIDATE_URLS if u != url_env]
    ok_backend = None
    for u in urls:
        try:
            prov = IQMProvider(u)
            bk = prov.get_backend()
            print(f"  OK   {u}")
            print(f"       백엔드 {bk.name},  큐빗 {bk.num_qubits}")
            ok_backend = (u, bk)
            break
        except Exception as exc:                          # noqa: BLE001
            msg = str(exc).splitlines()[0][:110]
            print(f"       실패  {u}\n             {type(exc).__name__}: {msg}")

    if ok_backend is None:
        print("\n" + "=" * 78)
        print(" 연결 실패 — 확인할 것")
        print("=" * 78)
        print("""
  1) 토큰이 IQM Resonance 대시보드에서 새로 발급한 것인지
     (구 형식 토큰은 401 을 반환합니다)
  2) 토큰에 공백이나 줄바꿈이 섞이지 않았는지
     ~/.iqm_emerald.sh 를 열어  export IQM_TOKEN=xxxx  형태인지 확인
     (따옴표는 있어도 없어도 됩니다. 값 뒤 공백은 안 됩니다)
  3) 계정에 해당 디바이스 사용 권한이 있는지
     Resonance 대시보드에서 Garnet / Emerald 접근 가능 여부 확인
  4) URL 형식
     신형:  export IQM_SERVER_URL=https://resonance.meetiqm.com
     구형:  export IQM_SERVER_URL=https://cocos.resonance.meetiqm.com/garnet
""")
        return 1

    u, bk = ok_backend
    cm = bk.coupling_map
    edges = sorted({tuple(sorted(e)) for e in cm.get_edges()})
    info = dict(url=u, backend=bk.name, num_qubits=int(bk.num_qubits),
                n_edges=len(edges), basis_gates=list(bk.operation_names),
                edges=[list(e) for e in edges])
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "iqm_live_backend.json"), "w") as fh:
        json.dump(info, fh, indent=2)

    print("\n" + "=" * 78)
    print(" 연결 성공")
    print("=" * 78)
    print(f"  URL      : {u}")
    print(f"  백엔드   : {bk.name}")
    print(f"  큐빗     : {bk.num_qubits}")
    print(f"  커플링   : 엣지 {len(edges)}개")
    print(f"  기본게이트: {bk.operation_names}")
    n16 = "가능" if bk.num_qubits >= 16 else "불가"
    n19 = "가능" if bk.num_qubits >= 19 else "불가"
    print(f"\n  coexp_ising_16  (16 큐빗) : {n16}")
    print(f"  string_maxcut_19(19 큐빗) : {n19}")
    print(f"\n저장 -> {OUT}/iqm_live_backend.json")
    print("\n다음 단계:  python scripts/08_iqm_submit.py --yes --only 1,2")
    return 0


if __name__ == "__main__":
    sys.exit(main())
