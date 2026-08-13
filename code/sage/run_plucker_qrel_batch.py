#!/usr/bin/env python3
"""Run Plucker q-relation probes over several finite-field primes.

Raw Groebner bases are written outside the repository by default.  The compact
q-relation and factor JSON files are written under ``code/sage/out``.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SAGE_DIR = ROOT / "code" / "sage"
OUT_DIR = SAGE_DIR / "out"


@dataclass(frozen=True)
class CaseSpec:
    name: str
    probe_script: str
    basis_prefix: str
    qrel_prefix: str


CASES = {
    "s7_78612": CaseSpec(
        name="s7_78612",
        probe_script="probe_overtie_plucker_locus.py",
        basis_prefix="plucker_locus_s7_78612",
        qrel_prefix="qrel_plucker_locus_s7_78612",
    ),
    "s7_78612_ansatz": CaseSpec(
        name="s7_78612",
        probe_script="probe_overtie_plucker_ansatz.py",
        basis_prefix="plucker_ansatz_s7_78612",
        qrel_prefix="qrel_plucker_ansatz_s7_78612",
    ),
    "s8_79656": CaseSpec(
        name="s8_79656",
        probe_script="probe_overtie_plucker_ansatz.py",
        basis_prefix="plucker_ansatz_s8_79656",
        qrel_prefix="qrel_plucker_ansatz_s8_79656",
    ),
}


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def parse_csv(text: str) -> list[str]:
    return [part.strip() for part in text.split(",") if part.strip()]


def run_command(cmd: list[str], env: dict[str, str], log_path: Path, timeout_seconds: int) -> int:
    started = time.time()
    actual_cmd = ["timeout", str(timeout_seconds), *cmd] if timeout_seconds > 0 else cmd
    print("RUN", " ".join(actual_cmd), flush=True)
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"$ {' '.join(actual_cmd)}\n")
        log.flush()
        proc = subprocess.Popen(
            actual_cmd,
            cwd=ROOT,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="", flush=True)
            log.write(line)
        code = proc.wait()
        elapsed = time.time() - started
        log.write(f"# exit={code} elapsed_seconds={elapsed:.3f}\n")
    print(f"DONE exit={code} elapsed_seconds={elapsed:.3f}", flush=True)
    return code


def append_jsonl(path: Path, record: dict) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", required=True, help="comma-separated primes")
    parser.add_argument(
        "--cases",
        default="s8_79656,s7_78612",
        help="comma-separated case names; default runs the faster s8 first",
    )
    parser.add_argument(
        "--python",
        default="/home/p.osinenko/miniforge3/envs/sage/bin/python",
        help="Sage Python executable",
    )
    parser.add_argument("--dot-sage", default="/tmp/gtz_sage_cache")
    parser.add_argument("--raw-dir", default="/tmp/gtz_plucker")
    parser.add_argument("--log", default="")
    parser.add_argument("--stop-after-hours", type=float, default=0.0)
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=0,
        help="optional timeout for each basis/qrel/factor stage",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--delete-raw-after-qrel", action="store_true")
    args = parser.parse_args()

    primes = [int(part) for part in parse_csv(args.primes)]
    nonprimes = [p for p in primes if not is_prime(p)]
    if nonprimes:
        raise SystemExit(f"not prime: {nonprimes}")

    case_names = parse_csv(args.cases)
    unknown = [name for name in case_names if name not in CASES]
    if unknown:
        raise SystemExit(f"unknown cases: {unknown}; known={sorted(CASES)}")

    raw_dir = Path(args.raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log) if args.log else OUT_DIR / f"plucker_qrel_batch_{int(time.time())}.log"
    jsonl_path = log_path.with_suffix(log_path.suffix + ".jsonl")

    env = os.environ.copy()
    env["DOT_SAGE"] = args.dot_sage

    started = time.time()
    total_failures = 0
    print("=" * 78)
    print("PLUCKER QREL BATCH")
    print(f"root: {ROOT}")
    print(f"primes: {primes}")
    print(f"cases: {case_names}")
    print(f"raw dir: {raw_dir}")
    print(f"log: {log_path}")
    print("=" * 78, flush=True)

    for case_name in case_names:
        spec = CASES[case_name]
        for prime in primes:
            if args.stop_after_hours and (time.time() - started) / 3600 >= args.stop_after_hours:
                print("time limit reached before next target", flush=True)
                return 0 if total_failures == 0 else 1

            basis_json = raw_dir / f"{spec.basis_prefix}_p{prime}.json"
            qrel_json = OUT_DIR / f"{spec.qrel_prefix}_p{prime}.json"
            factors_json = OUT_DIR / f"{spec.qrel_prefix}_p{prime}_factors.json"
            target = {
                "case": case_name,
                "prime": prime,
                "basis_json": str(basis_json),
                "qrel_json": str(qrel_json),
                "factors_json": str(factors_json),
                "started_at": time.time(),
            }
            append_jsonl(jsonl_path, {"event": "target_start", **target})

            if args.overwrite or not qrel_json.exists():
                if args.overwrite or not basis_json.exists():
                    code = run_command(
                        [
                            args.python,
                            str(SAGE_DIR / spec.probe_script),
                            "--case",
                            spec.name,
                            "--characteristic",
                            str(prime),
                            "--out",
                            str(basis_json),
                        ],
                        env,
                        log_path,
                        args.timeout_seconds,
                    )
                    if code != 0:
                        total_failures += 1
                        append_jsonl(jsonl_path, {"event": "target_failed", "stage": "basis", **target})
                        continue

                code = run_command(
                    [
                        args.python,
                        str(SAGE_DIR / "compute_q_power_relation.py"),
                        "--basis-json",
                        str(basis_json),
                        "--variable",
                        "q",
                        "--check-every",
                        "100",
                        "--out",
                        str(qrel_json),
                    ],
                    env,
                    log_path,
                    args.timeout_seconds,
                )
                if code != 0:
                    total_failures += 1
                    append_jsonl(jsonl_path, {"event": "target_failed", "stage": "qrel", **target})
                    continue
            else:
                print(f"SKIP existing {qrel_json}", flush=True)

            if args.overwrite or not factors_json.exists():
                code = run_command(
                    [
                        args.python,
                        str(SAGE_DIR / "factor_q_relation.py"),
                        "--input",
                        str(qrel_json),
                        "--out",
                        str(factors_json),
                    ],
                    env,
                    log_path,
                    args.timeout_seconds,
                )
                if code != 0:
                    total_failures += 1
                    append_jsonl(jsonl_path, {"event": "target_failed", "stage": "factor", **target})
                    continue
            else:
                print(f"SKIP existing {factors_json}", flush=True)

            if args.delete_raw_after_qrel and basis_json.exists():
                basis_json.unlink()
            append_jsonl(
                jsonl_path,
                {
                    "event": "target_done",
                    **target,
                    "elapsed_seconds": time.time() - target["started_at"],
                },
            )

    print("=" * 78)
    print(f"batch complete failures={total_failures}")
    print("=" * 78, flush=True)
    return 0 if total_failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
