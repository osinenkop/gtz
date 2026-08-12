#!/usr/bin/env python3
"""Periodically refresh CRT summaries for completed Plucker q-relations."""
from __future__ import annotations

import argparse
import glob
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "code" / "sage" / "out"
SAGE_PYTHON = "/home/p.osinenko/miniforge3/envs/sage/bin/python"


TARGETS = {
    "s7_78612": "qrel_plucker_locus_s7_78612_p*.json",
    "s8_79656": "qrel_plucker_ansatz_s8_79656_p*.json",
}


def paths_for(pattern: str) -> list[str]:
    paths = [
        path
        for path in glob.glob(str(OUT_DIR / pattern))
        if not path.endswith("_factors.json")
    ]
    return sorted(paths)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hours", type=float, default=5.5)
    parser.add_argument("--interval-seconds", type=int, default=600)
    parser.add_argument("--python", default=SAGE_PYTHON)
    parser.add_argument("--log", default=str(OUT_DIR / "qrel_crt_summary_watch.log"))
    args = parser.parse_args()

    deadline = time.time() + args.hours * 3600
    log_path = Path(args.log)
    while True:
        now = time.strftime("%Y-%m-%d %H:%M:%S")
        with log_path.open("a", encoding="utf-8") as log:
            log.write(f"=== {now} ===\n")
            for case, pattern in TARGETS.items():
                paths = paths_for(pattern)
                log.write(f"{case}: {len(paths)} qrel files\n")
                if len(paths) < 2:
                    continue
                primes = "_".join(Path(path).stem.rsplit("_p", 1)[1] for path in paths)
                out = OUT_DIR / f"qrel_plucker_{case}_crt_{primes}.json"
                cmd = [
                    args.python,
                    str(ROOT / "code" / "sage" / "summarize_qrel_crt.py"),
                    "--out",
                    str(out),
                    *paths,
                ]
                proc = subprocess.run(
                    cmd,
                    cwd=ROOT,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                log.write(proc.stdout)
                log.write(f"exit={proc.returncode}\n")
            log.flush()
        if time.time() >= deadline:
            break
        time.sleep(args.interval_seconds)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
