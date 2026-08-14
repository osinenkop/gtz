#!/usr/bin/env python3
"""Unattended relaxed-bound attack loop.

This runner repeatedly:

1. runs the reduced-obstruction SLSQP sweep ``v40`` at selected targets;
2. runs the relaxed KKT diagnostic ``v43`` on each best candidate;
3. canonicalizes active patterns with ``v42``;
4. updates a compact status JSON/Markdown file.

The goal is not to prove a theorem by itself.  It is an unattended search for
either a lower boundary than the parametric ansatz or a stable active-pattern
profile suitable for exact KKT/CAD certification.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def parse_target(text: str) -> Fraction:
    return Fraction(text)


def target_label(target: Fraction) -> str:
    if target.denominator == 1:
        return str(target.numerator)
    return f"{target.numerator}over{target.denominator}"


def predicted_ansatz_slack(target: Fraction) -> float | None:
    t = float(target)
    disc = 25.0 * t * t + 80.0 * t - 8.0
    if disc < 0:
        return None
    q = (2.0 + 5.0 * t - math.sqrt(disc)) / 6.0
    return q - t


def run_command(cmd: list[str], log_path: Path) -> tuple[int, float]:
    start = time.time()
    with open(log_path, "a") as log:
        log.write(f"\n[{utc_now()}] RUN {' '.join(cmd)}\n")
        log.flush()
        proc = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, text=True)
        elapsed = time.time() - start
        log.write(f"[{utc_now()}] EXIT {proc.returncode} elapsed={elapsed:.1f}s\n")
        log.flush()
    return proc.returncode, elapsed


def load_best(path: Path) -> dict:
    data = json.load(open(path))
    target = data["targets"][0]
    return target["best"]


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)


def write_markdown(path: Path, state: dict) -> None:
    lines = [
        "# Relaxed Long Attack Status",
        "",
        f"- started_utc: `{state['started_utc']}`",
        f"- updated_utc: `{state['updated_utc']}`",
        f"- hours_requested: `{state['hours_requested']}`",
        f"- cycles_completed: `{state['cycles_completed']}`",
        f"- current_cycle: `{state.get('current_cycle')}`",
        f"- current_target: `{state.get('current_target')}`",
        f"- status: `{state['status']}`",
        "",
        "## Best Rows",
        "",
        "| target | best s | ansatz s | gap best-ansatz | bound value | active obj | active core | KKT residual | source |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in state.get("best_by_target", []):
        ansatz = row.get("ansatz_s")
        gap = row.get("gap_to_ansatz")
        lines.append(
            f"| `{row['target']}` | `{row['s']:.12g}` | "
            f"`{ansatz:.12g}` | `{gap:.3e}` | `{row['bound_value']:.12g}` | "
            f"`{row.get('active_objective_count')}` | `{row.get('active_core_count')}` | "
            f"`{row.get('kkt_residual_norm'):.3e}` | `{row['path']}` |"
        )
    lines.extend([
        "",
        "## Events",
        "",
    ])
    for event in state.get("events", [])[-30:]:
        lines.append(f"- `{event['time']}` {event['message']}")
    path.write_text("\n".join(lines) + "\n")


def summarize_kkt(path: Path) -> dict:
    data = json.load(open(path))
    diag = data["diagnostic"]
    return {
        "active_objective_count": diag["active_objective_count"],
        "active_core_count": diag["active_core_count"],
        "kkt_residual_norm": diag["kkt_residual_norm"],
        "beta_min": diag["beta_min"],
        "mu_min": diag["mu_min"],
    }


def update_best(state: dict, target: Fraction, sweep_path: Path, kkt_path: Path) -> None:
    best = load_best(sweep_path)
    kkt = summarize_kkt(kkt_path)
    ansatz = predicted_ansatz_slack(target)
    row = {
        "target": str(target),
        "target_float": float(target),
        "s": float(best["s"]),
        "bound_value": float(best["bound_value"]),
        "F": float(best["F"]),
        "ansatz_s": ansatz,
        "gap_to_ansatz": float(best["s"]) - ansatz if ansatz is not None else None,
        "path": str(sweep_path),
        "kkt_path": str(kkt_path),
        **kkt,
    }

    rows = {item["target"]: item for item in state.get("best_by_target", [])}
    old = rows.get(str(target))
    if old is None or row["s"] < old["s"]:
        rows[str(target)] = row
        state.setdefault("events", []).append({
            "time": utc_now(),
            "message": (
                f"new best target={target} s={row['s']:.12g} "
                f"gap_to_ansatz={row['gap_to_ansatz']:.3e}"
            ),
        })
    state["best_by_target"] = [rows[key] for key in sorted(rows, key=lambda item: float(Fraction(item)))]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hours", type=float, default=8.25)
    parser.add_argument("--targets", default="1/8,0.126,0.127,16/125,0.1284,9/70")
    parser.add_argument("--random-starts", type=int, default=160)
    parser.add_argument("--maxiter", type=int, default=2500)
    parser.add_argument("--seed", type=int, default=2026081501)
    parser.add_argument("--out-dir", default="")
    args = parser.parse_args()

    started = time.time()
    started_stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out_dir = Path(args.out_dir) if args.out_dir else ROOT / "verify" / "out" / f"v44_relaxed_long_attack_{started_stamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "runner.log"
    state_path = out_dir / "state.json"
    md_path = out_dir / "STATUS.md"
    targets = [parse_target(item.strip()) for item in args.targets.split(",") if item.strip()]
    deadline = started + args.hours * 3600.0

    state = {
        "status": "running",
        "started_utc": utc_now(),
        "updated_utc": utc_now(),
        "hours_requested": args.hours,
        "targets": [str(t) for t in targets],
        "random_starts": args.random_starts,
        "maxiter": args.maxiter,
        "cycles_completed": 0,
        "events": [],
        "best_by_target": [],
        "out_dir": str(out_dir),
        "log_path": str(log_path),
    }
    write_json(state_path, state)
    write_markdown(md_path, state)

    cycle = 0
    while time.time() < deadline:
        cycle += 1
        state["current_cycle"] = cycle
        cycle_paths: list[Path] = []
        for idx, target in enumerate(targets):
            if time.time() >= deadline:
                break
            label = target_label(target)
            state["current_target"] = str(target)
            state["updated_utc"] = utc_now()
            write_json(state_path, state)
            write_markdown(md_path, state)

            sweep_path = out_dir / f"cycle{cycle:04d}_{label}_sweep.json"
            kkt_path = out_dir / f"cycle{cycle:04d}_{label}_kkt.json"
            seed = args.seed + 1000 * cycle + idx
            cmd = [
                sys.executable,
                "verify/v40_relaxed_target_sweep.py",
                "--targets",
                str(target),
                "--random-starts",
                str(args.random_starts),
                "--seed",
                str(seed),
                "--maxiter",
                str(args.maxiter),
                "--out",
                str(sweep_path),
            ]
            code, elapsed = run_command(cmd, log_path)
            state.setdefault("events", []).append({
                "time": utc_now(),
                "message": f"cycle={cycle} target={target} v40 exit={code} elapsed={elapsed:.1f}s",
            })
            if code != 0:
                state["updated_utc"] = utc_now()
                write_json(state_path, state)
                write_markdown(md_path, state)
                continue

            cycle_paths.append(sweep_path)
            cmd = [
                sys.executable,
                "verify/v43_relaxed_kkt_diagnostic.py",
                "--input",
                str(sweep_path),
                "--target-index",
                "0",
                "--out",
                str(kkt_path),
            ]
            code, elapsed = run_command(cmd, log_path)
            state.setdefault("events", []).append({
                "time": utc_now(),
                "message": f"cycle={cycle} target={target} v43 exit={code} elapsed={elapsed:.1f}s",
            })
            if code == 0:
                update_best(state, target, sweep_path, kkt_path)

            state["updated_utc"] = utc_now()
            write_json(state_path, state)
            write_markdown(md_path, state)

        if cycle_paths:
            summary_path = out_dir / f"cycle{cycle:04d}_patterns.json"
            cmd = [
                sys.executable,
                "verify/v42_relaxed_active_patterns.py",
                *[str(path) for path in cycle_paths],
                "--out",
                str(summary_path),
            ]
            code, elapsed = run_command(cmd, log_path)
            state.setdefault("events", []).append({
                "time": utc_now(),
                "message": f"cycle={cycle} v42 exit={code} elapsed={elapsed:.1f}s",
            })

        state["cycles_completed"] = cycle
        state["updated_utc"] = utc_now()
        write_json(state_path, state)
        write_markdown(md_path, state)

    state["status"] = "complete"
    state["current_target"] = None
    state["updated_utc"] = utc_now()
    write_json(state_path, state)
    write_markdown(md_path, state)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
