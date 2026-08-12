#!/usr/bin/env python3
"""Run GTZ D6/S6 reconstruction attempts as watched artifacts arrive."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


MAMBA = Path("/home/p.osinenko/miniforge3/bin/mamba")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def rel(path: Path, root: Path) -> str:
    return str(path.relative_to(root))


def tail(text: str, lines: int = 12) -> str:
    return "\n".join(text.splitlines()[-lines:])


def load_state(path: Path) -> dict:
    if not path.exists():
        return {"created_at": utc_now(), "attempts": {}}
    return json.loads(path.read_text())


def save_state(path: Path, state: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")


def case_specs(root: Path) -> list[dict]:
    out = root / "code/sage/out"
    return [
        {
            "name": "s7_support19_limited",
            "min_primes": 4,
            "root_file": out / "refine_overtie_s7_78612_t2_7_10_19_p1400.json",
            "prefix": out / "reconstruct_s7_78612_D6_S6_support19_limited",
            "inputs": [
                out / "local_separator_s7_78612_omit_active1_D6_S6_32003_support19.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32009.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32027_support19.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32029_support19.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32051_support19.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32057_support19.json",
            ],
        },
        {
            "name": "s7_support8_limited",
            "min_primes": 3,
            "root_file": out / "refine_overtie_s7_78612_t2_7_10_19_p1400.json",
            "prefix": out / "reconstruct_s7_78612_D6_S6_support8_limited",
            "inputs": [
                out / "local_separator_s7_78612_omit_active1_D6_S6_32003_support8.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32009_support8.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32027_support8.json",
            ],
        },
        {
            "name": "s8_support6_limited",
            "min_primes": 4,
            "root_file": out / "refine_overtie_s8_79656_t10_11_15_p1400.json",
            "prefix": out / "reconstruct_s8_79656_D6_S6_support6_limited",
            "inputs": [
                out / "local_separator_s8_79656_omit_active0_D6_S6_32003_support6.json",
                out / "local_separator_s8_79656_omit_active0_D6_S6_32009_support6.json",
                out / "local_separator_s8_79656_omit_active0_D6_S6_32027_support6.json",
                out / "local_separator_s8_79656_omit_active0_D6_S6_32029_support6.json",
                out / "local_separator_s8_79656_omit_active0_D6_S6_32051_support6.json",
                out / "local_separator_s8_79656_omit_active0_D6_S6_32057_support6.json",
            ],
        },
        {
            "name": "s7_support19_allcand",
            "min_primes": 3,
            "root_file": out / "refine_overtie_s7_78612_t2_7_10_19_p1400.json",
            "prefix": out / "reconstruct_s7_78612_D6_S6_support19_allcand",
            "short_prefix": out / "short_s7_78612_D6_S6_support19_allcand",
            "inputs": [
                out / "local_separator_s7_78612_omit_active1_D6_S6_32003_support19_allcand.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32009_support19_allcand.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32027_support19_allcand.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32029_support19_allcand.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32051_support19_allcand.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32057_support19_allcand.json",
            ],
        },
        {
            "name": "s7_support0_allcand",
            "min_primes": 3,
            "root_file": out / "refine_overtie_s7_78612_t2_7_10_19_p1400.json",
            "prefix": out / "reconstruct_s7_78612_D6_S6_support0_allcand",
            "short_prefix": out / "short_s7_78612_D6_S6_support0_allcand",
            "inputs": [
                out / "local_separator_s7_78612_omit_active1_D6_S6_32003_support0_allcand.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32009_support0_allcand.json",
                out / "local_separator_s7_78612_omit_active1_D6_S6_32027_support0_allcand.json",
            ],
        },
        {
            "name": "s8_support6_allcand",
            "min_primes": 3,
            "root_file": out / "refine_overtie_s8_79656_t10_11_15_p1400.json",
            "prefix": out / "reconstruct_s8_79656_D6_S6_support6_allcand",
            "short_prefix": out / "short_s8_79656_D6_S6_support6_allcand",
            "inputs": [
                out / "local_separator_s8_79656_omit_active0_D6_S6_32003_support6_allcand.json",
                out / "local_separator_s8_79656_omit_active0_D6_S6_32009_support6_allcand.json",
                out / "local_separator_s8_79656_omit_active0_D6_S6_32027_support6_allcand.json",
                out / "local_separator_s8_79656_omit_active0_D6_S6_32029_support6_allcand.json",
                out / "local_separator_s8_79656_omit_active0_D6_S6_32051_support6_allcand.json",
                out / "local_separator_s8_79656_omit_active0_D6_S6_32057_support6_allcand.json",
            ],
        },
    ]


def prime_from_path(path: Path) -> str:
    match = re.search(r"_(3[0-9]{4})(?:_|\.json)", path.name)
    if not match:
        raise ValueError(f"could not parse prime from {path}")
    return match.group(1)


def available_inputs(spec: dict) -> list[Path]:
    return [path for path in spec["inputs"] if path.exists()]


def prime_tag(paths: list[Path]) -> str:
    return "_".join(prime_from_path(path) for path in paths)


def attempt_key(spec: dict, paths: list[Path]) -> str:
    return f"{spec['name']}_{prime_tag(paths)}"


def short_attempt_key(spec: dict, paths: list[Path]) -> str:
    return f"{spec['name']}_short_{prime_tag(paths)}"


def reconstruction_command(root: Path, spec: dict, inputs_paths: list[Path]) -> list[str]:
    inputs = [rel(path, root) for path in inputs_paths]
    out_prefix = rel(Path(f"{spec['prefix']}_{prime_tag(inputs_paths)}_p1400"), root)
    return [
        str(MAMBA),
        "run",
        "-n",
        "sage",
        "python",
        "code/sage/reconstruct_separator_subspace.py",
        "--inputs",
        *inputs,
        "--root",
        rel(spec["root_file"], root),
        "--precision",
        "500",
        "--out-prefix",
        out_prefix,
    ]


def short_recovery_command(root: Path, spec: dict, inputs_paths: list[Path]) -> list[str]:
    inputs = [rel(path, root) for path in inputs_paths]
    out_prefix = rel(Path(f"{spec['short_prefix']}_p{prime_tag(inputs_paths)}_p1400"), root)
    return [
        str(MAMBA),
        "run",
        "-n",
        "sage",
        "python",
        "code/sage/recover_short_separator_lattice.py",
        "--inputs",
        *inputs,
        "--root",
        rel(spec["root_file"], root),
        "--precision",
        "1400",
        "--max-output",
        "80",
        "--out-prefix",
        out_prefix,
    ]


def run_command_attempt(
    root: Path,
    run_dir: Path,
    spec: dict,
    inputs_paths: list[Path],
    command: list[str],
    attempt_type: str,
    summary_json: Path,
) -> dict:
    started = utc_now()
    tag = prime_tag(inputs_paths)
    log_path = run_dir / "logs" / f"{spec['name']}_{attempt_type}_{tag}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        command,
        cwd=root,
        text=True,
        capture_output=True,
        check=False,
    )
    output = completed.stdout + completed.stderr
    log_path.write_text(output)
    return {
        "started_at": started,
        "finished_at": utc_now(),
        "returncode": completed.returncode,
        "spec_name": spec["name"],
        "attempt_type": attempt_type,
        "prime_tag": tag,
        "display_tag": f"{attempt_type}:{tag}",
        "prime_count": len(inputs_paths),
        "command": command,
        "log": str(log_path),
        "tail": tail(output),
        "summary_json": rel(summary_json, root),
    }


def run_attempt(root: Path, run_dir: Path, spec: dict, inputs_paths: list[Path]) -> dict:
    tag = prime_tag(inputs_paths)
    summary_json = Path(f"{spec['prefix']}_{tag}_p1400.json")
    return run_command_attempt(
        root,
        run_dir,
        spec,
        inputs_paths,
        reconstruction_command(root, spec, inputs_paths),
        "rref",
        summary_json,
    )


def run_short_attempt(root: Path, run_dir: Path, spec: dict, inputs_paths: list[Path]) -> dict:
    tag = prime_tag(inputs_paths)
    summary_json = Path(f"{spec['short_prefix']}_p{tag}_p1400.json")
    return run_command_attempt(
        root,
        run_dir,
        spec,
        inputs_paths,
        short_recovery_command(root, spec, inputs_paths),
        "short",
        summary_json,
    )


def write_summary(path: Path, state: dict, specs: list[dict]) -> None:
    lines = [
        "# GTZ D6/S6 Reconstruction Follow-up",
        "",
        f"- Updated: `{state.get('updated_at', '')}`",
        f"- Started: `{state.get('created_at', '')}`",
        f"- Poll count: `{state.get('poll_count', 0)}`",
        "",
        "| Case | Available | Attempted | Latest result |",
        "| --- | ---: | --- | --- |",
    ]
    attempts = state.get("attempts", {})
    for spec in specs:
        available = available_inputs(spec)
        available_label = prime_tag(available) if available else "-"
        spec_attempts = [
            (key, value)
            for key, value in attempts.items()
            if value.get("spec_name") == spec["name"] or key.startswith(f"{spec['name']}_")
        ]
        spec_attempts.sort(key=lambda item: item[1].get("finished_at", ""))
        tried = [
            value.get("display_tag") or value.get("prime_tag") or key.removeprefix(f"{spec['name']}_")
            for key, value in spec_attempts
        ]
        latest = spec_attempts[-1][1] if spec_attempts else None
        if latest is None:
            result = "not attempted"
        elif latest["returncode"] == 0:
            result = f"succeeded: `{latest['summary_json']}`"
        else:
            result = f"failed, see `{latest['log']}`"
        total_inputs = len(spec["inputs"])
        lines.append(
            f"| `{spec['name']}` | {len(available)}/{total_inputs} `{available_label}` | "
            f"{', '.join(str(item) for item in tried) if tried else '-'} | {result} |"
        )
    lines.append("")
    path.write_text("\n".join(lines))


def poll_once(root: Path, run_dir: Path, state: dict) -> None:
    specs = case_specs(root)
    state["poll_count"] = int(state.get("poll_count", 0)) + 1
    state["updated_at"] = utc_now()
    attempts = state.setdefault("attempts", {})

    for spec in specs:
        inputs_paths = available_inputs(spec)
        if len(inputs_paths) < spec["min_primes"]:
            continue
        key = attempt_key(spec, inputs_paths)
        if key in attempts:
            pass
        else:
            attempts[key] = run_attempt(root, run_dir, spec, inputs_paths)
            save_state(run_dir / "state.json", state)
        if "short_prefix" in spec:
            short_key = short_attempt_key(spec, inputs_paths)
            if short_key not in attempts:
                attempts[short_key] = run_short_attempt(root, run_dir, spec, inputs_paths)
                save_state(run_dir / "state.json", state)

    write_summary(run_dir / "latest_summary.md", state, specs)
    save_state(run_dir / "state.json", state)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default="/home/p.osinenko/Documents/gtz")
    parser.add_argument("--run-dir", default="")
    parser.add_argument("--interval-seconds", type=int, default=600)
    parser.add_argument("--hours", type=float, default=4.5)
    args = parser.parse_args()

    root = Path(args.root).resolve()
    run_dir = Path(args.run_dir).resolve() if args.run_dir else root / "code/mlcore/watch/d6_followup"
    run_dir.mkdir(parents=True, exist_ok=True)
    state_path = run_dir / "state.json"
    state = load_state(state_path)
    state["root"] = str(root)
    state["interval_seconds"] = args.interval_seconds
    state["max_hours"] = args.hours

    deadline = time.monotonic() + args.hours * 3600
    while True:
        poll_once(root, run_dir, state)
        if time.monotonic() >= deadline:
            state["stopped_at"] = utc_now()
            state["stop_reason"] = "deadline_reached"
            save_state(state_path, state)
            write_summary(run_dir / "final_summary.md", state, case_specs(root))
            write_summary(run_dir / "latest_summary.md", state, case_specs(root))
            return 0
        time.sleep(args.interval_seconds)


if __name__ == "__main__":
    raise SystemExit(main())
