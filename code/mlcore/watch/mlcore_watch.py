#!/usr/bin/env python3
"""Watch GTZ MLCore jobs and collect separator artifacts.

This watcher intentionally uses only the documented MLCore CLI surface:
`mlc job get`, `mlc job logs`, and `mlc job download artifacts`.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZipFile


DEFAULT_JOBS = [
    "gtz-sage-local-sep-o9fr7w",  # size-7 full extraction, p=32003
    "gtz-sage-local-sep-b321g6",  # size-8 full extraction, p=32003
    "gtz-sage-local-sep-uskwjl",  # size-7 rank-only, p=32009
    "gtz-sage-local-sep-6nm5vj",  # size-7 rank-only, p=32027
    "gtz-sage-local-sep-0l3ssu",  # size-8 rank-only, p=32009
    "gtz-sage-local-sep-lpvaxl",  # size-8 rank-only, p=32027
]


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def as_text(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def run_command(args: list[str], timeout: int) -> dict:
    started = utc_now()
    try:
        completed = subprocess.run(
            args,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
        return {
            "args": args,
            "started_at": started,
            "finished_at": utc_now(),
            "returncode": completed.returncode,
            "stdout": as_text(completed.stdout),
            "stderr": as_text(completed.stderr),
            "timed_out": False,
        }
    except subprocess.TimeoutExpired as exc:
        return {
            "args": args,
            "started_at": started,
            "finished_at": utc_now(),
            "returncode": None,
            "stdout": as_text(exc.stdout),
            "stderr": as_text(exc.stderr),
            "timed_out": True,
        }


def parse_scalar(text: str, key: str) -> str | None:
    match = re.search(rf"^{re.escape(key)}:\s*(.+)$", text, flags=re.MULTILINE)
    if not match:
        return None
    return match.group(1).strip().strip('"')


def parse_job_get(text: str) -> dict:
    return {
        "name": parse_scalar(text, "name"),
        "project": parse_scalar(text, "project"),
        "state": parse_scalar(text, "state"),
        "phase": parse_scalar(text, "phase"),
        "terminationReason": parse_scalar(text, "terminationReason"),
        "createTime": parse_scalar(text, "createTime"),
        "endTime": parse_scalar(text, "endTime"),
        "uiLink": parse_scalar(text, "uiLink"),
    }


def parse_log_tail(text: str) -> dict:
    parsed = {}
    patterns = {
        "matrix": r"matrix:\s*([0-9]+ x [0-9]+)",
        "rank": r"rank:\s*([0-9]+), generator rank:\s*([0-9]+), separator rank defect:\s*([0-9]+)",
        "kernel_dimension": r"kernel dimension:\s*([0-9]+)",
        "candidates": r"separator candidates:\s*([0-9]+), passing at root:\s*([0-9]+)",
        "largest_abs": r"largest \|s\(root\)\|:\s*(.+)",
        "first_passing": r"first passing basis index:\s*([0-9]+)",
        "passing_abs": r"\|s\(root\)\|:\s*(.+)",
        "wrote": r"wrote\s+(.+)",
    }
    matrix = re.search(patterns["matrix"], text)
    if matrix:
        parsed["matrix"] = matrix.group(1)
    rank = re.search(patterns["rank"], text)
    if rank:
        parsed["matrix_rank"] = int(rank.group(1))
        parsed["generator_rank"] = int(rank.group(2))
        parsed["separator_rank_defect"] = int(rank.group(3))
    kernel = re.search(patterns["kernel_dimension"], text)
    if kernel:
        parsed["kernel_dimension"] = int(kernel.group(1))
    candidates = re.search(patterns["candidates"], text)
    if candidates:
        parsed["candidate_count"] = int(candidates.group(1))
        parsed["passing_candidate_count"] = int(candidates.group(2))
    largest = re.search(patterns["largest_abs"], text)
    if largest:
        parsed["max_separator_abs_value"] = largest.group(1).strip()
    first_passing = re.search(patterns["first_passing"], text)
    if first_passing:
        parsed["first_passing_basis_index"] = int(first_passing.group(1))
    passing_abs = re.search(patterns["passing_abs"], text)
    if passing_abs:
        parsed["first_passing_abs_value"] = passing_abs.group(1).strip()
    wrote = re.findall(patterns["wrote"], text)
    if wrote:
        parsed["wrote"] = wrote[-1].strip()
    return parsed


def load_state(path: Path, jobs: list[str]) -> dict:
    if path.exists():
        data = json.loads(path.read_text())
    else:
        data = {"created_at": utc_now(), "jobs": {}}
    for job in jobs:
        data.setdefault("jobs", {}).setdefault(job, {})
    return data


def save_state(path: Path, state: dict) -> None:
    path.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")


def safe_extract_local_separator_jsons(zip_path: Path, out_dir: Path, wanted_basename: str | None = None) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    extracted = []
    with ZipFile(zip_path) as archive:
        for member in archive.namelist():
            basename = Path(member).name
            if wanted_basename and basename != wanted_basename:
                continue
            if not basename.startswith("local_separator") or not basename.endswith(".json"):
                continue
            target = out_dir / basename
            with archive.open(member) as source:
                target.write_bytes(source.read())
            extracted.append(target)
    return extracted


def summarize_separator_json(path: Path) -> dict:
    data = json.loads(path.read_text())
    fields = [
        "target_label",
        "multiplier_degree",
        "separator_degree",
        "characteristic",
        "matrix_rank",
        "generator_rank",
        "separator_rank_defect",
        "kernel_dimension",
        "rank_only",
        "candidate_count",
        "passing_candidate_count",
        "max_separator_abs_value",
        "max_separator_abs_basis_index",
    ]
    summary = {field: data.get(field) for field in fields if field in data}
    summary["file"] = str(path)
    return summary


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def make_markdown(state: dict, jobs: list[str], finished: bool, deadline_reached: bool) -> str:
    lines = [
        "# GTZ MLCore Watch Summary",
        "",
        f"- Updated: `{state.get('updated_at', '')}`",
        f"- Started: `{state.get('created_at', '')}`",
        f"- Poll count: `{state.get('poll_count', 0)}`",
        f"- Finished all jobs: `{finished}`",
        f"- Deadline reached: `{deadline_reached}`",
        "",
        "## Job Status",
        "",
        "| Job | State | Phase | Termination | Log signal |",
        "| --- | --- | --- | --- | --- |",
    ]
    for job in jobs:
        info = state["jobs"].get(job, {})
        status = info.get("status", {})
        signal = info.get("log_signal", {})
        signal_bits = []
        if "separator_rank_defect" in signal:
            signal_bits.append(f"defect {signal['separator_rank_defect']}")
        if "passing_candidate_count" in signal:
            signal_bits.append(f"passing {signal['passing_candidate_count']}")
        if "max_separator_abs_value" in signal:
            signal_bits.append(f"max |s| {signal['max_separator_abs_value']}")
        if "wrote" in signal:
            signal_bits.append(f"`{signal['wrote']}`")
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{job}`",
                    f"`{status.get('state', '')}`",
                    f"`{status.get('phase', '')}`",
                    f"`{status.get('terminationReason', '')}`",
                    "; ".join(signal_bits),
                ]
            )
            + " |"
        )

    separator_summaries = []
    seen_separator_files = set()
    for job in jobs:
        for summary in state["jobs"].get(job, {}).get("separator_jsons", []):
            filename = summary.get("file", "")
            if filename in seen_separator_files:
                continue
            seen_separator_files.add(filename)
            separator_summaries.append(summary)

    if separator_summaries:
        lines.extend(
            [
                "",
                "## Separator JSONs",
                "",
                "| File | p | D | S | rank | gen rank | defect | kernel | rank-only | candidates | passing | max |s(root)| |",
                "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | --- |",
            ]
        )
        for summary in separator_summaries:
            lines.append(
                "| "
                + " | ".join(
                    [
                        f"`{Path(summary.get('file', '')).name}`",
                        str(summary.get("characteristic", "")),
                        str(summary.get("multiplier_degree", "")),
                        str(summary.get("separator_degree", "")),
                        str(summary.get("matrix_rank", "")),
                        str(summary.get("generator_rank", "")),
                        str(summary.get("separator_rank_defect", "")),
                        str(summary.get("kernel_dimension", "")),
                        str(summary.get("rank_only", "")),
                        str(summary.get("candidate_count", "")),
                        str(summary.get("passing_candidate_count", "")),
                        str(summary.get("max_separator_abs_value", "")),
                    ]
                )
                + " |"
            )

    breakthrough = []
    completed_full = []
    zero_full = []
    for summary in separator_summaries:
        if summary.get("rank_only") is False:
            completed_full.append(summary)
            passing = summary.get("passing_candidate_count")
            if isinstance(passing, int) and passing > 0:
                breakthrough.append(summary)
            elif passing == 0:
                zero_full.append(summary)

    lines.extend(["", "## Conclusion", ""])
    if breakthrough:
        lines.append("Potential breakthrough: at least one full extraction has a nonzero-at-root separator candidate.")
    elif completed_full and len(zero_full) == len(completed_full):
        lines.append("All completed full extractions found zero passing kernel-basis separators so far.")
    elif deadline_reached and not finished:
        lines.append("Deadline reached before all jobs completed. Continue polling or inspect running jobs manually.")
    elif finished:
        lines.append("All watched jobs ended. No passing full-extraction candidate was observed in extracted summaries.")
    else:
        lines.append("Still waiting for running or pending jobs.")

    lines.append("")
    return "\n".join(lines)


def poll_once(args: argparse.Namespace, state: dict, jobs: list[str], root: Path) -> None:
    watch_dir = root / "code/mlcore/watch"
    raw_dir = watch_dir / "raw"
    log_dir = watch_dir / "logs"
    downloads_dir = root / "code/mlcore/downloads"
    sage_out_dir = root / "code/sage/out"

    state["poll_count"] = int(state.get("poll_count", 0)) + 1
    state["updated_at"] = utc_now()

    for job in jobs:
        job_state = state["jobs"].setdefault(job, {})
        job_raw_dir = raw_dir / job
        job_raw_dir.mkdir(parents=True, exist_ok=True)

        get_result = run_command(["mlc", "job", "get", job, "-p", args.project], args.command_timeout)
        write_text(job_raw_dir / "job_get_latest.txt", get_result["stdout"] + get_result["stderr"])
        job_state["last_get"] = {
            "returncode": get_result["returncode"],
            "timed_out": get_result["timed_out"],
            "finished_at": get_result["finished_at"],
        }
        if get_result["returncode"] == 0:
            job_state["status"] = parse_job_get(get_result["stdout"])

        logs_result = run_command(["mlc", "job", "logs", job, "-p", args.project], args.command_timeout)
        write_text(log_dir / f"{job}.log", logs_result["stdout"] + logs_result["stderr"])
        job_state["last_logs"] = {
            "returncode": logs_result["returncode"],
            "timed_out": logs_result["timed_out"],
            "finished_at": logs_result["finished_at"],
        }
        if logs_result["stdout"]:
            job_state["log_signal"] = parse_log_tail(logs_result["stdout"])

        status = job_state.get("status", {})
        if status.get("state") == "SUCCEEDED" and not job_state.get("artifacts_downloaded"):
            dst = downloads_dir / job
            dst.mkdir(parents=True, exist_ok=True)
            download_result = run_command(
                [
                    "mlc",
                    "job",
                    "download",
                    "artifacts",
                    job,
                    "-p",
                    args.project,
                    "--dst-folder",
                    str(dst),
                ],
                args.download_timeout,
            )
            write_text(job_raw_dir / "artifact_download_latest.txt", download_result["stdout"] + download_result["stderr"])
            job_state["last_artifact_download"] = {
                "returncode": download_result["returncode"],
                "timed_out": download_result["timed_out"],
                "finished_at": download_result["finished_at"],
            }
            if download_result["returncode"] == 0:
                job_state["artifacts_downloaded"] = True

        zip_path = downloads_dir / job / "GTZ_OUT.zip"
        if zip_path.exists() and not job_state.get("separator_jsons_extracted"):
            wrote = job_state.get("log_signal", {}).get("wrote", "")
            wanted_basename = Path(wrote).name if wrote else None
            extracted = safe_extract_local_separator_jsons(zip_path, sage_out_dir, wanted_basename)
            if wanted_basename and not extracted:
                extracted = safe_extract_local_separator_jsons(zip_path, sage_out_dir)
            job_state["extracted_files"] = [str(path) for path in extracted]
            job_state["separator_jsons"] = [summarize_separator_json(path) for path in extracted]
            job_state["separator_jsons_extracted"] = True


def all_jobs_ended(state: dict, jobs: list[str]) -> bool:
    for job in jobs:
        status = state["jobs"].get(job, {}).get("status", {})
        if status.get("phase") != "ENDED" and status.get("state") not in {"SUCCEEDED", "FAILED", "CANCELLED"}:
            return False
    return True


def has_breakthrough(state: dict, jobs: list[str]) -> bool:
    for job in jobs:
        for summary in state["jobs"].get(job, {}).get("separator_jsons", []):
            if summary.get("rank_only") is False and int(summary.get("passing_candidate_count") or 0) > 0:
                return True
    return False


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", default="aida")
    parser.add_argument("--jobs", nargs="*", default=DEFAULT_JOBS)
    parser.add_argument("--interval-seconds", type=int, default=600)
    parser.add_argument("--hours", type=float, default=10.0)
    parser.add_argument("--root", default=".")
    parser.add_argument(
        "--watch-dir",
        default="",
        help="directory for watcher state/logs; defaults to code/mlcore/watch under --root",
    )
    parser.add_argument(
        "--continue-after-breakthrough",
        action="store_true",
        help="keep polling after a nonzero full-extraction candidate is observed",
    )
    parser.add_argument(
        "--continue-until-deadline",
        action="store_true",
        help="keep polling until the time deadline even if all watched jobs have ended",
    )
    parser.add_argument("--command-timeout", type=int, default=180)
    parser.add_argument("--download-timeout", type=int, default=600)
    args = parser.parse_args()

    root = Path(args.root).resolve()
    watch_dir = Path(args.watch_dir).resolve() if args.watch_dir else root / "code/mlcore/watch"
    watch_dir.mkdir(parents=True, exist_ok=True)
    state_path = watch_dir / "state.json"
    latest_summary = watch_dir / "latest_summary.md"
    final_summary = watch_dir / "final_summary.md"
    jobs = list(args.jobs)

    state = load_state(state_path, jobs)
    state["project"] = args.project
    state["jobs_watched"] = jobs
    state["watch_started_at"] = state.get("watch_started_at") or utc_now()
    state["max_hours"] = args.hours
    state["interval_seconds"] = args.interval_seconds
    state.pop("stop_reason", None)
    state.pop("stopped_at", None)

    deadline = time.monotonic() + args.hours * 3600
    deadline_reached = False
    finished = False

    while True:
        poll_once(args, state, jobs, root)
        finished = all_jobs_ended(state, jobs)
        deadline_reached = time.monotonic() >= deadline
        markdown = make_markdown(state, jobs, finished, deadline_reached)
        write_text(latest_summary, markdown)
        save_state(state_path, state)

        if has_breakthrough(state, jobs) and not args.continue_after_breakthrough:
            state["stop_reason"] = "breakthrough"
            break
        if finished and not args.continue_until_deadline:
            state["stop_reason"] = "all_jobs_ended"
            break
        if deadline_reached:
            state["stop_reason"] = "deadline_reached"
            break
        time.sleep(args.interval_seconds)

    state["stopped_at"] = utc_now()
    final_markdown = make_markdown(state, jobs, finished, deadline_reached)
    write_text(latest_summary, final_markdown)
    write_text(final_summary, final_markdown)
    save_state(state_path, state)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
