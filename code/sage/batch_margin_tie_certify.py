#!/usr/bin/env python3
"""Batch high-precision refinement and Krawczyk certification of v30 roots.

This is a driver, not a new mathematical method.  For each saved
``verify/out/v30_margin_*`` point it:

1. recomputes the outside triples tied at the numerical margin maximum;
2. chooses subsets of those ties of size ``10-|A|`` so that the active-plus-tie
   determinant system is square in ``(z0,...,z8,q)``;
3. runs ``refine_margin_tie_root.py``;
4. runs ``certify_margin_tie_krawczyk.py``.

If there are more numerical ties than needed, the unselected ties are passed as
``--nonstrict-outside-indices``.  The certificate is only accepted if interval
arithmetic proves they still satisfy det(6P_TT-qI) <= 0.
"""
from __future__ import annotations

import argparse
import glob
import itertools
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from gtz63_semialgebraic import TRIPLES, parse_indices  # noqa: E402
from verify.v33_margin_table import summarise  # noqa: E402


def parse_ints(text: str) -> set[int]:
    if not text.strip():
        return set()
    return {int(part.strip()) for part in text.split(",") if part.strip()}


def tag(indices) -> str:
    return "_".join(str(i) for i in indices)


def tlabel(idx: int) -> str:
    return "".join(str(x) for x in TRIPLES[int(idx)])


def run_command(cmd, timeout, dry_run):
    if dry_run:
        return {
            "returncode": 0,
            "stdout": "",
            "stderr": "",
            "timeout": False,
            "cmd": cmd,
        }
    try:
        proc = subprocess.run(
            cmd,
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )
        return {
            "returncode": proc.returncode,
            "stdout": proc.stdout[-4000:],
            "stderr": proc.stderr[-4000:],
            "timeout": False,
            "cmd": cmd,
        }
    except subprocess.TimeoutExpired as exc:
        return {
            "returncode": None,
            "stdout": (exc.stdout or "")[-4000:] if isinstance(exc.stdout, str) else "",
            "stderr": (exc.stderr or "")[-4000:] if isinstance(exc.stderr, str) else "",
            "timeout": True,
            "cmd": cmd,
        }


def load_json(path):
    return json.loads(Path(path).read_text())


def maybe_run_refine(row, ties, args):
    tie_text = ",".join(str(i) for i in ties)
    base = f"s{row['size']}_{row['canon']}_t{tag(ties)}"
    out = Path(args.out_dir) / f"refine_tie_{base}_p{args.refine_precision}.json"
    cmd = [
        sys.executable,
        str(HERE / "refine_margin_tie_root.py"),
        "--input",
        row["file"],
        "--size",
        str(row["size"]),
        "--canon",
        str(row["canon"]),
        "--tie-indices",
        tie_text,
        "--precision",
        str(args.refine_precision),
        "--steps",
        str(args.refine_steps),
        "--tol-exp",
        str(args.refine_tol_exp),
        "--algdep-degree",
        str(args.algdep_degree),
        "--algdep-height",
        str(args.algdep_height),
        "--out",
        str(out),
    ]
    if out.exists() and not args.force:
        return out, {"skipped_existing": True, "cmd": cmd, "returncode": 0}
    result = run_command(cmd, args.command_timeout, args.dry_run)
    return out, result


def maybe_run_certify(row, ties, nonstrict_outside, refine_out, args):
    base = f"s{row['size']}_{row['canon']}_t{tag(ties)}"
    out = Path(args.out_dir) / f"krawczyk_tie_{base}.json"
    cmd = [
        sys.executable,
        str(HERE / "certify_margin_tie_krawczyk.py"),
        "--input",
        str(refine_out),
        "--precision",
        str(args.certify_precision),
        "--radius-exp-start",
        str(args.radius_exp_start),
        "--radius-exp-stop",
        str(args.radius_exp_stop),
        "--radius-exp-step",
        str(args.radius_exp_step),
        "--out",
        str(out),
    ]
    if nonstrict_outside:
        cmd.extend([
            "--nonstrict-outside-indices",
            ",".join(str(i) for i in nonstrict_outside),
        ])
    if out.exists() and not args.force:
        return out, {"skipped_existing": True, "cmd": cmd, "returncode": 0}
    result = run_command(cmd, args.command_timeout, args.dry_run)
    return out, result


def certificate_status(path):
    if not Path(path).exists():
        return {"exists": False, "certified": False}
    data = load_json(path)
    winner = data.get("winner") or {}
    branch = winner.get("branch_checks") or {}
    return {
        "exists": True,
        "krawczyk_certified": bool(data.get("krawczyk_certified", data.get("certified"))),
        "branch_certified": bool(data.get("branch_certified", data.get("certified"))),
        "certified": bool(data.get("certified")),
        "margin_lower": branch.get("margin_lower"),
        "margin_upper": branch.get("margin_upper"),
        "strict_outside_det_max_upper": branch.get("outside_shifted_det_max_upper"),
        "nonstrict_outside_det_max_upper": branch.get("nonstrict_outside_shifted_det_max_upper"),
        "branch_checks": {
            key: branch.get(key)
            for key in (
                "chart_d_positive",
                "q_gt_one",
                "active_minors_positive",
                "equation_tie_minors_positive",
                "outside_shifted_dets_negative",
                "nonstrict_outside_shifted_dets_nonpositive",
                "all_branch_checks",
            )
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--glob", default="verify/out/v30_margin_*.json")
    parser.add_argument("--sizes", default="6,7,8")
    parser.add_argument("--canons", default="")
    parser.add_argument("--tol", type=float, default=1e-6)
    parser.add_argument("--out-dir", default="code/sage/out")
    parser.add_argument("--summary-out", default="code/sage/out/batch_margin_tie_certificates.json")
    parser.add_argument("--refine-precision", type=int, default=500)
    parser.add_argument("--refine-steps", type=int, default=30)
    parser.add_argument("--refine-tol-exp", type=int, default=260)
    parser.add_argument("--certify-precision", type=int, default=750)
    parser.add_argument("--radius-exp-start", type=int, default=30)
    parser.add_argument("--radius-exp-stop", type=int, default=130)
    parser.add_argument("--radius-exp-step", type=int, default=10)
    parser.add_argument("--algdep-degree", type=int, default=0)
    parser.add_argument("--algdep-height", type=int, default=1000000)
    parser.add_argument("--command-timeout", type=int, default=180)
    parser.add_argument("--max-subsets", type=int, default=0)
    parser.add_argument("--all-subsets-after-success", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    sizes = parse_ints(args.sizes)
    canons = parse_ints(args.canons)
    rows = [summarise(path, args.tol) for path in sorted(glob.glob(args.glob))]
    rows = [row for row in rows if row and row["size"] in sizes]
    if canons:
        rows = [row for row in rows if row["canon"] in canons]

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    results = []
    print("=" * 78)
    print("BATCH MARGIN TIE CERTIFICATION")
    print(f"cases: {len(rows)}")
    print("=" * 78, flush=True)
    for row in rows:
        needed = 10 - len(row["active"])
        ties = list(row["max_tie_out"])
        print(
            f"({row['size']},{row['canon']}): active={len(row['active'])} "
            f"needed_ties={needed} numerical_ties={ties}",
            flush=True,
        )
        if needed <= 0 or len(ties) < needed:
            results.append({
                "size": row["size"],
                "canon": row["canon"],
                "skipped": True,
                "reason": "not a low-active square-subset target",
            })
            continue
        subsets = list(itertools.combinations(ties, needed))
        if args.max_subsets:
            subsets = subsets[: args.max_subsets]

        case_rows = []
        for subset in subsets:
            nonstrict = tuple(i for i in ties if i not in subset)
            label = ",".join(tlabel(i) for i in subset)
            extra = ",".join(tlabel(i) for i in nonstrict) or "-"
            print(f"  subset [{label}], nonstrict extras [{extra}]", flush=True)
            refine_out, refine_result = maybe_run_refine(row, subset, args)
            cert_out = None
            cert_result = {"returncode": None}
            status = {"exists": False, "certified": False}
            if refine_result.get("returncode") == 0 and not args.dry_run:
                cert_out, cert_result = maybe_run_certify(row, subset, nonstrict, refine_out, args)
                if cert_result.get("returncode") == 0:
                    status = certificate_status(cert_out)
            subset_row = {
                "ties": list(subset),
                "tie_triples": [list(TRIPLES[i]) for i in subset],
                "nonstrict_outside_indices": list(nonstrict),
                "nonstrict_outside_triples": [list(TRIPLES[i]) for i in nonstrict],
                "refine_out": str(refine_out),
                "certify_out": str(cert_out) if cert_out else None,
                "refine_result": refine_result,
                "certify_result": cert_result,
                "status": status,
            }
            case_rows.append(subset_row)
            print(
                "    certified="
                f"{status.get('certified')} margin=[{status.get('margin_lower')}, "
                f"{status.get('margin_upper')}]",
                flush=True,
            )
            if status.get("certified") and not args.all_subsets_after_success:
                break

        results.append({
            "size": row["size"],
            "canon": row["canon"],
            "file": row["file"],
            "active": row["active"],
            "active_triples": [list(TRIPLES[i]) for i in row["active"]],
            "numerical_ties": ties,
            "numerical_tie_triples": [list(TRIPLES[i]) for i in ties],
            "needed_ties": needed,
            "subsets": case_rows,
            "any_certified": any(r["status"].get("certified") for r in case_rows),
        })

    payload = {
        "inputs_glob": args.glob,
        "sizes": sorted(sizes),
        "canons": sorted(canons),
        "results": results,
    }
    Path(args.summary_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_out).write_text(json.dumps(payload, indent=1) + "\n")
    total = len(results)
    certified = sum(1 for row in results if row.get("any_certified"))
    print("=" * 78)
    print(f"certified cases: {certified}/{total}")
    print(f"wrote {args.summary_out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
