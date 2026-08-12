#!/usr/bin/env python3
"""Screen determinant active subsets by modular linear sections."""
from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gtz63_semialgebraic import (  # noqa: E402
    determinant_text_system,
    indices_from_mask,
    known_active_indices,
    parse_indices,
    text_system_summary,
)
from probe_determinant_ideals import (  # noqa: E402
    random_linear_forms,
    run_singular,
    write_probe_script,
)


def resolve_pool(args) -> tuple[int, ...]:
    sources = [
        args.active_indices is not None,
        args.mask is not None,
        args.known_label is not None,
    ]
    if sum(sources) != 1:
        raise SystemExit("provide exactly one of --active-indices, --mask, --known-label")
    if args.active_indices is not None:
        return parse_indices(args.active_indices)
    if args.mask is not None:
        return indices_from_mask(int(args.mask, 0))
    return known_active_indices(args.known_label)


def choose_subsets(pool: tuple[int, ...], subset_size: int, max_subsets: int, seed: int, randomize: bool):
    subsets = list(itertools.combinations(pool, subset_size))
    if randomize:
        rng = random.Random(seed)
        rng.shuffle(subsets)
    if max_subsets > 0:
        subsets = subsets[:max_subsets]
    return subsets


def top_dimension_from_section(section_dim: int | None, sections: int) -> int | None:
    if section_dim is None or section_dim < 0:
        return None
    return section_dim + sections


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--active-indices")
    parser.add_argument("--mask")
    parser.add_argument("--known-label")
    parser.add_argument("--subset-size", type=int, required=True)
    parser.add_argument("--invert-d", action="store_true")
    parser.add_argument("--max-subsets", type=int, default=25, help="0 means all")
    parser.add_argument("--randomize", action="store_true")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--method", default="slimgb", choices=["std", "slimgb", "facstd"])
    parser.add_argument("--characteristic", type=int, default=32003)
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument(
        "--linear-sections",
        type=int,
        default=None,
        help="default is 9 - subset_size, testing expected codimension subset_size",
    )
    parser.add_argument("--out-prefix", default="code/sage/out/screen_det_subsets")
    args = parser.parse_args()

    pool = resolve_pool(args)
    section_count = 9 - args.subset_size if args.linear_sections is None else args.linear_sections
    if section_count < 0:
        raise SystemExit("linear section count must be nonnegative")

    subsets = choose_subsets(
        pool,
        subset_size=args.subset_size,
        max_subsets=args.max_subsets,
        seed=args.seed,
        randomize=args.randomize,
    )
    out_prefix = Path(args.out_prefix)
    script_dir = out_prefix.with_suffix("")
    script_dir.mkdir(parents=True, exist_ok=True)
    records = []

    print(
        "subset,section_dim,inferred_top_dim,status,elapsed_seconds",
        flush=True,
    )
    for ordinal, subset in enumerate(subsets):
        system = determinant_text_system(subset, invert_d=args.invert_d)
        linear_forms = random_linear_forms(
            system.variables,
            count=section_count,
            seed=args.seed + ordinal,
            characteristic=args.characteristic,
            coefficient_bound=17,
        )
        subset_name = "_".join(str(i) for i in subset)
        script_path = script_dir / f"subset_{subset_name}_{args.method}.sing"
        write_probe_script(
            system,
            script_path,
            method=args.method,
            order="degrevlex",
            characteristic=args.characteristic,
            linear_forms=linear_forms,
        )
        run = run_singular(script_path, "Singular", args.timeout)
        parsed = run.get("parsed", {})
        section_dim = parsed.get("dimension") if isinstance(parsed, dict) else None
        section_dim = section_dim if isinstance(section_dim, int) else None
        inferred = top_dimension_from_section(section_dim, section_count)
        status = "timeout" if run["timed_out"] else f"exit {run['exit_code']}"
        record = {
            "subset": list(subset),
            "summary": text_system_summary(system),
            "linear_sections": section_count,
            "section_dim": section_dim,
            "inferred_top_dim": inferred,
            "status": status,
            "elapsed_seconds": run["elapsed_seconds"],
            "script": str(script_path),
            "parsed": parsed,
        }
        records.append(record)
        print(
            f"{list(subset)},{section_dim},{inferred},{status},{run['elapsed_seconds']:.3f}",
            flush=True,
        )

    result = {
        "pool": list(pool),
        "subset_size": args.subset_size,
        "n_tested": len(records),
        "method": args.method,
        "characteristic": args.characteristic,
        "invert_d": args.invert_d,
        "timeout_seconds": args.timeout,
        "linear_sections": section_count,
        "records": records,
    }
    result_path = out_prefix.with_suffix(".json")
    result_path.parent.mkdir(parents=True, exist_ok=True)
    result_path.write_text(json.dumps(result, indent=1) + "\n")
    print(f"wrote {result_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
