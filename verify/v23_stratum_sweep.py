#!/usr/bin/env python3
"""
v23_stratum_sweep.py -- batch driver that attacks the 124 low-active orbits from
v22 with the Sage/Singular ideal machinery in code/sage/.

GOAL AND CORRECT SUCCESS CRITERION.
  v22 reduced "does any extremal have |A| <= 9?" to 124 explicit S_6-orbits.
  Each orbit A defines a stratum
      E_A = { P : P^2=P, tr P=3,
              lambda_min(P_TT) = 1/6   for T in A,
              lambda_min(P_TT) < 1/6   for T not in A }.
  Finiteness (via the simple-active route) needs every one of these to be EMPTY.

  CRUCIAL: dimension 0 is NOT the success criterion here, and neither is
  "dimension > 0 means trouble".  What we need is EMPTINESS of the SEMIALGEBRAIC
  set, and the determinant ideal is only the EQUALITY part of it.  Concretely:

    * dim = -1 (empty ideal)              => stratum empty  => ORBIT KILLED, exactly.
    * dim >= 0 but no real point satisfies the PSD/inactive inequalities
                                          => stratum empty, but needs the
                                             inequality layer (lex slice / CAD /
                                             real root isolation) to certify.
    * dim >= 0 with a feasible real point => a genuine low-active extremal, which
                                             would BREAK the finiteness route and
                                             is the single most important thing
                                             this sweep could find.

  So this driver reports three buckets and never conflates them.  The cheap
  modular dimension pass is a SCREEN that finds the dim = -1 kills for free; the
  survivors are queued for the expensive exact inequality layer.

  Note the direction of the inference: over a prime field, dim = -1 does NOT by
  itself prove emptiness over Q (bad primes exist), so a modular -1 is recorded as
  "likely killed, needs characteristic-0 confirmation", and only a
  characteristic-0 dim = -1 is recorded as a proof.  This asymmetry is why the
  driver runs the modular pass first and re-runs survivors over QQ.

USAGE
  .venv/bin/python -u verify/v23_stratum_sweep.py --sizes 6,7 --characteristic 32003
  .venv/bin/python -u verify/v23_stratum_sweep.py --sizes 6 --characteristic 0
Deterministic; resumable via the JSON ledger.
"""
import argparse, itertools, json, os, subprocess, sys, time

TRIPLES = list(itertools.combinations(range(6), 3))
IDX = {t: i for i, t in enumerate(TRIPLES)}
SAGE = os.path.expanduser("~/miniforge3/bin/mamba")
LEDGER = "verify/out/v23_stratum_ledger.json"


def load_orbits(sizes):
    d = json.load(open("verify/out/v22_low_active.json"))
    out = []
    for x in d["full_pair_cover"]:
        if x["size"] in sizes:
            idx = sorted(IDX[tuple(t)] for t in x["triples"])
            out.append(dict(size=x["size"], canon=x["canon"],
                            triples=x["triples"], indices=idx))
    out.sort(key=lambda r: (r["size"], r["canon"]))
    return out


def run_probe(indices, char, timeout, tag, extra=()):
    prefix = f"code/sage/out/sweep_{tag}"
    cmd = [SAGE, "run", "-n", "sage", "python",
           "code/sage/probe_determinant_ideals.py",
           "--active-indices", ",".join(map(str, indices)),
           "--invert-d", "--methods", "slimgb",
           "--characteristic", str(char),
           "--timeout", str(timeout),
           "--out-prefix", prefix, *extra]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=timeout + 240)
        out = p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        return dict(status="driver-timeout", elapsed=time.time() - t0)
    el = time.time() - t0
    dim, deg = None, None
    for line in out.splitlines():
        s = line.strip().lower()
        if s.startswith("dimension") and dim is None:
            for tok in line.replace(":", " ").split():
                try:
                    dim = int(tok); break
                except ValueError:
                    pass
        if "degree" in s and deg is None:
            for tok in line.replace(":", " ").split():
                try:
                    deg = int(tok); break
                except ValueError:
                    pass
    js = prefix + ".json"
    if os.path.exists(js):
        try:
            j = json.load(open(js))
            dim = j.get("dimension", dim)
            deg = j.get("degree", deg)
        except Exception:
            pass
    return dict(status="ok" if dim is not None else "no-dimension",
                dimension=dim, degree=deg, elapsed=el, tail=out[-600:])


def classify(dim, char):
    if dim is None:
        return "unresolved"
    if dim < 0:
        return "EMPTY (proved)" if char == 0 else "empty mod p (needs char-0)"
    return f"dim {dim}: needs inequality layer"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", default="6")
    ap.add_argument("--characteristic", type=int, default=32003)
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    sizes = [int(x) for x in a.sizes.split(",")]

    orbits = load_orbits(sizes)
    if a.limit:
        orbits = orbits[:a.limit]
    print("=" * 78)
    print(f"STRATUM SWEEP over low-active orbits, sizes {sizes}")
    print(f"  orbits to process: {len(orbits)}   characteristic: {a.characteristic}")
    print("  success = stratum EMPTY.  dim>=0 is NOT failure: it means the")
    print("  determinant equalities alone do not decide, and the inequality")
    print("  layer (PSD + inactive) must be applied to that orbit.")
    print("=" * 78, flush=True)

    ledger = {}
    if os.path.exists(LEDGER):
        ledger = json.load(open(LEDGER))

    for n, orb in enumerate(orbits, 1):
        key = f"{orb['size']}:{orb['canon']}:{a.characteristic}"
        if key in ledger and ledger[key].get("status") == "ok":
            print(f"[{n}/{len(orbits)}] size {orb['size']} canon {orb['canon']} "
                  f"-- cached: {ledger[key]['verdict']}", flush=True)
            continue
        tag = f"s{orb['size']}_c{orb['canon']}_p{a.characteristic}"
        print(f"[{n}/{len(orbits)}] size {orb['size']} canon {orb['canon']} "
              f"indices {orb['indices']}", flush=True)
        r = run_probe(orb["indices"], a.characteristic, a.timeout, tag)
        r["verdict"] = classify(r.get("dimension"), a.characteristic)
        r["size"] = orb["size"]; r["canon"] = orb["canon"]
        r["triples"] = orb["triples"]
        ledger[key] = r
        print(f"      dim={r.get('dimension')} deg={r.get('degree')} "
              f"({r['elapsed']:.0f}s)  => {r['verdict']}", flush=True)
        os.makedirs("verify/out", exist_ok=True)
        json.dump(ledger, open(LEDGER, "w"), indent=1)

    # summary
    rows = [v for k, v in ledger.items() if v.get("size") in sizes]
    from collections import Counter
    print("\n" + "=" * 78)
    print("SUMMARY")
    for v, c in sorted(Counter(r["verdict"] for r in rows).items()):
        print(f"  {c:>4}  {v}")
    killed = [r for r in rows if r["verdict"].startswith(("EMPTY", "empty"))]
    live = [r for r in rows if r["verdict"].startswith("dim")]
    print(f"\n  killed (empty):                {len(killed)}")
    print(f"  need inequality layer:         {len(live)}")
    print(f"  unresolved/timeout:            {len(rows)-len(killed)-len(live)}")
    if live:
        print("\n  orbits needing the inequality layer (size, canon, dim):")
        for r in sorted(live, key=lambda z: (z["size"], -(z.get("dimension") or 0)))[:20]:
            print(f"    size {r['size']}  canon {r['canon']}  dim {r.get('dimension')}")
    print(f"\n  ledger: {LEDGER}")
    print("=" * 78)
