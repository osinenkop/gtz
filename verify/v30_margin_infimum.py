#!/usr/bin/env python3
"""
v30_margin_infimum.py -- test whether the "wide margin" of §3i.3/§3i.4 is real, by
DIRECTLY MINIMISING the overshoot on the active-equality locus.

WHY THIS RUN EXISTS.  §3i.3/§3i.4 claimed a uniform margin: on the active-equality
locus of a low-active orbit, some outside triple overshoots 1/6 by at least
~4.5e-2.  That claim came from 30-40 random starts per orbit.  v29, sampling 200
starts on the |A|=6 orbit, found a locus point whose best overshoot is only
2.2e-4 -- a factor 264 smaller.  A minimum that collapses when you sample more is
the signature of an INFIMUM APPROACHING ZERO, not of a uniform margin.

THE STAKES.  Define, on the locus V_A = {lambda_min(P_TT) = 1/6 for all T in A},
    g(P) := max_{T not in A} ( lambda_min(P_TT) - 1/6 ).
  * inf g > 0  => the active set always grows with a margin; the low-active strata
                  are empty and §3i.3's mechanism stands (with a smaller constant).
  * inf g = 0, attained  => there is a point with actual active set EXACTLY A.
                  Since |A| <= 9 < 10, that point is NON-SHARP, so it may sit on a
                  positive-dimensional piece of the extremal set.  THE FINITENESS
                  ROUTE OF proofs/finiteness-of-extremal-set.md WOULD BREAK.
  * inf g = 0, not attained => the strata are still empty but no margin exists, so
                  every certificate strategy premised on a margin (§3i.3-§3i.5) is
                  dead and the exact layer must handle a degenerate boundary.

Note this is NOT a counterexample to GTZ(6,3): such a point still has F = 1/6, not
below.  It attacks FINITENESS, not the hypothesis.  The GTZ tripwire is armed anyway.

METHOD.  Minimise g(P) with a hard penalty keeping P on the locus, from many starts,
and report the achieved (g, equality residual) pairs.  A point with g ~ 1e-12 and
residual ~ 1e-14 would be the decisive object; it is then re-validated by the
independent gap-rule active-set recomputation used everywhere else in this project.
"""
import argparse, itertools, json, os, sys
import numpy as np
from scipy.optimize import minimize

TRIPLES = list(itertools.combinations(range(6), 3))
IDX = {t: i for i, t in enumerate(TRIPLES)}
TH = 1.0 / 6.0


def retract(X):
    U, _, Vt = np.linalg.svd(X, full_matrices=False)
    return U @ Vt


def lam_all(P):
    return np.array([np.linalg.eigvalsh(P[np.ix_(T, T)])[0] for T in TRIPLES])


def parts(v, active, inact):
    A = retract(v.reshape(6, 3))
    P = A @ A.T
    lam = lam_all(P)
    r = float(sum((lam[i] - TH) ** 2 for i in active))
    g = float(max(lam[j] - TH for j in inact))
    return r, g, lam, P


def objective(v, active, inact, w):
    r, g, _, _ = parts(v, active, inact)
    return g + w * r


def gap_active(lam):
    dev = np.abs(lam - TH)
    order = np.argsort(dev)
    sd = dev[order]
    bk, br = 1, 0.0
    for k in range(1, 20):
        ratio = sd[k] / max(sd[k - 1], 1e-16)
        if ratio > br:
            bk, br = k, ratio
    return sorted(int(i) for i in order[:bk]), float(br)


def build_report(a, active, best, k_done, verdict=None, complete=False):
    best_sorted = sorted(best, key=lambda b: b["g"])
    payload = dict(canon=a.canon, size=a.size, active=active, starts=a.starts,
                   completed_starts=k_done, seed=a.seed, weight=a.weight,
                   maxiter=a.maxiter, complete=complete,
                   n_on_locus=len(best_sorted), verdict=verdict)
    if best_sorted:
        gs = np.array([b["g"] for b in best_sorted])
        payload.update(
            g_min=float(gs.min()),
            g_median=float(np.median(gs)),
            best=[{k: v for k, v in b.items() if k != "v"}
                  for b in best_sorted[:40]],
            best_point=best_sorted[0]["v"],
            best_points=[b["v"] for b in best_sorted[:a.keep_points]],
        )
    return payload


def write_report(path, payload):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(payload, f, indent=1)
    os.replace(tmp, path)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--canon", type=int, default=78593)
    ap.add_argument("--size", type=int, default=6)
    ap.add_argument("--starts", type=int, default=150)
    ap.add_argument("--seed", type=int, default=20260807)
    ap.add_argument("--weight", type=float, default=1e6)
    ap.add_argument("--maxiter", type=int, default=8000)
    ap.add_argument("--keep-points", type=int, default=8)
    ap.add_argument("--checkpoint-every", type=int, default=10)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    d = json.load(open("verify/out/v22_low_active.json"))
    orb = next(x for x in d["full_pair_cover"]
               if x["canon"] == a.canon and x["size"] == a.size)
    active = sorted(IDX[tuple(t)] for t in orb["triples"])
    inact = [i for i in range(20) if i not in active]
    os.makedirs("verify/out", exist_ok=True)
    out_path = a.out or f"verify/out/v30_margin_{a.size}_{a.canon}.json"

    print("=" * 78)
    print(f"MARGIN INFIMUM TEST  orbit size {a.size} canon {a.canon}")
    print(f"  active {active}")
    print("  minimising  g = max_{T outside} (lambda_min(P_TT) - 1/6)  ON the locus")
    print("  g -> 0 with small residual  =>  a genuine low-active extremal exists")
    print("                                  => FINITENESS ROUTE BREAKS")
    print("=" * 78, flush=True)

    rng = np.random.default_rng(a.seed)
    best = []
    tripwire = []
    for k in range(a.starts):
        v = rng.standard_normal(18)
        for w in (1e3, 1e5, a.weight):
            for _ in range(3):
                v = minimize(objective, v, args=(active, inact, w),
                             method="Nelder-Mead",
                             options=dict(maxiter=a.maxiter, maxfev=a.maxiter,
                                          xatol=1e-15, fatol=1e-18)).x
        r, g, lam, P = parts(v, active, inact)
        F = float(np.max(lam))
        if F < TH - 1e-6:
            tripwire.append((F, v.tolist()))
        if r < 1e-12:
            aa, gr = gap_active(lam)
            best.append(dict(g=g, resid=r, F=F, actual=aa, actual_size=len(aa),
                             gap_ratio=gr, matches=bool(aa == active),
                             v=v.tolist()))
        if a.checkpoint_every and (k + 1) % a.checkpoint_every == 0:
            write_report(out_path, build_report(a, active, best, k + 1,
                                                verdict="RUNNING",
                                                complete=False))
        if (k + 1) % 25 == 0:
            gs = [b["g"] for b in best]
            print(f"  [{k+1}/{a.starts}] on-locus {len(best)}  "
                  f"min g so far {min(gs) if gs else float('nan'):.4e}", flush=True)

    if tripwire:
        print("\n*** GTZ TRIPWIRE: F < 1/6 - 1e-6 ***")
        json.dump(tripwire, open("verify/out/v30_TRIPWIRE.json", "w"), indent=1)
        sys.exit(2)

    if not best:
        print("no on-locus points; increase --starts")
        sys.exit(1)

    best.sort(key=lambda b: b["g"])
    gs = np.array([b["g"] for b in best])
    print(f"\non-locus points: {len(best)}")
    print(f"  g:  min={gs.min():.6e}  p5={np.percentile(gs,5):.6e}  "
          f"median={np.median(gs):.6e}")
    print("\nten smallest g, with the RE-VALIDATED actual active set:")
    for b in best[:10]:
        print(f"   g={b['g']:.4e}  resid={b['resid']:.2e}  F={b['F']:.12f}  "
              f"|A|actual={b['actual_size']}  matches={b['matches']}  "
              f"gapratio={b['gap_ratio']:.1e}")

    tiny = [b for b in best if b["g"] < 1e-8]
    matched = [b for b in tiny if b["matches"]]
    print(f"\npoints with g < 1e-8:                 {len(tiny)}")
    print(f"   ... whose ACTUAL active set is A:  {len(matched)}")
    verdict = ("BREAKS-FINITENESS" if matched else
               "INFIMUM-NEAR-ZERO" if gs.min() < 1e-5 else "MARGIN-HOLDS")
    print(f"\nVERDICT: {verdict}")
    if verdict == "BREAKS-FINITENESS":
        print("  A genuine low-active extremal was found and re-validated.")
        print("  |A| <= 9 < 10 => non-sharp => the compact-discrete finiteness")
        print("  argument no longer applies.  This must be certified exactly.")
    elif verdict == "INFIMUM-NEAR-ZERO":
        print(f"  No exact zero, but g gets to {gs.min():.2e}: the '>= 4.5e-2 margin'")
        print("  of proofs/finiteness-of-extremal-set.md 3i.3/3i.4 is REFUTED.")
        print("  The strata may still be empty, but NOT with a uniform margin, so")
        print("  every margin-based certificate plan needs revising.")
    else:
        print(f"  g stays above {gs.min():.2e}; a margin survives this harder test.")

    write_report(out_path, build_report(a, active, best, a.starts,
                                        verdict=verdict, complete=True))
    print(f"\nwrote {out_path}")
    print("=" * 78)
