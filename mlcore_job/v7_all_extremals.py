#!/usr/bin/env python3
"""
v7_all_extremals.py -- enumerate the Nesterenko scaled-star extremal family for
(6,3) and run the EXACT sharpness certificate of v6 at every one of them.

BACKGROUND / WHAT WAS REVERSE-ENGINEERED.  boundary-obstruction.md cites "the nine
Nesterenko scaled-star matrices" but names only four graphs and never states the
weight rule.  Two facts were recovered here and verified exactly at the known case:

  * The construction is the weighted graphic matroid / resistor-network one:
    for a two-terminal series-parallel (TTSP) graph G on 6 edges with edge
    weights w, put Mtilde = W^{1/2} B_red^T and P = Mtilde (Mtilde^T Mtilde)^{-1} Mtilde^T.
  * The leverage of edge e is then exactly  ell_e = w_e * R_eff(e),  the classical
    identity (leverage = weight x effective resistance).  Verified at
    3xLR+path(1,1,1): R_eff(0,1) = 5/18 and (9/5) R_eff(0,2) = 13/18, matching
    Lemma 2.1's leverages exactly.

So an "extremal" is a TTSP graph on 6 edges together with a weight vector chosen
to make F = 1/6 exactly.  Rather than guess the corpus's nine, we ENUMERATE TTSP
graphs on 6 edges, and for each we SOLVE exactly for weights making the
non-singular triples all attain lambda_min = 1/6, then certify sharpness.

METHOD PER CANDIDATE (all exact, sympy over Q or Q(sqrt5)):
  1. build B_red for the graph; symbolic weights on each orbit of the graph's
     automorphism group (this keeps the solve small and matches the corpus's
     two-value weight patterns);
  2. impose  det(P_TT - I/6) = 0  for the non-singular triples and solve;
  3. keep real positive solutions with F(P) = 1/6 exactly and all leverages > 1/6
     (core region);
  4. run the v6 certificate: active set, simplicity of lambda_min, tangent basis,
     10-or-more x 9 functional matrix L, rank, KKT nullspace, STRICT positivity of
     all multipliers  =>  sharp minimum, kappa > 0.

Nothing is floating-point load-bearing.  Deterministic.
"""
import itertools, sys, json, os
import sympy as sp

SIX = sp.Rational(1, 6)
TRIPLES = list(itertools.combinations(range(6), 3))


def R(e):
    return sp.radsimp(sp.simplify(sp.expand(e)))


# --------------------------------------------------------------------------
# TTSP graph construction.  Represent a TTSP graph as a nested expression:
#   ('e',)                -> single edge
#   ('S', g1, g2, ...)    -> series composition
#   ('P', g1, g2, ...)    -> parallel composition
# Build the edge list (as vertex pairs) with fresh internal vertices.
# --------------------------------------------------------------------------
class Builder:
    def __init__(self):
        self.nv = 2                      # 0 = source, 1 = sink
        self.edges = []

    def new_vertex(self):
        v = self.nv
        self.nv += 1
        return v

    def build(self, g, s, t):
        if g == ('e',):
            self.edges.append((s, t))
            return
        kind, *kids = g
        if kind == 'P':
            for k in kids:
                self.build(k, s, t)
        elif kind == 'S':
            prev = s
            for i, k in enumerate(kids):
                nxt = t if i == len(kids) - 1 else self.new_vertex()
                self.build(k, prev, nxt)
                prev = nxt
        else:
            raise ValueError(kind)


E = ('e',)


def series(*ks):
    return ('S',) + tuple(ks)


def parallel(*ks):
    return ('P',) + tuple(ks)


def path(n):
    return E if n == 1 else series(*([E] * n))


def enumerate_ttsp(n_edges):
    """All TTSP expression trees with exactly n_edges edges, up to the obvious
    flattening (no S inside S, no P inside P), children sorted for canonicity."""
    memo = {}

    def gen(n, allow):
        key = (n, allow)
        if key in memo:
            return memo[key]
        out = []
        if n == 1:
            out.append(E)
        else:
            for kind in ('S', 'P'):
                if kind == allow:
                    continue
                # compositions of n into >=2 parts
                for parts in compositions(n):
                    if len(parts) < 2:
                        continue
                    child_opts = [gen(p, kind) for p in parts]
                    for combo in itertools.product(*child_opts):
                        out.append((kind,) + tuple(sorted(combo, key=repr)))
        res = list(dict.fromkeys(out))
        memo[key] = res
        return res

    def compositions(n):
        if n == 0:
            yield ()
            return
        for first in range(1, n + 1):
            for rest in compositions(n - first):
                yield (first,) + rest

    return list(dict.fromkeys(gen(n_edges, None)))


def name_of(g):
    if g == E:
        return "e"
    kind, *kids = g
    return ("S(" if kind == 'S' else "P(") + ",".join(name_of(k) for k in kids) + ")"


def reduced_incidence(edges, nv):
    """B_red: rows = vertices 1..nv-1 (drop vertex 0), cols = edges."""
    B = sp.zeros(nv - 1, len(edges))
    for j, (a, b) in enumerate(edges):
        if a != 0:
            B[a - 1, j] += 1
        if b != 0:
            B[b - 1, j] -= 1
    return B


def edge_orbits(edges, nv):
    """Crude but effective orbit detection: group edges by a structural signature
    (multiset of degrees of endpoints + whether both endpoints are terminals)."""
    deg = {}
    for a, b in edges:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    sig = {}
    for j, (a, b) in enumerate(edges):
        key = (tuple(sorted((deg[a], deg[b]))),
               tuple(sorted((a in (0, 1), b in (0, 1)))))
        sig.setdefault(key, []).append(j)
    return list(sig.values())


# --------------------------------------------------------------------------
def solve_weights(eqs, free):
    """Solve the exact weight equations robustly.

    sp.solve(..., dict=True) on several symbols raises "can only solve for one
    symbol at a time" for some of these systems, so try a ladder:
      1. groebner-based nonlinsolve (handles multivariate polynomial systems);
      2. sp.solve on the whole system;
      3. if there is a single free symbol, solve each equation separately and keep
         values that satisfy all of them.
    Returns a list of substitution dicts."""
    if not eqs:
        return [{}]
    out = []
    try:
        sset = sp.nonlinsolve(eqs, free)
        if sset is not sp.EmptySet and not isinstance(sset, sp.ConditionSet):
            for tup in sset:
                if len(tup) == len(free):
                    out.append(dict(zip(free, tup)))
    except Exception:
        pass
    if not out:
        try:
            s = sp.solve(eqs, free, dict=True)
            out.extend(s if isinstance(s, list) else [s])
        except Exception:
            pass
    if not out and len(free) == 1:
        x = free[0]
        cands = set()
        for e in eqs:
            try:
                for r in sp.solve(e, x):
                    cands.add(sp.nsimplify(r))
            except Exception:
                pass
        for r in cands:
            if all(sp.simplify(e.subs({x: r})) == 0 for e in eqs):
                out.append({x: r})
    # de-duplicate
    seen, uniq = set(), []
    for d in out:
        key = tuple(sorted((str(k), str(sp.nsimplify(v))) for k, v in d.items()))
        if key not in seen:
            seen.add(key)
            uniq.append(d)
    return uniq


def projector_from(B, wsym):
    W = sp.diag(*[sp.sqrt(x) for x in wsym])
    Mt = W * B.T
    if Mt.rank() != 3:
        return None
    return sp.simplify(Mt * (Mt.T * Mt).inv() * Mt.T)


def analyze(P, label, log):
    """Run the v6 sharpness certificate on an exact projector P. Returns a dict."""
    out = dict(label=label)
    P = P.applyfunc(R)
    if sp.simplify(sp.expand(P * P - P)) != sp.zeros(6, 6) or R(sp.trace(P)) != 3:
        out["ok"] = False
        out["why"] = "not a rank-3 projector"
        return out

    lev = [R(P[i, i]) for i in range(6)]
    out["leverages"] = [str(x) for x in lev]
    out["core"] = all(sp.simplify(x - SIX) > 0 for x in lev)

    active, singular, below = [], [], []
    for T in TRIPLES:
        M = sp.Matrix(3, 3, lambda a, b: P[T[a], T[b]]).applyfunc(R)
        if R(M.det()) == 0:
            singular.append(T)
            continue
        D = (M - SIX * sp.eye(3)).applyfunc(R)
        psd = True
        for k in (1, 2, 3):
            for idx in itertools.combinations(range(3), k):
                if not (sp.simplify(D[list(idx), list(idx)].det()) >= 0):
                    psd = False
        if psd and R(D.det()) == 0:
            active.append(T)
        elif not psd:
            below.append(T)
    out.update(n_active=len(active), n_singular=len(singular), n_below=len(below))
    # F = 1/6 exactly requires: at least one active, and nothing strictly between
    # 0 and 1/6 among the rest.
    strictly_between = []
    for T in singular + below:
        M = sp.Matrix(3, 3, lambda a, b: P[T[a], T[b]]).applyfunc(R)
        lm = min(sp.nsimplify(r) for r in M.eigenvals().keys())
        if sp.N(lm, 30) > sp.Float(10) ** -25 and sp.simplify(lm - SIX) < 0:
            strictly_between.append(T)
    out["F_is_1_6"] = (len(active) > 0 and not strictly_between)
    if not out["F_is_1_6"]:
        out["ok"] = False
        out["why"] = f"F != 1/6 ({len(strictly_between)} triples in (0,1/6))"
        return out

    # eigenvector simplicity
    Ev, simple = {}, True
    for T in active:
        M = sp.Matrix(3, 3, lambda a, b: P[T[a], T[b]]).applyfunc(R)
        ker = (M - SIX * sp.eye(3)).applyfunc(R).nullspace()
        if len(ker) != 1:
            simple = False
            break
        v = ker[0].applyfunc(R)
        Ev[T] = (v / sp.sqrt(R((v.T * v)[0]))).applyfunc(R)
    out["simple_lmin"] = simple
    if not simple:
        out["ok"] = False
        out["why"] = "lambda_min not simple at some active triple (needs SDP form)"
        return out

    # tangent space
    def gs(vs):
        o = []
        for v in vs:
            u = v
            for q in o:
                u = (u - R((q.T * v)[0]) * q).applyfunc(R)
            o.append((u / sp.sqrt(R((u.T * u)[0]))).applyfunc(R))
        return o

    U = gs([v.applyfunc(R) for v in P.columnspace()][:3])
    N = gs([v.applyfunc(R) for v in P.nullspace()][:3])
    TAN = [(U[a] * N[b].T + N[b] * U[a].T).applyfunc(R) for a in range(3) for b in range(3)]
    if len(TAN) != 9:
        out["ok"] = False
        out["why"] = "tangent space not 9-dimensional"
        return out

    L = sp.Matrix([[R((Ev[T].T * sp.Matrix(3, 3, lambda a, b: B[T[a], T[b]]) * Ev[T])[0])
                    for B in TAN] for T in active])
    rk = L.rank()
    ns = L.T.nullspace()
    out.update(L_shape=list(L.shape), rank_L=rk, dim_kkt=len(ns))

    if rk != 9:
        out["ok"] = False
        out["why"] = f"rank(L) = {rk} != 9: active gradients do not span"
        return out
    if len(ns) == 0:
        out["ok"] = False
        out["why"] = "no KKT multiplier vector (0 not in span) -> NOT a critical point"
        return out

    # need SOME nonneg combination that is strictly positive; with dim>1 search
    # the nullspace for a strictly-positive vector via LP-free reasoning:
    strict, lam_used = False, None
    if len(ns) == 1:
        ent = [R(ns[0][i]) for i in range(ns[0].rows)]
        sg = [sp.sign(sp.nsimplify(e)) for e in ent]
        if all(s == 1 for s in sg) or all(s == -1 for s in sg):
            strict = True
            lam_used = ent if sg[0] == 1 else [R(-e) for e in ent]
    else:
        # dim > 1: look for a strictly positive vector in the KKT nullspace.
        # NOTE: sp.solve on a list of strict inequalities in several unknowns is
        # unsupported (raises "can only solve for one symbol at a time") -- do not
        # call it here.  Try the basis vectors and a few nonneg combinations; that
        # suffices to CERTIFY when it succeeds, and we report honestly when it
        # does not (absence of a witness is not a proof of nonexistence).
        m = ns[0].rows
        cands = list(ns) + [sum(ns[1:], ns[0])]
        for k in range(len(ns)):
            cands.append(sum((ns[i] for i in range(len(ns)) if i != k), ns[k]))
        for cand in cands:
            for sgn in (1, -1):
                ent = [R(sgn * cand[i]) for i in range(m)]
                if all(sp.sign(sp.nsimplify(e)) == 1 for e in ent):
                    strict, lam_used = True, ent
                    break
            if strict:
                break
        if not strict:
            out["kkt_note"] = ("dim_kkt>1 and no strictly-positive witness found "
                               "among tried combinations; inconclusive, not negative")
    out["strict_multipliers"] = strict
    if lam_used is not None:
        tot = R(sum(lam_used))
        out["multipliers"] = [str(R(e / tot)) for e in lam_used]

    if not strict:
        out["ok"] = False
        out["why"] = "KKT multipliers not all strictly positive -> sharpness NOT certified"
        return out

    out["ok"] = True
    out["why"] = "SHARP: rank(L)=9 and all multipliers strictly positive => kappa>0"
    return out


# --------------------------------------------------------------------------
if __name__ == "__main__":
    log = []
    graphs = enumerate_ttsp(6)
    print(f"TTSP expression trees on 6 edges: {len(graphs)}", flush=True)

    found = []
    for gi, g in enumerate(graphs):
        b = Builder()
        b.build(g, 0, 1)
        if b.nv - 1 < 3:
            continue
        B = reduced_incidence(b.edges, b.nv)
        if B.rank() != 3:
            continue                      # need exactly a rank-3 (4-vertex) picture
        if B.rows != 3:
            continue
        orbits = edge_orbits(b.edges, b.nv)
        nm = name_of(g)
        # symbolic weight per orbit, normalized so the first orbit has weight 1
        ws = sp.symbols(f"w0:{len(orbits)}", positive=True)
        wvec = [None] * 6
        for oi, orb in enumerate(orbits):
            for j in orb:
                wvec[j] = sp.Integer(1) if oi == 0 else ws[oi]
        free = [ws[oi] for oi in range(1, len(orbits))]
        P = projector_from(B, wvec)
        if P is None:
            continue
        # impose det(P_TT - I/6) = 0 on all non-singular triples
        eqs = set()
        for T in TRIPLES:
            M = sp.Matrix(3, 3, lambda a, bb: P[T[a], T[bb]])
            d = sp.simplify(M.det())
            if d == 0:
                continue
            e = sp.simplify((M - SIX * sp.eye(3)).det())
            e = sp.simplify(sp.numer(sp.together(e)))
            if e != 0 and e.free_symbols:
                eqs.add(sp.factor(e))
        if not free:
            sols = [{}]
        else:
            sols = solve_weights(list(eqs), free)
        print(f"\n[{gi}] {nm}   edges={b.edges} orbits={[len(o) for o in orbits]} "
              f"free={len(free)} eqs={len(eqs)} sols={len(sols)}", flush=True)
        for sol in sols:
            if any((not v.is_real) or v.is_negative or v == 0
                   for v in sol.values() if hasattr(v, "is_real")):
                continue
            try:
                wv = [sp.nsimplify(sp.simplify(x.subs(sol) if hasattr(x, "subs") else x))
                      for x in wvec]
                if any((not sp.simplify(x).is_positive) for x in wv):
                    continue
                Pn = projector_from(B, wv)
                if Pn is None:
                    continue
                res = analyze(Pn, f"{nm} w={[str(x) for x in wv]}", log)
                print(f"     weights {[str(x) for x in wv]}", flush=True)
                print(f"     -> active={res.get('n_active')} sing={res.get('n_singular')} "
                      f"core={res.get('core')} F=1/6:{res.get('F_is_1_6')} "
                      f"rank={res.get('rank_L')} kkt={res.get('dim_kkt')} "
                      f"STRICT={res.get('strict_multipliers')}", flush=True)
                print(f"        {res.get('why')}", flush=True)
                if res.get("F_is_1_6"):
                    res["graph"] = nm
                    res["edges"] = b.edges
                    res["weights"] = [str(x) for x in wv]
                    found.append(res)
            except Exception as ex:
                print(f"     (skipped: {type(ex).__name__}: {ex})", flush=True)

    print("\n" + "=" * 78)
    print(f"EXTREMAL CANDIDATES WITH F = 1/6 EXACTLY: {len(found)}")
    sharp = [f for f in found if f.get("ok")]
    print(f"  of which SHARP (certified kappa > 0): {len(sharp)}")
    print(f"  not certified:                        {len(found) - len(sharp)}")
    for f in found:
        tag = "SHARP" if f.get("ok") else "NOT CERTIFIED"
        print(f"\n  [{tag}] {f['graph']}  w={f['weights']}")
        print(f"     leverages={f['leverages']}")
        print(f"     active={f['n_active']} singular={f['n_singular']} "
              f"rank(L)={f.get('rank_L')} dim_kkt={f.get('dim_kkt')}")
        if f.get("multipliers"):
            print(f"     multipliers={f['multipliers']}")
        if not f.get("ok"):
            print(f"     reason: {f['why']}")

    os.makedirs("verify/out", exist_ok=True)
    with open("verify/out/v7_extremals.json", "w") as fh:
        json.dump(found, fh, indent=1)
    print("\nwrote verify/out/v7_extremals.json")
    print("=" * 78)
