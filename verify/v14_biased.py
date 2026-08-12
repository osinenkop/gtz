#!/usr/bin/env python3
"""
v14_biased.py -- COMPLEMENTARY search for extremals using BIASED starts, to probe
regions uniform-random descent under-samples.

Rationale: v12 starts from N(0,1) Stiefel points, which concentrate near "generic"
leverage profiles.  Extremals with lopsided leverages, or with many active triples,
may sit in small basins that uniform starts rarely hit.  Here we seed from:
  (a) perturbations of each KNOWN extremal (probing for nearby distinct extremals,
      which is also a direct test of the certified radius r0);
  (b) leverage-targeted starts: build A with prescribed row norms across a grid of
      leverage patterns (k rows heavy, 6-k light), then descend;
  (c) high-symmetry starts: orbits of small groups, ETFs, and the icosahedral
      configuration mentioned in the corpus.

Reports any leverage pattern not in the known list.  NUMERICALLY SUPPORTED only.
"""
import itertools, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize

MASTER = 20260801
TRIPLES = list(itertools.combinations(range(6), 3))
TH = 1/6
CX = 1e-6
HIT = 1e-9
KNOWN_LEV = [sorted([5/18]*3+[13/18]*3), sorted([11/18]*4+[5/18]*2),
             sorted([7/18]*4+[13/18]*2),  sorted([5/14]*3+[9/14]*3)]
def is_known(key, tol=1e-6):
    k = np.sort(np.array(key))
    return any(np.max(np.abs(k - np.array(v))) < tol for v in KNOWN_LEV)

def retract(X):
    U,_,Vt=np.linalg.svd(X,full_matrices=False); return U@Vt
def lams(A):
    P=A@A.T
    return np.array([np.linalg.eigvalsh(P[np.ix_(T,T)])[0] for T in TRIPLES]),P
def F(v): return float(np.max(lams(retract(v.reshape(6,3)))[0]))
def softF(v,b):
    l=lams(retract(v.reshape(6,3)))[0]; m=float(np.max(l))
    return m+float(np.log(np.sum(np.exp(b*(l-m))))/b)

def run(args):
    kind,v0 = args
    v=np.array(v0,dtype=float)
    for b in (30.,200.,1500.):
        v=minimize(softF,v,args=(b,),method='Nelder-Mead',
                   options=dict(maxiter=5000,maxfev=5000,xatol=1e-12,fatol=1e-15)).x
    for _ in range(8):
        v=minimize(F,v,method='Nelder-Mead',
                   options=dict(maxiter=20000,maxfev=20000,xatol=1e-16,fatol=1e-18)).x
    f=F(v)
    if f-TH < -CX: return dict(kind='VIOLATION',src=kind,F=f,v=v.tolist())
    if abs(f-TH)>HIT: return dict(kind='miss',src=kind,F=f)
    A=retract(v.reshape(6,3)); lev=np.sort(np.diag(A@A.T))
    # NOTE: cluster at 1e-6, NOT 1e-9.  Biased starts converge less precisely than
    # uniform ones, and a 1e-9 key split ONE pattern into ~100 spurious "new"
    # patterns in an earlier run.  1e-6 is still far tighter than the gaps between
    # genuinely distinct leverage patterns (the closest pair differs by ~0.03).
    return dict(kind='hit',src=kind,F=f,v=v.tolist(),
                key=tuple(np.round(lev,6)),lev=lev.tolist())

def make_starts():
    rng=np.random.default_rng(MASTER); S=[]
    # (a) perturbations of known extremals (also tests the r0 radius)
    for f in sorted(os.listdir('verify/data')) if os.path.isdir('verify/data') else []:
        if not f.endswith('.npy'): continue
        P=np.load(f'verify/data/{f}')
        w,U=np.linalg.eigh(P); A=U[:,np.argsort(-w)[:3]]
        for scale in (0.02,0.08,0.25,0.6):
            for _ in range(30):
                d=rng.standard_normal((6,3)); d/=np.linalg.norm(d)
                S.append((f'perturb:{f}:{scale}',(A+scale*d).reshape(-1).tolist()))
    # (b) leverage-targeted starts
    for k in range(1,6):
        for heavy in (0.55,0.62,0.68,0.74,0.85):
            light=(3-k*heavy)/(6-k)
            if not (0.01<light<0.99): continue
            for _ in range(40):
                A=rng.standard_normal((6,3))
                nrm=np.linalg.norm(A,axis=1,keepdims=True)
                tgt=np.array([heavy]*k+[light]*(6-k))[:,None]
                A=A/nrm*np.sqrt(np.abs(tgt))
                S.append((f'lev:k{k}:h{heavy}',A.reshape(-1).tolist()))
    # (c) symmetric starts: icosahedral 6 lines, and small-group orbits
    phi=(1+5**0.5)/2
    ico=np.array([[0,1,phi],[0,1,-phi],[1,phi,0],[1,-phi,0],[phi,0,1],[-phi,0,1]],float)
    ico/=np.linalg.norm(ico,axis=1,keepdims=True)
    S.append(('ico',retract(ico).reshape(-1).tolist()))
    for _ in range(60):
        Om=np.linalg.qr(rng.standard_normal((3,3)))[0]
        S.append(('ico-rot',retract(ico@Om).reshape(-1).tolist()))
    return S

if __name__=='__main__':
    S=make_starts()
    ncpu=int(os.environ.get('GTZ_CPUS', max(1,(os.cpu_count() or 4)-2)))
    print(f"biased starts: {len(S)}   cores={ncpu}",flush=True)
    from collections import Counter
    print("by source:",dict(Counter(k.split(':')[0] for k,_ in S)),flush=True)
    with Pool(ncpu) as pool:
        res=pool.map(run,S,chunksize=2)
    viol=[r for r in res if r['kind']=='VIOLATION']
    if viol:
        print(f"\n*** {len(viol)} VIOLATION(S) -- HALT ***")
        json.dump(viol,open('verify/out/v14_VIOLATIONS.json','w'),indent=1); sys.exit(2)
    hits=[r for r in res if r['kind']=='hit']
    print(f"\nhits: {len(hits)} / {len(S)}")
    pats={}
    for r in hits: pats.setdefault(r['key'],[]).append(r)
    print(f"distinct patterns: {len(pats)}")
    new=[]
    for k,mem in sorted(pats.items(),key=lambda kv:-len(kv[1])):
        tag='KNOWN' if is_known(k) else '*** NEW ***'
        srcs=dict(__import__('collections').Counter(m['src'].split(':')[0] for m in mem))
        print(f"  {tag:<12} {np.round(np.array(k),9)}  n={len(mem)} src={srcs}")
        if not is_known(k):
            new.append(dict(key=list(k),n=len(mem),rep=mem[0]['v']))
            A=retract(np.array(mem[0]['v']).reshape(6,3))
            np.save(f"verify/data/v14_new_{k[0]:.6f}_{k[-1]:.6f}.npy",A@A.T)
    json.dump(dict(n_starts=len(S),n_hits=len(hits),n_patterns=len(pats),
                   new_patterns=new),open('verify/out/v14_biased.json','w'),indent=1)
    print(f"\nNEW patterns: {len(new)}")
    if new: print("*** certify each exactly with v11 machinery ***")
    else: print("No new pattern from biased starts either.")
    print("no counterexample throughout.")
