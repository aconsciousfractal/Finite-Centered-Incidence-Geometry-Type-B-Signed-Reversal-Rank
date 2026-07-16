import sys, os, json
from fractions import Fraction as F
from math import factorial, comb
# public repo layout: the shared exact engine lives in ../scripts/
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
    "..","scripts")))
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass
import p15_s3_exact_engine as s3
# ============================================================================
# G22: ODD PAIR-STANDARD Type-D positivity lam'^D>0 CLOSED via a DIRECT Type-D
# trace identity (mirror of P15-S9B on D_n), replacing the retracted PS-DOM/8^m route.
#
#   D-layer = elements of the (k,k) cell with an EVEN number of negative signs (w in D_n).
#   m_pair^D = sum of their permutation matrices (built by s3.channel_matrices).
#   FACTS (proved; finite-replay-verified n=5,7 exactly, the same standard as P15-S9B n=3,5,7):
#     (T1) trace(m_pair^D) = k|X^D| + NegFixSum^D            [each elt: k pos-fixed + neg-fixed]
#     (T2) trace(m_pair^D) = |X^D| + (m-1)lam'^D + mu^D      [orbital block form; rho in D_n, D_n normal]
#     (T3) mu^D <= |X^D|                                     [m_pair^D = |X^D| x doubly-stochastic; Perron]
#   => (m-1)lam'^D = (k-1)|X^D| + NegFixSum^D - mu^D >= (k-2)|X^D| + NegFixSum^D.
#   Hence for k>=2:  lam'^D>0  (k>=3: (k-2)|X^D|>0;  k=2: NegFixSum^D>0 parity-even witness).
#   k=1:  lam'^D(m,1)>0 by the S9C comparison K_m=(m-1)lam'^B(2m+1,1) > D_m=(m-1)|lam'^chi(m,1)|,
#         using the S9C-certified ratio rho_m=K_{m+1}/K_m>=4m^2 and D_{m+1}/D_m<4m^2 (m>=3).
# ============================================================================

def chi_even(eps): return sum(1 for e in eps if e==-1) % 2 == 0

def replay_typeD(n):
    rho=s3.reversal(n); center=(n-1)//2
    seen={center}; pairs=[]
    for i in range(n):
        if i in seen: continue
        pairs.append((i,rho[i])); seen.add(i); seen.add(rho[i])
    m=(n-1)//2; layers=s3.build_layers(n); rows=[]
    for k in range(1,m+2):
        full=layers.get((k,k),[])
        layer=[(pi,eps) for (pi,eps) in full if chi_even(eps)]
        if not layer: continue
        _,mp=s3.channel_matrices(layer,n)
        size=len(layer); gc=mp[center][center]
        mu=F(n*gc-size,n-1)
        trace=sum(mp[i][i] for i in range(n))
        negfix=sum(s3.neg_fixed_count(pi,eps) for pi,eps in layer)
        i0,ri0=pairs[0]; generic=[j for j in range(n) if j not in (i0,ri0,center)]
        lam=mp[i0][i0]+mp[ri0][i0]-2*mp[generic[0]][i0]
        # T2 support: lam' is the SAME for every non-center pair and every generic column
        lam_set=set()
        for (p0,rp0) in pairs:
            gg=[j for j in range(n) if j not in (p0,rp0,center)]
            for g in gg: lam_set.add(mp[p0][p0]+mp[rp0][p0]-2*mp[g][p0])
        lam_uniform=(len(lam_set)==1 and next(iter(lam_set))==lam)
        # character-split reference
        _,mpB=s3.channel_matrices(full,n)
        mpChi=[[0]*n for _ in range(n)]
        for pi,eps in full:
            w=1
            for e in eps:
                if e==-1: w=-w
            for i in range(n): mpChi[pi[i]][i]+=w
        lamB=mpB[i0][i0]+mpB[ri0][i0]-2*mpB[generic[0]][i0]
        lamChi=mpChi[i0][i0]+mpChi[ri0][i0]-2*mpChi[generic[0]][i0]
        rows.append(dict(n=n,k=k,size=size,mu=mu,lam=lam,negfix=negfix,trace=trace,
            T1=(trace==k*size+negfix),
            T2=(F(trace)==F(size)+F(m-1)*lam+mu),
            T2_lam_uniform_all_pairs=lam_uniform,
            T3_mu_le_size=(mu<=size),
            traceid=(F(m-1)*lam==F((k-1)*size+negfix)-mu),
            csplit_match=(F(lamB+lamChi,2)==lam),
            lamD_pos=(lam>0),
            k_ge2_bound=(F(m-1)*lam >= F(k-2)*size) if k>=2 else None,
            k2_witness=(negfix>0) if k==2 else None,
        ))
    return rows

# ---- k=1 S9C comparison (K_m vs D_m) ----
def L_poly_t(m):
    e=[2*m-3,-4*m*m-8*m+19,24*m*m+18*m-68,-56*m*m-24*m+122,56*m*m+24*m-92,
       -12*m*m-44*m+12,-16*m*m+56*m,8*m*m-24*m+8]
    v=[0]*9
    for d,c in enumerate(e): v[d+1]=2*c
    return v
def lmul(a,b):
    o=[0]*(len(a)+len(b)-1)
    for i,ci in enumerate(a):
        if ci:
            for j,cj in enumerate(b):
                if cj: o[i+j]+=ci*cj
    return o
def lpow(p,n):
    o=[1]; b=p[:]
    while n:
        if n&1: o=lmul(o,b)
        b=lmul(b,b); n>>=1
    return o
def phifac(poly,n): return sum(c*(2**(n-r))*factorial(n-r) for r,c in enumerate(poly) if c and 0<=r<=n)
def lamB1(m):
    if m==2: return 28
    if m==3: return 7320
    return phifac(lmul(lpow([1,-4,2],m-4),L_poly_t(m)),2*m)
def Km(m): return (m-1)*lamB1(m)
def Dm(m): return (m-1)*(2**m*(m*m-m-1))

res={}
# 1) finite replay of the Type-D trace identity (n=5,7)
replay=[]; ok=True
for n in (5,7):
    for row in replay_typeD(n):
        replay.append({kk:(str(vv) if isinstance(vv,F) else vv) for kk,vv in row.items()})
        for key in ("T1","T2","T2_lam_uniform_all_pairs","T3_mu_le_size","traceid","csplit_match","lamD_pos"):
            if row[key] is not True: ok=False
        if row["k"]>=2 and row["k_ge2_bound"] is not True: ok=False
        if row["k"]==2 and row["k2_witness"] is not True: ok=False
res["typeD_trace_replay_ok_n5_7"]=ok
res["replay_rows"]=replay
# 2) k=1 S9C comparison induction
k1=[]; k1ok=True
for m in range(2,40):
    km=Km(m); dm=Dm(m); rho=F(Km(m+1),km); dr=F(Dm(m+1),dm)
    row_ok = km>dm and (rho>=4*m*m) and (m<3 or 4*m*m>dr)
    if not (km>dm): k1ok=False
    if m<8: k1.append(dict(m=m,Km=km,Dm=dm,Km_gt_Dm=km>dm,rho=str(rho),thr=4*m*m,
                           rho_ge_4m2=rho>=4*m*m,Dratio=str(dr),step_ok=rho>dr))
res["k1_Km_gt_Dm_all_m2_39"]=k1ok
res["k1_base_rows"]=k1
res["k1_step_rho_gt_Dratio_all"]=all(F(Km(m+1),Km(m))>F(Dm(m+1),Dm(m)) for m in range(2,60))
res["k1_Dratio_lt_4m2_for_m_ge_3"]=all(2*m*(m*m+m-1)<4*m*m*(m-1)*(m*m-m-1) for m in range(3,4000))
# 3) EXPLICIT all-m parity-even k=2 witness family (n=2m+1, c=m):
#    pi swaps 1<->rho(1), fixes all others; eps=+1 on {0,c,1}, -1 elsewhere.
#    => positive id-fixed {0,c}=2, positive rev-hits {1,c}=2, negatives=2m-2 (even),
#       NegFixSum = 2m-3 > 0.  Verified below for a spread of n.
def k2_witness_family(n):
    m=(n-1)//2; c=m; rho=s3.reversal(n)
    pi=list(range(n)); pi[1]=rho[1]; pi[rho[1]]=1
    eps=[-1]*n
    for i in (0,c,1): eps[i]=1
    pi=tuple(pi); eps=tuple(eps)
    k,l=s3.positive_hit_counts(pi,eps,rho)
    negs=sum(1 for e in eps if e==-1); nf=s3.neg_fixed_count(pi,eps)
    return dict(n=n,k=k,l=l,even_negs=(negs%2==0),neg_fixed=nf,
                valid=(k==2 and l==2 and negs%2==0 and nf>0))
wit=[k2_witness_family(n) for n in (5,7,9,11,13,15,21,31,51)]
res["k2_parity_even_witness_family"]=wit
res["k2_witness_family_all_valid"]=all(w["valid"] for w in wit)

# 4) Type-D channel ranks by the paper's convention (rank_ref=rank(m_ref);
#    rank_pair-std=rank(P m_pair P), P=I-J/n).  Covers the n=3 base case explicitly.
import math
def typeD_ranks(n):
    layers=s3.build_layers(n); out=[]
    for k in range(1, n//2+2):
        D=[(pi,eps) for (pi,eps) in layers.get((k,k),[]) if chi_even(eps)]
        if not D: continue
        mref,mpair=s3.channel_matrices(D,n)
        rr=s3.rank_exact(mref); rp=s3.rank_exact(s3.standard_projected_matrix(mpair))
        ce=math.ceil(n/2)
        out.append(dict(n=n,k=k,size=len(D),rank_ref=rr,rank_pair_std=rp,ok=(rr==ce and rp==ce-1)))
    return out
rankrows=[]
for n in (3,5,7): rankrows+=typeD_ranks(n)
res["typeD_rank_rows_n3_5_7"]=rankrows
res["typeD_ranks_match_ceil_all"]=all(r["ok"] for r in rankrows)

print(json.dumps({k:v for k,v in res.items() if k not in ("replay_rows","k1_base_rows")},indent=2,default=str))
print("\n-- replay rows --")
for r in res["replay_rows"]: print(" ",r)
print("-- k1 base --")
for r in res["k1_base_rows"]: print(" ",r)
print("-- k2 witness family (n=2m+1: pi swaps 1<->rho(1); eps+ on {0,c,1}) --")
for r in res["k2_parity_even_witness_family"]: print(" ",r)
print("-- Type-D ranks (n=3,5,7; rank_ref, rank_pair-std vs ceil(n/2), ceil(n/2)-1) --")
for r in res["typeD_rank_rows_n3_5_7"]: print(" ",r)
allb=[v for k,v in res.items() if isinstance(v,bool)]
allpass=all(allb)
print("\nALL BOOLEAN CHECKS PASS:", allpass, f"({sum(allb)}/{len(allb)})")
res={"gate":"P15D-G22-PSDOM-TYPED-TRACE","status":"PASS" if allpass else "FAIL",
     "pass":allpass,"checks_passed":sum(allb),"checks_total":len(allb),**res}

import argparse
ap=argparse.ArgumentParser()
ap.add_argument("--no-write",action="store_true",help="skip writing the JSON certificate")
ap.add_argument("--out",default=None,help="output JSON path (default results/P15D_G22_PSDOM_TYPED_TRACE.json)")
args,_=ap.parse_known_args()
if not args.no_write:
    out=args.out or os.path.join(os.path.dirname(__file__),"..","results","P15D_G22_PSDOM_TYPED_TRACE.json")
    try:
        os.makedirs(os.path.dirname(out),exist_ok=True)
        with open(out,"w",encoding="utf-8") as f:
            json.dump(res,f,indent=2,default=str)
        print("JSON ->",out)
    except OSError as e:
        print(f"[warn] could not write JSON ({e}); results are in stdout above. Use --no-write to skip.")
sys.exit(0 if allpass else 1)
