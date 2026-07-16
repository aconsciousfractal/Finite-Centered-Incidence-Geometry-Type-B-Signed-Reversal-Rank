from math import comb
from fractions import Fraction as F
import sys, json, os
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass
# ============================================================================
# G20: OA-tail closed SYMBOLICALLY (upgrades the G15/G19 (OA) obligation).
#
# On the band 3k>2m the odd sign-scalar is
#     s^chi(m,k) = [A^k C^k] (1-AC)(1-A^2) F^{m-1},   F=(1+A)^2+(1+C)^2,
# and in the SD/RB/ORPH decomposition of G15
#     s^chi = sum_i [ t(i,k-1) - t(i,k) ],   t(i,k):=|absterm(m,i,k)| >= 0 (SD),
#   * RB   : common support => t(i,k) <= B t(i,k-1),  B=2(m-k+1)^2/k^2   (proved G15)
#   * ORPH : exactly one orphan i0=m+1-k with t(i0,k)!=0, t(i0,k-1)=0   (proved G15)
# The ONLY orphan cells are 3k in {2m+1, 2m+2}, i.e. the two families
#     A_r : m=3r+1, k=2r+1   (3k=2m+1),  i0=r+1
#     B_r : m=3r-1, k=2r     (3k=2m+2),  i0=r
# Retaining the single common index i1=i0+1 gives the master lower bound
#     s^chi >= (1-B) t(i1,k-1) - O ,    O:=t(i0,k),
# so s^chi>0 follows from ONE explicit binomial inequality R(r):=(1-B)t(i1,k-1)/O >= 1.
# The ratio has a PROVEN closed rational form (hand-derived; verified here exact):
#     R_A(r) = (2r^2-1)(2r-1)(2r+3) / (24(4r-1))            [valid r>=3]
#     R_B(r) = (2/3) r^2 (r-1)(2r+1) / (4r-3)               [valid r>=2]
# and R>=1 reduces to quartics qA,qB>=0 closed by derivative-positivity + one base.
# The single leftover cell (family A, r=2 = cell (7,5)) is checked directly (=840>0).
# NO non-cancellation / "|T(k-1)|>=max summand" claim is used anywhere.
# ============================================================================
def C(n,k):
    if k<0 or n<0 or k>n: return 0
    return comb(n,k)
def Diff(N,k): return C(2*N+1,k)-C(2*N+1,k-1)
def absterm(m,i,k): return -(C(m-1,i)*C(2*i,k)*Diff(m-1-i,k))   # t(i,k) >= 0 on band
def schi(m,k):
    T=lambda kk: sum(-absterm(m,i,kk) for i in range(m))
    return T(k)-T(k-1)
def B_of(m,k): return F(2*(m-k+1)**2,k*k)

def RA_closed(r): return F((2*r*r-1)*(2*r-1)*(2*r+3), 24*(4*r-1))
def RB_closed(r): return F(2,3)*F(r*r*(r-1)*(2*r+1), (4*r-3))
def R_actual(fam,r):
    m,k=(3*r+1,2*r+1) if fam=="A" else (3*r-1,2*r)
    i0=m+1-k; i1=i0+1
    O=absterm(m,i0,k); ak1=absterm(m,i1,k-1)
    return F((1-B_of(m,k))*ak1,O), O, ak1, (absterm(m,i1,k)!=0 and ak1!=0)
def orphans(m,k): return [i for i in range(m) if absterm(m,i,k)!=0 and absterm(m,i,k-1)==0]
def qA(r): return 8*r**4+8*r**3-10*r**2-100*r+27      # R_A>=1  <=>  qA>=0
def qB(r): return 4*r**4-2*r**3-2*r**2-12*r+9          # R_B>=1  <=>  qB>=0

res={}; RANGE=range(2,220)
# 1. hand-derived closed forms are exact
res["RA_closed_form_exact_r>=3"]=all(R_actual("A",r)[0]==RA_closed(r) for r in range(3,200))
res["RB_closed_form_exact_r>=2"]=all(R_actual("B",r)[0]==RB_closed(r) for r in range(2,200))
# 2. single orphan + i1 common support
res["ORPH_single_A"]=all(orphans(3*r+1,2*r+1)==[r+1] for r in RANGE)
res["ORPH_single_B"]=all(orphans(3*r-1,2*r)==[r] for r in RANGE)
res["i1_common_support_A_r>=3"]=all(R_actual("A",r)[3] for r in range(3,200))
res["i1_common_support_B_r>=2"]=all(R_actual("B",r)[3] for r in range(2,200))
# 3. R>=1 <=> quartics; quartics >=0 by base + monotonicity
res["qA_equiv_RA_ge_1"]=all((RA_closed(r)>=1)==(qA(r)>=0) for r in range(1,300))
res["qB_equiv_RB_ge_1"]=all((RB_closed(r)>=1)==(qB(r)>=0) for r in range(1,300))
res["qA(2)_neg_qA(3)_pos"]=(qA(2)<0 and qA(3)>0)   # r=2 is the base cell, r>=3 general
res["qB(2)_pos"]=(qB(2)>0)
# derivative positivity (closes the tail with no empirical cutoff):
# qA'(r)=32r^3+24r^2-20r-100>0 for r>=3, qA''=96r^2+48r-20>0 for r>=1 => qA up from r=3
# qB'(r)=16r^3-6r^2-4r-12>0 for r>=2, qB''=48r^2-12r-4>0  for r>=1 => qB up from r=2
res["qA'(3)_pos"]=(32*27+24*9-20*3-100>0)
res["qB'(2)_pos"]=(16*8-6*4-4*2-12>0)
# 4. master lower bound  s^chi >= O*(R-1)
def lower_ok(fam,r):
    m,k=(3*r+1,2*r+1) if fam=="A" else (3*r-1,2*r)
    Ract,O,_,_=R_actual(fam,r); return schi(m,k)>=O*(Ract-1)
res["masterbound_A_r>=3"]=all(lower_ok("A",r) for r in range(3,120))
res["masterbound_B_r>=2"]=all(lower_ok("B",r) for r in range(2,120))
# 5. final positivity on orphan cells (A r=2 = base case, =840)
res["schi_pos_A_r>=2"]=all(schi(3*r+1,2*r+1)>0 for r in range(2,200))
res["schi_pos_B_r>=2"]=all(schi(3*r-1,2*r)>0 for r in range(2,200))
res["A2_basecase_schi(7,5)"]=schi(7,5)
res["r1_exceptions"]={"schi(4,3)":schi(4,3),"schi(2,2)":schi(2,2)}
# 6. full (C): only band negatives for m<160 are the r=1 base cells (4,3),(2,2)
res["C_band_negatives_m<160"]=[(m,k) for m in range(2,160) for k in range(0,m+2)
                               if 3*k>2*m and schi(m,k)<0]

allb=[v for v in res.values() if isinstance(v,bool)]
allpass=all(allb)
# top-level status wrapper (uniform with the other certificates)
res={"gate":"P15D-G20-OATAIL","status":"PASS" if allpass else "FAIL",
     "pass":allpass,"checks_passed":sum(allb),"checks_total":len(allb),**res}
print(json.dumps({k:v for k,v in res.items() if k not in ("C_band_negatives_m<160","r1_exceptions")},indent=2,default=str))
print("\nALL BOOLEAN CHECKS PASS:", allpass, f"({sum(allb)}/{len(allb)})")

import argparse
ap=argparse.ArgumentParser()
ap.add_argument("--no-write",action="store_true",help="skip writing the JSON certificate")
ap.add_argument("--out",default=None,help="output JSON path")
args,_=ap.parse_known_args()
if not args.no_write:
    out=args.out or os.path.join(os.path.dirname(__file__),"..","results","P15D_G20_OATAIL_SYMBOLIC.json")
    try:
        os.makedirs(os.path.dirname(out),exist_ok=True)
        with open(out,"w",encoding="utf-8") as f:
            json.dump(res,f,indent=2,default=str)
        print("JSON ->",out)
    except OSError as e:
        print(f"[warn] could not write JSON ({e}); results are in stdout. Use --no-write to skip.")
sys.exit(0 if allpass else 1)
