from math import comb
from fractions import Fraction as F
import sys, json, os
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass
# ============================================================================
# G15: CLOSURE of (C) [odd sign-localization: s^chi_{2m+1,k}>=0 for 3k>2m]
# via the monotonicity of T(k)=S(k,k)-S(k-2,k)  (s^chi=T(k)-T(k-1)).
# Proof chain (all elementary):
#   (SD) support-incompatibility  => T(k)<=0 on band, sign-definite summand.
#   (RB) explicit per-i ratio bound ratio_i <= 2(m-k+1)^2/k^2 < 1 (AM-GM+factor2).
#   (ORPH) single orphan at i0=m+1-k, nonzero only for 3k<=2m+2 (<=2 lowest cells).
#   (OA) B + orphan/|T(k-1)| < 1  (B=2(m-k+1)^2/k^2), monotone ->1/2; finite base.
# ============================================================================
def C(n,k):
    if k<0 or n<0 or k>n: return 0
    return comb(n,k)
def Diff(N,k): return C(2*N+1,k)-C(2*N+1,k-1)              # Diff_i(k), i via N=m-1-i
def absterm(m,i,k): return -(C(m-1,i)*C(2*i,k)*Diff(m-1-i,k))
def Trep(m,k): return sum(-absterm(m,i,k) for i in range(m))
def Sfun(m,a,c):
    if a<0 or c<0: return 0
    return sum(C(m-1,i)*C(2*i,a)*C(2*(m-1-i),c) for i in range(m))
def schi(m,k): return Sfun(m,k,k)-Sfun(m,k-1,k-1)-Sfun(m,k-2,k)+Sfun(m,k-3,k-1)

R=range(4,140); res={}

# 0) master identity: s^chi = T(k)-T(k-1),  T(k)=S(k,k)-S(k-2,k)
id0=all(schi(m,k)==Trep(m,k)-Trep(m,k-1) for m in range(2,60) for k in range(1,m+2))
res["id_schi_eq_dT"]=id0

# 1) (SD): on band 3k>2m, every nonzero summand of T(k) has Diff_i(k)<=0 -> T(k)<=0
sd=True
for m in R:
    for k in range(0,m+2):
        if 3*k<=2*m: continue
        for i in range(m):
            if C(2*i,k)!=0 and Diff(m-1-i,k)>0: sd=False
        if Trep(m,k)>0: sd=False
res["SD_T_nonpos_on_band"]=sd

# 2) (RB): ratio_i <= B=2(m-k+1)^2/k^2 pointwise on common support, and B<1 on band
rb=True; Blt1=True
for m in R:
    for k in range(int(2*m//3)+1,m+2):
        B=F(2*(m-k+1)**2,k*k)
        if not (B<1): Blt1=False
        for i in range(m):
            ak,ak1=absterm(m,i,k),absterm(m,i,k-1)
            if ak==0 or ak1==0: continue
            if F(ak,ak1)>B: rb=False
res["RB_pointwise_ratio_le_B"]=rb
res["RB_B_lt_1_on_band_m>=4"]=Blt1

# 3) (ORPH): single orphan at i0=m+1-k; nonzero only when 3k<=2m+2
orph_single=True; orph_loc=True
for m in R:
    for k in range(int(2*m//3)+1,m+2):
        orph=[i for i in range(m) if absterm(m,i,k)!=0 and absterm(m,i,k-1)==0]
        if orph and orph!=[m+1-k]: orph_single=False
        if orph and not (3*k<=2*m+2): orph_loc=False
res["ORPH_single_at_m+1-k"]=orph_single
res["ORPH_nonzero_only_3k<=2m+2"]=orph_loc

# 4) (OA): B + orphan/|T(k-1)| < 1 on band, exceptions?
oa_exc=[]
for m in R:
    for k in range(int(2*m//3)+1,m+2):
        Tk1=abs(Trep(m,k-1))
        if Tk1==0:      # T(k-1)=0: (4,3)-type degenerate, handle by direct check
            if abs(Trep(m,k))>0 and abs(Trep(m,k))>Tk1 and (m,k)!=(4,3): oa_exc.append((m,k,"Tk1=0"))
            continue
        B=F(2*(m-k+1)**2,k*k); orph=absterm(m,m+1-k,k)
        if B+F(orph,Tk1)>=1: oa_exc.append((m,k,float(B+F(orph,Tk1))))
res["OA_exceptions_m4_139"]=oa_exc

# 5) FINAL: direct s^chi>=0 for 3k>2m (the theorem (C)) exceptions
final_exc=[(m,k) for m in range(2,160) for k in range(0,m+2)
           if 3*k>2*m and schi(m,k)<0]
res["C_direct_exceptions_m<160"]=final_exc

allb=[v for v in res.values() if isinstance(v,bool)]
allpass=all(allb)
# top-level status wrapper (uniform with the other certificates)
res={"gate":"P15D-G15-C-CLOSURE","status":"PASS" if allpass else "FAIL",
     "pass":allpass,"checks_passed":sum(allb),"checks_total":len(allb),**res}
print(json.dumps(res,indent=2,default=str))
print("\nALL BOOLEAN CHECKS PASS:", allpass, f"({sum(allb)}/{len(allb)})")

import argparse
ap=argparse.ArgumentParser()
ap.add_argument("--no-write",action="store_true",help="skip writing the JSON certificate")
ap.add_argument("--out",default=None,help="output JSON path")
args,_=ap.parse_known_args()
if not args.no_write:
    out=args.out or os.path.join(os.path.dirname(__file__),"..","results","P15D_G15_C_CLOSURE.json")
    try:
        os.makedirs(os.path.dirname(out),exist_ok=True)
        with open(out,"w",encoding="utf-8") as f:
            json.dump(res,f,indent=2,default=str)
        print("JSON ->",out)
    except OSError as e:
        print(f"[warn] could not write JSON ({e}); results are in stdout. Use --no-write to skip.")
sys.exit(0 if allpass else 1)
