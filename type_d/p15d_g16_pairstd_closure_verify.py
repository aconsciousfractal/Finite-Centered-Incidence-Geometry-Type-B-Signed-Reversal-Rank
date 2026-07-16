from fractions import Fraction as F
from math import comb
from itertools import combinations
import sys, json, os
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass
# ============================================================================
# G16: ODD PAIR-STANDARD channel CLOSED.  lam'^D=1/2(lam'^B+lam'^chi).
#   Closed form   lam'^chi_{2m+1,k}=[A^kC^k](AC-1)Q^m = S'(k-1,k-1)-S'(k,k)
#                 (= diagonal of Delta_odd),  S'(a,c)=[A^aC^c]((1+A)^2+(1+C)^2)^m.
#   |lam'^chi|<=2*8^m;   lam'^chi(m,1)=-2^m(m^2-m-1).
#   Domination reduces to ONE pure Type-B count inequality  lam'^B >= S'(k,k):
#     chi>=0 cells: lam'^D>=1/2 lam'^B>0 (P15-S9, free);
#     chi<0  cells: |lam'^chi|=S'(k,k)-S'(k-1,k-1)<S'(k,k)<=lam'^B => lam'^D>0.
# ============================================================================
def pmul(p,q):
    r={}
    for (a1,c1),v1 in p.items():
        for (a2,c2),v2 in q.items():
            r[(a1+a2,c1+c2)]=r.get((a1+a2,c1+c2),0)+v1*v2
    return r
def permanent_poly(M,n):
    total={}; full=tuple(range(n))
    for r in range(0,n+1):
        for S_ in combinations(full,r):
            prod={(0,0):1}; zero=False
            for i in range(n):
                rs={}
                for j in S_:
                    for mm,c in M[i][j].items(): rs[mm]=rs.get(mm,0)+c
                rs={mm:c for mm,c in rs.items() if c}
                if not rs: zero=True; break
                prod=pmul(prod,rs)
                if not prod: zero=True; break
            if zero: continue
            sgn=-1 if (r%2) else 1
            for mm,c in prod.items(): total[mm]=total.get(mm,0)+sgn*c
    if n%2==1: total={mm:-c for mm,c in total.items()}
    return {mm:c for mm,c in total.items() if c}
def body(n,chi):
    center=(n+1)//2 if n%2==1 else None; rho=lambda i:n+1-i; M=[]
    for i in range(1,n+1):
        row=[]
        for j in range(1,n+1):
            isd=(j==i); isa=(j==rho(i)); isc=(center is not None and i==center and j==center)
            if chi: e=({(1,1):1,(0,0):-1} if isc else {(1,0):1,(0,0):-1} if isd else {(0,1):1,(0,0):-1} if isa else {})
            else:   e=({(1,1):1,(0,0):1} if isc else {(1,0):1,(0,0):1} if isd else {(0,1):1,(0,0):1} if isa else {(0,0):2})
            row.append(e)
        M.append(row)
    return M
def forced(n,chi,col,entry):
    M=body(n,chi); nr=[{} for _ in range(n)]; nr[col-1]=entry; M[0]=nr; return M
def coeff(M,n,k): return permanent_poly(M,n).get((k,k),0)
def lamB(n,k):
    a=coeff(forced(n,False,1,{(1,0):1,(0,0):1}),n,k)
    b=coeff(forced(n,False,n,{(0,1):1,(0,0):1}),n,k)
    c=coeff(forced(n,False,2,{(0,0):2}),n,k)
    return a+b-2*c
def lamChi_machine(n,k):
    a=coeff(forced(n,True,1,{(1,0):1,(0,0):-1}),n,k)
    b=coeff(forced(n,True,n,{(0,1):1,(0,0):-1}),n,k)
    return a+b
def Sp(m,a,c):
    if a<0 or c<0: return 0
    return sum(comb(m,i)*comb(2*i,a)*comb(2*(m-i),c) for i in range(m+1))
def qm(m,a,c): return ((-1)**(a+c))*Sp(m,a,c) if a>=0 and c>=0 else 0

res={}; NS=[5,7,9,11,13]
# closed forms
res["lamChi_eq_ACminus1_Qm"]=all(lamChi_machine(n,k)==qm(n//2,k-1,k-1)-qm(n//2,k,k) for n in NS for k in range(1,n))
res["lamChi_eq_Sform"]=all(lamChi_machine(n,k)==Sp(n//2,k-1,k-1)-Sp(n//2,k,k) for n in NS for k in range(1,n))
res["lamChi_k1_closed"]=all((Sp(m,0,0)-Sp(m,1,1))==-(2**m)*(m*m-m-1) for m in range(2,40))
res["absLamChi_le_2_8m"]=all(abs(Sp(n//2,k-1,k-1)-Sp(n//2,k,k))<=2*8**(n//2) for n in NS for k in range(1,n))
# domination
res["lamB_ge_Spkk"]=all(lamB(n,k)>=Sp(n//2,k,k) for n in NS for k in range(1,n) if lamB(n,k)!=0 or Sp(n//2,k,k)!=0 or lamChi_machine(n,k)!=0)
res["lamB_pos"]=all(lamB(n,k)>0 for n in NS for k in range(1,n) if lamB(n,k)!=0 or lamChi_machine(n,k)!=0)
res["lamD_pos"]=all(F(lamB(n,k)+lamChi_machine(n,k),2)>0 for n in NS for k in range(1,n) if lamB(n,k)!=0 or lamChi_machine(n,k)!=0)
res["dom_on_negchi"]=all(lamB(n,k)>abs(lamChi_machine(n,k)) for n in NS for k in range(1,n) if lamChi_machine(n,k)<0)
# margins
worst=None
for n in NS:
    for k in range(1,n):
        lx=lamChi_machine(n,k)
        if lx<0:
            r=F(lamB(n,k),-lx)
            if worst is None or r<worst[0]: worst=(r,n,k)
res["worst_domination_margin"]=f"{float(worst[0]):.3f} at (n,k)={worst[1],worst[2]}"
res["worst_lamB_over_Spkk"]=f"{float(min(F(lamB(n,k),Sp(n//2,k,k)) for n in NS for k in range(1,n) if Sp(n//2,k,k)>0)):.3f}"
allb=[v for v in res.values() if isinstance(v,bool)]
allpass=all(allb)
# top-level status wrapper (uniform with the other certificates)
res={"gate":"P15D-G16-PAIRSTD-CLOSURE","status":"PASS" if allpass else "FAIL",
     "pass":allpass,"checks_passed":sum(allb),"checks_total":len(allb),**res}
print(json.dumps(res,indent=2,default=str))
print("\nALL BOOLEAN CHECKS PASS:", allpass, f"({sum(allb)}/{len(allb)})")

import argparse
ap=argparse.ArgumentParser()
ap.add_argument("--no-write",action="store_true",help="skip writing the JSON certificate")
ap.add_argument("--out",default=None,help="output JSON path")
args,_=ap.parse_known_args()
if not args.no_write:
    out=args.out or os.path.join(os.path.dirname(__file__),"..","results","P15D_G16_PAIRSTD_CLOSURE.json")
    try:
        os.makedirs(os.path.dirname(out),exist_ok=True)
        with open(out,"w",encoding="utf-8") as f:
            json.dump(res,f,indent=2,default=str)
        print("JSON ->",out)
    except OSError as e:
        print(f"[warn] could not write JSON ({e}); results are in stdout. Use --no-write to skip.")
sys.exit(0 if allpass else 1)
