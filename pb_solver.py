"""
Pridmore-Brown dispersion solver: fluid column 0<=y<=a bounded by TWO flexible plates,
carrying a Poiseuille mean flow U(y)=4 U0 eta(1-eta), eta=y/a, no-slip U(0)=U(a)=0.

Nondim PB (eta in [0,1]):
   p'' - (2 M xi f'/shat) p' + l^2 (shat^2 - xi^2) p = 0,   shat=Omega - M xi f(eta)
   f=4 eta(1-eta), f'=4(1-2 eta),  M=U0/c.
Plate-admittance (Robin) BCs (derived; U=0 at walls so no Ingard-Myers slip term):
   p'(0) + (eps/S1) p(0) = 0,   p'(1) - (eps/S2) p(1) = 0
   S1 = xi^4/Om^2 - 1,  S2 = rD xi^4/Om^2 - rm.
Solved by shooting in eta + root finding the wall residual over xi.
"""
import numpy as np
from scipy.integrate import solve_ivp

def f_prof(eta):  return 4*eta*(1-eta)
def fp_prof(eta): return 4*(1-2*eta)

def _rhs(eta,Y,Om,l,xi,M,shear=True):
    p,q=Y
    sh=Om - M*xi*f_prof(eta)
    dq=-(2*M*xi*fp_prof(eta)/sh if shear else 0.0)*q - l**2*(sh**2 - xi**2)*p
    return [q,dq]

def wall_residual(xi,Om,l,eps,M,rD=1.0,rm=1.0,shear=True):
    S1=xi**4/Om**2-1.0; S2=rD*xi**4/Om**2-rm
    Y0=[1.0, -(eps/S1)]                       # BC1: p(0)=1, p'(0)=-(eps/S1)
    sol=solve_ivp(_rhs,(0,1),Y0,args=(Om,l,xi,M,shear),rtol=1e-7,atol=1e-9,max_step=0.06)
    p1,q1=sol.y[0,-1],sol.y[1,-1]
    return q1 - (eps/S2)*p1                   # BC2 residual

def pb_roots(Om,l,eps,M,rD=1.0,rm=1.0,xmax=3.0,N=600,shear=True):
    xs=np.linspace(1e-2,xmax,N)
    vals=np.array([wall_residual(x,Om,l,eps,M,rD,rm,shear) for x in xs])
    r=[]
    for i in range(len(xs)-1):
        a,b=vals[i],vals[i+1]
        if np.isfinite(a) and np.isfinite(b) and a*b<0 and abs(a)<1e4 and abs(b)<1e4:
            lo,hi,flo=xs[i],xs[i+1],a
            for _ in range(45):
                mid=0.5*(lo+hi); fm=wall_residual(mid,Om,l,eps,M,rD,rm,shear)
                if flo*fm<=0: hi=mid
                else: lo=mid; flo=fm
            r.append(0.5*(lo+hi))
    return r

if __name__=="__main__":
    import num_solver as M0
    l,eps=3.0,0.25
    print("Validation at M=0 : PB shooting  vs  exact quiescent (full relation)")
    for Om in [0.7,1.3,1.8]:
        pb=sorted(pb_roots(Om,l,eps,0.0))
        ex=sorted(M0.roots(M0.res_general,Om,l,eps))
        print(f" Om={Om}: PB={[round(x,4) for x in pb]}")
        print(f"        exact={[round(x,4) for x in ex]}")
