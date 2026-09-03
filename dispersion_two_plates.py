"""Coupled dispersion of a fluid column bounded by TWO flexible plates.
Extension of Sarkar & Sonti (JSV 306, 2007). See report for derivation.

Residuals are written so that the propagating (xi<Omega) and evanescent (xi>Omega)
forms join with the SAME SIGN at xi=Omega, so the bracketing root-finder does not
see a spurious sign jump there (which would hide the evanescent coupled plane wave).
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def S_of(xi,Om): return xi**4/Om**2 - 1.0

def res_general(xi,Om,l,eps,rD=1.0,rm=1.0):
    """Full non-identical relation tan(tau)=(a1+a2)/(1-a1 a2), real residual."""
    xi=np.asarray(xi,float); d=Om**2-xi**2
    S1=S_of(xi,Om); S2=rD*xi**4/Om**2-rm
    out=np.empty_like(xi); p=d>=0; e=~p
    tau=l*np.sqrt(np.abs(d))
    a1=S1*tau/eps; a2=S2*tau/eps
    out[p]=(np.sin(tau)*(1-a1*a2)-np.cos(tau)*(a1+a2))[p]
    b1=S1*tau/eps; b2=S2*tau/eps                       # tau=i T  -> a_j=i b_j
    out[e]=(np.tanh(tau)*(1+b1*b2)-(b1+b2))[e]
    return out

def res_sym(xi,Om,l,eps):           # v=0 at centre ; alpha=-cot(tau/2)
    # residual = [S*tau*tan(tau/2)+eps]*cos(tau/2)  (sign-continuous across xi=Om)
    xi=np.asarray(xi,float); d=Om**2-xi**2; S=S_of(xi,Om); m=l*np.sqrt(np.abs(d))
    return np.where(d>=0, S*m*np.sin(m/2)+eps*np.cos(m/2),
                          eps*np.cosh(m/2)-S*m*np.sinh(m/2))   # evan tau=iT

def res_asym(xi,Om,l,eps):          # p=0 at centre ; alpha=+tan(tau/2)
    # residual = [S*tau*cot(tau/2)-eps]*sin(tau/2)
    xi=np.asarray(xi,float); d=Om**2-xi**2; S=S_of(xi,Om); m=l*np.sqrt(np.abs(d))
    return np.where(d>=0, S*m*np.cos(m/2)-eps*np.sin(m/2),
                          S*m*np.cosh(m/2)-eps*np.sinh(m/2))

def res_SS(xi,Om,l,eps):            # Sarkar-Sonti single plate ; alpha=-cot(tau)
    xi=np.asarray(xi,float); d=Om**2-xi**2; S=S_of(xi,Om); m=l*np.sqrt(np.abs(d))
    return np.where(d>=0, S*m*np.sin(m)+eps*np.cos(m),
                          eps*np.cosh(m)-S*m*np.sinh(m))       # evan tau=iT

def roots(fun,Om,l,eps,xmax=3.5,N=4000,**kw):
    xs=np.linspace(1e-3,xmax,N); v=fun(xs,Om,l,eps,**kw)
    sg=np.sign(v); idx=np.where((sg[:-1]*sg[1:]<0)&np.isfinite(v[:-1])&np.isfinite(v[1:]))[0]
    r=[]
    for i in idx:
        lo,hi=xs[i],xs[i+1]; flo=v[i]
        if abs(v[i])>1e3 and abs(v[i+1])>1e3:   # skip jumps across tan/cot poles
            continue
        for _ in range(60):
            mid=0.5*(lo+hi); fm=fun(np.array([mid]),Om,l,eps,**kw)[0]
            if flo*fm<=0: hi=mid
            else: lo=mid; flo=fm
        r.append(0.5*(lo+hi))
    return r

def sweep(fun,l,eps,Omv,**kw):
    pts=[(Om,x) for Om in Omv for x in roots(fun,Om,l,eps,**kw)]
    return np.array(pts) if pts else np.empty((0,2))

def uncoupled(Omv,l):
    return np.sqrt(Omv), Omv.copy()
