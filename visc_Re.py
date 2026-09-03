import numpy as np
from scipy.optimize import minimize
import pb_viscous as V
l,eps,Om=3.0,0.25,0.7
def refine(M,Re,seed,N):
    def obj(z): return V.sigmin(z[0]+1j*z[1],Om,l,eps,M,Re,N)
    r=minimize(obj,[seed.real,seed.imag],method='Nelder-Mead',
               options={'xatol':1e-8,'fatol':1e-14,'maxiter':600})
    return r.x[0]+1j*r.x[1]
def Nfor(Re): return 110 if Re<=3e3 else (140 if Re<=1.2e4 else 170)
print("Im(xi) [decay] vs Re:  below M_crit (M=0.70) vs above (M=0.92)")
print(f"{'Re':>8} {'Imxi(M=0.70)':>14} {'Imxi(M=0.92)':>14} {'CL part(0.92-0.70)':>20}")
seedb=complex(0.844,0.0003); seeda=complex(0.841,0.0002)
for Re in [1e3,3e3,1e4]:
    N=Nfor(Re)
    xb=refine(0.70,Re,seedb,N); xa=refine(0.92,Re,seeda,N); seedb,seeda=xb,xa
    print(f"{Re:>8.0e} {xb.imag:>+14.6f} {xa.imag:>+14.6f} {(xa.imag-xb.imag):>+20.6f}")
print("\nExpect: below-crit Im~Re^-1/2 ->0 ; above-crit retains a finite critical-layer absorption.")
# check Re^-1/2 scaling of the sub-critical (Stokes) damping
print("Re^-1/2 ratios (M=0.70): should be ~constant if Stokes-dominated")
