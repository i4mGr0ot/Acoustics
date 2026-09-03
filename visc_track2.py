import numpy as np
from scipy.optimize import minimize
import pb_viscous as V
l,eps,Om=3.0,0.25,0.7
def refine(M,Re,seed,N=120):
    def obj(z): return V.sigmin(z[0]+1j*z[1],Om,l,eps,M,Re,N)
    r=minimize(obj,[seed.real,seed.imag],method='Nelder-Mead',
               options={'xatol':1e-7,'fatol':1e-13,'maxiter':500})
    return r.x[0]+1j*r.x[1], r.fun
seed=complex(0.857,-0.001)
print("Track HIGH structural branch across M_crit (Re=3000, N=120). Decay = Im(xi)>0.")
print(f"{'M':>5} {'Re xi':>8} {'Im xi':>10} {'M*xi/Om':>8} {'regime':>12}")
prev=seed
for M in [0.60,0.72,0.80,0.84,0.88,0.92,0.95]:
    xi,res=refine(M,3000.0,prev); prev=xi
    cl=M*xi.real/Om
    print(f"{M:>5} {xi.real:>8.4f} {xi.imag:>+10.5f} {cl:>8.3f} {'CRITICAL' if cl>1 else 'subcrit':>12}  (r {res:.0e})")
