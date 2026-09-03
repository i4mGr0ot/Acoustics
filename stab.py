import numpy as np, warnings; warnings.filterwarnings('ignore')
from scipy.optimize import minimize
import pb_viscous as V, pb_solver as P
l,eps=3.0,0.25
def refine(Om,M,Re,seed,N=110):
    f=lambda z:V.sigmin(z[0]+1j*z[1],Om,l,eps,M,Re,N)
    r=minimize(f,[seed.real,seed.imag],method='Nelder-Mead',options={'xatol':1e-7,'fatol':1e-13,'maxiter':500})
    return r.x[0]+1j*r.x[1]
print("Convective stability: Im(xi) (>0 = decay/stable) for plane-wave branch Om=1.3, Re=3000")
r0=P.pb_roots(1.3,l,eps,0.0,N=240); seed=complex(min(r0,key=lambda z:abs(z-1.3267)),0.0)
prev=seed
for M in [0.0,0.2,0.4]:
    xi=refine(1.3,M,3000.0,prev); prev=xi
    print(f"  M={M}: xi={xi.real:.4f}{xi.imag:+.5f}i  -> {'STABLE (decay)' if xi.imag>0 else 'UNSTABLE (growth)'}")
