"""
Numerical validation of every asymptotic formula in Asymptotic_Catalogue.md
against exact roots of the coupled dispersion relation. Run: python3 validate_catalogue.py
"""
import numpy as np, dispersion_two_plates as M
l=3.0
taus=lambda O:l*np.sqrt(complex(O**2-O))
sym_struct =lambda O,e:np.sqrt(O)+(-1/(4*l*np.sqrt(complex(O-1))*np.tan(taus(O)/2))).real*e
anti_struct=lambda O,e:np.sqrt(O)+(np.tan(taus(O)/2)/(4*l*np.sqrt(complex(O-1)))).real*e
sym_plane  =lambda O,e:O+e/(O*l**2*(O**2-1))
def cuton(O,e,mm):
    x0=np.sqrt(O**2-(mm*np.pi/l)**2); S0=x0**4/O**2-1; return x0+2*e/(S0*l**2*x0)
def pr_mode(O,e,Mp):
    x0=np.sqrt(O**2-(Mp*np.pi/l)**2); S0=x0**4/O**2-1; return x0-2*(Mp*np.pi)**2*S0/(l**2*x0)/e
def nearest(res,O,g,e,**k):
    r=M.roots(res,O,l,e,**k); return min(r,key=lambda z:abs(z-g)) if r else np.nan
def two(res,O,g,e,**k):
    r=sorted(M.roots(res,O,l,e,**k),key=lambda z:abs(z-g))[:2]; return sorted(r)
def band(name,res,asy,Omv,e,**k):
    er=[abs(nearest(res,O,asy(O,e),e,**k)-asy(O,e))/abs(asy(O,e)) for O in Omv]
    er=np.array([x for x in er if np.isfinite(x)])
    print(f"{name:28s} eps={e}: max={er.max()*100:6.3f}%  mean={er.mean()*100:6.3f}%")

print("=== AWAY branches (eps=0.05) ===")
e=0.05
band("SYM struct  Om>1",M.res_sym,sym_struct,np.linspace(1.2,1.55,10),e)
band("ANTI struct Om<1",M.res_asym,anti_struct,np.linspace(0.4,0.9,10),e)
band("SYM plane wave",M.res_sym,sym_plane,np.linspace(1.3,2.3,10),e)
band("SYM cut-on m=2",M.res_sym,lambda O,e:cuton(O,e,2),np.linspace(2.25,2.5,8),e)

print("\n=== LOCAL sqrt(eps) (eps=0.03) ===")
e=0.03
ex=two(M.res_sym,1.0,1.0,e); print("SYM coincidence  pred",[round(1+s*np.sqrt(e)/(2*l),4) for s in(-1,1)],"exact",[round(x,4) for x in ex])
Oc=np.roots([l**2,-l**2,-(2*np.pi)**2]); Oc=Oc[Oc>0][0].real; j=np.sqrt(2*e)/(2*l)
ex=two(M.res_sym,Oc,np.sqrt(Oc),e); print(f"SYM cut-on Oc={Oc:.3f} pred",[round(np.sqrt(Oc)+s*j,4) for s in(-1,1)],"exact",[round(x,4) for x in ex])

print("\n=== LARGE eps pressure-release (eps=20) ===")
e=20.0
for Mp,res,O in [(1,M.res_sym,1.7),(2,M.res_asym,2.5)]:
    p=pr_mode(O,e,Mp); x=nearest(res,O,p,e); print(f"{'SYM' if Mp%2 else 'ANTI'} M={Mp} Om={O}: pred={p:.4f} exact={x:.4f} err={abs(x-p)/abs(x)*100:.2f}%")
