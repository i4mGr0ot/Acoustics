"""
Closed-form O(M) correction xi = xi0 + M*xi1 + O(M^2) for the Poiseuille two-plate waveguide,
derived by Fredholm solvability of the O(M) Pridmore-Brown problem (leading operator self-adjoint).

xi1 = [ (2/Om)(p0(1)^2+p0(0)^2) - (4/Om) I0 - l^2 Om If ]
      / [ l^2 I0 + (2 eps xi0^2/Om^2)( rD p0(1)^2/S20^2 + p0(0)^2/S10^2 ) ]
with I0=∫p0^2, If=∫ f p0^2, f=4η(1-η),  p0=cos(τ0 η)+ (B/A) sin(τ0 η), B/A=-1/alpha1.
"""
import numpy as np
from scipy.integrate import quad
import pb_solver as P, num_solver as NS

def xi1_closed(Om,l,eps,xi0,rD=1.0,rm=1.0):
    S10=xi0**4/Om**2-1.0; S20=rD*xi0**4/Om**2-rm
    tau0=l*np.sqrt(complex(Om**2-xi0**2))
    alpha1=S10*tau0/eps
    BoverA=(-1.0/alpha1)
    def p0(e):
        return (np.cos(tau0*e)+BoverA*np.sin(tau0*e)).real if abs(tau0.imag)<1e-12 \
               else (np.cosh((1j*tau0).real*e)+ (BoverA*np.sin(tau0*e))).real
    # robust: just use complex cos/sin (handles imaginary tau0 -> cosh/sinh), take real of p0^2
    def P0(e): return np.cos(tau0*e)+BoverA*np.sin(tau0*e)
    f=lambda e:4*e*(1-e)
    I0=quad(lambda e:(P0(e)**2).real,0,1)[0]
    If=quad(lambda e:(f(e)*P0(e)**2).real,0,1)[0]
    p0_0=(P0(0.0)**2).real; p0_1=(P0(1.0)**2).real
    num=-(2/Om)*(p0_1+p0_0)+(4/Om)*I0 - l**2*Om*If   # shear sign-corrected
    den=l**2*I0 + (2*eps*xi0**2/Om**2)*(rD*p0_1/S20**2 + p0_0/S10**2)
    return num/den

if __name__=="__main__":
    l,eps=3.0,0.25
    print("O(M) correction xi1 : closed-form vs numerical dξ/dM (central, M=±0.02)")
    print(f"{'Om':>4} {'xi0':>8} {'xi1_closed':>12} {'xi1_numeric':>12} {'rel.err':>9}")
    for Om in [1.3,1.8,2.2]:
        roots0=sorted(P.pb_roots(Om,l,eps,0.0))
        for xi0 in roots0:
            # numeric dxi/dM by matching nearest root at +-dM
            dM=0.02
            rp=P.pb_roots(Om,l,eps,+dM); rm_=P.pb_roots(Om,l,eps,-dM)
            xp=min(rp,key=lambda r:abs(r-xi0)); xm=min(rm_,key=lambda r:abs(r-xi0))
            xi1_num=(xp-xm)/(2*dM)
            xi1_cf=xi1_closed(Om,l,eps,xi0)
            err=abs(xi1_cf-xi1_num)/max(abs(xi1_num),1e-6)
            print(f"{Om:>4} {xi0:>8.4f} {xi1_cf:>12.4f} {xi1_num:>12.4f} {err*100:>8.2f}%")
