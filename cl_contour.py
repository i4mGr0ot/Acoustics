"""
Critical-layer eigenvalue by causal (Landau) contour deformation.
Integrate the Pridmore-Brown ODE along a complex-eta path that indents around the
critical layer(s) where sigma_hat=Om-M xi f=0.  Yields a COMPLEX wavenumber xi:
Im xi<0 => the critical layer absorbs acoustic energy (spatial decay).
"""
import numpy as np
l,eps=3.0,0.25
def f(e):  return 4*e*(1-e)
def fp(e): return 4*(1-2*e)

def contour(t,A,w):           # eta(t) dips below real axis near t=0.5
    g=-A*np.exp(-((t-0.5)/w)**2)
    gp=-A*np.exp(-((t-0.5)/w)**2)*(-2*(t-0.5)/w**2)
    return t+1j*g, 1.0+1j*gp

def integrate(xi,Om,M,A,w,Nt=1200):
    S1=xi**4/Om**2-1; 
    P=1.0+0j; Q=-(eps/S1)+0j
    ts=np.linspace(0,1,Nt); dt=ts[1]-ts[0]
    def deriv(t,P,Q):
        e,ep=contour(t,A,w); sg=Om-M*xi*f(e)
        dP=Q*ep
        dQ=(-(2*M*xi*fp(e)/sg)*Q - l**2*(sg**2-xi**2)*P)*ep
        return dP,dQ
    for i in range(Nt-1):
        t=ts[i]
        k1P,k1Q=deriv(t,P,Q)
        k2P,k2Q=deriv(t+dt/2,P+dt/2*k1P,Q+dt/2*k1Q)
        k3P,k3Q=deriv(t+dt/2,P+dt/2*k2P,Q+dt/2*k2Q)
        k4P,k4Q=deriv(t+dt,P+dt*k3P,Q+dt*k3Q)
        P=P+dt/6*(k1P+2*k2P+2*k3P+k4P); Q=Q+dt/6*(k1Q+2*k2Q+2*k3Q+k4Q)
    S2=xi**4/Om**2-1
    return Q-(eps/S2)*P     # wall residual at eta=1

def newton(Om,M,xi0,A,w,iters=40,tol=1e-11):
    xi=complex(xi0)
    for _ in range(iters):
        R=integrate(xi,Om,M,A,w); h=1e-7
        dR=(integrate(xi+h,Om,M,A,w)-integrate(xi-h,Om,M,A,w))/(2*h)
        step=R/dR; xi=xi-step
        if abs(step)<tol: break
    return xi

Om=0.7
# eigenvalue just BELOW M_crit (real) as seed
import pb_solver as P
M0=0.78
r=P.pb_roots(Om,l,eps,M0,N=220); xi_seed=min(r,key=lambda z:abs(z-np.sqrt(Om)))
print(f"Om={Om}: real structural xi at M={M0} (below M_crit) = {xi_seed:.4f}")
Mc=Om/xi_seed; print(f"approx M_crit=Om/xi={Mc:.3f}")

print("\nAbove M_crit: complex eigenvalue via contour deformation (A=dip depth, w=width)")
for M in [0.85,0.90]:
    xb=newton(Om,M,xi_seed,A=0.14,w=0.16)    # contour BELOW critical layers
    xa=newton(Om,M,xi_seed,A=-0.14,w=0.16)   # contour ABOVE (opposite indentation)
    print(f" M={M}: xi(below)= {xb.real:.4f}{xb.imag:+.4f}i   xi(above)= {xa.real:.4f}{xa.imag:+.4f}i")
    print(f"        jump (below-above) = {(xb-xa).real:+.4f}{(xb-xa).imag:+.4f}i   (residue/-i*pi signature)")

if __name__=="__main__" and True:
    print("\n=== Continuation through M_crit (Om=0.85, high structural branch) ===")
    import pb_solver as P
    Om=0.85
    r=sorted(P.pb_roots(Om,l,eps,0.7,N=300)); xi_seed=max(r)  # high structural ~0.93
    print(f"seed: M=0.7 real xi={xi_seed:.4f}, M_crit=Om/xi={Om/xi_seed:.3f}")
    xi=complex(xi_seed); prev=xi
    for M in [0.70,0.80,0.88,0.91,0.93,0.96]:
        xi=newton(Om,M,prev,A=0.13,w=0.18)
        prev=xi
        cl = (M*xi.real>Om)
        print(f" M={M}: xi={xi.real:.4f}{xi.imag:+.5f}i  {'(CRITICAL LAYER: Im xi<0 = absorption)' if cl else '(no CL)'}")
