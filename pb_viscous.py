"""
VISCOUS critical-layer eigensolver for the two-flexible-plate waveguide with Poiseuille flow.
Linearized compressible Navier-Stokes (viscous momentum, inviscid energy), primitive (u,v,p),
Chebyshev collocation, nonlinear eigenvalue in complex xi.  Re = c a / nu (acoustic Reynolds).

Equations (nondim; sigma=Om-M xi f, f=4 eta(1-eta)):
 (a) -i sigma u + (M/l) f' v + i xi p - (1/(l Re))(D2 - l^2 xi^2) u = 0
 (b) -i sigma v + (1/l) p_eta       - (1/(l Re))(D2 - l^2 xi^2) v = 0
 (c)  i xi u - i sigma p + (1/l) v_eta = 0
BCs:  u(0)=u(1)=0 (no-slip);  v - beta1 v_eta=0 at eta=0 ;  v + beta2 v_eta=0 at eta=1,
      beta_j = eps/(l^2 Om^2 S_j)  (plate admittance; reduces to inviscid PB Robin BC as Re->inf).
"""
import numpy as np

def cheb(N):
    x=np.cos(np.pi*np.arange(N+1)/N)
    c=np.hstack([2.,np.ones(N-1),2.])*(-1)**np.arange(N+1)
    X=np.tile(x,(N+1,1)).T; dX=X-X.T
    D=np.outer(c,1./c)/(dX+np.eye(N+1)); D-=np.diag(D.sum(1))
    return D,x

def build(xi,Om,l,eps,M,Re,N,rD=1.0,rm=1.0):
    D,x=cheb(N); eta=(x+1)/2; D1=2*D; D2=D1@D1
    n=N+1; I=np.eye(n); Z=np.zeros((n,n))
    f=4*eta*(1-eta); fp=4*(1-2*eta)
    sig=Om-M*xi*f
    Sg=np.diag(sig); F=np.diag(f); Fp=np.diag(fp)
    visc=(1/(l*Re))*(D2 - (l**2*xi**2)*I)
    # (a): -i Sg u + (M/l) Fp v + i xi p - visc u
    Aa=np.hstack([-1j*Sg - visc, (M/l)*Fp, 1j*xi*I])
    # (b): -i Sg v + (1/l) D1 p - visc v
    Ab=np.hstack([Z, -1j*Sg - visc, (1/l)*D1])
    # (c): i xi u - i Sg p + (1/l) D1 v
    Ac=np.hstack([1j*xi*I, (1/l)*D1, -1j*Sg])
    A=np.vstack([Aa,Ab,Ac]).astype(complex)
    S1=xi**4/Om**2-1; S2=rD*xi**4/Om**2-rm
    b1=eps/(l**2*Om**2*S1); b2=eps/(l**2*Om**2*S2)
    # node indexing: eta from x: x[0]=1 -> eta=1 (top), x[N]=0 -> eta=0 (bottom)
    top=0; bot=N
    # BC rows: u(top)=0, u(bot)=0  -> rows in (a) block (rows top, bot)
    A[top,:]=0; A[top,top]=1.0
    A[bot,:]=0; A[bot,bot]=1.0
    # v Robin: at bottom (eta=0): v - b1 v_eta =0 ; row in (b) block = n+bot
    rb=n+bot; A[rb,:]=0; A[rb,n:2*n]=-b1*D1[bot,:]; A[rb,n+bot]+=1.0
    # at top (eta=1): v + b2 v_eta=0 ; row n+top
    rt=n+top; A[rt,:]=0; A[rt,n:2*n]=b2*D1[top,:]; A[rt,n+top]+=1.0
    return A

def sigmin(xi,Om,l,eps,M,Re,N,**kw):
    return np.linalg.svd(build(xi,Om,l,eps,M,Re,N,**kw),compute_uv=False)[-1]

def solve_xi(Om,l,eps,M,Re,seed,N=90,iters=60,tol=1e-10,**kw):
    """Newton on smallest-singular-value via det-free: minimize sigmin using complex Newton on
       a scalar dispersion function f(xi)=log det(A) is unstable; use sigmin gradient descent + polish."""
    xi=complex(seed); h=1e-6
    for _ in range(iters):
        f0=sigmin(xi,Om,l,eps,M,Re,N,**kw)
        gx=(sigmin(xi+h,Om,l,eps,M,Re,N,**kw)-sigmin(xi-h,Om,l,eps,M,Re,N,**kw))/(2*h)
        gy=(sigmin(xi+1j*h,Om,l,eps,M,Re,N,**kw)-sigmin(xi-1j*h,Om,l,eps,M,Re,N,**kw))/(2*h)
        grad=gx+1j*gy   # d(sigmin)/d(conj?) approx; use as descent direction
        if abs(grad)<1e-14: break
        step=f0*np.conj(grad)/abs(grad)**2
        xi=xi-step
        if abs(step)<tol: break
    return xi, sigmin(xi,Om,l,eps,M,Re,N,**kw)

if __name__=="__main__":
    import pb_solver as P
    l,eps,Om=3.0,0.25,1.3
    # validation below M_crit: large Re viscous -> inviscid real root
    M=0.1
    r=sorted(P.pb_roots(Om,l,eps,M,N=200)); seed=min(r,key=lambda z:abs(z-np.sqrt(Om)))
    for Re in [2e3,1e4]:
        xi,res=solve_xi(Om,l,eps,M,Re,seed,N=90)
        print(f"Om={Om} M={M} Re={Re:.0e}: viscous xi={xi.real:.4f}{xi.imag:+.5f}i (resid {res:.1e}); inviscid {seed:.4f}")
