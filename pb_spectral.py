"""
INDEPENDENT validation of the Pridmore-Brown two-flexible-plate dispersion problem
by global Chebyshev spectral collocation (cf. the shooting solver pb_solver.py).

Nondim PB on eta in [0,1], plate-admittance Robin BCs:
   p'' - (2 M xi f'/shat) p' + l^2 (shat^2 - xi^2) p = 0,  shat=Om - M xi f, f=4 eta(1-eta)
   p'(0) + (eps/S1) p(0)=0 ,  p'(1) - (eps/S2) p(1)=0.
xi enters nonlinearly; we find xi where the collocation matrix A(xi) is singular
(smallest singular value -> 0), by scanning + bisection on a sign-tracked indicator.
"""
import numpy as np

def cheb(N):
    """Trefethen Chebyshev diff matrix on [-1,1], returns D,x (x[0]=1..x[N]=-1)."""
    if N==0: return np.array([[0.0]]),np.array([1.0])
    x=np.cos(np.pi*np.arange(N+1)/N)
    c=np.hstack([2.,np.ones(N-1),2.])*(-1)**np.arange(N+1)
    X=np.tile(x,(N+1,1)).T
    dX=X-X.T
    D=(np.outer(c,1./c))/(dX+np.eye(N+1))
    D-=np.diag(D.sum(axis=1))
    return D,x

def build_A(xi,Om,l,eps,M,rD,rm,D1,D2,eta):
    f=4*eta*(1-eta); fp=4*(1-2*eta)
    shat=Om-M*xi*f
    N=len(eta)-1
    # interior operator
    A=D2 + (np.diag(2*M*xi*fp/shat))@D1 + np.diag(l**2*(shat**2-xi**2))  # +shear (sign-corrected)
    A=A.astype(complex)
    S1=xi**4/Om**2-1.0; S2=rD*xi**4/Om**2-rm
    # eta ordered: eta[0]=1 (top wall y/a=1), eta[N]=0 (bottom wall). D1 acts in eta.
    # BC at eta=1 (row 0): p'(1) - (eps/S2) p(1) = 0
    A[0,:]=D1[0,:]; A[0,0]-=eps/S2
    # BC at eta=0 (row N): p'(0) + (eps/S1) p(0) = 0
    A[N,:]=D1[N,:]; A[N,N]+=eps/S1
    return A

def sigmin(xi,Om,l,eps,M,rD,rm,D1,D2,eta):
    A=build_A(xi,Om,l,eps,M,rD,rm,D1,D2,eta)
    return np.linalg.svd(A,compute_uv=False)[-1]

def spectral_roots(Om,l,eps,M,rD=1.0,rm=1.0,N=48,xmax=3.0,Ns=900):
    D1,x=cheb(N); eta=(x+1)/2.0; D1=2*D1; D2=D1@D1
    xs=np.linspace(0.05,xmax,Ns)
    sv=np.array([sigmin(X,Om,l,eps,M,rD,rm,D1,D2,eta) for X in xs])
    roots=[]
    for i in range(1,len(xs)-1):
        if sv[i]<sv[i-1] and sv[i]<sv[i+1]:
            lo,hi=xs[i-1],xs[i+1]
            for _ in range(50):
                m1=lo+0.382*(hi-lo); m2=lo+0.618*(hi-lo)
                if sigmin(m1,Om,l,eps,M,rD,rm,D1,D2,eta)<sigmin(m2,Om,l,eps,M,rD,rm,D1,D2,eta): hi=m2
                else: lo=m1
            xr=0.5*(lo+hi)
            if sigmin(xr,Om,l,eps,M,rD,rm,D1,D2,eta)<1e-4: roots.append(xr)
    out=[]
    for r in sorted(roots):
        if not out or abs(r-out[-1])>2e-3: out.append(r)
    return out

if __name__=="__main__":
    import pb_solver as SH, num_solver as NS
    l,eps=3.0,0.25
    print("INDEPENDENT spectral (Chebyshev) vs shooting vs analytic")
    for Om in [1.3,1.8]:
        for M in [0.0,0.15]:
            sp=sorted(spectral_roots(Om,l,eps,M))
            sh=sorted(SH.pb_roots(Om,l,eps,M,N=200))
            line=f"Om={Om} M={M}: spectral={[round(x,4) for x in sp]}"
            print(line)
            print(f"           shooting={[round(x,4) for x in sh]}")
            if M==0:
                an=sorted(NS.roots(NS.res_general,Om,l,eps))
                print(f"           analytic={[round(x,4) for x in an]}")
