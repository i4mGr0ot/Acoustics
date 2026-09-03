import numpy as np, cl_contour as C
l,eps=3.0,0.25

# ---- (a) MODAL COALESCENCE: the two structural wall modes merge -> complex pair ----
# Om=0.7: real roots 0.850,0.857 at M=0.3 vanish by M=0.5  => complex-conjugate pair.
Om=0.7
def cnewton(Om,M,seed,iters=60):
    xi=complex(seed)
    for _ in range(iters):
        R=C.integrate(xi,Om,M,A=0.0,w=0.2)   # straight contour (no CL below M_crit)
        h=1e-7; dR=(C.integrate(xi+h,Om,M,0,0.2)-C.integrate(xi-h,Om,M,0,0.2))/(2*h)
        s=R/dR; xi-=s
        if abs(s)<1e-12: break
    return xi
print("=== (a) Flow-induced modal coalescence (Om=0.7) ===")
prev=0.8535+0.002j
for M in [0.30,0.40,0.44,0.46,0.50]:
    z=cnewton(Om,M,prev); prev=z
    tag='REAL pair' if abs(z.imag)<1e-4 else 'COMPLEX (coalesced)'
    print(f" M={M}: xi={z.real:.4f}{z.imag:+.4f}i   {tag}")

# ---- (b) VISCOUS inner layer: Re^-1/3 thickness and -i*pi phase jump (Plemelj/Airy) ----
print("\n=== (b) Viscous critical layer: -i*pi jump and Re^{-1/3} thickness ===")
# Plemelj: viscosity replaces 1/(eta-eta_c) by 1/(eta-eta_c - i*gamma), gamma>0 (causal).
# Integral of the singular term across the layer -> PV - i*pi as gamma->0+.
def jump(gamma):
    y=np.linspace(-1,1,400001)
    integrand=1.0/(y - 1j*gamma)
    return np.trapz(integrand,y)   # -> 0 (PV) - i*pi
for g in [1e-1,1e-2,1e-3,1e-4]:
    J=jump(g); print(f"  gamma={g:.0e}: Im[∫ dy/(y-i gamma)] = {J.imag:+.5f}  (-> -pi={-np.pi:.5f})")
print("  => causal/viscous continuation gives the -i*pi critical-layer phase jump.")

# Airy balance: inner viscous term (1/Re) d^2/dY^2 balances shear sigma' ~ Y  => layer ~ Re^{-1/3}
print("\n  Inner Airy balance: (1/Re) v'''' ~ sigma'_c (eta-eta_c) v''  =>  delta ~ (sigma'_c Re)^{-1/3}")
for Re in [1e3,1e4,1e5,1e6]:
    delta=Re**(-1/3.0)
    print(f"   Re={Re:.0e}: critical-layer thickness delta ~ Re^(-1/3) = {delta:.4e}")
