import numpy as np, num_solver as NS

# ---------- (1) VEERING: two structural branches vs rD (rm=1), fixed Om ----------
def struct_roots(Om,l,eps,rD,rm=1.0):
    # roots of full relation near xi=sqrt(Om) (the two structural waves)
    rs=NS.roots(NS.res_general,Om,l,eps,rD=rD,rm=rm)
    s=np.sqrt(Om)
    near=sorted(rs,key=lambda r:abs(r-s))[:2]
    return sorted(near)
l,eps,Om=3.0,0.25,1.3
print("VEERING: two structural-branch wavenumbers vs r_D (Om=1.3, eps=0.25)")
print(f"{'rD':>5} {'xi_-':>8} {'xi_+':>8} {'gap':>8}")
gaps=[]
for rD in [0.6,0.8,0.9,1.0,1.1,1.2,1.4]:
    rr=struct_roots(Om,l,eps,rD)
    if len(rr)==2:
        g=rr[1]-rr[0]; gaps.append((rD,g)); print(f"{rD:>5} {rr[0]:>8.4f} {rr[1]:>8.4f} {g:>8.4f}")
# identical-plate gap formula
tau_s=l*np.sqrt(Om**2-Om)
gap_formula=eps/(2*l*np.sqrt(Om-1)*np.sin(tau_s))
print(f"identical-plate (rD=1) gap formula = {gap_formula:.4f}")

# gap scaling with eps at rD=1
print("\nGap at rD=1 vs eps (check linear scaling):")
for e in [0.05,0.1,0.2]:
    rr=struct_roots(Om,l,e,1.0); g=rr[1]-rr[0]
    print(f"  eps={e}: gap={g:.4f}  formula={e/(2*l*np.sqrt(Om-1)*np.sin(tau_s)):.4f}")

# ---------- (2) DAMPING: D->D(1+i eta) makes alpha complex -> complex xi ----------
# leading-order: xi = xi0 + (dxi/dS) * dS, with structural loss. Use plane-wave & struct asymptotics.
# Easiest: re-solve relation with complex S via Newton from real root.
def newton_complex(Om,l,eps,xi0,eta_loss,which='sym',iters=40):
    # complex dispersion fn for identical plates with loss factor eta_loss: S=(1+i eta)xi^4/Om^2 -1
    def F(xi):
        S=(1+1j*eta_loss)*xi**4/Om**2-1.0
        tau=l*np.sqrt(Om**2-xi**2+0j)
        if which=='sym':  return S*tau*np.tan(tau/2)+eps
        else:             return S*tau*np.cos(tau/2)/np.sin(tau/2)*np.sin(tau/2)/np.cos(tau/2)*0+S*tau/np.tan(tau/2)-eps
    xi=complex(xi0)
    for _ in range(iters):
        h=1e-6
        d=(F(xi+h)-F(xi-h))/(2*h)
        xi=xi-F(xi)/d
    return xi
print("\nDAMPING (sym branches, structural loss eta=0.02): complex xi = xi_r + i xi_i")
for xi0 in struct_roots(Om,l,eps,1.0):
    z=newton_complex(Om,l,eps,xi0,0.02,'sym')
    print(f"  xi0={xi0:.4f} -> {z.real:.4f} {z.imag:+.5f} i   (spatial decay rate -Im xi = {-z.imag:.5f})")

# ---------- (3) GROUP VELOCITY along branches (dOmega/dxi) ----------
print("\nGROUP VELOCITY  d(Omega)/d(xi)  along symmetric branches (finite diff):")
def branch_vg(res,Om,l,eps,xi_target,dO=0.01):
    rp=NS.roots(res,Om+dO,l,eps); rm=NS.roots(res,Om-dO,l,eps)
    xp=min(rp,key=lambda r:abs(r-xi_target)); xm=min(rm,key=lambda r:abs(r-xi_target))
    return (2*dO)/(xp-xm) if xp!=xm else np.inf
for xi0 in sorted(NS.roots(NS.res_sym,Om,l,eps)):
    vg=branch_vg(NS.res_sym,Om,l,eps,xi0)
    print(f"  xi={xi0:.4f}: vg=dOmega/dxi={vg:.3f}")
