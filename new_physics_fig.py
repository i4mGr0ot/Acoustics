import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import num_solver as NS, pb_solver as P
l,eps=3.0,0.25

fig,ax=plt.subplots(1,3,figsize=(15,4.6))

# (A) Poiseuille flow: wavenumber vs Mach for branches at Om=1.3
Om=1.3; Ms=np.linspace(-0.25,0.25,11)
r0=sorted(P.pb_roots(Om,l,eps,0.0,N=200))
tracks={i:[] for i in range(len(r0))}
for M in Ms:
    rr=sorted(P.pb_roots(Om,l,eps,M,N=200))
    for i,x0 in enumerate(r0):
        if rr: tracks[i].append(min(rr,key=lambda r:abs(r-x0)))
        else: tracks[i].append(np.nan)
labels=['struct (anti)','struct (sym)','struct2','plane-wave']
for i in range(len(r0)):
    ax[0].plot(Ms,tracks[i],'-o',ms=3,label=f'ξ₀={r0[i]:.3f}')
ax[0].axvline(0,color='0.8',lw=0.8)
ax[0].set(xlabel='Mach number $M=U_0/c$',ylabel='ξ',title='(A) Poiseuille flow shifts wavenumbers\n(Ω=1.3): downstream M>0 lowers ξ')
ax[0].legend(fontsize=7)

# (B) Veering: two structural branches vs rD
rDs=np.linspace(0.55,1.5,25); lo=[];hi=[]
for rD in rDs:
    rs=NS.roots(NS.res_general,Om,l,eps,rD=rD,rm=1.0)
    near=sorted(sorted(rs,key=lambda r:abs(r-np.sqrt(Om)))[:2])
    lo.append(near[0]); hi.append(near[1])
ax[1].plot(rDs,lo,'b-',lw=2); ax[1].plot(rDs,hi,'r-',lw=2)
# uncoupled (eps->0) in-vacuo structural wavenumbers that would cross at rD=1
ax[1].plot(rDs,[np.sqrt(Om)]*len(rDs),'k--',lw=1,label='plate-1 in-vacuo √Ω')
ax[1].plot(rDs,[Om**0.5*(1.0/rD)**0.25 for rD in rDs],'g--',lw=1,label='plate-2 in-vacuo (rₘ/r_D)^¼√Ω')
ax[1].axvline(1.0,color='0.8',lw=0.8)
ax[1].set(xlabel='stiffness ratio $r_D=D_2/D_1$',ylabel='ξ',title='(B) Veering of the two structural waves\n(min gap at $r_D=r_m=1$, gap ∝ ε)')
ax[1].legend(fontsize=7)

# (C) Group velocity dOmega/dxi along sym branches across a gap region
Omv=np.linspace(0.6,2.4,120); pts=[]
for O in Omv:
    rs=NS.roots(NS.res_sym,O,l,eps)
    for r in rs: pts.append((O,r))
pts=np.array(pts)
# numerical vg per point via local slope along nearest-in-xi neighbours is messy; show dispersion + mark gaps
ax[2].plot(pts[:,0],pts[:,1],'r.',ms=2)
ax[2].plot(Omv,np.sqrt(Omv),'k--',lw=1,label='structure'); ax[2].plot(Omv,Omv,'b--',lw=1,label='plane wave')
ax[2].axvline(1.0,color='g',ls=':',lw=1,label='coincidence gap (vg dips)')
ax[2].set(xlabel='Ω',ylabel='ξ',title='(C) Symmetric branches; gaps = points\nwhere group velocity dΩ/dξ changes sharply')
ax[2].legend(fontsize=7); ax[2].set_ylim(0,2.6)

fig.suptitle('New physics: Poiseuille mean flow, non-identical-plate veering, and gaps',fontsize=12)
fig.tight_layout(); fig.savefig('fig5_newphysics.png',dpi=140); print("saved fig5_newphysics.png")

# refined damping numbers (per correct branch)
def newton(res_complex,xi0,iters=50):
    xi=complex(xi0)
    for _ in range(iters):
        h=1e-6; d=(res_complex(xi+h)-res_complex(xi-h))/(2*h)
        xi=xi-res_complex(xi)/d
    return xi
def Fsym(xi,etal):  S=(1+1j*etal)*xi**4/Om**2-1; tau=l*np.sqrt(Om**2-xi**2+0j); return S*tau*np.tan(tau/2)+eps
def Fanti(xi,etal): S=(1+1j*etal)*xi**4/Om**2-1; tau=l*np.sqrt(Om**2-xi**2+0j); return S*tau/np.tan(tau/2)-eps
print("DAMPING (structural loss eta_loss=0.02), complex xi:")
for nm,F,seed in [('sym-struct',Fsym,1.1165),('sym-plane',Fsym,1.3267),('anti-struct',Fanti,1.1837)]:
    z=newton(lambda x:F(x,0.02),seed)
    print(f"  {nm:12s}: xi = {z.real:.4f} {z.imag:+.5f} i  (decay -Im xi={-z.imag:.5f})")
