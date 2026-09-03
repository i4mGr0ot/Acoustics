import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import dispersion_two_plates as M
l=3.0; eps=0.15
Omv=np.linspace(0.03,2.6,360)

# exact roots (faint reference)
ds=M.sweep(M.res_sym,l,eps,Omv); da=M.sweep(M.res_asym,l,eps,Omv)

# asymptotic closed forms
def taus(O): return l*np.sqrt(complex(O**2-O))
def sym_struct(O):  t=taus(O); return np.sqrt(O)+(-1/(4*l*np.sqrt(complex(O-1))*np.tan(t/2))).real*eps
def anti_struct(O): t=taus(O); return np.sqrt(O)+(np.tan(t/2)/(4*l*np.sqrt(complex(O-1)))).real*eps
def sym_plane(O):   return O+eps/(O*l**2*(O**2-1))
def cuton(O,mm):    xi0=np.sqrt(O**2-(mm*np.pi/l)**2); S0=xi0**4/O**2-1; return xi0+2*eps/(S0*l**2*xi0)

fig,axes=plt.subplots(1,2,figsize=(13,5.6),sharey=True)

# ---- SYMMETRIC panel ----
ax=axes[0]
ax.plot(ds[:,0],ds[:,1],'.',color='0.7',ms=3,label='exact (numerical)')
# structural away (avoid Om~1 and cut-on Oc~2.65)
O=Omv[(Omv>0.15)&(Omv<0.9)];  ax.plot(O,[sym_struct(x) for x in O],'r-',lw=2,label='struct (Eq S1)')
O=Omv[(Omv>1.15)&(Omv<2.45)]; ax.plot(O,[sym_struct(x) for x in O],'r-',lw=2)
# plane wave away
O=Omv[(Omv>0.1)&(Omv<0.85)];  ax.plot(O,[sym_plane(x) for x in O],'b-',lw=2,label='plane wave (Eq S2)')
O=Omv[(Omv>1.2)&(Omv<2.55)];  ax.plot(O,[sym_plane(x) for x in O],'b-',lw=2)
# coincidence local
ax.plot([1,1],[1-np.sqrt(eps)/(2*l),1+np.sqrt(eps)/(2*l)],'gs',ms=7,label='coincidence ±√ε/2l (Eq L1)')
# even cut-on m=2 (Oc where struct meets it ~2.65)
O=Omv[Omv>2.2]; ax.plot(O,[cuton(x,2) for x in O],'m-',lw=2,label='cut-on m=2 (Eq C)')
ax.set_title(f'SYMMETRIC family  (l={l}, ε={eps})'); ax.set_xlabel('Ω'); ax.set_ylabel('ξ')
ax.set_xlim(0,2.6); ax.set_ylim(0,2.6); ax.legend(fontsize=8,loc='upper left')

# ---- ANTISYMMETRIC panel ----
ax=axes[1]
ax.plot(da[:,0],da[:,1],'.',color='0.7',ms=3,label='exact (numerical)')
O=Omv[(Omv>0.1)&(Omv<1.55)]; ax.plot(O,[anti_struct(x) for x in O],'r-',lw=2,label='struct (Eq A1, regular at Ω=1)')
O=Omv[(Omv>1.75)&(Omv<2.55)];ax.plot(O,[anti_struct(x) for x in O],'r-',lw=2)
# odd cut-on m=1 (cuts on at Om=pi/l=1.047)
O=Omv[(Omv>1.15)&(Omv<1.55)]; ax.plot(O,[cuton(x,1) for x in O],'m-',lw=2,label='cut-on m=1 (Eq C)')
O=Omv[(Omv>1.75)&(Omv<2.55)]; ax.plot(O,[cuton(x,1) for x in O],'m-',lw=2)
# antisym cut-on crossing Oc (struct meets m=1) ~1.66
Oc=np.roots([1,-1,-(np.pi/l)**2]); Oc=Oc[Oc>0][0].real
ax.plot([Oc,Oc],[np.sqrt(Oc)-np.sqrt(2*eps)/(2*l),np.sqrt(Oc)+np.sqrt(2*eps)/(2*l)],'gs',ms=7,label='cut-on gap ±√(2ε)/2l (Eq L2)')
ax.set_title(f'ANTISYMMETRIC family  (l={l}, ε={eps})'); ax.set_xlabel('Ω')
ax.set_xlim(0,2.6); ax.set_ylim(0,2.6); ax.legend(fontsize=8,loc='upper left')

fig.suptitle('Fig.4  Full asymptotic catalogue (lines) vs exact numerical roots (grey) — two flexible plates',fontsize=12)
fig.tight_layout(); fig.savefig('fig4_catalogue.png',dpi=140)
print("saved fig4_catalogue.png")
