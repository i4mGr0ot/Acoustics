"""Regenerate fig1-3 (dispersion curves) with the fixed solver."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import dispersion_two_plates as M
l=3.0; Omv=np.linspace(0.03,2.5,240); st,pl=M.uncoupled(Omv,l)
def cut(ax,mult,style,lab):
    for i,n in enumerate(mult):
        val=Omv**2-(n*np.pi/l)**2
        ax.plot(Omv,np.where(val>0,np.sqrt(np.abs(val)),np.nan),style,lw=1,label=lab if i==0 else None)
# fig1
fig,ax=plt.subplots(figsize=(7,5)); ax.plot(Omv,st,'k--',lw=1,label='structure'); ax.plot(Omv,pl,'b--',lw=1,label='plane wave')
cut(ax,[1,2,3,4],'g:','rigid cut-ons'); d=M.sweep(M.res_SS,l,0.25,Omv); ax.plot(d[:,0],d[:,1],'r.',ms=2,label='coupled ε=0.25')
ax.set(xlim=(0,2.5),ylim=(0,3),xlabel='Ω',ylabel='ξ',title='Fig.1 Validation: single flexible + rigid wall'); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig('fig1_validation.png',dpi=140)
# fig2
fig,ax=plt.subplots(figsize=(7.5,5.5)); ax.plot(Omv,st,'k--',lw=1,label='structure'); ax.plot(Omv,pl,'b--',lw=1,label='plane wave')
cut(ax,[2,4,6],'g:','sym cut-ons 2nπ'); cut(ax,[1,3,5,7],'m:','anti cut-ons (2n+1)π')
ds=M.sweep(M.res_sym,l,0.25,Omv); da=M.sweep(M.res_asym,l,0.25,Omv)
ax.plot(ds[:,0],ds[:,1],'r.',ms=2.2,label='coupled SYM'); ax.plot(da[:,0],da[:,1],'.',color='darkorange',ms=2.2,label='coupled ANTI')
ax.set(xlim=(0,2.5),ylim=(0,3),xlabel='Ω',ylabel='ξ',title='Fig.2 Two identical flexible plates (l=3, ε=0.25)'); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig('fig2_two_identical.png',dpi=140)
# fig3
fig,ax=plt.subplots(figsize=(7.5,5.5)); ax.plot(Omv,st,'k--',lw=1,label='structure'); ax.plot(Omv,pl,'b--',lw=1,label='plane wave')
for rD,col,lab in [(1.0,'r','identical'),(4.0,'purple','plate 2 stiffer rD=4')]:
    d=M.sweep(M.res_general,l,0.25,Omv,rD=rD,rm=1.0)
    if len(d): ax.plot(d[:,0],d[:,1],'.',color=col,ms=2.0,label=lab)
ax.set(xlim=(0,2.5),ylim=(0,3),xlabel='Ω',ylabel='ξ',title='Fig.3 Non-identical plates (full relation)'); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig('fig3_nonidentical.png',dpi=140)
print("figs 1-3 regenerated")
