import numpy as np, pickle, warnings
warnings.filterwarnings('ignore')
import pb_solver as P
from pb_asym import xi1_closed
l,eps=3.0,0.25
def full_track(Om,xi0,Mvals):
    out=[]; prev=xi0
    for M in Mvals:
        r=P.pb_roots(Om,l,eps,M,N=200)
        if r:
            x=min(r,key=lambda z:abs(z-prev))
            if abs(x-prev)<0.3: prev=x; out.append(x)
            else: out.append(np.nan)
        else: out.append(np.nan)
    return np.array(out)
def coeffs(Om,xi0):
    h=0.04
    rp=P.pb_roots(Om,l,eps,h,N=200); rm=P.pb_roots(Om,l,eps,-h,N=200)
    xp=min(rp,key=lambda z:abs(z-xi0)); xm=min(rm,key=lambda z:abs(z-xi0))
    a1=(xp-xm)/(2*h); a2=(xp-2*xi0+xm)/h**2; return a1,a2
data={}
for tag,Om,seedguess in [('plane Ω=1.3',1.3,1.3267),('struct Ω=1.8',1.8,1.3469)]:
    r0=P.pb_roots(Om,l,eps,0.0,N=240); xi0=min(r0,key=lambda z:abs(z-seedguess))
    a1,a2=coeffs(Om,xi0)
    Mv=np.linspace(0,0.5,11); full=full_track(Om,xi0,Mv)
    oM=xi0+a1*Mv; oM2=xi0+a1*Mv+a2*Mv**2
    data[tag]=dict(Mv=Mv,full=full,oM=oM,oM2=oM2,xi0=xi0,a1=a1,a2=a2)
    eM=np.nanmax(np.abs(oM-full)/np.abs(full))*100; eM2=np.nanmax(np.abs(oM2-full)/np.abs(full))*100
    print(f"{tag}: xi0={xi0:.4f} a1={a1:.3f} a2={a2:.3f}  max%err O(M)={eM:.1f} O(M2)={eM2:.1f}")
pickle.dump(data,open('validity.pkl','wb'))
