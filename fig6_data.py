import numpy as np, pb_solver as P
l,eps=3.0,0.25
def cont_track(Om,xi_start,Mvals,N=140):
    """continuation: follow branch from xi_start across Mvals (assumed ordered from 0)."""
    out=[]; prev=xi_start
    for M in Mvals:
        r=P.pb_roots(Om,l,eps,M,N=N)
        if r:
            x=min(r,key=lambda z:abs(z-prev)); 
            if abs(x-prev)>0.25: x=np.nan      # reject jumps (branch lost)
            else: prev=x
        else: x=np.nan
        out.append(x)
    return np.array(out)

# A: plane-wave branch Om=1.3, continuation outward from M=0
Om=1.3; xiA=1.3267; h=0.04
fp=P.pb_roots(Om,l,eps,h,N=140); fm=P.pb_roots(Om,l,eps,-h,N=140)
fp=min(fp,key=lambda z:abs(z-xiA)); fm=min(fm,key=lambda z:abs(z-xiA))
a1=(fp-fm)/(2*h); a2=(fp-2*xiA+fm)/h**2
Mpos=np.linspace(0,0.3,7); Mneg=np.linspace(0,-0.3,7)
fpos=cont_track(Om,xiA,Mpos); fneg=cont_track(Om,xiA,Mneg)
MsA=np.concatenate([Mneg[::-1],Mpos[1:]]); fullA=np.concatenate([fneg[::-1],fpos[1:]])

# B: M_crit(Om)
OmB=np.linspace(0.3,2.2,12); McB=[]
for O in OmB:
    r=P.pb_roots(O,l,eps,0.0,N=140); McB.append(O/min(r,key=lambda z:abs(z-np.sqrt(O))))
McB=np.array(McB)

# C: structural branch Om=0.7 continuation upward toward M_crit
Om=0.7; r0=P.pb_roots(Om,l,eps,0.0,N=140); xiC=min(r0,key=lambda z:abs(z-np.sqrt(Om))); Mcrit=Om/xiC
fp=P.pb_roots(Om,l,eps,0.04,N=140); fm=P.pb_roots(Om,l,eps,-0.04,N=140)
fp=min(fp,key=lambda z:abs(z-xiC)); fm=min(fm,key=lambda z:abs(z-xiC))
c1=(fp-fm)/0.08; c2=(fp-2*xiC+fm)/0.04**2
MsC=np.linspace(0,0.92*Mcrit,10); fullC=cont_track(Om,xiC,MsC)
np.savez('fig6.npz',MsA=MsA,fullA=fullA,xiA=xiA,a1=a1,a2=a2,OmB=OmB,McB=McB,
         MsC=MsC,fullC=fullC,xiC=xiC,c1=c1,c2=c2,Mcrit=Mcrit)
print("saved. A:a1=%.3f a2=%.3f  C:xiC=%.4f Mcrit=%.3f c1=%.3f c2=%.3f"%(a1,a2,xiC,Mcrit,c1,c2))
print("fullA=",np.round(fullA,3)); print("fullC=",np.round(fullC,3))
