import numpy as np, pb_solver as P
l,eps=3.0,0.25

def track(Om,xi0,Ms,shear=True):
    out=[]
    for M in Ms:
        r=P.pb_roots(Om,l,eps,M,N=240,shear=shear)
        out.append(min(r,key=lambda z:abs(z-xi0)) if r else np.nan)
    return np.array(out)

# ---- O(M^2): Richardson extraction of xi0,xi1,xi2 for branches at Om=1.3 ----
print("=== O(M^2): xi(M) = xi0 + xi1 M + xi2 M^2 ; coefficients by central differences ===")
h=0.04
for xi0 in sorted(P.pb_roots(1.3,l,eps,0.0,N=240)):
    fp=track(1.3,xi0,[h])[0]; fm=track(1.3,xi0,[-h])[0]; f0=xi0
    xi1=(fp-fm)/(2*h); xi2=(fp-2*f0+fm)/h**2
    # shear-off (convection only)
    fps=track(1.3,xi0,[h],shear=False)[0]; fms=track(1.3,xi0,[-h],shear=False)[0]
    xi1_ns=(fps-fms)/(2*h); xi2_ns=(fps-2*f0+fms)/h**2
    print(f" xi0={xi0:.4f}: xi1={xi1:+.4f} xi2={xi2:+.4f} | conv-only xi2={xi2_ns:+.4f} | shear share of xi2={(xi2-xi2_ns):+.4f}")

# accuracy gain: at M=0.25 compare O(M) vs O(M^2) to full solver (plane-wave branch)
xi0=1.3267; M=0.25
full=track(1.3,xi0,[M])[0]
fp=track(1.3,xi0,[h])[0]; fm=track(1.3,xi0,[-h])[0]
xi1=(fp-fm)/(2*h); xi2=(fp-2*xi0+fm)/h**2
oM =xi0+xi1*M; oM2=xi0+xi1*M+xi2*M**2
print(f"\nAccuracy at M={M} (plane-wave branch): full={full:.4f}  O(M)={oM:.4f} (err {abs(oM-full)/full*100:.2f}%)  O(M^2)={oM2:.4f} (err {abs(oM2-full)/full*100:.2f}%)")

# ---- Critical layer: M_crit = Om/xi (phase Mach); structural branch below coincidence ----
print("\n=== Critical layer: sigma=Om-M xi f=0 reachable when M >= M_crit = Om/xi ===")
for Om in [0.5,0.7,0.9,1.3,1.8]:
    rs=P.pb_roots(Om,l,eps,0.0,N=240)
    # structural branch nearest sqrt(Om)
    xs=min(rs,key=lambda z:abs(z-np.sqrt(Om)))
    print(f" Om={Om}: structural xi0={xs:.4f}, phase speed c_p=Om/xi={Om/xs:.4f}=M_crit  ({'subsonic, CL at M<1' if Om/xs<1 else 'supersonic, no CL below M=1'})")
