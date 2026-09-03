import numpy as np, pickle, os, sys, warnings
warnings.filterwarnings('ignore')
import pb_solver as P
l,eps=3.0,0.25
M=float(sys.argv[1])
Omv=np.linspace(0.2,2.5,26)
pts=[]
for Om in Omv:
    r=P.pb_roots(Om,l,eps,M,N=160)
    for x in r:
        if 0.05<x<2.05: pts.append((Om,x))
d={} 
if os.path.exists('mapS.pkl'): d=pickle.load(open('mapS.pkl','rb'))
d[M]=np.array(pts); pickle.dump(d,open('mapS.pkl','wb'))
print(f"M={M}: {len(pts)} physical points")
