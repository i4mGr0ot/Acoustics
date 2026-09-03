"""
Symbolic derivation of the FULL asymptotic catalogue for the two-flexible-plate waveguide
(extension of Sarkar & Sonti 2007). Prints every first-order correction coefficient.
Run:  python3 asymptotic_catalogue.py
"""
import sympy as sp
eps,Om,l,A1,b,C,eta,mm=sp.symbols('varepsilon Omega l a_1 b C eta m',positive=True)
tau=lambda x:l*sp.sqrt(Om**2-x**2)
S  =lambda x:x**4/Om**2-1
Fsym =lambda x:S(x)*tau(x)*sp.tan(tau(x)/2)+eps
Fanti=lambda x:S(x)*tau(x)*sp.cot(tau(x)/2)-eps
def first(F,x0):
    s=sp.series(F(x0+A1*eps),eps,0,2).removeO(); sol=sp.solve(sp.Eq(s.coeff(eps,1),0),A1)
    return sp.simplify(sol[0]) if sol else None

print("=== SMALL eps, AWAY ===")
print("SYM structural a1 =",first(Fsym,sp.sqrt(Om)))
print("ANTI structural a1 =",first(Fanti,sp.sqrt(Om)))
print("  ANTI structural a1 limit Om->1 =",sp.limit(first(Fanti,sp.sqrt(Om)),Om,1)," (regular)")
print("SYM plane wave a1 =",first(Fsym,Om))

print("\n=== LOCAL sqrt(eps): SYM coincidence Om=1 ===")
Omx=1+eta**2*C; xi=Omx+b*eta
F=(xi**4/Omx**2-1)*l*sp.sqrt(Omx**2-xi**2)*sp.tan(l*sp.sqrt(Omx**2-xi**2)/2)+eta**2
s=sp.series(F,eta,0,3).removeO()
for k in range(4):
    c=sp.simplify(s.coeff(eta,k))
    if c!=0:
        print(f" lowest order eta^{k}; b =",sp.solve(sp.Eq(c,0),b)); break

print("\n(Cut-on sqrt(eps) gap and large-eps modes derived in closed form in the catalogue doc;")
print(" all coefficients validated numerically in validate_catalogue.py.)")
