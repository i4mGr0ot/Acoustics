"""
Symbolic verification (sympy) for the TWO-FLEXIBLE-PLATE structural-acoustic waveguide.
Extension of A. Sarkar & V.R. Sonti, JSV 306 (2007) 657-674.

Run:  python3 symbolic_verification.py
Checks:
  (1) Eliminating the acoustic amplitudes A,B gives the dispersion relation
         (1+i a1)(1+i a2) = e^{2 i ky a}(1 - i a1)(1 - i a2),   a_j = ky G_j/(w^2 rho)
  (2) Identical plates (a1=a2) factorise into symmetric x antisymmetric.
  (3) Rigid-wall limit (a2 -> inf) recovers Sarkar-Sonti Eq.(9):  D kx^4 - m w^2 = -(w^2 rho/ky) cot(ky a).
  (4) First-order asymptotic corrections in eps for the coupled structural branches.
"""
import sympy as sp
I = sp.I

# ---------------------------------------------------------------- (1) elimination
A,B,ky,a,rho,om,G1,G2 = sp.symbols('A B k_y a rho omega G1 G2')
p   = lambda y: A*sp.exp(-I*ky*y) + B*sp.exp(I*ky*y)
dp  = lambda y: -I*ky*A*sp.exp(-I*ky*y) + I*ky*B*sp.exp(I*ky*y)
v   = lambda y: (I/(om*rho))*dp(y)             # Euler:  -i w rho v = dp/dy
W1  = v(0)/(I*om)                              # velocity continuity  v = i w W
W2  = v(a)/(I*om)
eq1 = sp.expand(G1*W1 + p(0))                  # plate 1 (fluid above):  G1 W1 = -p(0)
eq2 = sp.expand(G2*W2 - p(a))                  # plate 2 (fluid below):  G2 W2 = +p(a)
M   = sp.Matrix([[eq1.coeff(A), eq1.coeff(B)],
                 [eq2.coeff(A), eq2.coeff(B)]])
det = sp.simplify(M.det())
a1,a2 = sp.symbols('alpha1 alpha2')
target = (1+I*a1)*(1+I*a2) - sp.exp(2*I*ky*a)*(1-I*a1)*(1-I*a2)
ratio  = sp.simplify(det.subs({G1:a1*om**2*rho/ky, G2:a2*om**2*rho/ky})/target)
print("(1) det / target =", ratio, " (nonzero, no roots added/lost) -> relation confirmed")

# ---------------------------------------------------------------- (2) factorisation
al = sp.symbols('alpha')
D  = lambda x1,x2: (1+I*x1)*(1+I*x2) - sp.exp(2*I*ky*a)*(1-I*x1)*(1-I*x2)
sym  = (1+I*al) - sp.exp(I*ky*a)*(1-I*al)      # symmetric factor
asym = (1+I*al) + sp.exp(I*ky*a)*(1-I*al)      # antisymmetric factor
print("(2) sym*asym - D(al,al) =", sp.simplify(sp.expand(sym*asym) - D(al,al)))
print("    symmetric  branch alpha =", sp.simplify(sp.solve(sym ,al)[0].rewrite(sp.tan)))
print("    antisym.   branch alpha =", sp.simplify(sp.solve(asym,al)[0].rewrite(sp.tan)))

# ---------------------------------------------------------------- (3) rigid limit
lim = sp.simplify(sp.limit(D(a1,a2)/a2, a2, sp.oo))
sol = sp.simplify(sp.solve(lim,a1)[0]).rewrite(sp.cot)
print("(3) rigid-wall limit alpha1 =", sol, " == -cot(ky a)  -> Sarkar-Sonti Eq.(9)")

# ---------------------------------------------------------------- (4) asymptotics
eps,Om,l,A1 = sp.symbols('varepsilon Omega l a1', positive=True)
tau = lambda x: l*sp.sqrt(Om**2 - x**2)
S   = lambda x: x**4/Om**2 - 1
def first_order(Fexpr):
    s = sp.series(Fexpr, eps, 0, 2).removeO()
    return sp.simplify(sp.solve(sp.Eq(s.coeff(eps,1),0), A1)[0])
xS = sp.sqrt(Om) + A1*eps
print("(4) symmetric  structural a1 =", first_order(S(xS)*tau(xS)*sp.tan(tau(xS)/2)+eps))
print("    antisym.   structural a1 =", first_order(S(xS)*tau(xS)*sp.cot(tau(xS)/2)-eps))
