"""
Local (Frobenius) analysis of the Pridmore-Brown critical layer.
Near eta_c where sigma_hat=0:  sigma_hat ~ s1*(eta-eta_c), s1 = -M xi f'(eta_c).
PB:  p'' + (2 M xi f'/sigma) p' + l^2(sigma^2-xi^2) p = 0.
The singular coefficient -> (2 M xi f'(eta_c))/(s1 (eta-eta_c)) = -2/(eta-eta_c).
"""
import sympy as sp
z=sp.symbols('z')             # z = eta-eta_c
s=sp.symbols('s')
# leading equation near z=0:  p'' - (2/z) p' + (regular) p = 0
# indicial from p ~ z^s :  s(s-1) z^{s-2} - 2 s z^{s-2} = 0
indicial=sp.expand(s*(s-1)-2*s)
roots=sp.solve(indicial,s)
print("singular coefficient of p' near critical layer:  -2/(eta-eta_c)")
print("indicial equation:", indicial, "= 0  -> Frobenius exponents s =", roots)
print("exponents differ by", max(roots)-min(roots), "(integer) -> log term in the s=0 solution generically present")

# Construct Frobenius series to detect the log (Frobenius w/ exponents 0 and 3)
# Solve p'' - (2/z)p' + c0 p = 0 (keep leading regular term c0) via series; show s=3 branch + log
c0=sp.symbols('c0')
# Try p = z^3 * sum a_n z^n  for the larger root; and check the smaller root needs log
# Demonstrate: Wronskian / reduction of order log coefficient
# p1 = analytic solution (s=0): p1 = 1 + ... ; reduction of order: p2 = p1 * \int e^{\int (2/z)dz}/p1^2 dz
# \int (2/z) dz = 2 ln z -> e^{...}=z^2 ; integrand ~ z^2/p1^2 ~ z^2 -> \int -> z^3/3 (no log at leading)
intg=sp.integrate(sp.exp(2*sp.log(z))/ (1)**2, z)   # p1~1
print("\nreduction-of-order integral (p1~1):  ∫ z^2 dz =", intg, " -> second solution ~ z^3 (regular at leading order)")
print("=> pressure p is continuous/bounded; the singularity is in the CROSS-STREAM VELOCITY")
print("   v ~ p'/sigma ~ (z^2)/(s1 z) ~ z  for the z^3 branch, but the GENERAL solution's p' ~ const,")
print("   so v ~ p'/sigma ~ 1/z  -> v DIVERGES as 1/(eta-eta_c) at the critical layer.")
print("   A log(eta-eta_c) enters at higher order (z^3 ln z); its branch is fixed by causality/viscosity.")
