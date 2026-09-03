# Coupled Dispersion of a 1-D Structural–Acoustic Waveguide Bounded by **Two Flexible Plates**

### An asymptotic extension of Sarkar & Sonti (JSV 306, 2007)

*Working document — derivation, literature/gap review, asymptotic analysis, and numerically validated results. Every result below has been verified symbolically (`symbolic_verification.py`) and numerically (`dispersion_two_plates.py`).*

---

## 0. Executive summary

Sarkar & Sonti (2007) derived analytical (asymptotic) expressions for the fluid–structure
coupled wavenumbers of a fluid column bounded **below by a flexible plate and above by a rigid
cover plate**. This document generalises that work to the case where **both bounding walls are
flexible plates** — in general non-identical. The main results are:

1. **General dispersion relation.** Eliminating the two acoustic amplitudes gives

$$\boxed{\;(1+\mathrm{i}\alpha_1)(1+\mathrm{i}\alpha_2)=e^{\,2\mathrm{i}k_y a}\,(1-\mathrm{i}\alpha_1)(1-\mathrm{i}\alpha_2)\;}\qquad
\alpha_j=\frac{k_y\,G_j}{\omega^2\rho},\quad G_j=D_j k_x^4-m_j\omega^2 .$$

   For real (propagating) transverse wavenumber this collapses to the compact real form

$$\boxed{\;\tan(k_y a)=\frac{\alpha_1+\alpha_2}{1-\alpha_1\alpha_2}\;}$$

2. **Symmetry decoupling.** For *identical* plates the relation factorises **exactly** into a
   **symmetric** branch ($v=0$ at the mid-plane, $\alpha=-\cot(k_y a/2)$) and an
   **antisymmetric** branch ($p=0$ at the mid-plane, $\alpha=+\tan(k_y a/2)$).

3. **Reduction / validation.** The rigid-wall limit ($\alpha_2\to\infty$) recovers Sarkar & Sonti
   Eq. (9) identically; the symmetric branch is *exactly* the Sarkar–Sonti single-plate problem with
   the gap and the fluid-loading parameter both halved. Numerically the two agree to 4 decimals at
   every frequency tested. This is exactly Fahy's symmetric-mode model — now joined by the
   previously-untracked antisymmetric family.

4. **New physics.** The single rigid-duct cut-on ladder splits into an **even (symmetric)** family
   ($k_y a = 2n\pi$) and an **odd (antisymmetric)** family ($k_y a = (2n{+}1)\pi$); the plane wave
   belongs to the symmetric family; the antisymmetric family has **no plane wave** (its lowest
   member is a pressure-release cut-on). For non-identical plates the two symmetry families
   **re-couple**, producing avoided crossings (veering) between the in-phase and out-of-phase
   structural waves — a feature absent from both Fahy and Sarkar–Sonti.

---

## 1. Literature & gap review

### 1.1 The line of work this extends

The canonical reference is **Fahy, *Sound and Structural Vibration* (1989)**, which set up the
problem of a fluid layer between two identical flexible plates and analysed the **symmetric modes
only**, presenting *numerical* dispersion solutions that could not be traced back to their
uncoupled (in-vacuo / rigid-duct) parents.

**Sarkar & Sonti (2007)** replaced one of Fahy's flexible plates with a rigid cover plate — which,
under Fahy's symmetry assumption, is *equivalent to the symmetric half-space* — and made the
decisive methodological advance: instead of solving numerically, they wrote the coupled dispersion
relation as the *in-vacuo structural dispersion relation plus a fluid-loading correction*, and built
an **asymptotic series in a fluid-loading parameter** $\varepsilon=\rho a/m$. Setting
$\varepsilon=0$ recovers the uncoupled wavenumbers; expanding for small (and separately large)
$\varepsilon$ gives closed-form expressions that let each coupled branch be *continuously tracked*
from its uncoupled parent. They showed that (i) every branch transits from a **rigid-walled** to a
**pressure-release** wavenumber as $\varepsilon$ grows, and (ii) wherever two uncoupled curves
intersect, coupling opens a **gap** (no crossings).

Related numerical/experimental studies of the same family of problems include Cabelli (1985, square
duct, one flexible wall), Huang *et al.* (2000) and Choi & Kim (2002) (membrane ducts), Ko (1994,
mean flow), and Fuller & Fahy (1982) / Pavić (1990) (fluid-filled cylindrical shells). The
asymptotic / matched-expansion philosophy traces to Morse & Ingard (1968) and to Crighton's fluid-
loading work (1989, 1992), but those treat **infinite (external) fluid domains**, where branch cuts
appear; for the **bounded** column considered here there are no branch cuts but instead **cut-on**
phenomena.

### 1.2 The gap

| | Walls | Modes tracked | Method | Coupling of symmetry families |
|---|---|---|---|---|
| Fahy (1989) | 2 flexible (identical) | symmetric only | numerical | — |
| Sarkar–Sonti (2007) | 1 flexible + 1 rigid | = symmetric only | **asymptotic** | — |
| **This work** | **2 flexible (general)** | **symmetric + antisymmetric** | **asymptotic** | **yes (non-identical)** |

The specific, defensible novelty is therefore: a **unified asymptotic treatment of the full
two-flexible-plate waveguide**, including (a) the antisymmetric family that Fahy excluded and
Sarkar–Sonti could not see, and (b) the **non-identical-plate case**, where the symmetric and
antisymmetric problems re-couple and the asymptotic machinery must be applied to the full
$2\times2$ relation rather than to two scalar relations. Note also that the second uploaded paper —
Morab, Sharma & Murallidharan (*J. Comput. Phys.* 514, 2024) — is a *computational* flow-induced
acoustics paper; it is unrelated in method but is a useful pointer if a later goal is to validate the
analytical dispersion results against a finite-volume/finite-difference solver.

---

## 2. Model and governing equations

### 2.1 Geometry

A two-dimensional acoustic fluid occupies the strip $0\le y\le a$ and is unbounded in $x$ (all
quantities independent of $z$). A thin elastic plate lies at $y=0$ (plate 1) and a second thin
elastic plate lies at $y=a$ (plate 2). Time dependence $e^{\mathrm{i}\omega t}$ is suppressed
throughout. Plate $j$ has flexural rigidity $D_j=E_jI_j$ with
$I_j=h_j^3/[12(1-\nu_j^2)]$, mass per unit area $m_j=\rho_{p,j}h_j$, and transverse displacement
$w_j(x)$ measured **positive in $+y$**.

> **Reduction map.** Plate 2 rigid $\Rightarrow$ Sarkar–Sonti (2007). Plate 2 rigid *and*
> $D_1\to\infty$ $\Rightarrow$ a rigid duct. $D_j\to0$ $\Rightarrow$ pressure-release walls.

### 2.2 Acoustic field

The pressure satisfies the 2-D Helmholtz equation $p_{xx}+p_{yy}+k^2p=0$, $k=\omega/c$, with the
general $x$-travelling solution

$$p(x,y)=\big(A\,e^{-\mathrm{i}k_y y}+B\,e^{+\mathrm{i}k_y y}\big)\,e^{-\mathrm{i}k_x x},
\qquad k_x^2+k_y^2=k^2 .$$

$k_x$ is the (sought) axial wavenumber and $k_y=\sqrt{k^2-k_x^2}$ the transverse wavenumber. The
linearised $y$-momentum (Euler) equation $-\mathrm{i}\omega\rho\,v=\partial p/\partial y$ gives the
fluid transverse velocity

$$v(x,y)=\frac{k_y}{\omega\rho}\big(A\,e^{-\mathrm{i}k_y y}-B\,e^{+\mathrm{i}k_y y}\big)e^{-\mathrm{i}k_x x}.$$

### 2.3 Plate equations and interface conditions

Each plate obeys the thin-plate (Euler–Bernoulli/Kirchhoff) equation
$D_j\,\partial_x^4 w_j-m_j\omega^2 w_j=q_j$, where $q_j$ is the net transverse load per area in $+y$.
With $w_j\propto e^{-\mathrm{i}k_x x}$ the operator becomes $G_j\equiv D_j k_x^4-m_j\omega^2$.

* **Velocity continuity** (fluid cannot separate from plate): $v=\mathrm{i}\omega w_j$ at each wall.
* **Loading sense.** The fluid is **above** plate 1, so it pushes it in $-y$: $q_1=-p(x,0)$. The
  fluid is **below** plate 2, so it pushes it in $+y$: $q_2=+p(x,a)$.

Substituting the field expressions:

$$G_1 W_1=-(A+B),\qquad G_2 W_2 = A e^{-\mathrm{i}k_y a}+B e^{+\mathrm{i}k_y a},$$
$$W_1=-\frac{\mathrm{i}k_y}{\omega^2\rho}(A-B),\qquad
W_2=-\frac{\mathrm{i}k_y}{\omega^2\rho}\big(A e^{-\mathrm{i}k_y a}-B e^{+\mathrm{i}k_y a}\big).$$

---

## 3. The coupled dispersion relation

### 3.1 General (non-identical plates)

Writing the two interface equations as a homogeneous linear system in $(A,B)$ and setting the
determinant to zero (so a non-trivial acoustic field exists) gives — after introducing the
dimensionless **fluid-loading admittance of each plate**

$$\alpha_j\equiv\frac{k_y\,G_j}{\omega^2\rho}=\frac{k_y\,(D_j k_x^4-m_j\omega^2)}{\omega^2\rho}\;:$$

$$(1+\mathrm{i}\alpha_1)(1+\mathrm{i}\alpha_2)=e^{\,2\mathrm{i}k_y a}\,(1-\mathrm{i}\alpha_1)(1-\mathrm{i}\alpha_2).\tag{$\star$}$$

Because $(1+\mathrm{i}\alpha)/(1-\mathrm{i}\alpha)=e^{\,2\mathrm{i}\arctan\alpha}$, ($\star$) is
equivalent to the **phase-matching condition**

$$k_y a=\arctan\alpha_1+\arctan\alpha_2+n\pi
\quad\Longleftrightarrow\quad
\tan(k_y a)=\frac{\alpha_1+\alpha_2}{1-\alpha_1\alpha_2}.\tag{$\star\star$}$$

This is the cleanest working form: it is a single **real** transcendental equation in $k_x$ for any
combination of plate properties. *(Verified in `symbolic_verification.py`, check 1: the determinant
equals the left-minus-right of ($\star$) times $e^{-\mathrm{i}k_y a}$ — a nonzero factor, so no roots
are created or destroyed.)*

### 3.2 Rigid-wall limit → Sarkar & Sonti (2007)

Letting plate 2 become rigid, $D_2\to\infty\Rightarrow\alpha_2\to\infty$, ($\star\star$) gives
$\tan(k_y a)=-1/\alpha_1$, i.e. $\alpha_1=-\cot(k_y a)$, which is

$$D_1 k_x^4-m_1\omega^2=-\frac{\omega^2\rho}{k_y}\cot(k_y a),$$

**identical to Sarkar–Sonti Eq. (9).** *(Verified symbolically, check 3.)*

### 3.3 Identical plates → symmetric / antisymmetric decoupling

If the plates are identical ($\alpha_1=\alpha_2=\alpha$, $G_1=G_2=G$) then ($\star$) is a perfect
square and factorises:

$$\Big[(1+\mathrm{i}\alpha)-e^{\mathrm{i}k_y a}(1-\mathrm{i}\alpha)\Big]\cdot
   \Big[(1+\mathrm{i}\alpha)+e^{\mathrm{i}k_y a}(1-\mathrm{i}\alpha)\Big]=0,$$

giving two **decoupled** scalar dispersion relations *(verified, check 2)*:

| Family | Mid-plane condition | Plate motion | Dispersion relation |
|---|---|---|---|
| **Symmetric** | $v=0$ at $y=a/2$ (rigid-like) | breathing, $W_2=-W_1$ | $\;D k_x^4-m\omega^2=-\dfrac{\omega^2\rho}{k_y}\cot\!\dfrac{k_y a}{2}\;$ |
| **Antisymmetric** | $p=0$ at $y=a/2$ (release) | sloshing, $W_2=+W_1$ | $\;D k_x^4-m\omega^2=+\dfrac{\omega^2\rho}{k_y}\tan\!\dfrac{k_y a}{2}\;$ |

The **symmetric** relation is Sarkar–Sonti's *with $a\to a/2$* — i.e. a flexible plate over a rigid
wall at the mid-plane. This is Fahy's symmetric mode, now reproduced analytically. The
**antisymmetric** relation is the new branch.

---

## 4. Non-dimensionalisation

Following Sarkar & Sonti, scale on the **coincidence** condition of plate 1. With
$\omega_c=c^2\sqrt{m_1/(E_1I_1)}$, $k_c=\omega_c/c$:

$$\Omega=\frac{\omega}{\omega_c}=\frac{k}{k_c},\quad
\xi=\frac{k_x}{k_c},\quad
l=k_c a,\quad
\varepsilon=\frac{\rho a}{m_1},\quad
\tau\equiv k_y a=l\sqrt{\Omega^2-\xi^2}.$$

The in-vacuo plate-1 operator scales as $G_1=m_1\omega^2\big(\xi^4/\Omega^2-1\big)$, so

$$\alpha_1=\frac{1}{\varepsilon}\Big(\frac{\xi^4}{\Omega^2}-1\Big)\tau,\qquad
\alpha_2=\frac{1}{\varepsilon}\Big(r_D\frac{\xi^4}{\Omega^2}-r_m\Big)\tau,\qquad
r_D=\frac{D_2}{D_1},\; r_m=\frac{m_2}{m_1}.$$

The two identical-plate relations become

$$\textbf{(sym)}\quad\Big(\frac{\xi^4}{\Omega^2}-1\Big)\,\tau\,\tan\frac{\tau}{2}+\varepsilon=0,
\qquad
\textbf{(antisym)}\quad\Big(\frac{\xi^4}{\Omega^2}-1\Big)\,\tau\,\cot\frac{\tau}{2}-\varepsilon=0.$$

Compare Sarkar–Sonti Eq. (11): $\big(\xi^4/\Omega^2-1\big)\,\tau\tan\tau+\varepsilon=0$. The
symmetric relation is obtained by $\tau\to\tau$, $\tan\tau\to\tan(\tau/2)$ — exactly the half-gap
substitution, as it must be.

### 4.1 Uncoupled ($\varepsilon=0$) skeleton

Setting $\varepsilon=0$ exposes the parents that the coupled branches perturb:

* **Structural** (both families): $\xi^4/\Omega^2-1=0\Rightarrow\xi=\Omega^{1/2}$ — the in-vacuo
  flexural wave. It appears in *both* families (in-phase and out-of-phase plate flexure).
* **Symmetric acoustic:** $\tau\tan(\tau/2)=0\Rightarrow$ either $\tau=0$ (the **plane wave**
  $\xi=\Omega$) or $\tau=2n\pi$ (the **even** rigid-duct cut-ons).
* **Antisymmetric acoustic:** $\tau\cot(\tau/2)=0\Rightarrow\tau=(2n{+}1)\pi$ (the **odd**, i.e.
  pressure-release-at-centre cut-ons). There is **no plane wave** — $\tau\cot(\tau/2)\to2\ne0$ as
  $\tau\to0$.

So the single cut-on ladder of the one-plate problem splits into two interleaved ladders, one per
symmetry family — the central qualitative new feature.

---

## 5. Asymptotic analysis (small $\varepsilon$)

The Sarkar–Sonti programme carries over branch-by-branch. As an illustration, perturb the structural
branch $\xi=\Omega^{1/2}+a_1\varepsilon+O(\varepsilon^2)$ and balance at $O(\varepsilon)$. Symbolic
expansion (`symbolic_verification.py`, check 4) gives the **first-order corrections**

$$\textbf{(sym)}\quad
\xi_s(\varepsilon)=\Omega^{1/2}-\frac{\varepsilon}{4\,l\sqrt{\Omega-1}}\,
\cot\!\Big(\tfrac{l}{2}\sqrt{\Omega^2-\Omega}\Big)+O(\varepsilon^2),$$

$$\textbf{(antisym)}\quad
\xi_a(\varepsilon)=\Omega^{1/2}+\frac{\varepsilon}{4\,l\sqrt{\Omega-1}}\,
\tan\!\Big(\tfrac{l}{2}\sqrt{\Omega^2-\Omega}\Big)+O(\varepsilon^2).$$

The symmetric expression is **identical to Sarkar–Sonti Eq. (14) under $l\to l/2$** — the half-gap
substitution again — which is a strong analytic consistency check. The antisymmetric expression is
new. As in the parent paper, these first-order forms break down near $\Omega\approx1$ (coincidence)
and near each cut-on $\Omega_c^n$, where the correction grows; there the appropriate local
expansions ($\Omega=1+\varepsilon C$, $\Omega=\Omega_c^n+\varepsilon C$, etc.) must be used, and at
each such point the coupling opens an $O(\sqrt\varepsilon)$ **gap** instead of a crossing — now
*two* interleaved sets of gaps, one per symmetry family. The large-$\varepsilon$ analysis proceeds
identically with $\varepsilon'=1/\varepsilon$ and yields branches that migrate from rigid-walled to
pressure-release wavenumbers.

> **Implementation note.** Higher-order corrections, the coincidence/cut-on local expansions, and the
> large-$\varepsilon$ expansions are all mechanical to generate with the supplied sympy script — just
> change the ansatz and the local variable. This is exactly the "simple to implement with a symbolic
> package" spirit the original authors emphasise.

---

## 6. The non-identical case — genuinely new behaviour

When $r_D\ne1$ or $r_m\ne1$ the factorisation fails and one must use the full ($\star\star$). The two
structural waves (in-phase and out-of-phase plate flexure) now have *different* in-vacuo
wavenumbers, and where they would cross they instead **veer** (avoided crossing), exchanging mode
shape — the classic signature of two coupled oscillators. Physically, asymmetry in stiffness/mass
breaks the mirror symmetry of the channel, so symmetric and antisymmetric motions are no longer
eigenstates and exchange energy. This regime has no counterpart in Fahy or Sarkar–Sonti and is the
most promising direction for a novel contribution; the asymptotic treatment proceeds by perturbing
($\star\star$) about the *coupled* identical-plate solutions, with $r_D-1$ and $r_m-1$ as additional
small parameters.

---

## 7. Numerical results and validation

All figures are produced by `dispersion_two_plates.py` (parameters $l=3$, $\varepsilon=0.25$, matching
Sarkar–Sonti's Fig. 5).

**`fig1_validation.png` — Validation.** The single-flexible-plate + rigid-wall solver reproduces the
Sarkar–Sonti coupled dispersion curves: coupled branches sit just off the in-vacuo flexural, plane-
wave and rigid cut-on parents, with gaps at coincidence and at $\Omega_c^n$.

**`fig2_two_identical.png` — Two identical flexible plates.** The symmetric (red) and antisymmetric
(orange) coupled families are overlaid on their uncoupled parents. The symmetric family hugs the
plane wave and even cut-ons; the antisymmetric family hugs the odd cut-ons and has no low-frequency
plane wave. Both perturb the common structural wave and open gaps at every intersection.

**`fig3_nonidentical.png` — Non-identical plates.** Comparison of the identical case ($r_D=1$) with a
stiffer second plate ($r_D=4$) via the full relation, showing the shift and the onset of veering
between the structural waves.

### Numerical cross-checks (all pass)

| Check | Result |
|---|---|
| Symmetric branch ($l=3,\varepsilon=0.25$) vs. Sarkar–Sonti single plate ($l=1.5,\varepsilon=0.125$) | identical to 4 d.p. at $\Omega=0.7,1.3,1.8$ |
| Full relation (identical plates) vs. (symmetric $\cup$ antisymmetric) | identical root sets (plane-wave $\xi=\Omega$ is a removable tangency) |
| Rigid-wall limit of ($\star\star$) | $\alpha_1=-\cot(k_y a)$ = Sarkar–Sonti Eq. (9) |

---

## 8. Suggested roadmap to a publishable / production-grade result

1. **Complete the asymptotic catalogue** for the antisymmetric family (coincidence, each cut-on,
   $O(\varepsilon^2)$, large-$\varepsilon$) — mechanical with the supplied sympy script.
2. **Non-identical plates:** derive the perturbation of ($\star\star$) in $(r_D-1, r_m-1)$ and obtain
   closed-form veering/gap widths between the structural waves.
3. **Damping & complex wavenumbers:** add structural loss factor $D_j\to D_j(1+\mathrm{i}\eta_j)$ and
   fluid absorption; track complex $\xi$ (the asymptotic series extends directly).
4. **Energy / group-velocity interpretation** of the gaps (Fuller & Fahy give the group-velocity
   no-crossing argument; reproduce it per symmetry family).
5. **Independent validation** against a numerical eigen-solver (FEM, or the FVM–FDM CFiA framework of
   Morab *et al.* 2024) to confirm the asymptotic curves beyond first order.

---

## 9. Reproducibility

| File | Purpose |
|---|---|
| `symbolic_verification.py` | sympy proof of ($\star$), the symmetric/antisymmetric factorisation, the rigid limit, and the first-order asymptotics |
| `dispersion_two_plates.py` | numerical root-tracker for all branches + figure generation + cross-checks |
| `fig1_validation.png` | reproduction of Sarkar–Sonti (validation) |
| `fig2_two_identical.png` | two identical flexible plates: symmetric + antisymmetric branches |
| `fig3_nonidentical.png` | non-identical plates via the full relation |

### Key symbols

| Symbol | Meaning |
|---|---|
| $D_j=E_jI_j$, $I_j=h_j^3/12(1-\nu_j^2)$ | flexural rigidity of plate $j$ |
| $m_j=\rho_{p,j}h_j$ | plate mass per unit area |
| $\rho,\,c$ | fluid density, sound speed |
| $k=\omega/c$, $k_x$, $k_y=\sqrt{k^2-k_x^2}$ | acoustic, axial, transverse wavenumbers |
| $G_j=D_jk_x^4-m_j\omega^2$ | in-vacuo plate operator |
| $\alpha_j=k_yG_j/\omega^2\rho$ | dimensionless plate fluid-loading admittance |
| $\Omega,\xi,l,\varepsilon,\tau$ | nondim. frequency, axial wavenumber, gap, fluid-loading, $k_ya$ |
| $r_D=D_2/D_1$, $r_m=m_2/m_1$ | plate stiffness / mass ratios |

*Prepared as a working extension document. The analytical results are exact reductions of, and
consistent with, Sarkar & Sonti (2007); the antisymmetric family and the non-identical-plate
coupling are the new contributions.*
