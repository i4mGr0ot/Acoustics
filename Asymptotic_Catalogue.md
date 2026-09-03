# Full Asymptotic Catalogue — Two-Flexible-Plate Structural–Acoustic Waveguide

### Complete small-ε and large-ε expansions for the symmetric and antisymmetric coupled wavenumbers

*Companion to `Two_Flexible_Plates_Extension.md`. Every formula below is derived in the
Sarkar–Sonti (2007) style and has been checked symbolically and validated against exact numerical
roots (errors quoted). Conventions: $\Omega=\omega/\omega_c$, $\xi=k_x/k_c$, $l=k_c a$,
$\varepsilon=\rho a/m$, $\tau=k_y a=l\sqrt{\Omega^2-\xi^2}$,
$S(\xi)=\xi^4/\Omega^2-1$, $\tau_s=l\sqrt{\Omega^2-\Omega}=l\sqrt{\Omega}\sqrt{\Omega-1}$.*

The two identical-plate dispersion relations are

$$\textbf{(SYM, }v=0\text{ at centre)}\quad S(\xi)\,\tau\,\tan\tfrac{\tau}{2}+\varepsilon=0,
\qquad
\textbf{(ANTI, }p=0\text{ at centre)}\quad S(\xi)\,\tau\,\cot\tfrac{\tau}{2}-\varepsilon=0.$$

---

## A. Small fluid-loading parameter ($0<\varepsilon\ll1$)

Each coupled wavenumber is its uncoupled parent plus an $O(\varepsilon)$ correction, **except** in
narrow $O(\varepsilon)$-wide windows around coincidence ($\Omega=1$) and each cut-on, where a local
$\sqrt\varepsilon$ expansion is needed (§C).

### A1. Coupled structural wave $\xi_s(\varepsilon)$ — both families

$$\boxed{\;\xi_s^{\text{sym}}(\varepsilon)=\Omega^{1/2}
-\frac{\varepsilon}{\,4\,l\sqrt{\Omega-1}\,\tan(\tau_s/2)}\;},\qquad
\boxed{\;\xi_s^{\text{anti}}(\varepsilon)=\Omega^{1/2}
+\frac{\varepsilon\,\tan(\tau_s/2)}{\,4\,l\sqrt{\Omega-1}\,}\;}.$$

For $\Omega<1$ (subsonic structure) write $\sqrt{\Omega-1}\to i\sqrt{1-\Omega}$ and
$\tan\to i\tanh$, giving real, positive corrections (mass loading: the wavenumber rises above its
in-vacuo value, as in Sarkar–Sonti §4.1).

* **Validity:** away from $\Omega=1$ and from any cut-on crossing.
* **Consistency:** the symmetric form is *exactly* Sarkar–Sonti Eq. (14) under the half-gap
  substitution $l\to l/2$ (their single-plate problem is our symmetric mode).
* **Validation (l=3, ε=0.05):** sym $\le0.08\%$ ($\Omega\in[1.2,1.55]$), $\le0.49\%$
  ($\Omega\in[0.4,0.85]$); anti $\le0.22\%$ ($\Omega\in[1.15,1.5]$), $\le0.009\%$
  ($\Omega\in[0.4,0.9]$).

> **Key physical contrast.** As $\Omega\to1$ the symmetric correction diverges
> ($\tan(\tau_s/2)\to0$) — it needs the local expansion C1. The antisymmetric correction is
> **finite**, $a_1\to\tfrac18$, so the **antisymmetric structural wave passes through coincidence
> smoothly** (no gap there). This is because the symmetric family contains the plane wave (which the
> structure collides with at $\Omega=1$) while the antisymmetric family does not.

### A2. Coupled acoustic plane wave $\xi_a(\varepsilon)$ — symmetric family only

$$\boxed{\;\xi_a^{\text{sym}}(\varepsilon)=\Omega+\frac{\varepsilon}{\,\Omega\,l^2(\Omega^2-1)\,}\;}.$$

The antisymmetric family has **no plane wave**. For $\Omega>1$ the correction places the coupled
plane wave slightly *above* $\Omega$, i.e. in the evanescent region $\xi>\Omega$.

* **Validity:** away from $\Omega=1$.
* **Validation (l=3, ε=0.05):** $\le0.016\%$ ($\Omega\in[1.3,2.3]$), $\le0.33\%$
  ($\Omega\in[0.3,0.8]$). Coefficient confirmed to ratio $1.000$ as $\varepsilon\to0$.

### A3. Coupled acoustic cut-on $\xi_a^{(m)}(\varepsilon)$ — both families (unified)

$$\boxed{\;\xi_a^{(m)}(\varepsilon)=\sqrt{\Omega^2-\Big(\tfrac{m\pi}{l}\Big)^2}
\;+\;\frac{2\,\varepsilon}{S_0\,l^2\,\xi_0}\;},\qquad
\xi_0=\sqrt{\Omega^2-\Big(\tfrac{m\pi}{l}\Big)^2},\;\;S_0=\frac{\xi_0^4}{\Omega^2}-1,$$

with $m$ **even** for the symmetric family ($k_y a=2n\pi$) and $m$ **odd** for the antisymmetric
family ($k_y a=(2n{+}1)\pi$). (Sarkar–Sonti's single-plate version, Eq. 21, has all integer $m$ and
the prefactor $1/S_0l^2\xi_0$; the factor $2$ here is the half-angle of the two-plate problem.)

* **Validity:** above the cut-on frequency, away from where the structural wave crosses it.
* **Validation (l=3, ε=0.05):** sym $m=2$: $\le0.075\%$; anti $m=1$ interior $\sim1\%$ (rising at the
  structural crossing, where C2 takes over).

---

## B. Large fluid-loading parameter ($\varepsilon\gg1$): pressure-release modes

Set $\varepsilon'=1/\varepsilon\ll1$. The parents become **pressure-release** cut-ons. The leading
coupled correction is, for both families (unified),

$$\boxed{\;\xi_{a}'^{(M)}(\varepsilon)=\sqrt{\Omega^2-\Big(\tfrac{M\pi}{l}\Big)^2}
\;-\;\frac{2\,(M\pi)^2\,S_0}{l^2\,\xi_0}\,\frac{1}{\varepsilon}\;},\qquad
\xi_0=\sqrt{\Omega^2-\Big(\tfrac{M\pi}{l}\Big)^2},\;\;S_0=\frac{\xi_0^4}{\Omega^2}-1,$$

with $M$ **odd** for the symmetric family and $M$ **even** for the antisymmetric family (including
$M=0$, a pressure-release "plane-wave-like" antisymmetric mode).

* **Validation (l=3, ε=20):** sym $M=1$: $0.14$–$1.8\%$; anti $M=2$: $2.8$–$3.9\%$ (interior),
  rising near its own cut-on — all consistent with a first-order $\varepsilon'=0.05$ expansion.

### The rigid → pressure-release migration, resolved per family

| | small $\varepsilon$ (rigid walls) | large $\varepsilon$ (pressure-release walls) |
|---|---|---|
| **Symmetric** | $k_y a = 2n\pi$ (even) + plane wave | $k_y a=(2m{+}1)\pi$ (odd) |
| **Antisymmetric** | $k_y a=(2n{+}1)\pi$ (odd) | $k_y a=2m\pi$ (even) + $M{=}0$ mode |

Each family migrates from the rigid-duct cut-on ladder to the pressure-release ladder as
$\varepsilon$ grows — a per-symmetry refinement of the single migration Sarkar–Sonti found. The two
ladders **interleave** in frequency, so the full two-plate spectrum is twice as dense as the
one-plate spectrum.

---

## C. Local $\sqrt\varepsilon$ expansions (transition windows)

Where two uncoupled curves intersect, coupling forbids the crossing and opens a gap of width
$O(\sqrt\varepsilon)$; the branches **swap identity** across it.

### C1. Symmetric coincidence ($\Omega=1$: structure meets plane wave)

$$\boxed{\;\xi(\varepsilon)=\Omega\pm\frac{\sqrt\varepsilon}{2l}\;}\qquad(\text{gap }=\sqrt\varepsilon/l).$$

The "+" branch above $\Omega=1$ continues the structural branch from below, and vice-versa (sign
fixed by continuity, exactly as Sarkar–Sonti §4.3).
*(Validation l=3, ε=0.03: predicted $\{0.9711,1.0289\}$ vs exact $\{0.9709,1.0287\}$.)*

### C2. Cut-on crossings ($\Omega=\Omega_c$: structure meets a cut-on) — both families

At the crossing $\Omega_c$ (where $l^2(\Omega_c^2-\Omega_c)=(m\pi)^2$, so $\xi_c=\sqrt{\Omega_c}$),

$$\boxed{\;\xi(\varepsilon)=\sqrt{\Omega_c}\pm\frac{\sqrt{2\varepsilon}}{2l}\;}\qquad
(\text{gap }=\sqrt{2\varepsilon}/l).$$

*(Validation l=3, ε=0.03: SYM $m{=}2$ at $\Omega_c{=}2.6533$, predicted $\{1.5881,1.6697\}$ vs exact
$\{1.5872,1.6689\}$; ANTI $m{=}1$ at $\Omega_c{=}3.6811$, predicted $\{1.8778,1.9595\}$ vs exact
$\{1.877,1.9587\}$.)*

### C3. Antisymmetric coincidence

No gap — see the boxed note after A1: the antisymmetric structural wave is regular at $\Omega=1$
with $\xi\approx\Omega^{1/2}(1+\varepsilon/8)$.

---

## D. Composite picture (Fig. 4, `fig4_catalogue.png`)

The figure overlays every closed-form expression above (coloured lines / green gap-markers) on the
exact numerical roots (grey). Left = symmetric family, right = antisymmetric family, $l=3$,
$\varepsilon=0.15$. The agreement is visually exact everywhere the relevant expansion is valid,
including the swap of branches across each $\sqrt\varepsilon$ gap. As in Sarkar–Sonti, a single
*physical* dispersion curve is a patchwork: e.g. on the symmetric side, plane-wave → (coincidence
gap) → structural → (cut-on gap) → cut-on.

---

## E. Summary table of the catalogue

| Region | Family | Formula | Order | Max error (validated) |
|---|---|---|---|---|
| Structural, away | sym | $\Omega^{1/2}-\varepsilon/[4l\sqrt{\Omega-1}\tan(\tau_s/2)]$ | $\varepsilon$ | 0.08–0.49% |
| Structural, away | anti | $\Omega^{1/2}+\varepsilon\tan(\tau_s/2)/[4l\sqrt{\Omega-1}]$ | $\varepsilon$ | 0.009–0.22% |
| Plane wave, away | sym | $\Omega+\varepsilon/[\Omega l^2(\Omega^2-1)]$ | $\varepsilon$ | 0.016–0.33% |
| Cut-on, away | both | $\xi_0+2\varepsilon/(S_0 l^2\xi_0)$, $m$ even/odd | $\varepsilon$ | 0.08–1% |
| Coincidence | sym | $\Omega\pm\sqrt\varepsilon/(2l)$ | $\sqrt\varepsilon$ | <0.03% |
| Coincidence | anti | regular, $\Omega^{1/2}(1+\varepsilon/8)$ | $\varepsilon$ | — |
| Cut-on crossing | both | $\sqrt{\Omega_c}\pm\sqrt{2\varepsilon}/(2l)$ | $\sqrt\varepsilon$ | <0.06% |
| Pressure-release (large ε) | both | $\xi_0-2(M\pi)^2 S_0/(l^2\xi_0)\cdot\varepsilon^{-1}$, $M$ odd/even | $\varepsilon^{-1}$ | 0.1–3.9% |

All derivations: `asymptotic_catalogue.py` (symbolic) and `validate_catalogue.py` (numerical).
Figure: `fig4_catalogue.png`.
