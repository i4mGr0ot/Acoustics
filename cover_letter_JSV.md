# Cover letter — Journal of Sound and Vibration

**To:** The Editors, *Journal of Sound and Vibration*

**Manuscript:** *Coupled dispersion of a two-flexible-plate structural–acoustic waveguide carrying a
Poiseuille mean flow: closed-form low-Mach corrections, critical-layer structure and wall-mode
survival*

**Author:** Ashray Saxena, Department of Mechanical Engineering, BITS Pilani (Pilani Campus)

---

Dear Editors,

I am pleased to submit the manuscript above for consideration as a research article in the *Journal of
Sound and Vibration*.

**What the paper does.** It asks what happens to the coupled structural–acoustic wavenumbers of a
fluid column bounded by two thin elastic plates when the column is set in motion by a fully developed
Poiseuille mean flow. With a sheared mean profile the pressure perturbation no longer satisfies the
Helmholtz equation but the Pridmore–Brown equation, and the paper derives that equation for this
geometry, solves it asymptotically in the centreline Mach number, and then examines what happens once
the flow becomes fast enough for a critical layer to appear inside the channel.

**The principal contributions are the following.**

1. **A closed-form first-order Mach correction to every coupled wavenumber.** The leading operator of
   the Mach expansion is shown to be self-adjoint under the plate-admittance boundary conditions, so a
   Fredholm solvability argument applies with the quiescent mode shape as its own adjoint null
   function. This yields a fully explicit expression for the correction in terms of the quiescent mode
   shape, its two endpoint values, and two elementary quadratures — no differential equation need be
   solved. The correction is proved negative on every branch, and it separates cleanly into
   boundary-shear, bulk-shear and Doppler-convection contributions. Agreement with numerically
   differentiated d*ξ*/d*M* is 0.1–0.5%.

2. **An exact structural simplification specific to this configuration.** Because a Poiseuille profile
   satisfies no-slip, the mean velocity vanishes precisely where the compliant-wall boundary condition
   is applied. The convective (Ingard–Myers) term is therefore identically zero rather than merely
   small, the plate-admittance conditions retain their quiescent Robin form at *every* Mach number, and
   the well-posedness difficulties that attend the Ingard–Myers condition in lined-duct acoustics do
   not arise. This is what makes the closed-form analysis possible.

3. **A critical-layer threshold controlled by a property of the structure.** The layer appears once
   *M* ≥ *M*_crit = Ω/ξ. Because the structural wave has ξ = Ω^(1/2), its threshold is √Ω, which is
   subsonic precisely below coincidence; every other branch requires sonic or supersonic flow.
   Coincidence is a plate property, so the accessibility of a critical layer in this waveguide is set
   by the structure and not by the fluid — a feature with no counterpart in a rigid or
   constant-impedance duct.

4. **The main physical result: the discrete wall modes are not absorbed by the critical layer.**
   Solving the viscous linearised Navier–Stokes problem shows the eigenvalue is regular through and
   above *M*_crit; the attenuation follows a clean *Re*^(−1/2) law, identifying it as wall Stokes-layer
   damping that vanishes inviscidly; and the *Re*-independent critical-layer absorption is zero to
   within 10⁻⁵. The reason is geometric: these modes are plate-driven and have small amplitude at the
   centreline, where the layer sits. This contrasts with the impedance-wall surface modes of Brambley,
   Darau & Rienstra (*JFM* 710, 2012), which the critical layer can absorb substantially, and suggests
   the general principle that what matters is the geometric overlap between layer and mode rather than
   the strength of the singularity.

**Relation to existing work, stated plainly.** The quiescent two-plate material in Sections 3 and 4 —
the coupled dispersion relation, its exact symmetric/antisymmetric factorisation, and the small- and
large-ε asymptotic catalogue — recovers, for this geometry, results consistent with the
asymptotic-tracking framework of Sarkar & Sonti (*JSV* 306, 2007), Sarkar, Kunte & Sonti (*CMES* 81,
2011) and Prakash & Sonti (*JSV* 373, 2016). I do not claim it as novel, and the manuscript says so
explicitly in Section 1.6 and again in the conclusions. It is included because it fixes the notation,
supplies the *M* = 0 eigenpairs that appear inside the closed-form mean-flow corrections, and makes the
paper self-contained. Its only incremental content is the antisymmetric family, which the
one-flexible-plate geometry cannot see, and the non-identical-plate veering. The new contribution of
the paper is the mean flow and everything that follows from it.

Equally, I want to be clear that sheared-flow duct acoustics and its critical layer are mature
subjects in their own right (Pridmore-Brown 1958; Maslowe 1986; Brambley, Darau & Rienstra 2012; King
*et al.* 2022). What is new here is the combination of that machinery with a plate-bounded waveguide,
the closed-form low-Mach corrections that the no-slip profile makes possible, and the wall-mode
conclusion.

**Fit with JSV.** The paper sits directly in the line of work JSV has published on asymptotic tracking
of coupled wavenumbers in flexible structural–acoustic waveguides, and extends it in the one direction
that line has not taken — a mean flow. The results are of the closed-form, directly usable kind that a
designer of a compliant duct or a flow-bearing elastic waveguide can apply, and the paper reports the
branch-dependent range of validity of every expression so that a reader knows when to stop trusting
them.

**Verification.** Three mutually independent numerical schemes — Runge–Kutta shooting, Chebyshev
spectral collocation, and a primitive-variable formulation reduced symbolically to Pridmore–Brown form
— agree on every branch to approximately 6 × 10⁻⁹. The algebra of the quiescent relation and the
equivalence of the inviscid and viscous wall conditions were verified with a symbolic-algebra package.

The manuscript is original, has not been published elsewhere, and is not under consideration by any
other journal. I have no competing interests to declare. I would be glad to suggest referees if that
is helpful.

Thank you for your consideration.

Yours sincerely,

**Ashray Saxena**
Department of Mechanical Engineering
Birla Institute of Technology & Science Pilani, Pilani Campus
Rajasthan 333031, India
f20231043@pilani.bits-pilani.ac.in

---

## Optional paragraph — prior submission history

*You are under no obligation to disclose a desk rejection, and many authors do not. Including a
sentence like the one below costs nothing and reads as good faith; delete it if you prefer. If you
instead decide to submit to JFS as a clean restart, this paragraph is **not** optional — the Editor
asked explicitly to be reminded of the correspondence — and the alternative wording is given after it.*

**For JSV (optional):**

> For completeness: an earlier and considerably shorter version of this work was declined without
> review by another journal, on the grounds that it presented its results without exposing the
> underlying technical detail. That assessment was fair. The present manuscript is a substantially
> different document: the governing equations are now stated and derived in full, the solvability
> argument and the critical-layer analysis are worked through step by step rather than summarised, the
> literature review has been rebuilt, and the contribution has been reframed so that the mean-flow
> results lead and the quiescent material is explicitly identified as recovered baseline.

**For JFS, if submitted there instead (required by the Editor):**

> This manuscript is submitted as a clean restart, not a revision. An earlier and much shorter version
> was declined without review by JFS, and the Editor-in-Chief kindly provided informal notes: that the
> paper hid too much technical detail, that the governing equations were not even stated, that critical
> layers cannot be assumed familiar to the general JFS reader, that the novelty was therefore hard to
> assess, and that the literature review was limited. He invited me to mention this correspondence if I
> chose to submit a new version. The present manuscript addresses each point directly — the linearised
> Euler system, the Pridmore–Brown derivation, the plate equations and the boundary conditions are now
> given in full; the critical-layer section is written to be self-contained; the literature review has
> been rebuilt; and the contribution is stated explicitly, with the quiescent material identified as
> recovered baseline rather than claimed as new. I understand that fit with JFS remains an open
> question and would welcome the Editor's view.

---

## Editor's criticisms → where the manuscript now answers them

| Criticism (JFS Editor) | Response in the rewritten manuscript |
|---|---|
| "You're hiding so much of the technical details" | Sections 2, 5 and 6 derive everything from first principles; five appendices carry the longer algebra. Length grew from ~12 to ~37 pages. |
| "We don't even know what the governing equations are" | §2.3 states the linearised Euler system (eqs. 12–14) explicitly; §5.2 derives the Pridmore–Brown equation in four labelled steps; §2.4 gives the plate equation and both interface conditions; §7.2 writes out the viscous LNS system. |
| "Critical layers aren't for the faint-hearted; you cannot assume the reader is familiar" | §8 is rewritten as a self-contained tutorial: what the singularity is physically, where it appears, the Frobenius structure, why causality selects the branch, and how viscosity regularises it — with the indicial equation, the Landau prescription and the Airy balance all shown. Appendix C carries the recurrence. |
| "The novelty of the contribution is difficult to assess" | §1.6 lists seven contributions with the status of each stated explicitly as *new*, *classical*, or *recovered baseline*. The abstract, introduction and conclusions all say which parts are not novel. |
| "The literature review is quite limited" | §1 is now six subsections covering fluid loading, asymptotic tracking, sheared-flow duct acoustics, the impedance boundary-condition controversy, critical layers, and flow over compliant walls. The bibliography went from 11 to 50 entries. |
| "I'm not sure JFS would be the right journal" | Retargeted to JSV, which published the asymptotic-tracking line this work extends. |

---

## Before you submit — please check

- [ ] **Verify every bibliography entry** against the actual paper. The list was assembled from
      literature searches and author recollection; volume, page and year details should be confirmed
      individually, particularly refs. 10 (Prakash & Sonti 2019, page numbers omitted), 17 (Raviprolu
      — confirm co-author), 24–28 (the older sheared-flow duct papers) and 14 (Choy & Huang).
- [ ] Compile `manuscript_JSV.tex` in the `Acoustics` folder so the ten figures resolve
      (`pdflatex` × 3 for the table of contents and cross-references).
- [ ] Decide whether to keep the table of contents and line numbers — both are useful for review but
      JSV may want them removed for the final version.
- [ ] Confirm the Table 4 shear-contribution figures against a fresh run of `om2_crit.py`.
- [ ] Consider adding a schematic figure of the geometry (two plates, fluid column, Poiseuille
      profile) as Fig. 1 — reviewers of a geometry-specific paper usually expect one, and none of the
      ten existing figures shows the configuration itself.
