# Mean Flow, Veering, Damping, O(M²) & Critical Layer — Toward a Publishable Result

### Poiseuille (sheared) mean flow via the Pridmore–Brown equation, with second-order Mach analysis, a critical-layer study, and independent validation

*Companion to `Two_Flexible_Plates_Extension.md` and `Asymptotic_Catalogue.md`. All results validated
against a numerical Pridmore–Brown BVP solver (`pb_solver.py`), an independent Chebyshev spectral
eigensolver (`pb_spectral.py`), and a symbolic primitive-variable reduction.*

> **Correction notice.** An earlier draft used the Pridmore–Brown shear term with the wrong sign. An
> independent primitive-variable (linearized-Euler) derivation showed the correct middle term is
> **`+2Mξf′/σ̂`** in the `e^{i(kx x − ωt)}`, `σ̂ = Ω − Mξ f` convention. All solvers, the O(M) formula,
> and the numbers below are now sign-corrected. One earlier claim — that the antisymmetric structural
> wave is "flow-insensitive" — was an artefact of the wrong sign and is **withdrawn** (its true
> ξ₁ ≈ −0.11).

---

## 1. Pridmore–Brown problem (sign-corrected)

Poiseuille mean flow `U(y) = 4U₀ (y/a)(1−y/a) = U₀ f(η)`, `f = 4η(1−η)`, `η = y/a`, no-slip
(`U=0` at both walls), Mach `M = U₀/c`. Pressure perturbations obey

```
p̂_ηη + (2Mξ f′/σ̂) p̂_η + l²(σ̂² − ξ²) p̂ = 0,   σ̂ = Ω − Mξ f,   f′ = 4(1−2η)
```

with plate-admittance (Robin) boundary conditions `p̂_η(0) + (ε/S₁)p̂(0) = 0`,
`p̂_η(1) − (ε/S₂)p̂(1) = 0`. No-slip kills the Ingard–Myers wall correction, so the flow enters only
through the bulk operator. At M=0 this reduces exactly to the quiescent relation
`tan(k_y a) = (α₁+α₂)/(1−α₁α₂)`.

## 2. Closed-form O(M) correction (sign-corrected)

```
        −(2/Ω)(p₀(1)²+p₀(0)²) + (4/Ω)I₀ − l²Ω I_f
ξ₁ = ──────────────────────────────────────────────────── ,   I₀=∫p₀², I_f=∫f p₀²
        l²I₀ + (2εξ₀²/Ω²)(r_D p₀(1)²/S₂² + p₀(0)²/S₁²)
```

Validated to 0.1–0.5% vs numerical dξ/dM. Corrected table (l=3, ε=0.25):

| Ω | branch (ξ₀) | type | ξ₁ (closed) | ξ₁ (numerical) |
|---|---|---|---|---|
| 1.3 | 1.1165 | sym structural | −0.1283 | −0.1290 |
| 1.3 | 1.1837 | anti structural | −0.1065 | −0.1068 |
| 1.3 | 1.3267 | plane wave | −0.7471 | −0.7468 |
| 1.8 | 1.8069 | plane wave | −1.1743 | −1.1750 |

**Physical readings:** ξ₁ < 0 on every branch (downstream flow lowers the axial wavenumber, upstream
raises it — the dispersion relation loses kx → −kx symmetry); the plane wave is the most
flow-sensitive (it samples the full sheared profile); the two structural branches are comparable and
more modest.

## 3. Second-order Mach correction + convection/shear split

ξ(M) = ξ₀ + ξ₁M + ξ₂M² + …  ; ξ₂ from Richardson-extrapolated central differences, decomposed into a
convection part (only σ̂ in the σ̂²−ξ² term) and a shear part (the U′ term). At Ω=1.3:

| branch (ξ₀) | ξ₁ | ξ₂ (full) | ξ₂ (convection) | shear share |
|---|---|---|---|---|
| 0.685 (low) | −0.851 | +0.562 | +0.597 | −0.034 |
| 1.1165 (sym struct) | −0.130 | −0.652 | −0.811 | +0.160 |
| 1.1837 (anti struct) | −0.107 | +0.416 | +0.178 | +0.239 |
| 1.3267 (plane wave) | −0.747 | +1.647 | +2.264 | −0.617 |

The shear (U′) contribution to ξ₂ reaches ~37% (plane-wave branch), so a uniform-flow (plug) model
would mispredict the second-order term — the full Poiseuille profile matters. ξ₂ is symmetric in M (a
net, direction-independent effect), unlike the antisymmetric O(M) Doppler shift, and including it
extends the valid M-range (Fig. 6A).

## 4. Critical-layer analysis

The PB operator is singular where σ̂(η) = Ω − Mξf = 0 (wave phase speed = local flow speed). Since
f peaks at the centreline, a critical layer first appears there when

```
M ≥ M_crit = Ω/ξ = c_ph     (structural wave: M_crit ≈ √Ω)
```

So the regime splits at coincidence: **below** coincidence (Ω<1) the structural wave is subsonic
(ξ>Ω), M_crit<1, and a critical layer is reachable at **subsonic** flow; **above** coincidence it is
supersonic and no critical layer forms below sonic flow. For M ≪ M_crit the low-Mach expansion is
regular and accurate; as M → M_crit the inviscid eigenfunction develops a critical-layer singularity
at the centre and a viscous/contour treatment is needed (next analytical step). The low-Mach, low-ε
wall modes here sit safely in the regular regime. (Fig. 6B–C.)

## 5. Veering, damping, group velocity (quiescent — unaffected by the flow-sign fix)

- **Veering:** the two structural waves avoid-cross as plates differ; minimum gap at r_D=r_m=1,
  `Δξ_veer = ε/(2l√(Ω−1) sin τ_s)`, validated ∝ ε.
- **Damping:** structural loss D→D(1+iη) gives complex ξ with spatial decay; structural branches
  attenuate most, plane wave least.
- **Group velocity:** dΩ/dξ changes sharply across each √ε gap (energy switches carrier between
  structure and fluid); becomes direction-dependent with flow.

## 6. Independent validation (three routes)

1. **Shooting** (RK) and **Chebyshev spectral collocation** eigensolvers agree to ~6×10⁻⁹ on every
   branch up to M=0.2.
2. Both reduce at M=0 to the symbolically-verified closed-form quiescent relation.
3. A **primitive-variable** (3-field u′,v′,p′ linearized-Euler) formulation, eliminated symbolically,
   reproduces the Pridmore–Brown operator exactly — and is what exposed the shear-term sign error.

## 7. Roadmap status

| Item | Status |
|---|---|
| O(M) convective correction (sign-corrected) | done, validated 0.1–0.5% |
| O(M²) + convection/shear split | done |
| Critical-layer threshold M_crit = Ω/ξ | done |
| Non-identical-plate veering; damping; group velocity | done |
| Independent spectral + primitive-variable validation | done (~10⁻⁹) |
| Viscous/nonlinear critical layer (M ≳ M_crit) | open (future) |

### Reproducibility
`pb_solver.py` (shooting BVP, shear toggle), `pb_spectral.py` (independent Chebyshev eigensolver),
`pb_asym.py` (closed-form O(M)), `om2_crit.py` (O(M²) + shear split + M_crit), `fig6_data.py` (Fig. 6
data), `veer_damp_vg.py` (veering/damping/group velocity).
