# Gate-2 ledger: widening to integral ⟨ξ,α^∨⟩

Goal: support covers with ⟨ξ,α^∨⟩ ∈ ℤ (not just 0), keeping (1−δ)ξ ∈ X*.
Motivation: δ-fixed representatives (e.g. GL(2,R) det/2 with ξ'=(−1/2,1/2))
live here; with them twist is honest and hermitian forms computable.
Invariant tests after EVERY batch: xi=0 baseline identical; cover-smoke,
tensor-bijection, hermitian-cover outputs identical ((1D) covers unaffected).

Derivation rule: storage = λ−ρ−ξ. Any site reconstructing a ⟨λ,α^∨⟩-class
from storage + ⟨ρ,α^∨⟩-constants needs +⟨ξ,α^∨⟩ (helper Rep_context::xi_level,
integral under gate). Affine reflection pivots anchored at basepoint get
offset 1+⟨ξ,α_s^∨⟩; pivots anchored at ρ_r/γ (basepoint-free) unchanged.
(1−θ_x)ξ ∈ X* for all x: since ⟨ξ,coroots⟩∈ℤ, (w−1)ξ ∈ root lattice, and
(1−δ)ξ∈X* by gate. (1+θ)ξ = 2ξ−(1−θ)ξ ∈ X* ✓ theta_plus_1_rho_xi stays Weight.

## Batches
- [ ] A: constructor widening + xi_level helper + colevel sites
      (repr.cpp reducibility_points; K_repr theta_plus_1_eval, is_standard,
      is_nonzero)
- [ ] B: parity mod-2/4 sites (repr: is_parity_at_0, is_semifinal, is_final,
      finals_for eval_lr, deformation_unit lambda_rho_real2, any_Cayley parity
      vector, Cayley rho_r_corr parity tests, common_context::is_parity/
      down_Cayley/up_Cayley)
- [ ] C: affine offsets (repr: make_dominant, complex_crosses, deform_readjust,
      finals_for; K_repr: make_dominant, make_theta_stable,
      to_canonical_involution, normalise, finals_for) — derive real-case pivots
      per site
- [ ] D: height parity asserts (repr height, K_repr heights) — audit
      ⟨(1+θ)ξ,2ρ^∨⟩ parity, floor or relax
- [ ] E: tests — GL(2,R) presentation-independence (ξ=(1/2,1/2) vs δ-fixed
      ξ'=(−1/2,1/2): same λ ⟹ storage shift by (1,0); all invariants equal);
      hermitian forms via κ=0 presentation; R^×-cover unchanged

## Site log
(appended as done)

- [x] A (done): constructor now requires integral ⟨ξ,simple coroots⟩ (was: zero);
  Rep_context::xi_level helper; corrections at reducibility_points lam_alpha,
  K_repr theta_plus_1_eval (alpha and beta terms), is_standard, is_nonzero.
  All batteries pass (only intended error-message diffs; SL(2,R) xi=1/2 still
  rejected: pairing 1/2 non-integral).

- [x] B: parity sites (commit 754b5f2).
- [x] C: affine pivots: make_dominant offsets complex 1+xl / real xl (class
  mod 2 suffices, y_pack renormalises); all 12 simple_reflect(s,lr,1) sites
  -> 1+xl; K_repr KGP real cross offset xl; integral make_dominant: xi term
  of V=2(lambda-xi)-2rho_r stays unreflected (correction after rd.reflect).
- [x] D: height evenness assert weakened on covers (floor division).
- [x] E: gate2.in. RESULTS: presentation-independence of blocks/KL/deform
  EXACT (xi=(1/2,1/2) vs delta-fixed (-1/2,1/2)); singular-gamma finality/
  block-membership consistency restored in both presentations.
  FINDING: h = normal*twist fails on ANY shifted coset (even kappa=0):
  the absorption move's parity is controlled by half-integrality of lambda
  (2(b,a)+(-1,1) parity computation), i.e. the identity is calibrated to
  rho+X^* specifically. Therefore: same-x h for ALL covers; forms allowed
  on equal-rank covers only (twist=id there), else gated pending
  extended-parameter theory on covers.
