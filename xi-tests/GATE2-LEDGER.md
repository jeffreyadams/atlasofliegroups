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

## induction.at (2026-07-11, commit dc2f514)
Forward direction cover-aware: theta_induce shift = rho_xi(G)-rho_xi(L);
rho_u_cover(L,G); Aq_genuine relaxing the integrality REQUIRE. Verified:
classic battery identical to stock binary; zero-shift consistency; integral
collapse; half-integral GL(2,R) Aq with exact tensoring dictionary.
DEFERRED: inverse maps (theta_stable_data, theta_stable_quasi_data, line-892
inverse) for cover-of-G inputs — need covered-Levi construction inside;
real_induce on matching-xi covers (expected to work verbatim — untested);
goodness classification on covers (intrinsic rho_u, likely fine — untested);
Aq_genuine cases where lambda-rho(u) has half-integral G-pairings are
correctly rejected (general-xi regime, v2).

## Endoscopy lifting port (E0 cont.) — 2026-07-11, scope map
Structure file DONE+verified (476b2a2). Lifting (endoscopic_lifting.at, the
L-homomorphism ε* transfer) sits on a substrate master & branch evolved
DIFFERENTLY. Reconciliation status:
- structure_constants.at: DONE. endoscopy_sc_support.at adds branch-only
  pieces (orientation type, m_pm int-overloads, orientation_function_vector);
  master's orientation_function/m_pm(vec,vec) VERIFIED to agree (same Geck /
  2-colouring), so signs are safe.
- tits.at: 161 diff-lines vs branch (tits_simple_reflection + more branch-only
  Tits fns; master's tits.at diverged, used by other scripts -> must
  supplement not clobber). THE BIG remaining piece; Tits sign conventions
  delicate.
- character_table_reps.at: 16 diff-lines. weak_packets.at: 32.
- copied branch files (homomorphism/L_hom/embed/param_pair_hash) carry
  branch-vs-master SYNTAX drift to fix (param_pair_hash.at:96 type error).
Copied-but-unwired; all.at unaffected.

## Endoscopy lifting port — tits.at RE-ASSESSED + cascade status (2026-07-11)
KEY CORRECTION: the "161 delicate diff-lines" in tits.at were a FALSE ALARM.
Reality: master defines NOTHING the branch lacks; branch tits.at = master's +
13 additive Tits-element fns; shared fns (multiply, inverse, left/right,
conjugate) FUNCTIONALLY IDENTICAL (only whitespace/comments/assert-syntax
differ; master is actually newer on @: asserts). No sign-convention risk.
=> clean additive supplement, like structure_constants.
Layers DONE (all clean additive supplements, verified no math divergence):
  structure_constants (endoscopy_sc_support.at), tits (endoscopy_tits_support.at
  - 13 fns + left/right ratvec overloads), synthetic2.at, to_dominant.
Cascade remaining (converging; each = a small branch-only helper/file):
  in_coweight_lattice (maybe_KGB.at dep) -> strong_real_forms -> ... then
  param_pair_hash.at:96 type-drift (embed.at), weak_packets (32-line diff).
RECOMMENDATION to finish efficiently: stop error-at-a-time; compute the FULL
transitive closure of branch-new/modified files the epstar lifting needs and
batch-port, rather than reactive one-by-one. Mechanical, low math-risk.
