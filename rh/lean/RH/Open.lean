/-
RH/Open.lean -- the ladder bookkeeping of the two former open edges
(r273 RETYPE: the open statements themselves moved to RH/Window.lean).

Provenance: the post-r259 endgame state, ledger rows
  * PRIME.PORT.COUPLEDTAU.TERMINAL_CROSSRATIO.01  [O]  (edge A, fiber)
  * PRIME.PORT.FULLSOURCE.PREFIX_RESUMMATION.01   [O]  (edge B, base)
registered by verification/v959_coupledtau_terminal_dictionary.py; campaign
state r260-r267 (probes terminal_crossratio_probe.py, terminal_triangle_
probe.py, cancellation_adjudication_probe.py, quenched_opening_probe.py,
s_monotonicity_probe.py, border_resolvent_identity_probe.py,
ranktrace_adjudication_probe.py; charter SHA 7b9e751d, r264).

r273 RETYPE (reviewer repair).  This file formerly stated the two edges
as universally quantified `sorry` theorems over the bare `WindowLadder`
(`terminal_crossratio_cofinal`) and over an arbitrary source predicate
(`prefix_resummation_positivity`).  Both statements are REFUTABLE in
that form -- the machine-checked counterexamples are permanent guards in
RH/Counterexamples.lean (`terminal_crossratio_not_universal`,
`prefix_resummation_not_universal`): the ladder fields carry no coupling
between `F_k` and `h_k`, and an arbitrary predicate carries no
positivity.  The truth-capable forms now live in RH/Window.lean and
(since the r305 reconstruction) RH/Closure.lean:
  * master theorem `augmented_prefix_positive` (conditioned on the
    honest opaque predicate `MainWindow`; since r305 PROVED in
    RH/Closure.lean from the two true holes `lstar_subordination` +
    `terminal_positive_main` via `lstar_terminal_implies_master`),
  * former edge A  => corollary `terminal_crossratio_main` (PROVED from
    the master theorem, RH/Closure.lean),
  * former edge B  => corollary `prefix_chain_positive_main` (PROVED
    from the master theorem, RH/Closure.lean).
This file keeps the `WindowLadder` bookkeeping structure (it correctly
records the measured structural facts and feeds the counterexample
guards) and the campaign kill lists.  It states no open theorem anymore.

KILL LIST (measured no-gos any future proof of edge A must avoid):
  - mass majorants (sup / Cauchy-Schwarz / triangle on masses): 0/42
    coverage, median x490 excess, break locus 0.97N            (r258)
  - bordered spectral bound f^T (I-Q)^{-1} f < 1: measured FALSE,
    g(1) = 6.06..15.41 on 42/42                                 (r257/r265)
  - naive s-coordinate: exact restatement of the endpoint       (r265)
  - IIKS s-family: definite dynamics, WRONG endpoint criterion
    (DEFINITENESS_WALL_EQUIVALENT)                              (r265)
  - eventual triangle (exceptions vanish for large N): refuted --
    the deepest rung kz52 (N = 878) is itself an exception      (r263)
  - any root-scale bound on the border drive |t_{N-2}|: fires the
    PAIRCORR detector (demand 0.66..2.31 dec)                   (r262/r263)

KILL LIST (edge B):
  - saddle dominance / level crossing: branches cross 5-16 degrees
    too early; MAIN crosses 45x without a flip                  (r259)
  - k-swap truncation / one-swap statistics: gap distributions
    fully overlap between MAIN and controls                     (r259)
  - subtraction-free positivity (LGV): no-go at depth k = 1     (r261)
  - global block involution: 4 unpaired negative configurations
    at n = 2 already                                            (r261)
  - transfer cone: best certificate IS the terminal budget
    positivity -- wall-equivalent, no new mechanism             (r261)

Claim boundary: research documentation.  NOT evidence for or against the
Riemann Hypothesis in either direction.  NO RH CLAIM.
-/
import RH.Basic
import RH.Recursion

namespace RH

/-- An abstract cofinal window ladder: for every rung `i` a chain datum
`data i` with terminal degree `N i`, a budget in the r243 form
`B = S_{N-1} + m` with one fixed positive floor scalar `m` (the corpus
imports `m = 5/7`, typed FLOOR_IMPORTED -- its source-pure derivation is
part of edge A), and terminal degrees growing cofinally. -/
structure WindowLadder (K : Type*) [Field K] [LinearOrder K]
    [IsStrictOrderedRing K] where
  data : ℕ → ChainData K
  N : ℕ → ℕ
  /-- the floor scalar (`5/7` in the corpus, FLOOR_IMPORTED per r258). -/
  m : K
  m_pos : 0 < m
  /-- terminal degrees are positive... -/
  N_pos : ∀ i, 0 < N i
  /-- ...and cofinal in the ladder. -/
  cofinal : ∀ n : ℕ, ∃ i, n ≤ N i
  /-- the r243 budget form `B = S_{N-1} + m` on every rung. -/
  budget_form : ∀ i, (data i).B = (data i).S (N i - 1) + m
  /-- the wall through the free window: the whole h-prefix is positive
  (the MAIN measurement; for the ladder this is itself only known
  finitely -- 42/42 rungs, N = 142..878). -/
  wall : ∀ i, ∀ k < N i, 0 < (data i).h k
  /-- the Chebyshev coefficients never vanish (by construction). -/
  cheb : ∀ i, ∀ k < N i, (data i).c k ≠ 0

variable {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- The terminal cross-ratio of rung `i`:
`q_N = rho_{N-1} / m` (r258 TERMINAL_Q_LAW). -/
def WindowLadder.q (L : WindowLadder K) (i : ℕ) : K :=
  (L.data i).rho (L.N i - 1) / L.m

/-!
## Where the two edges went (r273)

**EDGE A (fiber): the terminal cross-ratio inequality.**
Ledger: PRIME.PORT.COUPLEDTAU.TERMINAL_CROSSRATIO.01 [O].
Formerly stated here as `terminal_crossratio_cofinal : ∀ i, L.q i < 1`
over the bare ladder -- REFUTED as a universal by
`Counterexamples.terminal_crossratio_not_universal` (the model
`h = c = 1, F = 2, m = 1` satisfies every field and has `q_N = 4`).
The truth-capable form is `terminal_crossratio_main` in RH/Window.lean,
PROVED from the master theorem `augmented_prefix_positive`.
Measured support unchanged: 42/42 rungs, min margin 0.0139, trend FLAT
(r258); two-branch split r263: 35/42 rungs + both mains close by the
triangle certificate alone (`two_branch_cheap_strict`), the 7-rung
exception branch is the S6 requirement `|Z^RHP| + |E| < sqrt(m)`.

**EDGE B (base/orientation): globally oriented prefix resummation.**
Ledger: PRIME.PORT.FULLSOURCE.PREFIX_RESUMMATION.01 [O].
Formerly stated here as `prefix_resummation_positivity` over an
arbitrary source predicate -- REFUTED as a universal by
`Counterexamples.prefix_resummation_not_universal` (`MainSource := True`
and `h_0 = -1` give `tau_1 = -1`).  The truth-capable form is
`prefix_chain_positive_main` in RH/Window.lean (the wall derived from
the master positivity, not assumed); the missing DEFINITION of the main
source is exactly the opaque predicate `MainWindow` (r261: the object IS
the coherent sum; cancellation is ~5.7 decades deep at n = 6 already).

## r397 census (exact selected domain)

The C1 holes `lstar_canonical` / `terminal_q_canonical` are DEGRADED
to the alternative (rational-certificate) route: `CanonicalWindow`
may be empty after exact real transcription, and mesh-tolerance is
uncoupled from the L* margin.  They are kept as typed `sorry`s, not
deleted.  The load-bearing open kernel is the named Prop
`selected_augDualResolvent_gt_half` (`RH/Selected.lean`):
`∀ᶠ k, (R†(W^ℝ_k) − ½·1).PosDef` on the constructed real windows
of the selected sequence `a_k = 2^k`, `m_k = k·2^{⌊√k⌋}−1`.
New named Props (not sorrys): `SelectedMasterImpliesPlainReads`,
`ExactArchAgreesWithArchRead`.  Sorry census unchanged at 5.

## r417 census (source Schur sign)

The edge scalar `sch < 0` now has two closed formulae
(Woodbury `sch = den-2 + sᵀ(A₀+Ucd Ucdᵀ)⁻¹s` and the
unnormalized Sylvester chart).  Normalized `τ → 0` is a
census, not a theorem; the rate is `RATE_OPEN` (Uvarov / fold
identities do not force the observed slope).  Cofinal `sch < 0`
on the selected sequence is not proved and is not a new
`sorry` here -- it remains a named open of the selected
mincut `selected_augDualResolvent_gt_half`.
NO RH CLAIM.

## r418 census (phi_bb sign)

The border-border entry has the closed split
`φ_bb = (den-2) + sᵀ A₀⁻¹ s` and the r407 form
`A₀⁻¹ = 2(C+I)(C-I)⁻¹`.  Uniform `φ_bb < 0` on the
named living census is refuted (6/14 vacuous core
windows overflow).  Pole-dominance on living P1 is
refuted (3/28).  Two named mechanisms; selected
`a_k = 2^k` remains a negative census, not a new
`sorry`.  Mincut unchanged.
NO RH CLAIM.

## r419 census (vacuous overflow)

The r417/r418 tension is resolved as a named
reduction, not a theorem: H3 (vacuous family
empties) is refuted (VAC grows 3/14 → 7/14);
the razor-pole account of overflow is refuted
(s ⊥ C_min-mode); the six live by τ²>φ_bb at
finite N, but vacuous τ_un is the vanishing
one so a cofinal τ-floor fails.  Cofinal VAC
`sch < 0` reduces to `φ_bb < 0` on large-N
VAC (EXT/selected census).  No new `sorry`.
Mincut unchanged.
NO RH CLAIM.

## r420 census (c_J vs Sigma)

`den = 1 + ||b||²/B_w - v_t·s` is exact
(`B_w = S_{N-2}+5/7`).  `den<2` has an O(1) budget
gap, but the balance reserve `-φ_bb` is O(0.03) and
shrinks along `a_k=2^k`.  The naive occupied-max
`Σ`-bound is loose (200–24000×).  Cofinal
`c_J < -Σ` on large-N VAC is a census, not a new
`sorry`.  The whole edge is one `c_J`-vs-self-energy
balance.  Mincut unchanged.
NO RH CLAIM.

## r421 census (reserve limit)

Selected $a_k=2^k$ prefers a floor
`R_∞ ≈ 0.030` for `R = -φ_bb` (AIC vs `k`);
`R ~ c/log k` is killed.  `k=10` border-fails.
The T0–C_min linearisation is the r411
dictionary (ratio 1.0000 on vacuous).
Three saturation scales, not one exponent.
Cofinal `R_∞ > 0` is a census, not a new
`sorry`.  Mincut unchanged.
NO RH CLAIM.
-/

end RH
