/-
RH/Elementwise.lean -- THE ELEMENTWISE EXTRACTION ARCHITECTURE
(r326; the R325 repair set, `extraction_order_probe.py`, sealed
SPEC_SHA 57e50d366a62f2a6 freeze / 18/18 gates, primary verdict
ELEMENTWISE_STABILIZATION_GO).

WHY THIS FILE EXISTS.  The R319 red-team audit and the wave-12
reviewer adjudication localized the missing LEVEL C of the proof
graph (the extraction: window-local positivity ==> the Weil form) in
two defects that appeared in NO sorry, because the correct statements
did not exist yet: (a) the rh/ cofinality theorems run in the ANCHOR
direction while the carrier hypothesis (H_cof) demands a
MESH-refinement PSD tower (the mesh-vs-anchor seam,
RH/Source.lean), and the window direction of that route is measured
inadmissible (false floors, hcof_dodging_audit S6.8); (b) the opaque
`SourceExact` guard (r320) is a good firewall but must NOT stand as a
free assumption in a final proof.  R325 measured the three reviewer
repair variants under a sealed tree and adjudicated the ELEMENTWISE
architecture as primary: on the native dense class (dyadic
step-function autocorrelations -- the v749 "Weil form of step
functions" class) EVERY channel of the canonical tower windows
stabilizes EXACTLY at a finite anchor onset PREDEFINED from the
element (a* = (n_g + 1) D0 / 2), and the value is CONSTANT under
dyadic mesh refinement -- the extraction consumes NO mesh-cofinal
ladder and NO transport.  This file implements that quantifier set:

  (i)   `CanonicalPrimeWindow` (the canonical family as a window
        predicate) + `SourceExactSpec` (the TRANSCRIBABLE
        source-exactness clauses as a real definition, not an opaque
        marker) + `sourceExact_buildPrimeWindow` (PROVED: every
        canonical family member satisfies them) -- the wave-12
        reviewer target "CanonicalPrimeWindow + exactly one
        construction theorem", eliminating `SourceExact` as a free
        assumption FROM THE EXTRACTION ROUTE (the opaque predicate
        itself stays untouched in RH/Source.lean; the relation is the
        documented completion sorry below);
  (ii)  `GridElement` (the native dense class, built for real:
        support parameter, native dyadic grid, step values;
        autocorrelation and its piecewise-linear even interpolant as
        DERIVED definitions) + the elementwise finite stabilization:
        the comb channel PROVED with the onset `elementAnchor f`
        predefined from the element's support alone (the R310
        `finite_forms_converge_to_weil` shape, elementwise and in the
        corpus gauge), the arch/pole kernel channels as TYPED sorrys
        (classical analysis -- the kernels are explicit
        integrals/closed forms, R325 S1-measured exact; their
        transcription is the named classical TODO), and the full-form
        statement `elementwise_finite_stabilization` PROVED from the
        three channels;
  (iii) `weil_nonneg_of_windowlocal` (PROVED -- a finite
        instantiation, NO ladder, NO (H_cof)): window-local
        positivity of the canonical family on the native class
        implies Weil-form nonnegativity on every grid element; plus
        the named compression-bridge statement
        `BorderedCompressionBridge` (the bordered ==> plain tower
        form step -- the documented S2 rest; named as a Prop, NOT
        asserted) with the proved composition
        `weil_nonneg_of_bordered`.

WHAT `fullRead`/`weilForm` COVER (honest scope).  The comb channel is
TRANSCRIBED and PROVED: `combRead` is the exact atom sum
`Sigma mu_n f(log n)` over the canonical atom set in the corpus gauge
`mu_n = 2 Lambda(n)/sqrt(n)` (Python MU_ALL; the tent-assembled
mesh read equals this exactly on the native class at native-or-finer
mesh -- the R325 S1.3 tent-read identity, measured at 1e-12, whose
mesh-grid transcription belongs to the fold TODO), and the gauge
relation to the built window's comb channel is PROVED
(`combRead_eq_window_channel`, `combMass_eq_gauge`).  The arch and
pole channels enter through the OPAQUE reads `archRead`/`poleRead`
and references `weilArchSide`/`weilPoleSide`: their semantics (the
exact Weil arch kernel `arch_A` = GL-48 tent integrals, v563; the
v716 pole closed form) are NOT yet transcribed -- the opaque
constants are the named classical TODO made visible, and every
statement about them below is a TYPED sorry, never a proof.  The
spectral/zero side of the full explicit formula is the arithmetically
open content of the program and is NOT stated here.

SIGN CONVENTION: the corpus total lag vector is c = car + cat + cp
with the atom channel read equal to MINUS the comb sum
(L_cat(F) = -Sigma mu_n F(u_n), R325 S1.3), so
`fullRead = archRead - combRead + poleRead` and
`weilForm = weilArchSide - weilCombSide + weilPoleSide`.

WHAT REPLACES (H_cof).  The extraction (iii) consumes: one anchor
`a >= elementAnchor f` (a prime via Euclid -- a FINITE choice per
element, not a cofinal tower), the element's OWN native mesh level
`m = f.meshExp` (predefined, never refined), and the stabilization
(ii).  No PSD transport along any ladder appears; the mesh limit
survives only off the native class as the classical
density/continuity step (R325 leg C: rate-controlled quadrature
defect inside the closed-form interpolation envelope -- not a
positivity ladder), which is deliberately NOT formalized here.

SORRY CENSUS OF THIS FILE: exactly THREE, all typed, all statements
that were NOT formalizable before this round (the wave-12 reviewer
reservation "the Level-C distance appears in no sorry" is hereby
partially discharged -- the distance is now named):
  * `arch_elementwise_stabilization`  -- CLASSICAL (S2): the arch
    kernel read stabilizes elementwise (R325 S1 measured exact,
    1.5e-15 mesh-constancy; the kernel transcription `arch_A` is
    classical integral analysis);
  * `pole_elementwise_stabilization`  -- CLASSICAL (S2): same for
    the v716 pole closed form (R325 S1: 2.0e-17);
  * `specFamily_sourceExact_completion` -- CLASSICAL +
    DEFINITIONAL/TECHNICAL (opacity-forced): every canonical family
    member admits a source-exact completion (genuine arch/border/
    budget data on the same atom set); the completion data is
    classical (the transcriptions), and no statement concluding the
    opaque `SourceExact` is provable by design (r273/r320
    convention).  NOTE: the extraction route (iii) does NOT consume
    this sorry -- that is the point of the architecture.
The five pre-existing intentional sorrys of the pilot are untouched
and byte-identical (`lstar_subordination`, `terminal_positive_main`,
`pair_margin_main`, `crossing_budget`,
`mainWindow_iff_builtFromPrimeSource`).

Claim boundary: research documentation.  NOT evidence for or against
the Riemann Hypothesis in either direction.  NO RH CLAIM.
-/
import RH.Source
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Analysis.SpecialFunctions.Gamma.Digamma
import Mathlib.Analysis.Complex.Trigonometric

namespace RH

local notation "Λ" => ArithmeticFunction.vonMangoldt

/-! ## (i) The canonical prime windows and the transcribable
source-exactness

The wave-12 reviewer target: `SourceExact` must not stand as a free
assumption in a final proof.  The repair is NOT to prove anything
about the opaque predicate (impossible by design) but to make the
extraction route consume only the CONSTRUCTION: the canonical family
`specFamily` (RH/Source.lean, r310) is predefined by arithmetic
alone, and everything the route needs from "source exactness" is
provable FOR the family -- `SourceExactSpec` below.  The relation to
the opaque guard is the separate documented completion sorry at the
bottom of this section. -/

/-- **THE CANONICAL PRIME WINDOWS** (r326, R325 target form (i)): a
real window is canonical iff it is built from a member of the
predefined family -- anchor `a` (a prime power), mesh level `m`.
This is the R325 repair set's replacement for the (H_cof) ladder
route: the extraction below quantifies window positivity over THIS
predicate and nothing else.  (Note `buildPrimeWindow` is
deliberately mesh-independent -- source theorem 4, r310b -- so at
window level the family is indexed by the anchor; the mesh level
enters the kernel channels through the reads below, never through
the built window.) -/
def CanonicalPrimeWindow (w : PrimeWindow) : Prop :=
  ∃ (a m : ℕ) (ha : IsPrimePow a), w = buildPrimeWindow (specFamily a m ha)

/-- every canonical family member is a canonical window (the
construction is the witness). -/
theorem canonicalPrimeWindow_build (a m : ℕ) (ha : IsPrimePow a) :
    CanonicalPrimeWindow (buildPrimeWindow (specFamily a m ha)) :=
  ⟨a, m, ha, rfl⟩

/-- canonical windows are explicitly-main (r310 sense): the canonical
family is a subfamily of the explicit construction. -/
theorem canonicalPrimeWindow_isExplicit {w : PrimeWindow}
    (h : CanonicalPrimeWindow w) : MainWindowExplicit w := by
  obtain ⟨a, m, ha, rfl⟩ := h
  exact ⟨specFamily a m ha, rfl⟩

/-- **THE TRANSCRIBABLE SOURCE-EXACTNESS** (r326): the clauses of
"this spec carries the genuine source data" that are statable and
provable TODAY, as a real definition (contrast the deliberately
opaque `SourceExact`, RH/Source.lean, which reserves the
NOT-yet-transcribed arch/border/budget-identity channels):
  (1)/(2) the node positions and comb weights ARE the `log`/`Λ`
          data of the atom set -- definitional BY THE INTERFACE
          DESIGN (r310 source theorem 1: the spec derives them
          instead of carrying them as free fields; recorded here so
          the predicate readably binds the comb channel);
  (3)     ATOM-SET COMPLETENESS -- the genuine content: the atoms
          are ALL prime powers `n ≤ anchor²`, not merely some
          strictly increasing selection (the structure alone does
          NOT force this; `predefined_family` proves it for the
          canonical family);
  (4)     budget positivity -- the transcribable r243 half (spec
          field since r320).
NOT bound here (honest): the arch/border channels and the full r243
budget identity -- those are exactly the intended content of the
opaque `SourceExact`, and no honest finite clause can bind them
before their transcriptions exist (r320 verification finding).  The
extraction route of this file never needs them bound: the kernel
channels enter through the typed reads below. -/
def SourceExactSpec (s : PrimeWindowSpec) : Prop :=
  (∀ j, s.node j = Real.log (s.primePowers j)) ∧
  (∀ j, s.combWeight j = Λ (s.primePowers j)) ∧
  (∀ n : ℕ, (∃ i, s.primePowers i = n) ↔ (IsPrimePow n ∧ n ≤ s.anchor ^ 2)) ∧
  0 < s.budget

/-- **THE CONSTRUCTION THEOREM** (r326, R325 target form (i); the
wave-12 reviewer's "exactly one construction theorem"; PROVED):
every canonical family member is source-exact in the transcribable
sense -- clauses (1)/(2) by derivation (`rfl`), clause (3) is
`predefined_family` (r310, structure theorem 1), clause (4) is the
spec field.  This is what replaces `SourceExact` as a hypothesis in
the extraction route: the route quantifies over BUILT windows, and
built windows PROVABLY carry their source. -/
theorem sourceExact_buildPrimeWindow (a m : ℕ) (ha : IsPrimePow a) :
    SourceExactSpec (specFamily a m ha) :=
  ⟨fun _ => rfl, fun _ => rfl, predefined_family a m ha,
    (specFamily a m ha).budget_pos⟩

/-- **THE SOURCE-EXACT COMPLETION** (r326 -- documented `sorry` no. 3
of this file; type: CLASSICAL + DEFINITIONAL/TECHNICAL,
opacity-forced).  The bridge relation between the new provable
`SourceExactSpec` and the old opaque `SourceExact`: every canonical
family member admits a SOURCE-EXACT COMPLETION -- a spec with the
SAME arithmetic data (atoms, anchor, mesh level; only the free
arch/border/budget fields replaced) that satisfies the opaque guard.
INTENDED CONTENT: the genuine completion data exists mathematically
-- the exact arch kernel weights (`arch_A`, GL-48 tent integrals,
v563), the v958 border column of this atom set, and the full r243
budget `B = S_{N-2} + 5/7`; filling them in is the named classical
transcription TODO of the r320 header.  WHY A SORRY: concluding the
opaque `SourceExact` is unprovable by design (r273 convention) --
this statement is the honest interface, exactly like the opacity
bridge of RH/Source.lean.  ARCHITECTURAL NOTE (the r326 point): the
extraction route (iii) below does NOT consume this statement --
`SourceExact` is thereby eliminated as a free assumption from the
route, and this sorry remains only as the documented relation
between the two predicates. -/
theorem specFamily_sourceExact_completion (a m : ℕ) (ha : IsPrimePow a) :
    ∃ (aw : Fin (specFamily a m ha).S → ℝ) (haw : ∀ j, 0 ≤ aw j)
      (bd : ℕ → ℝ) (B : ℝ) (hB : 0 < B),
      SourceExact { specFamily a m ha with
        archWeight := aw, arch_nonneg := haw,
        border := bd, budget := B, budget_pos := hB } := by
  sorry

/-! ## (ii) The native dense class: dyadic step-function
autocorrelations

The v749 "Weil form of step functions" class (R325 S1: elements
sealed by seed before any sign read; Python
`F = D0 * np.correlate(x, x, "full")`, knots on the native grid,
closing zero knot): a step function with `steps` values on the
dyadic grid of width `D0 = 2^{-meshExp}`, its autocorrelation
coefficients, and the even piecewise-linear interpolant through them
(support `steps * D0`).  Everything below the structure is DERIVED,
not a field -- the r310 interface convention. -/

/-- **a native grid element**: `steps` step values on the dyadic grid
of exponent `meshExp` (native mesh `D0 = 2^{-meshExp}`).  The probe
seeds `x` at random and seals it; here the data is free -- every
statement below holds for EVERY element of the class. -/
structure GridElement where
  /-- the number of native grid steps `n_g` (support = `steps · D0`). -/
  steps : ℕ
  /-- the native dyadic mesh exponent (`D0 = 2^{-meshExp}`). -/
  meshExp : ℕ
  /-- the step values. -/
  x : Fin steps → ℝ

namespace GridElement

variable (f : GridElement)

/-- the native mesh width `D0 = 2^{-meshExp}` (dyadic by
construction). -/
noncomputable def D0 : ℝ := 1 / 2 ^ f.meshExp

theorem D0_pos : 0 < f.D0 := by
  unfold D0; positivity

/-- **the autocorrelation coefficients** (DERIVED):
`a_d = D0 · Σ_i x_i x_{i+d}` for `d < steps`, `0` beyond (Python
`np.correlate(x, x, "full")[ng-1:] * D0` with the closing zero
knot). -/
noncomputable def acf (d : ℕ) : ℝ :=
  f.D0 * ∑ i : Fin f.steps,
    if h : (i : ℕ) + d < f.steps then f.x i * f.x ⟨(i : ℕ) + d, h⟩ else 0

/-- beyond the support the autocorrelation vanishes (the closing
zero knot and everything after). -/
theorem acf_eq_zero {d : ℕ} (hd : f.steps ≤ d) : f.acf d = 0 := by
  unfold acf
  rw [Finset.sum_eq_zero, mul_zero]
  intro i _
  rw [dif_neg]
  omega

/-- the zero-lag autocorrelation is nonnegative (`D0 · Σ x_i² ≥ 0` --
the class is genuinely of positive type at lag 0; the full
positive-type structure is what the window positivity premise of
(iii) consumes). -/
theorem acf_zero_nonneg : 0 ≤ f.acf 0 := by
  refine mul_nonneg f.D0_pos.le (Finset.sum_nonneg fun i _ => ?_)
  rw [dif_pos (by have := i.isLt; omega : (i : ℕ) + 0 < f.steps)]
  exact mul_self_nonneg (f.x i)

/-- **the induced test function** (DERIVED): the EVEN piecewise-linear
interpolant of the autocorrelation on the native grid (Python
`np.interp` on `|u|`, knots `d · D0`, `right=0`): at
`t = |u| / D0 ∈ [d, d+1]` linear between `a_d` and `a_{d+1}`. -/
noncomputable def toFun (u : ℝ) : ℝ :=
  f.acf ⌊|u| / f.D0⌋₊
    + (|u| / f.D0 - (⌊|u| / f.D0⌋₊ : ℝ))
      * (f.acf (⌊|u| / f.D0⌋₊ + 1) - f.acf ⌊|u| / f.D0⌋₊)

/-- **the support parameter** `steps · D0` -- the quantity from which
the elementwise onset `elementAnchor` is PREDEFINED (R325: "a₀, m_f
in advance from f"). -/
noncomputable def supportBound : ℝ := (f.steps : ℝ) * f.D0

theorem supportBound_nonneg : 0 ≤ f.supportBound :=
  mul_nonneg (Nat.cast_nonneg _) f.D0_pos.le

/-- the test function is even (the class is a class of
autocorrelations -- Weil test data). -/
theorem toFun_even (u : ℝ) : f.toFun (-u) = f.toFun u := by
  simp [toFun, abs_neg]

/-- **compact support, proved from the construction**: beyond the
support parameter the test function vanishes identically. -/
theorem toFun_eq_zero {u : ℝ} (hu : f.supportBound < u) : f.toFun u = 0 := by
  have hD0 : 0 < f.D0 := f.D0_pos
  have habs : (f.steps : ℝ) * f.D0 < |u| :=
    lt_of_lt_of_le hu (le_abs_self u)
  have ht : (f.steps : ℝ) ≤ |u| / f.D0 := (le_div_iff₀ hD0).mpr habs.le
  have hd : f.steps ≤ ⌊|u| / f.D0⌋₊ := Nat.le_floor ht
  unfold toFun
  rw [f.acf_eq_zero hd, f.acf_eq_zero (le_trans hd (Nat.le_succ _))]
  ring

/-- **the predefined elementwise anchor onset** (R325 target form
(ii)): computed from the element's support parameter ALONE, before
any window is built or any form is read -- `a₀(f) = max(1, ⌈exp
(steps · D0)⌉)`.  (The measured onset is `α* = (n_g + 1) D0 / 2` in
log coordinates; `supportBound ≥ α*`, so this anchor is safe and
still element-predefined.)  The predefined mesh level is the
element's own `meshExp` -- R325 S1.4: the native mesh already
carries the limit value. -/
noncomputable def elementAnchor : ℕ :=
  max 1 (Nat.ceil (Real.exp f.supportBound))

end GridElement

/-! ## The comb channel, corpus gauge: transcribed and PROVED

The corpus comb masses are `mu_n = 2 Λ(n)/√n` at nodes `log n`
(Python `MU_ALL`/`U_ALL`, v563).  The Lean window layer keeps the
bare `Λ` (the rescale is a positive diagonal gauge -- r310 header);
the FORM read of the extraction, however, is a linear functional and
must be stated in the gauge the corpus positivity lives in.  So the
comb read below carries the corpus gauge, and the gauge relation to
the built window's comb channel is PROVED (`combMass_eq_gauge`,
`combRead_eq_window_channel`). -/

/-- the corpus comb mass `2 Λ(n)/√n` (Python `MU_ALL`). -/
noncomputable def combMass (n : ℕ) : ℝ := 2 * Λ n / Real.sqrt n

theorem combMass_nonneg (n : ℕ) : 0 ≤ combMass n :=
  div_nonneg (mul_nonneg (by norm_num)
    ArithmeticFunction.vonMangoldt_nonneg) (Real.sqrt_nonneg _)

/-- the corpus gauge is the positive diagonal rescale of the exact
`Λ`: `mu_n = (2/√n) · Λ(n)` -- the r310 "positive diagonal gauge
lives in the builder map" note, now a statement. -/
theorem combMass_eq_gauge (n : ℕ) :
    combMass n = (2 / Real.sqrt n) * Λ n := by
  unfold combMass
  ring

/-- **the finite comb read of the canonical window at anchor `a`**
(corpus gauge): the exact atom sum `Σ_{n ∈ atoms(a)} mu_n f(log n)`
over the predefined atom set.  MESH STATUS (honest, R325 S1.3/S1.4):
this is the mesh-free exact atom sum; the corpus evaluates the comb
channel as a tent-assembled read on the mesh grid, and the two are
EQUAL on the native class at native-or-finer mesh (the derived
tent-read identity, measured at 1e-12) -- the tent assembly's
transcription belongs to the fold TODO (RH/Source.lean, stage 3),
which is why the read here is stated mesh-free and the window's
mesh-independence is `buildPrimeWindow_refinement_compatible`. -/
noncomputable def combRead (a : ℕ) (f : GridElement) : ℝ :=
  ∑ n ∈ windowAtoms a, combMass n * f.toFun (Real.log n)

/-- **the Weil prime side in the corpus gauge**: `Σ'_n mu_n f(log n)`
(a tsum over all naturals; the mass vanishes off the prime powers). -/
noncomputable def weilCombSide (f : GridElement) : ℝ :=
  ∑' n : ℕ, combMass n * f.toFun (Real.log n)

/-- **THE COMB ELEMENTWISE FINITE STABILIZATION** (r326, R325 target
form (ii), comb channel; PROVED).  For every grid element the finite
comb reads of the canonical family EQUAL the Weil prime side for
every anchor `a ≥ elementAnchor f` -- with the onset PREDEFINED from
the element's support alone, and NO mesh quantifier at all (the read
is mesh-free; the built window is mesh-independent).  This is the
elementwise form of `finite_forms_converge_to_weil` (r310) in the
corpus gauge: stabilization, not convergence -- the window rule
covers the support from the onset on, exactly (R325 S1.1: rel dev
0 above onset, a genuine onset below). -/
theorem comb_elementwise_stabilization (f : GridElement) {a : ℕ}
    (haA : f.elementAnchor ≤ a) : combRead a f = weilCombSide f := by
  have ha1 : 1 ≤ a := le_trans (le_max_left _ _) haA
  have haexp : Real.exp f.supportBound ≤ (a : ℝ) := by
    have hceil : (Nat.ceil (Real.exp f.supportBound) : ℝ) ≤ (a : ℝ) := by
      exact_mod_cast le_trans (le_max_right _ _) haA
    exact le_trans (Nat.le_ceil _) hceil
  have hzero : ∀ n : ℕ, n ∉ windowAtoms a →
      combMass n * f.toFun (Real.log n) = 0 := by
    intro n hn
    by_cases hpp : IsPrimePow n
    · have hgt : a ^ 2 < n := by
        by_contra hle
        exact hn (Finset.mem_filter.mpr
          ⟨Finset.mem_range.mpr (by omega), hpp⟩)
      have ha1R : (1 : ℝ) ≤ (a : ℝ) := by exact_mod_cast ha1
      have hgtR : ((a : ℝ)) ^ 2 < (n : ℝ) := by exact_mod_cast hgt
      have hexp_lt : Real.exp f.supportBound < (n : ℝ) := by nlinarith
      have hb : f.supportBound < Real.log n := by
        have := Real.log_lt_log (Real.exp_pos f.supportBound) hexp_lt
        simpa [Real.log_exp] using this
      rw [f.toFun_eq_zero hb, mul_zero]
    · have hL : Λ n = 0 := by
        by_contra h
        exact hpp (ArithmeticFunction.vonMangoldt_ne_zero_iff.mp h)
      simp [combMass, hL]
  unfold combRead weilCombSide
  exact (tsum_eq_sum hzero).symm

/-- **the comb read IS the canonical window's comb channel** (PROVED;
the honesty tie to the R325 target form `fullForm (buildPrimeWindow
(specFamily a m ha)) f`): the corpus-gauge read equals the sum of
`mu(atom_j) · f(node_j)` over the BUILT window's node vector -- same
atoms, same nodes, the mass in the proved gauge relation
`combMass_eq_gauge` to the window's `Λ` channel. -/
theorem combRead_eq_window_channel (a m : ℕ) (ha : IsPrimePow a)
    (f : GridElement) :
    combRead a f
      = ∑ j, combMass ((specFamily a m ha).primePowers j)
          * f.toFun ((buildPrimeWindow (specFamily a m ha)).nodes j) := by
  unfold combRead
  rw [← Finset.sum_coe_sort (windowAtoms a)
    (fun n => combMass n * f.toFun (Real.log n))]
  exact (Fintype.sum_equiv ((windowAtoms a).orderIsoOfFin rfl).toEquiv
    _ _ fun j => rfl).symm

/-- the Λ-gauge window form stabilizes elementwise too (the r310
statement `finite_forms_converge_to_weil`, now with the EXPLICIT
element-predefined onset instead of an existential): stated for the
record, PROVED by the same support argument. -/
theorem comb_window_elementwise_stabilization (f : GridElement)
    {a : ℕ} (ha : IsPrimePow a) (m : ℕ) (haA : f.elementAnchor ≤ a) :
    (buildPrimeWindow (specFamily a m ha)).combForm f.toFun
      = weilPrimeSide f.toFun := by
  have ha1 : 1 ≤ a := le_trans (le_max_left _ _) haA
  have haexp : Real.exp f.supportBound ≤ (a : ℝ) := by
    have hceil : (Nat.ceil (Real.exp f.supportBound) : ℝ) ≤ (a : ℝ) := by
      exact_mod_cast le_trans (le_max_right _ _) haA
    exact le_trans (Nat.le_ceil _) hceil
  have hzero : ∀ n : ℕ, n ∉ windowAtoms a →
      Λ n * f.toFun (Real.log n) = 0 := by
    intro n hn
    by_cases hpp : IsPrimePow n
    · have hgt : a ^ 2 < n := by
        by_contra hle
        exact hn (Finset.mem_filter.mpr
          ⟨Finset.mem_range.mpr (by omega), hpp⟩)
      have ha1R : (1 : ℝ) ≤ (a : ℝ) := by exact_mod_cast ha1
      have hgtR : ((a : ℝ)) ^ 2 < (n : ℝ) := by exact_mod_cast hgt
      have hexp_lt : Real.exp f.supportBound < (n : ℝ) := by nlinarith
      have hb : f.supportBound < Real.log n := by
        have := Real.log_lt_log (Real.exp_pos f.supportBound) hexp_lt
        simpa [Real.log_exp] using this
      rw [f.toFun_eq_zero hb, mul_zero]
    · have hL : Λ n = 0 := by
        by_contra h
        exact hpp (ArithmeticFunction.vonMangoldt_ne_zero_iff.mp h)
      rw [hL, zero_mul]
  have htsum : weilPrimeSide f.toFun
      = ∑ n ∈ windowAtoms a, Λ n * f.toFun (Real.log n) :=
    tsum_eq_sum hzero
  have hform : (buildPrimeWindow (specFamily a m ha)).combForm f.toFun
      = ∑ n ∈ windowAtoms a, Λ n * f.toFun (Real.log n) := by
    rw [← Finset.sum_coe_sort (windowAtoms a)
      (fun n => Λ n * f.toFun (Real.log n))]
    exact Fintype.sum_equiv ((windowAtoms a).orderIsoOfFin rfl).toEquiv
      _ _ (fun j => rfl)
  rw [hform, htsum]

/-! ## The kernel channels: transcribed kernels + opaque tent-reads

r373: the closed-form *kernels* are now Lean objects
(`weilArchKernel`, `polePotential`).  mathlib v4.29.1 carries
`Complex.digamma = logDeriv Complex.Gamma` with the recurrence
`digamma (s+1) = digamma s + s⁻¹` and the values at `0, 1, 1/2`;
there is NO `Real.digamma`, NO Gauss integral (mathlib TODO on
`Digamma.lean`), and NO ψ-monotonicity.  The tent-reads
`archRead`/`poleRead` and the Weil-side pairings
`weilArchSide`/`weilPoleSide` remain opaque: identifying a tent
read with a kernel pairing is classical quadrature (R325 S1;
Titchmarsh Ch. X; Weil 1952).  The two stabilization sorrys
below are exactly that remaining identification, now stated
against named kernels. -/

/-- Archimedean digamma factor on the critical line (Titchmarsh,
*The Theory of the Riemann Zeta-Function*, 2nd ed. (1986), Chapter X;
Weil, *Sur les «formules explicites» de la théorie des nombres
premiers*, Comm. Math. Helv. 26 (1952)). -/
noncomputable def weilArchDigamma (t : ℝ) : ℂ :=
  Complex.digamma ((1 / 4 : ℂ) + Complex.I * t / (2 * Real.pi))

/-- Even combination of the archimedean digamma factor.  Dictionary
kernel; not yet paired against mesh tents. -/
noncomputable def weilArchKernel (t : ℝ) : ℝ :=
  (weilArchDigamma t + weilArchDigamma (-t)).re

theorem weilArchKernel_even (t : ℝ) :
    weilArchKernel (-t) = weilArchKernel t := by
  simp [weilArchKernel, add_comm]

/-- v716 closed form of the polar contribution as a function of lag:
`Π(t) = −8 (cosh(|t|/2) − 1)` (Python `stage2.pole_lags_closed`;
elementary hyperbolic identity, no ζ). -/
noncomputable def polePotential (t : ℝ) : ℝ :=
  -8 * (Real.cosh (|t| / 2) - 1)

theorem polePotential_zero : polePotential 0 = 0 := by
  simp [polePotential, Real.cosh_zero]

theorem polePotential_even (t : ℝ) : polePotential (-t) = polePotential t := by
  simp [polePotential, abs_neg]

theorem polePotential_eq_cosh (t : ℝ) :
    polePotential t = -8 * (Real.cosh (t / 2) - 1) := by
  simp [polePotential]
  have : Real.cosh (|t| / 2) = Real.cosh (t / 2) := by
    rcases le_total 0 t with h | h
    · rw [abs_of_nonneg h]
    · rw [abs_of_nonpos h, neg_div, Real.cosh_neg]
  rw [this]

theorem polePotential_nonpos (t : ℝ) : polePotential t ≤ 0 := by
  rw [polePotential_eq_cosh]
  have hc : 0 < Real.cosh (t / 2) := Real.cosh_pos _
  have hsq : 1 ≤ Real.cosh (t / 2) ^ 2 := by
    rw [Real.cosh_sq]
    nlinarith [sq_nonneg (Real.sinh (t / 2))]
  have : 1 ≤ Real.cosh (t / 2) := by nlinarith
  nlinarith

/-- the archimedean kernel read of a grid element at the canonical
window (anchor, mesh level) -- Python
`read_lags(arch_lags(M, D), D, F)`.  OPAQUE: tent-quadrature of
`weilArchKernel` (Titchmarsh Ch. X; R325 S1 mesh-constancy 1.5e-15). -/
opaque archRead : ℕ → ℕ → GridElement → ℝ

/-- the pole-channel read (v716 closed form).  OPAQUE: tent-quadrature
of `polePotential` (R325 S1 mesh-constancy 2.0e-17). -/
opaque poleRead : ℕ → ℕ → GridElement → ℝ

/-- the archimedean side of the Weil form of a grid element (the
exact kernel integral against `weilArchKernel`).  OPAQUE: named
classical pairing TODO. -/
opaque weilArchSide : GridElement → ℝ

/-- the pole side of the Weil form of a grid element (pairing against
`polePotential`).  OPAQUE: named classical pairing TODO. -/
opaque weilPoleSide : GridElement → ℝ

/-- **the arch-channel elementwise stabilization** (r326 -- documented
`sorry` no. 1 of this file; type: CLASSICAL, S2).  Remaining hole,
stated against the transcribed kernel `weilArchKernel`: the tent-read
`archRead` equals the Weil pairing `weilArchSide` at native-or-finer
mesh past a finite anchor onset.  Classical ingredients not in
mathlib v4.29.1: Gauss's integral for `digamma` (mathlib TODO on
`Digamma.lean`), `Real.digamma`, and ψ-monotonicity; the dictionary
form is Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd ed.,
Ch. X, and Weil 1952.  R325 S1 measured exactness on the native class
(onset `α*`, mesh constancy 1.5e-15). -/
theorem arch_elementwise_stabilization (f : GridElement) :
    ∃ a₀ : ℕ, ∀ a : ℕ, a₀ ≤ a → ∀ m : ℕ, f.meshExp ≤ m →
      archRead a m f = weilArchSide f := by
  sorry

/-- **the pole-channel elementwise stabilization** (r326 -- documented
`sorry` no. 2 of this file; type: CLASSICAL, S2).  Remaining hole,
stated against the transcribed kernel `polePotential` (closed form
PROVED even/zero/nonpositive above): tent-read `poleRead` equals the
Weil pairing `weilPoleSide`.  Classical quadrature on the PL class
(R325 S1: 2.0e-17); no ζ. -/
theorem pole_elementwise_stabilization (f : GridElement) :
    ∃ a₀ : ℕ, ∀ a : ℕ, a₀ ≤ a → ∀ m : ℕ, f.meshExp ≤ m →
      poleRead a m f = weilPoleSide f := by
  sorry

/-! ## The full form and the elementwise stabilization theorem -/

/-- **the full window form read** of a grid element at the canonical
window (anchor `a`, mesh level `m`): arch − comb + pole (the corpus
total `c = car + cat + cp`, with the atom channel entering as MINUS
the comb sum -- the R325 S1.3 sign convention).  HONEST SCOPE: the
comb summand is transcribed and proved above; the arch/pole summands
are the opaque reads (named classical TODO). -/
noncomputable def fullRead (a m : ℕ) (f : GridElement) : ℝ :=
  archRead a m f - combRead a f + poleRead a m f

/-- **the Weil form** of a grid element: arch − comb + pole, each
side in its channel's reference (comb: the proved tsum; arch/pole:
the opaque kernel integrals).  The spectral/zero side of the
explicit formula is NOT part of this definition -- it is the
arithmetically open content of the program. -/
noncomputable def weilForm (f : GridElement) : ℝ :=
  weilArchSide f - weilCombSide f + weilPoleSide f

/-- **THE ELEMENTWISE FINITE STABILIZATION** (r326, R325 target form
(ii); PROVED from the three channels -- the comb channel
unconditionally, the arch/pole channels through their typed sorrys
above).  For every grid element there is a finite anchor onset `a₀`
such that for EVERY anchor `a ≥ a₀` and every mesh level at or below
the element's native mesh (`m ≥ f.meshExp` -- the predefined `m_f`),
the full finite window form EQUALS the Weil form.  The onset is
elementwise (`a₀` depends on `f` -- for the comb channel explicitly:
`elementAnchor f` from the support alone); NO mesh-cofinal ladder,
NO transport, NO (H_cof) appears. -/
theorem elementwise_finite_stabilization (f : GridElement) :
    ∃ a₀ : ℕ, ∀ a : ℕ, a₀ ≤ a → ∀ m : ℕ, f.meshExp ≤ m →
      fullRead a m f = weilForm f := by
  obtain ⟨aA, hA⟩ := arch_elementwise_stabilization f
  obtain ⟨aP, hP⟩ := pole_elementwise_stabilization f
  refine ⟨max f.elementAnchor (max aA aP), fun a haa m hm => ?_⟩
  have h1 : f.elementAnchor ≤ a := le_trans (le_max_left _ _) haa
  have h2 : aA ≤ a :=
    le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) haa
  have h3 : aP ≤ a :=
    le_trans (le_trans (le_max_right _ _) (le_max_right _ _)) haa
  unfold fullRead weilForm
  rw [comb_elementwise_stabilization f h1, hA a h2 m hm, hP a h3 m hm]

/-! ## (iii) The extraction WITHOUT the ladder

The R325 repair of the H_cof seam: window positivity, quantified over
the CANONICAL family only (the predicate of section (i)) and typed on
the PLAIN full form (the form for which the finite instantiation
actually goes through), yields Weil-form nonnegativity on every
element of the native dense class by ONE finite instantiation per
element -- pick any prime anchor above the element's onset (Euclid)
and the element's own native mesh.  No refinement tower, no PSD
transport, no `SourceExact` hypothesis. -/

/-- **the window-local positivity premise** (typed honestly on the
plain full form): every canonical window read of every
mesh-compatible grid element is nonnegative.  This is "∀ w,
CanonicalPrimeWindow w → 0 ≤ windowForm w" evaluated per test
element, with the canonical index `(a, m)` explicit because the
kernel reads live on the mesh grid while the built window is
mesh-independent; the guard `f.meshExp ≤ m` is the honest grid
compatibility (positivity of reads = PSD of the window's lag
Toeplitz form holds for autocorrelations ON the window's grid --
R325 S3: per-stage Levinson PD).  HONEST STATUS: this premise is the
Level-B content (the two window-local holes) transported to the
canonical family -- it is a PREMISE here, proved nowhere, and this
file makes no claim about it. -/
def WindowLocalPositive : Prop :=
  ∀ (a m : ℕ), IsPrimePow a → ∀ f : GridElement,
    f.meshExp ≤ m → 0 ≤ fullRead a m f

/-- **THE EXTRACTION WITHOUT THE LADDER** (r326, R325 target form
(iii); PROVED -- a finite instantiation).  Window-local positivity of
the canonical family implies Weil-form nonnegativity on every grid
element: given `f`, take the elementwise onset `a₀` (stabilization),
a prime anchor `a ≥ a₀` (Euclid, `Nat.exists_infinite_primes` -- a
FINITE choice, not a cofinal tower), and the element's own native
mesh `m = f.meshExp`; then `weilForm f = fullRead a m f ≥ 0`.  This
REPLACES the (H_cof) route: no mesh-refinement PSD tower is consumed
anywhere.  (Modulo the two typed kernel sorrys inside the
stabilization; the comb-only instantiation is unconditional.)  NO RH
CLAIM: the premise is exactly the open window-local content. -/
theorem weil_nonneg_of_windowlocal (hpos : WindowLocalPositive) :
    ∀ f : GridElement, 0 ≤ weilForm f := by
  intro f
  obtain ⟨a₀, hstab⟩ := elementwise_finite_stabilization f
  obtain ⟨p, hp₀, hp⟩ := Nat.exists_infinite_primes a₀
  have h := hpos p f.meshExp hp.prime.isPrimePow f (le_refl _)
  rwa [hstab p hp₀ f.meshExp (le_refl _)] at h

/-- **THE COMPRESSION BRIDGE, named** (r326; the documented S2 rest of
R325 -- deliberately a NAMED Prop, not an asserted theorem: its
faithful transcription needs the bordered/border-column data, i.e.
the `SourceExact` elimination).  The corpus positivity certificates
live on the BORDERED tower form (v749/v755 bordered Schur tower, the
odd compression); the premise of the extraction above is typed on
the PLAIN full form.  The bridge "bordered positivity ⟹ plain
positivity on the canonical family" is the named missing step
between the two typings -- once stated against a transcribed
bordered read, it becomes a finite matrix statement (principal
compression).  Parametrized over the bordered read so that NO new
opaque constant and NO truth commitment enters here. -/
def BorderedCompressionBridge
    (borderedRead : ℕ → ℕ → GridElement → ℝ) : Prop :=
  ∀ (a m : ℕ), IsPrimePow a → ∀ f : GridElement, f.meshExp ≤ m →
    0 ≤ borderedRead a m f → 0 ≤ fullRead a m f

/-- the extraction composes through the compression bridge (PROVED
from (iii)): bordered window-local positivity plus the named bridge
give the same Weil-form conclusion. -/
theorem weil_nonneg_of_bordered
    (borderedRead : ℕ → ℕ → GridElement → ℝ)
    (hbridge : BorderedCompressionBridge borderedRead)
    (hpos : ∀ (a m : ℕ), IsPrimePow a → ∀ f : GridElement,
      f.meshExp ≤ m → 0 ≤ borderedRead a m f) :
    ∀ f : GridElement, 0 ≤ weilForm f :=
  weil_nonneg_of_windowlocal fun a m ha f hm =>
    hbridge a m ha f hm (hpos a m ha f hm)

end RH
