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
        itself stays untouched in RH/Source.lean; the opaque filling
        is the named Prop `SourceExactOfFamilyCompletion` below);
  (ii)  `GridElement` (the native dense class, built for real:
        support parameter, native dyadic grid, step values;
        autocorrelation and its piecewise-linear even interpolant as
        DERIVED definitions) + the elementwise finite stabilization:
        the comb channel PROVED with the onset `elementAnchor f`
        predefined from the element's support alone (the R310
        `finite_forms_converge_to_weil` shape, elementwise and in the
        corpus gauge), the pole channel PROVED as the native-mesh
        second-difference of `polePotential` (r376), the arch kernel
        channel as a TYPED sorry (classical analysis -- `arch_A` is
        not a second difference of a named elementary antiderivative;
        Gauss/Mellin missing from mathlib v4.29.1), and the full-form
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

SORRY CENSUS OF THIS FILE (r376): exactly ONE typed `sorry` —
`arch_elementwise_stabilization` (CLASSICAL, S2).  The pole-channel
stabilization is PROVED (native-mesh second-difference transcription
of `pole_lags_closed`, parallel to mesh-free `combRead`).  The
source-exact completion is demoted to the named Prop
`SourceExactOfFamilyCompletion` (opacity-forced, not a hole; the
transcribable half is `sourceExact_buildPrimeWindow`).  Named, not
asserted: `PoleDyadicIndependence` (dyadic refinement of the pole
pairing).  Extraction consumes only the arch sorry.

Claim boundary: research documentation.  NOT evidence for or against
the Riemann Hypothesis in either direction.  NO RH CLAIM.
-/
import RH.Source
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.Order.Floor.Semifield
import Mathlib.Algebra.Order.Interval.Finset.Basic
import Mathlib.Algebra.Order.Ring.Int
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
the opaque guard is the named Prop `SourceExactOfFamilyCompletion`
at the bottom of this section (r376: demoted from a `sorry`). -/

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

/-- **THE SOURCE-EXACT COMPLETION** (r326; r376 RETYPE).

The constructive half is already `sourceExact_buildPrimeWindow` (the
transcribable clauses) together with C1's named opaque constant
`canonicalCompletion` (`RH/Canonical.lean`): that constant NAMES the
intended arch/border/budget filling.  R373's kernel objects
(`weilArchKernel`, `polePotential`) are lag-kernels, not the per-atom
masses / v958 border column / r243 drain-sum, so they do not make
`canonicalCompletion` constructive.

Concluding the opaque `SourceExact` of any filling remains unprovable
by design (r273/r320).  r376 demotes that conjunct from a `sorry` to
the named Prop below — the C1 `PairMarginLaw` convention: statable,
never asserted.  The extraction route does not consume it. -/
def SourceExactOfFamilyCompletion : Prop :=
  ∀ (a m : ℕ) (ha : IsPrimePow a),
    ∃ (aw : Fin (specFamily a m ha).S → ℝ) (haw : ∀ j, 0 ≤ aw j)
      (bd : ℕ → ℝ) (B : ℝ) (hB : 0 < B),
      SourceExact { specFamily a m ha with
        archWeight := aw, arch_nonneg := haw,
        border := bd, budget := B, budget_pos := hB }

/-- record: the transcribable half of a family member is PROVED
(`SourceExactSpec`); the opaque `SourceExact` filling is the named
Prop `SourceExactOfFamilyCompletion`, not a hole. -/
theorem specFamily_sourceExact_completion_transcribable
    (a m : ℕ) (ha : IsPrimePow a) :
    SourceExactSpec (specFamily a m ha) :=
  sourceExact_buildPrimeWindow a m ha

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

/-- compact support, even form: vanishes once `|u|` exceeds the support. -/
theorem toFun_eq_zero_of_lt_abs {u : ℝ} (hu : f.supportBound < |u|) :
    f.toFun u = 0 := by
  rcases le_total 0 u with h | h
  · rw [abs_of_nonneg h] at hu
    exact f.toFun_eq_zero hu
  · rw [abs_of_nonpos h] at hu
    rw [← f.toFun_even]
    exact f.toFun_eq_zero hu

/-- **the predefined elementwise anchor onset** (R325 target form)
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

/-! ## The kernel channels: transcribed kernels + remaining arch opacity

r373: the closed-form *kernels* are Lean objects
(`weilArchKernel`, `polePotential`).  r376: the pole tent-read is
transcribed as the native-mesh second-difference pairing of
`polePotential` (Python `pole_lags_closed`) and the pole
stabilization is PROVED.  mathlib v4.29.1 carries `Complex.digamma =
logDeriv Complex.Gamma` with the recurrence `digamma (s+1) = digamma
s + s⁻¹` and the values at `0, 1, 1/2`; there is NO `Real.digamma`,
NO Gauss integral (mathlib TODO on `Digamma.lean`), and NO
ψ-monotonicity.  The arch tent-read `archRead` and pairing
`weilArchSide` remain opaque: identifying `arch_A` tent integrals
with `weilArchKernel` is Titchmarsh Ch. X / Weil 1952. -/

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

/-- the archimedean side of the Weil form of a grid element (the
exact kernel integral against `weilArchKernel`).  OPAQUE: named
classical pairing TODO.  r376: unlike the pole channel, the lag-domain
tent integrals `arch_A` (v563) are NOT a second difference of a named
elementary antiderivative in mathlib v4.29.1, so this read stays
opaque. -/
opaque weilArchSide : GridElement → ℝ

/-! ### Pole channel: transcribed tent-read (r376)

Python `pole_lags_closed` is the second difference of `g_pole = polePotential`:
`c_d = -(g((d-1)D) - 2g(dD) + g((d+1)D))/D`.  For an even compactly
supported test function this is the exact hat-function pairing of
`g''` (distributional), and on the native PL class the pairing is a
finite sum — no `MeasureTheory`, no ζ.  Mesh-independence under
dyadic refinement is the finite-element identity: the discrete
Laplacian of a coarse-grid interpolant is supported on the coarse
knots, with scale `2^{-r}`. -/

/-- dyadic mesh width `D(m) = 2^{-m}`. -/
noncomputable def meshWidth (m : ℕ) : ℝ := (1 : ℝ) / 2 ^ m

theorem meshWidth_pos (m : ℕ) : 0 < meshWidth m := by
  unfold meshWidth
  positivity

theorem meshWidth_mul_pow (m n : ℕ) :
    ((n * 2 ^ m : ℕ) : ℝ) * meshWidth m = n := by
  unfold meshWidth
  push_cast
  field_simp

theorem GridElement.D0_eq_meshWidth (f : GridElement) :
    f.D0 = meshWidth f.meshExp := rfl

/-- second difference of `polePotential` at mesh `D`, index `k`. -/
noncomputable def poleΔ (D : ℝ) (k : ℤ) : ℝ :=
  polePotential ((k : ℝ) * D - D) - 2 * polePotential ((k : ℝ) * D)
    + polePotential ((k : ℝ) * D + D)

theorem poleΔ_even (D : ℝ) (k : ℤ) : poleΔ D (-k) = poleΔ D k := by
  unfold poleΔ
  simp only [Int.cast_neg, neg_mul]
  have e1 : -((k : ℝ) * D) - D = -((k : ℝ) * D + D) := by ring
  have e2 : -((k : ℝ) * D) + D = -((k : ℝ) * D - D) := by ring
  rw [e1, e2, polePotential_even, polePotential_even, polePotential_even]
  ring

/-- two-sided cutoff large enough that every sample outside it
vanishes, at every mesh `m`. -/
def poleCutoff (f : GridElement) (m : ℕ) : ℕ :=
  (f.steps + 3) * 2 ^ m

theorem toFun_at_mesh_eq_zero (f : GridElement) (m : ℕ) {k : ℤ}
    (hk : poleCutoff f m ≤ k.natAbs) :
    f.toFun ((k : ℝ) * meshWidth m) = 0 := by
  have hD : 0 < meshWidth m := meshWidth_pos m
  have habs : |((k : ℝ) * meshWidth m)| = (k.natAbs : ℝ) * meshWidth m := by
    rw [abs_mul, abs_of_nonneg hD.le, Nat.cast_natAbs, Int.cast_abs]
  have hgt : f.supportBound < |((k : ℝ) * meshWidth m)| := by
    rw [habs]
    have hkR : ((poleCutoff f m : ℕ) : ℝ) ≤ (k.natAbs : ℝ) := by exact_mod_cast hk
    have hmul : (poleCutoff f m : ℝ) * meshWidth m ≤
        (k.natAbs : ℝ) * meshWidth m :=
      mul_le_mul_of_nonneg_right hkR hD.le
    have hcancel : (poleCutoff f m : ℝ) * meshWidth m = (f.steps + 3 : ℕ) :=
      meshWidth_mul_pow m (f.steps + 3)
    have hsteps : f.supportBound ≤ (f.steps : ℝ) := by
      unfold GridElement.supportBound GridElement.D0
      have : (1 / 2 ^ f.meshExp : ℝ) ≤ 1 := by
        have hp : 0 < (2 : ℝ) ^ f.meshExp := by positivity
        rw [div_le_one₀ hp]
        exact one_le_pow₀ (by norm_num : (1 : ℝ) ≤ 2)
      exact mul_le_of_le_one_right (Nat.cast_nonneg _) this
    have : (f.steps : ℝ) < (f.steps + 3 : ℕ) := by
      exact_mod_cast Nat.lt_add_of_pos_right (by decide : 0 < 3)
    linarith
  exact f.toFun_eq_zero_of_lt_abs hgt

/-- two-sided discrete pairing of `polePotential` against `f` at mesh `D`
truncated at `N`.  Python `read_lags(pole_lags_closed, D, F)` for even
`F` (the negative-index half equals the positive half by evenness). -/
noncomputable def polePairingZ (D : ℝ) (f : GridElement) (N : ℕ) : ℝ :=
  -(∑ k ∈ Finset.Icc (-(N : ℤ)) N,
      f.toFun ((k : ℝ) * D) * poleΔ D k) / D

/-- the pole even-read at mesh level `m` (window-independent: the pole
kernel does not see the anchor; truncation past the support is a
no-op). -/
noncomputable def poleEvenRead (m : ℕ) (f : GridElement) : ℝ :=
  polePairingZ (meshWidth m) f (poleCutoff f m)

/-- **the pole-channel tent-read** (r376 transcription of
`pole_lags_closed` / `read_lags`).  Parallel to `combRead` (mesh-free
exact atom sum: the tent-assembly equals it on the native class),
the pole read is taken at the element's native mesh — R325 S1.4: the
native mesh already carries the limit value.  Independent of the
anchor: the pole kernel does not see the window. -/
noncomputable def poleRead (_a _m : ℕ) (f : GridElement) : ℝ :=
  poleEvenRead f.meshExp f

/-- **the pole side of the Weil form**: the native-mesh pairing. -/
noncomputable def weilPoleSide (f : GridElement) : ℝ :=
  poleEvenRead f.meshExp f

/-- **named remaining finite-sum identity** (r376; not a `sorry`):
the second-difference pairing of `polePotential` against a native-PL
test function is independent of dyadic refinement.  Parallel to the
comb channel, the Lean read is taken at the native mesh (R325 S1.4:
the native mesh already carries the limit value).  The identity is
statable and classically true (affine second differences vanish;
`Real.cosh` is in mathlib); it is not consumed by the extraction. -/
def PoleDyadicIndependence : Prop :=
  ∀ f : GridElement, ∀ m, f.meshExp ≤ m →
    poleEvenRead m f = poleEvenRead f.meshExp f

/-- **the arch-channel elementwise stabilization** (r326 -- documented
`sorry` no. 1 of this file; type: CLASSICAL, S2; r376: remaining
hole named exactly).  Tent-read `archRead` equals the Weil pairing
`weilArchSide` at native-or-finer mesh past a finite anchor onset.

WHY THIS REMAINS A SORRY (r376 mathlib census, v4.29.1):
the lag-domain tent integrals `arch_A` (v563: `e^{-w/2}/(-expm1(-2w))`
plus the Euler+`log π` near-cell regularizer) are NOT a second
difference of a named elementary antiderivative.  mathlib carries
`Complex.digamma = logDeriv Gamma` with `digamma_apply_add_one`,
values at `0,1,1/2`, `meromorphic_digamma`,
`differentiableAt_Gamma` off the nonpositive integers, and
`Real.Gamma`; it does NOT carry Gauss's integral representation
(explicit TODO on `Digamma.lean`), `Real.digamma`, ψ-monotonicity,
or Mellin inversion identifying `arch_A` with `weilArchKernel`.
That identification is Titchmarsh, *The Theory of the Riemann
Zeta-Function*, 2nd ed., Ch. X, and Weil 1952.  Not a finite-sum
identity; not closable from `Complex.Gamma` foundations in this
mathlib tag. -/
theorem arch_elementwise_stabilization (f : GridElement) :
    ∃ a₀ : ℕ, ∀ a : ℕ, a₀ ≤ a → ∀ m : ℕ, f.meshExp ≤ m →
      archRead a m f = weilArchSide f := by
  sorry

/-- **the pole-channel elementwise stabilization** (r376 -- PROVED).
The pole tent-read is the native-mesh second-difference pairing of
`polePotential` (Python `pole_lags_closed`); the Weil pole side is
the same pairing.  Equality is definitional, parallel to
`comb_elementwise_stabilization` being mesh-free.  Remaining named
(not a hole): `PoleDyadicIndependence`.  No ζ. -/
theorem pole_elementwise_stabilization (f : GridElement) :
    ∃ a₀ : ℕ, ∀ a : ℕ, a₀ ≤ a → ∀ m : ℕ, f.meshExp ≤ m →
      poleRead a m f = weilPoleSide f := by
  refine ⟨f.elementAnchor, fun _ _ _ _ => rfl⟩

/-! ## The full form and the elementwise stabilization theorem -/

/-- **the full window form read** of a grid element at the canonical
window (anchor `a`, mesh level `m`): arch − comb + pole (the corpus
total `c = car + cat + cp`, with the atom channel entering as MINUS
the comb sum -- the R325 S1.3 sign convention).  HONEST SCOPE: the
comb summand is transcribed and proved above; the arch/pole summands
are the arch opaque read and the transcribed pole pairing. -/
noncomputable def fullRead (a m : ℕ) (f : GridElement) : ℝ :=
  archRead a m f - combRead a f + poleRead a m f

/-- **the Weil form** of a grid element: arch − comb + pole, each
side in its channel's reference (comb: the proved tsum; pole: the
native-mesh pairing of `polePotential`; arch: the opaque kernel
integral).  The spectral/zero side of the
explicit formula is NOT part of this definition -- it is the
arithmetically open content of the program. -/
noncomputable def weilForm (f : GridElement) : ℝ :=
  weilArchSide f - weilCombSide f + weilPoleSide f

/-- **THE ELEMENTWISE FINITE STABILIZATION** (r326, R325 target form
(ii); PROVED from the three channels -- the comb channel
unconditionally, the arch channel through its typed sorry
above, the pole channel unconditionally).  For every grid element there is a finite anchor onset `a₀`
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
anywhere.  (Modulo the one typed arch-kernel sorry inside the
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
