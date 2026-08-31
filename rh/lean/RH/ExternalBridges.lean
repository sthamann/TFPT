/-
RH/ExternalBridges.lean -- THE THREE EXTERNAL BRIDGES MISSING AT r462
(r463, PRIME.RDAGGER.LEAN_FIDELITY_REPAIR.01).

These declarations deliberately increase the sorry census.  They
replace three absent arrows by typed, auditable obligations:

  internal GridElement positivity
    -> full Weil test-class positivity
    -> standard explicit-formula positivity
    -> Mathlib.NumberTheory.LSeries.RiemannZeta.RiemannHypothesis.

Nothing in this file is an RH claim.  The endpoint in
`FrequentlySelected.lean` is named `_internal` until all three arrows
are proved without `sorry`.
-/
import RH.Elementwise
import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.Topology.Order.Basic

namespace RH

/-- Carrier for the fixed-support real Weil autocorrelation class.

The test is continuous, even, compactly supported in
`[-supportRadius, supportRadius]`, and carries an explicit real
autorrelation witness `g(u)=∫ h(t)h(t+u)dt`.  Further strip/transform
conditions needed by a chosen Guinand--Weil theorem remain in the
`admissible` field; they are no longer conflated with density. -/
structure FullWeilTest where
  toFun : ℝ → ℝ
  supportRadius : ℝ
  supportRadius_nonneg : 0 ≤ supportRadius
  continuous_toFun : Continuous toFun
  even_toFun : Function.Even toFun
  support_toFun : ∀ u : ℝ, supportRadius < |u| → toFun u = 0
  autocorrelation : ∃ h : ℝ → ℝ, MeasureTheory.Integrable h ∧
    ∀ u : ℝ, toFun u =
      ∫ t : ℝ, h t * h (t + u) ∂MeasureTheory.volume
  admissible : Prop

/-- Full-class archimedean channel.  Its concrete integral
normalization is the r475 `weilArchSide` formula with `F.toFun`;
factoring it here lets density be audited channel by channel. -/
opaque fullWeilArchSide : FullWeilTest → ℝ

/-- Full-class prime channel.  Fixed support makes the von-Mangoldt
pairing a finite sum, although the finite-support reduction is one of
the continuity lemmas below. -/
opaque fullWeilCombSide : FullWeilTest → ℝ

/-- Full-class pole channel.  This is the rank-two/polar functional in
the standard normalization. -/
opaque fullWeilPoleSide : FullWeilTest → ℝ

/-- The continuation of the internal three-channel pairing to the full
Weil test class, now visibly decomposed by channel. -/
noncomputable def fullWeilForm (F : FullWeilTest) : ℝ :=
  fullWeilArchSide F - fullWeilCombSide F + fullWeilPoleSide F

/-- The standard Weil explicit-formula quadratic form, independently
of the TFPT finite-window construction. -/
opaque standardExplicitFormula : FullWeilTest → ℝ

/-- Fixed-support approximation data.  `grid n` never exceeds the
target support, converges pointwise there, and converges in `L¹` on
that fixed compact interval.  These are the inputs needed for the
arch integral; the finite prime and polar channels require less. -/
def FullWeilTest.FixedSupportGridApproximation
    (F : FullWeilTest) (grid : ℕ → GridElement) : Prop :=
  (∀ n, (grid n).supportBound ≤ F.supportRadius) ∧
  (∀ u : ℝ, |u| ≤ F.supportRadius →
    Filter.Tendsto (fun n => (grid n).toFun u)
      Filter.atTop (nhds (F.toFun u))) ∧
  Filter.Tendsto
    (fun n => intervalIntegral
      (fun u : ℝ => |(grid n).toFun u - F.toFun u|)
      (-F.supportRadius) F.supportRadius MeasureTheory.volume)
    Filter.atTop (nhds 0)

/-- Density brick: every admissible fixed-support autocorrelation has
a dyadic `GridElement` autocorrelation approximation with the same
support bound.  Classical route: compactly supported `L²` step
approximation of the autocorrelation witness, then dyadic PL
interpolation. -/
def FullWeilFixedSupportGridDensity : Prop :=
  ∀ F : FullWeilTest, F.admissible →
    ∃ grid : ℕ → GridElement, F.FixedSupportGridApproximation grid

/-- Channel-continuity brick along the fixed-support approximation.
The arch component is dominated/L¹ convergence (r475 supplies the
u-space kernel and selected `Δ²` model); the comb component is a
finite sum of point evaluations; the pole component is finite-rank. -/
def FullWeilChannelContinuity : Prop :=
  ∀ (F : FullWeilTest) (grid : ℕ → GridElement), F.admissible →
    F.FixedSupportGridApproximation grid →
      Filter.Tendsto (fun n => weilArchSide (grid n))
        Filter.atTop (nhds (fullWeilArchSide F)) ∧
      Filter.Tendsto (fun n => weilCombSide (grid n))
        Filter.atTop (nhds (fullWeilCombSide F)) ∧
      Filter.Tendsto (fun n => weilPoleSide (grid n))
        Filter.atTop (nhds (fullWeilPoleSide F))

/-- Form convergence assembled from the three channel limits. -/
theorem fullWeilForm_tendsto_of_channels
    {F : FullWeilTest} {grid : ℕ → GridElement}
    (hchannels :
      Filter.Tendsto (fun n => weilArchSide (grid n))
        Filter.atTop (nhds (fullWeilArchSide F)) ∧
      Filter.Tendsto (fun n => weilCombSide (grid n))
        Filter.atTop (nhds (fullWeilCombSide F)) ∧
      Filter.Tendsto (fun n => weilPoleSide (grid n))
        Filter.atTop (nhds (fullWeilPoleSide F))) :
    Filter.Tendsto (fun n => weilForm (grid n))
      Filter.atTop (nhds (fullWeilForm F)) := by
  unfold weilForm fullWeilForm
  exact hchannels.1.sub hchannels.2.1 |>.add hchannels.2.2

/-- Positivity is closed under the form limit. -/
theorem fullWeilForm_nonneg_of_tendsto
    {F : FullWeilTest} {grid : ℕ → GridElement}
    (hpos : ∀ f : GridElement, 0 ≤ weilForm f)
    (hlim : Filter.Tendsto (fun n => weilForm (grid n))
      Filter.atTop (nhds (fullWeilForm F))) :
    0 ≤ fullWeilForm F :=
  ge_of_tendsto hlim (Filter.Eventually.of_forall fun n => hpos (grid n))

/-- Missing bridge 1: density/continuity from the native dyadic
`GridElement` class to the complete Weil test class. -/
def GridDenseExtension : Prop :=
  (∀ f : GridElement, 0 ≤ weilForm f) →
    ∀ F : FullWeilTest, F.admissible → 0 ≤ fullWeilForm F

/-- The logical dense-extension bridge is proved once the two precise
fixed-support bricks are supplied. -/
theorem grid_dense_extension_of_fixedSupport
    (hdense : FullWeilFixedSupportGridDensity)
    (hcontinuous : FullWeilChannelContinuity) :
    GridDenseExtension := by
  intro hpos F hF
  obtain ⟨grid, hgrid⟩ := hdense F hF
  exact fullWeilForm_nonneg_of_tendsto hpos
    (fullWeilForm_tendsto_of_channels (hcontinuous F grid hF hgrid))

/-- OPEN CLASSICAL BRICK 1a (r487): dyadic step approximation of the
compactly supported autocorrelation witness. -/
theorem fullWeil_fixedSupport_grid_density :
    FullWeilFixedSupportGridDensity := by
  sorry

/-- OPEN CLASSICAL BRICK 1b (r487): continuity of the three concrete
channels along the fixed-support approximation. -/
theorem fullWeil_channel_continuity :
    FullWeilChannelContinuity := by
  sorry

/-- Bridge 1 now consumes exactly the two named bricks above; its
positivity-transfer algebra is sorry-free. -/
theorem grid_dense_extension : GridDenseExtension :=
  grid_dense_extension_of_fixedSupport
    fullWeil_fixedSupport_grid_density fullWeil_channel_continuity

/-- Missing bridge 2: identify the continued custom three-channel form
with the standard Weil explicit formula. -/
def StandardExplicitFormulaIdentification : Prop :=
  ∀ F : FullWeilTest, F.admissible →
    fullWeilForm F = standardExplicitFormula F

/-- OPEN CLASSICAL BRIDGE 2 (r463): this must prove the prime,
archimedean, pole, and spectral/zero dictionaries with matching
normalizations. -/
theorem standard_explicit_formula_identification :
    StandardExplicitFormulaIdentification := by
  sorry

/-- A nontrivial zero in exactly Mathlib's sense: zero of
`riemannZeta`, not a negative even integer, and not the pole at `1`. -/
def IsNontrivialRiemannZetaZero (s : ℂ) : Prop :=
  riemannZeta s = 0 ∧
    (¬∃ n : ℕ, s = -2 * (n + 1)) ∧
    s ≠ 1

/-- Guinand--Weil explicit formula landing site (design, r487).

Intended theorem: for every admissible `F`, the arithmetic
prime+arch+pole form equals the multiplicity-weighted spectral sum
over nontrivial zeros of `riemannZeta`.  Mathlib v4.29.1 has no
`ZetaZero` enumeration or multiplicity API, so the spectral sum and
its local finiteness must be developed before this Prop is proved.

References: A. Weil (1952), E. Bombieri, ``The Riemann Hypothesis'',
Clay Mathematics Institute (2000), explicit-formula/Weil-criterion
sections. -/
def GuinandWeilExplicitFormula : Prop :=
  ∀ F : FullWeilTest, F.admissible →
    fullWeilForm F = standardExplicitFormula F

/-- Exact criterion brick: every off-critical nontrivial zeta zero can
be separated by an admissible autocorrelation whose standard explicit
formula is negative.  This packages the classical Weil-criterion
construction without inventing a nonexistent Mathlib zero API. -/
def FullWeilSeparatesOffCriticalZeros : Prop :=
  ∀ s : ℂ, IsNontrivialRiemannZetaZero s → s.re ≠ 1 / 2 →
    ∃ F : FullWeilTest, F.admissible ∧ standardExplicitFormula F < 0

/-- Missing bridge 3: the standard Weil criterion, stated against
mathlib's actual zeta interface `riemannZeta` and its formal
`RiemannHypothesis` predicate. -/
def StandardWeilCriterionToMathlibRH : Prop :=
  (∀ F : FullWeilTest, F.admissible →
      0 ≤ standardExplicitFormula F) →
    RiemannHypothesis

/-- The final Mathlib interface is pure logic once off-critical
separation is available. -/
theorem standard_weil_criterion_to_mathlib_rh_of_separation
    (hseparate : FullWeilSeparatesOffCriticalZeros) :
    StandardWeilCriterionToMathlibRH := by
  intro hpos s hz htrivial hpole
  by_contra hcritical
  obtain ⟨F, hF, hneg⟩ :=
    hseparate s ⟨hz, htrivial, hpole⟩ hcritical
  exact (not_lt_of_ge (hpos F hF)) hneg

/-- OPEN CLASSICAL BRICK 3 (r487): Weil's off-critical separation
construction.  The logical conversion to Mathlib `RiemannHypothesis`
is proved above. -/
theorem fullWeil_separates_offCritical_zeros :
    FullWeilSeparatesOffCriticalZeros := by
  sorry

/-- Bridge 3 is now a proved wrapper around the single named
separation theorem. -/
theorem standard_weil_criterion_to_mathlib_rh :
    StandardWeilCriterionToMathlibRH :=
  standard_weil_criterion_to_mathlib_rh_of_separation
    fullWeil_separates_offCritical_zeros

end RH
