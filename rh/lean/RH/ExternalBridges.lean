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
import Mathlib.Analysis.SpecialFunctions.Trigonometric.DerivHyp
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
  autocorrelation : ∃ h : ℝ → ℝ,
    MeasureTheory.MemLp h 2 MeasureTheory.volume ∧
    HasCompactSupport h ∧
    (∃ K : NNReal, LipschitzWith K h) ∧
    (∃ a : ℝ, ∀ t : ℝ, t < a ∨ a + supportRadius < t → h t = 0) ∧
    ∀ u : ℝ, toFun u =
      ∫ t : ℝ, h t * h (t + u) ∂MeasureTheory.volume
  admissible : Prop

/-- The actual half-open step function represented by a
`GridElement`: value `x i` on the dyadic cell `[iD,(i+1)D)` and zero
off the finite cell range. -/
noncomputable def GridElement.toStepFun (f : GridElement) (t : ℝ) : ℝ :=
  if ht : 0 ≤ t ∧ ⌊t / f.D0⌋₊ < f.steps then
    f.x ⟨⌊t / f.D0⌋₊, ht.2⟩
  else 0

/-- Explicit left-sampled dyadic grid element for a witness supported
in `[a,a+R]`.  Values are arbitrary reals; no value-grid rounding is
part of `GridElement`. -/
noncomputable def dyadicSampleGrid
    (h : ℝ → ℝ) (a R : ℝ) (m : ℕ) : GridElement where
  steps := Nat.floor (R * (2 : ℝ) ^ m)
  meshExp := m
  x i := h (a + (i : ℕ) * meshWidth m)

/-- Flooring the number of cells makes the sampled grid support fit
inside the declared target support. -/
theorem dyadicSampleGrid_supportBound_le
    (h : ℝ → ℝ) (a R : ℝ) (m : ℕ) (hR : 0 ≤ R) :
    (dyadicSampleGrid h a R m).supportBound ≤ R := by
  unfold dyadicSampleGrid GridElement.supportBound GridElement.D0
  dsimp
  have hpow : 0 < (2 : ℝ) ^ m := by positivity
  have hfloor :
      ((Nat.floor (R * (2 : ℝ) ^ m) : ℕ) : ℝ) ≤
        R * (2 : ℝ) ^ m :=
    Nat.floor_le (mul_nonneg hR hpow.le)
  calc
    ((Nat.floor (R * (2 : ℝ) ^ m) : ℕ) : ℝ) *
        (1 / (2 : ℝ) ^ m) ≤
      (R * (2 : ℝ) ^ m) * (1 / (2 : ℝ) ^ m) :=
        mul_le_mul_of_nonneg_right hfloor (by positivity)
    _ = R := by field_simp

/-- Full-class regularized u-space arch integrand, with the same
normalization as r475 `weilArchUIntegrand`. -/
noncomputable def fullWeilArchUIntegrand
    (F : FullWeilTest) (u : ℝ) : ℝ :=
  if u = 0 then 0
  else weilArchUWeight u *
    (Real.exp (-(3 / 2 : ℝ) * u) * F.toFun 0 - F.toFun u)

/-- Full-class archimedean channel, concretely the r475 u-space
pairing at the declared support radius. -/
noncomputable def fullWeilArchSide (F : FullWeilTest) : ℝ :=
  if 0 < F.supportRadius then
    let b := F.supportRadius
    let g0 := F.toFun 0
    let Cb := -(Real.eulerMascheroniConstant + Real.log Real.pi)
      - Real.log (1 - Real.exp (-2 * b))
    Cb * g0 +
      intervalIntegral (fullWeilArchUIntegrand F)
        0 b MeasureTheory.volume
  else 0

/-- Common finite anchor determined by the target support. -/
noncomputable def FullWeilTest.fullAnchor (F : FullWeilTest) : ℕ :=
  max 1 (Nat.ceil (Real.exp F.supportRadius))

/-- Full-class prime channel.  The support anchor makes the
von-Mangoldt pairing visibly finite. -/
noncomputable def fullWeilCombSide (F : FullWeilTest) : ℝ :=
  ∑ n ∈ windowAtoms F.fullAnchor, combMass n * F.toFun (Real.log n)

/-- Standard polar weight.  Since `Π''(u) = -2 cosh(u/2)` for the
r376 potential `Π = polePotential` and `polePairingZ` has an outer
minus sign, the integral weight is `+2 cosh(u/2)`.  This is
equivalently the two evaluations of the bilateral transform at
`±1/2`. -/
noncomputable def fullWeilPoleWeight (u : ℝ) : ℝ :=
  2 * Real.cosh (u / 2)

/-- Full-class pole channel as the bounded-weight integral on the
fixed support. -/
noncomputable def fullWeilPoleSide (F : FullWeilTest) : ℝ :=
  intervalIntegral (fun u : ℝ => fullWeilPoleWeight u * F.toFun u)
    (-F.supportRadius) F.supportRadius MeasureTheory.volume

/-- The continuation of the internal three-channel pairing to the full
Weil test class, now visibly decomposed by channel. -/
noncomputable def fullWeilForm (F : FullWeilTest) : ℝ :=
  fullWeilArchSide F - fullWeilCombSide F + fullWeilPoleSide F

/-- The standard Weil explicit-formula quadratic form, independently
of the TFPT finite-window construction. -/
opaque standardExplicitFormula : FullWeilTest → ℝ

/-- Fixed-support approximation data.  `grid n` never exceeds the
target support, converges uniformly there, and converges in `L¹` on
that fixed compact interval.  Uniform convergence supplies every
point evaluation in the finite prime channel; `L¹` controls bounded
integral weights such as the standard polar `cosh` functional. -/
def FullWeilTest.FixedSupportGridApproximation
    (F : FullWeilTest) (grid : ℕ → GridElement) : Prop :=
  (∀ n, (grid n).supportBound ≤ F.supportRadius) ∧
  (∀ ε : ℝ, 0 < ε → ∀ᶠ n in Filter.atTop,
    ∀ u : ℝ, |u| ≤ F.supportRadius →
      |(grid n).toFun u - F.toFun u| < ε) ∧
  (∀ n, IntervalIntegrable (grid n).toFun MeasureTheory.volume
    (-F.supportRadius) F.supportRadius) ∧
  Filter.Tendsto
    (fun n => intervalIntegral
      (fun u : ℝ => |(grid n).toFun u - F.toFun u|)
      (-F.supportRadius) F.supportRadius MeasureTheory.volume)
    Filter.atTop (nhds 0)

/-- Uniform-on-support convergence plus the common support bound gives
pointwise convergence on all of `ℝ` (both functions vanish outside). -/
theorem FullWeilTest.FixedSupportGridApproximation.tendsto_toFun
    {F : FullWeilTest} {grid : ℕ → GridElement}
    (happrox : F.FixedSupportGridApproximation grid) (u : ℝ) :
    Filter.Tendsto (fun n => (grid n).toFun u)
      Filter.atTop (nhds (F.toFun u)) := by
  rw [Metric.tendsto_atTop]
  intro ε hε
  obtain ⟨N, hN⟩ :=
    Filter.eventually_atTop.1 (happrox.2.1 ε hε)
  refine ⟨N, fun n hn => ?_⟩
  rw [Real.dist_eq]
  by_cases hu : |u| ≤ F.supportRadius
  · exact hN n hn u hu
  · have hRu : F.supportRadius < |u| := lt_of_not_ge hu
    rw [F.support_toFun u hRu]
    rw [(grid n).toFun_eq_zero_of_lt_abs
      (lt_of_le_of_lt (happrox.1 n) hRu)]
    simpa using hε

/-- Every approximating grid support fits below the target's common
finite prime anchor. -/
theorem GridElement.elementAnchor_le_fullAnchor
    {F : FullWeilTest} {f : GridElement}
    (hsupport : f.supportBound ≤ F.supportRadius) :
    f.elementAnchor ≤ F.fullAnchor := by
  unfold GridElement.elementAnchor FullWeilTest.fullAnchor
  exact max_le_max_left 1
    (Nat.ceil_mono ((Real.exp_le_exp).2 hsupport))

/-- The native prime side equals the target-anchor finite sum whenever
the native support lies in the target support. -/
theorem weilCombSide_eq_fullAnchor
    {F : FullWeilTest} {f : GridElement}
    (hsupport : f.supportBound ≤ F.supportRadius) :
    weilCombSide f =
      ∑ n ∈ windowAtoms F.fullAnchor,
        combMass n * f.toFun (Real.log n) := by
  symm
  exact comb_elementwise_stabilization f
    (f.elementAnchor_le_fullAnchor hsupport)

/-- **Comb continuity (r489, PROVED).**  Fixed support turns the prime
channel into one common finite sum, and uniform convergence supplies
the finitely many point evaluations. -/
theorem fullWeilCombSide_tendsto
    {F : FullWeilTest} {grid : ℕ → GridElement}
    (happrox : F.FixedSupportGridApproximation grid) :
    Filter.Tendsto (fun n => weilCombSide (grid n))
      Filter.atTop (nhds (fullWeilCombSide F)) := by
  unfold fullWeilCombSide
  have heq : (fun k => weilCombSide (grid k)) =
      (fun k => ∑ n ∈ windowAtoms F.fullAnchor,
        combMass n * (grid k).toFun (Real.log n)) := by
    funext k
    exact weilCombSide_eq_fullAnchor (happrox.1 k)
  rw [heq]
  exact tendsto_finset_sum (windowAtoms F.fullAnchor) fun n _ =>
    (happrox.tendsto_toFun (Real.log n)).const_mul (combMass n)

/-- The standard polar integral is continuous for the chosen
fixed-support topology.  This is the analytic rank-two continuity
statement; identifying r376's native-mesh second-difference read with
this integral is recorded separately below. -/
theorem fullWeilPoleIntegral_tendsto
    {F : FullWeilTest} {grid : ℕ → GridElement}
    (happrox : F.FixedSupportGridApproximation grid) :
    Filter.Tendsto
      (fun n => intervalIntegral
        (fun u : ℝ => fullWeilPoleWeight u * (grid n).toFun u)
        (-F.supportRadius) F.supportRadius MeasureTheory.volume)
      Filter.atTop (nhds (fullWeilPoleSide F)) := by
  rw [tendsto_iff_norm_sub_tendsto_zero]
  unfold fullWeilPoleSide
  let C : ℝ := 2 * Real.cosh (F.supportRadius / 2)
  let err : ℕ → ℝ := fun n => intervalIntegral
    (fun u : ℝ => |(grid n).toFun u - F.toFun u|)
    (-F.supportRadius) F.supportRadius MeasureTheory.volume
  have hFint : IntervalIntegrable F.toFun MeasureTheory.volume
      (-F.supportRadius) F.supportRadius :=
    F.continuous_toFun.intervalIntegrable _ _
  have hweight : Continuous fullWeilPoleWeight := by
    unfold fullWeilPoleWeight
    fun_prop
  have hC : Filter.Tendsto (fun n => C * err n)
      Filter.atTop (nhds 0) := by
    simpa only [err, mul_zero] using happrox.2.2.2.const_mul C
  refine squeeze_zero' (Filter.Eventually.of_forall fun _ => norm_nonneg _) ?_ hC
  filter_upwards [] with n
  have hgridInt := happrox.2.2.1 n
  have hgridWeighted :=
    hgridInt.continuousOn_mul hweight.continuousOn
  have hFWeighted :=
    hFint.continuousOn_mul hweight.continuousOn
  rw [← intervalIntegral.integral_sub hgridWeighted hFWeighted]
  have hintegrand :
      (fun u : ℝ => fullWeilPoleWeight u * (grid n).toFun u -
        fullWeilPoleWeight u * F.toFun u) =
      (fun u : ℝ => fullWeilPoleWeight u *
        ((grid n).toFun u - F.toFun u)) := by
    funext u
    ring
  rw [hintegrand]
  change ‖intervalIntegral
    (fun u : ℝ => fullWeilPoleWeight u *
      ((grid n).toFun u - F.toFun u))
    (-F.supportRadius) F.supportRadius MeasureTheory.volume‖ ≤ C * err n
  calc
    _ ≤ intervalIntegral
        (fun u : ℝ => C * |(grid n).toFun u - F.toFun u|)
        (-F.supportRadius) F.supportRadius MeasureTheory.volume := by
      apply intervalIntegral.norm_integral_le_of_norm_le
        (neg_le_self F.supportRadius_nonneg)
      · filter_upwards [] with u hu
        have huabs : |u| ≤ F.supportRadius := by
          exact abs_le.mpr ⟨hu.1.le, hu.2⟩
        have hhalf : |u / 2| ≤ |F.supportRadius / 2| := by
          rw [abs_div, abs_div, abs_of_nonneg F.supportRadius_nonneg]
          exact div_le_div_of_nonneg_right huabs (by norm_num)
        have hcosh : Real.cosh (u / 2) ≤
            Real.cosh (F.supportRadius / 2) := by
          exact (Real.cosh_le_cosh).2 hhalf
        simp only [norm_mul, Real.norm_eq_abs, fullWeilPoleWeight]
        rw [abs_of_pos (by norm_num : (0 : ℝ) < 2),
          abs_of_pos (Real.cosh_pos _)]
        dsimp [C]
        exact mul_le_mul_of_nonneg_right
          (mul_le_mul_of_nonneg_left hcosh (by norm_num))
          (abs_nonneg ((grid n).toFun u - F.toFun u))
      · exact ((hgridInt.sub hFint).norm.const_mul C)
    _ = C * err n := by
      rw [intervalIntegral.integral_const_mul]

/-- Density brick: every admissible fixed-support autocorrelation has
a dyadic `GridElement` autocorrelation approximation with the same
support bound.  Classical route: compactly supported `L²` step
approximation of the autocorrelation witness, then dyadic PL
interpolation. -/
def FullWeilFixedSupportGridDensity : Prop :=
  ∀ F : FullWeilTest, F.admissible →
    ∃ grid : ℕ → GridElement, F.FixedSupportGridApproximation grid

/-- Exact dyadic-density target after r491: the approximating sequence
is no longer existentially anonymous, but the explicit left-sampled
`dyadicSampleGrid h a R`.  The remaining proof must establish its
uniform autocorrelation and `L¹` convergence from the support and
Lipschitz contracts. -/
def FullWeilDyadicSampleConvergence : Prop :=
  ∀ F : FullWeilTest, F.admissible →
    ∃ (h : ℝ → ℝ) (K : NNReal) (a : ℝ),
      MeasureTheory.MemLp h 2 MeasureTheory.volume ∧
      HasCompactSupport h ∧
      LipschitzWith K h ∧
      (∀ t : ℝ, t < a ∨ a + F.supportRadius < t → h t = 0) ∧
      (∀ u : ℝ, F.toFun u =
        ∫ t : ℝ, h t * h (t + u) ∂MeasureTheory.volume) ∧
      F.FixedSupportGridApproximation
        (fun m => dyadicSampleGrid h a F.supportRadius m)

/-- The explicit sampled-grid convergence target implies the
existential fixed-support density interface. -/
theorem fullWeilFixedSupportGridDensity_of_dyadicSample
    (hdyadic : FullWeilDyadicSampleConvergence) :
    FullWeilFixedSupportGridDensity := by
  intro F hF
  obtain ⟨h, K, a, hhLp, hhcompact, hhLip, hhsupport, hac, happrox⟩ :=
    hdyadic F hF
  exact ⟨fun m => dyadicSampleGrid h a F.supportRadius m, happrox⟩

/-- Existing r376 finite-element seam, now placed exactly: the
native-mesh second-difference pole read equals the compact-support
`+2 cosh(u/2)` integral.  Once this dictionary holds, pole continuity
is the proved `fullWeilPoleIntegral_tendsto`. -/
def GridPoleIntegralIdentification : Prop :=
  ∀ (F : FullWeilTest) (grid : ℕ → GridElement), F.admissible →
    F.FixedSupportGridApproximation grid →
      ∀ n, weilPoleSide (grid n) =
        intervalIntegral
          (fun u : ℝ => fullWeilPoleWeight u * (grid n).toFun u)
          (-F.supportRadius) F.supportRadius MeasureTheory.volume

/-- Exact finite-element identity needed for the r376 pole seam.
Each `GridElement.toFun` is the finite hat expansion of its ACF
samples; integrating a hat against `-polePotential'' = 2 cosh(u/2)`
gives the second difference in `polePairingZ`.  The larger interval
causes no change because `toFun` vanishes outside its support. -/
def GridPoleHatIntegralIdentity : Prop :=
  ∀ (f : GridElement) (R : ℝ), f.supportBound ≤ R →
    weilPoleSide f =
      intervalIntegral
        (fun u : ℝ => fullWeilPoleWeight u * f.toFun u)
        (-R) R MeasureTheory.volume

/-- r493c1: evenness split plus dyadic cell decomposition.  The two-sided
pole-density integral equals twice the sum of the nonnegative cell
integrals.  Each cell is the r493b affine integrand on
`[d D0, (d+1) D0]` via `toFun_eq_affine_on_nonneg_cell`. -/
theorem gridPoleIntegral_eq_two_mul_sum_cell
    (f : GridElement) {R : ℝ} (hR : f.supportBound ≤ R) :
    intervalIntegral
        (fun u : ℝ => fullWeilPoleWeight u * f.toFun u)
        (-R) R MeasureTheory.volume =
      2 * ∑ d ∈ Finset.range f.steps,
        intervalIntegral
          (fun u : ℝ => fullWeilPoleWeight u * f.toFun u)
          ((d : ℝ) * f.D0) (((d : ℝ) + 1) * f.D0)
          MeasureTheory.volume := by
  have hswap (a b : ℝ) :
      intervalIntegral
          (fun u : ℝ => fullWeilPoleWeight u * f.toFun u)
          a b MeasureTheory.volume =
        intervalIntegral
          (fun u : ℝ => f.toFun u * (2 * Real.cosh (u / 2)))
          a b MeasureTheory.volume := by
    refine intervalIntegral.integral_congr fun u _ => ?_
    unfold fullWeilPoleWeight
    ring
  simp_rw [hswap]
  exact f.intervalIntegral_toFun_mul_two_cosh_eq_two_mul_sum_cell hR

/-- r493c2: the native-mesh pole pairing equals the two-sided
`2 cosh(u/2)` integral.  Each cell evaluates to an
`affineCoshPrimitive` increment; the increments telescope against
`poleΔ` after the even one-sided rewrite of `polePairingZ`. -/
theorem gridPoleHatIntegralIdentity : GridPoleHatIntegralIdentity := by
  intro f R hR
  have hswap (a b : ℝ) :
      intervalIntegral
          (fun u : ℝ => fullWeilPoleWeight u * f.toFun u)
          a b MeasureTheory.volume =
        intervalIntegral
          (fun u : ℝ => f.toFun u * (2 * Real.cosh (u / 2)))
          a b MeasureTheory.volume := by
    refine intervalIntegral.integral_congr fun u _ => ?_
    unfold fullWeilPoleWeight
    ring
  rw [hswap, f.intervalIntegral_toFun_mul_two_cosh_eq_two_mul_sum_increment hR]
  exact weilPoleSide_eq_two_mul_sum_cellCoshIncrement f

/-- The pointwise hat identity supplies the sequence-level pole
dictionary used by channel continuity. -/
theorem gridPoleIntegralIdentification_of_hat
    (hhat : GridPoleHatIntegralIdentity) :
    GridPoleIntegralIdentification := by
  intro F grid _ happrox n
  exact hhat (grid n) F.supportRadius (happrox.1 n)

/-- Arch continuity component.  This is continuity of r475's concrete
u-space pairing along the fixed-support topology; Gauss--digamma is
needed later for the standard explicit-formula dictionary, not for
the topological positivity-limit argument itself. -/
def FullWeilArchContinuity : Prop :=
  ∀ (F : FullWeilTest) (grid : ℕ → GridElement), F.admissible →
    F.FixedSupportGridApproximation grid →
      Filter.Tendsto (fun n => weilArchSide (grid n))
        Filter.atTop (nhds (fullWeilArchSide F))

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

/-- The three channel limits follow from the remaining arch component
and the native-grid pole dictionary.  Comb continuity and standard
polar-integral continuity are proved above. -/
theorem fullWeilChannelContinuity_of_components
    (harch : FullWeilArchContinuity)
    (hpole : GridPoleIntegralIdentification) :
    FullWeilChannelContinuity := by
  intro F grid hF happrox
  refine ⟨harch F grid hF happrox, fullWeilCombSide_tendsto happrox, ?_⟩
  have heq : (fun n => weilPoleSide (grid n)) =
      (fun n => intervalIntegral
        (fun u : ℝ => fullWeilPoleWeight u * (grid n).toFun u)
        (-F.supportRadius) F.supportRadius MeasureTheory.volume) := by
    funext n
    exact hpole F grid hF happrox n
  rw [heq]
  exact fullWeilPoleIntegral_tendsto happrox

/-- **Single remaining dense-completion package (r489; r493c2 shrink).**

The r376 native-grid pole hat dictionary is now a proved theorem
(`gridPoleHatIntegralIdentity`).  Remaining components:
compactly-supported dyadic `L²` step density and autocorrelation
uniform convergence; r475 u-space / Dini arch continuity. -/
def FullWeilFixedSupportCompletion : Prop :=
  FullWeilDyadicSampleConvergence ∧
    FullWeilArchContinuity

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

/-- OPEN CLASSICAL BRICK 1 (r489): the one remaining completion
package.  No infinite-window or mesh-PSD hypothesis is included. -/
theorem fullWeil_fixedSupport_completion :
    FullWeilFixedSupportCompletion := by
  sorry

/-- Grid density is the first component of the single completion
package. -/
theorem fullWeil_fixedSupport_grid_density :
    FullWeilFixedSupportGridDensity :=
  fullWeilFixedSupportGridDensity_of_dyadicSample
    fullWeil_fixedSupport_completion.1

/-- All channel limits follow from the remaining arch completion
component, the proved pole hat identity, and the proved comb and
polar-integral continuity theorems. -/
theorem fullWeil_channel_continuity :
    FullWeilChannelContinuity :=
  fullWeilChannelContinuity_of_components
    fullWeil_fixedSupport_completion.2
    (gridPoleIntegralIdentification_of_hat gridPoleHatIntegralIdentity)

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
