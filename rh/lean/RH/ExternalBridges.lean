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
import Mathlib.Analysis.Analytic.Order
import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.Analysis.Normed.Module.Connected
import Mathlib.LinearAlgebra.Complex.FiniteDimensional
import Mathlib.Topology.Order.Basic
import Mathlib.MeasureTheory.Integral.Bochner.Set
import Mathlib.MeasureTheory.Function.LocallyIntegrable
import Mathlib.Analysis.Calculus.LogDeriv
import Mathlib.Analysis.Calculus.FDeriv.Analytic
import Mathlib.Analysis.Meromorphic.Order
import Mathlib.Topology.DiscreteSubset

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

/-- `D(m) = 2^{-m} → 0`.  Mesh of the r491 left-sampled grid. -/
lemma meshWidth_tendsto_zero :
    Filter.Tendsto meshWidth Filter.atTop (nhds 0) := by
  have h : Filter.Tendsto (fun m : ℕ => ((1 : ℝ) / 2) ^ m)
      Filter.atTop (nhds 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one (by positivity) (by norm_num)
  convert h using 1
  funext m
  unfold meshWidth
  exact (one_div_pow (2 : ℝ) m).symm

/-- Compactly supported continuous functions on `ℝ` are uniformly
continuous (Heine on a 1-neighbourhood of the support, both sides
vanish outside). -/
lemma FullWeilTest.uniformContinuous_toFun (F : FullWeilTest) :
    UniformContinuous F.toFun := by
  rw [Metric.uniformContinuous_iff]
  intro ε hε
  have hcpt :
      IsCompact (Set.Icc (-(F.supportRadius + 1)) (F.supportRadius + 1)) :=
    isCompact_Icc
  have hUC :=
    hcpt.uniformContinuousOn_of_continuous F.continuous_toFun.continuousOn
  rw [Metric.uniformContinuousOn_iff] at hUC
  obtain ⟨δ', hδ', hδu⟩ := hUC ε hε
  refine ⟨min δ' 1, lt_min hδ' (by norm_num), fun x y hdxy => ?_⟩
  have hd1 : dist x y < 1 := lt_of_lt_of_le hdxy (min_le_right _ _)
  have hdδ : dist x y < δ' := lt_of_lt_of_le hdxy (min_le_left _ _)
  have hmem {z : ℝ} (hz : |z| ≤ F.supportRadius + 1) :
      z ∈ Set.Icc (-(F.supportRadius + 1)) (F.supportRadius + 1) :=
    abs_le.mp hz
  by_cases hx : |x| ≤ F.supportRadius + 1
  · by_cases hy : |y| ≤ F.supportRadius + 1
    · exact hδu x (hmem hx) y (hmem hy) hdδ
    · have hyR : F.supportRadius < |y| :=
        lt_of_lt_of_le' (not_le.mp hy)
          (le_add_of_nonneg_right (by norm_num : (0 : ℝ) ≤ 1))
      have hxR : F.supportRadius < |x| := by
        have habs : |(|y| - |x|)| ≤ dist y x := abs_abs_sub_abs_le_abs_sub y x
        have : |(|y| - |x|)| < 1 :=
          lt_of_le_of_lt habs (by rwa [dist_comm, Real.dist_eq] at hd1)
        rw [abs_lt] at this
        linarith
      rw [F.support_toFun x hxR, F.support_toFun y hyR, dist_self]
      exact hε
  · have hxR : F.supportRadius < |x| :=
      lt_of_lt_of_le' (not_le.mp hx)
        (le_add_of_nonneg_right (by norm_num : (0 : ℝ) ≤ 1))
    have hyR : F.supportRadius < |y| := by
      have habs : |(|x| - |y|)| ≤ dist x y := abs_abs_sub_abs_le_abs_sub x y
      have : |(|x| - |y|)| < 1 :=
        lt_of_le_of_lt habs (by rwa [Real.dist_eq] at hd1)
      rw [abs_lt] at this
      linarith
    rw [F.support_toFun x hxR, F.support_toFun y hyR, dist_self]
    exact hε

lemma FullWeilTest.toFun_zero_at_radius (F : FullWeilTest) :
    F.toFun F.supportRadius = 0 := by
  have hs : Filter.Tendsto (fun n : ℕ => F.supportRadius + 1 / (n + 1 : ℝ))
      Filter.atTop (nhds F.supportRadius) := by
    simpa using
      (tendsto_const_nhds (x := F.supportRadius)).add
        tendsto_one_div_add_atTop_nhds_zero_nat
  have hval : ∀ n : ℕ, F.toFun (F.supportRadius + 1 / (n + 1 : ℝ)) = 0 :=
    fun n => by
      apply F.support_toFun
      rw [abs_of_nonneg
        (add_nonneg F.supportRadius_nonneg (by positivity))]
      exact lt_add_of_pos_right _ (by positivity)
  have hseq : Filter.Tendsto
      (fun n : ℕ => F.toFun (F.supportRadius + 1 / (n + 1 : ℝ)))
      Filter.atTop (nhds 0) := by
    have : (fun n : ℕ => F.toFun (F.supportRadius + 1 / (n + 1 : ℝ))) =
        fun _ => (0 : ℝ) := funext hval
    rw [this]
    exact tendsto_const_nhds
  exact tendsto_nhds_unique
    ((F.continuous_toFun.tendsto F.supportRadius).comp hs) hseq

lemma lipschitz_witness_zero_left {h : ℝ → ℝ} {K : NNReal} {a R : ℝ}
    (hK : LipschitzWith K h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + R < t → h t = 0) :
    h a = 0 := by
  have hs : Filter.Tendsto (fun n : ℕ => a - 1 / (n + 1 : ℝ))
      Filter.atTop (nhds a) := by
    simpa using
      (tendsto_const_nhds (x := a)).sub tendsto_one_div_add_atTop_nhds_zero_nat
  have hval : ∀ n : ℕ, h (a - 1 / (n + 1 : ℝ)) = 0 := fun n =>
    hsupp _ (Or.inl (sub_lt_self a (by positivity)))
  have hseq : Filter.Tendsto (fun n : ℕ => h (a - 1 / (n + 1 : ℝ)))
      Filter.atTop (nhds 0) := by
    have : (fun n : ℕ => h (a - 1 / (n + 1 : ℝ))) = fun _ => (0 : ℝ) :=
      funext hval
    rw [this]
    exact tendsto_const_nhds
  exact tendsto_nhds_unique ((hK.continuous.tendsto a).comp hs) hseq

lemma lipschitz_witness_zero_right {h : ℝ → ℝ} {K : NNReal} {a R : ℝ}
    (hK : LipschitzWith K h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + R < t → h t = 0) :
    h (a + R) = 0 := by
  have hs : Filter.Tendsto (fun n : ℕ => a + R + 1 / (n + 1 : ℝ))
      Filter.atTop (nhds (a + R)) := by
    simpa using
      (tendsto_const_nhds (x := a + R)).add tendsto_one_div_add_atTop_nhds_zero_nat
  have hval : ∀ n : ℕ, h (a + R + 1 / (n + 1 : ℝ)) = 0 := fun n =>
    hsupp _ (Or.inr (lt_add_of_pos_right _ (by positivity)))
  have hseq : Filter.Tendsto (fun n : ℕ => h (a + R + 1 / (n + 1 : ℝ)))
      Filter.atTop (nhds 0) := by
    have : (fun n : ℕ => h (a + R + 1 / (n + 1 : ℝ))) = fun _ => (0 : ℝ) :=
      funext hval
    rw [this]
    exact tendsto_const_nhds
  exact tendsto_nhds_unique ((hK.continuous.tendsto (a + R)).comp hs) hseq

lemma lipschitz_witness_abs_le {h : ℝ → ℝ} {K : NNReal} {a R : ℝ}
    (hK : LipschitzWith K h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + R < t → h t = 0)
    (hR : 0 ≤ R) (t : ℝ) :
    |h t| ≤ K * R := by
  rcases lt_or_ge t a with hlt | hge
  · rw [hsupp t (Or.inl hlt), abs_zero]
    exact mul_nonneg (NNReal.coe_nonneg _) hR
  rcases lt_or_ge (a + R) t with hgt | _hle
  · rw [hsupp t (Or.inr hgt), abs_zero]
    exact mul_nonneg (NNReal.coe_nonneg _) hR
  · have ha : h a = 0 := lipschitz_witness_zero_left hK hsupp
    have hdist := hK.dist_le_mul t a
    rw [Real.dist_eq, Real.dist_eq, ha, sub_zero] at hdist
    have ht : t - a ≤ R := by linarith
    have ht0 : 0 ≤ t - a := sub_nonneg.mpr hge
    rw [abs_of_nonneg ht0] at hdist
    exact le_trans hdist (mul_le_mul_of_nonneg_left ht (NNReal.coe_nonneg _))

lemma autocorrelation_eq_intervalIntegral
    {h : ℝ → ℝ} {a R u : ℝ}
    (hc : Continuous h) (hs : HasCompactSupport h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + R < t → h t = 0)
    (hR : 0 ≤ R) :
    ∫ t : ℝ, h t * h (t + u) ∂MeasureTheory.volume =
      intervalIntegral (fun t : ℝ => h t * h (t + u))
        a (a + R) MeasureTheory.volume := by
  set φ : ℝ → ℝ := fun t => h t * h (t + u)
  have _hint :
      MeasureTheory.Integrable φ MeasureTheory.volume :=
    (hc.mul (hc.comp (by fun_prop))).integrable_of_hasCompactSupport
      (by simpa [φ, Pi.mul_apply] using
        hs.mul_right (f' := fun t => h (t + u)))
  have h0 : ∀ t : ℝ, t ∉ Set.Icc a (a + R) → φ t = 0 := fun t ht => by
    have hout : t < a ∨ a + R < t := by
      rw [Set.mem_Icc, not_and_or] at ht
      rcases ht with hta | htR
      · exact Or.inl (lt_of_not_ge hta)
      · exact Or.inr (lt_of_not_ge htR)
    change h t * h (t + u) = 0
    rw [hsupp t hout, zero_mul]
  have hab : a ≤ a + R := le_add_of_nonneg_right hR
  rw [← MeasureTheory.setIntegral_eq_integral_of_forall_compl_eq_zero h0]
  rw [MeasureTheory.integral_Icc_eq_integral_Ioc]
  exact (intervalIntegral.integral_of_le hab).symm

lemma autocorrelation_integrand_abs_le
    {h : ℝ → ℝ} {K : NNReal} {a R : ℝ}
    (hK : LipschitzWith K h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + R < t → h t = 0)
    (hR : 0 ≤ R) (t u : ℝ) :
    |h t * h (t + u)| ≤ (K * R) ^ 2 := by
  have hb := lipschitz_witness_abs_le hK hsupp hR
  have hKR : 0 ≤ (K : ℝ) * R :=
    mul_nonneg (NNReal.coe_nonneg _) hR
  calc
    |h t * h (t + u)| = |h t| * |h (t + u)| := abs_mul _ _
    _ ≤ (K * R) * (K * R) :=
      mul_le_mul (hb t) (hb (t + u)) (abs_nonneg _) hKR
    _ = (K * R) ^ 2 := by ring

lemma autocorrelation_integrand_lip
    {h : ℝ → ℝ} {K : NNReal} {a R : ℝ}
    (hK : LipschitzWith K h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + R < t → h t = 0)
    (hR : 0 ≤ R) (u x y : ℝ) :
    |h x * h (x + u) - h y * h (y + u)| ≤
      (2 * (K : ℝ) * K * R) * |x - y| := by
  have hb := lipschitz_witness_abs_le hK hsupp hR
  have hdx := hK.dist_le_mul x y
  have hdu := hK.dist_le_mul (x + u) (y + u)
  rw [Real.dist_eq, Real.dist_eq] at hdx
  rw [Real.dist_eq, Real.dist_eq, add_sub_add_right_eq_sub] at hdu
  have hKR : 0 ≤ (K : ℝ) * R :=
    mul_nonneg (NNReal.coe_nonneg _) hR
  have hK0 : 0 ≤ (K : ℝ) := NNReal.coe_nonneg _
  calc
    |h x * h (x + u) - h y * h (y + u)| =
        |(h x - h y) * h (x + u) + h y * (h (x + u) - h (y + u))| := by
      congr 1
      ring
    _ ≤ |(h x - h y) * h (x + u)| + |h y * (h (x + u) - h (y + u))| :=
      abs_add_le _ _
    _ = |h x - h y| * |h (x + u)| + |h y| * |h (x + u) - h (y + u)| := by
      rw [abs_mul, abs_mul]
    _ ≤ (K * |x - y|) * (K * R) + (K * R) * (K * |x - y|) :=
      add_le_add
        (mul_le_mul hdx (hb (x + u)) (abs_nonneg _) (mul_nonneg hK0 (abs_nonneg _)))
        (mul_le_mul (hb y) hdu (abs_nonneg _) hKR)
    _ = (2 * (K : ℝ) * K * R) * |x - y| := by ring

lemma intervalIntegral_lipschitz_left_riemann
    {φ : ℝ → ℝ} {L D α : ℝ}
    (hL : ∀ x y : ℝ, |φ x - φ y| ≤ L * |x - y|)
    (hL0 : 0 ≤ L) (hD : 0 ≤ D)
    (hint : IntervalIntegrable φ MeasureTheory.volume α (α + D)) :
    |intervalIntegral φ α (α + D) MeasureTheory.volume - φ α * D| ≤
      L * D ^ 2 := by
  have hα : α ≤ α + D := le_add_of_nonneg_right hD
  have hconst : IntervalIntegrable (fun _ : ℝ => φ α)
      MeasureTheory.volume α (α + D) :=
    intervalIntegrable_const
  have hsub := intervalIntegral.integral_sub hint hconst
  have hφconst :
      intervalIntegral (fun _ : ℝ => φ α) α (α + D)
        MeasureTheory.volume = φ α * D := by
    rw [intervalIntegral.integral_const, smul_eq_mul]
    ring
  have hrew :
      intervalIntegral φ α (α + D) MeasureTheory.volume - φ α * D =
        intervalIntegral (fun t : ℝ => φ t - φ α) α (α + D)
          MeasureTheory.volume := by
    rw [← hφconst, hsub]
  rw [hrew]
  have habs :=
    intervalIntegral.abs_integral_le_integral_abs
      (f := fun t : ℝ => φ t - φ α) (μ := MeasureTheory.volume) hα
  refine le_trans habs ?_
  have hbound : ∀ t ∈ Set.uIcc α (α + D),
      |(φ t - φ α)| ≤ L * D := by
    intro t ht
    rw [Set.uIcc_of_le hα] at ht
    have htα : |t - α| ≤ D := by
      have h1 : α ≤ t := ht.1
      have h2 : t ≤ α + D := ht.2
      rw [abs_of_nonneg (sub_nonneg.mpr h1)]
      linarith
    exact le_trans (hL t α) (mul_le_mul_of_nonneg_left htα hL0)
  have hintabs : IntervalIntegrable (fun t : ℝ => |φ t - φ α|)
      MeasureTheory.volume α (α + D) :=
    (hint.sub hconst).abs
  have hmono :=
    intervalIntegral.integral_mono_on hα hintabs
      (intervalIntegrable_const (c := L * D))
      (fun t ht => by
        have : t ∈ Set.uIcc α (α + D) := by
          rw [Set.uIcc_of_le hα]
          exact ht
        exact hbound t this)
  refine le_trans hmono ?_
  rw [intervalIntegral.integral_const, smul_eq_mul]
  ring_nf
  nlinarith

lemma dyadicSampleGrid_D0_eq (h : ℝ → ℝ) (a R : ℝ) (m : ℕ) :
    (dyadicSampleGrid h a R m).D0 = meshWidth m := rfl

lemma dyadicSampleGrid_supportBound_tail_lt
    (h : ℝ → ℝ) (a R : ℝ) (m : ℕ) (hR : 0 ≤ R) :
    R - (dyadicSampleGrid h a R m).supportBound < meshWidth m := by
  unfold dyadicSampleGrid GridElement.supportBound GridElement.D0
  dsimp
  have hpow : 0 < (2 : ℝ) ^ m := by positivity
  have hlt : R * (2 : ℝ) ^ m <
      ((Nat.floor (R * (2 : ℝ) ^ m) : ℕ) : ℝ) + 1 :=
    Nat.lt_floor_add_one _
  have : R - ((Nat.floor (R * (2 : ℝ) ^ m) : ℕ) : ℝ) * (1 / (2 : ℝ) ^ m) <
      1 / (2 : ℝ) ^ m := by
    have hgoal : R < ((Nat.floor (R * (2 : ℝ) ^ m) : ℕ) : ℝ) *
        (1 / (2 : ℝ) ^ m) + 1 / (2 : ℝ) ^ m := by
      field_simp
      linarith
    linarith
  simpa [meshWidth] using this

lemma lipschitz_witness_abs_le_mesh
    {h : ℝ → ℝ} {K : NNReal} {a R : ℝ} {m i : ℕ}
    (hK : LipschitzWith K h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + R < t → h t = 0)
    (hR : 0 ≤ R)
    (hi : (dyadicSampleGrid h a R m).steps ≤ i) :
    |h (a + (i : ℝ) * meshWidth m)| ≤ K * meshWidth m := by
  set D := meshWidth m
  set n := (dyadicSampleGrid h a R m).steps
  have hD : 0 < D := meshWidth_pos m
  have hnR : (n : ℝ) * D ≤ R := by
    simpa [GridElement.supportBound, dyadicSampleGrid_D0_eq] using
      dyadicSampleGrid_supportBound_le h a R m hR
  have htail : R - (n : ℝ) * D ≤ D :=
    le_of_lt (by
      simpa [GridElement.supportBound, dyadicSampleGrid_D0_eq] using
        dyadicSampleGrid_supportBound_tail_lt h a R m hR)
  have haR : h (a + R) = 0 := lipschitz_witness_zero_right hK hsupp
  rcases le_or_gt ((i : ℝ) * D) R with hle | hgt
  · have hdist := hK.dist_le_mul (a + (i : ℝ) * D) (a + R)
    rw [Real.dist_eq, Real.dist_eq, haR, sub_zero, add_sub_add_left_eq_sub] at hdist
    have hid : 0 ≤ R - (i : ℝ) * D := sub_nonneg.mpr hle
    rw [abs_sub_comm, abs_of_nonneg hid] at hdist
    have hile : (n : ℝ) * D ≤ (i : ℝ) * D :=
      mul_le_mul_of_nonneg_right (Nat.cast_le.mpr hi) hD.le
    have hrest : R - (i : ℝ) * D ≤ R - (n : ℝ) * D := by linarith
    exact le_trans hdist
      (mul_le_mul_of_nonneg_left (le_trans hrest htail) (NNReal.coe_nonneg _))
  · rw [hsupp _ (Or.inr (by linarith)), abs_zero]
    exact mul_nonneg (NNReal.coe_nonneg _) hD.le

lemma dyadicSampleGrid_acf_eq
    (h : ℝ → ℝ) (a R : ℝ) (m d : ℕ) :
    (dyadicSampleGrid h a R m).acf d =
      meshWidth m *
        ∑ i : Fin (dyadicSampleGrid h a R m).steps,
          if hi : (i : ℕ) + d < (dyadicSampleGrid h a R m).steps then
            h (a + (i : ℕ) * meshWidth m) *
              h (a + (i : ℕ) * meshWidth m + (d : ℝ) * meshWidth m)
          else 0 := by
  unfold GridElement.acf
  rw [dyadicSampleGrid_D0_eq]
  congr 1
  refine Finset.sum_congr rfl fun i _ => ?_
  dsimp [dyadicSampleGrid]
  split_ifs with hi
  · have harg : a + ((↑i + d : ℕ) : ℝ) * meshWidth m =
        a + (i : ℕ) * meshWidth m + (d : ℝ) * meshWidth m := by
      have hcast : ((↑i + d : ℕ) : ℝ) = (i : ℕ) + (d : ℝ) :=
        Nat.cast_add (i : ℕ) d
      rw [hcast]
      ring
    rw [harg]
  · rfl

lemma dyadicSampleGrid_acf_riemann_gap
    {h : ℝ → ℝ} {K : NNReal} {a R : ℝ}
    (hK : LipschitzWith K h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + R < t → h t = 0)
    (hR : 0 ≤ R) (m d : ℕ) :
    let D := meshWidth m
    let n := (dyadicSampleGrid h a R m).steps
    let φ : ℝ → ℝ := fun t => h t * h (t + (d : ℝ) * D)
    |D * ∑ i : Fin n,
          (if hi : (i : ℕ) + d < n then φ (a + (i : ℕ) * D) else 0) -
        D * ∑ i : Fin n, φ (a + (i : ℕ) * D)| ≤
      ((K : ℝ) * R) ^ 2 * D := by
  intro D n φ
  have hD : 0 < D := meshWidth_pos m
  have hK0 : 0 ≤ (K : ℝ) := NNReal.coe_nonneg _
  have hKR : 0 ≤ (K : ℝ) * R := mul_nonneg hK0 hR
  have hnR : (n : ℝ) * D ≤ R := by
    simpa [GridElement.supportBound, dyadicSampleGrid_D0_eq] using
      dyadicSampleGrid_supportBound_le h a R m hR
  have hfac :
      |D * ∑ i : Fin n,
            (if hi : (i : ℕ) + d < n then φ (a + (i : ℕ) * D) else 0) -
          D * ∑ i : Fin n, φ (a + (i : ℕ) * D)| =
        D * |∑ i : Fin n,
            ((if hi : (i : ℕ) + d < n then φ (a + (i : ℕ) * D) else 0) -
              φ (a + (i : ℕ) * D))| := by
    rw [← mul_sub, abs_mul, abs_of_pos hD]
    congr 2
    exact (Finset.sum_sub_distrib (s := Finset.univ)
      (f := fun i : Fin n =>
        if hi : (i : ℕ) + d < n then φ (a + (i : ℕ) * D) else 0)
      (g := fun i : Fin n => φ (a + (i : ℕ) * D))).symm
  rw [hfac]
  have hterm (i : Fin n) :
      |(if hi : (i : ℕ) + d < n then φ (a + (i : ℕ) * D) else 0) -
          φ (a + (i : ℕ) * D)| ≤ (K * R) * (K * D) := by
    by_cases hi : (i : ℕ) + d < n
    · rw [dif_pos hi, sub_self, abs_zero]
      exact mul_nonneg hKR (mul_nonneg hK0 hD.le)
    · rw [dif_neg hi, zero_sub, abs_neg]
      have hge : n ≤ (i : ℕ) + d := Nat.le_of_not_gt hi
      have htailv :=
        lipschitz_witness_abs_le_mesh (m := m) (i := (i : ℕ) + d)
          hK hsupp hR (by simpa [n] using hge)
      have hb := lipschitz_witness_abs_le hK hsupp hR (a + (i : ℕ) * D)
      have hφ : |φ (a + (i : ℕ) * D)| =
          |h (a + (i : ℕ) * D)| *
            |h (a + (i : ℕ) * D + (d : ℝ) * D)| := by
        simp [φ, abs_mul]
      have harg : a + (i : ℕ) * D + (d : ℝ) * D =
          a + ((↑i + d : ℕ) : ℝ) * D := by
        have hcast : ((↑i + d : ℕ) : ℝ) = (i : ℕ) + (d : ℝ) :=
          Nat.cast_add (i : ℕ) d
        rw [hcast]
        ring
      rw [hφ, harg]
      exact mul_le_mul hb htailv (abs_nonneg _) hKR
  have hsumle :=
    Finset.abs_sum_le_sum_abs
      (s := (Finset.univ : Finset (Fin n)))
      (f := fun i : Fin n =>
        (if hi : (i : ℕ) + d < n then φ (a + (i : ℕ) * D) else 0) -
          φ (a + (i : ℕ) * D))
  refine le_trans (mul_le_mul_of_nonneg_left hsumle hD.le) ?_
  have hsumb :
      ∑ i : Fin n,
          |(if hi : (i : ℕ) + d < n then φ (a + (i : ℕ) * D) else 0) -
            φ (a + (i : ℕ) * D)| ≤
        (n : ℝ) * ((K * R) * (K * D)) := by
    refine le_trans (Finset.sum_le_sum fun i _ => hterm i) ?_
    simp [Finset.sum_const, nsmul_eq_mul]
  refine le_trans (mul_le_mul_of_nonneg_left hsumb hD.le) ?_
  calc
    D * ((n : ℝ) * ((K * R) * (K * D))) =
        ((n : ℝ) * D) * ((K * R) * (K * D)) := by ring
    _ ≤ R * ((K * R) * (K * D)) :=
      mul_le_mul_of_nonneg_right hnR (mul_nonneg hKR (mul_nonneg hK0 hD.le))
    _ = ((K : ℝ) * R) ^ 2 * D := by ring

lemma dyadicSampleGrid_riemann_cell_error
    {h : ℝ → ℝ} {K : NNReal} {a R : ℝ}
    (hK : LipschitzWith K h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + R < t → h t = 0)
    (hR : 0 ≤ R) (m d : ℕ) :
    let D := meshWidth m
    let n := (dyadicSampleGrid h a R m).steps
    let φ : ℝ → ℝ := fun t => h t * h (t + (d : ℝ) * D)
    |D * ∑ k ∈ Finset.range n, φ (a + (k : ℝ) * D) -
        intervalIntegral φ a (a + (n : ℝ) * D)
          MeasureTheory.volume| ≤
      2 * ((K : ℝ) * R) ^ 2 * D := by
  intro D n φ
  have hD : 0 < D := meshWidth_pos m
  have hK0 : 0 ≤ (K : ℝ) := NNReal.coe_nonneg _
  have hL0 : 0 ≤ 2 * (K : ℝ) * K * R := by positivity
  have hnR : (n : ℝ) * D ≤ R := by
    simpa [GridElement.supportBound, dyadicSampleGrid_D0_eq] using
      dyadicSampleGrid_supportBound_le h a R m hR
  have hφc : Continuous φ :=
    hK.continuous.mul (hK.continuous.comp (by fun_prop))
  have hφint (b c : ℝ) : IntervalIntegrable φ MeasureTheory.volume b c :=
    hφc.intervalIntegrable _ _
  have hφlip (x y : ℝ) :
      |φ x - φ y| ≤ (2 * (K : ℝ) * K * R) * |x - y| :=
    autocorrelation_integrand_lip hK hsupp hR ((d : ℝ) * D) x y
  have hadj :=
    intervalIntegral.sum_integral_adjacent_intervals
      (f := φ) (μ := MeasureTheory.volume) (n := n)
      (a := fun k : ℕ => a + (k : ℝ) * D)
      (fun k _ => hφint _ _)
  have hadj' :
      ∑ k ∈ Finset.range n,
          intervalIntegral φ (a + (k : ℝ) * D)
            (a + ((k + 1 : ℕ) : ℝ) * D) MeasureTheory.volume =
        intervalIntegral φ a (a + (n : ℝ) * D)
          MeasureTheory.volume := by
    convert hadj
    · simp
  have hcell (k : ℕ) (_hk : k ∈ Finset.range n) :
      |intervalIntegral φ (a + (k : ℝ) * D)
          (a + ((k + 1 : ℕ) : ℝ) * D) MeasureTheory.volume -
        φ (a + (k : ℝ) * D) * D| ≤
        (2 * (K : ℝ) * K * R) * D ^ 2 := by
    have hcell' : a + ((k + 1 : ℕ) : ℝ) * D = (a + (k : ℝ) * D) + D := by
      rw [Nat.cast_succ]
      ring
    rw [hcell']
    exact intervalIntegral_lipschitz_left_riemann hφlip hL0 hD.le (hφint _ _)
  have hsumcells :
      |∑ k ∈ Finset.range n,
            intervalIntegral φ (a + (k : ℝ) * D)
              (a + ((k + 1 : ℕ) : ℝ) * D) MeasureTheory.volume -
          ∑ k ∈ Finset.range n, φ (a + (k : ℝ) * D) * D| ≤
        ∑ k ∈ Finset.range n, (2 * (K : ℝ) * K * R) * D ^ 2 := by
    have hsub :
        ∑ k ∈ Finset.range n,
              intervalIntegral φ (a + (k : ℝ) * D)
                (a + ((k + 1 : ℕ) : ℝ) * D) MeasureTheory.volume -
            ∑ k ∈ Finset.range n, φ (a + (k : ℝ) * D) * D =
          ∑ k ∈ Finset.range n,
            (intervalIntegral φ (a + (k : ℝ) * D)
              (a + ((k + 1 : ℕ) : ℝ) * D) MeasureTheory.volume -
              φ (a + (k : ℝ) * D) * D) :=
      (Finset.sum_sub_distrib
        (s := Finset.range n)
        (fun k => intervalIntegral φ (a + (k : ℝ) * D)
          (a + ((k + 1 : ℕ) : ℝ) * D) MeasureTheory.volume)
        (fun k => φ (a + (k : ℝ) * D) * D)).symm
    rw [hsub]
    exact le_trans (Finset.abs_sum_le_sum_abs _ _)
      (Finset.sum_le_sum fun k hk => hcell k hk)
  have hrewsum :
      ∑ k ∈ Finset.range n, φ (a + (k : ℝ) * D) * D =
        D * ∑ k ∈ Finset.range n, φ (a + (k : ℝ) * D) := by
    simp [mul_comm, Finset.mul_sum]
  rw [hadj', hrewsum, abs_sub_comm] at hsumcells
  refine le_trans hsumcells ?_
  have hcard :
      ∑ _k ∈ Finset.range n, (2 * (K : ℝ) * K * R) * D ^ 2 =
        (n : ℝ) * ((2 * (K : ℝ) * K * R) * D ^ 2) := by
    simp [Finset.sum_const, nsmul_eq_mul]
  rw [hcard]
  calc
    (n : ℝ) * ((2 * (K : ℝ) * K * R) * D ^ 2) =
        ((n : ℝ) * D) * (2 * (K : ℝ) * K * R) * D := by ring
    _ ≤ R * (2 * (K : ℝ) * K * R) * D :=
      mul_le_mul_of_nonneg_right
        (mul_le_mul_of_nonneg_right hnR (by positivity)) hD.le
    _ = 2 * ((K : ℝ) * R) ^ 2 * D := by ring

lemma dyadicSampleGrid_tail_integral_le
    {h : ℝ → ℝ} {K : NNReal} {a R : ℝ}
    (hK : LipschitzWith K h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + R < t → h t = 0)
    (hR : 0 ≤ R) (m d : ℕ) :
    let D := meshWidth m
    let n := (dyadicSampleGrid h a R m).steps
    let φ : ℝ → ℝ := fun t => h t * h (t + (d : ℝ) * D)
    |intervalIntegral φ a (a + (n : ℝ) * D) MeasureTheory.volume -
        intervalIntegral φ a (a + R) MeasureTheory.volume| ≤
      ((K : ℝ) * R) ^ 2 * D := by
  intro D n φ
  have hD : 0 < D := meshWidth_pos m
  have hnR : (n : ℝ) * D ≤ R := by
    simpa [GridElement.supportBound, dyadicSampleGrid_D0_eq] using
      dyadicSampleGrid_supportBound_le h a R m hR
  have htail : R - (n : ℝ) * D < D := by
    simpa [GridElement.supportBound, dyadicSampleGrid_D0_eq] using
      dyadicSampleGrid_supportBound_tail_lt h a R m hR
  have hφc : Continuous φ :=
    hK.continuous.mul (hK.continuous.comp (by fun_prop))
  have hφint (b c : ℝ) : IntervalIntegrable φ MeasureTheory.volume b c :=
    hφc.intervalIntegrable _ _
  have hsplit :=
    intervalIntegral.integral_add_adjacent_intervals
      (hφint a (a + (n : ℝ) * D)) (hφint (a + (n : ℝ) * D) (a + R))
  have hrew :
      intervalIntegral φ a (a + R) MeasureTheory.volume -
          intervalIntegral φ a (a + (n : ℝ) * D) MeasureTheory.volume =
        intervalIntegral φ (a + (n : ℝ) * D) (a + R)
          MeasureTheory.volume := by
    linarith [hsplit]
  rw [abs_sub_comm, hrew]
  have hle : a + (n : ℝ) * D ≤ a + R := by linarith
  have habs :=
    intervalIntegral.abs_integral_le_integral_abs
      (f := φ) (μ := MeasureTheory.volume) hle
  refine le_trans habs ?_
  have hmono :=
    intervalIntegral.integral_mono_on hle
      (hφint (a + (n : ℝ) * D) (a + R)).abs
      (intervalIntegrable_const (c := (K * R) ^ 2))
      (fun t _ht =>
        autocorrelation_integrand_abs_le hK hsupp hR t ((d : ℝ) * D))
  refine le_trans hmono ?_
  rw [intervalIntegral.integral_const, smul_eq_mul]
  have hlen : (a + R) - (a + (n : ℝ) * D) = R - (n : ℝ) * D := by ring
  rw [hlen, mul_comm]
  exact mul_le_mul_of_nonneg_left (le_of_lt htail) (sq_nonneg _)

lemma dyadicSampleGrid_sum_sample_eq
    (h : ℝ → ℝ) (a : ℝ) (m d n : ℕ) :
    let D := meshWidth m
    let φ : ℝ → ℝ := fun t => h t * h (t + (d : ℝ) * D)
    ∑ i : Fin n, φ (a + (i : ℕ) * D) =
      ∑ k ∈ Finset.range n, φ (a + (k : ℝ) * D) := by
  intro D φ
  rw [Finset.sum_fin_eq_sum_range]
  refine Finset.sum_congr rfl fun k hk => ?_
  rw [dif_pos (Finset.mem_range.mp hk)]

lemma dyadicSampleGrid_acf_sub_toFun_le
    {F : FullWeilTest} {h : ℝ → ℝ} {K : NNReal} {a : ℝ}
    (hK : LipschitzWith K h) (hhKs : HasCompactSupport h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + F.supportRadius < t → h t = 0)
    (hac : ∀ u : ℝ, F.toFun u =
      ∫ t : ℝ, h t * h (t + u) ∂MeasureTheory.volume)
    (m d : ℕ) :
    |(dyadicSampleGrid h a F.supportRadius m).acf d -
        F.toFun ((d : ℝ) * meshWidth m)| ≤
      4 * ((K : ℝ) * F.supportRadius) ^ 2 * meshWidth m := by
  set R := F.supportRadius
  set D := meshWidth m
  set n := (dyadicSampleGrid h a R m).steps
  set φ : ℝ → ℝ := fun t => h t * h (t + (d : ℝ) * D)
  have hR : 0 ≤ R := F.supportRadius_nonneg
  have hD : 0 < D := meshWidth_pos m
  rw [dyadicSampleGrid_acf_eq, hac,
    autocorrelation_eq_intervalIntegral hK.continuous hhKs hsupp hR]
  have hgap := dyadicSampleGrid_acf_riemann_gap hK hsupp hR m d
  have hsum := dyadicSampleGrid_sum_sample_eq h a m d n
  have hcells := dyadicSampleGrid_riemann_cell_error hK hsupp hR m d
  have htailI := dyadicSampleGrid_tail_integral_le hK hsupp hR m d
  have hmid :
      |D * ∑ i : Fin n, φ (a + (i : ℕ) * D) -
          intervalIntegral φ a (a + (n : ℝ) * D)
            MeasureTheory.volume| ≤
        2 * ((K : ℝ) * R) ^ 2 * D := by
    simp only [D, φ] at hsum hcells ⊢
    rw [hsum]
    exact hcells
  refine le_trans
    (abs_sub_le (D * ∑ i : Fin n,
        (if hi : (i : ℕ) + d < n then φ (a + (i : ℕ) * D) else 0))
      (D * ∑ i : Fin n, φ (a + (i : ℕ) * D))
      (intervalIntegral φ a (a + R) MeasureTheory.volume)) ?_
  refine le_trans (add_le_add hgap
      (le_trans
        (abs_sub_le (D * ∑ i : Fin n, φ (a + (i : ℕ) * D))
          (intervalIntegral φ a (a + (n : ℝ) * D) MeasureTheory.volume)
          (intervalIntegral φ a (a + R) MeasureTheory.volume))
        (add_le_add hmid htailI))) ?_
  ring_nf
  nlinarith [sq_nonneg ((K : ℝ) * R), hD.le]

lemma dyadicSampleGrid_toFun_uniform
    {F : FullWeilTest} {h : ℝ → ℝ} {K : NNReal} {a : ℝ}
    (hK : LipschitzWith K h) (hhKs : HasCompactSupport h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + F.supportRadius < t → h t = 0)
    (hac : ∀ u : ℝ, F.toFun u =
      ∫ t : ℝ, h t * h (t + u) ∂MeasureTheory.volume)
    {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ m in Filter.atTop, ∀ u : ℝ, |u| ≤ F.supportRadius →
      |(dyadicSampleGrid h a F.supportRadius m).toFun u - F.toFun u| < ε := by
  have hUC := F.uniformContinuous_toFun
  rw [Metric.uniformContinuous_iff] at hUC
  obtain ⟨δ, hδ, hδu⟩ := hUC (ε / 6) (by positivity)
  have hDδ : ∀ᶠ m in Filter.atTop, meshWidth m < δ := by
    have h := meshWidth_tendsto_zero
    rw [Metric.tendsto_atTop] at h
    obtain ⟨N, hN⟩ := h δ hδ
    refine Filter.eventually_atTop.2 ⟨N, fun m hm => ?_⟩
    have := hN m hm
    rwa [Real.dist_eq, sub_zero, abs_of_nonneg (meshWidth_pos m).le] at this
  have hknot : ∀ᶠ m in Filter.atTop,
      4 * ((K : ℝ) * F.supportRadius) ^ 2 * meshWidth m < ε / 6 := by
    have htend :
        Filter.Tendsto
          (fun m : ℕ => 4 * ((K : ℝ) * F.supportRadius) ^ 2 * meshWidth m)
          Filter.atTop (nhds 0) := by
      simpa using meshWidth_tendsto_zero.const_mul
        (4 * ((K : ℝ) * F.supportRadius) ^ 2)
    rw [Metric.tendsto_atTop] at htend
    obtain ⟨N, hN⟩ := htend (ε / 6) (by positivity)
    refine Filter.eventually_atTop.2 ⟨N, fun m hm => ?_⟩
    have := hN m hm
    have hnn : 0 ≤ 4 * ((K : ℝ) * F.supportRadius) ^ 2 * meshWidth m :=
      mul_nonneg (mul_nonneg (by norm_num) (sq_nonneg _))
        (meshWidth_pos m).le
    rwa [Real.dist_eq, sub_zero, abs_of_nonneg hnn] at this
  filter_upwards [hDδ, hknot] with m hDm hkm
  intro u hu
  set R := F.supportRadius
  set D := meshWidth m
  set f := dyadicSampleGrid h a R m
  have hD : 0 < D := meshWidth_pos m
  have hR : 0 ≤ R := F.supportRadius_nonneg
  have hFeven : F.toFun u = F.toFun |u| := by
    rcases le_total 0 u with hu0 | hu0
    · rw [abs_of_nonneg hu0]
    · rw [abs_of_nonpos hu0, F.even_toFun]
  have hω {x y : ℝ} (hxy : dist x y < δ) :
      |F.toFun x - F.toFun y| < ε / 6 := by
    have := hδu hxy
    rwa [Real.dist_eq] at this
  have hknotu (k : ℕ) :
      |f.acf k - F.toFun ((k : ℝ) * D)| < ε / 6 :=
    lt_of_le_of_lt (dyadicSampleGrid_acf_sub_toFun_le hK hhKs hsupp hac m k) hkm
  set q : ℕ := ⌊|u| / D⌋₊
  have hqle : (q : ℝ) ≤ |u| / D := Nat.floor_le (div_nonneg (abs_nonneg _) hD.le)
  have hqt : |u| / D < (q : ℝ) + 1 := Nat.lt_floor_add_one _
  have hqD : (q : ℝ) * D ≤ |u| := (le_div_iff₀ hD).mp hqle
  have hqD' : |u| < ((q + 1 : ℕ) : ℝ) * D := by
    have : |u| / D < ((q + 1 : ℕ) : ℝ) := by
      rw [Nat.cast_succ]
      exact hqt
    exact (div_lt_iff₀ hD).mp this
  have hqd : dist ((q : ℝ) * D) |u| < δ := by
    rw [Real.dist_eq, abs_sub_comm, abs_of_nonneg (sub_nonneg.mpr hqD)]
    have : |u| - (q : ℝ) * D < D := by
      have hsucc : ((q + 1 : ℕ) : ℝ) = (q : ℝ) + 1 := Nat.cast_succ q
      have hqD'' := hqD'
      rw [hsucc] at hqD''
      linarith
    exact lt_trans this hDm
  have hFe : |F.toFun ((q : ℝ) * D) - F.toFun u| < ε / 6 := by
    rw [hFeven]
    exact hω hqd
  by_cases hq : q < f.steps
  · have hθ : 0 ≤ |u| / D - q := sub_nonneg.mpr hqle
    have hθ1 : |u| / D - q < 1 := by linarith
    have hD0 : f.D0 = D := dyadicSampleGrid_D0_eq h a R m
    have hto : f.toFun u =
        f.acf q + (|u| / D - q) * (f.acf (q + 1) - f.acf q) := by
      unfold GridElement.toFun
      rw [hD0]
    rw [hto]
    have hconv :
        f.acf q + (|u| / D - q) * (f.acf (q + 1) - f.acf q) =
          (1 - (|u| / D - q)) * f.acf q + (|u| / D - q) * f.acf (q + 1) := by
      ring
    rw [hconv]
    have hqd1 : dist (((q + 1 : ℕ) : ℝ) * D) |u| < δ := by
      rw [Real.dist_eq, abs_of_nonneg (sub_nonneg.mpr hqD'.le)]
      have : ((q + 1 : ℕ) : ℝ) * D - |u| ≤ D := by
        have : |u| ≥ (q : ℝ) * D := hqD
        have hsucc : ((q + 1 : ℕ) : ℝ) = (q : ℝ) + 1 := Nat.cast_succ q
        rw [hsucc]
        linarith
      exact lt_of_le_of_lt this hDm
    have hFe1 : |F.toFun (((q + 1 : ℕ) : ℝ) * D) - F.toFun u| < ε / 6 := by
      rw [hFeven]
      exact hω hqd1
    set θ := |u| / D - (q : ℝ)
    have hθ01 : 0 ≤ θ := hθ
    have hθle : θ ≤ 1 := le_of_lt hθ1
    have habs :
        |(1 - θ) * f.acf q + θ * f.acf (q + 1) - F.toFun u| ≤
          (1 - θ) * |f.acf q - F.toFun u| + θ * |f.acf (q + 1) - F.toFun u| := by
      have : (1 - θ) * f.acf q + θ * f.acf (q + 1) - F.toFun u =
          (1 - θ) * (f.acf q - F.toFun u) + θ * (f.acf (q + 1) - F.toFun u) := by
        ring
      rw [this]
      refine le_trans (abs_add_le _ _) ?_
      rw [abs_mul, abs_mul, abs_of_nonneg (sub_nonneg.mpr hθle), abs_of_nonneg hθ01]
    refine lt_of_le_of_lt habs ?_
    have h0 : |f.acf q - F.toFun u| < ε / 3 :=
      lt_of_le_of_lt
        (abs_sub_le (f.acf q) (F.toFun ((q : ℝ) * D)) (F.toFun u))
        (by linarith [hknotu q, hFe])
    have h1 : |f.acf (q + 1) - F.toFun u| < ε / 3 :=
      lt_of_le_of_lt
        (abs_sub_le (f.acf (q + 1))
          (F.toFun (((q + 1 : ℕ) : ℝ) * D)) (F.toFun u))
        (by linarith [hknotu (q + 1), hFe1])
    have hcomb :
        (1 - θ) * |f.acf q - F.toFun u| +
            θ * |f.acf (q + 1) - F.toFun u| ≤ ε / 3 := by
      have := add_le_add
        (mul_le_mul_of_nonneg_left (le_of_lt h0) (sub_nonneg.mpr hθle))
        (mul_le_mul_of_nonneg_left (le_of_lt h1) hθ01)
      refine le_trans this ?_
      ring_nf
      linarith
    exact lt_of_le_of_lt hcomb (by linarith)
  · have hsup : f.supportBound ≤ |u| := by
      have hqge : f.steps ≤ q := Nat.le_of_not_gt hq
      have : (f.steps : ℝ) * D ≤ (q : ℝ) * D :=
        mul_le_mul_of_nonneg_right (Nat.cast_le.mpr hqge) hD.le
      unfold GridElement.supportBound
      rw [dyadicSampleGrid_D0_eq]
      exact le_trans this hqD
    rw [f.toFun_eq_zero_of_supportBound_le hsup, zero_sub, abs_neg]
    have hFR : F.toFun R = 0 := F.toFun_zero_at_radius
    have hRu : dist |u| R < δ := by
      rw [Real.dist_eq, abs_sub_comm, abs_of_nonneg (sub_nonneg.mpr hu)]
      have htail : R - f.supportBound < D :=
        dyadicSampleGrid_supportBound_tail_lt h a R m hR
      have : R - |u| ≤ R - f.supportBound := by linarith
      exact lt_of_le_of_lt this (lt_trans htail hDm)
    have : |F.toFun u| < ε / 6 := by
      have hωu : |F.toFun |u| - F.toFun R| < ε / 6 := hω hRu
      rw [hFR, sub_zero] at hωu
      rwa [hFeven]
    linarith

lemma dyadicSampleGrid_l1_tendsto
    {F : FullWeilTest} {h : ℝ → ℝ} {K : NNReal} {a : ℝ}
    (hK : LipschitzWith K h) (hhKs : HasCompactSupport h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + F.supportRadius < t → h t = 0)
    (hac : ∀ u : ℝ, F.toFun u =
      ∫ t : ℝ, h t * h (t + u) ∂MeasureTheory.volume) :
    Filter.Tendsto
      (fun m => intervalIntegral
        (fun u : ℝ =>
          |(dyadicSampleGrid h a F.supportRadius m).toFun u - F.toFun u|)
        (-F.supportRadius) F.supportRadius MeasureTheory.volume)
      Filter.atTop (nhds 0) := by
  rw [Metric.tendsto_atTop]
  intro ε hε
  have hR : 0 ≤ F.supportRadius := F.supportRadius_nonneg
  have hε' : 0 < ε / (2 * F.supportRadius + 1) := by positivity
  have huni :=
    dyadicSampleGrid_toFun_uniform hK hhKs hsupp hac hε'
  obtain ⟨N, hN⟩ := Filter.eventually_atTop.1 huni
  refine ⟨N, fun m hm => ?_⟩
  have hpt := hN m hm
  have hgrid :
      IntervalIntegrable
        (dyadicSampleGrid h a F.supportRadius m).toFun
        MeasureTheory.volume (-F.supportRadius) F.supportRadius :=
    (dyadicSampleGrid h a F.supportRadius m).intervalIntegrable_toFun _ _
  have hFint : IntervalIntegrable F.toFun MeasureTheory.volume
      (-F.supportRadius) F.supportRadius :=
    F.continuous_toFun.intervalIntegrable _ _
  have habs :
      IntervalIntegrable
        (fun u =>
          |(dyadicSampleGrid h a F.supportRadius m).toFun u - F.toFun u|)
        MeasureTheory.volume (-F.supportRadius) F.supportRadius :=
    (hgrid.sub hFint).abs
  have hle : -F.supportRadius ≤ F.supportRadius := neg_le_self hR
  have hmono :=
    intervalIntegral.integral_mono_on hle habs
      (intervalIntegrable_const (c := ε / (2 * F.supportRadius + 1)))
      (fun u hu => by
        have huabs : |u| ≤ F.supportRadius := abs_le.mpr ⟨hu.1, hu.2⟩
        exact le_of_lt (hpt u huabs))
  have hbound :
      intervalIntegral
          (fun _ : ℝ => ε / (2 * F.supportRadius + 1))
          (-F.supportRadius) F.supportRadius MeasureTheory.volume =
        (2 * F.supportRadius) * (ε / (2 * F.supportRadius + 1)) := by
    rw [intervalIntegral.integral_const, smul_eq_mul]
    ring
  rw [Real.dist_eq, sub_zero]
  have hnn : 0 ≤ intervalIntegral
      (fun u : ℝ =>
        |(dyadicSampleGrid h a F.supportRadius m).toFun u - F.toFun u|)
      (-F.supportRadius) F.supportRadius MeasureTheory.volume :=
    intervalIntegral.integral_nonneg hle fun u hu => abs_nonneg _
  rw [abs_of_nonneg hnn]
  refine lt_of_le_of_lt (le_trans hmono (le_of_eq hbound)) ?_
  have hden : 0 < 2 * F.supportRadius + 1 := by positivity
  have : (2 * F.supportRadius) * (ε / (2 * F.supportRadius + 1)) < ε := by
    rw [mul_div_assoc']
    exact (div_lt_iff₀ hden).2 (by linarith)
  exact this

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
integral weights such as the standard polar `cosh` functional.
The last conjunct is the common Dini / Lipschitz majorant required
by the `1/u` arch weight: every approximating PL interpolant shares
one Lipschitz constant (dyadic ACF knots give this uniformly). -/
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
    Filter.atTop (nhds 0) ∧
  ∃ Kgrid : NNReal, ∀ n, LipschitzWith Kgrid (grid n).toFun

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
    simpa only [err, mul_zero] using happrox.2.2.2.1.const_mul C
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

/-- Cell slopes vanish past the last ACF knot. -/
lemma GridElement.cellSlope_eq_zero_of_steps_le (f : GridElement) {d : ℕ}
    (hd : f.steps ≤ d) : f.cellSlope d = 0 := by
  unfold GridElement.cellSlope
  rw [f.acf_eq_zero hd, f.acf_eq_zero (le_trans hd (Nat.le_succ _)),
    sub_zero, zero_div]

/-- Even PL interpolant depends only on `|u|`. -/
lemma GridElement.toFun_eq_toFun_abs (f : GridElement) (u : ℝ) :
    f.toFun u = f.toFun |u| := by
  rcases le_or_gt 0 u with hu | hu
  · rw [abs_of_nonneg hu]
  · rw [abs_of_neg hu, f.toFun_even]

/-- Same-cell increment of the nonnegative PL interpolant. -/
lemma GridElement.toFun_sub_le_of_mem_cell (f : GridElement) (d : ℕ)
    {L x y : ℝ} (hL : ∀ k : ℕ, |f.cellSlope k| ≤ L)
    (hx : (d : ℝ) * f.D0 ≤ x) (hxy : x ≤ y)
    (hy : y ≤ ((d : ℝ) + 1) * f.D0) :
    |f.toFun y - f.toFun x| ≤ L * (y - x) := by
  have hylo : (d : ℝ) * f.D0 ≤ y := le_trans hx hxy
  have hx0 : 0 ≤ x :=
    le_trans (mul_nonneg (Nat.cast_nonneg _) f.D0_pos.le) hx
  have hy0 : 0 ≤ y := le_trans hx0 hxy
  rw [f.toFun_eq_affine_on_nonneg_cell d hx (le_trans hxy hy),
    f.toFun_eq_affine_on_nonneg_cell d hylo hy]
  have hdiff :
      f.cellIntercept d + f.cellSlope d * y -
        (f.cellIntercept d + f.cellSlope d * x) =
      f.cellSlope d * (y - x) := by ring
  rw [hdiff, abs_mul, abs_of_nonneg (sub_nonneg.mpr hxy)]
  exact mul_le_mul_of_nonneg_right (hL d) (sub_nonneg.mpr hxy)

/-- Telescoping ACF increments from knot `i` across `n` cells. -/
lemma GridElement.acf_sub_eq_sum_cell (f : GridElement) (i n : ℕ) :
    f.acf (i + n) - f.acf i =
      ∑ k ∈ Finset.range n, (f.acf (i + k + 1) - f.acf (i + k)) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_range_succ, ← ih]
    ring

/-- Nonnegative restriction of the PL interpolant is Lipschitz. -/
lemma GridElement.toFun_nonneg_lip (f : GridElement) {L s t : ℝ}
    (hL : ∀ d : ℕ, |f.cellSlope d| ≤ L) (hs : 0 ≤ s) (hst : s ≤ t) :
    |f.toFun t - f.toFun s| ≤ L * (t - s) := by
  have hs0 : 0 ≤ s / f.D0 := div_nonneg hs f.D0_pos.le
  have htD : 0 ≤ t / f.D0 := div_nonneg (le_trans hs hst) f.D0_pos.le
  let i : ℕ := ⌊s / f.D0⌋₊
  let j : ℕ := ⌊t / f.D0⌋₊
  have hi : i = ⌊s / f.D0⌋₊ := rfl
  have hj : j = ⌊t / f.D0⌋₊ := rfl
  have hilo : (i : ℝ) * f.D0 ≤ s :=
    (le_div_iff₀ f.D0_pos).mp (Nat.floor_le hs0)
  have hjlo : (j : ℝ) * f.D0 ≤ t :=
    (le_div_iff₀ f.D0_pos).mp (Nat.floor_le htD)
  have hihi : s ≤ ((i : ℝ) + 1) * f.D0 := by
    have := Nat.lt_succ_floor (s / f.D0)
    rw [← hi, Nat.cast_succ] at this
    exact (div_le_iff₀ f.D0_pos).mp (le_of_lt this)
  have hjhi : t ≤ ((j : ℝ) + 1) * f.D0 := by
    have := Nat.lt_succ_floor (t / f.D0)
    rw [← hj, Nat.cast_succ] at this
    exact (div_le_iff₀ f.D0_pos).mp (le_of_lt this)
  have hij : i ≤ j := Nat.floor_le_floor
    (div_le_div_of_nonneg_right hst f.D0_pos.le)
  rcases eq_or_lt_of_le hij with hijeq | hijlt
  · rw [hijeq] at hilo hihi
    simpa [abs_sub_comm, hijeq] using
      f.toFun_sub_le_of_mem_cell (x := s) (y := t) j hL hilo hst hjhi
  · have hle : i + 1 ≤ j := Nat.succ_le_of_lt hijlt
    have h1 : |f.toFun t - f.toFun ((j : ℝ) * f.D0)| ≤
        L * (t - (j : ℝ) * f.D0) :=
      f.toFun_sub_le_of_mem_cell (x := (j : ℝ) * f.D0) (y := t) j hL
        le_rfl hjlo hjhi
    have h3 : |f.toFun (((i : ℝ) + 1) * f.D0) - f.toFun s| ≤
        L * (((i : ℝ) + 1) * f.D0 - s) := by
      simpa [abs_sub_comm] using
        f.toFun_sub_le_of_mem_cell (x := s) (y := ((i : ℝ) + 1) * f.D0)
          i hL hilo hihi le_rfl
    have h2 : |f.toFun ((j : ℝ) * f.D0) -
        f.toFun (((i : ℝ) + 1) * f.D0)| ≤
          L * ((j : ℝ) * f.D0 - ((i : ℝ) + 1) * f.D0) := by
      have hjk : f.toFun ((j : ℝ) * f.D0) = f.acf j :=
        f.toFun_nat_mul_D0 j
      have hik : f.toFun (((i : ℝ) + 1) * f.D0) = f.acf (i + 1) := by
        simpa [Nat.cast_succ] using f.toFun_nat_mul_D0 (i + 1)
      rw [hjk, hik]
      have hsum := f.acf_sub_eq_sum_cell (i + 1) (j - (i + 1))
      have hji : i + 1 + (j - (i + 1)) = j := by omega
      rw [hji] at hsum
      rw [hsum]
      refine le_trans (Finset.abs_sum_le_sum_abs _ _) ?_
      have hterm : ∀ k ∈ Finset.range (j - (i + 1)),
          |f.acf (i + 1 + k + 1) - f.acf (i + 1 + k)| ≤ L * f.D0 := by
        intro k _
        have hslope := hL (i + 1 + k)
        unfold GridElement.cellSlope at hslope
        rw [abs_div, abs_of_pos f.D0_pos] at hslope
        exact (div_le_iff₀ f.D0_pos).mp hslope
      refine le_trans (Finset.sum_le_sum hterm) ?_
      simp only [Finset.sum_const, nsmul_eq_mul, Finset.card_range]
      have hcast :
          ((j : ℝ) - ((i : ℝ) + 1)) * f.D0 =
            ((j - (i + 1) : ℕ) : ℝ) * f.D0 := by
        have : (j : ℝ) - ((i : ℝ) + 1) = ((j - (i + 1) : ℕ) : ℝ) := by
          rw [Nat.cast_sub hle, Nat.cast_add, Nat.cast_one]
        rw [this]
      have hspan :
          (j : ℝ) * f.D0 - ((i : ℝ) + 1) * f.D0 =
            ((j : ℝ) - ((i : ℝ) + 1)) * f.D0 := by ring
      rw [hspan, hcast]
      linarith
    have htri :=
      abs_sub_le (f.toFun t) (f.toFun ((j : ℝ) * f.D0)) (f.toFun s)
    have htri2 :=
      abs_sub_le (f.toFun ((j : ℝ) * f.D0))
        (f.toFun (((i : ℝ) + 1) * f.D0)) (f.toFun s)
    nlinarith [h1, h2, h3, htri, htri2]

/-- A uniform cell-slope bound is a Lipschitz constant for `toFun`. -/
lemma GridElement.lipschitz_toFun_of_slope_bound (f : GridElement) {L : ℝ}
    (hL0 : 0 ≤ L) (hL : ∀ d : ℕ, |f.cellSlope d| ≤ L) :
    LipschitzWith ⟨L, hL0⟩ f.toFun := by
  refine LipschitzWith.of_dist_le_mul fun x y => ?_
  rw [Real.dist_eq, Real.dist_eq]
  change |f.toFun x - f.toFun y| ≤ L * |x - y|
  rw [f.toFun_eq_toFun_abs x, f.toFun_eq_toFun_abs y]
  have habs : |(|x| - |y|)| ≤ |x - y| := abs_abs_sub_abs_le_abs_sub _ _
  refine le_trans ?_ (mul_le_mul_of_nonneg_left habs hL0)
  rcases le_total |x| |y| with hst | hts
  · have hlip := f.toFun_nonneg_lip hL (abs_nonneg x) hst
    rw [abs_sub_comm] at hlip
    convert hlip using 1
    rw [abs_sub_comm (|x|) (|y|), abs_of_nonneg (sub_nonneg.mpr hst)]
  · have hlip := f.toFun_nonneg_lip hL (abs_nonneg y) hts
    convert hlip using 1
    rw [abs_of_nonneg (sub_nonneg.mpr hts)]

/-- Autocorrelation of a Lipschitz compactly-supported witness is
itself Lipschitz, with constant `K² R²`. -/
lemma FullWeilTest.lipschitz_toFun (F : FullWeilTest) :
    ∃ K : NNReal, LipschitzWith K F.toFun := by
  obtain ⟨h, _, hhKs, hhLip, hhSupp, hac⟩ := F.autocorrelation
  obtain ⟨K, hK⟩ := hhLip
  obtain ⟨a, hsupp⟩ := hhSupp
  set R := F.supportRadius
  have hR : 0 ≤ R := F.supportRadius_nonneg
  have hc : Continuous h := hK.continuous
  refine ⟨⟨(K : ℝ) * K * R * R, by positivity⟩, ?_⟩
  refine LipschitzWith.of_dist_le_mul fun u v => ?_
  rw [Real.dist_eq, Real.dist_eq, hac u, hac v]
  rw [autocorrelation_eq_intervalIntegral hc hhKs hsupp hR,
    autocorrelation_eq_intervalIntegral hc hhKs hsupp hR]
  have hφu :
      IntervalIntegrable (fun t : ℝ => h t * h (t + u))
        MeasureTheory.volume a (a + R) :=
    (hc.mul (hc.comp (continuous_add_right u))).intervalIntegrable _ _
  have hφv :
      IntervalIntegrable (fun t : ℝ => h t * h (t + v))
        MeasureTheory.volume a (a + R) :=
    (hc.mul (hc.comp (continuous_add_right v))).intervalIntegrable _ _
  rw [← intervalIntegral.integral_sub hφu hφv]
  have hrew :
      (fun t : ℝ => h t * h (t + u) - h t * h (t + v)) =
        fun t : ℝ => h t * (h (t + u) - h (t + v)) := by
    funext t
    ring
  rw [hrew]
  have hab : a ≤ a + R := le_add_of_nonneg_right hR
  have hpt : ∀ t : ℝ,
      |h t * (h (t + u) - h (t + v))| ≤
        ((K : ℝ) * R) * (K : ℝ) * |u - v| := fun t => by
    rw [abs_mul]
    have hb := lipschitz_witness_abs_le hK hsupp hR t
    have hdu := hK.dist_le_mul (t + u) (t + v)
    rw [Real.dist_eq, Real.dist_eq] at hdu
    have hspan : |t + u - (t + v)| = |u - v| := by
      congr 1
      ring
    rw [hspan] at hdu
    have h1 : |h t| * |h (t + u) - h (t + v)| ≤
        ((K : ℝ) * R) * |h (t + u) - h (t + v)| :=
      mul_le_mul_of_nonneg_right hb (abs_nonneg _)
    have h2 : ((K : ℝ) * R) * |h (t + u) - h (t + v)| ≤
        ((K : ℝ) * R) * ((K : ℝ) * |u - v|) :=
      mul_le_mul_of_nonneg_left hdu (by positivity)
    linarith
  have hconst :
      intervalIntegral
          (fun _ : ℝ => ((K : ℝ) * R) * (K : ℝ) * |u - v|)
          a (a + R) MeasureTheory.volume =
        ((K : ℝ) * K * R * R) * |u - v| := by
    rw [intervalIntegral.integral_const, smul_eq_mul]
    ring
  rw [← Real.norm_eq_abs]
  refine le_trans
    (intervalIntegral.norm_integral_le_of_norm_le hab
      (by
        filter_upwards [] with t _ht
        simpa [Real.norm_eq_abs] using hpt t)
      intervalIntegrable_const)
    (le_of_eq hconst)

/-- `w(u) e^{-3u/2} = 2 e^{-2u}/(1-e^{-2u})` for `u ≠ 0`. -/
lemma weilArchUWeight_mul_exp {u : ℝ} (hu : u ≠ 0) :
    weilArchUWeight u * Real.exp (-(3 / 2 : ℝ) * u) =
      2 * Real.exp (-2 * u) / (1 - Real.exp (-2 * u)) := by
  unfold weilArchUWeight
  have hden : 1 - Real.exp (-2 * u) ≠ 0 := by
    intro h
    have hexp : Real.exp (-2 * u) = 1 := (sub_eq_zero.mp h).symm
    have hu2 : -2 * u = 0 := by
      rwa [← Real.exp_eq_one_iff]
    exact hu (by nlinarith [hu2])
  field_simp [hden]
  rw [← Real.exp_add]
  ring_nf

/-- The regularized weight is positive on `(0, ∞)`. -/
lemma weilArchUWeight_pos {u : ℝ} (hu : 0 < u) : 0 < weilArchUWeight u := by
  unfold weilArchUWeight
  have hnum : 0 < 2 * Real.exp (-u / 2) := by positivity
  have hden : 0 < 1 - Real.exp (-2 * u) :=
    sub_pos.mpr (Real.exp_lt_one_iff.2 (by linarith))
  positivity

/-- Dini majorant: `w(u)·u ≤ e^{3u/2}` on `(0, ∞)`. -/
lemma weilArchUWeight_mul_id_le {u : ℝ} (hu : 0 < u) :
    weilArchUWeight u * u ≤ Real.exp ((3 / 2 : ℝ) * u) := by
  have hden : 0 < 1 - Real.exp (-2 * u) :=
    sub_pos.mpr (Real.exp_lt_one_iff.2 (by linarith))
  have hge : (2 * u) * Real.exp (-2 * u) ≤ 1 - Real.exp (-2 * u) := by
    have hexp := Real.add_one_le_exp (2 * u)
    have hmul := mul_le_mul_of_nonneg_right hexp (Real.exp_pos (-2 * u)).le
    rw [add_mul, ← Real.exp_add, one_mul] at hmul
    have hsum : (2 * u) + (-2 * u) = 0 := by ring
    rw [hsum, Real.exp_zero] at hmul
    linarith
  unfold weilArchUWeight
  have hpos : 0 < 2 * u * Real.exp (-2 * u) := by positivity
  have hnum0 : 0 ≤ 2 * Real.exp (-u / 2) * u := by positivity
  have hw :
      2 * Real.exp (-u / 2) / (1 - Real.exp (-2 * u)) * u =
        2 * Real.exp (-u / 2) * u / (1 - Real.exp (-2 * u)) := by
    ring
  rw [hw]
  have hquot :
      2 * Real.exp (-u / 2) * u / (1 - Real.exp (-2 * u)) ≤
        2 * Real.exp (-u / 2) * u / (2 * u * Real.exp (-2 * u)) :=
    div_le_div_of_nonneg_left hnum0 hpos hge
  refine le_trans hquot ?_
  have hu0 : u ≠ 0 := hu.ne'
  have hne : Real.exp (-2 * u) ≠ 0 := (Real.exp_pos _).ne'
  have hsimp :
      2 * Real.exp (-u / 2) * u / (2 * u * Real.exp (-2 * u)) =
        Real.exp ((3 / 2 : ℝ) * u) := by
    field_simp [hu0, hne]
    rw [← Real.exp_add]
    congr 1
    ring
  exact le_of_eq hsimp

/-- Lipschitz increment against the tilted constant term. -/
lemma archDiff_le_of_lipschitz {g : ℝ → ℝ} {K : NNReal} {u : ℝ}
    (hK : LipschitzWith K g) (hu : 0 ≤ u) :
    |Real.exp (-(3 / 2 : ℝ) * u) * g 0 - g u| ≤
      ((3 / 2 : ℝ) * |g 0| + (K : ℝ)) * u := by
  have hexp : |Real.exp (-(3 / 2 : ℝ) * u) - 1| ≤ (3 / 2 : ℝ) * u := by
    have ha : 0 ≤ (3 / 2 : ℝ) * u := by positivity
    have hle : -((3 / 2 : ℝ) * u) + 1 ≤ Real.exp (-((3 / 2 : ℝ) * u)) :=
      Real.add_one_le_exp (-((3 / 2 : ℝ) * u))
    have hnp : Real.exp (-(3 / 2 : ℝ) * u) ≤ 1 := by
      refine Real.exp_le_one_iff.2 ?_
      nlinarith [ha]
    have hrew : Real.exp (-(3 / 2 : ℝ) * u) = Real.exp (-((3 / 2 : ℝ) * u)) := by
      ring_nf
    rw [hrew] at hnp ⊢
    rw [abs_of_nonpos (sub_nonpos.mpr hnp), neg_sub]
    linarith
  have hLip := hK.dist_le_mul 0 u
  rw [Real.dist_eq, Real.dist_eq, zero_sub, abs_neg, abs_of_nonneg hu] at hLip
  have hrew :
      Real.exp (-(3 / 2 : ℝ) * u) * g 0 - g u =
        (Real.exp (-(3 / 2 : ℝ) * u) - 1) * g 0 + (g 0 - g u) := by
    ring
  rw [hrew]
  refine le_trans (abs_add_le _ _) ?_
  rw [abs_mul]
  have h1 : |Real.exp (-(3 / 2 : ℝ) * u) - 1| * |g 0| ≤
      ((3 / 2 : ℝ) * u) * |g 0| :=
    mul_le_mul_of_nonneg_right hexp (abs_nonneg _)
  linarith [h1, hLip]

/-- Pointwise Dini bound for a Lipschitz test against the r475 weight. -/
lemma archIntegrand_dini_bound {g : ℝ → ℝ} {K : NNReal} {u b : ℝ}
    (hK : LipschitzWith K g) (hu : 0 < u) (hub : u ≤ b) :
    |weilArchUWeight u * (Real.exp (-(3 / 2 : ℝ) * u) * g 0 - g u)| ≤
      Real.exp ((3 / 2 : ℝ) * b) *
        ((3 / 2 : ℝ) * |g 0| + (K : ℝ)) := by
  have hw0 : 0 ≤ weilArchUWeight u := (weilArchUWeight_pos hu).le
  have hwu := weilArchUWeight_mul_id_le hu
  have hdiff := archDiff_le_of_lipschitz hK hu.le
  have hexp : Real.exp ((3 / 2 : ℝ) * u) ≤ Real.exp ((3 / 2 : ℝ) * b) :=
    Real.exp_le_exp.2 (by linarith)
  calc
    |weilArchUWeight u * (Real.exp (-(3 / 2 : ℝ) * u) * g 0 - g u)|
        = weilArchUWeight u *
            |Real.exp (-(3 / 2 : ℝ) * u) * g 0 - g u| := by
          rw [abs_mul, abs_of_nonneg hw0]
      _ ≤ weilArchUWeight u * (((3 / 2 : ℝ) * |g 0| + (K : ℝ)) * u) :=
          mul_le_mul_of_nonneg_left hdiff hw0
      _ = (weilArchUWeight u * u) * ((3 / 2 : ℝ) * |g 0| + (K : ℝ)) := by
          ring
      _ ≤ Real.exp ((3 / 2 : ℝ) * u) *
            ((3 / 2 : ℝ) * |g 0| + (K : ℝ)) :=
          mul_le_mul_of_nonneg_right hwu (by positivity)
      _ ≤ Real.exp ((3 / 2 : ℝ) * b) *
            ((3 / 2 : ℝ) * |g 0| + (K : ℝ)) :=
          mul_le_mul_of_nonneg_right hexp (by positivity)

/-- Dyadic sampled ACF interpolants share one Lipschitz constant. -/
lemma dyadicSampleGrid_uniform_lipschitz
    {F : FullWeilTest} {h : ℝ → ℝ} {K : NNReal} {a : ℝ}
    (hK : LipschitzWith K h) (hhKs : HasCompactSupport h)
    (hsupp : ∀ t : ℝ, t < a ∨ a + F.supportRadius < t → h t = 0)
    (hac : ∀ u : ℝ, F.toFun u =
      ∫ t : ℝ, h t * h (t + u) ∂MeasureTheory.volume) :
    ∃ Kgrid : NNReal, ∀ m,
      LipschitzWith Kgrid
        (dyadicSampleGrid h a F.supportRadius m).toFun := by
  obtain ⟨KF, hKF⟩ := F.lipschitz_toFun
  set R := F.supportRadius
  let L : ℝ := (KF : ℝ) + 8 * ((K : ℝ) * R) ^ 2
  have hL0 : 0 ≤ L := by positivity
  refine ⟨⟨L, hL0⟩, fun m => ?_⟩
  refine GridElement.lipschitz_toFun_of_slope_bound _ hL0 fun d => ?_
  set f := dyadicSampleGrid h a R m
  have hD : 0 < f.D0 := f.D0_pos
  by_cases hd : f.steps ≤ d
  · rw [f.cellSlope_eq_zero_of_steps_le hd, abs_zero]
    exact hL0
  · have h1 := dyadicSampleGrid_acf_sub_toFun_le hK hhKs hsupp hac m d
    have h2 := dyadicSampleGrid_acf_sub_toFun_le hK hhKs hsupp hac m (d + 1)
    have hF := hKF.dist_le_mul (((d : ℝ) + 1) * f.D0) ((d : ℝ) * f.D0)
    rw [Real.dist_eq, Real.dist_eq] at hF
    have hspan : |(((d : ℝ) + 1) * f.D0) - ((d : ℝ) * f.D0)| = f.D0 := by
      rw [add_mul, one_mul, add_sub_cancel_left, abs_of_pos hD]
    rw [hspan] at hF
    unfold GridElement.cellSlope
    rw [abs_div, abs_of_pos hD]
    have hDmesh : f.D0 = meshWidth m := rfl
    have h1' : |f.acf d - F.toFun ((d : ℝ) * f.D0)| ≤
        4 * ((K : ℝ) * R) ^ 2 * f.D0 := by
      simpa [f, hDmesh, R] using h1
    have h2' : |f.acf (d + 1) - F.toFun (((d : ℝ) + 1) * f.D0)| ≤
        4 * ((K : ℝ) * R) ^ 2 * f.D0 := by
      simpa [f, hDmesh, R, Nat.cast_succ] using h2
    have h1c : |F.toFun ((d : ℝ) * f.D0) - f.acf d| ≤
        4 * ((K : ℝ) * R) ^ 2 * f.D0 := by
      rw [abs_sub_comm]
      exact h1'
    have htri1 :=
      abs_sub_le (f.acf (d + 1)) (F.toFun (((d : ℝ) + 1) * f.D0)) (f.acf d)
    have htri2 :=
      abs_sub_le (F.toFun (((d : ℝ) + 1) * f.D0))
        (F.toFun ((d : ℝ) * f.D0)) (f.acf d)
    have hdiff : |f.acf (d + 1) - f.acf d| ≤ L * f.D0 := by
      have hsum :
          |f.acf (d + 1) - f.acf d| ≤
            |f.acf (d + 1) - F.toFun (((d : ℝ) + 1) * f.D0)| +
              |F.toFun (((d : ℝ) + 1) * f.D0) - F.toFun ((d : ℝ) * f.D0)| +
                |F.toFun ((d : ℝ) * f.D0) - f.acf d| := by
        linarith [htri1, htri2]
      have hnum :
          |f.acf (d + 1) - F.toFun (((d : ℝ) + 1) * f.D0)| +
            |F.toFun (((d : ℝ) + 1) * f.D0) - F.toFun ((d : ℝ) * f.D0)| +
              |F.toFun ((d : ℝ) * f.D0) - f.acf d| ≤
            8 * ((K : ℝ) * R) ^ 2 * f.D0 + (KF : ℝ) * f.D0 := by
        linarith [h2', hF, h1c]
      have hLrw : 8 * ((K : ℝ) * R) ^ 2 * f.D0 + (KF : ℝ) * f.D0 = L * f.D0 := by
        simp only [L]
        ring
      linarith [hsum, hnum, hLrw]
    exact (div_le_iff₀ hD).mpr hdiff

/-- r493d: the explicit left-sampled grid converges uniformly in the
autocorrelation and in `L¹` on the fixed support.  The witness is
already Lipschitz with interval support, so Heine plus a Lipschitz
Riemann-sum comparison of the ACF knots close both limits. -/
theorem fullWeil_dyadic_sample_convergence :
    FullWeilDyadicSampleConvergence := by
  intro F _hF
  obtain ⟨h, hhLp, hhKs, hhLip, hhSupp, hac⟩ := F.autocorrelation
  obtain ⟨K, hK⟩ := hhLip
  obtain ⟨a, hsupp⟩ := hhSupp
  refine ⟨h, K, a, hhLp, hhKs, hK, hsupp, hac, ?_⟩
  refine ⟨fun m =>
      dyadicSampleGrid_supportBound_le h a F.supportRadius m
        F.supportRadius_nonneg,
    ?unif, fun m =>
      (dyadicSampleGrid h a F.supportRadius m).intervalIntegrable_toFun
        (-F.supportRadius) F.supportRadius,
    ?l1, ?lip⟩
  · intro ε hε
    exact dyadicSampleGrid_toFun_uniform hK hhKs hsupp hac hε
  · exact dyadicSampleGrid_l1_tendsto hK hhKs hsupp hac
  · exact dyadicSampleGrid_uniform_lipschitz hK hhKs hsupp hac

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

/-- `w` is continuous on any set that avoids the origin. -/
lemma continuousOn_weilArchUWeight {s : Set ℝ}
    (hs : ∀ x ∈ s, x ≠ 0) : ContinuousOn weilArchUWeight s := by
  unfold weilArchUWeight
  refine ContinuousOn.div (by fun_prop) (by fun_prop) ?_
  intro x hx
  have hx0 : x ≠ 0 := hs x hx
  have hexp : Real.exp (-2 * x) ≠ 1 := by
    intro h
    have : -2 * x = 0 := by
      rwa [← Real.exp_eq_one_iff]
    exact hx0 (by nlinarith [this])
  exact sub_ne_zero.2 hexp.symm

/-- The tilted weight is continuous away from the origin. -/
lemma continuousOn_weilArchUWeight_mul_exp {s : Set ℝ}
    (hs : ∀ x ∈ s, x ≠ 0) :
    ContinuousOn
      (fun u : ℝ => weilArchUWeight u * Real.exp (-(3 / 2 : ℝ) * u)) s :=
  (continuousOn_weilArchUWeight hs).mul (by fun_prop)

/-- Measurability of the regularized integrand of a measurable test. -/
lemma measurable_weilArchUIntegrand_of_measurable {f : ℝ → ℝ}
    (hf : Measurable f) :
    Measurable (fun u : ℝ =>
      if u = 0 then (0 : ℝ)
      else weilArchUWeight u *
        (Real.exp (-(3 / 2 : ℝ) * u) * f 0 - f u)) := by
  refine Measurable.ite (measurableSet_singleton (0 : ℝ))
    measurable_const ?_
  unfold weilArchUWeight
  fun_prop

/-- Derivative of `log(1-e^{-2u})` is the tilted r475 weight. -/
lemma hasDerivAt_log_one_sub_exp_neg_two {x : ℝ} (hx : 0 < x) :
    HasDerivAt (fun t : ℝ => Real.log (1 - Real.exp (-2 * t)))
      (weilArchUWeight x * Real.exp (-(3 / 2 : ℝ) * x)) x := by
  have harg : 0 < 1 - Real.exp (-2 * x) :=
    sub_pos.mpr (Real.exp_lt_one_iff.2 (by linarith))
  have hder : HasDerivAt (fun t : ℝ => 1 - Real.exp (-2 * t))
      (2 * Real.exp (-2 * x)) x := by
    convert (hasDerivAt_const x (1 : ℝ)).sub
      ((Real.hasDerivAt_exp (-2 * x)).comp x
        ((hasDerivAt_id x).const_mul (-2))) using 1
    ring
  have hlog := (Real.hasDerivAt_log (ne_of_gt harg)).comp x hder
  convert hlog using 1
  · rw [weilArchUWeight_mul_exp hx.ne']
    field_simp [harg.ne']

/-- Antiderivative of the tilted weight on a compact interval away from `0`. -/
lemma weilArch_weight_exp_integral {b R : ℝ} (hb : 0 < b) (hle : b ≤ R) :
    intervalIntegral
        (fun u : ℝ =>
          weilArchUWeight u * Real.exp (-(3 / 2 : ℝ) * u))
        b R MeasureTheory.volume =
      Real.log (1 - Real.exp (-2 * R)) -
        Real.log (1 - Real.exp (-2 * b)) := by
  have hcont :
      ContinuousOn
        (fun u : ℝ =>
          weilArchUWeight u * Real.exp (-(3 / 2 : ℝ) * u))
        (Set.uIcc b R) := by
    refine continuousOn_weilArchUWeight_mul_exp ?_
    intro x hx
    rw [Set.uIcc_of_le hle] at hx
    exact (lt_of_lt_of_le hb hx.1).ne'
  apply intervalIntegral.integral_eq_sub_of_hasDerivAt
  · intro x hx
    have hxpos : 0 < x := by
      rw [Set.uIcc_of_le hle] at hx
      exact lt_of_lt_of_le hb hx.1
    exact hasDerivAt_log_one_sub_exp_neg_two hxpos
  · exact hcont.intervalIntegrable

/-- The r475 integrand of a Lipschitz test is interval-integrable
on any compact subinterval of `[0, ∞)`. -/
lemma intervalIntegrable_weilArchUIntegrand
    {f : ℝ → ℝ} {K : NNReal} {a b : ℝ}
    (hf : LipschitzWith K f) (ha : 0 ≤ a) (hab : a ≤ b) :
    IntervalIntegrable
      (fun u : ℝ =>
        if u = 0 then (0 : ℝ)
        else weilArchUWeight u *
          (Real.exp (-(3 / 2 : ℝ) * u) * f 0 - f u))
      MeasureTheory.volume a b := by
  let g : ℝ → ℝ := fun _ =>
    Real.exp ((3 / 2 : ℝ) * b) * ((3 / 2 : ℝ) * |f 0| + (K : ℝ))
  refine IntervalIntegrable.mono_fun' (g := g) intervalIntegrable_const ?meas ?bound
  · exact (measurable_weilArchUIntegrand_of_measurable
      hf.continuous.measurable).aestronglyMeasurable
  · refine (MeasureTheory.ae_restrict_iff' measurableSet_uIoc).2 ?_
    refine Filter.Eventually.of_forall fun u hu => ?_
    have huI : u ∈ Set.uIcc a b := Set.uIoc_subset_uIcc hu
    have hu0 : 0 ≤ u := by
      rw [Set.uIcc_of_le hab] at huI
      exact le_trans ha huI.1
    have hub : u ≤ b := by
      rw [Set.uIcc_of_le hab] at huI
      exact huI.2
    by_cases h0 : u = 0
    · simp [h0, g, Real.norm_eq_abs]
      positivity
    · have hupos : 0 < u := lt_of_le_of_ne hu0 (Ne.symm h0)
      simpa [h0, g, Real.norm_eq_abs] using
        archIntegrand_dini_bound hf hupos hub

/-- Cutoff invariance: if `toFun` vanishes past `supportBound ≤ R`,
the pairing at `supportBound` equals the pairing at `R`. -/
lemma weilArchSide_eq_cutoff (f : GridElement) {R : ℝ} {K : NNReal}
    (hR : 0 < R) (hsupp : f.supportBound ≤ R)
    (hK : LipschitzWith K f.toFun) :
    weilArchSide f =
      (-(Real.eulerMascheroniConstant + Real.log Real.pi)
        - Real.log (1 - Real.exp (-2 * R))) * f.toFun 0 +
        intervalIntegral (weilArchUIntegrand f) 0 R
          MeasureTheory.volume := by
  by_cases hb0 : 0 < f.supportBound
  · set b := f.supportBound
    have hbR : b ≤ R := hsupp
    have hI0b : IntervalIntegrable (weilArchUIntegrand f)
        MeasureTheory.volume 0 b :=
      intervalIntegrable_weilArchUIntegrand hK le_rfl f.supportBound_nonneg
    have hIbR : IntervalIntegrable (weilArchUIntegrand f)
        MeasureTheory.volume b R :=
      intervalIntegrable_weilArchUIntegrand hK hb0.le hbR
    have hadd :
        intervalIntegral (weilArchUIntegrand f) 0 R
            MeasureTheory.volume =
          intervalIntegral (weilArchUIntegrand f) 0 b
            MeasureTheory.volume +
          intervalIntegral (weilArchUIntegrand f) b R
            MeasureTheory.volume :=
      (intervalIntegral.integral_add_adjacent_intervals hI0b hIbR).symm
    have htail :
        intervalIntegral (weilArchUIntegrand f) b R
            MeasureTheory.volume =
          f.toFun 0 *
            intervalIntegral
              (fun u : ℝ =>
                weilArchUWeight u * Real.exp (-(3 / 2 : ℝ) * u))
              b R MeasureTheory.volume := by
      have heq : Set.EqOn (weilArchUIntegrand f)
          (fun u : ℝ =>
            weilArchUWeight u * Real.exp (-(3 / 2 : ℝ) * u) * f.toFun 0)
          (Set.uIcc b R) := by
        intro u hu
        have hu0 : u ≠ 0 := by
          rw [Set.uIcc_of_le hbR] at hu
          exact (lt_of_lt_of_le hb0 hu.1).ne'
        unfold weilArchUIntegrand
        rw [if_neg hu0]
        have htu : f.toFun u = 0 := by
          rw [Set.uIcc_of_le hbR] at hu
          exact f.toFun_eq_zero_of_supportBound_le
            (le_trans hu.1 (le_abs_self u))
        rw [htu, sub_zero]
        ring
      have hcongr :=
        intervalIntegral.integral_congr (μ := MeasureTheory.volume) heq
      have hmul :
          (fun u : ℝ =>
            weilArchUWeight u * Real.exp (-(3 / 2 : ℝ) * u) * f.toFun 0) =
            fun u : ℝ =>
              f.toFun 0 *
                (weilArchUWeight u * Real.exp (-(3 / 2 : ℝ) * u)) := by
        funext u
        ring
      rw [hcongr, hmul, intervalIntegral.integral_const_mul]
    unfold weilArchSide
    rw [if_pos (by simpa [b] using hb0)]
    rw [hadd, htail, weilArch_weight_exp_integral hb0 hbR]
    ring
  · have hb : f.supportBound = 0 :=
      le_antisymm (not_lt.mp hb0) f.supportBound_nonneg
    unfold weilArchSide
    rw [if_neg hb0]
    have hto0 : f.toFun 0 = 0 := by
      exact f.toFun_eq_zero_of_supportBound_le (by simpa [hb] using abs_nonneg (0 : ℝ))
    have hI0 : intervalIntegral (weilArchUIntegrand f) 0 R
        MeasureTheory.volume = 0 := by
      have heq : Set.EqOn (weilArchUIntegrand f) (fun _ => (0 : ℝ))
          (Set.uIcc 0 R) := by
        intro u hu
        unfold weilArchUIntegrand
        split_ifs with h0
        · rfl
        · have htu : f.toFun u = 0 :=
            f.toFun_eq_zero_of_supportBound_le (by simpa [hb] using abs_nonneg u)
          have ht00 : f.toFun 0 = 0 := hto0
          rw [ht00, htu]
          ring
      rw [intervalIntegral.integral_congr (μ := MeasureTheory.volume) heq,
        intervalIntegral.integral_zero]
    rw [hto0, hI0]
    ring

/-- r493e: the r475 u-space pairing is continuous along a fixed-support
approximation that carries a common Lipschitz / Dini majorant. -/
theorem fullWeil_arch_side_tendsto
    (F : FullWeilTest) (grid : ℕ → GridElement) :
    F.admissible → F.FixedSupportGridApproximation grid →
      Filter.Tendsto (fun n => weilArchSide (grid n))
        Filter.atTop (nhds (fullWeilArchSide F)) := by
  intro _hF happrox
  obtain ⟨KF, hKF⟩ := F.lipschitz_toFun
  obtain ⟨Kgrid, hKgrid⟩ := happrox.2.2.2.2
  set R := F.supportRadius
  have hR0 : 0 ≤ R := F.supportRadius_nonneg
  by_cases hRz : R = 0
  · have hF0 : fullWeilArchSide F = 0 := by
      unfold fullWeilArchSide
      have : ¬0 < F.supportRadius := by
        simpa [hRz, R] using (lt_irrefl (0 : ℝ))
      simp [this]
    have hn0 : ∀ n, weilArchSide (grid n) = 0 := fun n => by
      unfold weilArchSide
      have : ¬0 < (grid n).supportBound :=
        not_lt.mpr (le_trans (happrox.1 n) (le_of_eq hRz))
      simp [this]
    simpa [hF0, funext hn0] using tendsto_const_nhds (x := (0 : ℝ))
  have hRpos : 0 < R := lt_of_le_of_ne hR0 (Ne.symm hRz)
  have hFeq :
      fullWeilArchSide F =
        (-(Real.eulerMascheroniConstant + Real.log Real.pi)
          - Real.log (1 - Real.exp (-2 * R))) * F.toFun 0 +
          intervalIntegral (fullWeilArchUIntegrand F) 0 R
            MeasureTheory.volume := by
    unfold fullWeilArchSide
    rw [if_pos (by simpa [R] using hRpos)]
  have hgeq : ∀ n,
      weilArchSide (grid n) =
        (-(Real.eulerMascheroniConstant + Real.log Real.pi)
          - Real.log (1 - Real.exp (-2 * R))) * (grid n).toFun 0 +
          intervalIntegral (weilArchUIntegrand (grid n)) 0 R
            MeasureTheory.volume :=
    fun n => weilArchSide_eq_cutoff (grid n) hRpos (happrox.1 n) (hKgrid n)
  have hfun :
      (fun n => weilArchSide (grid n)) =
        fun n =>
          (-(Real.eulerMascheroniConstant + Real.log Real.pi)
            - Real.log (1 - Real.exp (-2 * R))) * (grid n).toFun 0 +
            intervalIntegral (weilArchUIntegrand (grid n)) 0 R
              MeasureTheory.volume := funext hgeq
  rw [hfun, hFeq]
  refine Filter.Tendsto.add ?_ ?_
  · exact tendsto_const_nhds.mul (happrox.tendsto_toFun 0)
  · let C : ℝ :=
      Real.exp ((3 / 2 : ℝ) * R) *
        ((3 / 2 : ℝ) * (|F.toFun 0| + 1) + max (KF : ℝ) (Kgrid : ℝ))
    let bound : ℝ → ℝ := fun _ => C
    have hg0 : ∀ᶠ n in Filter.atTop,
        |(grid n).toFun 0| ≤ |F.toFun 0| + 1 := by
      have htend := Metric.tendsto_nhds.1 (happrox.tendsto_toFun 0) 1
        (by norm_num)
      filter_upwards [htend] with n hn
      have hlt : |(grid n).toFun 0 - F.toFun 0| < 1 := by
        simpa [Real.dist_eq] using hn
      have htri : |(grid n).toFun 0| ≤
          |F.toFun 0| + |(grid n).toFun 0 - F.toFun 0| := by
        have := abs_add_le ((grid n).toFun 0 - F.toFun 0) (F.toFun 0)
        convert this using 1 <;> ring
      exact le_trans htri (add_le_add_right hlt.le _)
    have hIgrid (n : ℕ) :
        IntervalIntegrable (weilArchUIntegrand (grid n))
          MeasureTheory.volume 0 R :=
      intervalIntegrable_weilArchUIntegrand (hKgrid n) le_rfl hR0
    have hIF :
        IntervalIntegrable (fullWeilArchUIntegrand F)
          MeasureTheory.volume 0 R :=
      intervalIntegrable_weilArchUIntegrand hKF le_rfl hR0
    have heq_int (φ : ℝ → ℝ) :
        intervalIntegral φ 0 R MeasureTheory.volume =
          ∫ u in Set.Ioc (0 : ℝ) R, φ u :=
      intervalIntegral.integral_of_le hR0
    simp_rw [heq_int]
    have hbound_int :
        MeasureTheory.Integrable bound
          (MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) R)) :=
      (intervalIntegrable_const (c := C)).1
    refine MeasureTheory.tendsto_integral_filter_of_dominated_convergence
      (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) R))
      bound ?meas ?dom hbound_int ?lim
    · filter_upwards with n
      exact (measurable_weilArchUIntegrand_of_measurable
        (hKgrid n).continuous.measurable).aestronglyMeasurable
    · filter_upwards [hg0] with n hn
      rw [MeasureTheory.ae_restrict_iff' measurableSet_Ioc]
      refine Filter.Eventually.of_forall fun u hu => ?_
      have hupos : 0 < u := hu.1
      have hub : u ≤ R := hu.2
      rw [Real.norm_eq_abs]
      have hpt := archIntegrand_dini_bound (hKgrid n) hupos hub
      have hIeq : weilArchUIntegrand (grid n) u =
          weilArchUWeight u *
            (Real.exp (-(3 / 2 : ℝ) * u) * (grid n).toFun 0 -
              (grid n).toFun u) := by
        unfold weilArchUIntegrand
        rw [if_neg hupos.ne']
      rw [hIeq]
      refine le_trans hpt ?_
      refine mul_le_mul_of_nonneg_left ?_ (Real.exp_pos _).le
      refine add_le_add ?_ (le_max_right _ _)
      exact mul_le_mul_of_nonneg_left hn (by norm_num)
    · rw [MeasureTheory.ae_restrict_iff' measurableSet_Ioc]
      refine Filter.Eventually.of_forall fun u hu => ?_
      have hu0 : u ≠ 0 := hu.1.ne'
      simp only [weilArchUIntegrand, fullWeilArchUIntegrand, hu0, ↓reduceIte]
      exact ((tendsto_const_nhds.mul (happrox.tendsto_toFun 0)).sub
        (happrox.tendsto_toFun u)).const_mul _

def FullWeilArchContinuity : Prop :=
  ∀ (F : FullWeilTest) (grid : ℕ → GridElement), F.admissible →
    F.FixedSupportGridApproximation grid →
      Filter.Tendsto (fun n => weilArchSide (grid n))
        Filter.atTop (nhds (fullWeilArchSide F))

theorem fullWeil_arch_continuity : FullWeilArchContinuity :=
  fullWeil_arch_side_tendsto

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

/-- **Single remaining dense-completion package (r493d shrink).**

The r376 pole hat dictionary and the r493d dyadic sampled-grid
`L²`/ACF transfer are proved theorems.  The remaining component is
r475 u-space / Dini arch continuity. -/
def FullWeilFixedSupportCompletion : Prop :=
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

/-- r493e: the remaining dense-completion package is the proved
Dini / u-space arch continuity. -/
theorem fullWeil_fixedSupport_completion :
    FullWeilFixedSupportCompletion :=
  fullWeil_arch_continuity

/-- Grid density is the proved r493d dyadic sampled-grid transfer.
Formerly the first projection `fullWeil_fixedSupport_completion.1`. -/
theorem fullWeil_fixedSupport_grid_density :
    FullWeilFixedSupportGridDensity :=
  fullWeilFixedSupportGridDensity_of_dyadicSample
    fullWeil_dyadic_sample_convergence

/-- All channel limits follow from the remaining arch completion
component, the proved pole hat identity, and the proved comb and
polar-integral continuity theorems. -/
theorem fullWeil_channel_continuity :
    FullWeilChannelContinuity :=
  fullWeilChannelContinuity_of_components
    fullWeil_fixedSupport_completion
    (gridPoleIntegralIdentification_of_hat gridPoleHatIntegralIdentity)

/-- Bridge 1 now consumes exactly the two named bricks above; its
positivity-transfer algebra is sorry-free. -/
theorem grid_dense_extension : GridDenseExtension :=
  grid_dense_extension_of_fixedSupport
    fullWeil_fixedSupport_grid_density fullWeil_channel_continuity

/-! ### r497: multiplicity API for `riemannZeta` (bridge [2], first brick)

Lane decision: prove the explicit-formula identification **before**
off-critical separation.  `FullWeilSeparatesOffCriticalZeros` is
stated against the opaque `standardExplicitFormula`.  That value
cannot be shown negative without an identification (bridge [2]) or a
spectral-side definition, which is the same work.  Mathlib's
`RiemannHypothesis` is the quantified zero statement
`ζ(s)=0 ∧ nontrivial ∧ s≠1 → Re s = 1/2`; the logical wrapper
`standard_weil_criterion_to_mathlib_rh_of_separation` is already
proved.  A contour/Mellin argument for one test family is a
special case of Guinand--Weil, not a bypass of [2].

Named remainder of [2] after this brick:
  [2b] meromorphic `ζ'/ζ` (simple poles, residue = order);
  [2c] local finiteness of zeros and absolute convergence of the
      spectral pairing for admissible tests;
  [2d] contour evaluation = the three-channel form
      (`standard_explicit_formula_identification`).
-/

/-- `ζ` is holomorphic on `ℂ \ {1}`. -/
lemma differentiableOn_riemannZeta_compl_one :
    DifferentiableOn ℂ riemannZeta ({1}ᶜ) :=
  fun _z hz => (differentiableAt_riemannZeta hz).differentiableWithinAt

/-- Analyticity on a neighbourhood of every point of `ℂ \ {1}`. -/
lemma analyticOnNhd_riemannZeta_compl_one :
    AnalyticOnNhd ℂ riemannZeta ({1}ᶜ) :=
  differentiableOn_riemannZeta_compl_one.analyticOnNhd isOpen_compl_singleton

/-- Pointwise analyticity of `ζ` away from the pole. -/
lemma analyticAt_riemannZeta {s : ℂ} (hs : s ≠ 1) :
    AnalyticAt ℂ riemannZeta s :=
  analyticOnNhd_riemannZeta_compl_one s hs

/-- `ℂ \ {1}` is connected (real rank of `ℂ` is 2). -/
lemma isConnected_compl_one : IsConnected ({1}ᶜ : Set ℂ) :=
  isConnected_compl_singleton_of_one_lt_rank
      (Complex.rank_real_complex ▸ Nat.one_lt_ofNat) (1 : ℂ)

/-- `ℂ \ {1}` is preconnected (real rank of `ℂ` is 2). -/
lemma isPreconnected_compl_one : IsPreconnected ({1}ᶜ : Set ℂ) :=
  isConnected_compl_one.isPreconnected

/-- Identity theorem: `ζ` is not eventually zero at any non-polar
point, because `ζ(0) = -1/2 ≠ 0` and `ℂ \ {1}` is preconnected. -/
lemma analyticOrderAt_riemannZeta_ne_top {s : ℂ} (hs : s ≠ 1) :
    analyticOrderAt riemannZeta s ≠ ⊤ := by
  have h0 : (0 : ℂ) ≠ 1 := by norm_num
  have hord0 : analyticOrderAt riemannZeta 0 ≠ ⊤ := by
    have hz0 : riemannZeta 0 ≠ 0 := by
      rw [riemannZeta_zero]
      norm_num
    have hord :
        analyticOrderAt riemannZeta 0 = 0 :=
      (analyticAt_riemannZeta h0).analyticOrderAt_eq_zero.2 hz0
    simp [hord]
  exact analyticOnNhd_riemannZeta_compl_one.analyticOrderAt_ne_top_of_isPreconnected
    isPreconnected_compl_one h0 hs hord0

/-- A genuine zero (away from the pole) has positive order. -/
lemma analyticOrderAt_riemannZeta_ne_zero {s : ℂ}
    (hz : riemannZeta s = 0) (hs : s ≠ 1) :
    analyticOrderAt riemannZeta s ≠ 0 :=
  (analyticAt_riemannZeta hs).analyticOrderAt_ne_zero.2 hz

/-- **r497 brick [2a].**  At every non-polar zero the analytic order
of `riemannZeta` is a positive finite `ℕ∞` (neither `0` nor `⊤`).
This is the multiplicity number that a future spectral sum must
carry.  Mathlib v4.29.1 has `analyticOrderAt` but no `ZetaZero`
enumeration. -/
lemma riemannZeta_analyticOrderAt_finite_pos {s : ℂ}
    (hz : riemannZeta s = 0) (hs : s ≠ 1) :
    analyticOrderAt riemannZeta s ≠ 0 ∧
      analyticOrderAt riemannZeta s ≠ ⊤ :=
  ⟨analyticOrderAt_riemannZeta_ne_zero hz hs,
    analyticOrderAt_riemannZeta_ne_top hs⟩

/-- Local multiplicity as a natural number (junk `0` at the pole). -/
noncomputable def riemannZetaMultiplicity (s : ℂ) : ℕ :=
  analyticOrderNatAt riemannZeta s

lemma riemannZetaMultiplicity_pos {s : ℂ}
    (hz : riemannZeta s = 0) (hs : s ≠ 1) :
    0 < riemannZetaMultiplicity s := by
  have hfin := riemannZeta_analyticOrderAt_finite_pos hz hs
  unfold riemannZetaMultiplicity analyticOrderNatAt
  refine Nat.pos_of_ne_zero ?_
  intro h0
  exact (not_or.mpr ⟨hfin.1, hfin.2⟩) (ENat.toNat_eq_zero.mp h0)

/-! ### r498: meromorphic `ζ'/ζ` at non-polar zeros (bridge [2b], first cut)

Standard local calculus, specialised to `riemannZeta`.  If `f` is
analytic of finite order `m` at `s`, then
`f'/f = m/(z-s) + analytic` on a punctured neighbourhood.  At a
genuine non-polar zero this is a simple pole whose residue equals
`riemannZetaMultiplicity`.  The polar point `s = 1` (intended
residue `-1`) is not in this cut.  Meromorphy on `ℂ \ {1}` is the
carrier a later contour ([2d]) will restrict to a strip rectangle.
-/

section LogDerivZeta

open Filter

open scoped Topology

/-- Logarithmic derivative of a translated monomial.  Used only on
`z ≠ s`; the identity is algebraic. -/
lemma logDeriv_sub_const_pow (s : ℂ) (n : ℕ) (z : ℂ) :
    logDeriv (fun w : ℂ => (w - s) ^ n) z = n / (z - s) := by
  have hcomp :
      (fun w : ℂ => (w - s) ^ n) =
        (fun y : ℂ => y ^ n) ∘ fun w : ℂ => w - s := rfl
  rw [hcomp, logDeriv_comp (by fun_prop) (by fun_prop)]
  simp [logDeriv_pow]

/-- Local factorization of the logarithmic derivative of a
non-identically-vanishing analytic germ: finite order `m` gives
`f'/f = m/(z-s) + h` with `h` analytic at `s`.  Harmless at regular
points (`m = 0`). -/
lemma exists_analytic_logDeriv_eq_order_div
    {f : ℂ → ℂ} {s : ℂ} (hf : AnalyticAt ℂ f s)
    (hfin : analyticOrderAt f s ≠ ⊤) :
    ∃ h : ℂ → ℂ, AnalyticAt ℂ h s ∧
      ∀ᶠ z in 𝓝[≠] s,
        logDeriv f z = (analyticOrderNatAt f s : ℂ) / (z - s) + h z := by
  obtain ⟨g, hg_an, hg_ne, hg_eq⟩ := hf.analyticOrderAt_ne_top.mp hfin
  set m := analyticOrderNatAt f s
  have hg_log : AnalyticAt ℂ (logDeriv g) s :=
    (hg_an.deriv).div hg_an hg_ne
  refine ⟨logDeriv g, hg_log, ?_⟩
  have hg_nz : ∀ᶠ z in 𝓝 s, g z ≠ 0 :=
    hg_an.continuousAt.preimage_mem_nhds
      (isOpen_compl_singleton.mem_nhds hg_ne)
  have hprod : f =ᶠ[𝓝 s] fun z => (z - s) ^ m * g z := by
    filter_upwards [hg_eq] with z hz
    simpa [smul_eq_mul] using hz
  obtain ⟨r, hr, hball⟩ := hg_an.exists_ball_analyticOnNhd
  rw [eventually_nhdsWithin_iff]
  filter_upwards [hprod.eventually_nhds, hg_nz,
    Metric.ball_mem_nhds s hr] with z hf_nhds hgz hzball hzne
  have hz : z ≠ s := hzne
  have hf_eq : f =ᶠ[𝓝 z] fun w => (w - s) ^ m * g w := hf_nhds
  have hpow_ne : (z - s) ^ m ≠ 0 :=
    pow_ne_zero _ (sub_ne_zero.mpr hz)
  have hlog :
      logDeriv f z = logDeriv (fun w => (w - s) ^ m * g w) z := by
    simp [logDeriv_apply, hf_eq.deriv_eq, hf_eq.eq_of_nhds]
  have hdf : DifferentiableAt ℂ (fun w : ℂ => (w - s) ^ m) z := by
    fun_prop
  have hdg : DifferentiableAt ℂ g z :=
    (hball z hzball).differentiableAt
  rw [hlog, logDeriv_mul z hpow_ne hgz hdf hdg,
    logDeriv_sub_const_pow s m z]

/-- **r498 brick [2b], lemma (1).**  At a non-polar zero,
`ζ'/ζ = m/(z-s) + analytic`, with `m = riemannZetaMultiplicity s`. -/
lemma logDeriv_riemannZeta_eq_multiplicity_div_add_analytic {s : ℂ}
    (_hz : riemannZeta s = 0) (hs : s ≠ 1) :
    ∃ h : ℂ → ℂ, AnalyticAt ℂ h s ∧
      ∀ᶠ z in 𝓝[≠] s,
        logDeriv riemannZeta z =
          (riemannZetaMultiplicity s : ℂ) / (z - s) + h z :=
  exists_analytic_logDeriv_eq_order_div
    (analyticAt_riemannZeta hs) (analyticOrderAt_riemannZeta_ne_top hs)

/-- The logarithmic derivative is meromorphic at every non-polar
point (zeros and regular points of `ζ`). -/
lemma meromorphicAt_logDeriv_riemannZeta {s : ℂ} (hs : s ≠ 1) :
    MeromorphicAt (logDeriv riemannZeta) s := by
  obtain ⟨h, hh, heq⟩ :=
    exists_analytic_logDeriv_eq_order_div
      (analyticAt_riemannZeta hs) (analyticOrderAt_riemannZeta_ne_top hs)
  refine MeromorphicAt.iff_eventuallyEq_zpow_smul_analyticAt.mpr ?_
  refine ⟨-1,
    fun z => (analyticOrderNatAt riemannZeta s : ℂ) + (z - s) * h z,
    by fun_prop, ?_⟩
  filter_upwards [heq, self_mem_nhdsWithin] with z hz hzne
  have hz0 : z - s ≠ 0 := sub_ne_zero.mpr hzne
  rw [hz, zpow_neg_one, smul_eq_mul]
  field_simp [hz0]

/-- Carrier for a later strip-rectangle contour: `ζ'/ζ` is meromorphic
on `ℂ \ {1}`. -/
lemma meromorphicOn_logDeriv_riemannZeta_compl_one :
    MeromorphicOn (logDeriv riemannZeta) ({1}ᶜ) :=
  fun _s hs => meromorphicAt_logDeriv_riemannZeta hs

/-- At a genuine non-polar zero the meromorphic order of `ζ'/ζ` is
exactly `-1` (simple pole). -/
lemma meromorphicOrderAt_logDeriv_riemannZeta {s : ℂ}
    (hz : riemannZeta s = 0) (hs : s ≠ 1) :
    meromorphicOrderAt (logDeriv riemannZeta) s = (-1 : ℤ) := by
  have hm : 0 < riemannZetaMultiplicity s :=
    riemannZetaMultiplicity_pos hz hs
  obtain ⟨h, hh, heq⟩ :=
    logDeriv_riemannZeta_eq_multiplicity_div_add_analytic hz hs
  set G : ℂ → ℂ :=
    fun z => (riemannZetaMultiplicity s : ℂ) + (z - s) * h z
  have hGan : AnalyticAt ℂ G s := by fun_prop
  have hGne : G s ≠ 0 := by
    simp [G]
    exact_mod_cast hm.ne'
  have hfmer : MeromorphicAt (logDeriv riemannZeta) s :=
    meromorphicAt_logDeriv_riemannZeta hs
  rw [meromorphicOrderAt_eq_int_iff hfmer]
  refine ⟨G, hGan, hGne, ?_⟩
  filter_upwards [heq, self_mem_nhdsWithin] with z hz hzne
  have hz0 : z - s ≠ 0 := sub_ne_zero.mpr hzne
  rw [hz, zpow_neg_one, smul_eq_mul]
  field_simp [hz0]
  rfl

/-- Residue form of lemma (1): `(z-s)·ζ'/ζ → m` at a non-polar zero. -/
lemma tendsto_mul_logDeriv_riemannZeta {s : ℂ}
    (hz : riemannZeta s = 0) (hs : s ≠ 1) :
    Tendsto (fun z => (z - s) * logDeriv riemannZeta z)
      (𝓝[≠] s) (𝓝 (riemannZetaMultiplicity s : ℂ)) := by
  obtain ⟨h, hh, heq⟩ :=
    logDeriv_riemannZeta_eq_multiplicity_div_add_analytic hz hs
  set G : ℂ → ℂ :=
    fun z => (riemannZetaMultiplicity s : ℂ) + (z - s) * h z
  have hagree :
      (fun z => (z - s) * logDeriv riemannZeta z) =ᶠ[𝓝[≠] s] G := by
    filter_upwards [heq, self_mem_nhdsWithin] with z hz hzne
    have hz0 : z - s ≠ 0 := sub_ne_zero.mpr hzne
    rw [hz]
    field_simp [hz0]
    rfl
  have hcont : ContinuousAt G s := by
    unfold G
    exact continuousAt_const.add
      ((continuousAt_id.sub continuousAt_const).mul hh.continuousAt)
  have hGs : G s = (riemannZetaMultiplicity s : ℂ) := by simp [G]
  rw [← hGs]
  exact Tendsto.congr' (EventuallyEq.symm hagree)
    (hcont.tendsto.mono_left nhdsWithin_le_nhds)

end LogDerivZeta

/-! ### r499: local finiteness of `ζ`-zeros (bridge [2c], first cut)

The identity theorem ([2a]) isolates zeros on `ℂ \ {1}`.  Combined
with Mathlib's `codiscreteWithin` / locally-finite complement, every
compact `K ⊆ ℂ \ {1}` meets the zero set finitely often.  The residue
at `s = 1` excludes zeros from a punctured neighbourhood of the pole,
so the same finiteness holds on every compact in `ℂ` (a closed strip
rectangle may contain `1`).  A residue sum of `ζ'/ζ` over such a
rectangle is therefore a finite sum, using the r498 residue lemma.
The absolute-convergence half of [2c] is not in this cut.
-/

section ZetaZeroFiniteness

open Filter Set

open scoped Topology

/-- Isolation: away from the pole, `ζ` is not identically zero, so it
is non-vanishing on a punctured neighbourhood. -/
lemma riemannZeta_eventually_ne_zero_punctured {s : ℂ} (hs : s ≠ 1) :
    ∀ᶠ z in 𝓝[≠] s, riemannZeta z ≠ 0 :=
  (analyticAt_riemannZeta hs).eventually_eq_zero_or_eventually_ne_zero.resolve_left
    fun h => analyticOrderAt_riemannZeta_ne_top hs (analyticOrderAt_eq_top.mpr h)

/-- The residue ` (s-1)ζ(s) → 1 ` forbids zeros in a punctured
neighbourhood of the pole. -/
lemma riemannZeta_eventually_ne_zero_nhdsNE_one :
    ∀ᶠ s in 𝓝[≠] 1, riemannZeta s ≠ 0 := by
  have h := riemannZeta_residue_one.eventually
    (isOpen_compl_singleton.mem_nhds (by norm_num : (1 : ℂ) ≠ 0))
  filter_upwards [h] with s hs
  exact right_ne_zero_of_mul hs

/-- Non-polar zeros of `ζ` are a discrete subset of `ℂ`. -/
lemma isDiscrete_riemannZeta_zeros_compl_one :
    IsDiscrete (riemannZeta ⁻¹' {0} ∩ ({1}ᶜ : Set ℂ)) :=
  isDiscrete_of_codiscreteWithin (s := riemannZeta ⁻¹' {0})
    (analyticOnNhd_riemannZeta_compl_one.preimage_zero_mem_codiscreteWithin
      (by
        rw [riemannZeta_zero]
        norm_num)
      (by norm_num : (0 : ℂ) ≠ 1) isConnected_compl_one)

/-- Local finiteness on `ℂ \ {1}`: every point has a neighbourhood
meeting only finitely many zeros. -/
lemma riemannZeta_zeros_locallyFinite_compl_one :
    ∀ z ∈ ({1}ᶜ : Set ℂ), ∃ t ∈ 𝓝 z,
      (t ∩ {w : ℂ | w ≠ 1 ∧ riemannZeta w = 0}).Finite := by
  have hcod :
      riemannZeta ⁻¹' {0}ᶜ ∈ codiscreteWithin ({1}ᶜ : Set ℂ) :=
    analyticOnNhd_riemannZeta_compl_one.preimage_zero_mem_codiscreteWithin
      (by
        rw [riemannZeta_zero]
        norm_num)
      (by norm_num : (0 : ℂ) ≠ 1) isConnected_compl_one
  rw [codiscreteWithin_iff_locallyFiniteComplementWithin] at hcod
  intro z hz
  obtain ⟨t, ht, hfin⟩ := hcod z hz
  refine ⟨t, ht, hfin.subset fun w hw => ?_⟩
  simp only [mem_inter_iff, mem_diff, mem_preimage, mem_compl_iff,
    mem_singleton_iff, mem_setOf_eq] at hw ⊢
  exact ⟨hw.1, hw.2.1, not_not.mpr hw.2.2⟩

/-- **r499 brick [2c], lemma (1).**  On every compact `K ⊆ ℂ \ {1}`
the zero set of `ζ` is finite. -/
lemma finite_riemannZeta_zeros_of_isCompact {K : Set ℂ}
    (hK : IsCompact K) (hK1 : K ⊆ ({1}ᶜ : Set ℂ)) :
    {z ∈ K | riemannZeta z = 0}.Finite := by
  choose! t ht_nhds ht_fin using fun z (hz : z ∈ K) =>
    riemannZeta_zeros_locallyFinite_compl_one z (hK1 hz)
  obtain ⟨F, hFmem, hFcov⟩ := hK.elim_nhds_subcover t (fun z hz => ht_nhds z hz)
  have hunion :
      {z ∈ K | riemannZeta z = 0} ⊆
        ⋃ z ∈ F, t z ∩ {w : ℂ | w ≠ 1 ∧ riemannZeta w = 0} := by
    intro w hw
    obtain ⟨z, hzF, hwt⟩ := mem_iUnion₂.mp (hFcov hw.1)
    exact mem_iUnion₂.mpr ⟨z, hzF, hwt, hK1 hw.1, hw.2⟩
  exact ((F.finite_toSet).biUnion fun z hz =>
    ht_fin z (hFmem z (Finset.mem_coe.mp hz))).subset hunion

/-- The same finiteness on every compact in `ℂ`: zeros cannot
accumulate at the pole. -/
lemma finite_riemannZeta_zeros_of_isCompact_ne_one {K : Set ℂ}
    (hK : IsCompact K) :
    {z ∈ K | riemannZeta z = 0 ∧ z ≠ 1}.Finite := by
  have hnear := riemannZeta_eventually_ne_zero_nhdsNE_one
  rw [eventually_nhdsWithin_iff, eventually_nhds_iff] at hnear
  obtain ⟨U, hU, hUopen, h1U⟩ := hnear
  have hKc : IsCompact (K \ U) := hK.diff hUopen
  have hK1 : K \ U ⊆ ({1}ᶜ : Set ℂ) := by
    intro z hz h1
    exact hz.2 (h1 ▸ h1U)
  have hfin := finite_riemannZeta_zeros_of_isCompact hKc hK1
  refine hfin.subset ?_
  intro z hz
  refine ⟨⟨hz.1, fun hzU => ?_⟩, hz.2.1⟩
  exact hU z hzU hz.2.2 hz.2.1

/-- Closed strip rectangle `[σ₁, σ₂] × [-T, T]`. -/
def zetaClosedRect (σ₁ σ₂ T : ℝ) : Set ℂ :=
  Icc σ₁ σ₂ ×ℂ Icc (-T) T

lemma isCompact_zetaClosedRect (σ₁ σ₂ T : ℝ) :
    IsCompact (zetaClosedRect σ₁ σ₂ T) :=
  isCompact_Icc.reProdIm isCompact_Icc

/-- **r499 brick [2c], lemma (2).**  Only finitely many zeros of `ζ`
lie in a closed strip rectangle (the [2d] contour support). -/
lemma finite_riemannZeta_zeros_on_closedRect (σ₁ σ₂ T : ℝ) :
    {z ∈ zetaClosedRect σ₁ σ₂ T | riemannZeta z = 0 ∧ z ≠ 1}.Finite :=
  finite_riemannZeta_zeros_of_isCompact_ne_one (isCompact_zetaClosedRect σ₁ σ₂ T)

end ZetaZeroFiniteness

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
