/-
RH/FrequentlySelected.lean -- SEMIDEFINITE + FREQUENTLY SELECTED
(r430, PRIME.RH.SEMIDEFINITE_FREQUENT_SELECTED.01; reviewer
adjudication DCCXCVII).

TWO QUANTIFIER WINS, formalized against the existing extraction
(RH/Elementwise.lean, RH/Selected.lean).  No positivity transport,
no mesh ladder -- the architecture is built for a finite
instantiation per test element.

  (1) SEMIDEFINITE SUFFICES.  Elementwise extraction needs
      `fullRead ≥ 0`, not a strict window margin.  On the R†
      layer the Loewner identity
        `R† ⪰ ½I  ⟺  I − G† ⪰ 0`
      is `Rdagger_ge_half_iff_augmented_posSemidef`
      (RH/DualResolvent.lean; same spectral comparison as A3
      with `PosSemidef` in place of `PosDef`).  The graph-
      resolvent face is `graphResolvent_ge_half_iff`.

  (2) INFINITELY OFTEN SUFFICES (FREQ).  For each `GridElement f`
      the existing `elementwise_finite_stabilization` supplies a
      finite onset `a₀(f)` (comb proved, pole proved, arch the
      existing typed sorry).  `selected_covers` then supplies
      eventual `a_k ≥ a₀(f)` and `m_k ≥ f.meshExp`.  A single
      good index `k` past that onset with `R_k† ⪰ ½I` yields
      `weilForm f = fullRead(W_k, f) ≥ 0`.  Eventually-good
      (`∀ᶠ`) is strictly stronger than frequently-good (`∃ᶠ`);
      the extraction never needed a tail.

HONEST HYPOTHESES of `weil_nonneg_of_frequently_selected`:

  * FREQ of the selected semidefinite cone
    (`∀ K, ∃ k ≥ K, selectedWindowConeSemidef`), equivalently
    `∃ᶠ k in atTop` via `frequently_atTop`;
  * the named Prop `SelectedSemidefImpliesPlainReads`
    (R† ⪰ ½I on window k ⇒ plain `fullRead ≥ 0` for
    mesh-compatible grid elements -- the PSD face of
    `SelectedMasterImpliesPlainReads`; not a theorem here,
    same reason: the bordered read is still typed on the
    rational `VonMangoldtWindow`);
  * onset + mesh coverage PROVED (`selected_covers`);
  * arch-channel stabilization CONSUMED, not re-asserted
    (`elementwise_finite_stabilization`).

THE NEW MINCUT is the named Prop
`frequently_selected_augDualResolvent_ge_half`.  The r397
Prop `selected_augDualResolvent_gt_half` (`∀ᶠ`, strict
`PosDef`) is KEPT as the stronger alternative form; it
implies the new mincut (`frequently_selected_of_eventually_gt_half`).

ALSO IN THIS FILE (sorry-free arithmetic):

  * positive lower density ⇒ frequently (decidable proxy;
    `liminf`-density > 0 is the eventual lower bound);
  * the mean-value trick: a nonnegative integer index `κ`
    with block mean `< 1` yields a zero in the block
    (fallback for an averaged Potapov index).

SORRY CENSUS OF THIS FILE: ZERO.  New openness is named Props
only.  The FREQ extraction consumes the existing classical
sorry `arch_elementwise_stabilization` through
`elementwise_finite_stabilization`; it introduces no new
`sorry`.

Claim boundary: research documentation.  NOT evidence for or
against the Riemann Hypothesis in either direction.  NO RH CLAIM.
-/
import RH.Selected
import RH.GraphResolvent
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.Order.Archimedean.Basic
import Mathlib.Data.Real.Basic

namespace RH

open Filter Matrix Finset
open scoped Topology BigOperators

/-! ## Semidefinite selected-window cone

Same μ-ONB transcription as `selectedWindowCone`; the R†
margin is Loewner (`PosSemidef`) rather than strict
(`PosDef`).  `dualZ ≻ 0` stays -- R† is still an inverse. -/

/-- R† ⪰ ½I for a μ-ONB transcription of one selected real
window. -/
def selectedWindowConeSemidef (k : ℕ) (hk : 0 < k) : Prop :=
  let W := (selectedRealWindow k hk).toPrimeWindow
  ∃ (E : Matrix (Fin W.cap) (Fin W.cap) ℝ) (v : Fin W.cap → ℝ) (γ : ℝ),
    RepresentsLEnsembleReal W W.cap E v γ ∧
    (dualZ E v γ).PosDef ∧
    (augDualResolvent E v γ
      - (1 / 2 : ℝ) •
        (1 : Matrix (Fin W.cap ⊕ Unit) (Fin W.cap ⊕ Unit) ℝ)).PosSemidef

theorem selectedWindowConeSemidef_of_posDef {k : ℕ} (hk : 0 < k)
    (h : selectedWindowCone k hk) :
    selectedWindowConeSemidef k hk := by
  obtain ⟨E, v, γ, hrep, hZ, hR⟩ := h
  exact ⟨E, v, γ, hrep, hZ, hR.posSemidef⟩

/-- **THE NEW MINCUT** (r430): infinitely often the selected
real windows satisfy `R†(W^ℝ_k) ⪰ ½ I`.  Named open kernel.
The r397 Prop `selected_augDualResolvent_gt_half` is the
stronger (`∀ᶠ`, strict) alternative and is kept.  NO RH CLAIM. -/
def frequently_selected_augDualResolvent_ge_half : Prop :=
  ∃ᶠ k in atTop, ∃ hk : 0 < k, selectedWindowConeSemidef k hk

theorem frequently_selected_iff_forall_exists :
    frequently_selected_augDualResolvent_ge_half ↔
      ∀ K, ∃ k, K ≤ k ∧ ∃ hk : 0 < k, selectedWindowConeSemidef k hk :=
  frequently_atTop

/-- The r397 eventually-strict cone implies the r430 mincut. -/
theorem frequently_selected_of_eventually_gt_half
    (h : selected_augDualResolvent_gt_half) :
    frequently_selected_augDualResolvent_ge_half :=
  h.frequently.mono fun _ ⟨hk, hcone⟩ =>
    ⟨hk, selectedWindowConeSemidef_of_posDef hk hcone⟩

/-! ## FREQ plain reads -- the extraction the architecture actually
needs (onset + mesh + one nonnegative `fullRead`). -/

/-- selected-sequence window-local positivity, infinitely often
(plain `fullRead`, the type the extraction consumes). -/
def FrequentlySelectedWindowLocalPositive : Prop :=
  ∃ᶠ k in atTop,
    0 < k ∧ ∀ f : GridElement, f.meshExp ≤ selectedMesh k →
      0 ≤ fullRead (selectedAnchor k) (selectedMesh k) f

/-- **THE FREQ EXTRACTION** (r430; PROVED as a finite
instantiation per element).  Frequently-along-the-sequence
plain-form positivity plus the existing elementwise
stabilization imply Weil-form nonnegativity on every grid
element: coverage is *eventual* (`selected_covers`), positivity
is only *frequent*, and `Eventually.and_frequently` yields one
good covering index.

Honest hypotheses: `hpos` is FREQ of (4) in the r397 list;
onset and mesh coverage are proved; stabilization is the
existing `elementwise_finite_stabilization` (arch sorry
consumed).  NO RH CLAIM. -/
theorem weil_nonneg_of_frequently_plain
    (hpos : FrequentlySelectedWindowLocalPositive) :
    ∀ f : GridElement, 0 ≤ weilForm f := by
  intro f
  obtain ⟨a₀, hstab⟩ := elementwise_finite_stabilization f
  have hcov := selected_covers a₀ f.meshExp
  have hboth := hcov.and_frequently hpos
  obtain ⟨k, ⟨⟨_hkpos, ha, hm⟩, ⟨_, hread⟩⟩⟩ := hboth.exists
  have hrd := hread f hm
  have heq := hstab (selectedAnchor k) ha (selectedMesh k) hm
  rwa [heq] at hrd

/-- Named (never asserted): a semidefinite selected-window cone
`R† ⪰ ½I` implies the plain `fullRead` of every mesh-compatible
grid element at that window.  PSD face of
`SelectedMasterImpliesPlainReads`: the Loewner identity
`Rdagger_ge_half_iff_augmented_posSemidef` gives a nonnegative
augmented form, and the existing `BorderedCompressionBridge`
takes bordered → plain.  Not a theorem in this round: L† is
still typed on the rational `VonMangoldtWindow`. -/
def SelectedSemidefImpliesPlainReads : Prop :=
  ∀ (k : ℕ) (hk : 0 < k),
    selectedWindowConeSemidef k hk →
      ∀ f : GridElement, f.meshExp ≤ selectedMesh k →
        0 ≤ fullRead (selectedAnchor k) (selectedMesh k) f

theorem frequently_plain_of_frequently_selected
    (hbridge : SelectedSemidefImpliesPlainReads)
    (hgood : frequently_selected_augDualResolvent_ge_half) :
    FrequentlySelectedWindowLocalPositive :=
  hgood.mono fun _ ⟨hk, hcone⟩ => ⟨hk, hbridge _ hk hcone⟩

/-- **FREQ of R† ⪰ ½I ⇒ Weil ≥ 0** (r430; PROVED as a function
of the named bridge, never asserting it).  The `hgood`
quantifier is the reviewer's `∀ K, ∃ k ≥ K`.  NO RH CLAIM. -/
theorem weil_nonneg_of_frequently_selected
    (hbridge : SelectedSemidefImpliesPlainReads)
    (hgood : ∀ K, ∃ k, K ≤ k ∧ ∃ hk : 0 < k,
      selectedWindowConeSemidef k hk) :
    ∀ f : GridElement, 0 ≤ weilForm f :=
  weil_nonneg_of_frequently_plain
    (frequently_plain_of_frequently_selected hbridge
      (frequently_atTop.mpr hgood))

/-- Composition to the existing RH-pilot interface:
`∀ f, 0 ≤ weilForm f`.  The spectral/zero side of the explicit
formula is not formalized (RH/Elementwise: not part of
`weilForm`).  Named mincut and named bridge consumed, never
asserted.  NO RH CLAIM. -/
theorem rh_of_frequently_selected
    (hmincut : frequently_selected_augDualResolvent_ge_half)
    (hbridge : SelectedSemidefImpliesPlainReads) :
    ∀ f : GridElement, 0 ≤ weilForm f :=
  weil_nonneg_of_frequently_selected hbridge (frequently_atTop.mp hmincut)

/-! ## Positive lower density ⇒ frequently

A decidable proxy `p` with eventual density `≥ ε > 0` is
frequently true (if `p` died out the density would tend to 0).
`liminf` of the density being positive is this eventual lower
bound.  Applied to the mincut via a proxy that implies the
cone. -/

lemma filter_range_eq_of_not_from {p : ℕ → Prop} [DecidablePred p]
    {N₀ N : ℕ} (hN : N₀ ≤ N) (hp : ∀ n, N₀ ≤ n → ¬ p n) :
    (Finset.range N).filter p = (Finset.range N₀).filter p := by
  ext n
  simp only [mem_filter, mem_range]
  constructor
  · intro h
    refine ⟨?_, h.2⟩
    by_contra hn
    exact hp n (le_of_not_gt hn) h.2
  · intro h
    exact ⟨lt_of_lt_of_le h.1 hN, h.2⟩

/-- Positive lower density of a decidable predicate implies
it holds infinitely often. -/
theorem frequently_of_pos_lower_density {p : ℕ → Prop}
    [DecidablePred p] {ε : ℝ} (hε : 0 < ε)
    (h : ∀ᶠ N : ℕ in atTop,
      ε ≤ (((Finset.range N).filter p).card : ℝ) / (N : ℝ)) :
    ∃ᶠ n in atTop, p n := by
  by_contra hnot
  rw [not_frequently, eventually_atTop] at hnot
  obtain ⟨N₀, hN₀⟩ := hnot
  rw [eventually_atTop] at h
  obtain ⟨N₁, hN₁⟩ := h
  obtain ⟨Nraw, hraw⟩ : ∃ n : ℕ, (N₀ : ℝ) / ε < n := exists_nat_gt _
  let N := N₀ + N₁ + Nraw + 1
  have hN₀le : N₀ ≤ N := by omega
  have hN₁le : N₁ ≤ N := by omega
  have hNrawle : Nraw ≤ N := by omega
  have hNpos : (0 : ℕ) < N := Nat.succ_pos _
  have hdens := hN₁ N hN₁le
  have hcard : (((Finset.range N).filter p).card : ℝ) ≤ (N₀ : ℝ) := by
    have heq := filter_range_eq_of_not_from (N₀ := N₀) (N := N) hN₀le
      (fun n hn => hN₀ n hn)
    have hle : ((Finset.range N).filter p).card ≤ N₀ := by
      rw [heq]
      exact (card_filter_le _ _).trans (by simp [card_range])
    exact_mod_cast hle
  have hNposR : (0 : ℝ) < N := Nat.cast_pos.mpr hNpos
  have hlt : (((Finset.range N).filter p).card : ℝ) / (N : ℝ) < ε := by
    have hle : (((Finset.range N).filter p).card : ℝ) / (N : ℝ) ≤
        (N₀ : ℝ) / (N : ℝ) :=
      div_le_div_of_nonneg_right hcard hNposR.le
    have hN₀lt : (N₀ : ℝ) / (N : ℝ) < ε := by
      have hgt : (N₀ : ℝ) / ε < N :=
        hraw.trans_le (Nat.cast_le.mpr hNrawle)
      have hmul : (N₀ : ℝ) < ε * N := by
        have := (div_lt_iff₀ hε).mp hgt
        rwa [mul_comm] at this
      exact (div_lt_iff₀ hNposR).mpr hmul
    exact lt_of_le_of_lt hle hN₀lt
  exact (not_le_of_gt hlt) hdens

/-- Positive lower density of a decidable proxy that implies
the semidefinite cone yields the r430 mincut. -/
theorem frequently_selected_of_pos_lower_density {p : ℕ → Prop}
    [DecidablePred p]
    (himp : ∀ k, p k → ∃ hk : 0 < k, selectedWindowConeSemidef k hk)
    {ε : ℝ} (hε : 0 < ε)
    (hdens : ∀ᶠ N : ℕ in atTop,
      ε ≤ (((Finset.range N).filter p).card : ℝ) / (N : ℝ)) :
    frequently_selected_augDualResolvent_ge_half :=
  (frequently_of_pos_lower_density hε hdens).mono himp

/-! ## Mean-value trick (integer Potapov index)

If `κ : ℕ → ℕ` (nonnegative) has block mean `< 1` on
`[K, K+N)`, some index in the block is `0`.  Pure arithmetic;
the fallback route "averaged `κ(Θ†) < 1` on blocks ⇒ an
index-0 window exists". -/

/-- Block mean of a nonnegative integer index `< 1` yields a
zero in the block. -/
theorem exists_index_zero_of_block_mean_lt_one (κ : ℕ → ℕ)
    (K N : ℕ) (hN : 0 < N)
    (hmean : (∑ i ∈ Finset.range N, (κ (K + i) : ℝ)) / (N : ℝ) < 1) :
    ∃ i < N, κ (K + i) = 0 := by
  by_contra h
  push Not at h
  have hge : ∀ i ∈ Finset.range N, 1 ≤ κ (K + i) := by
    intro i hi
    exact Nat.one_le_iff_ne_zero.mpr (h i (mem_range.mp hi))
  have hsum : N ≤ ∑ i ∈ Finset.range N, κ (K + i) := by
    calc
      N = ∑ i ∈ Finset.range N, 1 := by simp [sum_const, card_range]
      _ ≤ ∑ i ∈ Finset.range N, κ (K + i) := sum_le_sum hge
  have hsumR : (N : ℝ) ≤ ∑ i ∈ Finset.range N, (κ (K + i) : ℝ) := by
    exact_mod_cast hsum
  have : 1 ≤ (∑ i ∈ Finset.range N, (κ (K + i) : ℝ)) / (N : ℝ) := by
    rw [le_div_iff₀ (Nat.cast_pos.mpr hN)]
    linarith
  linarith

end RH
