/-
RH/InnerBridges.lean -- FINITE INNER-BRIDGE ALGEBRA
(r464, PRIME.RDAGGER.INNER_BRIDGES.01).

The selected A-cap implication has two logically separate parts:

  (1) CHANNEL TRANSCRIPTION: every mesh-compatible GridElement read is
      the quadratic value of a direction in the augmented cap space;
  (2) FINITE PSD: A_cap.PosSemidef makes that quadratic value nonnegative.

Part (2) is proved below.  Part (1) is the one named remainder
`SelectedReadQuadraticRepresentation`.  It includes the real
quantifier issue: `GridElement.steps` is not bounded by the cap, so
the old bridge cannot be obtained from a degree-<cap polynomial
identification without an additional high-degree reduction theorem.

r470 seals the exact finite obstruction to part (1): the original
quantifier admits `steps > cap` (proved), and a mesh-and-onset-
compatible native Hessian at k=5 is indefinite of rank larger
than `cap+1` while `A_cap` is PD (numerical seal).  No new `sorry`.

Claim boundary: research documentation.  NO RH CLAIM.
-/
import RH.Selected

namespace RH

open Matrix

/-- The augmented quadratic read of a selected production window. -/
noncomputable def selectedReadQuadratic
    (k : ℕ) (hk : 0 < k)
    (z : Fin (selectedRealWindow k hk).cap ⊕ Unit → ℝ) : ℝ :=
  z ⬝ᵥ ((selectedRealWindow k hk).toPrimeWindow.A
    (selectedRealWindow k hk).cap *ᵥ z)

/-- **The exact remaining channel transcription.**

For every selected window and every mesh-compatible native test
element, construct the augmented cap-space direction whose quadratic
value is the concrete `fullRead`.  This is where the Fourier/Chebyshev
coordinates, whitening, border coordinate, and pole correction must
be identified.  It is a named Prop, not asserted.

The proposition deliberately retains the original unbounded `steps`
quantifier; no hidden `steps ≤ cap` weakening is introduced. -/
def SelectedReadQuadraticRepresentation : Prop :=
  ∀ (k : ℕ) (hk : 0 < k) (f : GridElement),
    f.meshExp ≤ selectedMesh k →
      ∃ z : Fin (selectedRealWindow k hk).cap ⊕ Unit → ℝ,
        fullRead (selectedAnchor k) (selectedMesh k) f =
          selectedReadQuadratic k hk z

/-- The original r434 target, retained byte-for-byte in mathematical
content: selected A-cap PSD implies every mesh-compatible plain read
is nonnegative. -/
def SelectedACapPsdImpliesPlainReads : Prop :=
  ∀ (k : ℕ) (hk : 0 < k),
    ((selectedRealWindow k hk).toPrimeWindow.A
      (selectedRealWindow k hk).cap).PosSemidef →
      ∀ f : GridElement, f.meshExp ≤ selectedMesh k →
        0 ≤ fullRead (selectedAnchor k) (selectedMesh k) f

/-- **FINITE PSD HALF CLOSED (r464).**  Once the read direction is
transcribed, A-cap semidefiniteness proves the original plain-read
bridge by the defining quadratic inequality. -/
theorem selectedACapPsdImpliesPlainReads_of_representation
    (hrep : SelectedReadQuadraticRepresentation) :
    SelectedACapPsdImpliesPlainReads := by
  intro k hk hA f hm
  obtain ⟨z, hz⟩ := hrep k hk f hm
  rw [hz]
  exact hA.dotProduct_mulVec_nonneg z

/-! ## r470: exact finite obstruction to the channel map

The original representation quantifies over every mesh-compatible
`GridElement`, with `steps` unbounded.  The selected cap is
`(S+1)/2` with `S = selectedMesh k`, so a degree-`< cap`
identification cannot cover the original quantifier.  The
arithmetic witness is proved below.

A stronger finite obstruction is numerical (probe
`quadrep_probe.py`): even after the onset guard
`elementAnchor f ≤ selectedAnchor k`, the native Hessian of
`fullRead` at `k=5` is indefinite of rank `24 > cap+1 = 11`,
while `A_cap` is positive definite.  A negative native read
cannot equal a quadratic value of a PSD matrix.  The lemma
`quadraticRepresentation_refuted_of_negative_read` records
that finite implication; the signed witness itself is the
named Prop `SelectedOnsetCompatibleNegativeRead` (not a
`sorry`).  NO RH CLAIM. -/

lemma selectedRoot_five : selectedRoot 5 = 2 := by
  unfold selectedRoot
  rw [Nat.floor_eq_iff (Real.sqrt_nonneg _)]
  constructor
  · rw [Real.le_sqrt (Nat.cast_nonneg 2) (Nat.cast_nonneg 5)]
    norm_num
  · have h3 : ((2 : ℕ) : ℝ) + 1 = 3 := by norm_num
    rw [h3, Real.sqrt_lt (Nat.cast_nonneg 5) (by norm_num : (0 : ℝ) ≤ 3)]
    norm_num

lemma selectedMesh_five : selectedMesh 5 = 19 := by
  unfold selectedMesh
  rw [selectedRoot_five]
  native_decide

lemma selectedWindow_S_eq_mesh (k : ℕ) (hk : 0 < k) :
    (selectedRealWindow k hk).toPrimeWindow.S = selectedMesh k :=
  rfl

lemma selected_cap_five (hk : 0 < 5) :
    (selectedRealWindow 5 hk).toPrimeWindow.cap = 10 := by
  unfold PrimeWindow.cap
  rw [selectedWindow_S_eq_mesh 5 hk, selectedMesh_five]

/-- Concrete mesh-compatible native element with `steps > cap`
at the r464 k=5 arithmetic witness. -/
def selectedCapStepsWitness : GridElement where
  steps := 11
  meshExp := 19
  x := fun _ => 0

/-- **r470, proved:** the original mesh-only guard admits a
native element whose step count strictly exceeds the selected
cap.  This is the quantifier obstruction to a degree-`< cap`
channel, not a statement about infinitely many `k`. -/
theorem exists_mesh_compatible_steps_gt_cap :
    ∃ (k : ℕ) (hk : 0 < k) (f : GridElement),
      f.meshExp ≤ selectedMesh k ∧
      (selectedRealWindow k hk).toPrimeWindow.cap < f.steps := by
  refine ⟨5, Nat.succ_pos 4, selectedCapStepsWitness, ?_, ?_⟩
  · simp [selectedCapStepsWitness, selectedMesh_five]
  · simp [selectedCapStepsWitness, selected_cap_five (Nat.succ_pos 4)]

/-- **r470 named obstruction.**  A mesh-and-onset-compatible
native element with strictly negative `fullRead`.  On a PSD
cap this refutes `SelectedReadQuadraticRepresentation` by
`quadraticRepresentation_refuted_of_negative_read`.  Sealed
numerically at `k=5` (`meshExp=3`, `steps=24`); not a `sorry`. -/
def SelectedOnsetCompatibleNegativeRead : Prop :=
  ∃ (k : ℕ) (f : GridElement),
    0 < k ∧
    f.meshExp ≤ selectedMesh k ∧
    f.elementAnchor ≤ selectedAnchor k ∧
    fullRead (selectedAnchor k) (selectedMesh k) f < 0

/-- Finite implication: a negative mesh-compatible read on a
PSD selected cap refutes the r464 channel representation.
No new `sorry`.  NO RH CLAIM. -/
theorem quadraticRepresentation_refuted_of_negative_read
    {k : ℕ} (hk : 0 < k) (f : GridElement)
    (hm : f.meshExp ≤ selectedMesh k)
    (hneg : fullRead (selectedAnchor k) (selectedMesh k) f < 0)
    (hA : ((selectedRealWindow k hk).toPrimeWindow.A
      (selectedRealWindow k hk).cap).PosSemidef) :
    ¬ SelectedReadQuadraticRepresentation := by
  intro hrep
  obtain ⟨z, hz⟩ := hrep k hk f hm
  have hnn : 0 ≤ fullRead (selectedAnchor k) (selectedMesh k) f := by
    rw [hz]
    exact hA.dotProduct_mulVec_nonneg z
  exact (not_le_of_gt hneg) hnn

end RH
