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

end RH
