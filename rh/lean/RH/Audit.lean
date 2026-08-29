/-
RH/Audit.lean -- THE FINAL AXIOM AUDIT (C1, the reviewer final-domain
contract, goal 4): `#print axioms` on the main theorem chain, run at
every build.

Results are recorded verbatim below after each command (and mirrored
in rh/lean/README.md).  The audit target: the master theorem and its
corollaries may consume ONLY the two canonical arithmetic window
holes (`lstar_canonical`, `terminal_q_canonical` -- both entering as
`sorryAx`), the named dictionary/classical imports (also `sorryAx`,
each individually typed at its statement), and the three standard
mathlib axioms `propext`, `Classical.choice`, `Quot.sound`.  The
sorry-free layer (the reconstruction theorem, the L† equivalences,
`crossing_budget`, the r373 Haynsworth theorems, the r373 dual-resolvent
bridge) must show the three standard axioms ONLY -- in
particular NO `sorryAx`: machine-checkable freedom from hidden holes.

NOTE on `#print axioms` granularity: Lean reports `sorryAx` once,
without naming which sorried declarations feed it.  The reduction "the
master consumes exactly the two canonical holes" is therefore checked
in two machine-verifiable halves: (a) the hypothesis-form
reconstruction `lstar_terminal_implies_master` is sorry-free (audited
below), and (b) the canonical master is literally its application to
`lstar_canonical`/`terminal_canonical` (one-line proof term,
RH/Canonical.lean) -- no other sorried name occurs in the proof.

Claim boundary: research documentation.  NOT evidence for or against
the Riemann Hypothesis in either direction.  NO RH CLAIM.
-/
import RH.Canonical
import RH.DualResolvent
import RH.Haynsworth
import RH.Elementwise
import RH.PivotCoordinate
import RH.FlankEntry
import RH.Selected
import RH.OneDefect
import RH.GraphResolvent

namespace RH

/-! ## (a) The sorry-free load-bearing layer
Expected AND MEASURED (C1 record, verbatim build output):
`[propext, Classical.choice, Quot.sound]` -- NO `sorryAx`:
  * `RH.lstar_terminal_implies_master`
  * `RH.augmentedSubordination_iff_lstar_and_terminal`
  * `RH.crossing_budget`
  * `RH.canonical_budget_pos` -/

-- the r305 reconstruction theorem (hypothesis form)
#print axioms lstar_terminal_implies_master
-- the r332 L† equivalence (hypothesis form)
#print axioms augmentedSubordination_iff_lstar_and_terminal
-- the Jacobi/Sylvester crossing budget -- PROVED since C1
#print axioms crossing_budget
-- the canonical budget positivity (the proved half of the old
-- terminal statement)
#print axioms canonical_budget_pos

/-! ## (b) The two canonical arithmetic holes (the final domain)
Expected AND MEASURED: `[propext, sorryAx, Classical.choice,
Quot.sound]` -- each IS a typed `sorry`. -/

#print axioms lstar_canonical
#print axioms terminal_q_canonical

/-! ## (c) The main theorem chain on the final domain
Expected AND MEASURED: `[propext, sorryAx, Classical.choice,
Quot.sound]`, the `sorryAx` entering ONLY through (b) -- see the
granularity note in the header: the proof terms apply the sorry-free
layer (a) to the two holes (b) and to nothing else (the pair closure
additionally consumes the ONE named dictionary import
`pair_terminal_dictionary`, typed transcription-blocked). -/

#print axioms augmented_prefix_positive
#print axioms augmentedSubordination_main
#print axioms terminal_crossratio_main
#print axioms prefix_chain_positive_main
#print axioms free_window_positivity
#print axioms pair_closes_main

/-! ## (d) The Level-C extraction (classical kernel imports)
Expected AND MEASURED (r376): `[propext, sorryAx, Classical.choice,
Quot.sound]` -- honest: the ladder-free extraction consumes the ONE
remaining typed CLASSICAL kernel-channel sorry
(`arch_elementwise_stabilization`; the comb channel is unconditional
and the pole channel is PROVED as the native-mesh second-difference
pairing of `polePotential`), never the two arithmetic holes and never
`SourceExact`/`MainWindow`. -/

#print axioms weil_nonneg_of_windowlocal
#print axioms pole_elementwise_stabilization

/-! ## (e) The r362/r373 dual-resolvent layer (reviewer priority 2)
Sorry-free finite algebra AND the window bridge, expected AND MEASURED:
`[propext, Classical.choice, Quot.sound]` -- NO `sorryAx`:
  * `RH.posDef_one_sub_iff_dualResolvent_gt_half` (r356-A / A2)
  * `RH.posDef_one_sub_borderedGram_iff_augDualResolvent` (A3)
  * `RH.posDef_one_sub_borderedGram_iff_qDagger` (A5)
  * `RH.augDualResolvent_fromBlocks` (A4 Sherman–Morrison)
  * `RH.augDualResolvent_gt_smul_implies_dualResolvent` (A7-min)
  * `RH.augmentedSubordination_iff_dualResolvent` -- r373: the μ-ONB
    Gram transcription `RepresentsLEnsemble` is a real Prop (whitening
    equations, R319); the iff is congruence of `A_cap` onto `I − G†`
    then A3.  Census 8 → 7.  NO RH CLAIM. -/

#print axioms posDef_one_sub_iff_dualResolvent_gt_half
#print axioms posDef_one_sub_borderedGram_iff_augDualResolvent
#print axioms posDef_one_sub_borderedGram_iff_qDagger
#print axioms augDualResolvent_fromBlocks
#print axioms augDualResolvent_gt_smul_implies_dualResolvent
#print axioms augmentedSubordination_iff_dualResolvent

/-! ## (f) The r373 Haynsworth layer (reviewer goal 1)
Expected AND MEASURED: `[propext, Classical.choice, Quot.sound]` --
NO `sorryAx`.  Finite real algebra; does not assert r367 P1/P2 on
any window. -/

#print axioms haynsworth_two_rank
#print axioms haynsworth_mixed
#print axioms haynsworth_sigNeg_₁₁
#print axioms haynsworth_sigNeg_₂₂

/-! ## (g) r373/r376 kernel objects
Expected AND MEASURED: `[propext, Classical.choice, Quot.sound]` --
NO `sorryAx` on the transcribed closed forms, and NO `sorryAx` on
`pole_elementwise_stabilization` (r376).  The arch stabilization
sorry remains in (d). -/

#print axioms polePotential_even
#print axioms polePotential_eq_cosh
#print axioms weilArchKernel_even

/-! ## (h) Historical opacity bridge (Alt-Last, not a load-bearing hole)
`mainWindow_iff_builtFromPrimeSource` (RH/Source.lean) is the r310/r320
opacity bridge.  Since C1 the load-bearing chain quantifies over
`CanonicalWindow`, not `MainWindow`.  The bridge is consumed only by
`mainWindow_explicit_bridge` (itself a corollary of the bridge).  It
is an Alt-Last of the historical interface, not an arithmetic or
classical hole; it is NOT deleted (the r273 opaque marker and the
U1–U3 guards still refer to it).  Expected: `[propext, sorryAx,
Classical.choice, Quot.sound]`. -/

#print axioms mainWindow_iff_builtFromPrimeSource

/-! ## (i) The r380 pivot-coordinate layer (DCCXXXVII finite algebra)
Expected AND MEASURED: `[propext, Classical.choice, Quot.sound]` --
NO `sorryAx` on the proved faces.  Named Props
(`ComplementaryDualHankelInertia`, `DPPIdentity`,
`SignedBorodinComplement`, `K2EqHankelRatio`, `P1EqCapInertia`,
`P2EqPostcapAlternation`) are hypotheses, not holes.  Census
unchanged at 5.  Does not assert (P1)/(P2) on any window.
NO RH CLAIM. -/

#print axioms rankOne_inertia_antitone
#print axioms adaptive_band_from_entry
#print axioms sigNeg_full_hankel_eq_sigNeg_weights
#print axioms postcap_pivot_ratio_eq_h_form
#print axioms indNeg_hankel_eq_neg_pivot_count
#print axioms p1_p2_iff_cap_posDef

/-! ## (j) The r384 flank-entry layer
Expected AND MEASURED: `[propext, Classical.choice, Quot.sound]` --
NO `sorryAx` on the proved faces.  Named Props (`FlankEntryPrefix`,
`ChristoffelPivotBound`) are hypotheses, not holes.  Census
unchanged at 5.  Does not assert (P1)/(P2) or L* on any window.
NO RH CLAIM. -/

#print axioms pair_energy_identity
#print axioms h_zero
#print axioms h_one
#print axioms h0_pos_of_mass
#print axioms h1_pos_of_pairEnergy
#print axioms three_atom_mass_pos
#print axioms threeAtom_flank_pivots
#print axioms clusterRun3_H3
#print axioms fiveAtom_energy
#print axioms christoffel_bound_k0
#print axioms christoffel_bound_k1
#print axioms threeAtom_christoffel_k1
#print axioms indNeg_entry_of_flank
#print axioms adaptive_band_from_flank_entry

/-! ## (k) The r397 exact selected domain
Expected AND MEASURED: the sequence identities and cofinality
theorems are sorry-free
(`[propext, Classical.choice, Quot.sound]` -- NO `sorryAx`).
`weil_nonneg_of_selected_windows` consumes the existing classical
arch-channel sorry through `elementwise_finite_stabilization`, so
`[propext, sorryAx, Classical.choice, Quot.sound]` -- the same
`sorryAx` as `weil_nonneg_of_windowlocal`, no new hole.
Named Props (`selected_augDualResolvent_gt_half`,
`SelectedMasterImpliesPlainReads`, `ExactArchAgreesWithArchRead`)
are hypotheses, not holes.  Census of `sorry` declarations
unchanged at 5; the two C1 holes are degraded to the alternative
route.  NO RH CLAIM. -/

#print axioms selectedDelta_eq
#print axioms selectedAnchor_tendsto
#print axioms selectedRoot_tendsto
#print axioms selectedMesh_tendsto
#print axioms selectedDelta_tendsto_zero
#print axioms selected_covers
#print axioms realCanonicalWindow_B_pos
#print axioms ExactFold_B
#print axioms weil_nonneg_of_selected_windows
#print axioms weil_nonneg_of_selected_master

/-! ## (l) The r406 one-defect finite algebra
Expected AND MEASURED: `[propext, Classical.choice, Quot.sound]` --
NO `sorryAx`.  Pure finite real matrix algebra (Schur rank-1
criterion, Woodbury identity, Cauchy--Schwarz Gram lift).
Independent of the source-side R404/R405 probes.  Census of
`sorry` declarations unchanged at 5.  NO RH CLAIM. -/

#print axioms indNeg_eq_zero_of_posDef
#print axioms indNeg_sub_rankOne_le_one
#print axioms posDef_sub_rankOne_iff
#print axioms woodbury_inv
#print axioms oneDefect_update_posDef_iff
#print axioms posDef_of_contractive_lift
#print axioms cMin_normSq
#print axioms posDef_gram_sub_rankOne_iff

/-! ## (m) The r412 graph-resolvent finite algebra
Expected AND MEASURED: `[propext, Classical.choice, Quot.sound]` --
NO `sorryAx`.  Pure finite real matrix algebra (spectral dictionary
of `R = C(I+C)⁻¹`, Möbius inertia, energy-split contraction).
Named Prop `GraphResolventIsLEnsembleInv` is a hypothesis, not a
hole.  Census of `sorry` declarations unchanged at 5.  NO RH CLAIM. -/

#print axioms graphResolvent_eq_one_sub_inv
#print axioms graphResolvent_eq_dualResolvent_inv
#print axioms graphResolvent_lt_one
#print axioms graphResolvent_posDef
#print axioms graphResolvent_sub_half_eq
#print axioms indNeg_one_add_inv_congruence
#print axioms indNeg_inv_congruence
#print axioms indNeg_graphResolvent_sub_half
#print axioms indNeg_mobius
#print axioms posDef_one_sub_inv_iff
#print axioms posSemidef_one_sub_inv_iff
#print axioms graphResolvent_gt_half_iff
#print axioms contractive_iff_gram_le_one
#print axioms energy_split_contractive
#print axioms energy_split_at_most_one
#print axioms p1_coord_graphResolvent
#print axioms augDualResolvent_gt_half_of_C_gt_one

end RH
