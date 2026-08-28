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

end RH
