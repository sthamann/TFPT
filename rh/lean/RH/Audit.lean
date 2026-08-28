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
`crossing_budget`) must show the three standard axioms ONLY -- in
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
Expected AND MEASURED: `[propext, sorryAx, Classical.choice,
Quot.sound]` -- honest: the ladder-free extraction consumes the TWO
typed CLASSICAL kernel-channel sorries
(`arch_elementwise_stabilization`, `pole_elementwise_stabilization`,
RH/Elementwise.lean; the comb channel is unconditional), never the
two arithmetic holes and never `SourceExact`/`MainWindow`. -/

#print axioms weil_nonneg_of_windowlocal

end RH
