/**
 * Live update feed for /prime-front ("The Prime Front").
 *
 * HOW TO POST A NEW ENTRY after a completed agent run:
 *   1. Open this file.
 *   2. Prepend a new object to `PRIME_FRONT_UPDATES` (newest first).
 *   3. Fill every field below; keep `summary` to 1–2 sentences.
 *   4. Rebuild / redeploy the website — no other files required for a feed post.
 *
 * Do not invent numbers. Copy verdicts and anchors from experiments/next.txt
 * (or the promoted verification/vN_*.py). Firewall: no "RH evidence" language.
 */

export type PrimeFrontBadge = "sandbox" | "machine-verified" | "running";

export type PrimeFrontVerdict =
  | "EXACT"
  | "PARTIAL"
  | "DEFLATED"
  | "DEAD"
  | "KILLED-AS-NAIVE"
  | "MIXED"
  | "BENCHMARK"
  | "REPAIRED"
  | "HARDENED"
  | "FOUNDED"
  | "PROMOTED"
  | "MACHINE-VERIFIED"
  | "CLOSED"
  | "TERRAIN-MAPPED"
  | "SCALABLE-AND-SIGNED"
  | "BRIDGE-FOUND"
  | "FUNCTOR-WIRED"
  | "RUNNING"
  | "TYPED"
  | "DENSITIES-CLOSED"
  | "OUTSIDE-KOHNEN-SCOPE"
  | "CENTRAL-VALUES-WIRED"
  | "IMPL-OPEN"
  | "UV-NO-SIGNAL"
  | "SCALE-TORSOR"
  | "FAMILY-AGGREGATED"
  | "MINIATURE-RUNS"
  | "INFINITE-CARRIER"
  | "MARGINAL-WEIGHT"
  | "PACKED"
  | "TRANSFORM-REQUIRED"
  | "HALF-STABLE"
  | "STATUS"
  | "DIRAC-SQRT-EXACT"
  | "PI-FRONT-CLOSED"
  | "DELETION-UNIVERSAL"
  | "LINEAR-PLUS-COMBINATION"
  | "COMPLEMENTARY-PAIR"
  | "DOORS-FURNISHED"
  | "RECIPE-UNIVERSAL-ON-BATTERY"
  | "LEMMA-CLASSICAL-SHAPED"
  | "WINDOW-PROVED"
  | "LEDGER-CLOSES"
  | "RESERVE-PARTIAL"
  | "AVOIDANCE-FAILS"
  | "ARCH-INTERNAL"
  | "INVARIANT-NULL"
  | "LIFT-WORKS-UNANCHORED"
  | "LEMMA-CLOSES-LAMBDA"
  | "LEMMA-FULLY-CLOSED"
  | "DICTIONARY-EXACT-CORE"
  | "TIGHT-SET-PARAMETRIZED"
  | "CROSSING-MAPPED"
  | "CORE-DISSECTED"
  | "BAND-PARTIAL"
  | "T-SKELETON"
  | "MIXED"
  | "BLIND-100"
  | "T-CONTINUUM-NUMERIC"
  | "RELAY-CONFIRMED"
  | "ALIGNMENT-ONLY"
  | "LAW-CONFIRMED-MECHANISM-OPEN"
  | "DECAY-LAW-FOUND"
  | "REMAINDER-CLOSES-ZONES"
  | "CROWDING-TRENDS"
  | "MECHANISM-IDENTIFIED"
  | "INSTRUMENT-IMPROVED"
  | "CHAIN-PARTIAL"
  | "ONE-OF-TWO"
  | "DENSITY-MAPPED"
  | "SCALAR-TRACTABLE"
  | "EPSILON-IDENTITY"
  | "BOUNDARY-CERTIFIED"
  | "MARGIN-PROPAGATES"
  | "CROSSING-CONFIRMED"
  | "SCALING-PARTIAL"
  | "SUBSTANCE-CONFIRMED"
  | "WALL-DISSOLVES"
  | "TRANSPORT-BLOCKED"
  | "RICCATI-PARTIAL"
  | "THEOREM-SHAPED"
  | "TWO-OF-THREE";

export type PrimeFrontUpdate = {
  /** ISO date (YYYY-MM-DD) of the agent run. */
  date: string;
  /**
   * Diary part number (Teil N).
   * Use `0` for meta / promotion announcements that span multiple diary parts
   * (string anchors are not supported by this field).
   */
  part: number;
  title: string;
  verdict: PrimeFrontVerdict;
  /** 1–2 sentences; plain English; scientifically honest. */
  summary: string;
  badge: PrimeFrontBadge;
  /** Optional probe / verification script basename. */
  script?: string;
};

/**
 * Newest first. Future agent runs: prepend here.
 */
export const PRIME_FRONT_UPDATES: readonly PrimeFrontUpdate[] = [
  {
    date: "2026-07-27",
    part: 118,
    title:
      "Two of three lemmas stand — saturation is an identity here, the Schur floor is rescued by the CBS shift onto the oscillation Gram, and the corner lemma corrects itself; one genuinely new analytic statement remains",
    verdict: "TWO-OF-THREE",
    summary:
      "T118 (TWO-OF-THREE, 36/36) — contract SYMBOL.LEMMAS: measure the three missing lemmata of the T117 theorem candidate one by one against their classical addresses. (Q1) Lemma B, attempt one — REFUTED, WITH A REASON: the exact two-grid symbol of the oscillation Schur complement S is the WEIGHTED HARMONIC MEAN of the aliasing pair, 1/σ_S = cos²(θ/2)/f(θ+π) + sin²(θ/2)/f(θ) (derived, verified on exactly positive model symbols). A harmonic mean needs BOTH entries positive, and the window's comb dips make inf σ_S = 0 on 14/14 (zone, level) rows — every pointwise route through σ_S is VACUOUS on the real symbol, not merely weak (a perturbative split recovers at most 2.6e-5 of λ_min(S); the comb is non-perturbative). The classical reason is computed, not asserted: λ_min at finite M is a Fejér–Riesz MOMENT problem, not a pointwise infimum — the ground state of S CONCENTRATES on the negative set (3.3–43× the uniform share) and positivity is a cancellation of relative size 10²–10⁴; the one-cell-averaged symbol is positive on 14/14 windows (a heuristic, declared as such). GROWTH: on an FFT-only lever of up to 128× in D the LOG model beats the POWER model (b = 0.44–0.93 per log(1/D)) — a logarithm, as the log-singular kernel demands (Kac–Murdock–Szegő 1953 / Widom 1974); and the pole mass in the oscillation space is ≤ 2.1e-5 (an identity via the DTFT factor |sin(θ/2)|) — which is why the finite-section effect is absent there. (Q2) Lemma A, the corner — THE ADDRESS DOES NOT TRANSFER: the edge-exponent estimator is calibrated (model class |2sin(θ/2)|^{2s} exact, p = s; the KMS 1953 closed-form inverse as the s = 0 control), but on the real symbol the corner exponent DRIFTS (+0.238 per halving of D, 3–7.4σ) — log(1/|θ|) is the boundary of all powers, so there is no asymptotic (s, κ) to cite. (The old T108/T109 24-cell exponent was DOUBLY mis-estimated; its apparent stability was the cancellation of two errors.) The (H3) correction: the single-cell form is mis-dimensioned — the concentration factor grows ~h^0.99, the slope mass sits on a fixed FRACTION of cells; the T117 chain is REPAIRED at its last step: the rate lives on ‖y‖² (stage-2 bound D^1.751 against ε ~ D^1.770 — lossless; the single-cell stage costs D^2.74). (Q3) Lemma C — AND THE RESCUE OF LEMMA B: saturation is not an assumption for a nested PWC pair but an IDENTITY, ε_c − ε_f = yᵀSy (rel 4.9e-11), so β is COMPUTED: β ∈ [0.252, 0.336] over 20 triples, monotone in the favourable direction — what remains of Lemma C is the one hypothesis Bank–Smith 1993 themselves assume. Their CBS constant γ (Yserentant 1986) is the lever that rescues Lemma B: the identity λ_min(L_z⁻¹SL_z⁻ᵀ) = 1−γ² (rel 2.1e-15) moves the symbol question off S and onto the oscillation Gram with σ_z = sin²(θ/2)f(θ) + cos²(θ/2)f(θ+π) — the ARITHMETIC mean of the same aliasing pair, the low-frequency negativity suppressed QUADRATICALLY by sin²(θ/2) instead of poisoning a harmonic mean: inf σ_z > 0 on 7/14 windows (systematically: every window with α ≤ 1.30 and h ≥ 200) with up to 52% sharpness; on ALL 14 windows λ_min(S) ≥ (1−γ²)λ_min(T_h(σ_z)) at 32–76% — Lemma B is reduced to the positivity of a plain Toeplitz section of an explicit symbol; and on the long FFT lever (M up to 15,680, no matrix ever formed) the certified floor RISES logarithmically and CROSSES ZERO on 3/5 zones — the failing windows were UNDER-RESOLVED, not obstructed. (Q4) Theorem V2: five unconditional links (identities + isometry + Grenander–Szegő + certified Choleskys) + (H1) γ-uniformity (measured 1−γ² ≥ 0.181) + (H2) ‖y‖² ≥ c·D^1.75 — THE one genuinely new analytic statement the theorem needs, not importable from the KMS/Widom power-law corner theory because this symbol is log-singular with atoms. Per-lemma missing lists: L-B two points (both ARITHMETIC: inf σ_z > 0 for D < D₀(α) — comb depth against Nyquist-band mass, an atom-counting question — and narrow-dip positivity of finite sections), L-C one point (γ-uniformity, (H1)), L-A two points ((H3′) mean-square form; the Fisher–Hartwig extension of the corner theory to log-singular symbols with atoms). T119 (OSCILLATION.MASS, oscillation_mass_probe.py) is running: the arithmetic half of L-B (the D₀(α) formula, deep verification), the (H2) section, the (H3′) mean-square form, a theorem candidate V3. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "symbol_lemmas_probe.py",
  },
  {
    date: "2026-07-27",
    part: 117,
    title:
      "Theorem-shaped — the one inequality is now an identity plus a certified lower-bound chain that loses no power of D; the jumps have closed forms, and T116's factor-120 drops were a sweep artifact",
    verdict: "THEOREM-SHAPED",
    summary:
      "T117 (THEOREM-SHAPED, 23/23) — contract EPSILON.THEOREM: make the one remaining inequality proof-shaped. (O1) The Galerkin identity, exact — four identities on 6 windows: t̃ is the cell functional of the D-INDEPENDENT function 2sinh(x/2) (independent quadrature, rel 5.2e-16); the family is EXACTLY NESTED (T_c = PᵀT_fP, t̃_c = Pᵀt̃_f, rel 2e-14) — so ε is a Galerkin best-approximation error of ONE bilinear form and its monotonicity under refinement is a THEOREM, not a fit; ε = min_v [1 − 2t̃ᵀv + vᵀTv] with the two-level Pythagoras ε_c − ε_f = ‖u_f − Pu_c‖²_T (rel 6.2e-11); and the dual-norm residual form gives ε > 0 ⟺ t̃ ∉ T(V_c) — positivity is a NON-MEMBERSHIP statement. (O2) The lower bound, honestly (the direction trap handled in the open: Céa/Aubin–Nitsche/Bramble–Hilbert give UPPER bounds; the missing lower-bound ingredient has a classical name — the saturation assumption, Bank–Smith 1993 / Dörfler–Nochetto 2002): the psd-minorant route is DEAD — criticality forces relative sharpness 8.3e-4 in the u-direction, which is circular; the TWO-LEVEL CHAIN survives: ε ≥ ε_c − ε_f = yᵀSy ≥ λ_min(S)‖y‖² ≥ λ_min(S)·D_f³·max|slope|²/2 (S the T-Schur complement onto the oscillation space; the last step Payne–Weinberger 1960 in its lower form) — CERTIFIED on 19/19 pairs, bound/ε ∈ [0.111, 0.185], at rate θ' = 1.74 against θ = 1.76: NO POWER OF D LOST (λ_min(S) GROWS under refinement, +0.20 per halving of D — the chain costs a constant, spread 1.195). The one-cell version costs 0.93 powers — the classical price of replacing a sum by its largest term. (O3) The jumps, exactly — and a correction to T116: adding one cell per end is a Levinson bordering, ε(h+1) = ε(h) − r₀²/s₀ (rel 3.7e-12); a prime-power entry is a corner update of rank EXACTLY k₀+1 whose Woodbury closed form reproduces all 23 measured entries (rel 2.1e-11). The factor-120 'jumps' were a SWEEP ARTIFACT: all 23 entries move ε UP (factors 1.0003–1.0506, share −0.3% of the log-drop); the falls are the accumulated smooth bordering product between sweep points. The jump side-condition on the ansatz is REMOVED; what remains is a smooth α^(−6.07±0.03) decay (fit, short lever). (O4) The theorem candidate is written out (H1 T ≻ 0 certified; H2 λ_min(S) ≥ κ certified per window; H3 a non-vanishing discrete slope on one cell; conclusion the four-line chain; line-by-line attribution exact/Cholesky/classical/fit). The honest missing list: three named analytic lemmas about ONE symbol, each with a classical address — (1) corner asymptotics of T⁻¹ (Kac–Murdock–Szegő 1953 / Widom 1974), (2) a lower bound for λ_min(S) at ONE level (a symbol question at the Nyquist frequency π/D), (3) the saturation constant (measured [0.675, 0.744], classical, costs no rate); two of the three are CONSTANTS, not rates. None of this touches T116's [SUPPORT] comb obstruction (a separate state-compression question). T118 (SYMBOL.LEMMAS, symbol_lemmas_probe.py) is running: the three lemmas against their classical addresses (the S-symbol at Nyquist, KMS/Widom corners, Bank–Smith saturation). Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "epsilon_theorem_probe.py",
  },
  {
    date: "2026-07-27",
    part: 116,
    title:
      "Riccati partial — the induction step IS a boundary process: the pole rides exactly in a 12×12 state and the march runs 903× deeper at flat cost; what refuses compression is the prime comb itself",
    verdict: "RICCATI-PARTIAL",
    summary:
      "T116 (RICCATI-PARTIAL, 33/33) — contract BOUNDARY.FORMULATION: never represent the old interior at all — run the induction as a pure boundary process. (N1) The boundary recursion stands: the state (Σ, r, σ) is the partial minimisation of the pole-free form over the interior (Haynsworth 1968); the identity tᵀT⁻¹t = rᵀΣ⁻¹r + σ holds on 6 zones at rel 6.7e-10; the Haynsworth double-Albert equivalence certifies 6/6 with sharp negative controls; and the Riccati chaining is exact (synthetic 1.4e-16, bit-exact on real windows). POLE BOOKKEEPING SOLVED, FOR FREE: the pole is genuinely global (38.6% boundary value at half depth vs 0.2% for the arch part), but (r, σ) carries it EXACTLY via bordered elimination — no truncation anywhere. TAIL BOOKKEEPING BREAKS, unexpectedly: the archimedean part decays cleanly (exp(−1.73ω)), but the FULL symbol does not decay at all (≥ 93% of its maximum beyond any width) — 100% of the off-band mass sits in the stripes of the two PRIME COMBS, carried by the reflection comb: every incoming cell couples via a Hankel term to the interior cell at α − log n_j, for every prime power √n < n_j < n. The truncation budget is 1.3e2…1.1e6 × the margin ε it must respect and grows like n^2.17 (μ_j ~ n^−1/4 vs ε ~ n^−1.8); 11/18 real truncations break outright; a comb-aware state keeps 83.6% of the window and still only reaches rel 20; the ceiling of the whole idea is a Θ(log²n) gain. (N2) The deep boundary run: 169,236 prepends to h = 1,354,088 cells = 903× the old cap, at FLAT cost (76 µs/step, first/last-decile ratio 1.00), with the entire state 12×12+12+1 numbers. Declared honestly: a cost-geometry demonstration, NOT a Weil certificate (at D = 1.9e-6 the band contains no atom; the stopper is ε < 0 at h = 1.35e6 — the band model loses positivity, not the machine); atom reach stays n = 173. (N3) The one inequality, sharpened: ε ~ 8.34·D^1.790·α^(−6.04) (jackknife errors ≤ 0.041, rms ≤ 0.097; φ from a fixed-D sweep ± 0.08); ε is NOT smooth in α (factor-120 jumps at prime-power entries — a smooth ansatz is wrong); the continuum form is EXACTLY CRITICAL (t̃ᵀT⁻¹t̃ up to 0.999975745; Szegő–Kolmogorov: the pole is representable in the limit); all measurements sit ≥ 1.1e5 above the Cholesky floor. THE CLASSICAL HIT: the Galerkin reading — t̃ represents the D-independent function 2sinh(x/2), ε is the error of a piecewise-constant Galerkin method, and Aubin–Nitsche duality predicts D^{2s}: measured θ = 1.79 ⟹ s ≈ 0.90 (log singularity + window edges). TARGET INEQUALITY formulated: ε(α,D) ≥ c₀·D^θ·α^φ across the jumps. (N4) The remaining list: WON — exact pole, exact Riccati chaining, flat costs; BROKEN — [SUPPORT]: the comb, no sparse faithful state, the stripes drift with fill-in; [P1] the one inequality, now with φ and a jump side-condition; [P2] transport DEVALUED (a single march needs no regrid); the basis is uncritical. T117 (EPSILON.THEOREM, epsilon_theorem_probe.py) is running: the exact Galerkin identity, an Aubin–Nitsche lower-bound chain with constants, exact jump formulas (Sherman–Morrison), a theorem candidate. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "boundary_formulation_probe.py",
  },
  {
    date: "2026-07-27",
    part: 115,
    title:
      "Transport blocked, cap broken — compression certifies a margin-free step at n = 155,921 and chains of 10 steps; the remaining list is three points, only one an inequality",
    verdict: "TRANSPORT-BLOCKED",
    summary:
      "T115 (TRANSPORT-BLOCKED with a large compression gain, 26/26) — contract SCHUR.TRANSPORT: the two last bricks of the margin-free programme — transport the exact Schur complement between the ladder's non-nested grids, and break the computational cap by multi-resolution compression. (M1) The transport: Haynsworth's partial minimisation (yᵀSy = min_z [y;z]ᵀQ′[y;z]) makes λ_min(S) inverse-free and margin-free — a two-sided transport bracket on 11 real regrid pairs (ρ = 1.023–3.425). The error bracket is an IDENTITY (bookkeeping); the evidence sits in η and in the lower end. Operator control: the exact cell-overlap projection beats three deliberately wrong transfer maps 9/9 (by up to a factor 270). THE LOWER END IS POSITIVE ON 7/11 — a clean split: every pair with ρ ≤ 2.061 certifies, every pair with ρ ≥ 2.291 does not; synthetic threshold ρ* = 1.83, median real ratio 2.001 (39% of the ladder's 7686 refinement pairs covered). The reason is PRINCIPLED, not a bound problem: on nested ladders, where the transport error is exactly zero, λ_min(S) itself falls like ρ^(−1.73…−1.65) — the drop is real, so no bound can undo it; transport certifies only for mild refinement. Non-nestability is sharper than T114 thought: the criterion is integrality (not dyadicity), and 0/14807 ratios are integer (closest miss 3.4e-5). η measured 2.7e-2…8.1e-2; the a-priori Céa/Strang surrogate is 10³–10⁶× too coarse. (M2) The compression: the two-scale form (interior merged into blocks, boundary kept fine, the merge anchored at the centre) keeps X_mixed = Q_old,mixed EXACTLY (rel 0.0) — the compressed step is still margin-free; Albert certifies 66/66 (zone, q) combinations up to q = 64; the compression error is certified ONE-SIDED and second order in the projection defect (Rayleigh–Ritz plus stationarity of the fine minimiser, Céa/Strang, 66/66); savings down to m/h = 0.043. (M3) The deep run: THE CAP BREAKS — the margin-free step certifies at n = 155,921 (117× T114's 1331), on a fine lattice of h = 93,470 cells (62× the hard cap h ≤ 1500) compressed to m = 1490; re-coarsening is free (Rayleigh–Ritz, rel ≤ 2.1e-12); the LONGEST CHAIN IS 10 STEPS (T114: 4), 33 certified steps over 4 chains; the stopper is the cost cap on 3/4 chains and a failing step on 0 — every certificate sits 10⁵–10¹¹× above the Cholesky backward-error floor. Honest: λ_min(S_mixed) ≈ 5.1 is 52× the uniform-grid size (the coarser space measures a larger λ_min), flat over depth. (M4) The remaining list — the shortest ever: three points, only ONE of which is an inequality. (1) [P1] is ONE scalar inequality: ε = 1 − t̃ᵀT⁻¹t̃ ≥ c(D) > 0 — the classical Szegő–Levinson prediction-error bound for one symbol (Haynsworth double-Albert: T ≻ 0 ∧ ε ≥ 0 ⟺ Q ⪰ 0); measured ε ~ ρ^−1.71, falling faster than λ_min(S) ~ ρ^−1.69. (2) An a-priori η bound = regularity of z* (would make transport certified rather than measured for ρ ≤ 1.8). (3) A BOUNDARY FORMULATION (suggested by the exact-zero lemma): never represent the old interior at all — that would remove the cap AND the regrid. Contiguous reach (what an induction needs): n ≤ 125 → n ≤ 5437 (factor 43.5) — the n·log n barrier is DIVIDED, not broken. Deleted from T114's list: the margin, the (R) demand, floor/ε transport, the compression risk. T116 (BOUNDARY.FORMULATION, boundary_formulation_probe.py) is running: the induction as a pure boundary process (Riccati-type Schur recursion, pole via Woodbury, kernel-tail truncation), a deep boundary run, and the ε(D) law with the classical comparison. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "schur_transport_probe.py",
  },
  {
    date: "2026-07-27",
    part: 114,
    title:
      "Wall dissolves — the margin-free chain runs to n = 1331 and stops only at compute cost; the wall was an artifact of the requirement",
    verdict: "WALL-DISSOLVES",
    summary:
      "T114 (WALL-DISSOLVES, 22/22) — contract MARGIN.FREE.STEP: rebuild the chain without the margin division — does the margin wall disappear as a requirement artefact? It does. (L1) The margin-free step: the sub-block of the grown window is bit-exactly Q_old (the atom entry is the exact zero matrix) — margin-freeness IS that property. Route (i), Albert 1969 (the generalized Schur complement with the Moore–Penrose pseudo-inverse; the range condition is Douglas 1966), certifies 27/27 ladder steps, ELEVEN of them beyond the old wall (467, 479, 512, 529, 547, 577, 661, 773, 887, 1129, 1331; deepest h = 1495). The exact 4×4 Schur complement is an O(0.1) object with NO cancellation: λ_min(S) = 0.068–0.154, i.e. 42–67% of the block scale. Route (ii), the T110 graded minorant at x_in = 0, dies structurally 27/27 — the T110 margin was not paying for positivity but for an abandoned direction. Route (iii), unshifted Cholesky, passes 27/27. Negative controls are sharp (5/5); frame invariance holds at ν = 8 (3/3). THE WALL IN ONE LINE: the same quantity via the norm bound (Weyl: λ_min(A) − ‖C‖²/λ_min(X)) is NEGATIVE by a factor 2.4e5–9.6e7 — an O(1) numerator divided by a 1e-6 artifact floor — so every norm route had to fail, while the exact Schur complement passes through none of it. All seven zones where the T109 chain tore (including the wall zone n = 449 itself, need109/m = 1.241) are certified by the margin-free step. (L2) The ratio closure: ε ~ D^1.77–2.01 and κ ~ D^3.64–4.84 — the preregistered 'same power' is refuted IN THE FAVOURABLE DIRECTION: κ falls faster, so r = κ/ε ~ D^1.6–3.1 shrinks under refinement; the ratio closure r = 2e-4…0.103 holds on 27/27 steps (11/11 beyond 463). (L3) The circle [base PSD + margin-free step + ratio closure] closes 27/27; ten multi-step chains run (21 certified steps, longest 4); every chain ends at the cap h ≤ 1500 — a window COST — never at a step. Regrid honesty: 0/464 gap ratios are dyadic, so the non-nested transport object stays. (L4) NEW, [P5]: the (R) demand is GRID-BOUND — holding (μ/2)P₋ physically fixed and refining, ω grows monotonically (5/5) and crosses 1 (4/5): (R) is false in the continuum and is to be DELETED, not transported; the exact Schur complement needs no demand surrogate. Cancellation robustness: the 1e-7 cancellation has 5.3–8.8 decimal digits of headroom above the Cholesky backward-error floor and sits almost entirely in the rank-1 pole — whose positivity is exactly the T108 Szegő identity ε. THE NEW CORE is two objects: [P1] 'ε > 0 with controlled size relative to κ' (the exact Szegő identity), and [P2, sharpened] transporting the exact O(0.1) Schur complement between non-nested grids. T115 (SCHUR.TRANSPORT, schur_transport_probe.py) is running: a transport certificate plus multi-resolution compression, and a deep run beyond 1331. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "margin_free_step_probe.py",
  },
  {
    date: "2026-07-27",
    part: 113,
    title:
      "Substance confirmed — the margin wall is real in every currency, but it measures the discretization, not the spectrum",
    verdict: "SUBSTANCE-CONFIRMED",
    summary:
      "T113 (SUBSTANCE-CONFIRMED with a fundamental reinterpretation, 27/27) — contract LIMIT.OPERATOR: is the margin fall in the scaled frame a normalisation artefact, and what is the limit operator? (K1) The currency table: the floor/requirement ratio carries THE SAME exponent −1.168 ± 0.259 in ALL five currencies (raw, /λ_max, /trace density, /D, /D²) — structurally irremovable: the PWC basis is L2-orthonormal (Gram = I, no free D-power to spend), and the zone matrices are restrictions of ONE form to subspaces of the same L2 space (cross-grid form identity 1.2e-14). Floor homogeneity is exactly 1.0, need109 homogeneity 1.0 under legitimate rescaling. The wall is currency-invariant — substance, not an artefact of bookkeeping. (K2) BUT the substance is different from what was expected — two surprises. First, the limit object is NOT 'arch + pole with atoms as a perturbation': the atom-free form is deeply negative (λ_min = −3.08…−19.36, growing ~n^+0.60), the atom block is equally large (norm 3.21…19.44), and the positive floor survives only as a CANCELLATION of relative size 1.3e-7…9.7e-5 (on the floor vector: −11.64 against +11.64); the Weyl bound is saturated — norm perturbation theory is five orders of magnitude too coarse. Second, NO PLATEAU under refinement (fixed window, nested grids, 8 levels, a factor 128): λ₁ ~ D^1.83 and λ₂ ~ D^1.76 — the SAME power, λ₂/λ₁ flat. THE CONTINUUM WINDOW FORM HAS NO GAP; m_k measures the gap of the DISCRETIZATION. This resolves the T112 tension exactly (λ₂/λ₁ constant while the floor falls: both eigenvalues carry the same D-power). (K3) The atom part: the entering (deepest) atom contributes almost nothing to the coupling (8.3e-11…1.3e-6) — which is why handover steps certify at retention 1.000000 while the margin tears. The floor is almost pure geometry (log m = const + 1.87·log D − 2.07·log α; arithmetic residual a factor 1.13, no Λ(n) dependence). θ→0 stress: the operator norm does NOT diverge (θ^−0.005); an extra atom at the window edge shifts the floor only 0.2–0.4%, the same atom in the interior 10⁷ times more — the 1/θ explosion is purely a COST statement (h ~ νn/θ). (K4) Regrid: the raw rate (D′/D)^3.85 is mostly the floor's own D-power (a normalisation jump, not an instability); the reserve certifies f_crit = 9.5e-2 (better than the T112 estimate); after subtracting the D-power the reserve buys 15 steps (raw jump: 0). THE NEW HARDNESS BALANCE: [P1] semidefiniteness of the BALANCED form without a gap (strict limit-floor formulations demand what refinement refuses); [P2] certify the D-power itself (Grenander–Szegő type); [P3] NEW AND LEADING: need109 carries D^0.63 against the floor's D^2.38 — the T109 requirement chain divides by an artifact margin; the repair is a margin-free step certificate (Schur/Loewner from semidefiniteness alone); [P4] θ→0 is only a cost-uniformity demand. T114 (MARGIN.FREE.STEP, margin_free_step_probe.py) is running: rebuild the chain without the margin division (Albert/pseudo-inverse Schur, ratio closure) — does the margin wall disappear as a requirement artefact? Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "limit_operator_probe.py",
  },
  {
    date: "2026-07-27",
    part: 112,
    title:
      "Scaling partial — gap-coupled scaling fells two of three walls structurally; the margin wall is frame-invariant at exponent −0.974",
    verdict: "SCALING-PARTIAL",
    summary:
      "T112 (SCALING-PARTIAL, 20/20) — contract ADAPTIVE.SCALING: stop freezing the cell width D and couple it to the local prime gap — do the three T111 walls survive the n-adaptively scaled frame? (J1) Two couplings built: frame A gap-coupled, D_k = g_k/(2ν), and frame B mean-field, D_k = log n/(2νn). FRAME B HAS NO LADDER (only 38% of pairs nest; first failure already at n = 3) — the ladder lives exclusively in the gap coupling. T105 admissibility sets a resolution floor: ν ≥ 4 (at ν = 1, 3/10 zones are inadmissible). The flatness verdict: NOT flat — but the exponent is FRAME-INVARIANT: the components shift massively (m_k: n^−1.93 → n^−0.95; need109: n^−0.98 → n^+0.03, practically n-independent in the scaled frame!), yet the difference stays −0.974. (J2) The three walls, scaled: the LADDER WALL IS STRUCTURALLY GONE — the nested step has exactly ν cells per end BY CONSTRUCTION (twin pairs included), the atom entry is the exact zero matrix by construction instead of an arithmetic accident, and the T111 killer pair 461→463 now certifies spectrally (h = 1418, retention 1.000000) with the reserve OPENING: f_crit = 8.1e-4…4.2e-3 instead of 1.00 — roughly a factor 1000 (521→523 would need h = 1634 > the cap: arithmetically co-certified, spectrally declared honestly as not run). The ω WALL IS NO LONGER A DEPTH WALL: the depth coefficient is n^−0.079 ± 0.121, compatible with zero — ω is driven by resolution and the local gap, with a residual requirement-vacuum tail of 5 zones on the deepest ladder, not booked as gone. The MARGIN WALL STAYS, practically unmoved: the first sub-1 zone sits at n = 449 against the frozen frame's 461–463 (factor 1.07) — it is not geometry but the substance of the T109 requirement. (J3) The limit object: the arch (0.657) and pole (0.558) parts contract Cauchy-like with amplitude laws D^−0.26 / D^+0.37; the atom part does NOT contract (0.935 — it IS the prime-gap statistics in scaled coordinates); λ₂/λ₁ does not drift (n^−0.029 ± 0.052). Formulation: [a deterministic limit shape: arch + rank-1 pole] + [atoms as a controlled, non-convergent perturbation]. New and uncomfortable: the REGRID OBJECT — the zone frames are not nested, Rayleigh–Ritz transfers nothing between them; measured rate (D′/D)^+2.93; in the frozen frame this object does not exist. (J4) The trade: three arithmetic walls exchanged for two analytic demands — [P1] positivity of the limit operator (ONE operator instead of a matrix sequence) and [P2] a convergence rate below the step reserve. The prime-gap inputs are disclosed: Bertrand–Chebyshev 1852 (used, verified) and the trivial even-gap bound suffice upward; θ_k = g_k·n/log n = 0.18…4.12 is the only remaining arithmetic parameter, and downward θ is provably unboundedly small (Zhang 2014 / Maynard 2015: infinitely many bounded gaps force θ → 0 on subsequences) — any gap-coupled construction must be uniform in θ → 0, at the price of the 1/θ cost explosion. T113 (LIMIT.OPERATOR, limit_operator_probe.py) is running: the currency question (is the margin fall in the scaled frame a normalisation artefact? λ₂/λ₁ constant vs m ~ n^−0.95), the limit operator made explicit, the atom perturbation as a gap functional, regrid vs reserve. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "adaptive_scaling_probe.py",
  },
  {
    date: "2026-07-27",
    part: 111,
    title:
      "Crossing confirmed — the deep ladder measures the wall at n* ≈ 462 and splits it into three; the mechanism itself never fails",
    verdict: "CROSSING-CONFIRMED",
    summary:
      "T111 (CROSSING-CONFIRMED, 23/23) — contract DEEP.ZONE.STRESS: stop extrapolating and MEASURE — does T110's induction circle survive deep zones? (I1) The deep ladder: the T110 grid is reproduced bit-exactly, then the ladder samples ALL prime-power zones to 1000 plus a thin tail to 3001 (199 zones; chain ladder n = 2…521 with 117 handovers). The ratio m_k/need109 falls 179.59 → 0.0174, and THE CROSSING IS MEASURED, NOT EXTRAPOLATED: the last zone ≥ 1 is n = 461, the first below is n = 463 (0.935) — T110's n ≈ 170 was a fit artefact, but the wall is real. The failure follows the PRIME BRANCH (prime powers 512/529/625/729 with small μ survive longer). Fit families compared honestly (best by AIC: broken scale, n* = 501; jackknife band of the pure power [298, 352]); the depth-preserving resolution study pushes the ratio DOWN under refinement — n* ≈ 462 is an UPPER bound. (I2) The uniformity object: nsoft* = 1…8 across the 117 handovers (no longer flat at 1), but retention stays 1.000000 at every step (largest single loss 5.96e-08), and the atom entry is constructively free 117/117 (exact zero matrix). New and harder: the LADDER WALL — the nested step needs a lag gap > 2D, and the twin-prime pair 521→523 (gap 0.003831 < 2D = 0.005617) kills the ladder arithmetically; at frozen cell width n = 521 is the deepest reachable zone — the bound is set by twin primes, not spectra. An arithmetic scan to 40000 finds 4052/4284 pairs endangered; a co-moving D ~ 1/n forces h ~ n·log n/4 (horizon n ≈ 690 at h ≤ 1500). (I3) Reserve and circle: the no-reserve finding softens — f_crit ~ n^−0.39 (fit), the reserve OPENS with depth (the bottleneck stays the first step 2→3); the circle runs from the certified base case (factor 179.6) and BREAKS FIRST AT k = 109, n = 463 — exactly where I1 measures the crossing. IT TEARS AT THE RATIO, NOT AT A STEP: all 117 handovers certify with retention 1.000000, even beyond the break. (I4) The balance: the driver is 1/κ (+0.86) against the margin (−0.68); THREE SEPARATE WALLS — the margin wall n ≈ 462 (measured, upper bound), the ladder wall 521→523 (twin primes, purely arithmetic), and the requirement wall n = 727 (ω_cert ≥ 1, need109 vacuous on 46 zones). The reinterpretation: the operating variable is DEPTH, not n — doubling the cheeks moves the wall from 463 to 47; a proof must beat the exponent gap −0.974 = −0.681 − 0.293 and couple δ₀ to the local prime gap. T112 (ADAPTIVE.SCALING, adaptive_scaling_probe.py) is running: the construction in the n-adaptively scaled frame (D ~ local gap) — do all three walls fall at once? Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "deep_zone_stress_probe.py",
  },
  {
    date: "2026-07-27",
    part: 110,
    title:
      "Margin propagates — the induction circle closes end-to-end on the measured zones; three sharp gaps remain",
    verdict: "MARGIN-PROPAGATES",
    summary:
      "T110 (MARGIN-PROPAGATES with a precision statement, 28/28) — contract MARGIN.PROPAGATION: does the one strict-margin input that T109 left — λ_min(Q|odd) above an explicit need109 — propagate itself through the induction? On the measured zones the circle closes END-TO-END: certified base case, 15 certified handover steps, atom entry structurally free. (H1) The margin map on a common grid (16 zones, δ₀ = 0.00284): λ_min(Q|odd) falls monotonically 6.95e-2 (n = 2) → 1.66e-5 (n = 29) against need109 = 3.09e-7…3.87e-4 — m_k ≥ need109 on 16/16 (ratio 2.09–179.6); dilution during pure window growth is weak (factor 1.02–2.05 per atom gap = 0.993–0.997 per cell, a 1 − O(1/M) loss; the fits are honestly poor, power law wins 7/16). THE ATOM ENTRY COSTS NOTHING: it RAISES λ_min on 15/15 (the Kato first-order term has helping sign); at n = 29 the atom-free window is not even positive (−9.2e-3) — the atom makes the window co-positive. The T106 killer argument has no analogue for the scalar floor. (H2) The step law: the growth step splits exactly (embedding error ≤ 2.6e-14); max|μN| = 0.0 on 15/15 — the new atom lies outside the old lag reach (gap 0.0606 vs δ₀ = 0.00284), the restriction to the old window is the exact zero matrix; everything is decided by the bordering. Scalar routes fail decisively (bordered Weyl 0/15; the Schur/Friedrichs cap is vacuous 15/15). What passes certified is the GRADED LOEWNER MINORANT (its validity is exactly the induction hypothesis, then bisection through Cholesky/Sylvester): nsoft = 1 holds 15/15 with retention 1.0000; strictly scalar (nsoft = dim) 0/15. (H3) Base case CERTIFIED: n = 2 by explicit Cholesky at three resolutions (λ_min ≥ 6.93e-2, factor 88–180 over need). BUT the trend is unfavourable: m_k ~ n^−1.93 against need109 ~ n^−0.98, ratio ~ n^−0.96 (fit, rms 0.735) — extrapolated crossing at n ≈ 170: propagation gets HARDER with k. (H4) The circle test: strictly scalar induction breaks at k = 3 (the Schur factor compounds to 2.7e-36); the GRADED chain runs through COMPLETELY from the certified base value — 16/16 zones above need109, each step one Cholesky, with the propagated floor as input. THE THREE GAPS, sharp: (1) no reserve — f_crit at the first handover is 1.00 and the chain lives on retention 1 − 2e-7; a factor-2 loss would break it within two steps — what is needed is a step law WITHOUT loss; (2) a scalar step law is structurally excluded (boundary layer: the new cells attach at the edge where the pole source sits); (3) k-uniformity is missing (every certificate is one finite Cholesky; flat-in-k is a measurement, not a theorem — the n^−0.96 warning). T111 (DEEP.ZONE.STRESS, deep_zone_stress_probe.py) is running: the deep ladder to n ~ 200, testing the n ≈ 170 crossing directly. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "margin_propagation_probe.py",
  },
  {
    date: "2026-07-27",
    part: 109,
    title:
      "Boundary certified — both remaining scalars certified; the chain for (R) closes 16/16 on exactly one strict-margin input",
    verdict: "BOUNDARY-CERTIFIED",
    summary:
      "T109 (BOUNDARY-CERTIFIED, 29/29) — contract BOUNDARY.DECAY: both remaining scalars of (R) are now CERTIFIED — ω unconditionally via a graded matrix cap, the boundary value via a residual certificate that carries the cancellation — and the whole chain closes 16/16 conditional on exactly one strict-margin input that is 10²–10⁶ weaker than the conclusion. (G1) The mechanism is NOT decay but CANCELLATION: support separation is geometrically excluded (t̃ has its maximum at r = 0 — the source sits ON the boundary); the Combes–Thomas hypothesis is satisfied (Jaffard class, local rate ρ(s) = D(1/2 + 2/(e^{2s}−1)), factor 0.85–1.00 as a fit) but its conclusion is empty: the Green sum for x₀ cancels by 44…4.6e4 on 15/16 zones (n = 2, the known outlier, does not cancel), with a nearly linear boundary zero r^0.77…1.28 (fit) — a boundary condition, not a decay length. (G2) Four x₀ certificate candidates: (i) Combes–Thomas REFUTED — tested with the exact Green row, so no constant can rescue it (1/16); (ii) T-metric Cauchy–Schwarz at the edge 1/16; (iii) Szegő/Levinson exact (the two-point evaluation exhibits the cancellation) but yields no inequality; (iv) the residual / goal-oriented certificate (Prager–Synge, second order in the residuals) carries the cancellation instead of bounding it away — sharpness 1.0000–2.4937, 16/16 within budget. (G3) ω CRACKED, UNCONDITIONALLY: a graded matrix cap in PSD order (one level per soft direction, Cholesky-verified, ntop = 17–512) gives ω_cert = 0.2561–0.7569 against measured 0.2366–0.7531 (factor 1.004–1.124), WITHOUT any hypothesis input; the compression Schur distance (the T108 blockade) shrinks from 0.18–0.52 to 0.0024–0.039; the ungraded variant is vacuous everywhere — the grading is the whole difference. (G4) Both scalars come from ONE object (Γ_cert = S_cert⁻¹); the chain closes 16/16 with margin 1.3–10.4 and is M-stable (44/48 ladder points up to M = 3000; the certification price 1.5–28.6 does not grow). WHAT REMAINS is exactly ONE input: the STRICT positivity margin of the odd channel, λ_min(Q|odd) ≥ 1.7e-6…7.5e-3 (measured 1.4e-5…5.1e-2) — 10²–10⁶ BELOW μ/2, strictly weaker than the conclusion; the induction so far delivers only the non-strict ⪰ 0. Caveats: trial data numerically generated (validity is trial-independent via identity + Cholesky + Loewner), floating point not audited, λ_min via Rayleigh–Ritz. T110 (MARGIN.PROPAGATION, margin_propagation_probe.py) is running. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "boundary_decay_probe.py",
  },
  {
    date: "2026-07-27",
    part: 108,
    title:
      "Epsilon identity — ε's positivity is an exact identity; (R) drops to two scalars, literally one boundary value",
    verdict: "EPSILON-IDENTITY",
    summary:
      "T108 (EPSILON-IDENTITY, 44/44) — contract RATIO.CERTIFICATE: epsilon's positivity is now an EXACT identity that coincides with the induction positivity itself. (F1) Everything hangs on the vector x := T_odd⁻¹t̃ (the pole response). Three EXACT identities (verified to 10⁻¹¹…10⁻¹³): ε = xᵀQ|odd x/τ — ε IS a Q|odd energy; 1/ε = 1 + t̃ᵀQ|odd⁻¹t̃; and the overlap of ε and κ is exactly ONE rank-1 term hhᵀ/ε, so (R) ⟺ a generalized Rayleigh quotient r̂ = (μ/2)λ_max(Γ + hhᵀ/ε) ≤ 1. T107's Cauchy–Schwarz chain turns out to be exactly the Weyl bound (slack 1.0000). The soft-edge picture is CORRECTED: τ sits to 99.66% on modes λ ≥ 1, and ‖h‖² is small through sign cancellation (Σ|k| = 401 vs Σk = 1). r̂ is flat over M (6.1× more stable than ε); r, r̂ < 1 at all 92 (zone, M) points up to M = 3000 (tightest 0.9648). (F2) ω = 0.2366…0.7531 measured 16/16; the T105 structure survives the parity split literally (VᵀS|odd V = −½I at 4.9e-14); ω_cert is blocked solely by the compression Schur distance (0.18…0.52). Three ‖h‖² certificate candidates die quantified (cancellation 401, Bessel 10⁷, pole-CS 10⁶). (F3) The ε route (i) carries: x is practically the ground mode of Q|odd (Rayleigh 0.9248…0.9993); the fully explicit form ε ≥ λ_ind‖t̃‖²/(t̃ᵀT_odd t̃) costs only 1.51…9.81. Route (iii) delivers the classical structure: t̃ is exactly a two-term geometry, τ a 2×2 Christoffel–Darboux form, and ε is EXACTLY the square of the last Cholesky pivot (the Szegő/Levinson prediction error) — ε's positivity coincides with the induction positivity in the odd channel. (F4) The chain (R) ⟸ [ω < 1] ∧ [λ_min(Q|odd) ≥ (μ/2)τθ/(1−ω)]: C1 closes 16/16 (margin 5.5–178; ≥ 1 at all 48 ladder points: 2.7–343), C2 15/16, C3 0/16 — and for the FIRST TIME the chain is NON-VACUOUS: the required induction margin is explicit at 9.3e-8…3.7e-3, a factor 10²–10⁶ BELOW μ/2. REDUCTION: (R) drops from an M/2-dimensional Loewner statement to TWO scalars — ω < 1 and the statement at the vector x. Exactly open: a certified upper bound on the avoidance norm ‖VᵀT_odd⁻¹t̃‖² — on 8 zones literally the SINGLE BOUNDARY VALUE x₀ (suppressed: |x₀|/max|x| = 2.9e-3…5.7e-3, linear edge profile r^1.01 as a fit) — a boundary-decay statement about an explicit vector, the object class of T105's support separation. T109 (BOUNDARY.DECAY, boundary_decay_probe.py) is running. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "ratio_certificate_probe.py",
  },
  {
    date: "2026-07-27",
    part: 107,
    title:
      "Scalar tractable — the last matrix statement is exactly one scalar ratio, with two orders of magnitude of room",
    verdict: "SCALAR-TRACTABLE",
    summary:
      "T107 (SCALAR-TRACTABLE, 30/30) — contract ODD.CHANNEL.CLOSURE: the last matrix statement is now exactly ONE scalar ratio. (E1) The Sherman–Morrison reduction stands EXACTLY: (R) ⟺ s* ≤ 1 (identity 3.4e-11; equivalence confirmed against Rayleigh–Ritz 16/16; precondition G ≻ 0 via Cholesky 16/16). BUT s* is SATURATED: 0.99885…0.99999 (only n = 2: 0.2569) — any uniform s* chain is hopeless. The Woodbury split rescues it: s* = τ + κ, and (R) ⟺ r := κ/ε ≤ 1 with ε = 1 − τ; measured r = 0.0053…0.1814 — TWO ORDERS OF MAGNITUDE of room instead of five decimal places. s* is resolution-stable (spread ≤ 2.60% over M = 600/900/1200). (E2) The symbol route is structurally dead: the T_odd lower bound via the symbol is certified but NEGATIVE (λ_cert = −2.54…−0.81) — measured λ_min(T_odd) = 4.7e-5…6.4e-2 against a symbol infimum of −22…−6 (five orders of magnitude plus a sign flip). The soft edge of T_odd is a pure finite-section effect that no symbol sees — Grenander–Szegő cannot deliver in principle. (Side finding: the completion freedom is a real degree of freedom; continuum vs Dirichlet wins 6/16 vs 10/16.) (E3) t is exact in closed form (‖t‖² = 0.0154…1.84). Coarse certified chains (Cauchy–Schwarz against λ_min) close 0/16 — a category error, not a precision problem (they fail even with the measured λ_min). t follows the measure (ratio 1.21 between Fourier mass and circle measure on the deep levels) — no symbol-level avoidance is exploitable; the real avoidance lives between t and the demand space in the T_odd⁻¹ metric. The Cauchy–Schwarz bound on κ is practically lossless (factor 1.00079). Deep sweep to M = 3000: budget factor q = 1.73…2.54 (drift declared as a fit); ε ~ M^(−1.99…−1.69). (E4) Certified closure 0/16, measured 16/16 (q = 1.31…4.16 at M = 1200 — replacing T106's 0.747…1.156, a factor 1.7–2.5 gain). The bonus β₀ = 1.31 changes nothing about the uniform chain (the obstruction is the soft edge); r_CS ≤ 1 holds 16/16 with and without it. EXACTLY ONE OBJECT is missing: a certified positive lower bound for ε = 1 − t̃ᵀT_odd⁻¹t̃ (how far the pole vector is from exhausting the odd channel) — or, since ε and κ both fall with M, a certificate for the RATIO r itself. The named path runs through ω (a T105-type wing statement) plus the avoidance norm ‖h‖², not through Szegő. T108 (RATIO.CERTIFICATE) is running. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "odd_channel_closure_probe.py",
  },
  {
    date: "2026-07-26",
    part: 106,
    title:
      "Density mapped — the parity split localizes all remaining hardness in the odd channel",
    verdict: "DENSITY-MAPPED",
    summary:
      "T106 (DENSITY-MAPPED, 32/32) — contract FRIEDRICHS.ANGLE: two routes honestly killed, and a parity breakthrough that halves the remaining object. (D1) The density chain is dead for good: Landau–Widom is the WRONG classical anchor (rms 0.50–1.55); the right one is Grenander–Szegő on the Planck-coarsened symbol (averaging over the mode-spacing cell π/M; hits 0.92–2.43, better on 67/72 pairs). It explains why a negative symbol (dip −21.3…−3.4) still gives Q ⪰ 0: the atoms oscillate exactly at the window's Planck scale. Inertia: the Toeplitz part has exactly ONE negative eigendirection and the rank-2 Weil pole lifts it. A density certificate would need a hard gap w₀ ≥ 0.55–2.44; measured 5.3e-6…1.1e-3 — a factor 5·10²–3.6·10⁵ short. (D2) Invariant amplification is a NO-GO with a measured mechanism: the wing demand gives β_wing = 2.21–14.43 > 1 on 16/16, BUT inherited demands give β_int ~ 10⁻² on 15/15 — the accumulated invariant with absolute anchoring is FALSE. The reason is transport: in the larger window an old wing pair becomes an interior pair and draws 99.7–99.95% of its demand from the 64 softest modes (coupling 751–2653× stronger). No self-propagation (β_on < 1 on 15/15, loss 13–320×, θ_k = 0.9866–0.9998): the wing margin is an EDGE property of the window, not a property of the pair. What survives: the window-wise family σ_k(δ₀) ≥ 1.31·(μ_k/2), uniform over M (drift ≤ 2.01×). (D3) The breakthrough: the Weil pole splits EXACTLY as abᵀ+baᵀ = ssᵀ − ttᵀ (10⁻¹⁸; J a = b) — a positive rank-1 lift in the J-even channel and a negative rank-1 pressure in the J-odd channel. The one negative Toeplitz direction is J-even and is removed by ssᵀ exactly there; the odd channel never needed the pole. The soft edge belongs to the even channel 16/16 (3.7–100× softer). The DANGEROUS channel is the odd one (ρ₋ = 0.557–0.826 vs ρ₊ = 0.076–0.464) — the one with the better density: not the density but the avoidance decides the angle. T105's bare bound sits 16/16 exactly in the hard channel; ρ₋ is resolution-stable (≤ 17.6%). (D4) Channel-wise angle chain: even 16/16 (0.107–0.583× budget), odd 9/16 (0.747–1.156), unsplit only 3/16; all 16 handoffs close at operating depth in both channels under the exact Schur criterion; certified b0 = 89–93%. The remaining object is ONE Loewner statement on half the dimensions: (R) Q_full(α_k)|_{J=−1} ⪰ (μ_k/2)·P₋|_{J=−1} on ⌈p/2⌉ dimensions — Toeplitz part POSITIVE, soft edge 3.7–100× higher, bare certified. Certified: bare, parity split, pole splitting, channel inertia, avoidance, growth inequality, accumulation bounds. Measured: ρ₋. Refuted: accumulated invariant + propagation, density chain, Landau–Widom anchor. T107 (ODD.CHANNEL.CLOSURE) is running. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "friedrichs_angle_probe.py",
  },
  {
    date: "2026-07-26",
    part: 105,
    title:
      "One of two — the bare bound is certified in closed form; the avoidance law becomes a theorem",
    verdict: "ONE-OF-TWO",
    summary:
      "T105 (ONE-OF-TWO, 28/28) — contract BARE.AVOIDANCE.CORE: the chain's first certified closed-form ingredient. The T104 arm discrepancy dissolves — bare depends only on (u_k, δ) (the window centre α cancels; ≤ 0.74% drift under 8× refinement); arm A (4.39–26.77) is the corners of the (zone × γ) surface, arm B (2.06–4.60) its γ = 1 column; one currency, exact split bare = μ_k/2 + b0. bare is CERTIFIED: the exact E₋ structure Q_full|E₋ = (μ_k/2)Id + Arch + Pole (all foreign atoms vanish identically by support separation, slack ≥ 0.042), and three classical steps (Bessel; Legendre/level-set optimized at k(π/δ); Cauchy–Schwarz at the pole pair) collapse to the closed lower bound b0 ≥ (δ/π)∫₀^{π/δ} k(t)dt − δ·K_max − 4(cosh(u/2)−1)sinh(δ/2) — the band mean of the archimedean kernel over the band conjugate to the wing, ~log(1/2δ)−1. Positive 16/16 (1.27–3.54) at 81.1–92.7% of the measured value (median 92%), sharper the flatter the wing; no eigenvalue or induction data needed. The AVOIDANCE LAW is upgraded to a theorem: Q_full ⪰ 0 forces ‖Bᵀv_i‖² ≤ w_i·Λ₋ (semi-inner-product argument; coupling proportional to the mode's eigenvalue; verified at every mode of every zone, worst ratio 0.377) — it explains the measured 0.00% avoidance exactly. Plus an exact parity superselection: JQJ = Q (5.7e-15), JB = −BR (3.4e-15) — two channels that never mix; the softest mode is J-even 16/16. Euclidean angles are NOT small (cos up to 1.0) — the avoidance lives in the form metric: ρ_s is a Q-Friedrichs angle, resolution-stable (M-sweep: the margin falls 4×, ρ moves ≤ 8.5%). Synthesis: with certified bare the chain closes 15–16/16 in the handoff window (certifying costs zero zones; reach 0.44–0.80), and the whole remaining hardness is ONE Loewner statement: Q_full(α) ⪰ (μ_k/2)P₋ ⟺ ρ ≤ 1 − (μ_k/2)/bare_k — a Friedrichs-angle / spectral-density condition, not a margin condition. The free cap L ≤ Λ₋ = 5.63–7.52 diverges like log(1/D): the certificate must come from the spectral density near the soft end — L stays a measurement. T106 (FRIEDRICHS.ANGLE) is running. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "bare_wing_bound_probe.py",
  },
  {
    date: "2026-07-26",
    part: 104,
    title:
      "Chain partial — two independent arms; the margin route is dead, exact spectral splits close 16/16",
    verdict: "CHAIN-PARTIAL",
    summary:
      "T104 (CHAIN-PARTIAL, 21/21 + 47/47) — contract SCHUR.PROFILE.BOUND, run TWICE: a file collision between two parallel workers was resolved by keeping both versions, so the same contract has two independent implementations (arm B: schur_profile_bound_probe.py; arm A: schur_profile_chain_probe.py, “chain variant; independent twin”) — an unplanned independent replication, and all core findings converge. The naive margin route (M ≥ m·Id) is dead: 0/16 zones; the margin m = 2.4e-6…1.5e-3 vanishes in the continuum like M^(−1.7) (fit); coupling mass over margin is 10³…10⁶. The MECHANISM is avoidance: the coupling B₋ decouples from the bulk's soft modes (the softest mode carries 0.00% of the coupling mass for n ≥ 3; the mass sits in mid/high frequency bands, lowest band only 2–12%) — structural avoidance saves positivity, not a margin. Exact spectral-split chains close 16/16: arm A via the exact block-inverse identity W = m₀A⁻¹m₀ᵀ + RΣ⁻¹Rᵀ with matrix caps in PSD order (headroom H = 1.15…1.44, reach γ_max = 0.50…0.90, headroom trend flat −0.018 ± 0.003); arm B via a threshold split (O(1) threshold w = 2.47…5.90, resolution-stable — but needing r_min = 64…1024 soft modes ∝ M, so the continuum statement is a spectral-density condition). The dressing is nearly rank-1 (effective rank 1.07…2.7; near-null direction 99.7–100% in ONE prolate mode), BUT the prolate wing basis loses to the raw basis in both arms — the T96/T103 recommendation is tested and rejected; the two-scalar form (w, L) closes only 5/16. The hard core MOVES: a lower bound on bare_k = λ_min(Q_full|E₋), the soft dressing scalar L, finite induction data, and a one-sided edge estimate instead of a fit — the chain itself is then classical (Schur complement, exact block inverse, Parseval/Bessel, Weyl) without fits. T105 (BARE.AVOIDANCE.CORE) is running. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "schur_profile_bound_probe.py + schur_profile_chain_probe.py",
  },
  {
    date: "2026-07-26",
    part: 103,
    title:
      "Instrument improved — the race slope halves; closure jumps to 44/64",
    verdict: "INSTRUMENT-IMPROVED",
    summary:
      "T103 (INSTRUMENT-IMPROVED, 29/29) — door C, tool-building: the T101 race curve is exactly reproduced (slope −0.1622 ± 0.0562), then re-run with a θ-weighted band sum (certified weights, chain ρ ≤ b_band ≤ b_tail ≤ b_t99 at 64/64 samples) and full m(Λ) exploitation — Λ_ok stays bounded across all 16 zones (0.771…3.640 instead of 2.3…376; demand reduction 3.0–103×; honest price: explicit modes grow 2 → 232). New race slope −0.0748 ± 0.0116 (2.2× flatter, fit); r_k falls only 9.33 → 2.70 and never leaves the spectrum; the closure map jumps 7/64 → 44/64 with one fixed k-uniform instrument (Λ₀ = 3, r = 2). Measured verdict: the remaining loss is NOT in the bulk (θ-weighting and finite rank exhausted; the bulk is not low-rank — effective rank up to 0.579·dim E₋) but in the wing slack S = 1 − ρ (pencil nearly saturated; S falls 0.2091 → 0.0392) — next lever: a wing-adapted prolate/Slepian basis or a Fredholm shape of the equality argument. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "instrument_probe.py",
  },
  {
    date: "2026-07-26",
    part: 102,
    title:
      "Mechanism identified — the onset is manufactured by Schur dressing; C/g was a proxy",
    verdict: "MECHANISM-IDENTIFIED",
    summary:
      "T102 (MECHANISM-IDENTIFIED, 42/42) — door A, attack on the handoff law's lower bound: the mechanism is exact — the k-th atom acts on E₋/E₀/E₊ as diag(−1/2, 0, +1/2), so μ_k enters exactly once, linearly; a two-sided sandwich over the Schur profile σ_k(δ) brackets the handoff (2w_k between the crossings in 16/16 zones; the onset is anchored at δ_c = 2w_k, R² 0.968 — the T96 essential-singularity reading is compatible but not singled out). The binding constraint FLIPS: no concentration condition binds — the bare E₋ form is strongly positive (4–14× μ_k/2) and the classical ceilings (Cauchy–Schwarz, Landau–Pollak/prolate) are saturated near 97%; the onset is manufactured entirely by the Schur dressing against E₀⊕E₊ (35.7%…97.3% of the bare eigenvalue), i.e. by the coupling to the induction hypothesis. The decomposition test is triply negative: g_k is causally impossible, statistically dispensable, arithmetically only a ceiling (the extrapolation violates it from k = 69; checked over 18 120 prime-power atoms to n = 200 000) — the T101 law C/g was a proxy. The hard core localizes to one scalar per zone: a lower bound on σ_k(δ_ref) just above atom entry — probe-level typing of where the hardness sits, not progress on it. T103's wing slack S = 1 − ρ and T102's Schur dressing are the same object from two sides; T104 (SCHUR.PROFILE.BOUND) is running. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "arithmetic_bound_probe.py",
  },
  {
    date: "2026-07-26",
    part: 101,
    title:
      "Crowding trends — but the race is lost by the instrument, not the math",
    verdict: "CROWDING-TRENDS",
    summary:
      "T101 (CROWDING-TRENDS with location, 31/31) — the fork: 16 zones resolved (n = 2..29). HEADLINE: the collapsed law w_k = 0.0838·(atom gap)/μ_k (fit) — “the handoff window is the atom spacing divided by the atom strength” (scatter 8.85× → 2.66×, residual trend null). The PRIMITIVES are FLAT across 16 zones (w/g flat; relative D_k margin flat; D_k ≤ μ_k/2 holds 64/64 and never fails); ONLY the closing instrument loses (r ~ exp(−0.16k), fit/extrapolation). Core sentence: “The crowding sits in the proof family, not in the mathematics it is trying to prove — the most hopeful version of the verdict.” An asymptotics theorem would need (A) the arithmetic lower bound of the collapsed law [THE localized hardness], (B) uniform relative margin, (C) a better bulk instrument, (D) a finite check. All laws marked as fits; k > 16 as extrapolation. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "k_asymptotics_probe.py",
  },
  {
    date: "2026-07-26",
    part: 100,
    title:
      "The 100th probe — remainder closes zones 2–4; zone-5 tip is an equality",
    verdict: "REMAINDER-CLOSES-ZONES",
    summary:
      "T100 (REMAINDER-CLOSES-ZONES, 27/27) — the 100th probe: closure 11/24 → 24/24 (6/6 in every zone); the drift was a lattice artefact; one lever gained 1.7–69× (“the Bessel step threw the induction data away twice”); zones 2–4 fully closed; zone-5 tip typed as an equality problem (Fredholm shape, simple degeneration). Classical: Bessel/Parseval, Slepian, Schur test, Fredholm alternative. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "bulk_remainder_probe.py",
  },
  {
    date: "2026-07-26",
    part: 99,
    title:
      "Decay law found — parity selection rule; recursion terminates arithmetically",
    verdict: "DECAY-LAW-FOUND",
    summary:
      "T99 (DECAY-LAW-FOUND, 23/23): exact parity selection rule (J₋Q₋₀J₀ = −Q₋₀) — the fragile near-null mode is excluded from the binding channel by symmetry; recursive inequality with only 1.01–1.20× loss; termination is arithmetic (240/240 in ≤ 4 steps to the classical zone). Sandbox; not RH evidence. Classical: Schur, parity/reflection identities.",
    badge: "sandbox",
    script: "dk_recursion_probe.py",
  },
  {
    date: "2026-07-26",
    part: 98,
    title:
      "Law confirmed, mechanism open — circular lemma replaced by exact scalar identity",
    verdict: "LAW-CONFIRMED-MECHANISM-OPEN",
    summary:
      "T98 (LAW-CONFIRMED-MECHANISM-OPEN with REPLACED TARGET, 44/44) — the third self-correction of the weekend, the most valuable: the conjectured one-vector lemma was circular (Douglas range inclusion — “the law is forced by positivity itself and carries no independent information”); three T97 premises honestly refuted; REPLACEMENT: the exact scalar inequality D_k(α) ≤ μ_k/2 (identity, no constant, no vector) — holds in all four zones, lattice-stable, saturates exactly at the zone tips where the next atom takes over. PLUS both certificate upgrades: E₋ from 43% to 93% mean (the WHOLE zone in 3 of 4 — via the probability-measure identity: “the archimedean term on the wings is a MEAN of k against a probability measure; Slepian forces narrow wings to high frequencies where k is positive”); E₊ certified for the first time. Skeleton stand: 8 pieces PROVED, 2 certificates, 3 refuted, 3 open (the D_k bound, E₊ cross-blocks, certificate completion). Classical: Douglas 1966, Schur, Slepian–Pollak–Landau, Rellich, Rayleigh–Ritz, Richardson. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "cross_lemma_probe.py",
  },
  {
    date: "2026-07-26",
    part: 97,
    title:
      "Alignment-only with certified half-step — induction takes proof shape",
    verdict: "ALIGNMENT-ONLY",
    summary:
      "T97 (ALIGNMENT-ONLY with certified half-step, 105/105): the relay induction step gets proof shape — alignment is sharp (sign alignment ⟺ coupling window nonempty, without exception); the t=0 killer loss on the anti-bump space is PROVED (k_eff = (1−cos(tu))k(t), gain ×2–4.8); structure pearl: the E₀ block is literally the same form on the smaller window — “the induction hypothesis appears inside its own decomposition as self-similarity” (7e-14). Sandbox; not RH evidence. Classical: Schur complement, Paley–Wiener, Prolate.",
    badge: "sandbox",
    script: "relay_induction_probe.py",
  },
  {
    date: "2026-07-26",
    part: 96,
    title:
      "Edge artifact + relay confirmed — each prime arrives before it is needed",
    verdict: "RELAY-CONFIRMED",
    summary:
      "T96 (EDGE-ARTIFACT + RELAY-CONFIRMED, 21/21): (i) the T95 “edge” at α* was a map artefact — λ_min stays positive on all of [0.38, 0.86] (value exactly reproduced, interpretation corrected; “the margin never crosses zero — it collapses exponentially, λ ~ exp(−49α)”). (ii) THE RELAY IS REAL: without the log3-atom λ_min crashes to −0.445, the loser is exactly the anti-double-bump at distance log 3 (alignment −0.99), rescue identity to 5e-15 — and the handover windows are all positive: “every prime atom arrives strictly BEFORE it is needed (+0.025/+0.009/+0.011/+0.007 for the first four atoms).” (iii) Strategy shift: not an edge problem but a margin problem; numerics as witness exhausted past α ≈ 0.55; the counterfactual is the proof target (O(0.1)-sizes instead of 1e-6). Second self-correction of the weekend, same anchor discipline. Classical: Paley–Wiener, Prolate, Galerkin/Richardson. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "relay_test_probe.py",
  },
  {
    date: "2026-07-26",
    part: 95,
    title:
      "Continuum attempt — C1 fully proved; atom-extremal directions are safe",
    verdict: "T-CONTINUUM-NUMERIC",
    summary:
      "T95 (T-CONTINUUM-NUMERIC, 28/28): the continuum proof attempt — the C1 chain is FULLY PROVED (unconditional: |h_f(log2)| ≤ 1/2 via disjoint support intervals in the band; ‖S‖ = 1/2 exact with characterised eigenspace; the atom-extremal directions satisfy the target inequality WITH margin — “the directions that maximize the atom cost are provably safe”); the continuum margin curve stays positive everywhere; the extremizer is NOT the two-bump structure — the binding mechanism is atom rescue (confirmed in corrected coordinates); lower bound open, missing instrument structurally named. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "continuum_extension_probe.py",
  },
  {
    date: "2026-07-26",
    part: 94,
    title:
      "Blind demo — 753/753 primes predicted with zero errors, no division",
    verdict: "BLIND-100",
    summary:
      "T94 (BLIND-100, 16/16): the compiler predicted all 753 primes of a never-seen window [1,000,001–1,010,000] with zero errors — from pure lattice counting (odd n prime ⟺ r₄(n) = 8(n+1), Jacobi 1834 on the rank-|μ₄| theta tower), with an AST-enforced prediction path containing no divisibility tests (no isprime, no sieve, no %, //, / operators); predictions committed by MD5 before any truth was computed. Honesty: a simple sieve is ~820× faster — this is structure completeness, not algorithmic progress; spectral prediction (zeros → primes) remains bound to I5/RH. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "prime_blind_demo_probe.py",
  },
  {
    date: "2026-07-26",
    part: 93,
    title:
      "Reconciliation — the self-check found the checker's bug; a_neg recalibrated",
    verdict: "MIXED",
    summary:
      "T93 (MIXED, 41/41): T92 reported two calibration discrepancies against the band map; an independent third implementation found the T89/T91 map SURVIVES — 6.98 was only a stale prose value (the probes always computed 6.2898); the second discrepancy was pure factor-2 convention plus a 1% interpolation artefact (a_neg precision-improved: 0.7410 → 0.7486); the only real bug was in the checker T92 itself (constants imported untranslated — it accidentally certified a harder four-atom region). Survival list: one-atom zone, atom rescue, a* = 0.9253, uncertainty constants, REST vectors — all stand. Sentence: “The self-check found the checker's bug — and precision-improved one constant by 1%. The anchor discipline works in both directions.” Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "band_reconciliation_probe.py",
  },
  {
    date: "2026-07-26",
    part: 92,
    title:
      "Zone-extension skeleton — finite blocks are the wrong instrument",
    verdict: "T-SKELETON",
    summary:
      "T92 (T-SKELETON, 36/36): machine-assisted proof attempt of the zone extension — certified is Q ≥ 6.7e-12 > 0 on the 8-dim window subspace over the whole scanned region (error certificate, 1057-point covering); the full extension does NOT stand: λ_min collapses geometrically (factor ~11 per mode), the complement would need ~5000 modes — “the finite-block route is structurally the wrong instrument.” Pearl: k(0) = −γ − 3log2 − π/2 − log π exact. (Calibration flags against T89/T91 were later resolved by T93 as the checker's bug.) Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "zone_extension_proof_probe.py",
  },
  {
    date: "2026-07-26",
    part: 91,
    title:
      "The band's law — the first prime rescues positivity; named target inequality (T)",
    verdict: "BAND-PARTIAL",
    summary:
      "T91 (BAND-PARTIAL, 19/19, 3/4 closed): the band is the one-atom zone log 2 < a ≤ 0.9253 with inner edge a_neg = 0.7486 (recalibrated by T93; was 0.7410) — and the geography's cleanest finding: beyond a_neg the prime-free margin changes sign and THE PRIME ATOM BECOMES LOAD-BEARING — the first prime rescues positivity where the archimedean margin is exhausted (the T89 balance point is exactly this sign change). Atom turn-on law exact (k = 2m+1, Beta integrals); uncertainty constants decided: a·t_rms → π (Wirtinger) and a·t_cent → 2Si(π) − 4/π = 2.4306 (new closed constant); band and tight curves are two orbit regions of the same functional with shared exact scale ∫k_ζ = 2θ_RS (Riemann–Siegel). Named target inequality (T): (P_pole + A_arch)(f) ≥ √2·log2·h_f(log2) on the band — provable-shaped as a ZONE EXTENSION beyond Bombieri's log 2 (a self-standing classical target!); RH ⇒ (T), (T) ⇏ RH. Honest open: the super-exponential λ_pf rate remains empirical. Classical: Wirtinger/Rayleigh, Beta integrals, Lambert-W, Si integral, Bombieri, θ_RS. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "thin_band_analytic_probe.py",
  },
  {
    date: "2026-07-26",
    part: 90,
    title:
      "The residual dissected — explicit Gaussian modes, and 3 of 10 vectors under no control",
    verdict: "CORE-DISSECTED",
    summary:
      "T90 (CORE-DISSECTED, 17/17): the I5 core is now explicit — the 1–4 residual vectors are n-stable (≤ 2.1°) and have closed shape (Gauss×cos / Gauss×sin modes, 99+% capture); coverage matrix: 3/10 vectors controlled by no known structure (neither certificate extension nor CC pole-vanishing nor Sonin projection); the core question decided: the core is NOT only pole coupling — from a ≈ 0.85 genuine pole-free atom↔arch content remains, plus an odd atom-coupled mode at a = 1.2. Final requirement line: an I5 idea must deliver positivity for an explicit family of Gaussian-modulated modes in the thin band — and it is not reducible to pole cleanup. Classical cited: Bombieri, Connes–Consani vanishing, Slepian. Sandbox; no spectral identification; I5 remains ⟺ RH; not RH evidence.",
    badge: "sandbox",
    script: "residual_dissection_probe.py",
  },
  {
    date: "2026-07-26",
    part: 89,
    title:
      "Crossing mapped — no sharp wall at log 2; a thin band and a 1–4-dim residual",
    verdict: "CROSSING-MAPPED",
    summary:
      "T89 (CROSSING-MAPPED, 19/19): window compression of the Weil form is exactly Bombieri's classical object (Bombieri 2000 — unconditional positivity is classically known precisely up to support width log 2: Yoshida 1992, Bombieri 2000, Connes–Consani 2021). (i) Proven zone reproduced on our objects (prime side identically zero below log 2, λ_min ≥ 0 everywhere). (ii) The boundary is NOT sharp — margin falls smoothly (classical band-limitation), no collapse at log 2; atom turn-on is (a−log2)³-soft. (iii) The minimal object: the residual subspace controlled by neither program grows softly 0 → 4 dimensions (of 32) across the crossing and carries the global margin minimum — “the I5 core, for the first time, as a concrete finite-dimensional object per window width.” (iv) The real attackable content is the atom↔arch balance in the thin band log 2 < a ≲ 1.0; the full-form minimum lives exactly in the pole-coupled DC direction that Connes–Consani exclude by vanishing conditions. Honest self-correction: BOUNDARY-SHARP retyped (“small at the boundary” ≠ “collapses at the boundary”). Fence: margins beyond the proven zone are measurements; no spectral identification; I5 remains ⟺ RH. Classical: Weil 1952, Guinand, Yoshida 1992, Bombieri 2000/2003, Connes–Consani 2021/2023, Connes–Consani–Moscovici, Suzuki. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "sonin_crossing_probe.py",
  },
  {
    date: "2026-07-26",
    part: 88,
    title:
      "I5 tightness landscape — eight parametrized curves carry the RH content",
    verdict: "TIGHT-SET-PARAMETRIZED",
    summary:
      "T88 (TIGHT-SET-PARAMETRIZED, 27/27): the I5 tightness landscape — band-limited zero plateau (support below log 2, consistent with the T87 sector boundary), safe zones, and 8 nearly vertical tight curves whose spacing follows the Γ-side density law (ratio 1.01 ± 0.08 — described via smooth Γ-density, zero-free; no spectral identification). Zero negatives on true autocorrelations (earlier T76 negatives were missing arch/p=2 bookkeeping, as the ledger predicted). Validation at machine precision (null test 4.3e-13, ledger 1.1e-14). Conclusion: the RH content of I5 concentrates on low-dimensional, parametrizable curves. Milestone: 2001/2001 sandbox checks. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "i5_tightness_probe.py",
  },
  {
    date: "2026-07-26",
    part: 87,
    title:
      "Connes dictionary exact at the core — two positivity programs are complementary",
    verdict: "DICTIONARY-EXACT-CORE",
    summary:
      "T87 (DICTIONARY-EXACT-CORE, 22/22): the Connes dictionary is exact at the core — Q_cert atoms = Connes' finite orbit terms (rel 6e-16); A_fam − A_shift = Connes' archimedean W_∞ term including principal-value constant (rel 3e-23); the internal kernel k_ζ = 2θ'_RS (the Riemann–Siegel phase derivative, rel 0.0); the dihedral shift = the scaling group R*₊ (Weyl commutation exact). Sharpest finding: the two positivity programs are complementary — Connes–Consani prove positivity BEFORE the primes (prime-free Sonin window, boundary exactly at the first prime atom u = log 2); the compiler AFTER the primes (atom certificates); the sectors touch, do not overlap; the I5 coupling lives in the crossing region outside both. Transferable-shaped: Sonin compression (same kernel, same group — verified). Classical: Connes 1999, Connes–Consani 2021, Weil, Guinand, Bombieri, Burnol, Yoshida. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "i5_connes_dictionary_probe.py",
  },
  {
    date: "2026-07-25",
    part: 86,
    title:
      "Q-pairing closes the non-coherent tail — matching lemma on all atom classes",
    verdict: "LEMMA-FULLY-CLOSED",
    summary:
      "T86 (LEMMA-FULLY-CLOSED, 29/29): the Q-pairing involution (d ↔ qd at the least 3-mod-4 prime factor) closes the non-coherent tail in a provably-shaped way — flip law 0/146 987 violations; weight control 4|ψ| = E·Θ exact; paired certificate max 0.836 < 1 (1.97M pairs). Structure pearl: the unpaired remainder is exactly the Z[i]-norm sublattice = the g_λ support — “the pairing hands its remainder exactly to the λ-channel.” Class cover exact: no clash atom uncovered. Honest: paired envelope better-divergent (not convergent); divergence survives only in coherent directions — exactly what the λ-channel closes. Fence: window-format, named classics, not RH evidence.",
    badge: "sandbox",
    script: "correlated_cancellation_probe.py",
  },
  {
    date: "2026-07-25",
    part: 0,
    title:
      "Promoted: matching lemma and transport ledger package (v541) — two named limits as content",
    verdict: "MACHINE-VERIFIED",
    summary:
      "The T78–T85 proof package is now load-bearing (33 checks, ~10 s; every core check recomputed, not cited). (A) The matching lemma is PROVED exact-integer on [4, 10⁶]: 7S < 40A at every atom (0 violations over 939 870 clash atoms; exact margin X = 0.082159, ρ_crit = 1.144; structure laws — 8-law, seed towers, Cohen seeds −48/−80/−8 — at 0 tolerance). (B) The transport ledger closes exactly: Q_Weil = Q_cert + Δ_arch + Δ₂ with Δ_pole ≡ Δ_conv ≡ 0 proven; the odd-prime side of the Weil functional equals the certified plus combination exactly. (C) The signed envelope is character-exact: −ψ = (χ₋₄ + ¼χ₈ + ¼χ₋₈)·Θ; confinement is a set equality. (D) The archimedean term is internal: Δ_arch is exactly the Γ-difference via Legendre duplication ((2π)^{−s}Γ(s) = ½Γ_R(s)Γ_R(s+1)); the kernel identity k_ζ = K_fam − K_shift holds pointwise. (E) The coherent class (= Z[i]-norms, set equality) is closed by the λ-equivariant channel (CM carrier g_λ exact in two independent routes; canonical phase mean μ₁ ∈ [−1,1]; λ-certificate 0.065 vs 0.236 unlifted). THE TWO NAMED LIMITS ARE THE CONTENT: (i) one classically-shaped OPEN lemma remains — the correlated-cancellation lemma on the non-coherent uniform tail (provably-shaped ≠ formal proof); (ii) I5 in one-family form (Q_cert + Δ₂ + A_fam − A_shift ≥ 0) is the single remaining TFPT-specific object — by the closed ledger EQUIVALENT to Weil positivity ⟺ RH (an equivalence typing, no progress claim). NOT 'almost RH'; ZETA.HP.CARRIER untouched; no marker moves.",
    badge: "machine-verified",
    script: "v541_matching_lemma_ledger.py",
  },
  {
    date: "2026-07-25",
    part: 85,
    title:
      "λ-equivariant design closes the coherent class — lemma end-status",
    verdict: "LEMMA-CLOSES-LAMBDA",
    summary:
      "T85 (LEMMA-CLOSES-LAMBDA, 24/24): the λ-equivariant design closes the coherent class in a provably-shaped way — CM carrier g_λ exact (two independent routes, Hecke–CM rules, support = Z[i]-norms as set equality); canonical phase mean μ₁ ∈ [−1,1] exact; λ-window certificate 0.065 vs 0.236 unlifted (3.6× margin); λ-ledger 90/90 including all 7 earlier gap rows; lifted chain never crosses. Lemma end-status: window proved (10⁶) + signed structure exact + demand conditioning exact + coherent class closed via λ-channel — remaining ingredients: proven classics (Hecke 1918/1920, L(1,λ)≠0, Cornacchia, Dirichlet, Mertens-AP, Landau, Gronwall/Robin, Cohen, Alaoglu–Erdős) + ONE classically-shaped open lemma (correlated cancellation, non-coherent uniform tail). Fence: provably-shaped, not a formal proof; I5 untouched; not RH evidence.",
    badge: "sandbox",
    script: "lambda_equivariant_design_probe.py",
  },
  {
    date: "2026-07-25",
    part: 84,
    title:
      "Gaussian lift works unanchored — last gap sits in the compiler's Z[i] home",
    verdict: "LIFT-WORKS-UNANCHORED",
    summary:
      "T84 (LIFT-WORKS-UNANCHORED, 29/29): the coherent class equals primitive Z[i]-norms exactly; Grossencharacter phases replace Mertens divergence by L(1,λ)-convergence — the lifted chain never crosses; frontier jumps from ~10²³ to ~10^{5.9·10¹²}. Missing: λ-equivariant certificate design (embedding obstruction exact: ψ all-minus vs λ mixed). Circle closes: the last gap sits in the Z[i] sector — the series' origin object (μ₄-glue) — and exactly there the compiler's own character structure supplies the control that is provably impossible over ℚ. Sandbox; not RH evidence. Classical: Hecke 1918, Grossencharacter L-functions, Landau.",
    badge: "sandbox",
    script: "gaussian_lift_probe.py",
  },
  {
    date: "2026-07-25",
    part: 83,
    title:
      "FE symmetrization is absorbed — the wall is FE-transversal, not positional",
    verdict: "INVARIANT-NULL",
    summary:
      "T83 (INVARIANT-NULL, 27/27): FE symmetrization is fully absorbed — 5/24 and λ* already are the numbers of the symmetric sector (the test region was right). Depth find: the product of the two FE reflections J₁∘J_{1/2} = e^{±u} is the unit-line shift — the centre delta is literally the transport operator; the value side is invariant under the whole infinite-dihedral ladder, the spectral cone only under its own reflection — “the wall is FE-transversal, not FE-positional.” Bonus: explicit-formula null test ≤ 2e-12 validates the whole convention. Sandbox; not RH evidence. Classical: Mellin involutions, infinite dihedral group, explicit formula.",
    badge: "sandbox",
    script: "fe_symmetrization_probe.py",
  },
  {
    date: "2026-07-25",
    part: 82,
    title:
      "The archimedean term was never external — I5 becomes one-family self-consistency",
    verdict: "ARCH-INTERNAL",
    summary:
      "T82 (ARCH-INTERNAL, 22/22): Δ_arch is exactly the internal Γ-difference via Legendre duplication ((2π)^{−s}Γ(s) = ½Γ_R(s)Γ_R(s+1)); battery 18/18 rel 6.6e-15. The family carries its Γ factor as the Mellin signature of its heat-sum nature. Consequence: I5 changes type — new form Q_cert(h) + Δ₂(h) + A_fam(h) − A_shift(h) ≥ 0: self-consistency of one heat family (atom expansion vs Mellin signature of the same theta objects). Nearest classical relative: Connes 1999 semi-local trace-formula positivity (context, not used). Fence: type change ≠ proof; I5 remains ⟺ RH. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "heat_arch_probe.py",
  },
  {
    date: "2026-07-25",
    part: 81,
    title:
      "Avoidance fails: coherent targets require coherent rescalings — theorem-shaped",
    verdict: "AVOIDANCE-FAILS",
    summary:
      "T81 (AVOIDANCE-FAILS, 31/31): the same multiplicativity that makes the lever exact proves the counter-lever: coherent targets are reachable only by coherent rescalings — the freedom does not exist on coherent demand (coh(m·k) = coh(m) ∧ coh(k), 0 exceptions to 10⁶). Salvage exact: the T76 recipe was already minimally coherent (0 unforced keys, 100/100); forced coherent clash closed per certificate (83/90, worst ratio 0.18); on avoidant demand the class vanishes identically. The probe chain ends here honestly: what remains is classical open problems or RH itself. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "recipe_coherent_avoidance_probe.py",
  },
  {
    date: "2026-07-25",
    part: 80,
    title:
      "Signed tail to ~10²³ — last gap confined to the χ₋₄-coherent class",
    verdict: "RESERVE-PARTIAL",
    summary:
      "T80 (RESERVE-PARTIAL, 29/29): the signed envelope is character-exact (−ψ = (χ₋₄ + ¼χ₈ + ¼χ₋₈)·Θ); signed certificate with 8.9× margin; tail closed to N* ≈ 10²³. Last gap exactly confined to the χ₋₄-coherent atom class (only 1-mod-4 prime factors — provably zero cancellation; Landau class density classical); on the superabundants the credits even dominate. Fence: Matching Lemma not yet complete modulo classical; I5 untouched; not RH evidence. Classical: Pólya–Vinogradov, Mertens in AP, Landau.",
    badge: "sandbox",
    script: "tail_correlation_lemma_probe.py",
  },
  {
    date: "2026-07-25",
    part: 79,
    title:
      "The transport ledger closes: the wall is one inequality, provably ⟺ RH",
    verdict: "LEDGER-CLOSES",
    summary:
      "T79 (LEDGER-CLOSES, 22/22): the transport ledger closes exactly — Q_Weil = Q_cert + Δ_arch + Δ₂ (identity ~7e-16 on 100/100); Δ_pole ≡ Δ_conv ≡ 0 proved; the odd-prime side equals the certified combination exactly. The wall is decomposed: 2 exact identities + 2 classically provable-shaped bounds + exactly one coupled core inequality I5 (the prime↔archimedean coupling), which by identity is equivalent to Weil positivity ⟺ RH. The wall is no longer a metaphor: it is one named inequality — equivalence typing, not progress toward proving it. The 2/6 negative T76 Q-values are shadows of missing Arch/p=2 bookkeeping. Sandbox; not RH evidence. Classical: Weil 1952, Guinand, digamma terms.",
    badge: "sandbox",
    script: "transport_ledger_probe.py",
  },
  {
    date: "2026-07-25",
    part: 78,
    title:
      "Matching Lemma machine-proved on [4, 10⁶] — tail ingredient named",
    verdict: "WINDOW-PROVED",
    summary:
      "T78 (WINDOW-PROVED, 25/25): the Matching Lemma is machine-proved on [4, 10⁶] — exact-integer inequality chain, full enumeration over 939 870 clash atoms, 0 violations, exact margin 0.082159; four structure laws at 0 tolerance. Tail honestly open: Robin 1983 + constants miss by factor 6.16 — missing residual ingredient named: a correlation lemma (thinning × cancellation). Sandbox; not RH evidence. Classical: Gronwall, Robin 1983, Cohen seeds.",
    badge: "sandbox",
    script: "lemma_window_proof_probe.py",
  },
  {
    date: "2026-07-25",
    part: 77,
    title:
      "Matching lemma is classically shaped — Gronwall/Robin suffice; transport wall untouched",
    verdict: "LEMMA-CLASSICAL-SHAPED",
    summary:
      "T77 (LEMMA-CLASSICAL-SHAPED, 21/21): the Matching Lemma reduces to a restricted σ₋₁ divisor bound (s = 1, Gronwall regime), h-free up to thinning; true margin stable ≤ 0.58 over windows to 8000; worst-case atoms are the classical superabundant suspects (5040, 7920 — Alaoglu–Erdős). Proof ingredients fully classical: 4-periodic sign law, weight-5/2 coefficient bounds, Gronwall 1913 / Robin 1983 unconditional (the RH-equivalent sharpening is not needed), Abel summation. Two-stage residue honest (raw envelope closes to ~10⁹; with thinning ~10^3131, then typed cancellation reserve). Fence: proof ingredients listed, lemma not proven; even a proven lemma yields only value-side representability — the transport wall remains the open wall. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "matching_lemma_structure_probe.py",
  },
  {
    date: "2026-07-25",
    part: 76,
    title: "A universal certificate recipe — and its named lemma",
    verdict: "RECIPE-UNIVERSAL-ON-BATTERY",
    summary:
      "T76 (RECIPE-UNIVERSAL-ON-BATTERY, 22/22): the hybrid recipe certifies 91/91 nontrivial adversarial Weil directions (100%); cost is polynomial / window-extensive (λ_m ~ m^{5/2} = Eisenstein law); lattice discreteness is the floor (δ_h → 0 does not break S1). Conjecture form with named core lemma: the MATCHING-LEMMA on the log lattice (a Diophantine divisor-sum problem). Implication architecture typed, NOT claimed: Lemma ⇒ value-side representability ⇒ [transport wall open] ⇒ Weil positivity. “Any RH content would relocate into a universality proof plus the transport wall.” Sandbox; not RH evidence. Classical: Weyl, Euler–Maclaurin, Farkas.",
    badge: "sandbox",
    script: "hybrid_universality_probe.py",
  },
  {
    date: "2026-07-25",
    part: 75,
    title:
      "Doors A and B furnished: two no-go results and a λ* functional equation",
    verdict: "DOORS-FURNISHED",
    summary:
      "T73–T75: uniform cone route closed even in the continuum (sign-constancy lemma; window pin [0, δ_h)), yet per-direction hybrid cones exist for all 19 uncovered directions (HYBRID-GAINS). Door A (VACUUM-STRUCTURE): spectral Dirac phase carries the metaplectic sign datum (~93%) but has a Hecke defect; L2 no-go — the spectral world is provably sign-blind; Dirac vacuum sits on the minimal Waldspurger-mass atom; final R1–R5 (Krein quotient: ♭ as null sector, not summand). Door B (LAMBDA-STRUCTURED): λ* has closed form, own FE (orbit invariant tanh(σ²ω²/2) under multipliers ⋊ dilations), critical width σ_c = √2, convexity; target inequality measured both sides (safety 273 → < 1 at ω = 4.2). Sandbox; not RH evidence. Classical: Farkas, Dirac sea / Krein, Fejér / convex analysis.",
    badge: "sandbox",
    script: "door_b_lambda_transport_probe.py",
  },
  {
    date: "2026-07-25",
    part: 0,
    title:
      "Promoted: amplitude route and positive linear carrier (v540) — the open boundary λ* is part of the claim",
    verdict: "MACHINE-VERIFIED",
    summary:
      "The T67–T72 amplitude chain is now load-bearing (34 checks, ~3 s). (A) The amplitude Dirac D = [[0,V],[Vᵀ,0]] with D² = family kernel exists exactly and is Hecke-equivariant. (B) The geometric polarisation b = N₊ − N₋ is exact; Θ = N₊ + N₋ is a pure Siegel–Weil Eisenstein σ₃-eigenform with the Cohen seed Θ(d) = −48·L(−1,χ_d) (exact-rational, 159 live d). (C) Every coefficient bilinear form inherits the even-k deletion (Cauchy–Littlewood, five channels) — and the deletion is exactly the square-class double counting of the towers. (D) The positive linear carrier ℓ²(d, 48|L(−1,χ_d)|·|d|^{−a}) carries full weights and the plus balance Q = Q_ζ(g₋) + Q_ζ(g₊). (E) The completed FE Λ_Θ(s) = 8^{1−s}·Λ_Θ†(5/2−s) holds exactly (Fricke closed, rel ~1e-40); the guaranteed cone is FE-self-dual. HONEST BOUNDARY, inside the claim: the positivity is Euler-region positivity (edge L-values, NOT the central line); the residual distance to the Weil cone is the FE-covariant gap functional λ* on n ≡ 6 mod 8 — named, measured, Farkas-certified, and not removable by any finite signed theta library. NOT 'almost RH'; ZETA.HP.CARRIER untouched; no marker moves.",
    badge: "machine-verified",
    script: "v540_amplitude_linear_carrier.py",
  },
  {
    date: "2026-07-25",
    part: 72,
    title:
      "The cone library saturates: the gap becomes one functional",
    verdict: "COMPLEMENTARY-PAIR",
    summary:
      "T72 (COMPLEMENTARY-PAIR, 34/34): the cone library is surveyed — twist absorbs the sign class n ≡ 0,1 mod 4 (gap −26% mean, −90% max) but coverage saturates at 5/24; h(0) > 0 pins every Weil element against every twist. Final compression: the entire residual distance to the Weil cone is the FE-covariant gap functional λ* on the atoms n ≡ 6 mod 8 — no finite theta library can remove it. The transport problem, compressed to a single measurable object. Sandbox; this is not RH evidence; Euler-region positivity is not a central-line statement. Classical named: Farkas/LP duality, Weil 1952.",
    badge: "sandbox",
    script: "cone_enlargement_probe.py",
  },
  {
    date: "2026-07-25",
    part: 71,
    title: "A positive linear carrier with plus-only ζ-content",
    verdict: "LINEAR-PLUS-COMBINATION",
    summary:
      "T70–T71: the positive linear carrier stands — Θ(d) = −48·L(−1,χ_d) exact (Cohen 1975); full weights [1,1,1,1]; Weil balance Q = Q_ζ(g₋) + Q_ζ(g₊) plus-only (rel ~1e-15); ζ(2s) only as the squarefree sieve of the carrier. FE exact: Λ_Θ(s) = 8^{1−s}Λ_Θ†(5/2−s) (rel ~1e-40, Fricke closed); plus survives reflection; mirror family has a rigid sign law; guaranteed cone is FE-selfdual (Weil cone is not); overlap 5/24, violations exactly at the first spectral node. Sandbox; Euler-region positivity ≠ central-line claim. Classical: Cohen, Shintani, Jacobi/Fricke, Weil 1952.",
    badge: "sandbox",
    script: "linear_measure_weil_probe.py",
  },
  {
    date: "2026-07-25",
    part: 69,
    title: "The square level closes: every coefficient square deletes",
    verdict: "DELETION-UNIVERSAL",
    summary:
      "T67–T69: amplitude Dirac D = [[0,V],[Vᵀ,0]] with D² = family kernel is exact and Hecke-equivariant; signs of b(d) are genuine metaplectic residue (52% mixed fibres). Geometric polarisation b = N₊ − N₋ is exact; Θ = N₊+N₋ is pure Siegel–Weil Eisenstein (σ₃ eigenvalues, null cusp); b² = Θ² − 4N₊N₋ — but the minus is polarisation-invariant. Cauchy–Littlewood: every coefficient bilinear form inherits even-k deletion (theorem-like, five channels); the minus is exactly the square-class double-counting of the towers (inclusion–exclusion). Square level closed; amplitude/linear route is the remaining path. Sandbox; not RH evidence. Classical: Siegel–Weil, Cauchy–Littlewood, Rankin.",
    badge: "sandbox",
    script: "mixed_channel_full_weight_probe.py",
  },
  {
    date: "2026-07-25",
    part: 67,
    title:
      "Amplitude-level Dirac candidate: D² = family kernel, exact and Hecke-equivariant",
    verdict: "DIRAC-SQRT-EXACT",
    summary:
      "T67 (DIRAC-SQRT-EXACT, 27/27): the Dirac candidate D = [[0,V],[Vᵀ,0]] with D² = family kernel exists exactly and is Hecke-equivariant; amplitude towers carry degree-2 Euler factors without the doubling factor (the minus candidate is structurally absent at amplitude level). Signs of b(d) are a genuine metaplectic residue beyond all characters (52% mixed fibres, conditional entropy 0.73 bit — Kohnen depth, classical); the categorical seesaw is quantified (amplitude: no ζ-core, no minus; square: ζ-core + minus, via Clebsch–Gordan). Surprise: the phase lives in the d-side structure, not the m-side kernel. Open: the canonical polarisation that recovers the ζ-core with positive sign upon squaring. Sandbox; not RH evidence.",
    badge: "sandbox",
    script: "amplitude_dirac_sqrt_probe.py",
  },
  {
    date: "2026-07-25",
    part: 66,
    title:
      "π-front closed: digits are uniform noise; primes are arithmetic",
    verdict: "PI-FRONT-CLOSED",
    summary:
      "T65–T66 (PI-NULL-ARITHMETIC-DISTINCT 16/16; FOUR-LEVEL-NULL 23/23): π digits at prime places are uniform; density-detrended cross-correlations null; full placebo battery (e, √2, log 2, crypto, p±1) plus blind windows — no replication. Phase plane classical: the π spike is the continued-fraction convergent 355/113; compiler constant 1/(8π) is not an outlier (z = −0.90). Contrast: π-driven Cramér randomness reproduces prime density perfectly, but only true primes track Hardy–Littlewood pair correlation (corr +0.81 vs flat) — what makes primes primes is arithmetic, not randomness. Sandbox; π-digit route closed; not RH evidence.",
    badge: "sandbox",
    script: "pi_prime_four_level_probe.py",
  },
  {
    date: "2026-07-25",
    part: 0,
    title:
      "Program status: one categorical obstruction remains",
    verdict: "STATUS",
    summary:
      "Consolidated stand: TFPT generates an infinite stabilized arithmetic trace structure that agrees with the Weil framework up to one explicitly identified categorical obstruction — the passage from a quadratic automorphic family to its canonical metaplectic or linear square root. Progression: v539 [machine-verified] isolates two obstructions (minus + Corr); sandbox T64 resolves Corr as the exact det/det₂ Hilbert–Carleman Jacobian (~1.5e-16); one categorical minus remains. Not RH evidence; no almost-RH. (π-digit route later closed at T65–T66; amplitude Dirac at T67 — see newer feed entries.)",
    badge: "sandbox",
  },
  {
    date: "2026-07-25",
    part: 0,
    title: "Promoted: Weil structure of the compiler family (v539) — two isolated obstructions",
    verdict: "MACHINE-VERIFIED",
    summary:
      "The compiler family now has a fully identified Weil structure up to TWO EXPLICITLY ISOLATED OBSTRUCTIONS (obstructions are the claim, not a footnote). (A) GNS = direct integral over the Gelfand spectrum: fibre twist-mix ≤5.037e-16. (B) Trivial Sato–Tate isotype = GL(1) core G₀=(1+Y)/(1−Y)=ζ_p(w−3)²/ζ_p(2w−6)=Σ 2^ω(n) n^{−u}. (C) Q_fam = 2Q_ζ(g) − 2Q_ζ(g♭) + Arch + Corr (Prime_F rel ≤9e-16). Obstruction 1: doubling enters with a MINUS — family positivity does NOT imply Q_ζ≥0. Obstruction 2: Corr is non-automorphic (e^{−Σ p^{−u}}-type). Finite-class Q_fam∈[4.369,11.486] measured; dense-class positivity open/RH-adjacent — not claimed. NOT almost-RH; ZETA.HP.CARRIER untouched; no marker moves.",
    badge: "machine-verified",
    script: "v539_weil_structure_family.py",
  },
  {
    date: "2026-07-25",
    part: 0,
    title: "Promoted: the compiler relative-trace identity (v538)",
    verdict: "MACHINE-VERIFIED",
    summary:
      "v535, v536 and v537 are now load-bearing as three projections of one finite relative-trace identity: geometric Tr_V(nu_p) from the glue census equals spectral Tr_V from the eta product for every (p, projector) at p = 3,5,7 (anchors e.g. 727680 full / 19840 on f8 at p = 7); the Eichler split is the RTF orbit dictionary (lambda_Eis + a_p^2); and the Waldspurger period side is literally quaternary lattice counting with R = 23.1873585645 (rel ~1e-12). Verdict ONE-FORMULA. Honest fence: the identity is finite — the infinite RTF (all p,d + archimedean term) stays open; classical Jacquet/Eichler/Waldspurger named classical; no RH claim.",
    badge: "machine-verified",
    script: "v538_relative_trace_identity.py",
  },
  {
    date: "2026-07-25",
    part: 64,
    title: "The final relation and its two obstructions",
    verdict: "HALF-STABLE",
    summary:
      "T63–T64: the global core identity ζ(u)²/ζ(2u) = Σ 2^ω(n) n^{-u} is exact, and the linear relation Q_fam = 2Q_ζ − 2Q_ζ(♭) + Arch + Corr closes to machine precision. Obstruction 2 is resolved in sandbox: the extra term is exactly the det/det₂ Hilbert–Carleman Jacobian (rel 1.5e-16) — 'the determinant of a forgotten compact operator'; det₂ is necessary on the critical line. Obstruction 1 remains fundamental: the minus sign (32/32 endoscopy characters; Krein analysis) is the square/sym² signature of the family. Not RH evidence; the open question is categorical — what is the metaplectic square root of the family?",
    badge: "sandbox",
    script: "rtf_stabilization_probe.py",
  },
  {
    date: "2026-07-25",
    part: 62,
    title: "From GNS to Weil: the gap gets a closed form",
    verdict: "TRANSFORM-REQUIRED",
    summary:
      "T59–T62: the gap between the positive GNS form and the classical Weil quadratic (Weil 1952) has closed form — SU(2) characters Φ_k at Satake angles. Prime averages contract Φ_k onto Sato–Tate moments (1,0,2,0); k=1 unconditional via Rankin–Selberg. The trivial isotype is exactly the GL(1) core G₀ = ζ_p(w−3)²/ζ_p(2w−6), and fibre decomposition kills twist-mix (6e-16 = null). Sandbox structure map — not a positivity proof and not RH evidence.",
    badge: "sandbox",
    script: "weil_gns_identification_probe.py",
  },
  {
    date: "2026-07-25",
    part: 58,
    title: "Both infinities packed — and no new xi-sector",
    verdict: "PACKED",
    summary:
      "T57–T58: both infinities (p-towers with closed Euler factors + d-family) sit in one bilaterally verified double series Z(s,w) (errors down to 8e-8; classical Goldfeld–Hoffstein multiple Dirichlet series named classical). Residual/diagonal/Weyl/Eisenstein/GNS probes find only classical GL(1) shadows — Eisenstein floor is an exact ζ-product; no new ξ-sector. Route A (degeneration → ξ) closed at packing level. This is not RH evidence.",
    badge: "sandbox",
    script: "infinite_rtf_packing_probe.py",
  },
  {
    date: "2026-07-25",
    part: 56,
    title: "The canonical measure question — GNS builds its own Hilbert space",
    verdict: "MARGINAL-WEIGHT",
    summary:
      "T54–T56: |d|^{-5/2} is the critical line of the family measure, not a canonical measure (kills fired); the RTF forces weight |d|^{-1}. The positive pairing at that weight is cutoff-form-independent, Hecke-equivariant, and generates ℓ²(d, b²/|d|) by GNS. Sharpening to D=60000 / 6032 live d closes the factor chain to 0.09% — the old σ_eff 'tension' was a category error (slope vs density). Sandbox; central values at GL(2) centre 2, not the ξ-line.",
    badge: "sandbox",
    script: "rtf_period_pairing_probe.py",
  },
  {
    date: "2026-07-25",
    part: 53,
    title: "The reframe: family + trace formula instead of one operator",
    verdict: "INFINITE-CARRIER",
    summary:
      "T51–T53: the Waldspurger family kernel K_D is the series' first non-collapsing infinite carrier (rank grows 8→192; spectrum = central-value family). The Hecke path space realises L-functions as fibre Fredholm determinants — global det ≡ 1 honestly discarded. Review conjecture confirmed: v535/v536/v537 are three projections of one finite relative-trace identity (later promoted as v538); g is the quaternary lattice form n=(x²+y²)/2+2z²+u²+2w². Classical Waldspurger/Jacquet named classical; no RH claim.",
    badge: "sandbox",
    script: "waldspurger_family_kernel_probe.py",
  },
  {
    date: "2026-07-25",
    part: 50,
    title: "The central-value family aggregates exactly — but away from the xi-line",
    verdict: "FAMILY-AGGREGATED",
    summary:
      "The self-convolution of the bridge object g factorises exactly through the Shimura relation (10000/10000 coefficients), the Waldspurger constant extends to 46 discriminants with spread 7.7e-16, and the aggregate identity D_gg = R * sum |d|^{3/2-s} L(f8 x chi_d, 2) E_d(s) + non-fundamental part closes to 1e-12. Honest boundary: the aggregate lives at abscissa 5/2 — farther from the xi-centre 1/2 than the pointwise values; aggregation works but does not move toward the centre.",
    badge: "sandbox",
    script: "half_integral_selfconvolution_probe.py",
  },
  {
    date: "2026-07-25",
    part: 49,
    title:
      "The only mechanism that ever proved an RH now runs inside the suite — with one cell vacant",
    verdict: "MINIATURE-RUNS",
    summary:
      "The classical Rankin/Landau positivity argument is rebuilt as a machine-verified inference: nonnegativity plus the measured abscissa 4.0026 proves |a_p| = O(p^1.55) without checking values individually, and the same trick gives alpha_tau = 11.9997 after confirming Delta = (theta2 theta3 theta4 / 2)^8 lies in the compiler theta monoid (the Ramanujan tau world included). The strategic map: on the eigenvalue side all three mechanism roles (family, positivity, pole) are occupied in-suite; on the zeros side over Z exactly one cell is vacant — the family role (the classical F1 problem). What is missing for RH is now an empty table cell next to demonstrably working ones; the miniature proves eigenvalue magnitudes, not zero locations.",
    badge: "sandbox",
    script: "rankin_positivity_miniature_probe.py",
  },
  {
    date: "2026-07-25",
    part: 48,
    title: "The scale-uniqueness kill test: route 1's explanatory claim is dead",
    verdict: "SCALE-TORSOR",
    summary:
      "A zero-free kill test (no Riemann zeros loaded, enforced by the probe itself) enumerated all compiler-native scale candidates by a mechanical generation rule and checked every verified selection principle: Poisson self-duality selects a scale only under a free root-norm convention (the {pi, 2pi} torsor remains), the unimodularity gate is scale-blind, and KMS/BW was already deflated. A residual two-element family means the compiler does not fix the prolate scale before looking at zeros — the explanatory claim of stage-4 route 1 is killed, and the construction survives only as external mathematics. All three stage-4 routes are now either killed-as-explanation or honestly typed as long-term theory projects; the series rests fully mapped.",
    badge: "sandbox",
    script: "lambda_scale_uniqueness_probe.py",
  },
  {
    date: "2026-07-25",
    part: 47,
    title:
      "UV/Dirac prototype: no pointwise signal — and preregistration caught a near-false-positive",
    verdict: "UV-NO-SIGNAL",
    summary:
      "The numerically clean track shows no pointwise Dirac-to-zeros match out of sample (mean rel 0.34, null controls comparable), while the exterior-coupling track looked spectacular (2.5% train, 4% out-of-sample) but failed the preregistered convergence kill — a near-false-positive correctly killed by the probe's own discipline. Connes-Moscovici 2022 remains an external asymptotic counting anchor, not a few-percent prototype; the stage-4 self-duality class stays R1-certified with the remaining levers being genuine theory projects, not quick probes. The prime-front series (parts 11-47) rests here at its natural saturation point: three promoted modules, a three-channel functor map, and an honestly mapped stage-4 frontier.",
    badge: "sandbox",
    script: "stage4_prolate_uv_dirac_probe.py",
  },
  {
    date: "2026-07-25",
    part: 46,
    title:
      "Stage-4 prototype: the self-duality operator class passes R1, the IR shortcut does not",
    verdict: "IMPL-OPEN",
    summary:
      "The prolate spheroidal operator — the canonical object of Fourier self-duality (Slepian, classical; Connes-Moscovici 2022 as external anchor) — is implemented and validated against literature eigenvalues to 4.5e-4, giving the stage-4 map its first R1-certified unbounded self-adjoint prototype. The naive IR assignment |mu_n| ~ gamma_n^2 + 1/4 fails preregistered criteria (~90% off, null controls confirm no signal there): the published statement concerns the UV counting / Dirac square root, not interior IR eigenvalues. Path not dead — the next probe targets the UV/Dirac construction; no TFPT result and no RH progress is claimed either way.",
    badge: "sandbox",
    script: "stage4_prolate_prototype_probe.py",
  },
  {
    date: "2026-07-24",
    part: 0,
    title: "Promoted: the half-integral bridge (v537)",
    verdict: "MACHINE-VERIFIED",
    summary:
      "The half-integral bridge is now load-bearing: v537 verifies in one module that the compiler theta monoid contains a unique weight-5/2 Shimura preimage of f8 (Sh(g) = -8 f8 to n = 120), that g is an exact T(p^2)-eigenform with eigenvalues a_p (p = 3..13), that the Kohnen-plus route is structurally out of scope (an honest negative built into the claim), and that the Waldspurger quotient is constant to thirteen decimal places across ten discriminants (R = 23.1873585645..., d = 5 mod 8 vanishing structurally). Audit OK, Wolfram mirrors 603/603. Classical theorems named as classical; no RH claim.",
    badge: "machine-verified",
    script: "v537_halfintegral_bridge.py",
  },
  {
    date: "2026-07-24",
    part: 45,
    title: "Central values wired: Waldspurger constancy to thirteen decimal places",
    verdict: "CENTRAL-VALUES-WIRED",
    summary:
      "With functional-equation-based twisted L-values (the part-44 outlier d = 89 was purely numerical and now sits exactly on the constant), the Waldspurger quotient R(d) = b(d)^2 / (d^{3/2} L(f8 x chi_d, 2)) is constant to 2.7e-13 across ten fundamental discriminants d = 1 mod 8, while the d = 5 mod 8 class vanishes structurally via root number -1. The central L-value family is realised as the coefficient squares of the compiler object g — the third functor channel is wired (classical Waldspurger/Baruch-Mao; GL(2) centre s = 2, not the xi-line; no RH content).",
    badge: "sandbox",
    script: "baruch_mao_metaplectic_probe.py",
  },
  {
    date: "2026-07-24",
    part: 44,
    title: "The Kohnen-plus route is structurally closed",
    verdict: "OUTSIDE-KOHNEN-SCOPE",
    summary:
      "Three independent reasons close the central-value route via the Kohnen plus space: the compiler object g has exactly zero mass on n = 0 mod 4 (so the spectral projector sees no plus component), the support-cut is not a Shimura-Kohnen preimage, and the level 32 = 4*8 with even M = 8 lies outside Kohnen's 1982 theorem. Exploratory Waldspurger quotients nevertheless cluster around ~21-23 with one badly-converged outlier — motivating one final probe of this line via the Baruch-Mao metaplectic formulas for non-plus forms (classical).",
    badge: "sandbox",
    script: "kohnen_operator_projection_probe.py",
  },
  {
    date: "2026-07-24",
    part: 0,
    title:
      "Promoted: Hecke from Geometry (v535) and the Eichler trace layer (v536)",
    verdict: "MACHINE-VERIFIED",
    summary:
      "Two modules from this series are now part of the load-bearing verification suite: v535 (Kneser correspondence carries the census Hecke action; the frozen neighbour operator outputs the a_p; census redundancy is purely 2-adic oldform structure) and v536 (the Eichler trace layer: closed Witt densities, the two-sided identity lambda_geom = lambda_Eis + a_p^2, and the signed O(p^3) extraction of a_p). Full ledger/paper/website sync, bash build.sh audit -> AUDIT OK, Wolfram mirrors 598/598. No RH claim; the weight-4 -> GL(1) transport stays an open research contract.",
    badge: "machine-verified",
    script: "v535_hecke_from_geometry.py",
  },
  {
    date: "2026-07-24",
    part: 43,
    title: "Kohnen plus space: eigenform found, isomorphism preimage not yet",
    verdict: "PARTIAL",
    summary:
      "A plus-space object with the correct T(p^2) eigenvalues (-4, -2, 24) does live in the q^4-extended compiler monoid, but it is not the Kohnen-isomorphism preimage of f8 (its Shimura lift vanishes at t = 2), so the Waldspurger central-value quotient scatters (~294%) instead of being constant. The central-value wire stays open; the named next lever is a true modular plus projection (operator, not a support cut) or the metaplectic extension. Classical ingredients (Kohnen 1980/82, Waldspurger) named as classical.",
    badge: "sandbox",
    script: "kohnen_plus_waldspurger_probe.py",
  },
  {
    date: "2026-07-24",
    part: 42,
    title: "Local densities closed: the Eichler layer is promotion-ready",
    verdict: "DENSITIES-CLOSED",
    summary:
      "The separating invariant is found — Type-A isotropic cosets are exactly those representing norm p, and the closed densities N_A = min(240(1+p^3), #iso-1) predict p = 7 correctly (82560/743040, verified live). The O(1) assembler lambda_geom = lambda_Eis + a_p^2 now runs to p = 31 with no enumeration, and the checklist for promoting the Eichler trace layer (candidate v536) is fully green. Classical ingredients (Witt, Eichler, Siegel densities) named as classical; the claim is the two-sided in-suite identity.",
    badge: "sandbox",
    script: "e8_local_densities_probe.py",
  },
  {
    date: "2026-07-24",
    part: 41,
    title: "The bridge is Hecke-equivariant: g carries the prime fingerprints as eigenvalues",
    verdict: "FUNCTOR-WIRED",
    summary:
      "The weight-5/2 compiler object g from part 38 is an exact T(p^2)-eigenform with eigenvalues precisely a_p(f8) for p = 3, 5, 7, 11, 13 (-4, -2, 24, -44, 22) — the Shimura bridge is a mechanism, not a coefficient coincidence (classical Shimura correspondence). A trial functional equation completes the half-integral side around centre 5/4. Honest gap: g is not in the Kohnen plus space, so the Waldspurger central-value window fails for now — the Kohnen projection is the named next lever.",
    badge: "sandbox",
    script: "functor_wiring_halfintegral_probe.py",
  },
  {
    date: "2026-07-24",
    part: 38,
    title: "The cuspidal bridge exists: a Shimura preimage inside the compiler",
    verdict: "BRIDGE-FOUND",
    summary:
      "Among 70 compiler-native theta candidates, exactly one weight-5/2 combination g = theta2(q^2)^2 * theta3(q^2) * theta4 * theta4(q^2) Shimura-lifts to -8*f8, verified coefficient-by-coefficient to n = 120 — the prime fingerprints a_p are now reachable from half-integral weight, the world where the zeta centre 1/2 lives (classical Shimura correspondence; honest caveat: half-integral L-series generally lack Euler products). The cuspidal self-convolution also places zeta exactly ON the line after weight normalisation (offset 0, vs offset ±1 for the abelian route).",
    badge: "sandbox",
    script: "cuspidal_bridge_probe.py",
  },
  {
    date: "2026-07-24",
    part: 37,
    title: "Signed and scalable: the prime fingerprint from small shells",
    verdict: "SCALABLE-AND-SIGNED",
    summary:
      "The naive invariant classification died at p = 5 (preregistered kill) and was repaired via the Shell(p)-image refinement; the assembler identity R(p) = σ₃ − 1 − c(p²)/8 reproduces the cusp residues 16 and 4 live and predicts R(7) = 576. The sign of a_p now comes from Shell(p) (order p³ work instead of p⁷): a_p = −c(p)/8 gives −4, −2, +24 with correct signs. A single Salié sum is not a_p (kill fired); promotion needs closed local densities first.",
    badge: "sandbox",
    script: "e8_salie_signed_cusp_probe.py",
  },
  {
    date: "2026-07-24",
    part: 40,
    title: "Stage-4 terrain map",
    verdict: "TERRAIN-MAPPED",
    summary:
      "Today's operator algebra has a two-point spectrum; a Hilbert–Pólya carrier needs infinitely many. Only two candidate classes remain (seam modular flow; adelic BC completion), each with preregistered kills. Distance to RH stated honestly: large.",
    badge: "sandbox",
    script: "weil_recovery_compiler_sector_probe.py",
  },
  {
    date: "2026-07-24",
    part: 39,
    title: "Centre atlas — abelian weight drop closed",
    verdict: "CLOSED",
    summary:
      "The ξ-line (centre 1/2) is reached exactly by weight ≤ 1 theta factors: Mellin(θ₃) → ζ(2s) and ζ_{ℚ(i)} = ζ(s)L(s,χ₄). Abelian drop = factorisation + Mellin, typed closed. Cuspidal channel remains at centre 2 — no spectral / RH content claimed.",
    badge: "sandbox",
    script: "zeta_center_atlas_probe.py",
  },
  {
    date: "2026-07-24",
    part: 36,
    title: "λ_geom evaluator — two-sided Eichler to p ≤ 5",
    verdict: "HARDENED",
    summary:
      "Independent f₈-free geometry: λ_geom = λ_i + λ_ii hits anchors λ_geom(3)=352 and λ_geom(5)=3784; residual R(p)=a_p² exactly. Two-sided confirmation (mod-p geometry ↔ eta product) at p ≤ 5; p ≥ 7 gated by shell size.",
    badge: "sandbox",
    script: "e8_lambda_charsum_evaluator_probe.py",
  },
  {
    date: "2026-07-24",
    part: 35,
    title: "Rankin–Selberg translates only the abelian shadow",
    verdict: "TYPED",
    summary:
      "RS-TRANSLATES-EISENSTEIN: multiplicative convolution of compiler factors yields exact ζ(s)L(s,χ₄)ζ(s−2)L(s−2,χ₄) products — but after weight normalisation no GL(1) factor sits on the ξ-line, and f₈ is invisible. Cuspidal channel does not cross this bridge.",
    badge: "sandbox",
    script: "rankin_selberg_functor_probe.py",
  },
  {
    date: "2026-07-24",
    part: 34,
    title: "Compiler functor founded with anchor",
    verdict: "FOUNDED",
    summary:
      "In-suite RH centre 1/2 is defined via Mellin((θ₃−1)/2); unique factorisation among 45 monomials picks only (2,6,0) for the signed glue theta. Bare Shimura on θ₃ sees a_p not at all — contract ZETA.COMPILER.FUNCTOR [O] founded with kills K1–K3.",
    badge: "sandbox",
    script: "compiler_functor_theta_probe.py",
  },
  {
    date: "2026-07-24",
    part: 33,
    title: "Many-prime hardening — Eichler identity",
    verdict: "HARDENED",
    summary:
      "HARDENED-TO-p≤100: geometric count splits as elementary Witt part + exactly a_p². Closed forms for all 24 odd primes ≤ 100; Ramanujan bound and cusp rule a_p = b − σ₃ hold. Smooth background + coherent interference.",
    badge: "sandbox",
    script: "e8_nu_rule_many_primes_probe.py",
  },
  {
    date: "2026-07-24",
    part: 32,
    title: "PROMOTED — Hecke from geometry (v535)",
    verdict: "PROMOTED",
    summary:
      "T27+T31+T32 consolidated into verification/v535_hecke_from_geometry.py (HECKE.GEOM.01, 25/25, AUDIT OK). The one load-bearing result of this series: Kneser + affine ν_p + 2-adic oldforms. No RH claim; weight-4→GL(1) stays open.",
    badge: "machine-verified",
    script: "v535_hecke_from_geometry.py",
  },
  {
    date: "2026-07-24",
    part: 32,
    title: "Census redundancy is 2-adic oldform structure",
    verdict: "EXACT",
    summary:
      "REDUNDANCY-IS-OLDFORM: dim V = 7 = 5 (E₄ copies, levels 1..16) + 2 (f₈ copies, levels 8,16). Recovery projectors π_cusp = (28−T₃)/32, π_Eis = (T₃+4)/32; new spaces each dim 1. Level census purely 2-adic.",
    badge: "sandbox",
    script: "census_newform_recovery_probe.py",
  },
  {
    date: "2026-07-24",
    part: 31,
    title: "Neighbour operator = affine Hecke element",
    verdict: "REPAIRED",
    summary:
      "REPAIRED-BY-DICTIONARY: ν_p = a·Id + b·T_p with b = σ₃ + a_p — one (a,b) for all shells/channels. Frozen geometry outputs a₃=−4, a₅=−2, |a₇|=24. T30's broken reading was a normalisation artifact.",
    badge: "sandbox",
    script: "e8_neighbor_operator_decomposition_probe.py",
  },
  {
    date: "2026-07-24",
    part: 30,
    title: "Uniformity test broken — p=3 coincidence",
    verdict: "DEAD",
    summary:
      "The Teil-28 reading a_p = #P³ − λ_odd/8 fails at p=5 (λ_odd=3784 ≠ 1264). Transport structure itself holds; only the eigenvalue reading was a coincidence. Named repair path → T31.",
    badge: "sandbox",
    script: "e8_transport_cusp_uniformity_probe.py",
  },
  {
    date: "2026-07-24",
    part: 29,
    title: "Primon transport category — relations verified",
    verdict: "EXACT",
    summary:
      "RELATIONS-VERIFIED: coprimality, ladder, and constructive decompositions of T_n hold on compiler channels. Stable 4×4 on (Tot, Spinor, Eis, f₈): T₃ = diag(28,28,28,−4), T₅ = diag(126,126,126,−2).",
    badge: "sandbox",
    script: "primon_transport_category_probe.py",
  },
  {
    date: "2026-07-24",
    part: 28,
    title: "Cusp coefficient visible from Kneser transport",
    verdict: "PARTIAL",
    summary:
      "At p=3, full enumeration of 1120 neighbour lattices yields a₃ = −4 exactly from mod-3 geometry — no modular-form import. Honest fence: uniformity untested (later broken at T30).",
    badge: "sandbox",
    script: "e8_kneser_transport_holonomy_probe.py",
  },
  {
    date: "2026-07-24",
    part: 27,
    title: "Kneser p-neighbours carry Hecke (Eisenstein half)",
    verdict: "PARTIAL",
    summary:
      "#iso_lines = σ₃(p)·#P³ enumerated at p=2,3,5,7 (135/1120/19656/137600). Eisenstein eigenvalue is lattice-native; μ4-refined line counts do not yet recover a_p(f₈). Later consolidated into v535.",
    badge: "sandbox",
    script: "e8_pneighbor_hecke_probe.py",
  },
  {
    date: "2026-07-24",
    part: 26,
    title: "Shell quotienting dead",
    verdict: "DEAD",
    summary:
      "Weyl-orbit and W(D5)×W(A3)×glue quotients break multiplicativity. Only classical primitive stratification β = r_prim/240 survives — still weight 4, no GL(1). Full W(E8) does not preserve μ4 glue class.",
    badge: "sandbox",
    script: "shell_multiplicity_one_quotient_probe.py",
  },
  {
    date: "2026-07-24",
    part: 25,
    title: "Scattering phase dead — arch route closed",
    verdict: "DEAD",
    summary:
      "Last non-counting observable fails: seam overlaps lack the thermal Unruh magnitude; discrete sector clock (π/4, π/2) remains. Arch-from-seam route closed in four preregistered stages. Seam = discrete μ4 clock, not a hidden Gamma.",
    badge: "sandbox",
    script: "seam_tate_phase_probe.py",
  },
  {
    date: "2026-07-24",
    part: 24,
    title: "Truncated-trace / RvM dictionary dead",
    verdict: "DEAD",
    summary:
      "Under the only allowed self-dual rule L=√T, Peschel counting leads 1/(2π²) not 1/(2π). The fit that would hit RvM is algebraically forbidden. π-factor misses in every allowed counting version.",
    badge: "sandbox",
    script: "seam_truncated_trace_rvm_probe.py",
  },
  {
    date: "2026-07-24",
    part: 23,
    title: "Boost–arch dictionary dead on three lines",
    verdict: "DEAD",
    summary:
      "Chord-vs-arc 2/π hypothesis fails; coupled cutoffs miss; digamma-in-the-seam finds no named ψ/log candidate. Seam carries the boost, not the archimedean slot of the explicit formula.",
    badge: "sandbox",
    script: "seam_boost_arch_dictionary_probe.py",
  },
  {
    date: "2026-07-24",
    part: 22,
    title: "Interval cut carries the missing log — dictionary still fails",
    verdict: "PARTIAL",
    summary:
      "Entanglement Hamiltonian of a seam half-cut has 1/ln L spacing and Peschel universal form (K1/K2). One-constant dictionary ln L ↔ ln(T/2π) misses by classical factor 2/π (K3).",
    badge: "sandbox",
    script: "seam_halfcut_dilation_arch_probe.py",
  },
  {
    date: "2026-07-24",
    part: 21,
    title: "Can TFPT predict primes? Three honest channels",
    verdict: "TYPED",
    summary:
      "Exact geometric primality: n prime ⟺ shell(2n) has 240(1+n³) vectors (0 errors to 10⁴). Glue characters predict two-squares type at 100%. Positional prediction needs the zero spectrum — that operator does not exist.",
    badge: "sandbox",
    script: "prime_prediction_mechanics_probe.py",
  },
  {
    date: "2026-07-24",
    part: 20,
    title: "Arch kernel killed as naive",
    verdict: "KILLED-AS-NAIVE",
    summary:
      "Raw seam mode density is flat/falling O(1); classical digamma arch kernel grows as log. Archimedean term does not live in free seam DOS. Route narrowed, not absolutely dead — until T25 closed it.",
    badge: "sandbox",
    script: "seam_archimedean_kernel_probe.py",
  },
  {
    date: "2026-07-24",
    part: 19,
    title: "Prime thinning ↔ seam truncation",
    verdict: "MIXED",
    summary:
      "At suite temperatures the small prime tower carries almost all of log ζ(β) ({2,3,5}: 99.962% at β=2π). Post-hoc 'exactly {2,3,5} is the 99% set' is dead. TFPT tones are not distinguished on the prime clock.",
    badge: "sandbox",
    script: "prime_thinning_seam_truncation_probe.py",
  },
  {
    date: "2026-07-24",
    part: 18,
    title: "Guinand–Weil benchmark standing; arch gap typed",
    verdict: "BENCHMARK",
    summary:
      "Classical Guinand–Weil balance verified as executable reference (500 zeros). Weil↔v221 dictionary for pole+prime; archimedean digamma term has no seam counterpart — the precise gap of ZETA.WEIL.RECOVERY.",
    badge: "sandbox",
    script: "zeta_weil_recovery_trace_probe.py",
  },
  {
    date: "2026-07-24",
    part: 17,
    title: "Functional-equation duality on the census",
    verdict: "PARTIAL",
    summary:
      "Poisson involution Λ_E8(s)=Λ_E8(4−s) to ≥20 digits; character channels form an exact dual pair with factor 1. Weight transport 4→½ remains open with kill preregistered.",
    badge: "sandbox",
    script: "zeta_funceq_duality_probe.py",
  },
  {
    date: "2026-07-24",
    part: 16,
    title: "One Hecke polynomial covers three character channels",
    verdict: "PARTIAL",
    summary:
      "Weight-4 Hecke polynomial 1−a_p x+p³x² covers trivial/μ4/μ3 channels; Eisenstein factors as (1−x)(1−p³x); f₈ Deligne-pure. Seam derivation of the Hecke operation and GL(1) transport stay open.",
    badge: "sandbox",
    script: "zeta_local_euler_probe.py",
  },
  {
    date: "2026-07-24",
    part: 15,
    title: "Arithmetic completion — lattice Hamiltonian",
    verdict: "PARTIAL",
    summary:
      "E8 shell Hamiltonian H|x⟩=log(|x|²/2) yields Tr e^{−sH}=240ζ(s)ζ(s−3) — all primes, but multiplicity 240σ₃(n)≠1. Naive all-N completion killed; contract narrowed to multiplicity-1 reduction.",
    badge: "sandbox",
    script: "zeta_arith_completion_probe.py",
  },
  {
    date: "2026-07-24",
    part: 14,
    title: "Self-dual width slogan deflated",
    verdict: "DEFLATED",
    summary:
      "v526 measures β_steps=N (kernel detailed balance); angle 2π is universal BW/Unruh conversion — 'measured 2π is the self-dual width' overreached. Only E8 among atlas lattices admits a Poisson-self-dual width. Kill chain begins here.",
    badge: "sandbox",
    script: "seam_selfdual_width_probe.py",
  },
  {
    date: "2026-07-24",
    part: 13,
    title: "Glue character atlas — μ4 is not alone",
    verdict: "EXACT",
    summary:
      "Nine maximal root sub-splits yield cyclic glue orders {2,3,4,5,6}. μ3 mirrors the Teil-11 pattern (Q(ω), η(3τ)⁸); μ5 is χ₅-blind on that resolution. Character content is glue-order-limited and split-selective.",
    badge: "sandbox",
    script: "e8_glue_character_atlas_probe.py",
  },
  {
    date: "2026-07-24",
    part: 12,
    title: "Surprise bridges — Apéry / ζ(3) and L-series",
    verdict: "EXACT",
    summary:
      "Cusp form f₈=η(2τ)⁴η(4τ)⁴ is the Beukers/Ahlgren–Ono form: Apéry numbers A((p−1)/2) ≡ a_p mod p² for all odd p≤97. Signed channel is the unique pole-free L-series. Self-dual-width slogan later deflated at T14.",
    badge: "sandbox",
    script: "e8_glue_lseries_selfdual_probe.py",
  },
  {
    date: "2026-07-24",
    part: 11,
    title: "Signed census — θ₃²·θ₄⁶ is a tensor factor",
    verdict: "EXACT",
    summary:
      "Counting E8 shell points by glue colour and taking the signed difference yields exactly θ₃²·θ₄⁶. The Gaussian-integer theta (π via L(1,χ₄)=π/4) is a literal tensor factor of the μ4-glue census. Three channels: total, signed, spinor.",
    badge: "sandbox",
    script: "e8_glue_chi4_signed_theta_probe.py",
  },
];

/** Display order for the narrative sections (anchor ids). */
export const PRIME_FRONT_SECTIONS = [
  { id: "hook", label: "Hook" },
  { id: "compiler", label: "Compiler" },
  { id: "census", label: "Signed census" },
  { id: "bridges", label: "Surprise bridges" },
  { id: "kill-chain", label: "Kill chain" },
  { id: "predict", label: "Prime prediction" },
  { id: "hecke", label: "Hecke from geometry" },
  { id: "eichler", label: "Eichler layer" },
  { id: "weight-drop", label: "Weight drop" },
  { id: "stage-4", label: "Stage-4 map" },
  { id: "july-25-arc", label: "July 25 arc" },
  { id: "program-status", label: "Where it stands" },
  { id: "amplitude-route", label: "Amplitude route" },
  { id: "doors-furnished", label: "Doors furnished" },
  { id: "three-perspectives", label: "Three perspectives" },
  { id: "i5-geography", label: "I5 geography" },
  { id: "mechanism", label: "The mechanism" },
  { id: "instrument-race", label: "The race" },
  { id: "two-doors", label: "Two doors" },
  { id: "meaning", label: "What it would mean" },
  { id: "updates", label: "Live updates" },
] as const;
