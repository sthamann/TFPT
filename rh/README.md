# rh/ — the consolidated RH-program workspace

> **Claim boundary (non-negotiable).** This workspace is **research
> documentation**. Nothing in it is evidence **for or against the Riemann
> Hypothesis, in either direction. NO RH CLAIM.** Status markers follow the
> repository ledger (`verification/status_ledger.csv`, which wins on any
> disagreement): **[E]** = certified finite theorems (machine-checked exact
> identities), **[O]** = honestly open. Mincut **base 4 / refined 5** —
> unchanged by anything in this tree.

`rh/` is the single home for the TFPT prime-front RH program: the certified
finite reduction chain, a focused paper, a Lean 4 formalization pilot, and an
RH-specific verification suite. It **references** the canonical sources in the
main tree (verification modules, sealed discovery probes, ledger rows) with
pinned SHA-256 hashes — it **never duplicates** them.

## What the RH program is

The program studies the *window object*: the finite Galerkin compression of
the localized Weil operator onto a window of the von-Mangoldt double comb
(the "MAIN world"), together with scramble / Epstein / smooth control worlds.
Positivity of this finite object through the whole window — "the wall" — is
a *finite, machine-decidable* statement per window. Two years of rounds have
compressed the entire question into an exact dictionary (tau functions,
bordered Riemann–Hilbert problems, a coupled two-term recursion) whose
**certified finite part is [E]** and whose **remaining arithmetic edge is two
precisely named open statements [O]**.

## The reduction chain (current state)

| Stage | Statement | Status | Module / round |
|---|---|---|---|
| 1 | Window object: Galerkin of the localized Weil operator; wall = quasi-definiteness of the signed defect measure `mu - nu` (`WALL_EQUIVALENT`) | [E] | `v955` (r224–r226) |
| 2 | Half-filling law `N_w = ceil(S/2)`: MAIN survives its entire free moment window, controls die at 0.11–0.15 `N_w`; four exact source-pure coordinates; dual-FIK L-gauge | [E] | `v956` (r228–r231) |
| 3 | Bordered-RHP dictionary: bordered PSD = wall + budget; Schlesinger rank-1 insertion; centering congruence; three exact error formulas; contour `R_1`; PSD base theorem | [E] | `v958` (r244–r253) |
| 4 | Coupled-tau recursion `tau_{n+1} = a_n tau_n`, `tau^aug_{n+1} = a_n tau^aug_n + b_n tau_n` ⇒ Riccati drain `D_{n+1} = D_n - rho_n`; bilinear identity; Schur budget telescope `D_N = B - sum rho`; terminal reduction `margin > 0 ⟺ q_N < 1` | [E] | `v959` (r256–r259) |
| 5 | Terminal-surface closure: two-branch theorem (`Z_w = r_{N-1}`, cheap branch 35/42 + both mains **without cancellation**, exception set 7 rungs) + sealed phase bounds (edge truncation, block alternation) + kz15 interval certificate (dps 640, width 1.5e-92, margin +0.0268) + universal pair theorem H1–H4 (Lean-proved) — **the 42-rung census `q_N < 1` certified on all 42 rungs** (41 mechanism target-blind + kz15 exact-finite); census + certificate, NOT cofinal | [E] | `v960` (r260–r275) |
| 6 | Midpoint-orientation dictionary: base Casoratian `= h_n` (`c' = 1`), h-free midpoint form = node polynomial, augmented telescope `D_{n+1} = B − W^aug/W^base`; Maslov census GO (rule R2 = Jacobi interlacing/reality, blind 42/42, controls at flip+1, one-way; raw atom-Sturm census REFUTED as winding quantity); metric-firewall measurements (graded continuum, exact Hellmann–Feynman gradients, u-profile predictive, stability law perturbative-only) | [E] | `v961` (r274–r278) |
| 7 | **Half-filling pinning theory** (wave 5): T1 moment counting (free pivots exactly `h_0..h_{N_w−1}` — half-filling is the end of the free moment space, "why half-filling" answered by counting) + T2 crossing budget (`#(h<0) = S_−`, Jacobi/Sylvester, world-blind) + T3 two-sided parity (h-blind; census bilanz + 87376-case exhaustion) + T4 main window reduction (**the entire open statement is `minC ≥ N_w` ⟺ `∀ n < N_w : h_n > 0`**) — plus the four refutations (no universal O(1) pinning, no extremality, no generic Maslov obstruction, no simple offset law) | [E] | `v962` (r279–r281) |
| 8 | **The L* reduction dictionary** (wave 6): the r283 A2 chain exact (mu-frame congruence `minor_k(D_mu − G) == D_k(mutilde)`; contraction `h > 0` through the window ⟺ `lambda_max(E_m) < 1`; crossing exactly at `minC + 1`; pigeonhole `minC ≤ S_+`; capacity-as-counting refuted by the rank-one pair) ⇒ **L\* is the canonical reduction** (`∫p²dnu < ∫p²dmu` for `deg p < N_w`); the r285 bookkeeping `lambda = maxdiag·(1+assist)` exact with sign-exact budget equivalence; the r282 four-language elimination as named negative gates (`CONTEST_ALL_DEAD` — every classical language forces positivity exactly for the positive measure class); typed measurements: one-atom wall, (D) separates, sub-classical `p = 0.38`, ensemble LOW_OUTLIER, the first two MAIN_SEPARATING detectors, the honest margin decay | [E] | `v963` (r282–r285) |
| 9 | **Edge A (fiber, cofinal front):** hypothesis H5 / lemma L2 — pair-sum decay via a non-adjacent/global generic mechanism (address c3, `delta' > 0.21` of available 0.45, TRUTH_ALLOWS); the r287 framework census (consumed in wave 7, `v964`) measures the generic half CLASSICAL: the exact van der Corput bound at `H = ceil(sqrt(m))` delivers `delta' +0.31` world-blindly and certifies 6/7 exceptions — only the kz15 razor stays open (0.07–0.12 dec on every surviving frame; exact-finite per r270) — the vdC form is registered as the named classical candidate: prove it as a *chain* statement (the remaining input: the chain origin of the P-variance scaling; r297 froze the target inequality — sigma ≤ −0.516, truth −0.714, margin 0.198 — and adjudicated the three classical provenance routes `CHAIN_PROVENANCE_PARTIAL(B3)`: the chain is exactly orthogonal on the *window* measure, the named rest is the bordered-form measure transfer + the signed-to-absolute gap; r298 delivered that transfer **exactly and sign-preservingly** — `S_F = B(omega,omega) + B(Delta,omega+beta)` on the frozen positional Fejér kernel, adjudicated `TRANSFER_DOMINANT`: the window main term is empty, `S_F ≈ B(PDelta,PDelta)` — the vdC input IS the Fejér energy of the difference measure; the wave-9 gap object is *prove the Δ-energy decay*; r299 decomposed that energy exactly into its frequency spectrum and adjudicated `DECAY_SPLIT`: the energy is LOWPASS (main-lobe share 0.93), Δ is a *pure c-value difference measure* (the unions share every position, 42/42), and the live edge is the pair split — the diagonal `sum PDelta_j²` already falls fast enough (slope −0.571 ≤ σ* −0.516) with a *falling* pair ratio B/D (−0.168): the frozen rest pair is `DIAG_TARGET` + ratio flatness (plus `CVALUE_TARGET` on the B3 side); the ET/Abel composition fails loudly (+1.948); r300 anatomized the diagonal exactly — `sl_D = 2·sl_L1 − sl_neff` with the effective-carrier count `n_eff = L1²/D` growing at +0.963 against an L1 mass at +0.196: the diagonal decay IS a participation statement; both magnitude factorizations fail (B1 −0.384, B2 −0.346 vs σ* −0.516; the fill decay −0.225 is the non-factorizable part; the |dw| sum-rule attach BREAKS), the **ratio half is settled structurally** (`RATIO_BOUNDED_STRUCTURAL`: the exact kernel-envelope majorant R_env med 1.61 falls at −0.122), and the refined freeze is ONE inequality — `NEFF_TARGET`: prove `slope(n_eff) ≥ +0.908` (measured +0.963, margin 0.055), with the falling D_rank (sp −0.81) as the named bridge to provable terrain; **rounds r296–r300 consumed in wave 9, `v966` = `PRIME.L2.REDUCTION_CHAIN.01` [E]** — the chain stands reduced to exactly this inequality, registered `PRIME.L2.NEFF_TARGET.01` [O]; r301 executed NEFF_TARGET (`neff_target_probe.py`, 32/32, `NEFF_SPLIT(B1|B3)`): the exact identity `n_eff = n_act/(1 + CV²)` with the PERFECT count link `n_act == m` on 42/42 — the carrier count IS the constructive level-2 block count (+1.002) — compresses the rest into **`UNIF_TARGET`**: slope(1 + CV²) ≤ +0.094 (measured +0.039, margin 0.055, jackknife-stable 0/42 below NEED); the weighted-discrepancy/max route is closed honestly (bound slope +0.074, the |dc|-weighted discrepancy is near-flat), the local sum-rule census breaks 24/24; experiments-side, consumed in wave 10, `v967`; r302 executed UNIF_TARGET (`unif_target_probe.py`, 30/30, **`UNIF_DERIVED(B1+A2)`** — the first sealed DERIVED verdict of the lane, every hypothesis typed MEASURED): the normalized block profile of \|PDelta\| is **stationary** (pooled KS(G1,G3) 0.043 over a 3× depth range, an exact-1/N transient onto a finite second moment m2_inf = 1.973 approached from above — the r299 pointwise-cconv negative was the wrong convergence notion), the exact coherence identity `1 + CV² = n_act·χ/(surv²·n_eff_atom)` (χ = D/Q the in-block coherence factor, med 0.63 DESTRUCTIVE and non-growing −0.060, the r301 lag-1 anti-correlation ρ₁ = −0.22 its atom-level mechanism) decomposes the uniformity head exactly; the r301 deep-half flag is carried ENTIRELY by the CV² head (+0.228), not the count (+1.030), and exactly this excess is what the 1/N transient explains; the ONE remaining measured growth statement is **`ATOM_TARGET`**: slope(n_eff_atom) ≥ +0.888 (measured +0.942, margin 0.055 — the same thin margin one level down); the local-pattern route (within-share 0.685) and the recursion-gain damping clause (g/N med 1.079 vs 1.0) closed honestly; experiments-side, consumed in wave 10, `v967`; r303 adjudicated the thrice-identical 0.055 margin (`atom_target_probe.py`, 26/26, **`REGRESS_CONFIRMED`**): the four target margins m_D/m_NEFF/m_UNIF/m_ATOM are ONE algebraic number `S = σ* − sl_D = +0.0547` (invariance ≤ 9e-16; Fractions re-proof; the ½-conversion conjecture refuted) — the cascade r297→r302 is an exact *dictionary* around one measured core, **closed**; the first non-tautological mechanism test (sealed synthetic dc rearrangements, exact marginal, steered ρ₁, 1008 builds) fired **`MIXING_INSUFFICIENT`** with a monotone ladder: flipping the mixing sign kills the inequality (margin −0.044), matching ρ₁ reproduces slopes + margin but not the χ level (0.764 vs 0.630) — the within-block structure beyond lag 1 is the named mechanism gap; `n_eff_atom` is a pure marginal functional (invariance 1.0e-15): stationarity and mixing are complementary; ρ₁-sign census honest 41/42; experiments-side, consumed in wave 10, `v967`; r304 executed the reviewer-mandated short-range-law test (`shortrange_law_probe.py`, 37/37, **`LONGRANGE_STRUCTURE` + `LAW_LONGRANGE`**): the lag profile ρ₁..₁₆ of the dc profile is a *stable period-4 comb* (strong lags at k = 4m and 4m+2 up to 16, halves-stable, world-specific), NOT a decaying tail — the sealed short-range rule finds no k₀ and the reviewer stop case fires; the reviewer condition *splits*: net covariance NC(16) = 0.712 < 1 HOLDS while summability-with-small-tail FAILS (SUM(16) = 1.563); the exact χ lag decomposition (χ = 1 + 2Σ T_k/Q, ward 4.7e-16) shows the χ-relevant within-block structure IS short-range (shares die out by k ~ 4; the r303 miss 0.134 sits at k ≤ 3, dominated by lag 2) — the comb is invisible to χ; lag-8 matching reproduces the χ level (0.652 vs 0.630) but breaks the slopes by 0.028/0.027 — the mechanism is TWO-SCALE, no lag-matching family reproduces both; graduated ρ₁ ladder monotone in χ, NOT in the margin (absolute levels do not kill the inequality — the r303 kill came from per-rung flipped targets); sign pattern −/−/+/+ Fractions-exact on kz18/kz23; the lane's global-profile mixing route is CLOSED (documented stop); experiments-side, consumed in wave 10, `v967`; **rounds r301–r304 consumed in wave 10, `v967` = `PRIME.L2.CASCADE_CLOSURE.01` [E] — the cascade closure and the documented lane stop**: the cascade r297→r302 stands RETYPED as an exact reduction dictionary around one measured core `S = σ* − sl_D = +0.0547` (r303 `REGRESS_CONFIRMED`, invariance ≤ 9e-16; the ½-conversion refuted; hard rule: a round counts only with NEW information); the first causal coordinate is real (ρ₁ ladder monotone, the flip kills at −0.044) but two-scale-insufficient; the global-profile mixing route is documented CLOSED (`LAW_LONGRANGE`: stable world-specific period-4 comb, no k₀; NC(16) = 0.712 < 1 holds, summability fails) — the honest stop state: **L2 generic ⟺ anti-concentration of the explicit dc block field with long-range period-4 structure; return with new tools**; `PRIME.L2.NEFF_TARGET.01` [O] retyped to the documented stop state) | **[O, documented stop]** | `PRIME.PORT.COUPLEDTAU.TERMINAL_CROSSRATIO.01` → `PRIME.PORT.L2.VDC_LEMMA.01` (wave 7, `v964`) → `PRIME.L2.NEFF_TARGET.01` (wave 9, `v966`) → `PRIME.L2.CASCADE_CLOSURE.01` (wave 10, `v967` — documented stop) |
| 10 | **Edge B (base/orientation):** the oriented midpoint theorem — rounds 279–281 executed and consumed (wave 5): the theorem share promoted into stage 7, the generic two-sided obstruction refuted exactly, the center precisified to the budget localization; the successor representation contest is executed and consumed (wave 6, stage 8): `CONTEST_ALL_DEAD` — no classical representation carries the positivity, certified in `v963` with the common deep reason | **[O]** | `PRIME.PORT.FULLSOURCE.PREFIX_RESUMMATION.01` → `PRIME.PORT.RHP.MIDPOINT.ORIENTED_THEOREM.01` → `PRIME.PORT.REPRESENTATION.CONTEST.01` (consumed) |
| 11 | **THE NORTH STAR (the open center, canonical form since wave 6): lemma L\*** — for every admissible anchor `z` and every real polynomial `p ≠ 0` with `deg p < N_w`: `∫p²dnu_z < ∫p²dmu_z`, equivalently `lambda_max(E_{N_w}) < 1` (the nu-dressed mu-CD kernel is a strict contraction at half-filling depth). Via the v963 A-chain this ONE two-measure subordination inequality IS the full free-window quasi-definiteness (`∀ n < N_w : h_n > 0`, the reinstform north star of wave 5), hence via T4 the entire open statement of the wall. State: 42/42 measured, and since r286 **57/57 sign-safe measured** — the sealed extension rule (window cap lifted) added 15 new anchors `N_w` 942–1218, **all margins positive** (min `+1.81e-8` at z = 529, every sub-1e-5 margin mp-verified dps 30/45): **no counterexample**; the honest margin decay is now a measured law (POWER, `alpha ≈ 3.05` vs S, flattening at large S) with the harmless reading quantified (the margin falls *slower* than the local eigenvalue-loading speed at a fixed O(1) census offset; the r258-flat terminal `q_N` re-measured as a different object; r286–r289 consumed in wave 7, `v964` = `PRIME.LSTAR.COHERENCE_CENSUS.01` [E]). Anatomy: (D) per-atom (discrete, sub-classical `p = 0.38`) + (C) coherence (MAIN-destructive, ensemble LOW_OUTLIER, two MAIN_SEPARATING detectors; r288 carrier map: the wall destructivity = antiphase next-nearest (3–4 fold) ARCH–ARCH interference, sign-reversed on every control; plain phase dispersion world-blind — no source-pure separator yet, and the destructive sign flips already at jitter dose 0.005 while the phase field barely rotates; r289 twin adjudication `METRIC_ONLY`: replacing every `log n` by a small-denominator rational multiple of the grid spacing — position cost 2e-9, all exact log-relations destroyed — changes *nothing* (mp-confirmed at the 1.7e-4 margin), while metrically matched fraction shuffles lose exactly as the plain jitter ladder predicts; the coherence reads the tent-split fraction profile at ~1e-3 of the local gap, indifferent to number-theoretic exactness — Baker-class nonresonance input unnecessary on all tested scales; **r331 audit repair (r328C finding C2): the twin RESOLUTION is measured** — at the r289 precision the twin footprint sits ≥ 3.4 orders below every margin down to 1.8e-8 (11/11 anchors RESOLVE, mp pair-verified), the tight anchors are twin-stable from 1e-6 gap, and a gap-preserving genuinely diophantine detune keeps the wall at w9 AND kz119: `METRIC_ONLY_CONFIRMED_AT_RESOLUTION` — the switch carries a resolution certificate on the tight range; **r290–r295 consumed in wave 8, `v965` = `PRIME.LSTAR.CURVATURE_ARC.01` [E]**: the working-set geometry is now *measured* — a soft-shouldered anisotropic tube (killfrac 0.38/0.62/1.00 at 5e-4/1e-3/2e-3 gap), a real one-degree ridge (plateau minC 185, first-order budget threshold in (1.280, 1.291], MAIN-specific), a rank-1 DENS curvature valley (92.5 pct, lam_top −0.418) and simple h_184 flip zeros (alpha = 1; MSTAR_NO_LAW) — while the closed-functional search is sealed honestly NEGATIVE over five class families (r294 F10_FRAGILE + r295 F10_SP_MAJORITY 14/20: F10 unpromoted, R293_LUCK); the successor question is `PRIME.LSTAR.CLOSED_FUNCTIONAL.01` [O]; **r296 consumed in wave 9, `v966`: the DENS-identity fork on the rank-2 DENS core is closed — `DENS_WORLD_BLIND`, coupling number +0.394 < 0.40, the lane routed to L2**). Lean: `lstar_subordination` + `free_window_positivity` (`rh/lean/RH/Window.lean`; L\* ⇒ free-window PROVED); standalone problem statement `rh/problem/lstar_problem.tex`. Predecessor contract `PRIME.PORT.RHP.FULLSOURCE.QUASIDEFINITENESS.01` consumed by `v963` | **[O]** | `PRIME.LSTAR.SUBORDINATION.01` |
| 12 | **The metric firewall (global form):** wall lemma stable under small-dose permutations while consuming the metric data exactly (certified state: perturbative-only) | **[O]** | `PRIME.PORT.WALL.METRIC_FIREWALL.01` |

**The wall.** Every single-component attack class is measured dead
(compression, constants, band locality, mass majorants, level crossing,
one-swap, smooth dressing, bordered spectral bound, naive and IIKS
`s`-coordinates — `DEFINITENESS_WALL_EQUIVALENT`, r257 + r265). The obstacle
is coherent cancellation over the configuration ensemble; the wall contract is
`PRIME.CASE.PAIRCORR.CONTRACT.01` **[O]**. Since r267 the wall has an
**external anchor**: arXiv:2608.13637 (Alpöge–Furman, Lean-4-formalized,
unconditional ≥ 2/3 of zeros simple and on the line) proves a *ceiling*
(`p0 <= 0.6818` for all two-moment bandwidth-one certificates) whose binding
constraint — pair-correlation cancellation beyond diagonal control — is
**exactly the location of our measured wall** (`CEILING_IS_OUR_WALL`), while
being orthogonal in scope (`NO_IMPORT`, frame mismatch quantified: their
zero-atom near-tight frame vs our prime-comb half-filling regime, R ≈ 4/3 vs
73–126).

**The exception scalar.** After r263 the fiber edge is a two-branch statement:
the generic branch is a target-blind triangle certificate (no cancellation
needed); the whole difficulty is concentrated in one scalar `Z_w` on the
7-rung exception branch, with requirement **S6**: `|Z^RHP| + |E| < sqrt(5/7)`
on the branch `g_w < 0`. The r264 charter froze the campaign (SHA
`7b9e751d`); r264 delivered the exact RHP dictionary
`Z^RHP = (border-column readout)/(Y1 root)` and the data signature
(`S1_TERMDRIVE`, separation 0.959); r265 closed the `s`-coordinate class
(no admissible `s`-coordinate carries both the right endpoint and definite
dynamics); r266 searches for a border-dressed resolvent identity; r267
adjudicated the external position. r268 built the local drive asymptotics:
the exception branch is a **cancellation-depth anomaly** at fixed comb mass
(DEPTH power 0.812; mass/edge/participation all blind), the L-gauge frame
provably adds no new locality at the terminal (node collapse exact, dual
degree `N`), and the edge-region triangle bound (outer 20 % of the hull)
delivered the **first source-pure certification of an exception rung**
(kz22, margin +0.04), leaving 6 rungs open by 0.18–0.44 dec
(`DRIVE_LOCAL_PARTIAL`). r269 attacked exactly that gap with sealed
**phase-aware bulk bounds** (block alternation over adjacent sign-run pairs,
exact triangle theorems, chain-envelope variant honestly typed too coarse):
`PHASE_BOUND_PARTIAL(c3F0.30, cert 6/7)` — kz20/22/36/38/39/52 are now
certified source-purely (37/42 on the full ladder), no detector fired; the
razor kz15 misses by 0.12 dec, and its measured perfect-pairing potential is
0.49 dec against the needed 0.44: the remaining kz15 gap sits at the
**measured limit of the local class** at this window split (sub-paircorr).
r270 fought the kz15 boss under the binding reviewer adjudication (no F0.20
level-1 re-optimisation; three sealed attacks): the level-2 anatomy
**confirmed the reviewer's second-cancellation-level hypothesis in
substance** (pair-sum potential pot2 = 0.22 dec at kz15 vs the needed 0.18 —
headroom exists but only +0.04 dec) while measuring its limit: the pair-sum
sign sequence is coin-like (alternation fraction med 0.39), so the sealed
blind level-2 block triangle certifies 5/7 and leaves kz15 at 0.06 dec
(kz39 at 0.002 dec); the a-priori 2/3-absolute-mass split failed loudly
(`SPLIT_RULE_FAILED`, L2 4/7); and the third attack **closed the surface
honestly**: a certified outward-rounded interval evaluation of the full
bordered chain (mpmath.iv, dps 640, every f64 source atom exactly
converted, ~548 digits consumed by interval dependency over 202 steps)
encloses `Z_kz15` with width 1.5e-92 and proves `|Z| < sqrt(5/7)` **strictly,
margin +0.0268** — typed `EXACT_FINITE_CERTIFICATE`, a finite enclosure
statement (EnclOK kind), *not* a mechanism and *not* target-blind class
evidence. Verdict `KZ15_EXACT_ONLY`: the 42-rung surface is closed as
**41/42 mechanism-certified target-blind + kz15 exact-finite**; the open
mechanism edge is the kz15 level-2 sign structure beyond blind pairing
(0.06 dec, sub-paircorr).
r271 promoted the *fixed-form* r269 winner c2PAIR to a **theorem in pure
form** (hypothesis list H1–H5; H5 — the margin `|Z_local| + eps < sqrt(5/7)`
— is the only window-dependent input) with three sealed parameter-free
refinements: the exact boundary group (b1), the level-2 double-pair triangle
(b2, radius 4), the sharp envelope constant (b3, honestly coarser).
Verdict `UNIVERSAL_STILL_PARTIAL(c2PAIR, 5/7)` under the sealed tie-break;
the honest edge: **b2LEVEL2 leaves kz39 at 0.002 dec and kz15 at 0.06 dec**
(full ladder 31/42), no detector fired, the regression ward held. The
N-scaling census typed the cofinal entry door (lemma list L1–L5); the
measured pair-sum trend `sp(N, eps) = +0.67` says lemma L2 (pair-sum decay)
must supply a mechanism the plain pair triangle does not have. The finite
pair algebra is **proved in Lean** (`rh/lean/RH/PairBound.lean`); H5 stays
the documented `sorry`.
r272 dissected that alarm in a pure anatomy round (no new certificate, no
bound modification): the sealed c2PAIR bound factorizes **exactly** as
`eps = A x B x C + tail` (pair count x amplitude x bound-side cancellation,
all ratios to the demand `sqrt(5/7)`), and the measured exponents are
**alpha(A) = +1.01, beta(B/M) = -0.81, gamma(C) = -0.01**: the +0.67 trend
is carried by count x amplitude (`~N^0.2`) while the pairing's cancellation
capture is exactly **flat** (~53 % of the bulk mass at every depth). The
truth side moves the **other way**: the true bulk rest falls
(`sp(N, |R|/M) = -0.31`), the true cancellation deepens
(`gamma_true = +0.45`), and the truth margin `M - |Z|` **rises** with N
(sp +0.36, halves +0.26 → +0.41) — sealed adjudication
**`BOUND_COARSENESS_CONFIRMED`** (possibility B: H5 is not contradicted by
any measured truth trend; the bound loses cancellation, 0.35 → 0.63 dec).
The loss decomposes exactly as `loss = gap12 + slack2`; the growing part is
**c3 slack2** (beyond blind level-2 pairing, 0.22 → 0.50 dec, coin-like
P signs at every lag) — the quantified flip condition
(`FLIP_CONDITION`, typed task, no claim): the L2 lemma needs a
**non-adjacent/global mechanism** supplying `N^{-delta'}` with
`delta' > 0.21`, of the measured available truth decay `0.45`
(`TRUTH_ALLOWS`, ~0.25 exponent margin).
r273 ran the reviewer-adjudicated **mechanism round**: a sealed Euler
perturbation ladder (four graduated source surgeries on the comb — weight
assignment at reach `n^theta`, support jitter at local atom-spacing scale,
Euler-family (SAMEP) decoupling via the r254 integer-root labels, segment
shuffle at preserved local length `n^(1-theta)`; 3 thetas x 3 pinned
replicates x the full 42-rung ladder = 1512 perturbed worlds, conservation
identities exact on all) against `gamma_true`. Sealed outcome
**`PERTURBATION_INSENSITIVE`**: *no* surgery collapses the truth-side
cancellation at any strength — every stage median stays ≥ +0.20 and mostly
at or above the MAIN +0.45 (free permutations cancel even *deeper*): the
`N^{-0.45}` deepening is **near the generic root-scale baseline, not an
arithmetic fingerprint**. The sharp finding is the **firewall map**: every
tested surgery destroys free-prefix positivity (flip degrees 33–43 for
assignment/jitter/window, 66–101 for the family surgery — the r243 flip
census reproduced at ladder scale, pp = 0.00) — **the wall, not the
cancellation rate, is what is arithmetically special about MAIN**. The
resulting L2 hypothesis (typed task, no claim): the missing
`delta' > 0.21` is *generic* level-2 sign equidistribution; the lemma task
is to capture generic cancellation source-purely without crossing the
PAIRCORR firewall — mechanism candidates that need SAMEP coherence, exact
support or assignment structure are measured to aim at the wrong carrier.
In the same round the **Lean pilot was repaired** under the reviewer
adjudication: the three open statements were universally quantified over
pure bookkeeping structures and hence **refutable by trivial models**
(`q_N = 4` ladders, `tau_1 = -1` chains, `2 < 1` margins) — the
refutations are now permanent machine-checked guards
(`RH/Counterexamples.lean`), the window is a concrete structure
(`VonMangoldtWindow`: nodes, signed comb/archimedean weights, window
rule, derived half-filling cap, border vector, budget), and the two
former edges are consolidated into **one master target**:
`augmented_prefix_positive` — positivity of the augmented matrix
`A_{w,n} = [[H_n, u_n], [u_n^T, B]]` through the half-filling cap on a
MAIN window (`MainWindow` = the honest opaque source predicate). The
derivation chain is *proved* in Lean (upper block ⇒ Sylvester prefix
chain `h > 0` = former edge B; Schur complement ⇒ `D_n > 0` ⇒ `q_N < 1`
= former edge A; `D_n = B − u^T H_n^{-1} u` proved via mathlib
`det_fromBlocks₁₁`); the 5/7 floor falls out as a consequence instead of
entering as an axiom. Sorry count 7 → 6, every one truth-capable.
r274 built the reviewer-demanded **smallest hard test** before the oriented
midpoint theorem: the exact **Wronskian dictionary** for the two-sided
midpoint geometry (`wronskian_dictionary_probe.py`, 32/32). The right
solution exists in **three exact constructions** (residue route; *from the
right* with the exact Dirichlet boundary `q_S == 0` — the node polynomial
is the right boundary condition; r230 mirrored dual chain transported
through the r231 L-gauge), and the dictionary closes in both variants:
**base** `h_n = pihat_{n+1} q_n − pihat_n q_{n+1}` with constant `c' = 1`
(exact in rationals, symbolic in z, and at every degree on the mains and
the 42-rung ladder), with the **h-free midpoint form**
`pihat_{n+1} pihat#_{S−1−n} − beta_{n+1} pihat_n pihat#_{S−2−n} = L(z)` —
the r231 sign-blindness as an identity (the gauge-normalized midpoint
Wronskian is the node polynomial, carrying *no* orientation); and
**augmented** `D_{n+1} = B − W^aug_n / W^base_n` with the border-paired
Casoratian `W^aug_n = ∫ (pihat_{n+1} qs_n − pihat_n qs_{n+1}) dsigma =
h_n S_n` (exact against the independent determinant route on the r263
instances, at mp to the w9 terminal `D_N = +0.561250` dev 4.2e−12, and at
the kz15 razor terminal dev 3.5e−11). The honest c-typology: `c_n =
1/W^base_n = 1/(h_0 prod beta_k)` — a source-pure chain-coefficient
product whose positivity on the free prefix *is* the prefix positivity;
the orientation content of the dictionary sits **entirely in c_n** (r231
G12 at dictionary grade), and c is neither the wall nor the target
(fingerprints 0.74/0.36 < 0.9, scopes clean). Verdict
**`WRONSKIAN_DICTIONARY_GO + ORIENTATION_PREVIEW`**: the Prüfer band
`0 < Theta^L − Theta^R < pi` equals `h_n > 0` exactly, MAIN stays in-band
through the whole free window, the controls exit **degree-exactly at
25/21/27**, and the w9 winding census (mp, full continuation) measures
262/366 in-band with first exit at 184 = N_w — the **Maslov round is
freed**, the target `W_{w,n} > 0` (orientation, not magnitude) stays OPEN.
r275 ran the reviewer-authorized **track 2 of the double attack** — the
`l^1 → l^2` currency reform in state-space form (`kyp_memory_probe.py`,
24/24): the exact level-2 block decomposition of the bulk rest (run pairs,
offset 0, odd tail its own block) as a dimension-2 state system
`x_j = (acc, 1)`, `T_j = [[1, b_j], [0, 1]]`, `C_j = (0, b_j)` (pure exact
bookkeeping, `acc_J == Z` rational dev 0 on 47/47 worlds), with the
binding KYP memory ladder constant → even/odd → period 4 → Riccati
adjudicated **exactly** (rational LMI checks, no floats in any gate). The
outcome is a **two-sided no-go**: (i) every uniform quadratic memory is
**INFEASIBLE by structural algebra** — the homogenization fiber forces
`p11 = 0` (killing terminal usefulness) and the descent windows
`sum b^2 > 0` kill the LMI chain outright, machine-verified on 46/46
worlds including the trivial SMOOTH anchor (a structure theorem about
additive accumulation of sign-free drives, not a cancellation issue);
(ii) the unique admissible source-dependent memory (backward Riccati) is
**FORCED to carry the exact signed tail sums** — its budget is
`Z^2 + sum b^2`, the inverted terminal readout itself
(`KYP_IS_TARGET_INVERSE`, fingerprint sp +1.00; it "certifies" the
broken SCRAMBLE world — the friendliness signature — while failing EPST
exactly where its wall kips at the terminal); and (iii) the naive generic
root-scale `l^2` repackaging is QM-AM-coarser than the sealed `l^1`
pairing on every world (overhead med +0.15 dec, generic census 8/42).
Reading: the covariance a quadratic form would need is exactly what r273
typed GENERIC — per-window data, not source law; quadratic storage
without a sign/covariance law either dies (uniform) or reads the answer
(Riccati). The L2 task keeps its r272 address (non-adjacent/global
mechanism, `delta' > 0.21`) — with the state-space vehicle now excluded
alongside `l^1` refinement.
r276 measured the reviewer-adjudicated **dose-response curve of the
arithmetic** (`minimal_firewall_probe.py`, 25/25): the r273 firewall map
(every surgery kills the wall at binary strength) refined from binary to
graduated to **minimal** — five sealed surgeries (neighbor weight swap,
support jitter at `theta` × the local atom gap, family decoupling, and two
new grid-density surgeries separating the **sign field** from the
**magnitude field**) at doses SINGLE + `theta` 0.02–0.25 on three windows
(360 worlds, conservation and dose accounting exact, every borderline flip
mp-counter-checked at dps 40 — including the v956 boundary flip **at**
`N_w = 184`). Sealed outcome **`FIREWALL_LAW` + `CONTINUUM`**: the wall is
**not** an all-or-nothing firewall — single minimal operations typically
leave survival depth 0.88–1.00 of the free window (position-dependent
lethal exceptions at degrees 51–152), the v956 control level (0.11–0.15)
is reached only inside the graduated ladder, and the deficit laws are
shallow and saturating (`D ~ theta^b`, `b` +0.04…+1.09 — immediate onset,
no critical dose). The property ranking is the round's surprise:
**support exactness is the most wall-critical property** (jitter at 2 % of
the local atom gap already costs 3/4 of the depth; tolerance ranking
P2_JIT 0.343 < B2_MAG 0.389 < B1_SIGN 0.536 < P1_SWAP 0.621 < P3_FAM
0.700) — the wall reads the **metric placement** of the signed source
(positions, magnitude-position pairing) more exactly than the Euler-family
bookkeeping; the map transports across windows (sp +0.84/+0.86). Typed
MEASUREMENT_ONLY: the firewall hypothesis (wall carried by exact metric
placement, graded loss law) is a quantified falsifiable task, not a
mechanism claim.
r277 ran the reviewer-demanded **blind Prüfer/Maslov census**
(`maslov_census_probe.py`, 29/29) — predict the flip degrees from
phase/counting data instead of restating `h_n > 0`. The sealed protocol
(training = w9 + the 4 smallest-N rungs + their scramble(seed 2)
variants; three sealed candidate classes with fixed priority; blind =
the remaining 37 rungs + w12/w13/w26 + EPSTEIN/SCRAMBLE/SMOOTH) selected
**R2 = interlacing/reality of the Jacobi zeros** and passed the blind
GO: **37/37 rungs + all mains SAFE at full depth, controls fire at
26/22/28 = flips 25/21/27 + 1** (inside the sealed ±1), r259
energy-branch separation 17/6/9, and **not h-equivalent** (79 pattern
mismatches against the 78 h re-entry pivots on the w9 continuation —
the rule is a one-way break detector, not the sign chain). The round's
central finding answers the r274 winding warning: the **raw atom-counted
Sturm census is NOT the winding quantity** — on MAIN it breaks at
n = 56/48 with h positive throughout (zeros escape the atom hull, then
pair up inside single atom gaps: the positive-measure separation theorem
genuinely fails for the signed comb), while the interlacing/reality
structure provably holds on every positive prefix (β > 0 ⇒
symmetrizable ⇒ real + Cauchy-strict) and breaks at exactly `nf + 1`
everywhere (w9 continuation fire 185 = N_w + 1).
`STURM_CHAIN_VERIFIED` is honestly **not** awarded (the c1 atom-count
expectation is refuted — that refutation is the finding); the Sturm
rotation `#eig(J_n) < x* == n − #signchanges` is machine-exact on both
mains, and the oriented midpoint theorem (Maslov index through the
interlacing break) is the named next round.
r278 turned the r276 continuum into the reviewer's demanded **local
stability law** (`metric_stability_probe.py`, 31/31): the **exact
first-order sensitivities** of the wall against atom positions, derived
Hellmann–Feynman-style (`d log h_n/du_j = ⟨wdot_j, pihat_n²⟩/h_n` — the
polynomial variation drops by orthogonality) and via the border/CD kernel
(`dF_n/du_j = −⟨wdot_j, pihat_n B_n⟩`, `B_n = Σ_{k<n} F_k pihat_k/h_k`),
pushed through the sealed tent-grid channel (piecewise linear, lever
`m_j/(2D)`) and machine-gated against Richardson FD and mp (dps 40) to
1e−7. Structural disclosure: on w9 `alpha = log 16` makes the grid
commensurate with `log 2` — the **entire 2-power family sits exactly on
tent nodes**, so the gradient is carried as a one-sided pair with
side-selected predictions. Sealed outcome **`METRIC_STABILITY_LAW` +
`GRADIENT_EXPLAINS_DOSE(PARTIAL)`**: the local law exists (curvature
TAME, consistency 1.02), the stability inequality
`margin(c) ≥ margin_MAIN − L_D(w)·theta − O(theta²)` holds with exact
measured `L_D(w)`, **but its strict window `theta*` ~ 8e−4 (w9) / 1.2e−4
(kz55) sits 25–170× below the smallest r276 dose**: the measured 0.02
collapse is already a nonlinear cascade (first order under-predicts the
losses ~2.4×, ratio 0.41, and reaches the flip criterion on only 2/9
replicates). The **u-profile** (the r276 follow-up) is closed and
surprising twice: the gap-weighted sensitivity is **bottom-loaded — the
small primes 2, 3, 5 carry it** (hull correlation −0.80…−0.83, weight
+0.74…+0.82), it **predicts the single-op lethality** (deterministic
±gap singles, 66/70 flip, sp −0.82 — the r276 lethal exceptions located),
and the terminal-margin map shares the same geometry (sp +0.83). The
N-scaling says the margin sensitivity does **not** grow with N (halves
slope −1.09, non-monotone, decade 17.3): no fragility razor, but also
**no uniform stability constant** — the honest `MainWindow` consequence
is typed **PERTURBATIVE_ONLY**: `MetricNear(c, MAIN, theta) :=
max_j |u_j − u_j^MAIN|/g_j ≤ theta ≤ theta*(w)` inherits wall positivity
from MAIN positivity per window with window-specific `L_D(w)`; the MAIN
positivity itself remains the open center (H5 untouched). Control
contrast: at low degrees the arithmetic combs (MAIN, EPSTEIN) are
position-**blind** (~1e−14) while the metrically wrong worlds (SCRAMBLE,
SMOOTH) already carry O(1e−2…1) sensitivity — the low-degree gradient
separates metric randomness from arithmetic structure; the wall depth
separates the arithmetic combs among themselves.
**Wave 4 (2026-08-25)** promoted rounds r260–r278 into the load-bearing
suite as **`v960_terminal_surface_closure.py`** (r260/262/263/268–271
theorem set + r272/273/275 typed measurements; ledger
`PRIME.PORT.COUPLEDTAU.SURFACE_CLOSURE.01` [E]) and
**`v961_midpoint_orientation_dictionary.py`** (r274/277 theorem set +
r276/278 typed measurements; ledger
`PRIME.PORT.RHP.MIDPOINT.ORIENTATION_DICTIONARY.01` [E]) — both embed
their probes byte-exact and execute them in the sealed `--smoke` stage
(the full records stay sealed experiments-side and are re-verified by
this suite); two new open contracts registered:
`PRIME.PORT.RHP.MIDPOINT.ORIENTED_THEOREM.01` [O] (round 279 in flight)
and `PRIME.PORT.WALL.METRIC_FIREWALL.01` [O]. The rh/ pipeline gate
landed: `bash build.sh audit` now runs `run_rh.py --fast` as its own
"RH workspace" section.
r279 ran the **oriented-midpoint proof round** for exactly that first
open contract (`oriented_theorem_probe.py`, 32/32,
`PRIME.PORT.RHP.MIDPOINT.ORIENTED_THEOREM.01` — experiments-side, no
status move) — build the theorem or name the missing step-5 ingredient
exactly. The round constructs the **two-sided index bilanz** (left
chain degree `n` vs dual degree `S−1−n` against the S common nodes of
the gauge polynomial `L`) and delivers two new exact, world-blind
theorems: the **TWO_SIDED_PARITY_THEOREM** — from the r274 node
identity and the `L'` alternation, the union polynomial
`Q_n = pihat_n · pihat#_{S−1−n}` has `sign Q_n(x_j) =
e_j (−1)^{S−1−j} sign(h_n)`, so its real zeros are ODD-many in every
weight-agreement gap (≥ 1 forced, never empty) and EVEN-many in every
disagreement gap, at **every** degree and independent of the h sign
(exact in rationals on the toys; 246 212 sign gates + full-degree
occupancy sweeps on the real combs; the |D| weight-sign boundaries are
exactly the slack of the bilanz, and the census identity
`c + c# = (S−1) − |D| + 2·scD` is exhaustion-proved to k = 8); and the
**CROSSING_BUDGET_THEOREM** (Jacobi's minor-sign rule + Sylvester
congruence): the number of crossings over the FULL algebraic
continuation is **exactly the negative-atom count**,
`#{n < S : h_n < 0} = S_−`, machine-gated on every world (w9 104, w13
98, EPST/SCR/SMOOTH 141/94/6, kz15 121, kz52 551). The sealed
candidate adjudication is honest: **STEP5_OPEN** — no non-restatement
invariant kips at `N_w` (tight two-sidedness never fails anywhere = a
theorem, not a separator; union reality / shared gaps / left-in-D die
at degrees 1–4 world-blindly; the first crossing itself is the h
restatement), and the b2 obstruction implication is
**OBSTRUCTION_REFUTED by exact rational counterexamples** (the
two-sided compatibility state does not forbid a crossing). Two
measured surprises (disclosed amendments a1/a2): the localization
equality is **not universal** — w13 crosses first at `N_w + 2`, kz15
at `N_w + 1` (w9 and kz52 exactly at `N_w`) — and the r277 R2 anchor
is co-located with the **first crossing** (`min C + 1`), not with
`N_w`. Verdict **`MAIN_SPECIFICITY(LOCALIZATION_IS_THE_ARITHMETIC)`**:
every theorem of the bilanz holds on all controls before *and* after
their flips (25/21/27) — the two-sided machinery carries no
arithmetic; the **entire** arithmetic content of the wall is *where*
the world-blindly fixed budget `S_−` is spent. The precise
Lean-suitable gap statement (the honest step-5 sorry, matching the
`Window.lean` master theorem): for every MAIN window,
`∀ n < N_w : h_n > 0`, i.e. `min C ≥ N_w` — the localization of the
Jacobi-fixed crossing budget into the upper half of the continuation.
NO RH CLAIM.
r280 measured the **localization anatomy** of exactly that budget
(`budget_localization_probe.py`, 29/29,
`PRIME.PORT.WALL.BUDGET_LOCALIZATION.01` — experiments-side, no status
move): the **full 42-rung offset census** — `min C − N_w` takes the
values {+0: 18, +1: 10, +2: 6, +3: 6, +4: 1, +5: 1} across the entire
frame-A ladder (max **+5** at kz43; the v956/r279 anchors 0/2/2/3/1,
kz15 +1, kz52 +0 reproduced exactly; half-filling `N_w = ceil(S/2)` on
42/42; mp dps-40 arbitration exact), with `sp(offset, N) = +0.096` over
`N ∈ [142, 878]` — the offset stays O(1) with **no N-growth trend**, and
the r276 dose curve is the same coordinate (`min C == nf` definitional;
the P2_JIT medians 0.250/0.207 reproduced with the exact r276 seeds).
The round's core finding is the **extremality answer**: the w9 extremal
localization `min C = N_w` is **NOT a local maximum** of the
localization functional — the r278-gradient directions
(OPT/OPT_SAFE/SMALLPRIME) push the first crossing **past half-filling**
(`min C` 184 → 185 at `theta` 7.7e−5…1.6e−4, all mp-confirmed;
first-order census: every second random direction raises, because the
crossing is *shallow*, |h_184| rel margin 1.9e−3) — the b3 variational
hypothesis is **REFUTED in its local-maximum form**; the crossing
log-gradient is anti-aligned with the last free pivot on all three
tested worlds (`cos_w ≈ −0.97`, the h-sign flip of a raw-gradient
lockstep — `CRITICALITY_STRUCTURED` by the sealed cosine rule); the
offset windows kz15/w13 stay UNRESOLVED at the tested doses. The
**duality lens is exact**: `h_{S−m}·h#_{m−1} = 1` and the sign mirror
`sign h#_k = sign h_{S−1−k}` (toys in rationals, w9 at mp) translate the
localization statement exactly — `min C ≥ N_w` **is** the confinement of
all `S_− = 104` dual negative pivots to the lower dual half, and w9
**saturates** the bound (max negative dual pivot 182 = `S−1−N_w`);
typed `DUAL_RESTATEMENT` (same signature, no isolated carrier). The
moment view is honest: the natural mu/nu **Weyl minor-perturbation
argument dies at n = 3** (`X_weyl = 2` of `N_w = 184`, 1.1 %), the
raw-moment rest zone does not exist on any world, and the PAIRCORR
detector types **both** arguments WORLD_BLIND (SMOOTH exceeds MAIN at
Weyl level) — neither carries the localization; the r279-b3 gap
statement stands unchanged. NO RH CLAIM.
r281 ran the **pinning-anatomy round** on the reviewer's reframed
question — not maximality but *pinning*: why is half-filling the
natural pinning point? (`halffilling_pinning_probe.py`, 28/28,
`PRIME.PORT.WALL.HALFFILLING_PINNING.01` — experiments-side, no status
move). The round's exact core is the **counting theorem**: pivot `h_n`
consumes `m_0..m_{2n}`, an S-atom signed measure has exactly S free
moments (Vandermonde bijection, freedom demonstrated exactly: `dm =
e_{S−1}` moves the last free pivot alone), all `m_k, k ≥ S` are forced
by the node-polynomial L-recurrence (toys exact in rationals; w9 at mp
dps 200 with max rel residual 7.6e−157) — so the **free pivots are
EXACTLY `h_0..h_{N_w−1}`** and `h_{N_w}` is the first forced pivot:
"why half-filling" has the exact answer **because the free moments end
there**. The upper-side theorem attempt is adjudicated honestly:
**`UPPER_PINNING_MEASURED`** — the v956 boundary statement ("the wall
dies immediately past the cap") is two computational paths + an mp
ward on five windows, i.e. *measurement*; the measured pinning
constant is `C = 5` (42-rung census, minC ≥ N_w everywhere); the only
*proven* upper bound today is the pigeonhole ceiling `minC ≤ S − S_−`
(from the r279 budget theorem; w9: 184 ≤ 263 = p, 79 degrees above
`N_w`); and **world-blind O(1) upper pinning is REFUTED exactly** — a
single-tiny-negative toy measure has `minC = S−1` (offset `N_w − 2`,
unbounded in S): any O(1) upper side must consume the comb structure.
The offset anatomy is a clean negative: six sealed source-pure
tail-coupling candidates (last free margin, margin slope, first-forced
moment cancellation, negative-weight fraction, razor position) reach
max |spearman| = **0.273** vs offset over the 42 rungs —
`OFFSET_UNSTRUCTURED`, no predictive offset formula in this census.
The relay question closes honestly too: the per-degree condition
`B_n = [beta_n > 0]` is typed **RESTATEMENT** (no h-blind margin proxy
predicts the flip location on all five worlds), and the r280
anti-alignment `cos_w ≈ −0.97` is measured to be the **h-sign flip of
a raw-gradient lockstep** (`cos_raw ≈ +0.96..0.98` at every crossing;
all four dead-world crossings type CREEPING — crossing *type* does not
separate the worlds). Verdict **`PINNING_DECOMPOSED`**: pinning =
counting theorem (upper location, exact) + forced-pivot death within
O(1) (measured) + free-window survival (the lower side — the open
center, now with the sharpened reading "why does MAIN survive all its
free degrees"). NO RH CLAIM.
r282 ran the **representation contest** on the north-star question
"what IS `h_n`?" (`representation_contest_probe.py`, 30/30,
`PRIME.PORT.REPRESENTATION.CONTEST.01` — experiments-side, no status
move): four classical representation classes, each with a construction
attempt on exact toys and one brutal gate — the representation must
force MAIN positivity *structurally* AND provably reject an
early-dying control (calibrated on a positive-class toy, where the
classical structures MUST exist — and do). All four die with exact
eliminations: **K1 (RHP/monodromy)** is a RESTATEMENT — the Pontryagin
split `h_n = P_n − N_n` is exact everywhere, but the
`|m|²`/sum-of-squares structure exists iff the negative register is
empty (only the positive class; MAINLIKE has N/P = 0.538 at its last
free degree, w9 carries negfrac 0.498 at the wall), and the surviving
clause `N < P` is `h > 0` verbatim. **K2 (Hodge/Kasteleyn
orientation)** is NOT CONSTRUCTIBLE — full exhaustion over all `2^S`
gauges shows the coherent orientation class over the Cauchy–Binet
configuration space is *exactly* `{±sign w}` (the quadratic/edge class
collapses at the n = 3 cells), and that unique gauge changes the value
by exactly 2×(negative mass) (w9 defect 0.115): orientation exists iff
`S_− = 0`, and MAIN is signed — the r261 no-go upgraded from a pairing
failure to a **class obstruction**. **K3 (de Branges/canonical
system)** is a RESTATEMENT — Hamiltonian positivity ⟺ `h > 0` exactly,
and the only independent shadow (R2 interlacing) lags +1 on every
flipping world (toys exact, controls 26/22/28/26 = flip+1, w9 SAFE
through 183): it detects, it never forbids; forcing the Hamiltonian
positive (`|β|`) provably changes the measure. **K4 (operator/dual
pair)** is WORLD-BLIND by theorem — `h_n · h#_{S−1−n} = 1` exactly, so
the dual pair is sign-synchronized *identically* in every world: the
"joint positive structure" adds no second constraint, and global
positive structures are killed by the r279 budget theorem. Overall
verdict **`CONTEST_ALL_DEAD`** — a program finding: the pivot
positivity has *no representation in the four classical languages at
this resolution*; each language forces positivity exactly for the
positive measure class, and MAIN is signed. The lower side (why MAIN
survives its free window) stays the open center. NO RH CLAIM.
r283 ran the **big RHP attack** on the north star (full source →
monodromy property → quasi-definiteness of the complete free window;
`fullsource_quasidefiniteness_probe.py`, 26/26,
`PRIME.PORT.RHP.FULLSOURCE.QUASIDEFINITENESS.01` — experiments-side,
no status move): three sealed candidates. **A2 (monodromy
factorization over μ/ν)** is the carrier — verdict
**`MECHANISM_PARTIAL`**: in the μ-orthonormal frame the full Hankel
form of μ̃ is `I − M_n`, and the node Gram `E_n = B_n B_nᵀ` is
*exactly* the v881/r224 IIKS kernel (μ-CD kernel dressed by the ν
weights). Steps s1–s4 are EXACT (split; unit-triangular congruence
`minor_k(D_μ − G) == D_k(μ̃)` — rationals on five toys, slogdet dev
3.4e-9 on w9; frame contraction `h > 0` through degree m ⟺
`λ_max(E_m) < 1`; monotone rank-one loading ⇒ `λ_max` crosses 1
**exactly at minC+1**). World test exact on all seven worlds: w9 185
= minC+1, w13 171 = N+3, kz15 205 = N+2, controls 26/22/28/26 =
flip+1 — the *same* degrees as the r277 R2 fires: the contraction IS
the split-source form of the r277 reality/interlacing condition. The
capacity step is adjudicated both ways: the pigeonhole ceiling
`minC ≤ S₊` is *re-proved* in the frame (null-polynomial witness,
exact), and capacity-as-COUNTING is *refuted* exactly (rank-one pair
with identical rank/support and opposite window fate; the deciding
value is metric: 1.25e-13 vs 1.25e+2). What remains is ONE formal
scalar, the **missing lemma L\***: for MAIN, every real polynomial
`p` with `deg p < N_w` has `∫p²dν < ∫p²dμ`, i.e. `λ_max(E_{N_w}) < 1`
(measured margin 1.68e-4, top eigenvalues 0.99983/0.99874/0.99597 —
broadly tight, not a single-mode squeeze); given s1–s4 this single
inequality ⟺ full free-window quasi-definiteness. **A3
(free-window Fredholm determinant)** COLLAPSES onto A2 (zero-freeness
of `det(I − sE_{N_w})` on (0,1] is the contraction scalar; s* = 1/ρ
dictionary + degree-resolved slogdet co-location 185/26). **A1
(Schwarz/J-symmetry)** is DEAD at the residue-positivity step: the
reflection symmetry of the real jump data holds identically on every
world (world-blind, cannot carry), the Herglotz condition fails on
MAIN itself (S₋ = 104 negative residues), and the windowed η*
anatomy is world-blind by the sealed distance rule. Forbidden-list
bilanz: triangle/Gershgorin bounds die at n = 21 ≪ 185; no
target-inverse quantity, no truncation, no posthoc window (mutants
flagged). The RHP language holds exactly ONE invariant on this
question — the contraction scalar; the open center is L* itself.
NO RH CLAIM.
**Wave 5 (2026-08-25)** froze rounds r279–r281 as a **small mathematical
theory** (the reviewer demand: "not just as experiments") in
**`v962_halffilling_pinning_theory.py`** (ledger
`PRIME.WALL.HALFFILLING_PINNING_THEORY.01` [E]; the three probes embedded
byte-exact in the sealed `--smoke` stage, plus a module-own exact section
in pure Fractions): four named theorems — **T1 moment counting** (the
free pivots are exactly `h_0..h_{N_w−1}`; "why half-filling" = the free
moments end there), **T2 crossing budget** (`#(h<0) = S_−`,
Jacobi/Sylvester, world-blind), **T3 two-sided parity** (h-blind; census
bilanz + exhaustion), **T4 main window reduction** (the entire open
statement is `minC ≥ N_w` ⟺ `∀ n < N_w : h_n > 0`) — and four named
refutations (**no universal O(1) pinning** — exact, unbounded-offset
one-negative family; **no extremality** — the w9 crossing is liftable
past half-filling; **no generic Maslov obstruction** — 3 exact rational
counterexamples; **no simple offset law** — max |sp| 0.273, all
world-blind). Rounds 279–281 are consumed on
`PRIME.PORT.RHP.MIDPOINT.ORIENTED_THEOREM.01` [O] (theorem share
promoted, center precisified); the two successor lanes are registered as
`PRIME.PORT.REPRESENTATION.CONTEST.01` [O] and
`PRIME.PORT.RHP.FULLSOURCE.QUASIDEFINITENESS.01` [O] (the r282/r283
rounds are sealed experiments-side and NOT consumed by wave 5). In Lean
the hole is now fog-free: `free_window_positivity` is the one central
`sorry` in base coordinates (T1/T4 proved for real, T2 stated, the exact
one-negative instance a permanent counterexample guard; build green,
8 intentional `sorry`s). NO RH CLAIM.
r284 read the missing lemma L\* as a **two-measure subordination
problem** (`lstar_two_measure_probe.py`, 30/30,
`PRIME.PORT.LSTAR.TWO_MEASURE_ANATOMY.01` — experiments-side, no status
move): ν (104 atoms) strictly μ-subordinate (263 atoms) on polynomials
of degree < `N_w`, with the structurally motivated Nyquist/shielding
hypothesis (resolution `π/n` on the θ grid vs the ν→μ pairing scale)
adjudicated sealed. Outcome, honest on every leg: **`SHIELDING_BLIND`**
— the source-pure pairing geometry is a *builder-channel* property (on
w9 every ν atom sits exactly one fold step from a μ atom, 100/104 fully
μ-neighbored; all four dead controls mirror the same near-perfect
pairing; all four sealed geometry statistics world-blind);
**`NYQUIST_ORDERING` with the honesty clause** — the sealed crossing
formula `n_pred = L/(4s)` holds the factor-2 band on all 42 rungs but
**fails on all four controls by 1.5–3 octaves** (they die at 22–28
despite near-perfect pairing): the support-geometry form of the
resolution hypothesis is refuted where it matters (the +0.998 spearman
is carried by trivial `N_w` scaling); **`CHRISTOFFEL_RANGE`** — the
round's core: the exact sandwich `max_k v_k K_n(y_k) ≤ λ_max(E_n) ≤
trace(E_n)` brackets the wall, and on w9 the pure **single-atom
Christoffel bound crosses at `n_DIAG = 187`, only two degrees above the
true crossing 185** (coherent assist 3.1 %, aggregate ν coherence
destructive: slack 50×), while on all four controls `n_DIAG` is never
reached — **MAIN's crossing is a near-single-atom Christoffel event,
control death is a collective mode**; **`EXTREMAL_ANATOMY`** — the
near-breaking direction is 99.7 % two background atoms at the shallow-u
hull edge *below the first prime* (folds 2/4, `u` = 0.03/0.06, the two
weakest-shielded ν atoms; SMALLP_DEPLETED — the r278 small-prime
sensitivity profile and the near-null direction of `E` are different
coordinates), the top-3 eigenvectors are ONE oscillatory edge band, and
the controls' breaking modes are more diffuse (PR 5.0–9.7 vs 1.89).
Detector: all seven sealed stats world-blind — the arithmetic lives in
the spectral coherence of the weights, not in the support geometry.
mp dps-60 ward: `ρ_184 = 0.99983248 < 1 < 1.00003660 = ρ_185`.
NO RH CLAIM.
r285 adjudicated the **decomposition perspective on L\***
(`christoffel_decomposition_probe.py`, 33/33,
`PRIME.PORT.LSTAR.CHRISTOFFEL_DECOMPOSITION.01` — experiments-side, no
status move): L\* split candidate-wise into **(D)** the per-atom
diagonal condition `v_k K_{N_w}(y_k) < 1 − δ` and **(C)** the coherence
condition `assist = λ_max/maxdiag − 1 < δ'`, with **exact bookkeeping**
(multiplicative identity `λ_max = maxdiag·(1+assist)` and the budget
equivalence `margin > 0 ⟺ ρ < 1` gated at every degree; the additive
Weyl form measured coarse, slack 0.287 — the multiplicative form is the
honest one). Verdict **`DECOMPOSITION_EXACT`**: the ladder budgets are
measured on all 42 rungs (`maxdiag_{N_w}` 0.9343–0.9941 with
`sp(N, maxdiag) = +0.65` — **the per-atom budget tightens toward 1 with
N**; `assist_{N_w}` 0.0059–0.0702 with `sp = −0.65` — the coherent
assist shrinks: **no fixed δ/δ' budget is supported by the measured
trend**); the controls rip at (C) at their death degree (maxdiag
0.26–0.38, assist 1.69–2.99, off-diagonal share 0.63–0.75).
**`PERATOM_CLASSICAL`**: the binding atom is the shallow-u ARCH edge
below the first prime on 42/42 rungs (w9 leaders folds 2/4/6/8/10,
near-degenerate 0.9411–0.9700); the classical arcsine yardstick (typed
COMPARISON_ONLY) measures the binding atom's kernel growth at exponent
**0.38 — sub-classical** (bulk law 1): the classical isolation degree
`n*_bind = 58` misses `n_DIAG = 187` by 1.68 octaves and covers 0/42
rungs — **(D) survives *because* the discrete kernel grows
sub-classically at the binding atoms**; any provability route for (D)
must be a discrete bound, not continuum asymptotics. The round's
surprise is the **b3 answer**: (D) is **not** world-blind at window
scale — the controls' `n_DIAG` = 50–91 (≪ 184) vs MAIN's 187 > 185:
the per-atom half itself separates the worlds, refuting the clean
split "(D) metric-blind + (C) carries all arithmetic" *as stated*
(sealed verdict falls to `COHERENCE_UNSTRUCTURED` by the blind-D
clause), while at the **death degree** the split is clean and the
assist separates sharply (MAIN 0.0195 at its crossing vs controls
1.69–2.99, `MAIN_SEPARATING`; z-scores: MAIN −3.15 **destructive** vs
controls +4.95..+12.44 **constructive** — sign-clean, band-MIXED).
**`ENSEMBLE_POSITION`**: against both conservation-gated replicate
ensembles around MAIN (16 sign re-assignments, 12 position scrambles;
crossings 5–28), MAIN's **wall assist is the extreme low outlier
(pct 0.00)** while its early assist is scramble-generic (pct 0.75) —
the near-single-atom wall is MAIN-specific, the r273 genericity
finding measured from the other side. mp dps-60 ward incl.
`maxdiag_184`. NO RH CLAIM.
r286 measured the **decay law of the L\* margin and hunted the
large-S counterexample** (`lstar_margin_scaling_probe.py`, 29/29,
`PRIME.PORT.LSTAR.MARGIN_SCALING.01` — experiments-side, no status
move; the probe runs the **standalone document pipeline**
`rh/problem/verify_lstar_instance.py` verbatim and re-gates it
against the repository builders on the flagship *and* every new
anchor, devs ≤ 2.9e-15). **The law** (`MARGIN_LAW`): on the 42-rung
family the margin `1 − λ_max(E_{N_w})` follows a **power law in S**
(Theil–Sen `alpha = 3.05` vs S / 3.06 vs `N_w` / 2.39 vs z; MAD 0.61
vs 0.85 for the exp form; halves 3.67/2.07 — the law **flattens** at
large S, the honest anti-extrapolation flag), and the driver is
neither the diagonal nor the assist alone but the **cancellation
sharpening**: in the exact identity `margin = (1 − maxdiag)(1 − c_w)`
the factor `c_w = (λ − maxdiag)/(1 − maxdiag)` rises to 0.999991
while `1 − maxdiag` only decays with slope −0.69 — the r284/285
shallow-edge extremal band persists on all 42 rungs (PR median 1.81).
**The hunt** (`NO_CROSSING_EXTENDED`): the sealed extension rule (the
document admissibility rule with the window cap lifted, comb
completeness `z² ≤` table cap) admits 83 candidates; the first 15 by
`(N_w, kz)` — `N_w` 942–1218, S 1883–2435, z 101–547 — were built,
cross-checked and measured under the sealed sign-safety tier protocol
(f64 two-route equivalence via the r283 dictionary `margin > 0 ⟺
minC ≥ N_w`, plus mp dps-30/45 staggered verification on every
sub-1e-5 margin; 32/42 + 15/15 mp-verified, worst f64-vs-mp rel dev
6.6e-6): **all 15 margins strictly positive** (min `+1.806e-8` at
kz=119/z=529), the triple-verification counterexample protocol never
fired, and the 42-rung power fit **predicts** all 15 out-of-sample
margins inside one decade (`EXTRAP_CALIBRATED`, max dev 0.78 dec).
**The census** (`CENSUS_EXTENDED`): `minC − N_w` on the new anchors
is {0: 1, +1: 10, +2: 2, +3: 1, +4: 1}, max +4 — **the O(1) wall
offset survives beyond the family cap** (by the r283 equivalence
this is exactly why the margin *cannot* cross there). **The
reconciliation** (`QN_RECONCILED`): the r258 terminal bordered ratio
`q_N` measured on the same 57 rungs stays FLAT and O(1)-positive
(min `1 − q_N` = 0.019, spearman +0.35) while the spectral margin
falls 4.0 decades (spearman −0.87) — **two different objects**: the
exact Rayleigh inequality `λ_max ≥ ‖B[:, N_w−1]‖²` holds on 42/42
with the terminal direction far from the top (median 0.448), and the
extremal coefficient vector is **full-sector** (terminal-degree share
median 1.9e-8, degree PR median 188) — the terminal pivot cannot see
the sector the spectral margin lives in. The sealed HARMLESS/DANGER
adjudication (two disclosed calibration amendments, r270 precedent;
the original rule's outcome printed permanently) lands **HARMLESS**:
no crossing + census O(1) + the margin falls *slower* than the
pre-wall eigenvalue-loading speed (`margin/v8` 19–652, median 144 —
the benign direction) + `q_N` flat-positive: **the margin decay is
the spectral signature of the known O(1) wall offset (kinematic
eigenvalue-density shrinkage), not an approach to failure** — under
both sealed law forms the margin is positive at every finite S, and
any counterexample claim needs the measured census, not the fit.
NO L\* claim. NO RH CLAIM.
r287 ran the **standard-framework census for the L2 lemma** (the
30-pct lane of the reviewer split 70 L\* / 30 L2;
`l2_deterministic_cancellation_probe.py`, 25/25,
`PRIME.PORT.L2.DETERMINISTIC_CANCELLATION.01` — experiments-side, no
status move): WHICH classical deterministic cancellation framework
captures the level-2 block sums `P_j`? Four sealed frames, each an
exact per-rung bound on `|sum P_j| = |R|`, tested against the r272
flip condition `delta' > 0.21` and the 7-exception certification.
Verdict **`FRAMEWORK_PARTIAL`**: **F1** discrepancy/Erdős–Turán (exact
Abel form `D* × V`) dies loudly — both factors grow (`D*/sqrt(m)` med
1.90, no `N^{−β}`), delta' −0.21, cert 0/7, paircorr demand **FIRES**
(max +1.44 dec); **F2 van der Corput** at the frozen window
`H = ceil(sqrt(m))` — an **exact world-blind theorem**, no constant,
no fit, no arithmetic — delivers **delta' +0.31 > 0.21** and certifies
**6/7 exceptions** (kz20/22/36/38/39/52; ladder 38/42, demand max
+0.04); **F3** Salem–Zygmund sub-Gauss (measured constant `C` med 0.90
max 1.68, spread 1.86 STABLE, controls covered) delta' +0.38, 6/7;
**F4** MDS/Azuma against the sign-bin chain filtration (md_ratio 0.29
— the `P_j` are MDS-like; `C4` max 0.68 stable, covered) delta' +0.41,
6/7. The sealed `(cert, delta')` key adjudicates F4; the honest
research edge is **F2, the exact theorem at 6/7**. The kz15 razor
misses on every surviving frame by only **0.07–0.12 dec**. The
round's surprise, the **autocorrelation lag map**: only 5/44 rungs
(0 exceptions) have a positive root-scale cancellation deficit — on
39/44 worlds `(sum P)^2 > sum P^2`, the block sequence **reinforces**:
the c3 cancellation is already absorbed *into* the `P_j` magnitudes,
and the classical frames win by the `sqrt(m)`-count economy, not by
lag structure (large-lag mean `|rho|` 0.058). The L2 lemma candidate
after this round: **the vdC form as a provable chain statement**
(source-pure, window-independent) + the r269-known kz15 exact-finite
certificate for the razor. mp dps-60 wards; r272 regressions exact.
NO RH CLAIM.
r288 measured the **sign anatomy of the destructive coherence** (the
r285 `z_v = −3.15` address; `destructive_coherence_probe.py`, 31/31,
`PRIME.PORT.LSTAR.DESTRUCTIVE_COHERENCE.01` — experiments-side, no
status move): all E signs come from the mu-CD kernel `K_n` between nu
positions; phases measured against the zeros of the mu-orthonormal
`P_{N_w}` (Jacobi eigenvalues of the chain; Sturm counts, midpoint
alternation and mp dps-60 edge wards exact). **The carrier map**: the
sign field is band-structured — fold distance 1–2 is 87 pct
**positive** (in-phase sampling), distance 3–4 is 80 pct **negative**
(antiphase); the top mode at `N_w` rides the in-phase adjacent pairs
(band 1–2 = +0.68 of |T|, ARCH–ARCH +0.74) while the wall
destructivity in the source frame at the crossing (`X_v = −0.0517`,
`C_off = −0.105`) is carried by the **antiphase next-nearest (3–4
fold) ARCH–ARCH pairs**; the controls at their own crossing are ALL
constructive (`C_off` +0.30..+1.00) — a total sign reversal. The pure
sinc/CD sign formula explains the field only PARTIALLY (agreement
0.728 plain / 0.878 |E|-weighted, below the sealed 0.90 bar →
`SAMPLING_BLIND`); the sealed source-pure phase-dispersion separator
is **not found** (K_S1 world-blind, ensemble-generic; the plain
dispersion does not predict the wall, sp ≈ +0.2). The r287 object
comparison lands **`DIFFERENT_OBJECTS`** (the w9 block sequence
reinforces, canc = −4.5e-2 < 0, while the kernel-side wall
interference is destructive — two coordinates). The round's finding,
the **dose curve**: under conservation-faithful union jitter the
phase field barely rotates (turn rate 0.24/dose == the static
prediction; the chain zeros co-move) yet `z_v` flips destructive →
constructive already at dose 0.005 and the depth collapses — the
destructive coherence sits in hyper-fine alignments far below the
median phase resolution: **the r276 metric firewall is sharper than
any phase-rotation account** (measured; r276 P2_JIT quoted
COMPARISON_ONLY, different channel). Ensembles = the 28 r285
replicates, identical seeds, records re-gated; equidistribution
statistics typed MEASUREMENT_ONLY. NO L\* claim. NO RH CLAIM.
r289 ran the **arch-kernel diophantine relevance test** on exactly that
r288 finding (`arch_kernel_diophantine_probe.py`, 42/42,
`PRIME.PORT.LSTAR.ARCH_KERNEL_DIOPHANTINE_ANATOMY.01` —
experiments-side, no status move; strict reviewer contract: no proof
attack, **no Baker/Matveev application**, need documentation only) and
decided the sealed three-way question with the **arithmetic twin**
family. Verdict **`METRIC_ONLY`**: the **rational twin** — every tent
center `log n` replaced by a small-denominator rational multiple of
`Delta` (CF convergents, denominators ≤ 5.7e4, position cost ≤ 2.1e-9
= 1e-8 × the local gap; weights/cells/on-node family conserved
exactly) — **keeps the full signature identically** (minC 184,
crossing 185, `z_v` −3.149, carrier band 3–4 −0.105, AA −0.056; mp
dps-60 ward at the 1.7e-4 margin): every exact log-linear-form
relation is destroyed with **zero effect**; the **shuffle twin**
(tent-split fractions permuted, fraction *distribution*/cells/weights
preserved exactly, only the assignment destroyed) loses totally
(depth 0.21, `z_v` +4.4) but at effective metric dose 0.12 gap, where
the plain jitter ladder already loses — fully metric-explained. The
construction question is settled first: the lag vector factors
**exactly** through (cell, tent-split fraction, mass) — the fractions
`{log n/Delta} = {46 log2 n}` are the *only* sub-grid entry of the
source into the weights (completeness 4.2e-14; differences alone
underdetermine, LOUD; the folded weights are downstream-complete).
Q1 anatomy: `X_v` decomposes exactly over CD degrees; the
destructivity is carried by the **mid degree band 50–75 pct** (−3.17
z-units; carrier set −2.48), *not* the terminal degrees; at dose
0.005 with degree and frame fixed **no kernel band moves O(1)** — the
O(1) mover is the **crossing relocation itself** (pivot cascade 185 →
~131) while the phase stands still (0.0012). The measured **metric
threshold**: the coherence needs the fraction *profile* at ~1e-3 of
the local gap (ladder keeps at 1e-3, collapses at 3e-3) and is
indifferent to its number-theoretic exactness; the linear forms are
documented (28 exact 2-power resonances inside the window, min
nonzero `|{46 log2(n/m)}|` = 1.34e-6 at (181, 241), attribution tops
at the small primes 3, 7, 25, 8, 29, 17), and the need calculation
shows Baker/Matveev-class bounds ~7900× too weak in the exponent for
the unbounded family — **Baker is unnecessary on all tested scales**;
the sub-gap information the wall reads is metric, not diophantine.
NO L\* claim. NO RH CLAIM.
r290 mapped the **geometry of the working set in profile space**
(`profile_functional_probe.py`, 32/32, SPEC_SHA `f953dd71`,
`PRIME.PORT.LSTAR.PROFILE_FUNCTIONAL.01` — experiments-side, no
status move): profile = the signed lag density on the fixed grid
(the r289 completeness object; linear in the tent-split fractions
at fixed cells), distances = gap-equivalent doses `theta_eq` in the
exact lag coordinate (analytic reference `0.5 Σ m g/Δ`; amendment a1,
demanded by the sealed self-consistency gate). **The onset map**
(five sealed weight-profile interpolation paths, NEAR = 0.90 bar):
world-directed axes kill at `theta_eq` 4e-5..2.3e-4 — **5..50×
below** the 1e-3..3e-3 jitter threshold — the rim is strongly
**anisotropic**; every onset is **GRADUAL** (soft shoulder, s 0.95 →
0.77..0.89; curvature |kappa| ≤ 0.08 — no cliff). **Isotropy**: a
thin **TUBE** (killfrac 1.00 at 3e-3 over 16 conservation-gated
random directions; NEAR-radius ~5e-4..2e-3, consistent with the
r285 ensemble pct 0.00). The **SMOOTH (arithmetic) direction** is a
privileged killer axis (onset 1.36e-4) yet **orthogonal** to the
first-order wall gradient (cos −3e-5): its lethality is collective,
not gradient-linear. The **r280 ridge axis** is real geometry: the
rebuilt OPT direction (theta_up 3.873e-5 anchor exact) **lifts** the
wall (minC 185) at extension factors 1..8 = `theta_eq` up to 1.2e-3
— inside a tube where half the random directions already kill — and
only dies at 2.4e-3. The **functional contest** (four sealed
source-pure profile scalars: antiphase 3–4 autocorrelation = the
r288 carrier band, total variation, perturbative r278 gradient
alignment, mid-band 50–75 pct deviation energy = the r289 flip
carrier) ends **`ALL_FUNCTIONALS_BLIND`**: best spearman +0.263 over
the 187-point test set (bars 0.5/0.8), none separates the worlds
non-trivially — the coherence functional, if it exists as a closed
profile scalar, is **not among these classes**; the working set
remains implicitly characterized. Channel identity gated bitwise on
MAIN; r288/r289 records, the r289 jitter-ladder anchor (exact
seeds), the r285 control flips and all 28 replicates (ENS_SIGN as
exact density sign-flips) reproduced. NO L\* claim. NO RH CLAIM.
r291 dissected the **anatomy of the one known raising axis**
(`ridge_anatomy_probe.py`, 30/30, SPEC_SHA `bb512c17`,
`PRIME.PORT.LSTAR.RIDGE_ANATOMY.01` — experiments-side, no status
move; amendments a1–a3 disclosed: r290 pinned metric quadruple,
1e-9 conservation-gate headroom, matched-dose/global reporting
split). **Ridge section**: the lift is **not** a distinguished
sparse support, band, or hidden functional — near the threshold it
is a **first-order budget phenomenon**: a conservation-gated
sub-direction lifts exactly when its side-selected linearized wall
budget `-Σ c_j·dose` exceeds one sharp scalar threshold in
**(1.280, 1.291]** (~1.3× second-order resistance over the naive
flip level −1; perfect separation over all 18 matched-dose subset
cases: PRIME/HEAD/XIPOS lift, POWER/MID/TAIL/XINEG do not; k_min =
9 top-budget atoms n = 2, 3, 5, 13, 11, 4, 29, 7, 89 — small-prime
head-heavy, not one-atom), with exactly one **overdrive
retraction** deep in the over-driven few-atom regime (TOP6 at
factor 8, margin 9.09, back to minC 184 while fully alive).
**Fixpoint question**: the step-controlled r280 recipe saturates
at a **one-degree plateau** (minC 185 from iterate 1, never 186;
axis decoherence cos +0.92/+1.00/+0.43 → negative; every full step
rejected from iterate 3 on; final `theta_eq` 1.43e-3 = the tube
rim) — **`RIDGE_NO_FIXPOINT`**, no stationary better world; and
the lift is **`LIFT_MAIN_SPECIFIC`**: from EPSTEIN the same
iteration terminates at step 0 (first-order **flat** wall,
theta_up 5.03e10) and the matched-dose ladder 7.75e-5..0.1 never
lifts. **Non-local functional contest** (sealed, fresh 143-point
conservation-gated corpus, 103-point test split): ridge projection
(best, sp +0.471), rank-2 lethality Gram (train/test split-sealed,
+0.137) and antiphase deviation pair correlation (−0.463) all stay
below the trivial size baseline |sp| 0.881 —
**`NONLOCAL_BLIND`**, the working set remains implicitly
characterized. **SMOOTH anatomy**:
**`SMOOTH_COLLECTIVE_2ND_ORDER`** — the killer axis is delocalized
(PR 112 of dim 368) and its margin curvature is Richardson-stable
−23.3 vs −0.98 for a same-length random direction (ratio 23.7 ≥
10): a genuine **quadratic wall valley**, not an r259-style
higher-order resummation gap. NO L\* claim. NO RH CLAIM.
r292 measured the **curvature two-form** of the working set
(`curvature_form_probe.py`, 36/36, SPEC_SHA `050821ff`,
`PRIME.PORT.LSTAR.CURVATURE_FORM.01` — experiments-side, no status
move; amendments a1–a4 disclosed: Hessian step ladder halved at the
measured margin-valid scale, one-sided estimator for the
side-selected budget kink, fractional-dose gate pad 1e-7,
L2-unit-consistent ridge display). **Hessian spectroscopy** (29
sealed theta_eq-normalized directions, 3-step diagonals + the full
406-pair polarization matrix, expansion-identity crosscheck 15/15):
**`HESSIAN_LOW_RANK(1)`** in the sealed L2 spectral metric — one
eigendirection carries 92.5 pct of the curvature mass, a pure
on-support **DENS combination**, *not* the SMOOTH axis (|cos| 0.07)
and *not* the ridge (**`RIDGE_CURV_FLAT`**, L2-rank 28/29: the lift
axis sits in the flattest curvature sector, the r291 expectation);
all 29 diagonals negative — a genuine multi-directional quadratic
valley. **EPSTEIN contrast**: the first-order flat wall is *also*
second-order structureless (ratio 5.4e-15) — MAIN's wall has
curvature structure, EPSTEIN's has none. **The threshold object
m\***: the measured log|h_184| quadratic term along the 18 r291
matched-dose sub-directions is **negative on all 18** (one-sided
cubic-triple solve, linear identity vs the analytic budget 0.3 pct)
— a second-order *assist*, not resistance: the corrected budgets
still separate perfectly but the bracket midpoint moves to 1.761,
*away* from the naive flip level —
**`THRESHOLD_NOT_EXPLAINED`** + **`RETRACTION_HIGHER_ORDER`**
(B2/B3 over-predict the TOP6 factor-8 retraction): the ~1.3×
factor and the retraction are **near-flip nonlinearities invisible
to low-order jets**; m\* stays an emergent threshold object.
**Two-form functional contest** (split-sealed, disjoint 94-point
test corpus): F10 curvature energy sp **+0.884** vs the size
baseline 0.907 — **`CURVATURE_BLIND`**, the fourth honest negative
under the sealed rule, but the closest miss of the whole program
(r290 best 0.263, r291 best 0.471; AUC(dead) 0.097): the
working-set geometry is predominantly second-order visible, the
open edge is the theta_eq-vs-L2 metric reconciliation. NO L\*
claim. NO RH CLAIM.
r293 adjudicated that open edge — the **metric reconciliation**
(`metric_reconciliation_probe.py`, 43/43, SPEC_SHA `33c44cc6`,
`PRIME.PORT.LSTAR.METRIC_RECONCILIATION.01` — experiments-side, no
status move; amendment a1 disclosed: the Leg-C death watch redrawn
from the raw h_185 sign to the real criterion minC < 184, a
diagnostic display fix). **Leg A sealed three metrics before any
functional evaluation** (M1 theta_eq, M2 L2 on the signed density,
M3 the curvature metric \|delta\|² = x|H_tr|x with the
spectral-absolute-value rule; distortion factors max/min 11.0 /
59.5 / 123.5, largest-deviation sector DENS = the r292
top-eigenaxis sector). **The hammer landed**:
**`METRIC_RECONCILED(F10 in M2)`** — evaluated against the size
baseline *of its own metric*, the split-sealed curvature energy
F10 beats a baseline for the first time in four rounds (sp
**+0.884** vs sp(F0_M2) −0.860) *and* holds
**`FUNCTIONAL_BEYOND_BASELINE`** (partial spearman sp(F10, s |
F0_M2) = **+0.423** on the 94-point test split, sign-stable +0.826
on the train-side replica) — the first predictive closed
functional of the working set; **promotion candidate flagged** for
the next consolidation wave (nothing promoted).
**`MIX_IS_CAUSE`**: the r292 mixed margin −0.023 vs the reconciled
home margin +0.024 — the metric mix *was* the cause of the r292
near-miss. **m\* as a singularity object**: all 8 selected flips
are **simple zeros** of h_184 (alpha = 1.000 ± 0.003 both sides;
every flip lands at minC 185 with h_185 < 0) but
**`MSTAR_NO_LAW`**: the per-direction thresholds m\*_dir = b·f\*
range 1.13–1.39 — no fixed safety distance to the singularity in
budget, theta_eq, L2 or curvature units (spreads 3.0/3.0/3.8/10.2
vs bar 1.25): the r291 "one scalar m\*" bracket was a fixed-dose
artifact, the ~1.3 factor is typical, not universal. The TOP6
factor-8 retraction **is a second simple h_184 zero** at f_ret =
7.107 — on the side the monotone budget prognosis is structurally
blind to. **Baseline anatomy**: `F0_DIRECTIONAL` (median
within-band |sp| 0.370 > 0.3, dead-only −0.890 / alive-only
−0.005) — the baseline is not a pure radius artifact, the
r290–r292 bars stand as drawn. NO L\* claim. NO RH CLAIM.
r294 stress-tested exactly that promotion candidate **before** any
consolidation — the **F10 stability round**
(`f10_stability_probe.py`, 42/42, SPEC_SHA `88c6fd1e`,
`PRIME.PORT.LSTAR.F10_STABILITY.01` — experiments-side, no status
move). **Leg A, the sealed promotion bar** (fixed before any corpus
was built): five fresh, seed-disjoint 147-point corpora (same
conservation discipline and r293 test-family mixture; tag
disjointness gated pairwise and against the r292 training tags)
against the **unchanged, hash-sealed r293 `H_tr`** (re-training
mutant CAUGHT by the SHA gate). The **|sp| win replicates on every
corpus** (5/5: F10 +0.675..+0.787 vs F0_M2 −0.660..−0.714, margins
+0.003..+0.115) — F10 is *not* a split/corpus artifact as a
predictor — **but the partial-correlation channel does not
replicate at its r293 magnitude** (partials +0.248..+0.555, median
0.299; only 2/5 clear the sealed +0.3 bar, two miss by ≤ 0.003):
sealed verdict **`F10_FRAGILE`** (win 5/5, part 2/5, full 2/5) —
**promotion precondition NOT met**, recommendation NO for wave 8;
the r293 partial +0.423 was the *top* of the fresh distribution,
not the center. **Leg B robustness**: `TRAIN_ROBUST` (leave-3-out
jackknife sigma 0.0101 ≤ 0.012) + **`RANK2_CARRIES`** (the top DENS
eigenaxis alone reaches +0.855, two whitened eigendirections beat
the baseline; top-axis DENS share 0.989) + `STEP_STABLE` (half FD
ladder sp 0.9995; the double ladder MARGIN_INVALID — the r292-a1
NaN boundary reproduced and excluded by the sealed rule). **Leg C
window transport** (w7/w11, complete window-own chains: walls
59/63, extraction dose caps 0.157/0.068 after the measured w11
full-ridge collapse, 67-point window corpora): **knife-edge** —
w11 beats its own L2 baseline by +0.002, w7 loses by 0.009 →
**`WINDOW_PARTIAL(w11)`**; the solid transport finding is
**mechanism constancy**: the top curvature eigenaxis is a DENS
combination on w7/w9/w11 alike (shares 0.795/0.989/0.834). **Leg
D**: **`L2_NOT_DENS`** — the DENS-neutralized baseline (scales
2.2..75.6) is *weaker* than plain L2 (−0.828 vs −0.860) and the
partial against it rises to +0.586: the F10 information is not the
DENS-sector length distortion; *why* L2 is the right coordinate
stays open. NO promotion. NO L\* claim. NO RH CLAIM.
r295 executed the two r294 fronts — the **F10 sp-hardening round**
(`f10_sp_hardening_probe.py`, 44/44, SPEC_SHA `4d7d8095`,
`PRIME.PORT.LSTAR.F10_SP_HARDENING.01` — experiments-side, no
status move; no amendments after the freeze). **Leg A, the sealed
three-clause hardening bar** (fixed before any corpus was built,
untouchable): **twenty** fresh seed-disjoint 147-point corpora
(seeds 300000+1000k; the full 192-seed forbidden set of
r292/r293/r294 enumerated and gated) against the unchanged,
hash-gated `H_tr` (published seal `3447ed198a56` re-verified). The
result is an **honest negative**: F10 wins the |sp| contest only
on **14/20** corpora (margins −0.079..+0.125, median +0.028, six
losses) → sealed verdict **`F10_SP_MAJORITY`** — the HARDENED bar
(≥ 18/20 AND margin median ≥ +0.02 AND no loss beyond 0.02) fails
on two clauses: **the r294 5/5 was itself top-of-distribution**,
exactly the way the r293 partial was; the partial-free ranking
statement is a *documented regularity* (19/25 combined census),
not a theorem candidate — **promotion recommendation NO for wave
8**, the candidate statement (Leg C) VOID by the sealed rule.
PARTIAL20 median +0.279, IQR [+0.182, +0.401] (the r293 +0.423
confirmed top-of-distribution; **`PARTIAL_STD`** — the partial at
the fixed 147-point mix — sealed as the standardized statistic of
future rounds). **Leg B partial anatomy** (pooled 2921 points): no
direction family reaches the STRONG 0.4 bar (PATH +0.156 / WORLD
+0.245 / FRAC +0.067 NULL / DENS +0.104) — the beyond-distance
information is diffuse, strongest on structured deviations toward
dead worlds, near-null on random FRAC rays; the
r293-composition-matched subsample (the PATH-heavy r293 mix) lifts
the median partial to +0.346 — the sealed gain clause passes
(+0.067 ≥ +0.05) but the level clause fails by 0.004 →
**`R293_LUCK`**: composition explains *part* of the r293/r294
discrepancy, the rest was n=94 sampling noise. **Leg D**:
`L2_VIA_CONSERVATION` formally (the eta_0 projection leaves the
contest exactly invariant) but near-tautological (projected-out
share median 2.9e-30) — *why* L2 stays open in substance.
Must-fails m1 hash / m2 seed collision / m3 bar sharpness (every
clause has teeth, both directions) / m4 conservation all CAUGHT.
NO promotion. NO L\* claim. NO RH CLAIM.
**Wave 6 (2026-08-25)** froze the L\* arc (rounds r282–r285 plus the
standalone problem document) as
**`v963_lstar_reduction_dictionary.py`** (ledger
`PRIME.LSTAR.REDUCTION_DICTIONARY.01` [E]; the four probes embedded
byte-exact in the sealed `--smoke` stage, plus a module-own exact
section — Hankel determinants, frame congruence minors, full
Cauchy–Binet gauge exhaustion, dual chains in pure Fractions): the
r283 A2 chain as exact gates ⇒ **L\* is the canonical reduction of
the open center** (stage 8/11 of the chain above), the r282
four-language elimination as named negative gates
(`CONTEST_ALL_DEAD`), the r285 decomposition bookkeeping exact, and
the r284/r285 anatomy + the honest DCXX margin-decay finding as typed
measurement pins. The open center itself is re-registered in
canonical form as **`PRIME.LSTAR.SUBORDINATION.01` [O]** (the
predecessor `PRIME.PORT.RHP.FULLSOURCE.QUASIDEFINITENESS.01` and
`PRIME.PORT.REPRESENTATION.CONTEST.01` are consumed). In Lean the
canonical form is stated (`lstar_subordination`, the wave-6 `sorry`)
and the direction L\* ⇒ free-window positivity is **proved**
(`lstar_implies_free_window` via the quadratic-form dictionary
`hankel_quadform`; build green, 9 intentional `sorry`s). The r286
margin-scaling round and the r287 L2 census (above) were sealed
experiments-side at that cut and are consumed in wave 7 (below).
NO RH CLAIM.
**Wave 7 (2026-08-25)** froze the coherence arc (rounds r286–r289) as
**`v964_lstar_coherence_census.py`** (ledger
`PRIME.LSTAR.COHERENCE_CENSUS.01` [E]; the four probes embedded
byte-exact in the sealed `--smoke` stage — SPEC SHAs
`0a44ac4e`/`761d88fa`/`da46f7ee`/`91cdc2b1` — plus a module-own exact
section: the Fejér window-sum identity, the van der Corput inequality
with its Cauchy–Schwarz proof route, the Abel/Erdős–Turán bound, the
mutant must-fails, the 2-power resonance count C(8,2) = 28, in pure
Fractions): the DCXX margin warning **resolved harmless** (r286: 15
new anchors N_w 942–1218 all mp-sign-safe positive, O(1) offset
survives, flattening power law, no counterexample), the generic L2
half **classical by a named theorem** (r287: exact vdC at
H = ceil(sqrt(m)), delta' +0.309 > 0.21 world-blind, 6/7 + 38/42 —
new open contract **`PRIME.PORT.L2.VDC_LEMMA.01` [O]**: the chain
origin of the P-variance scaling), the destructive coherence with a
**named carrier class** (r288: antiphase next-nearest ARCH–ARCH
pairs, z_v = −3.149, total control reversal, finest alignments below
phase resolution), and the diophantine route **excluded** (r289
`METRIC_ONLY`: the rational twin keeps the full signature identical;
metric coherence threshold 1e-3..3e-3 of the local gap; Baker ~7900x
too weak and unnecessary) — the open front is the **profile
functional** (`PRIME.LSTAR.SUBORDINATION.01` [O] carried at the
wave-7 state; r290 above, in flight at this cut, probes exactly this
question).  In Lean no new statement is needed — the hole stays
`lstar_subordination` (docstring refreshed to the wave-7 state).
NO RH CLAIM.

**Wave 8 (2026-08-26)** froze the curvature arc (rounds r290–r295) as
**`v965_lstar_curvature_arc.py`** (ledger
`PRIME.LSTAR.CURVATURE_ARC.01` [E]; the six probes embedded byte-exact
in the sealed `--smoke` stage — SPEC SHAs
`f953dd71`/`bb512c17`/`050821ff`/`33c44cc6`/`88c6fd1e`/`4d7d8095` —
plus a module-own exact section S0: the r292 polarization identity
H(u,v) = [d2(u+v) − d2(u−v)]/4 with its wrong-prefactor must-fail, the
r291 budget telescope (exact partition additivity) with the
dropped-atom must-fail, the r293 alpha = 1 simple-zero doubling law
with the double-zero mutant, and the sealed r294/r295 decision bars as
exact Fractions logic with tipping must-fails).  The wave freezes the
**measured geometry of the working set** — a soft-shouldered
anisotropic tube (killfrac 0.38/0.62/1.00 at 5e-4/1e-3/2e-3
gap-equivalent; world axes kill 5–50× earlier; the SMOOTH axis a
privileged yet gradient-orthogonal killer), a real one-degree ridge
(plateau minC 185; the lift a first-order budget phenomenon with the
threshold in (1.280, 1.291]; MAIN-specific), a rank-1 DENS curvature
valley (92.5 pct, lam_top −0.418; not SMOOTH, not the ridge), and
simple h_184 flip zeros (alpha = 1; the TOP6 retraction a second
simple zero at f_ret 7.107) — **together with the honest negatives of
the closed-functional search, verbatim**: `ALL_FUNCTIONALS_BLIND`
over five sealed class families, r294 `F10_FRAGILE` (win 5/5, part
2/5) and r295 `F10_SP_MAJORITY` (14/20) — **F10 is NOT promoted**,
the sealed bars rejected it twice and the rejections are the frozen
content; `MSTAR_NO_LAW`; `R293_LUCK` (the r293 partial and the r294
5/5 were both top-of-distribution).  The characterization question is
re-registered as **`PRIME.LSTAR.CLOSED_FUNCTIONAL.01` [O]** (open
fronts: loss-corpus forensics K05/K07/K19, the rank-2 DENS core as a
sharper candidate, why L2 is substantially the right coordinate);
`PRIME.LSTAR.SUBORDINATION.01` [O] carried at the wave-8 state (no
round in flight at this cut).  In Lean no new statement is needed —
the wave carries no new formal target form; the hole stays
`lstar_subordination`.  NO RH CLAIM.

r296 ran the reviewer's **hard identity fork** on the r292 axis — the
**DENS-identity round** (`dens_identity_probe.py`, 39/39, SPEC_SHA
`ffb413c8`; `PRIME.PORT.LSTAR.DENS_IDENTITY.01` — experiments-side;
consumed in wave 9, `v966`): is the 92.5-pct DENS top
curvature eigendirection *mathematically identifiable* and *coupled to
the L\* extremizer*?  One hash-sealed candidate library, frozen before
any cosine: **T1 = grad lambda_max(E_184)** via Hellmann–Feynman from
the extremal eigenvector (FD ward 3.8e-7; lam2/lam3 gradients as
excluded controls), the von-Mangoldt family, local-metric family,
kernel diagonals, the 4-dim moment-gradient subspace, 8 nulls + decoy;
primary metric = the in-span cosine on the 29-direction span (e_top
exists only there), honesty columns raw/lag/capture.  **Sealed verdict:
`DENS_WORLD_BLIND`** (base `DENS_PARTIAL(T5_MOMENTS, 0.825)` downgraded
by the world-contrast rule: the moment-subspace share *holds* at
EPST/SCR 0.57/0.60 ≥ 0.40 — construction-adjacent low-moment
smoothness, with the disclosed tension that the absolute 0.40 collapse
bar is conservative for subspace members).  **The coupling number of
the round: cos(e_top, grad lambda_max) = +0.394** — below the sealed
0.40 bar (noise max 0.331), and the lam3 control couples *harder*
(−0.574): what coupling exists is a near-crossing eigenvalue-*band*
property, not top-specific.  Every arithmetic candidate stays ≤ 0.38;
raw-grid overlaps ≤ 0.08 at captures ≤ 0.36 — the axis is dominated by
structure no on-support arithmetic profile in the library represents.
Window transport w7/w9/w11 measured (window DENS shares 0.795/0.834 ==
r294); amendments a1–a4 disclosed pre-freeze.  NO L\* claim, NO RH
CLAIM.

**Routing note (r296 → r297, pre-adjudicated reviewer fork):** with
`DENS_WORLD_BLIND` the DENS-identity lane is **closed** (coupling
number +0.394 below the sealed bar, no arithmetic candidate above
0.38) — resources move to the **L2 front**, where the exact vdC
theorem (v964-S0) already waits with a 0.24 exponent pad (need
> 0.21, truth ~0.45): **L2 is the active front**.

r297 opened the L2 lane with the **vdC chain-provenance round**
(`vdc_chain_provenance_probe.py`, 27/27, SPEC_SHA `e42a76eb`;
`PRIME.PORT.L2.VDC_CHAIN_PROVENANCE.01` — experiments-side; consumed
in wave 9, `v966`): can the vdC *input* (the P-variance /
autocorrelation scaling) be derived source-purely from the chain?
**Leg A freezes the target inequality** — the weakest sufficient input
bound, all constants live: `sigma := slope(S_F/M²) ≤ sigma* =
2·(slope(eps_c2/M) − 0.21) − slope(pref) = −0.516`; measured sigma
−0.714 — the truth clears the target with 0.198 sigma margin (0.099 in
delta'): *what is missing is not room, it is provenance* (pad
arithmetic exact; the pad-dropped mutant composes to delta' = 0).
Three source-pure routes, every chain step exact on 47 worlds
(world-blind by the same algebra): **B1 pair-identity** delta' −0.112
(break: the MAX pair gap falls at slope −0.535 where ≤ −1.178 would
suffice — the Jacobi near-balance is a *mean* statement); **B2
node-density** delta' −0.097 (the run-boundary *positions* are
near-equidistributed with falling discrepancy D_rank med 0.024, slope
−0.42 — provable terrain — but the *mass* imbalance grows, slope
+0.244); **B3 Fejér/Parseval** delta' −0.097, with **the surprise of
the round: the chain is exactly orthogonal w.r.t. the WINDOW measure**
(cross/diag devs 0.000; BORDER and WDIFF break) — the Parseval sum
rule (proved symbolically in Fractions) *attaches* there; what remains
is the measure transfer (the vdC input lives on the border union) and
the signed-to-absolute gap (med 4.3, growing).  **Sealed verdict:
`CHAIN_PROVENANCE_PARTIAL(B3_PARSEVAL, slope miss 0.307)`** — the
common obstruction named: every universal majorization of A(0) through
magnitudes pays a growing max/mean imbalance factor; the *signed*
structure carrying sigma = −0.71 is what magnitude chains cannot see.
kz15 structurally outside all three routes (r270 exact-finite stays
the permanent closure).  The two assets left for the lane: the exact
window-measure orthogonality (a bordered-form transfer statement is
the concrete named object) and the falling rank-coordinate
equidistribution.  NO L2-lemma claim, NO RH CLAIM.

r298 executed the named asset-1 attachment point with the
**window-border transfer round** (`window_border_transfer_probe.py`,
28/28, SPEC_SHA `05e831be`;
`PRIME.PORT.L2.WINDOW_BORDER_TRANSFER.01` — experiments-side;
consumed in wave 9, `v966`): does an exact, *sign-
preserving* transfer statement exist from the window measure (where
the chain is exactly orthogonal) to the border form (where the vdC
input lives)?  **It does, and it is exact:** with `Delta := border ⊖
window` and the frozen positional Fejér block kernel, the vdC input is
the quadratic form `S_F = B(beta,beta)` and splits per-rung exactly
(dev ≤ 8.8e-16 of scale, 47 worlds) as **`S_F = B(omega,omega) +
B(Delta, omega+beta)` = MAIN + T(Δ)** — nothing passes through
magnitudes (the `|Delta|` mutant breaks the identity by a factor ~1.7
of scale, LOUD: the r297 magnitude error certified as a must-fail);
the *linear* corollary is measured at float precision: the window
image of the drive functional **vanishes** — the drive is 100 pct
transfer at the linear level.  **Sealed verdict: `TRANSFER_DOMINANT`**
— the window main term is *empty* (med −3.94 dec, two decades below
S_F, slope −1.386 collapsing) and the omega-Delta cross term
negligible (med −1.4e-4), so **`S_F ≈ B(PDelta, PDelta)`: the vdC
input IS the Fejér energy of the difference measure itself**; sigma =
−0.714 is the decay exponent of this one named source-pure quantity;
the in-T sign cancellation (med 1.69) grows at +0.207 — the r297
imbalance factor localized *inside* T.  Leg-D bycatch: the imbalance
growth is broad-based (no outlier carrier), and its carriers coincide
with the T carriers (sp med +0.69) — **the two r297 gap halves are one
structure**.  Candidate theorem VOID (not subdominant); the precisely
narrowed wave-9 gap object: *prove the Fejér-energy decay of Delta*.
World-blind (EPST/SCR same class); anchors bit-near; amendments a1/a2
disclosed.  NO L2-lemma claim, NO RH CLAIM.

r299 executed that gap object with the **Fejér decay round**
(`fejer_decay_probe.py`, 32/32, SPEC_SHA `f432e944`;
`PRIME.PORT.L2.FEJER_DECAY.01` — experiments-side; consumed in wave
9, `v966`): can the decay of `B(PDelta,PDelta)` (sigma ≤
sigma* = −0.516) be derived source-purely from position
equidistribution plus a controllable signed-mass statement?  The
module-own **exact spectral representation** `B = (1/L) sum_k
F_H(theta_k) |Dhat(k)|²` (recomposition/Parseval/majorization wards ≤
7.1e-16 on 47 worlds) answers the anatomy question: **LOWPASS** — the
Δ energy sits inside the Fejér main lobe (band shares 0.93/0.04/0.02,
q50 = 0.19 lobes): sigma is a low-frequency phenomenon.  The
structural surprise: the border and window unions share their **full
support on 42/42 rungs** — `Delta_fresh ≡ 0`, **Δ is a pure c-value
difference measure** on one shared node set, and the relative
c-difference does *not* fall (cconv med 0.86, slope +0.045): the decay
is aggregate cancellation of a non-converging difference profile.
**Sealed verdict: `DECAY_SPLIT`** — the Erdős–Turán/Abel composition
fails loudly (composed slope +1.948 vs sigma* −0.516; the
position/kernel factor grows +1.504: the r297 magnitude wall recurs at
the frequency level; `MASS_TARGET` missed by 2.46), but the **pair
split fires**: the diagonal `sum PDelta_j²` already falls fast enough
(slope −0.571 ≤ −0.516, margin 0.055) and the pair ratio B/D (med
1.29) *falls* (−0.168) — the frozen rest pair is **`DIAG_TARGET`** (a
magnitude-density bound on one positive quantity) **+ ratio flatness**
(the off-diagonal reinforcement does not grow), with `CVALUE_TARGET`
frozen on the B3 side; the conditional candidate theorem (Leg-B bound
+ rest hypotheses ⇒ sigma ≤ sigma* ⇒ r297 target ⇒ v964-S0 vdC ⇒
delta' > 0.21 generic half) is printed, wave-9 candidate NOT promoted.
Bycatch: the in-T cancellation gap is entirely low-band and *not*
carried by the mass blocks (sp −0.01 vs the +0.69 energy coupling);
and the O-sign class separates MAIN (pair field *reinforces*, O < 0 on
13/42 only) from both dead controls (O < 0) — the first
world-separating sign class of the L2 lane, disclosed as
`WORLD_SENSITIVE`.  Anchors bit-near; Fractions ET section exact with
the prefactor-mutant family CAUGHT (v964-S0 T4 anchor); amendments
a1/a2 disclosed.  NO L2-lemma claim, NO RH CLAIM.

r300 executed the frozen rest pair with the **diag-target round**
(`diag_target_probe.py`, 31/31, SPEC_SHA `55218b5d`;
`PRIME.PORT.L2.DIAG_TARGET.01` — experiments-side; consumed in wave
9, `v966`): can the diagonal decay `sl_D ≤ sigma* = −0.516`
(margin 0.055, thin) be derived source-purely — and the ratio flatness
with it?  **The anatomy relocates the target:** in the exact
participation coordinates of the one positive vector `|PDelta|`
(`n_eff = L1²/D`, `fill = D/(mx·L1)`), the diagonal decay decomposes
*exactly* (dev 6.7e-16) as **`sl_D = 2·sl_L1 − sl_neff = 2·(+0.196) −
(+0.963)`** — the effective-carrier count grows at ~N (sp(N, n_eff)
+0.96) against a mildly growing L1 mass: "many small instead of few
large" is now a measured identity.  Both magnitude derivation routes
fail exactly the way r297 predicted: the B1 chain-norm Cauchy–Schwarz
on the `|dw|` measure composes to −0.384 (misses σ* by 0.133) *and*
the `|dw|` identity census BREAKS (cross med 0.932) — the proven
window sum rule does **not** attach to the difference measure; the B2
`max × L1` factorization composes to −0.346 (misses by 0.170) even
though the max itself is fully atom-controlled (a single c-value
distance, `mx/maxatom` med 1.07): the irreducible loss is the fill
decay −0.225, invisible to any max × mass bound.  **Sealed verdicts:
`DIAG_SPLIT(B3/NEFF_TARGET)` + `RATIO_BOUNDED_STRUCTURAL`** — the
ratio half of the r299 rest pair is *settled at the structural level*:
the exact kernel envelope `F_H(θ) ≤ min(H, 1/(H sin²(θ/2)))` majorizes
`B/D ≤ R_env` (med 1.61) and R_env itself falls (−0.122); the r299
lobe-width heuristic is refuted (sp(B/D, q50) = −0.47).  The refined
freeze is smaller than r299's: **one inequality remains — `NEFF_TARGET`:
prove `slope(n_eff) ≥ 2·sl_L1 − sigma* = +0.908`** (measured +0.963,
margin 0.055; disclosed as an exact reparametrization of DIAG_TARGET),
with the falling boundary discrepancy (D_rank(Δ) slope −0.117,
sp(D_rank, n_eff) = −0.81) as the named bridge to provable
equidistribution terrain.  Bycatch: D itself is BROADBAND
(0.27/0.26/0.46) — the low-pass character of sigma is pure Fejér
weighting; and the FILL class separates MAIN (FILL_LOW) from both dead
controls (FILL_HIGH) — the second world-separating class of the lane,
disclosed as `WORLD_SENSITIVE`.  Anchors bit-near (r297 + r298 + full
r299 set); amendment a1 (scale floor on the degenerate SMOOTH Δ ≡ 0)
disclosed.  NO L2-lemma claim, NO RH CLAIM.

**Wave 9 (2026-08-26)** froze the DENS fork and the L2 reduction chain
(rounds r296–r300) as **`v966_l2_reduction_chain.py`** (ledger
`PRIME.L2.REDUCTION_CHAIN.01` [E]; the five probes embedded byte-exact
in the sealed `--smoke` stage — SPEC SHAs
`ffb413c8`/`e42a76eb`/`05e831be`/`f432e944`/`55218b5d` — plus a
module-own exact section S0: the r297 σ\*-composition arithmetic
(δ′(σ\*(pad)) == pad exactly, the pad-dropped mutant composes to 0 —
the r297-m2 sharpness as an identity), the r298 decomposition identity
`S_F = B(omega,omega) + B(Delta,omega+beta)` on the exact positional
Fejér block form (toy 3 = 1 + 2; the |Δ| and wrong-weight mutants
break by exactly 4 and 2/3), the r300 participation identity
`D·n_eff == L1²` with ratio multiplicativity and the machine-checked
DIAG↔NEFF equivalence (margin 0.055 exact in both coordinates; the
signed-sum mutant breaks by 8/3), the r300 kernel envelope
`F_H(s) ≤ min(H, 1/(H·s))` exact-rational with equality witnesses on
both branches, and the five sealed verdict bars as exact decision
logic with tipping mutants).  The wave freezes **two things**: the
DENS-identity lane is **closed honestly** (r296 `DENS_WORLD_BLIND`:
coupling number +0.394 below the sealed 0.40 bar, λ3 couples harder
−0.574, every arithmetic candidate ≤ 0.38, the moment subspace
construction-adjacent — the lane routed to L2 per the pre-adjudicated
reviewer fork), and the **L2 chain stands reduced to ONE inequality**:
`[NEFF_TARGET, open] ⟹ sl_D ≤ σ* (exact participation algebra) ⟹
with RATIO_BOUNDED_STRUCTURAL sigma ≤ σ* ⟹ the r297 target
(margin 0.198) ⟹ the v964-S0 vdC theorem ⟹ delta' > 0.21 on the
generic half` (exceptions: 6 via r287 F2, kz15 exact-finite via
r270).  The magnitude no-go catalog is frozen verbatim (r297 B1/B2,
the r299 ET/Abel wall +1.948, the r300 chain-norm and max×L1 closures
with the invisible fill decay).  The one remaining inequality is
registered as **`PRIME.L2.NEFF_TARGET.01` [O]** (slope(n_eff) ≥
+0.908, measured +0.963, margin 0.055; the D_rank bridge sp −0.81 is
correlational, not causal; **r301 in flight, not consumed**).  In
Lean **no new statement** is added: NEFF_TARGET is a measured
ladder-slope aggregate (a halves log-slope estimator over 42 rungs
with a threshold derived from measured slopes), not an exact finite
identity and not a matrix theorem, and a universally quantified
ladder form would be refutable exactly like the pre-r273 universal
edge forms (the machine-checked guard pattern of
`RH/Counterexamples.lean`) — the hole stays `lstar_subordination`.
NO RH CLAIM.

r301 executed the wave-9 open inequality with the **neff-target
round** (`neff_target_probe.py`, 32/32, SPEC_SHA `6f8cc404`;
`PRIME.PORT.L2.NEFF_TARGET.01` — experiments-side; consumed by `v967`, wave 10): *why* does `n_eff = L1²/D` grow ~N, and can
`slope(n_eff) ≥ +0.908` (measured +0.963, margin 0.055) be derived
source-purely?  **The anatomy:** the whole Rényi order family
`N_2/N_3/N_4/N_inf` (med 37.41/27.88/24.04/15.47) grows at nearly one
exponent (+0.963/+0.926/+0.894/+0.738, tail slope only −0.069) —
**echt anti-concentration**, not a broadening tail; the positional
footprint is flat.  **The stability:** the thin margin is
jackknife-STABLE (min +0.936 ≥ +0.908, 0/42 leave-one-out ladders
below NEED — no `MARGIN_FRAGILE`), while the half-ladder census shows
the deep half at +0.802: the honest anti-extrapolation flag stands.
**The structural finding (B1):** `n_eff = n_act/(1 + CV²)` EXACTLY
(module-own Fractions re-proof, ward ≤ 6.8e-16 on 47 worlds) with the
PERFECT count link **`n_act == m` on 42/42** — the effective-carrier
count IS the constructive level-2 block count, growing at +1.002; the
entire NEFF_TARGET compresses into **`UNIF_TARGET`: slope(1 + CV²) ≤
sl_nact − NEED = +0.094** (measured +0.039, margin 0.055) — a growth
statement became a *boundedness* statement about the block profile of
`|PDelta|`; `NEFF_DERIVED` missed only the sealed CV-flat clause
(+0.039 > 0).  **Closed honestly (B2):** the exact
weighted-discrepancy interval chain (`|w_j/W − λ_j| ≤ 2·δ_w`, a
star-discrepancy theorem, wards exact on 47 worlds) fails by 0.834 —
the `|dc|`-WEIGHTED discrepancy (med 0.167, slope −0.017) is 10× the
raw D_rank and near-flat: the r300 bridge correlation lived on the
*raw* positions, the max-based equidistribution shortcut is closed.
**B3:** the difference profile decorrelates at lag ONE (l_loc med
1.0) and the atom-level participation grows at +0.942 (sp +0.96 to
the block level), but the LOCAL sum-rule census breaks 24/24 — the
`|dw|` non-orthogonality is scale-free.  Sealed verdict
`NEFF_SPLIT(B1(UNIF_TARGET) | B3(MIX_TARGET))` + `REST_FROZEN`;
bycatch: the new UNIF class does *not* separate worlds (MAIN CV² 1.03
bracketed by EPST 0.72 / SCR 2.76) — unlike the O-sign and FILL
classes.  Must-fails: the Fractions identity mutants (9 and 7/24
exact), the halved discrepancy prefactor (0.3 exact), the transposed
aggregation matrix (37/55 exact + column-sum flag), overlapping
sub-ladders CAUGHT, the synthetic one-block collapse to `n_eff == 1`.
NO L2-lemma claim, NO RH CLAIM.

r302 executed the r301 frozen rest with the **unif-target round**
(`unif_target_probe.py`, 30/30, SPEC_SHA `36df9424`;
`PRIME.PORT.L2.UNIF_TARGET.01` — experiments-side; consumed by `v967`, wave 10): *why* is the block profile of `|PDelta|` quasi-uniform
(slope(1 + CV²) ≤ +0.094, measured +0.039), and can it be derived
source-purely?  **The anatomy (A1):** the mild CV² rise is
**structureless** — no sealed class scheme (sign / position band /
atom-count pattern) fires the carrier rule (between-shares
0.012/0.216/0.111).  **The stationarity (A2, mechanism iii):** the
per-rung normalized block values have ONE stable distribution along
the ladder — pooled KS(G1,G3) 0.043, KS(G2,G3) 0.021, far under the
sealed 0.125: `PROFILE_STATIONARY`.  CV² boundedness is the second
moment of a stationary normalized profile, and the r299 cconv
negative is explained: convergence in *distribution*, not pointwise.
**The depth honesty (A3):** the r301 half-ladder flag is resolved by
exact per-half additivity — the deep-half flattening of n_eff
(+0.982 → +0.802) is carried ENTIRELY by the CV² head (+0.228), not
the count (+1.030); the sealed finite-size model fires
`TRANSIENT_1_OVER_N` (m2 group medians 2.07/2.02/2.00 fall on an
exact 1/N law onto A = 1.973, held-out dev 0.002): `DEPTH_HALF_MISS`
stands as an honest flag, `DEPTH_CAVEAT` does **not** fire; the
sl_cv2p jackknife is 0/42 above UNIF_NEED.  **The structural finding
(B1):** the exact coherence identity `1 + CV² =
n_act·χ/(surv²·n_eff_atom)` with `χ = D/Q` the in-block coherence
factor (module-own Fractions re-proof, wards ≤ 6.0e-16 on 47 worlds,
per-block Cauchy–Schwarz cap χ ≤ kmax with slack 0.0) — χ med 0.63
is **destructive and non-growing** (−0.060; the r301 lag-1
anti-correlation ρ₁ = −0.22 is its atom-level mechanism), surv flat
(−0.020), composed slope +0.039 == sl_cv2p exactly (5.4e-16).
**Closed honestly:** B2 local-pattern (PAT12 within-share 0.685 >
0.5 — the node pattern does not determine `|PDelta|`; k-profile KS
0.181 > 0.125) and the B3 damping clause (localized-perturbation
gain g/N med 1.079 vs the sealed 1.0 — the recursion responds *at*
the polynomial-degree rate, miss 0.079).  Sealed verdict
**`UNIF_DERIVED(B1 + A2)`** — the first DERIVED verdict of the lane,
with every hypothesis typed MEASURED; Leg C printed the full
r297→r302 candidate theorem (wave-10 candidate, NOT promoted).  The
one remaining measured growth statement is **`ATOM_TARGET`**:
slope(n_eff_atom) ≥ sl_nact + sl_χ − 2·sl_surv − UNIF_NEED = +0.888
(measured +0.942, margin 0.055 — the same thin margin relocated one
level down, onto an atom-level anti-concentration statement about
the c-difference profile itself).  Must-fails: wrong class weight
(3/4 exact), lag-0 double counting (5/7 exact), unnormalized-KS
mutant (re-anchored to the relative form per amendment a1, original
0.335 < 0.5 printed permanently; CAUGHT 8×), the two-point synthetic
(CV² = 729/196, LOUD on every route), scope mutants flagged.
World-blind by algebra on 47 worlds; UNIF non-separation reproduced
(MAIN 1.03 bracketed by EPST 0.72 / SCR 2.76).  NO L2-lemma claim,
NO RH CLAIM.

r303 ran the **atom-target round with hard regress audit**
(`atom_target_probe.py`, 26/26, SPEC_SHA `375e9f2b`;
`PRIME.PORT.L2.ATOM_TARGET.01` — experiments-side; consumed by `v967`, wave 10), adjudicating the warning signal that r300/r301/r302 froze
three targets with the *identical* margin 0.055.  **The regress
audit (Leg A): `REGRESS_CONFIRMED`** — the four target margins
m_D/m_NEFF/m_UNIF/m_ATOM all equal **+0.0547** with invariance devs
≤ 9e-16 (module-own `margin_chain`, re-proved exact in Fractions:
all four margins 3/50 on the rational slope set; the frozen
halves-slope estimator is exactly log-additive over the
r300/r301/r302 product identities): the cascade is an exact
**dictionary** around ONE measured core `S = σ* − sl_D = +0.0547`
(the r299 DIAG margin), not three reductions — further "reduction"
rounds may **not** be counted as progress on the inequality; the
r297 σ level is *not* the core (σ* − σ = +0.1976 = S + the ratio
surplus +0.1429 exactly; only the surplus's boundedness is
structural), and the charter's ½-conversion conjecture is
**refuted** (+0.0988 ≠ +0.0547 — the factor 2 lives only in
`NEED = 2·sl_L1 − σ*` and cancels in every margin difference).
**The mixing mechanism test (Leg B, the first non-tautological
test of the lane):** sealed synthetic rearrangements of the dc
profile per rung (bitwise-exact marginal, ρ₁ steered to matched /
zero / flipped targets by seeded greedy swaps; 1008 builds,
convergence 1008/1008): sealed verdict **`MIXING_INSUFFICIENT`**
with a *monotone* mechanism ladder — χ real 0.630 → matched 0.764 →
zero 1.029 → flipped 1.342, end margin +0.055 → +0.057 → +0.032 →
**−0.044**: flipping the mixing sign *kills* the target inequality
and matching ρ₁ reproduces the slopes and the end margin, but not
the destructive-coherence level (χ miss 0.134 > 0.05) — the
within-block structure beyond lag 1 carries the rest.  Bycatch
theorem-grade fact: `n_eff_atom = L1a²/Q` is a pure **marginal**
functional (invariance 1.0e-15 on all 1008 builds) — ATOM_TARGET's
growth lives in the value distribution alone; stationarity (r302
A2) and mixing are complementary, not redundant.  The ρ-sign census
is honest: S1 < 0 on **41/42** only (one rung positive — the
adjacent-π-orthogonality connection cannot be a raw per-rung
theorem; Fractions-exact negative certificates on the two smallest
rungs).  Consequence map (ii): the reduction cascade r297→r302 is
**closed as a dictionary**; the honest next object is the joint
short-range law of the dc profile (ρ₁..ρ_k / within-block
coherence), not another coordinate change.  Must-fails: the
wrong-½-factor propagation breaks by exactly sl_L1 = 1/5
(Fractions), the unmatched-marginal family LOUD (5.9e-1), the seed
collision CAUGHT, the swapped-D/Q coherence mutant breaks by 48/35
exact, scope mutants flagged; NO amendment after freeze.  NO
L2-lemma claim, NO RH CLAIM.

r304 ran the **short-range-law round** (`shortrange_law_probe.py`,
37/37, SPEC_SHA `2cc5d23f`; `PRIME.PORT.L2.SHORTRANGE_LAW.01` —
experiments-side; consumed by `v967`, wave 10), executing the reviewer
contract verbatim: *is the anti-correlation genuinely short-range
and summable (control of 1 + 2Σρ_k), or can higher lags eat the
lag-1 advantage?*  **The law (Leg A): `LONGRANGE_STRUCTURE`** — the
centered lag profile ρ₁..₁₆ over all 42 rungs is a **stable
period-4 comb** (med −0.222/−0.140/+0.089/+0.130, then recurring
strong lags at every k = 4m (≈ +0.13..+0.16) and k = 4m+2
(≈ −0.14..−0.15) up to 16; halves-stable 3/3 in sign;
world-specific — EPST/SCR differ), *not* a decaying tail: the
sealed rule (|med ρ_j| ≤ 0.05 for all j > k₀ ≤ 8) finds **no k₀**.
The reviewer condition **splits**: the net covariance
NC(16) = 1 + 2Σ med ρ_k = **0.712 < 1 holds** (net-negative,
the sign half of the contract) while summability-with-small-tail
**fails** (SUM(16) = 1.563, no decay); the zero-sum tautology
(1 + 2Σ_all ρ_k ≡ 0 for mean-removed finite profiles) is disclosed
and gated live — the truncated NC is the content.  **The structural
counterpoint (exact):** the χ lag decomposition χ = D/Q =
1 + 2Σ_k T_k/Q over within-block position pairs (recomposition ward
4.7e-16 on every non-degenerate world) shows the *χ-relevant*
structure **is** short-range — the real shares die out by k ~ 4
(c₁ −0.345, c₂ −0.139, c₃ +0.068, c₄ +0.028, Σ(5..16) −0.026,
tail 0.000) and the r303 χ-level miss 0.134 is attributed to
k ≤ 3 (dominated by lag 2: +0.181): the long-range comb is
**invisible to χ**.  **The sufficiency test (Leg B):** lag-8
matching (252 builds, 252/252 converged at tol 0.02, marginal
bitwise) reproduces the destructive χ **level** (0.652 vs 0.630,
|d| 0.022 ≤ 0.05 — where r303's ρ₁-only family missed by 0.134)
but breaks the **slopes** (sl_cv2p/sl_D miss by 0.028/0.027 >
0.02): the exact *complement* of r303 — the mechanism is
**two-scale** (level from the short lags, slopes from the
rung-wise lag-1 trend); no lag-matching family reproduces both
within the sealed bands.  The graduated reviewer ladder (absolute
ρ₁ targets +0.2/0.0/−0.1/−0.2/−0.3, 840/840 converged): χ
1.361/1.019/0.889/0.765/0.632 **monotone** in ρ₁, but the end
margin (+0.036/+0.074/+0.019/+0.045/+0.059) is **not** — honest
negative: absolute-level targets do not kill the inequality; the
r303 kill (−0.044) came from per-rung *flipped* targets — the
margin responds to the rung-wise trend of ρ₁, not its level.  Sign
pattern: med signs of S₁..₄ = −/−/+/+ with Fractions-exact
certificates on both smallest rungs (kz18, kz23: 4/4).  Leg 0
anchors bit-near (margin chain +0.0547 ×4; the r303 three-family
ladder replicated with the same seeds: χ 0.764/1.029/1.342,
margins +0.057/+0.032/−0.044; RHO_SIGN 41/42).  Consequence map
(iii, the reviewer stop case): the lane's **global-profile mixing
route is closed** (documented stop) — L2 ⟺ anti-concentration of
an explicit block field with *long-range (period-4 comb)*
structure; what survives as honest structure: the exact
within-block short-range decomposition of χ, the two-scale
mechanism split, NC < 1 and the exact sign pattern.  Must-fails:
factor-1 net covariance (1/2 ≠ 0 exact) and factor-1 χ
decomposition (6/7 vs 5/7 exact), unmatched marginal LOUD
(5.9e-1), uncentered field (3/4 exact + real 0.028), the
ρ₁-only-family smuggle against the lag-k₀ ward (dev 0.332 > 0.02)
CAUGHT, scope mutants flagged; NO amendment after freeze (the
pre-record G32 print fix disclosed in the spec).  NO L2-lemma
claim, NO RH CLAIM.

**Wave 10 (2026-08-26)** froze the cascade closure and the documented
lane stop (rounds r301–r304) as **`v967_l2_cascade_closure.py`**
(ledger `PRIME.L2.CASCADE_CLOSURE.01` [E]; the four probes embedded
byte-exact in the sealed `--smoke` stage — SPEC SHAs
`6f8cc404`/`36df9424`/`375e9f2b`/`2cc5d23f` — plus a module-own exact
section S0: the r303 margin-invariance algebra (all four target
margins `m_D = m_NEFF = m_UNIF = m_ATOM = 3/50` EXACT on the rational
slope set; the halved-need mutant breaks by exactly `sl_L1 = 1/5`; the
record-decimal chain 0.055/0.908/0.094/0.888 exact; the σ
decomposition 0.198 = 0.055 + 0.143 exact; the ½-conjecture
refutation 0.099 ≠ 0.055), the r302 coherence identity
`(1 + CV²)·surv²·n_eff_atom == n_act·χ` on two rational toys (mutant
breaks 5/7 and 48/35 exact), the r304 χ-lag decomposition
`χ = 1 + 2Σ T_k/Q == D/Q` on both toys (5/7 and 5/19 exact; factor-1
mutant 1/7), the r304 zero-sum tautology as a gate (factor-1 and
uncentered mutants 1/2 and 3/4 exact; NC(1) = 0.556 exact, NC(16) =
0.712 < 1; the period-4 comb certificate with NO k₀), the −/−/+/+
sign-pattern certificate with the kz18/kz23 Fractions pins, the four
sealed verdict bars as exact decision logic with tipping mutants, and
the lane-stop composition gate).  **The wave freezes the rejections
and retypings as the content**: the cascade r297→r302 is RETYPED as
an exact reduction *dictionary* (coordinate finding around one
measured core `S = σ* − sl_D = +0.0547`, not six proof steps; the
hard reviewer rule binds — a round counts only with NEW information);
the first causal coordinate is real but insufficient (the ρ₁ flip
kills the inequality, the χ level needs the short lags: the mechanism
is TWO-SCALE); and the sealed reviewer stop case fired
(`LAW_LONGRANGE`): **the global-profile mixing route of the L2 lane
is documented CLOSED** — the honest stop state is *L2 generic ⟺
anti-concentration of an explicit block field with long-range
(period-4 comb) structure; return later with new tools*.  What stands
honestly: the exact identity dictionary (count, coherence, χ-lag),
the two-scale split, NC < 1, the −/−/+/+ sign pattern, the
marginal-functional invariance of `n_eff_atom` (1e-15 on all
2100 + 1008 builds).  `PRIME.L2.NEFF_TARGET.01` [O] is retyped to the
documented stop state; `PRIME.L2.REDUCTION_CHAIN.01` [E] keeps its
exact identities with the reading refined.  In Lean **no new
statement** is added: the stop carries no new exact target form (the
period-4 comb and the two-scale split are measured ladder aggregates,
not exact finite identities, and a universally quantified ladder form
would be refutable — the r273 guard pattern of
`RH/Counterexamples.lean`); the hole stays `lstar_subordination`.
NO RH CLAIM.

**Lane status after wave 10.**  **L2 (edge A, the generic half):
documented STOPPED** — the global-profile mixing route is closed
(`LAW_LONGRANGE`), the clean equivalence state is frozen (`L2 generic
⟺ anti-concentration of the explicit dc block field`, core slack
`S = +0.0547` measured), and per the hard rule no further coordinate
change counts as progress; return with new tools.  **L\* (the open
center): dormant since r296** — the last L\*-side round was the DENS
fork closure (`DENS_WORLD_BLIND`, consumed in wave 9); the contract
`PRIME.LSTAR.SUBORDINATION.01` [O] stands at 57/57 sign-safe measured
with the harmless margin law, and `PRIME.LSTAR.CLOSED_FUNCTIONAL.01`
[O] keeps its named open fronts (loss-corpus forensics, why-L2, the
window knife edge).  Both edges rest with a clean documented state;
the mincut (base 4 / refined 5) is unchanged.  NO RH CLAIM.

r306 ran the **Rényi-3 round** (`renyi3_probe.py`, 27/27, SPEC_SHA
`3bb365e1`; contract `PRIME.L2.ATOM.RENYI3.01` — experiments-side,
exploration only, no ledger row; reviewer plan §5, the sharpest new
fiber attack — explicitly sanctioned as *new information* after the
r304 documented stop because it is **pointwise** instead of a slope
target and uses the **cubic moment** instead of the closed
max-/discrepancy routes): prove or refute
`Σ_j q_j³ ≤ C·(log m)^A/m²` pointwise on every rung of the atomic
|PΔ| profile (q = |PΔ|/L1 over the m level-2 blocks), the family
A ∈ {0, 1, 2} preregistered and C frozen on the first 5 rungs of the
(N, kz)-sorted ladder.  **Verdict `RENYI3_GO(C = 1.069, A = 2)`** —
the inequality **holds pointwise on all 57 rungs** (42 core + the 15
r286 extension anchors N_w 942..1218, rebuilt with the full block
machinery: contribution ward 4.9e-13, count link n_act == m 15/15)
with **growing reserve** (trend −0.322 on both frozen estimators;
the whole extension holds at reserve 2.6..5.7).  The polylog
exponent is real: A = 0 fails (the third shape moment
M₃ = m²S₃ rises +0.153 before decaying — exactly the disclosed
algebraic prior from the r301 records), A = 1 fails on precisely the
two near-critical rungs kz53/kz67, which sit 1.4–1.9 % under C₂ —
polylog-**squared** is the honest sharp form, with kz53/kz67 as the
sharpness witnesses.  Via the exact Rényi/Hill chain (two Lagrange
witnesses in pure Fractions, dev 0; equality exactly on the uniform
vector; the bridge rational with all constants; Lean-ready statement
printed for the Lean worker) the GO fixes the **new pointwise proof
target of the fiber**: `n_eff = N₂ ≥ N₃ ≥ m/(1.034·log m)` —
asymptotically above every demanded power m^0.888.  The anatomy is
the friendliest possible: the cubic mass is a **narrow-band**
functional (97 % within a factor 4 of q_max) carried ever more
**broadly** (top-8 share 0.780 falling at −0.268; strict triples
92 % of the cube) — no heavy-few obstruction; the minimal sufficient
statement is named **SHAPE3_TARGET** (M₃ polylog-bounded, the
third-moment sibling of the r302-proven M₂ stationarity 1.973; M₂
med here 2.03).  World census with the right sign:
**WORLD_SENSITIVE** — EPSTEIN holds the bound, SCRAMBLE **breaks**
it by 1.67× (the law is recursion structure, destroyed by
scrambling, not block combinatorics — the first world-separating
pointwise law of the lane).  Leg 0 anchors bit-near (Rényi family
meds/slopes, core slack S = +0.0547 at dev 6.7e-16, count link
42/42, n_eff_atom 118.0/+0.942, marginal invariance ≤ 8.7e-16).
Honest negatives: C₂ was frozen on the *shallowest* rung and the
kz53/kz67 margin is 1.4–1.9 % — the GO fixes a proof **target**, it
proves nothing beyond the 57 measured rungs; the atom-level (dc)
sibling stays unformalized.  Must-fails: power-mean flip breaks by
1/64 exact, the preregistration-breach mutant CAUGHT (set ward +
real C_full/C_cal 2.10 at A = 0), unnormalized q LOUD (factor L1³),
one-atom collapse 212× LOUD, scope mutants flagged; the two
disclosed pre-record code fixes (print format; mutant reference
aligned to the spec's w9) moved no bar, band or rule.  NO L2-lemma
claim, NO RH CLAIM.

r307 ran the **fixed-head kill round** (`fixed_head_probe.py`, 28/28,
SPEC_SHA `ec2bb008`; reviewer adjudication §4 solution C — the
cheapest kill test of proof architecture C (fixed head + contractive
tail); contract `PRIME.LSTAR.FIXED_HEAD.01`, experiments-side,
exploration only, no ledger row — this reopens the L\*-side after the
r296 dormancy noted above): in the μ-orthonormal frame
`E = E_H + E_T` with H = the first r **flat archimedean edge atoms
below the first prime** (u < log 2; the r284 carriers), under two
sealed window-fixed orderings (ascending fold; the frozen w9 extremal
rank permutation, regated exactly), λ_max(E_T(r)) was measured for
r = 0..16 on **77 windows** (42 core + 15 r286 anchors + 20 *new*
deeper windows N_w 1218..1650).  **Verdict `HEAD_GROWS`** — the
kill test *splits* the reviewer's question: **no fixed r ≤ 8 exists**
(best candidate leaves min reserve 6.5e-6, four decades below the
sealed 1e-2 bar; a fixed r_min ≤ 16 exists on 15/77 windows only,
all N_w ≤ 434; the fixed-head reserve falls with S at sp −0.93,
TS slope −3.15 — the r286 margin-law class: a fixed head changes the
prefactor, not the decay class), **but** the *full* edge (13..96
atoms, growing with the window at sp(edge, N_w) = +0.98) restores a
macroscopic reserve on 77/77 (min 2.19e-3, median amplification
3.1e4): the near-one-atom anatomy is **real and proof-relevant as a
coordinate, refuted as a fixed dimension** — architecture C in its
fixed-head form is dead; the honest successor object is the
window-dependent head of size ~ log2/(2Δ_w), whose own reserve
decays slowly (TS −0.87, disclosed).  The Schur/Woodbury ledger is
exact: `det(I−E) = det(I−E_T)·det(S_H)` on 57/57 (rational-exact on
the JF9 toy), the border-augmented identity against the independent
r244/r258 route on 57/57 (q_N < 1 reproduced, min gap 0.0195);
λ_min(S_H(8)) > 0 on 57/57 **but equals the wall margin to 4 digits
and is anchor-monotone along the r286 margin decay (sp −0.92) — the
hoped-for simplification is empty: S_H carries the whole problem**;
the border is inert (S_aug positive, degradation ≥ 0.873 —
`BORDER_BREAKS_HEAD` did not fire).  Hardening: the r289 rational
twin is `TWIN_METRIC_OK` (identical r_min, max |Δreserve| 1e-7);
EPSTEIN/SCRAMBLE are `CONTROL_COLLECTIVE` (EPST 2191 → 2189 at
r = 16, SCR 1.4e7 → 3.2e4 — no small head rescues a dead world);
Gershgorin/Schur-test row sums stay 1.5..3.6 on the tail (classical
bounds dead, as r283 measured on the full E); mp spots 3/3 at dps 30
(worst 2.8e-14); must-fails m1 target-inverse ordering AST-flagged /
m2 wrong-sign Woodbury 0.163 loud / m3 monomial frame hides the toy
wall 0.805 loud / m4 r = 0 exact.  Four pre-spec sizing scratches
disclosed in the spec.  NO L\* claim, NO RH CLAIM.

r308 ran the **block-Green discovery round** (`block_green_probe.py`,
30/30, SPEC_SHA `d5147850`; reviewer adjudication §4 solution A /
plan round "305" — experiments-side, exploration only, no ledger
row): the automated search for a **matrix-valued discrete Green
identity on fourfold fold blocks** — the route the reviewer named
after the scalar profile functionals (r290–r295), the
absolute-value majorants (r297/r299/r300) and the diophantine
relations (r289 `METRIC_ONLY`) all failed: does the augmented form
`Q(p,t) = ∫p² d(μ−ν) + 2t·u(p) + B·t²` (the v958 bordered object,
`B = S_{N−2} + 5/7`) admit an exact decomposition
`Q = Σ_r ⟨G_r Δ_r, Δ_r⟩ + (5/7)t²` with sealed fourfold-block
coordinates (first/second differences, the r288 antiphase
`layer_i − layer_{i+2}` coordinates, the local gross-mass mean, the
border `t`) and source-pure blocks `G_r`?  **Verdict
`IDENTITY_EXISTS` + `BLOCK_INDEFINITE_MAIN` + `FEAS_DIAG` world
separation.**  The identity **exists everywhere** (exact Fractions
with entrywise reconstruction on the hand toy, four synthetic small
models and the frozen real MINI16; f64 residuals ≤ 8.1e-14 on w9,
the r289 rational twin, EPSTEIN/SCRAMBLE/SMOOTH and all 57 rungs at
deg ≤ 8; kernel exclusion 100 %) — existence is the *disclosed easy
half* (linear in G, hugely underdetermined, dof 15..7593).  The
sealed eigenvalue-free minimum-Frobenius selection is indefinite
**everywhere** — even on the positive-class calibrator SM0: the
*selection*, not the world, produces the indefiniteness — so
`BLOCK_INDEFINITE_MAIN` fires on the letter of the sealed rule.
**The round's central measurement** (sealed diagnosis clause,
calibration amendment a1: one uniform two-stage Dykstra schedule
200/2000 for every world): the SDP-like feasibility census
**separates the worlds** — w9 and the rational twin at deg ≤ 8
converge to *genuinely all-psd block families* (min eig rel
+6.6e-16 / +2.0e-17 at affine residual 3.7e-14, after 200 steps)
while EPSTEIN/SCRAMBLE stall at −0.45/−0.49 after 2000 steps: on
the restricted subspace, where **every** world's target is positive
definite and nothing is theorem-forced, MAIN and its
diophantine-trivialized twin admit an all-psd fourfold-block Green
decomposition and the two hard controls (at this schedule) do
not — the first block-level world discriminator of the L\* lane,
**one-sided** (non-convergence proves no infeasibility) and
diagnosis-grade only.  SMOOTH sits marginal (−1.3e-05; its crossing
28 equals the second cap DEG_B — the sealed boundary case);
w9/twin at DEG_B = 28 stay **open** (−4.2e-05 after 2000).  The
DEG_B control negativity (crossings 26/22/28) is *forced by
theorem* given existence — gated live on EPSTEIN/SCRAMBLE (held),
disclosed, never sold as discovery.  R282 demarcation held
explicit: the Kasteleyn/SOS refutations concern the *full* signed
configuration class (every Cauchy–Binet cell, `h_n` itself); this
identity lives on the restricted subspace deg ≤ 8/28 ≪ N_w plus
one border coordinate — a different language, no contradiction.
Anchors bit-near (w9 367/263/104, crossing 185; r288 carrier pin
z_v −3.149 destructive; r282 pin N₂(MAINLIKE) exact; B_w9 =
8.368649 with 182 nonnegative ρ_k).  Must-fails: cholesky-of-target
FLAGGED, omitted remainder breaks by exactly 5/7 at (t,t), posthoc
family FLAGGED, partial monomial coverage leaves an exact defect
4.5e+01 LOUD.  Honest: a feasibility census fixes a **target**
(find a sealed *constructive* psd G rule on MAIN, extend it past
DEG_B, prove the world split), it proves neither L\* nor any
cofinal statement.  NO RH CLAIM.

r309 ran the **paired-cone pilot round** (`paired_cone_probe.py`,
33/33, SPEC_SHA `f8d99877`; reviewer adjudication §4 solution B /
plan round 306 — experiments-side, exploration only, no ledger
row): the third full proof architecture of the L\* lane — the
positive and negative source processed **pairwise as rank-one
updates** of the same augmented form `A = [[H, u], [uᵀ, B]]` as
r308, base and border together.  The form is rebuilt *exactly* as
an ordered signed rank-one stream (sealed base = the imported 5/7
`t`-floor + reserved/unmatched μ atoms; one signed budget step
`|S_{N−2}|`; every ν atom greedy-matched to its nearest free μ
partner; every border atom split exactly into a ±(η ± e_t) pair),
and at every negative step the reviewer's **local reserve**
`r = 1 − L_bb + L_ab²/(1 + L_aa)` is measured — the
Sherman–Morrison cancellation as an *explicit positive square*,
with the determinant dictionary `det(M₁) = (1 + L_aa)·det(M)`,
`det(M₂) = r·det(M₁)` exact on the Fractions toys (hand toy PC4:
reserve 7/9, dets 9/14 → 27/14 → 3/2) and f64-warded on w9
(reconstruction 4.8e-16; mp ward dps 60, dev 3.0e-13).  **Verdict
`DECOMP_EXACT` + `MAIN_RESERVE_POSITIVE` + `CONE_PARTIAL` +
`WORLD_DISCRIMINATOR`(all sealed stats world-blind).**  MAIN
survives the ordered chain at **both** caps with macroscopic
margin (w9 min r +0.9011 at deg ≤ 8 / +0.3232 at deg ≤ 28; the
r289 rational twin bit-near; 57/57 ladder rungs all-positive, min
+0.7357@kz97), and w9's two tightest pair reserves sit *exactly*
on the r284 extremal folds {2, 4}.  The reviewer's decisive cone
question **splits**: the sealed mass-ratio cone
`r ≥ 1 − c₁·(v/w_a)` (c₁ = 8.5436, frozen on the five shallowest
rungs) is **sound** — zero violations on all 20714 pointwise test
steps (52 rungs + w9 + twin) — but certifies only 56 % of the
steps (positive predicted bounds on 11643/20714); the sharper
local shield cone `r ≥ 1 − c₂·(v/SM₂)` overshoots exactly where
it is positive (54/54 positive bounds violated, every violated
step has true r ≥ +0.901).  Reading: the reserve **is**
predictable in the sound direction from one source quantity, but
no sealed *local* φ is yet a certificate — the named successor
object is a φ with a base-depth term (the deg-8 base carries most
of the reserve), not retro-fitted.  The information test stays
silent by the sealed letter (pooled tightest decile best |sp| =
0.50 at pair distance; on w9 alone the mass-ratio/shield features
reach −0.61/−0.58, reported), and the reviewer-B4 O(log N)
fold-bit register is measured **irrelevant** at this cap (max
|sp| 0.32).  Honest world census: at deg ≤ 8 the ordered *pair*
reserve is world-blind (the controls' pair phases are all
positive); the *entire* EPSTEIN/SCRAMBLE indefiniteness there is
the negative **budget step** (S_{N−2} = −3.992/−5.237 past the
flip — the r243 budget coordinate reads the wall; amendment a1
disclosed: the design-time "every deg-8 target is PD" sentence
was corrected, spec-prose + reporting only, no bar or rule
moved); the deg-28 control breaks are theorem-forced and
disclosed.  Anti-circularity at r275 sharpness: the o1/o2 +
forced-tail-sum pin re-run exact; the target-inverse cone mutant
and the posthoc-pairing mutant FLAGGED; wrong-sign
Sherman–Morrison caught by the exact det dictionary (dev 8/9
exact); the order-break mutant shifts the reserve by exactly 4/9
while the reconstruction stays exact — the r-sequence, not the
sum, is the ordered object.  Honest: a passing cone fixes a
**target** for proof work; this round proves neither L\* nor any
cofinal statement.  NO RH CLAIM.

r311 ran the **block-Green nontriviality round**
(`blockgreen_nontriviality_probe.py`, 34/34, SPEC_SHA `fac7a8df`;
the reviewer's *mandatory* round before any constructive G-search
— experiments-side, exploration only, no ledger row): is the r308
psd block family stronger and more local than `Q_w ⪰ 0`, or a
sparse-Gram / chordal restatement whose SDP separates the worlds
by force?  **Verdict `STRICT_SOURCE_CONE` with the binding
separation-mechanism clause** — both reviewer scenarios are true
at once, on different layers.  (1) **The r308 world separation is
fully the budget sign**: the sparsity graph of the Δ coordinates
(width-3 atom band + universal border vertex `t`) is **chordal**
on every size including the full w9 S = 367 graph (MCS +
Tarjan–Yannakakis, maximal cliques == the 364 sealed blocks;
C4/C5 control correctly non-chordal), every block spans its full
local R⁵, so on *evaluation* space the image cone is exactly the
chordal pattern-psd cone (Agler et al. / Grone et al.); the
control targets are **indefinite** at deg ≤ 8 through the one
negative budget diagonal (EPST −3.992 / SCR −5.237; λ_min rel
−0.86/−1.00 vs w9/twin +8.4e-4 PD) — the r308 Dykstra separation
was theorem-forced; the trivial common dual `Y_t = e_t e_tᵀ`
(**one** rational entry — the r243 budget positivity
`S_{N−2} ≥ 0`, the *known* wall reading) closes the r308
one-sided stalls two-sidedly and predicts two **new blind**
scramble worlds 2/2 before their runs; the decisive ablations
seal it: EPST/SCR with `(t,t) ← |S|` become PD and the Dykstra
**converges**, w9 with the EPST budget transplanted **dies** —
flip one scalar and the separation inverts; the exact Schur-tail
lemma (`Schur(A_sys@d) = S_{N−2} − Σ_{k≤d} F_k²/h_k`, exact on
SM0/SM1) shows the wall sign enters the target mechanically as
the budget ρ-tail and decides all four r308 small-model OPEN
stalls as trivially target-indefinite (even the positive-class
SM0: `−(ρ_3..ρ_8) = −7.08e-2` exact); SMOOTH's marginality is a
*boundary* target (S = +3.48, λ_min ≈ 0).  (2) **But the cone is
genuinely strict on the compressed subspace** — the round's
second discovery: 0/16 sealed random in-span PD samples decompose
(w9A 0/6, MINI16 0/6, SM1 0/4; staged Dykstra up to 20000), every
MINI16/SM1 stall carries a valid polished-numeric Farkas dual,
and MINI#0 + SM1#0 carry **fully exact in-span rational
certificates** (exact Chebyshev basis, denominator 1e6,
`⟨Y,Q⟩ < 0` exact; MINI16 = the real-w9 miniature with exact
full span 55/55) — `C_w` is *strictly* smaller than
span ∩ S₊: block membership is a genuinely stronger property
than positivity, which MAIN + twin **have** at deg ≤ 8 and
generic positive forms **lack**.  The ambient span layer is exact
too (codims TOY4 0 / SM0 11 / SM1 11 / SM2 2 / SM3 1 / MINI16 0
/ w9A 4 / w9B 236; exact annihilators, SM1 rank-one ambient
counter-model with exact rational Farkas).  Protocol audits
clean (r308 `feas_diag` objective-free, start ablations);
w9@deg 28 stays honestly **open** (PD-thin, −1.78e-5 after
20000, dual diverges — one-sided).  Consequence: the planned
R312 rank-one-conductivity round premised on the r308 *world
discriminator* is dead in that form (the discriminator was the
wall sign); the re-scoped GO object is cone **membership** — why
the bordered-Hankel world targets sit inside the strict cone
while generic PD forms do not.  NO RH CLAIM.

r312 ran the **block-Green membership round**
(`blockgreen_membership_probe.py`, 32/32, SPEC_SHA `6c32f749`;
THE one construction round granted by the reviewer tree after
r311's `STRICT_SOURCE_CONE` — experiments-side, exploration only,
no ledger row): does a source-pure **constructive rule** (no SDP
optimization in the final statement) write down an explicit psd
block family `G_r = Σ_ℓ c_{r,ℓ} v_ℓ v_ℓᵀ` with `c ≥ 0` for the
world targets, over a *pre-sealed* 22-generator library in the Δ
basis (fold distances 1/2, antiphase 3/4, the r307 arch
coordinates, border-vs-interior pairs `D_a ± D6`, symmetric /
antisymmetric four-block modes), with an explicit coefficient
formula `c_{r,ℓ}`(source quantities)?  **Verdict
`COEFFICIENT_SIGN_WALL`** — the library **span** carries the w9
identity exactly (unconstrained rel 2.4e-13; abstract span 15/21
exact) but the **cone** does not: deterministic Lawson–Hanson
NNLS stalls at rel 4.7e-4 with a **verified I-polished Farkas
certificate** (δ = 1e-4, min normalized column product +9.7e-5,
`⟨y,q⟩ = −0.9975`), and the certificate mass localizes on the
**border/budget row** (t-row fraction 0.583 ≥ the sealed 0.5; at
deg 28 even 0.725) — the identity forces negative coefficients
exactly on the wall coordinates: *the sign question of the c IS
the wall, in the constructive language too*; the twin is
type-identical (METRIC_ONLY holds).  **Exact grade** on the
miniatures: SM1 and MINI16 carry **fully exact I-polished
rational Farkas certificates** (den 100, `⟨Y,q⟩ = −0.9988 /
−1.0000` exact rational < 0, all 154/286 exact column products
≥ 0), and on SM1 even the library *span* fails exactly
(`[M_lib | q]` rank 41, EXISTS FALSE); rung census 1/57
NNLS-feasible (worst rel 8e-4 at kz12) — the sign wall is uniform
across the ladder; the four sealed formula candidates
Φ_A..Φ_D (masses, tent fractions, budget share, per-type table;
calibrate-on-w9 / freeze / verify-on-57+twin) never fire (best
calibration rel 7.7e-3 ≫ 1e-8).  **Membership anatomy** (the
round's discovery): the converged r308 Dykstra family is
high-rank ({rank 3: 17, rank 4: 79, rank 5: 268}) and
library-*aligned* (top-eigvec alignment med 0.9973) with
library-cone share med **0.976** — 97.6 % of the psd family lives
in the sealed rank-one cone and the *whole* membership
obstruction sits in the remaining 2.4 % of SDD-forbidden negative
cross mixtures; the dual anatomy of the 10 r311 sample duals
shows the generic-missing directions are low/mid-degree
**interior** mass (|z_t| med 0.000), not border mass.  w9@deg 28
stays honestly **open**, two-sidedly hardened (staged anchor
−1.784e-5 == the r311 record; library NNLS@B infeasible, same
wall class; +80000 extension steps improve only to −7.66e-6, a
slow tail consistent with asymptotic feasibility, not decided;
the strengthened dual POCS runs but polishes invalid).  Leg 0
reproduced every r311 record bit-near on the identical rng stream
(0/16 strictness, 10/10 polished duals, both exact certificates).
Consequence (sealed reviewer tree): **not a GO** — lane A closes
as cone language with the membership mechanism *named*
(block-psd membership **without** rank-one-SDD membership,
obstructed at the budget/border sign); the resources move to the
fiber.  Two disclosed calibration amendments (f64 I-polish
certificate bookkeeping; the m3 harness probes the r289
rationalization variable u/Δ) — the four-way adjudication tree
never moved.  NO RH CLAIM.

r313 ran the **Rényi-3 proof-form fork round**
(`renyi3_proof_fork_probe.py`, 32/32, SPEC_SHA `6505dd10`;
reviewer plan: triple incidence vs four-step Floquet —
experiments-side, exploration only, no ledger row): two sealed
proof architectures for the r306 pointwise cubic law
`Σq³ ≤ C·(log m)²/m²`, executed side by side under the binding
reviewer instructions (do **not** chase the sharp 1.069; the
theorem must not depend on von-Mangoldt arithmetic — sought class
`RecursiveDifferenceProfile` with MAIN + EPSTEIN in, SCRAMBLE
out).  **Verdict `BOTH_PARTIAL` (named remnants).**  Route R3A
(triple incidence): the raw atomic presentation is **exact** —
every block value is the sum of local c-differences, atoms at one
position form a β/ω **fold pair**, and the cube splits exactly
(block ward 3.9e-16 on 61 live worlds) into T1 diagonal / T2
two-equal / T3 near (shared fold) / T4 far (fully separated);
med signed shares +0.38 / +2.15 / −1.14 / −0.42 — the types
**cancel against each other** (T2 alone carries 2.15× the cube),
they are not independently small.  Two proof hopes are banked:
the fold multiplicity is **uniformly 2** on 57/57 rungs (an exact
structural class bound, not a census) and the far triples
genuinely **telescope** (TC_far med 0.069 — 93 % of their abs
mass cancels, falling).  Per-type prereg at A = 2
(first-5-frozen constants): T2 and T3 HOLD 0/57 with falling
trends (−0.400/−0.299); T1 fails 2/57 (kz55 2.54×, kz53 1.70×)
and T4 fails 1/57 (kz55 1.84×), both with **falling** trends —
shallow-calibration artifacts of the r306 A ≤ 1 kind, but the
sealed rule stands: COMPOSITION REFUTED-IN-FORM.  World census:
SCRAMBLE violates exactly T2 (2.85×) but EPSTEIN also violates
T2/T3 while its *total* holds the r306 bound — MAIN-frozen
per-type constants are **not world-portable**.  Route R3B
(four-step Floquet): the local transfer step is **real and
exact** (`r_{k+1} = t_k + ap_k·r_k + bp_k·r_{k−1}`, ward 1.9e-12,
det dictionary 4.4e-16, SMOOTH degenerate-skipped), but the
four-step monodromy on the split cubic mode space **expands**
(med sr 1.0600, contracting rungs 0/57, product exponent +1.917
vs the required −2): after r304 (single step) and this round
(period-4 Floquet step), **no transfer-operator contraction form
of the cubic law remains open at the chain level**.  The class
finding (sharpest honest negative): the four-block comb property
P4 excludes SCRAMBLE but **also EPSTEIN**; every other sealed
candidate (fold share, recursion, multiplicity, mass, assignment,
telescope) is world-**blind** — no boolean property separates
worlds the way the r306 total bound does; the separation lives in
the *size* of the type terms, so the class needs a
**quantitative membership functional**, not a property list.  The
m4 shuffle instance is rejected by P6 (293/300 atoms outside
their interval) while mass conservation still holds on it.
Anchors bit-near r306/r304 throughout.  Honest: nothing here
proves or refutes the cofinal law.  NO RH CLAIM.

r314 ran the **signed-cubic-flux round, part 1**
(`signed_cubic_flux_probe.py`, 29/29, SPEC_SHA `841b3196`;
reviewer contract after the r313 adjudication — **exact algebra
only**, no size estimate, quantification deferred to R315/R316;
experiments-side, exploration only, no ledger row): with
`x_j = (PΔ)_j`, `σ_j = sign(x_j)`, expand `|x_j|³ = σ_j·x_j³`
over the r313 fold genealogy into the signed cubic tensor
`C_{αβγ} = Σ_j σ_j·x_{j,α}x_{j,β}x_{j,γ}` and decompose
`C_cube = ΔF + C_collision + C_boundary`.  **Verdict
`SIGNED_CUBE_IDENTITY` with a vanishing boundary** — the
decomposition exists exactly: on w9 it is decided in **exact
Fractions on the real data** (35 blocks / 150 fold groups; brute
tensor enumeration == Newton aggregates == flux telescope, every
deviation the rational number 0 — `Fraction(float)` is exact on
dyadic f64 atoms) and to 4.5e-17 in f64 on all 57 rungs + mains
+ EPST/SCR (bar 1e-13).  The divergence form is genuinely
**local**: far triples enter as telescoping edge fluxes
`dF = 3G(s1² − s2)` consuming only the transported prefix state
along the position-ordered fold groups (telescope Newton-vs-path
worst 4.1e-16), and the **opening-flux lemma**
`F₁ = G³ − 3G·G² + 2G³ == 0` is cubic algebra, so
`C_boundary == 0` — the sealed language class closes without
remainder (disclosed: the razor acts *before* the genealogy; an
unmasked presentation would re-introduce a boundary term).  The
collision space is exactly **countable** by the banked r313
multiplicity-2 asset: config counts full/pair/far ==
`n / 3n(n−1) / n(n−1)(n−2)` and atom-triple count
`3p₁p₂ − 2p₃ == 8(n + 3n(n−1))` exact on 61/61 live worlds.
Measured preview for R315/R316 (labeled, no bound claim): med
signed shares ΔF −0.4226 (falling −0.422), C_pair +0.5980,
C_full +0.8537 (rising +0.213) — the two collision terms
together (+1.45) carry the cube against the negative far flux;
flux cancellation FC med 0.629 falling −0.141; the world table
(w9/w13/EPSTEIN/SCRAMBLE) is sealed record data — the SCRAMBLE
localization census names the **full-collision** column, but
EPSTEIN sits far from MAIN in every share while holding the
r306 total bound: single shares do **not** separate worlds; the
R315 functional Φ₃ must be built from the divergence-form
*structure* and must not be calibrated on these tables (sealing
discipline, disclosed).  Anchors bit-near r313 (type shares
+0.3808/+2.1542/−1.1442/−0.4226, TC_far 0.069/−0.050, fold
multiplicity == 2 on 57/57).  Must-fails: σ dropped (toy break
2 exact, w9 1.8e0 LOUD), double-counted fold group (breaks 2
and 98 exact, CAUGHT), unordered support (per-edge break
(0, 108, 108) exact, the total disclosed permutation-blind,
LOUD), mult-blind collision count (break 56 exact, 8× wrong,
CAUGHT), scope mutants flagged.  Honest: part 1 is bookkeeping
only — nothing here bounds anything.  NO RH CLAIM.

r315 ran the **Φ₃-functional round, part 2**
(`phi3_functional_probe.py`, 29/29, SPEC_SHA `92d35a3a`; reviewer
contract part 2 — define the quantitative membership functional
Φ₃ *in advance* from the r314 divergence-form **structure**,
normalized to the Rényi-3 scale `NORM = m²/((log m)²·L³)`, and
run the blind world test; experiments-side, exploration only, no
ledger row).  Three candidates sealed before evaluation:
`Φ₃a = NORM·(|COLL| + |BND| + |DFLUX|)` (reviewer raw),
`Φ₃b = NORM·|COLL|` (the combined signed collision sum),
`Φ₃c = NORM·|COLL|·FCIX` (flux-corrected); the honesty conditions
are **machine-audited** (AST identifier scan: no Σq³ read-back;
sealed literal scan: no r314 world-table number in any builder —
the e1/e2 must-fail mutants prove both audits bite), and the
literal `|Σ|x|³ − ΔF|` form is adjudicated read-back-adjacent and
demoted to a diagnostic column (== Φ₃b to 1.6e-15).  **Verdict
`PHI3_ALL_BLIND` on the sealed letter** — no candidate passes all
four reviewer criteria, and the anatomy is sharp and
concentrated: the *entire* failure sits on **kz55 + kz67** (the
only violators of all three first-5-frozen C₀ = 2.63/1.51/0.94;
kz55 alone tops SCRAMBLE in the rank test, Φ₃a 4.61 > 2.39, and
blows the b/c bands 66/88 vs bar 30) — the same near-critical
family that killed r306 `A ≤ 1` and the r313 T1/T4 constants,
with every trend **falling** (−0.40..−0.77 on both estimators):
the known shallow-calibration artifact, but the mid-ladder rule
was not sealed and is not applied post hoc.  The **local cause is
named by the probe's own coordinates**: kz55/kz67 carry FCIX
0.955/0.915 against med 0.629 — exactly the rungs where the
intra-block flux cancellation *dies*; the FCIX → 1 stratum is the
R316-relevant obstruction, and FCIX is source-pure and computable
in advance.  What *holds* of the reviewer table: EPSTEIN holds
all three frozen bounds (Φ₃c separates it 7× downward via its
world FC 0.101), the rational twin w13 is **near-identical** to
w9 (factor 1.07–1.17 — the strongest twin result of the lane),
the SCRAMBLE cause is named **COLL** (component devs dflux 1.33 /
coll 3.69 / fcix 0.10, K4 fires), and candidate (a) passes its
band bar (22.85 ≤ 30).  Both reviewer-table contrasts fired: the
multiplicity-3 control (half-atom split, G1/blocks/mass preserved
*bitwise*) breaks the collision count bound 48762 vs 14448 (toy
216 vs 64, break 152 exact) and is rejected by the cap while Φ₃
is provably **blind** to it — the multiplicity cap is a necessary
class side condition, not a consequence of the functional; the
assignment shuffle (seed 315001) is rejected 289/300 with the
flux profile broken edgewise.  Must-fails: e1 Σq³ read-back
AST-CAUGHT, e2 r314-table calibration literal-CAUGHT, e3
mult-2-built control accepted as normal (sharpness), e4 wrong
normalization exponent LOUD (factor m = 35 exact), scope mutants
flagged.  Anchors bit-near r314/r306.  Honest: ALL_BLIND stands —
the divergence-form structure with first-5-frozen constants is
*not* the membership functional; the kz55/kz67 + FCIX-stratum
diagnosis is post-hoc census, no rule; nothing here bounds
anything.  NO RH CLAIM.

r316 ran the **two-regime-bound round, part 3**
(`two_regime_bound_probe.py`, 35/35, SPEC_SHA `5c28b12b`; reviewer
contract part 3 — formulate the pointwise Rényi-3 bound as a
coupled two-regime statement over the source-pure FCIX stratum,
with θ = 0.85 sealed from the r315 record gap and a **mid-ladder
calibration sealed in advance** (CAL_START = n//3, N_CAL = 5,
small-m ranks certified individually); experiments-side,
exploration only, no ledger row).  Four sealed regime majorants,
all ≥ ρ₂ by exact algebra warded live (worst slack 6.5e-16 on 69
worlds): `ΦL2 = NORM(|COLL|+|BND|+FE)` (flux-index form), `ΦL1 =
NORM(|COLL|+|BND|+|DFLUX|)` (r315 raw), `ΦH2 = NORM(CNT3+|BND|+
FE)` (the multiplicity-2 collision **counting** majorant), `ΦH1 =
(m·QMAX/log m)²` (concentration); purity machine-audited (AST +
record-literal scans; θ-post-hoc, cubic read-back and literal
calibration mutants all caught).  **Verdict `TWO_REGIME_DEAD`
(L fails; H fails)** on the sealed letter — and the anatomy is
the round's find: (1) **the FCIX stratum is not the near-critical
family** — the mid-ladder freeze is killed on the regime-L side
by kz53 (ρ₂ 1.049 with bulk-normal FCIX 0.654) and deep spikes
kz83/kz105: the obstruction family cuts *across* the
flux-cancellation stratum, so FCIX is the wrong (or an
incomplete) regime coordinate; (2) ρ₂ itself is **not below its
mid-ladder window** (kz53 1.049 / kz83 0.779 vs the cal-window
max 0.458): any tight pointwise majorant fails at kz53 under a
mid-ladder freeze — the r306 first-5 protocol was not a
shallow-calibration artifact but **load-bearing** (the profile
is non-monotone with recurring near-critical spikes); (3) the
**regime-H mechanics object** delivered as census: the FCIX
outliers kz55/kz67 are near-**one-block** worlds — top-1 cube
share 0.558/0.785 vs 0.18 regime-L med (dev 2.68, the named
discriminator; kz55 additionally 10× compensation κ 0.105);
concentration, not counting tightness (dev 0.53), distinguishes
them — but kz67 misses the mid-window H constants by 1.4×/2.6×
and kz55 fell into the small-m set (rank 20/65), so the H test
stratum was one rung.  What holds: anchors bit-near (r314
shares/FC/mult, r315 C₀ + FCIX outliers 0.955/0.915, r306 C₂
1.069 with 0/57); the EXT2 extension clean (8 deep anchors to
N_w 1393, all mult 2, **no new H member** — the FCIX stratum did
not grow with depth); the 21-rung small-m certificate table
(C_small 1.0694 = the r306 shallow maximum); the class machinery
intact (EPSTEIN/twin/main admitted with the strongest twin band
of the lane 1.04, SCRAMBLE rejected via COLL attribution 3.69 +
edgewise flux break 290/300, mult-3 control rejected via the cap
with all four majorants provably blind).  Honest: DEAD stands —
no theorem candidate printed; the regime-L failure is partly
majorant looseness (NORM|COLL| rises vs ρ₂, trend +0.086/+0.202
— exactly the r314 warning that C_full grows), so a sharper
flux-exploiting majorant is not excluded, but the kz53 obstacle
binds every tight majorant under a mid-ladder freeze; R317
direction (census-grade): the regime coordinate must capture the
near-critical family itself (TOP1/QMAX separates kz55/kz67; kz53
needs a second coordinate) or the spikes become a certified
exception family (r287-F2 pattern).  NO RH CLAIM.

r317 executed **reviewer fork (b) — the exception-family census**
(`exception_families_probe.py`, 38/38, SPEC_SHA `04fbe5c0`;
contract `PRIME.L2.RENYI3.EXCEPTION_FAMILIES.01`,
experiments-side, exploration only, no ledger row; coexisting
with the parallel r318/r319 rounds).  The statement form under
test: *all sufficiently large generic windows + explicitly
classified exception FAMILIES*, with the **hard gate sealed in
advance** — at most TWO source-pure families plus the generic
theorem; any uncovered violator, family growth or world leak
fires `WHAC_A_MOLE` by contract.  Two sealed source-pure class
functionals (AST + record-literal audited): `F_A` = rank-local
QMAX spike ratio (concentration), `F_B` = rank-local ΦL2 spike
ratio (divergence-mass spike), window W = 5, thresholds frozen by
the sealed largest-gap rule (k ≤ 6, gap ≥ 1.25) **before any
bound evaluation**.  **Verdict `WHAC_A_MOLE(kz53, kz83)`** on the
sealed letter — the hard gate fired exactly as the contract
demands: the sealed rule recovers **B = {kz55, kz67}** (the
r315/r316 FCIX pair, F_B 7.23/4.96 with a clean 1.78 gap; finite
certificates C_B = 1.0536 ≤ the r306 shallow max, ΦL2-warded,
non-growing, 0 EXT2 members), but class A stays EMPTY (best gap
1.233 misses the sealed 1.25 bar by 0.017), and on the 63-rung
complement the mid-ladder generic constant C_gen = 0.4579 is
violated by exactly the two uncovered rungs kz53 (ρ₂ 1.0493) and
kz83 (0.7791) — a third exception form would be needed; no third
class added, abort by contract, recommendation back to fork (a).
The census behind the letter is the round's find: the **F_A
top-3 are exactly kz53 (2.47), kz83 (2.39), kz67 (2.38)** — the
rank-local QMAX ratio ranks the complete mid/deep near-critical
family on top, refuting at census level the r316 conjecture that
kz53 needs a second coordinate; but the distribution below is a
**continuum** (1.93, 1.90, 1.74, …), not a gap-separated family,
and F_B over-ranks kz12 (2.79, ρ₂ 0.38 — a ΦL2 spike that is NOT
a ρ₂ spike) above kz53 (2.70): threshold classification on these
coordinates cannot cover kz53/kz83 without swallowing harmless
rungs — the exception-family FORM is the wrong statement shape
for a spike continuum.  What holds: anchors bit-near incl. the
complete r316 anatomy; the generic complement trend falls
(−0.170, reserve med 2.06); 21 small-m certificates; world
machinery intact (twin 1.04, SCRAMBLE rejected via COLL 3.69 +
edgewise shuffle break 288/300).  R318 direction (census-grade):
fork (a) with the QMAX local ratio as a *continuous* regime
coordinate (bound ρ₂ BY the coordinate, not classify by
threshold), or the generic constant at the r306 first-5 scale
where no exception family is needed on the measured ladder.
NO RH CLAIM.

**Wave 11 (2026-08-26)** froze the **architecture day** (rounds
r305–r316 incl. r310b) as **`v968_architecture_adjudication.py`**
(ledger `PRIME.LSTAR.ARCHITECTURE_ADJUDICATION.01` [E] +
`PRIME.SOURCE.RANKONE.UPDATE.IDENTITY.01` [E] +
`PRIME.SOURCE.TENSOR.MECHANISM.01` [O] +
`PRIME.L2.RENYI3.PROVENANCE.01` [O]; the TEN probes embedded
byte-exact in the sealed `--smoke` stage — SPEC SHAs
`3bb365e1`/`ec2bb008`/`d5147850`/`f8d99877`/`fac7a8df`/`6c32f749`/
`6505dd10`/`841b3196`/`92d35a3a`/`5c28b12b` — plus a module-own exact
section S0: the r306 Hill/Lagrange chain with both witnesses and the
Rényi bridge with all constants, the r314 opening-flux lemma
`F₁ == 0` with the exact flux telescope, the collision counting
identity `3p₁p₂ − 2p₃ == 8(n + 3n(n−1))` with the exact factor-8
mutant break (95 vs 760), an r311-class exact rational Farkas
mini-certificate — a PD target in the span with a *forced* negative
generator weight, strictly outside the cone, `⟨Y,Q⟩ = −2` exact —,
the TEN sealed verdict bars as exact decision logic with tipping
mutants at every live clause boundary, and the four-level
composition gate).  **The wave freezes the day in the reviewer's
binding FOUR-LEVEL STRUCTURE**: *Level 1* — real formal theorems
(Lean, already green, NO Lean change ships with the promotion:
`lstar_terminal_implies_master`, the four closed Inertia theorems,
the real PrimeWindow construction, source exactness + mass
conservation of any folding + the strong stabilization form, the
r309 rank-one update identities, the Hill monotonicity + Rényi-3⟹N₂
bridge; sorries 9 → 5, the two TRUE holes stand alone).  *Level 2* —
certified finite statements (Rényi-3 on 57 rungs at `C = 1.069`,
`A = 2`; the fixed-head census 77/77; the block-Green identity +
the strict cone with exact rational Farkas certificates; paired-cone
soundness on 20714 steps; fold multiplicity exactly 2; 93 %
far-triple cancellation; the signed cubic identity with vanishing
boundary term, Fractions-exact; the 21 small-m certificates).
*Level 3* — negative architecture decisions verbatim
(`FIXED_HEAD_DEAD`; `PAIRED_CONE_NO_INDUCTION` +
`PAIRED_CONE_WORLD_BLIND` — B does NOT overtake A; the r308
discriminator = the budget sign / a chordal restatement;
`COEFFICIENT_SIGN_WALL` — Lane A closed as cone language with the
mechanism *named*: block-psd membership without rank-one SDD,
obstructed at the budget/border sign; `FLOQUET_EXPANDING`;
`TRIPLE_TYPE_MAJORANTS_WRONG`; `PHI3_ALL_BLIND`; `TWO_REGIME_DEAD`
with the anatomy — the obstruction family cuts *across* the FCIX
stratum, the r306 first-5 protocol was load-bearing, kz55/kz67 are
near-one-block worlds).  *Level 4* — open mechanisms typed (the
constructive PSD-G rule as a question, the signed cubic flux
theorem, the quantitative source functional, the `crossing_budget`
mathlib gap, the definitional source bridge, L\*, the terminal).
**The claim split is binding and NOT one TRANSFER title**:
`RANKONE.UPDATE.IDENTITY` [E] carries only the exact update
identities (determinant dictionary, Sherman–Morrison chain,
±(η±e_t) border split, local Schur reserve, signed source updates);
`TENSOR.MECHANISM` [O] carries the open question, registered *after*
r313/r316 with the documented negatives.  Universality-class
framing: the base is fine-metric (`METRIC_ONLY`, twin invariance
everywhere), the fiber is recursive (EPSTEIN holds, SCRAMBLE breaks)
— shared algebra, not necessarily a shared mechanism; the reviewer
maturity raise 4 → 5 (formal architecture 9.5) is typed as a
**reviewer assessment**, not a machine-checked fact.  In Lean **no
change ships** with this promotion (r305/r310/r310b are already
merged and green).  NO RH CLAIM.

**Lane status after wave 11.**  **Lane A (block-Green / constructive
PSD-G): CLOSED as a cone language** — r311 `STRICT_SOURCE_CONE`
(exact rational in-span Farkas certificates; the r308 world
separation demarcated as the budget sign) and r312
`COEFFICIENT_SIGN_WALL` with the mechanism named (block-psd
membership without rank-one SDD, obstructed at the budget/border
sign); the PSD-G rule stays open as a *question*.  **Lane B
(paired-cone / rank-one stream): CLOSED, identities banked** — the
r309 B1–B4 adjudication (B does not overtake A) froze the exact
rank-one update identities per reviewer case 3
(`PRIME.SOURCE.RANKONE.UPDATE.IDENTITY.01` [E]); the mechanism
question lives in `PRIME.SOURCE.TENSOR.MECHANISM.01` [O].  **Lane C
(fixed head): PARKED** — `FIXED_HEAD_DEAD` (r307), the one-atom
anatomy is a coordinate, not a dimension.  **The fiber: the r306
Rényi-3 GO stands** (57/57, growing reserve), its **provenance is
open** (`PRIME.L2.RENYI3.PROVENANCE.01` [O]; trilogy: identity exact
(r314), functional blind (r315), two-regime dead (r316)); the
**R317 fork sits with the reviewer** (near-critical family
coordinate OR certified exception family) — r317 executed fork
(b) and sealed `WHAC_A_MOLE(kz53, kz83)`: the spike-family
census is real (B = {kz55, kz67} recovered blind; the F_A top-3
are exactly the near-critical family) but the tail is a
continuum — back to fork (a).  **Lean: 5 sorries, two
TRUE holes** (`lstar_subordination`, `terminal_positive_main`; plus
the `crossing_budget` mathlib gap, the definitional source bridge
and the pair-coordinate refinement).  The mincut (base 4 / refined
5) is unchanged.  NO RH CLAIM.

r318 ran the **indefinite-fork round — the new L\* idea class**
(`indefinite_fork_probe.py`, 25/25, SPEC_SHA `f2d98683`; reviewer
perspective shift after the r312 anatomy — "the signedness itself
is the mathematical object that is not yet understood": ONE
theoretical representation test, Pontryagin/Krein **index** vs
**sign regularity**, under the binding stop rule (both routes
world-blind/restatement ⇒ `FORK_DEAD`, lane stop);
experiments-side, exploration only, no ledger row).  **Verdict
`P2_MAIN_SPECIFIC` (fingerprint)**: the 2.4 % negative cross
mixtures of the converged r308 Dykstra family are **not
structureless** — on MAIN-class worlds the block residual lives
lawfully on the **antiphase pair (D3, D4) with fixed sign −1**
(modal share med 0.699, 12/12 sealed rungs, the r289 rational
twin exact — `METRIC_ONLY` holds), while the dead controls put
their residual on the **arch-mean × border pair (D5, D6)**
(shares 0.953/1.000, d6-class 0.962/1.000 vs MAIN 0.027) — the
fingerprint of the signedness is world-separating **in shape**;
honest caveat typed with the letter: the control fingerprints are
measured on non-psd-feasible 200-step iterates (the r308
non-convergence), labeled ITERATE.  **The P1 index route closes
as language**: the index bookkeeping is exact (spectral count ==
mp pivot count on all seven worlds, window + guarded deep depth;
window index defect 0/0/0 on w9/w13/twin vs 55/37/4/31 on
EPST/SCR/SMOOTH/HL2 — the control flips *are* negative directions
inside the window; ladder 42/42 defect 0), but the index
statement "n₊(window) == N_w" is a total **restatement** of L*
(equivalence with minC ≥ N_w on 3 exact toys + 7 real worlds,
both truth values realized), the global inequality n₊ ≥ N_w is
**vacuous** (carried by the μ-channel majority on every world),
and the negative-subspace invariants are all world-blind (and
not crossing-location proxies — no lever either way).  The SR
half of P2 is an honest negative: no census object is
MAIN-specifically sign-regular (orientation-sensitive minors are
coin-flip ~0.50 even on MAIN) ⇒ the variation-diminishing
implication chain to L* fails at its premise
(`PREMISE_FAILS_ON_MAIN`).  Dig site named per the stop rule: why
the antiphase (D3, D4) cross mixture carries a fixed negative
sign law on MAIN-class worlds.  Protocol disclosures: prefilled
placeholder record tables removed completely before the first run
(the r316 error class); one smoke-stage harness bugfix; one
reporting-only amendment a1 (dig-site label derived from the
measured modal pair).  NO RH CLAIM.

r321 ran **reviewer fork (a) in the r317-sharpened form — the
continuous-coordinate round** (`continuous_coordinate_probe.py`,
39/39, SPEC_SHA `e68883ad`; contract
`PRIME.L2.RENYI3.CONTINUOUS_COORDINATE.01`, experiments-side,
exploration only, no ledger row; coexisting with the parallel
r318/r320 rounds).  The statement-shape change the r317 abort
demanded: **bound ρ₂ BY the single source-pure coordinate**
F_A = rank-local QMAX ratio (r317 import, W = 5) — ρ₂(w) ≤
g(F_A(w)) with g explicit and monotone — instead of classifying
by threshold.  Sealed before evaluation: the exact concentration
bracket qmax·ΦH1 ≤ ρ₂ ≤ ΦH1 = (F_A·B)² with the source-pure
baseline B = medloc × m/log m (warded two-sided live on 69
worlds); a three-form monotone g-family with mid-ladder
calibration and derivation-strength precedence (G_SQ = b·F² with
b = (max cal B)², consuming **no target value**; G_TT = c₁ +
c₂·F³, the r317-named heuristic; G_LIN = a·F); a fit-free
upper-envelope rule on the declared test split; gain bar 1.5; a
four-branch verdict tree.  **Verdict `SLIDING_BOUND_GO(G_SQ)`**:
ρ₂ ≤ 1.3056·F_A² pointwise with **0/39 mid-ladder test
violations and all four named r316/r317 violators INSIDE**
(kz53/kz83/kz67/kz55 at reserves 7.0–9.6 — exactly what killed
every flat mid-ladder constant since r316 is absorbed by the
coordinate); reserve min/med 2.71/5.35, trend −0.341 falling;
the test envelope is **strictly monotone in F_A** (bin Spearman
+1.000, top bin 1.0536), bulk Spearman +0.84 — stronger than the
r317 continuum reading; gain 15.95.  Candidate theorem (sliding
cubic bound) printed: Σq³ ≤ 1.3056·F_A(w)²·(log m)²/m² for m ≥
73 + 21 small-m certificates (C_small 1.0694); corollary: F_A ≤
2.47 measured ⇒ uniform C_impl = 7.97 (disclosed 7.5× looser
than the r306 first-5 constant — the round buys FORM, not
sharpness).  Two honest structural finds: (1) the pure-algebra
transfer route does **not** close — B is not bounded by its cal
max (test max 1.4088 = 1.23×, trend +0.122 rising); the G_SQ
certificate holds on the measured ρ₂ directly because the
bracket's qmax slack falls faster than B² rises; (2) SCRAMBLE is
**not** rejected by the coordinate bound alone (F_ins 2.00,
covered at 5.21) — the rejection is carried by the r317 class
side condition (COLL 3.69, shuffle 294/300), honestly typed.
The remaining provenance question is now source-pure, local and
split in two: is F_A bounded (max 2.47), and what bounds the
qmax-share ρ₂/(F_A²·B²)?  NO RH CLAIM.

r322 ran the **antiphase-sign-law round — the dig after the r318
fork decision** (`antiphase_sign_law_probe.py`, 25/25, SPEC_SHA
`761b51d4`; contract `PRIME.LSTAR.ANTIPHASE_SIGN_LAW.01`,
experiments-side, exploration only, no ledger row; coexisting
with the parallel r320/r321 rounds).  The three dig questions:
hardening (is the (D3, D4) sign law of the block-Green
cross-mixture residual iteration/construction-invariant, and
what do *fair* control objects carry — the r318 caveat was that
control fingerprints were ITERATE-grade), the law itself (is it
identity-forced — exact Fractions forced-component analysis of
the linear identity constraints on the (D3, D4) sector), and the
theorem candidate ("psd cone part + antiphase-supported
fixed-sign residual").  **Verdict `ALGORITHM_ARTIFACT`** on the
sealed tree: the law breaks under the sealed random-start
variants — their near-feasible points (staged to 20 000 Dykstra
steps, min eig rel −8.5e-9/−1.2e-9) carry modal pair (0, 2) with
shares 0.31/0.26 and cone share 0.78: **the sign law and the
r312 97.6/2.4 cone anatomy are properties of the
least-norm-proximal Dykstra basin, not of the psd solution set
as a whole** (the canonical variants LSTART/LONG/ZERO/REV all
carry the law exactly and coincide within rel 2e-7 — exactly the
9/15 variant pairs involving the random starts are distinct).
The r318 caveat resolved on fair objects: the least-norm census
reproduces the r318 shape contrast construction-fairly, **but
the budget-ablated converged controls split** — EPST-abl carries
the MAIN law *stronger than MAIN* (share 0.742 vs 0.692), SCR-abl
breaks (0.379) ⇒ `FAIR_CONTROL_CARRIES` 1/2: the r318 pattern
dichotomy does not survive fair convergence.  The law is **not
identity-forced** (TOY4 calibrator machine-exact 4/7; SM1 exact
free fraction 0.477 with *positive* forced value; MINI16
rowspace membership false) — the negative sign law lives
entirely in the selection geometry, as the D3 = −(D2 + 2 D1)
in-block dependence predicted.  Anatomy: total sign consistency
(1.000), position-structured (mid/tail thirds 0.29/0.86/0.93),
scale-stable magnitudes; best coupling correlate: the
block-local antiphase mass pairing (+0.81, below the sealed 0.9
bar).  The candidate form certifies 57/57 rungs + twin censally
but is demoted by the basin finding and the sharpened L5 gap.
Protocol: draft record block removed before the first run; one
amendment a1 (gate typing only, letter identical across passes).
NO RH CLAIM.

**Wave 12 (2026-08-27)** froze the **red-team / fork morning**
(rounds r317–r322, incl. the r319 extraction-chain audit and the
r320 Lean repair) as **`v969_forks_and_redteam.py`** (ledger
`PRIME.REDTEAM.EXTRACTION_AUDIT.01` [E] +
`PRIME.L2.RENYI3.SLIDING_BOUND.01` [O] + updates to
`PRIME.L2.RENYI3.PROVENANCE.01` [O] and
`PRIME.LSTAR.SUBORDINATION.01` [O]; the FOUR probes embedded
byte-exact in the sealed `--smoke` stage — SPEC SHAs
`04fbe5c0` (r317 `exception_families`, 38/38) /
`f2d98683` (r318 `indefinite_fork`, 25/25) /
`e68883ad` (r321 `continuous_coordinate`, 39/39) /
`761b51d4` (r322 `antiphase_sign_law`, 25/25) — plus a module-own
exact section S0: the exact two-sided concentration bracket
`qmax·ΦH1 ≤ ρ₂ ≤ ΦH1` re-proved from scratch with both mutant
directions caught, the U1–U3 red-team witnesses in exact
arithmetic (U1: the `B = −1` window against the terminal type;
U2: the mesh-level-0 total node collision reduced to three
rational facts about `e` on the certified bracket
`685/252 ≤ e ≤ 685/252 + 1/35280`, kill `0 < 0` exact; U3: the
adversary `(|F|+1)² < 5/7` false for every F), the r320
separation-discipline and canonical-split certificates (witness
mesh level 4 via `2⁷ < 3⁵ < 2⁸`, reduction `1/8 < 2/5` exact;
level-0 failure `4/3 < 4`), six sealed verdict bars with tipping
mutants, and the wave-12 composition gate; NO Lean call — the
r319/r320 Lean rounds are consumed as REPORTS, their artifacts in
`rh/lean/` re-verified by `lake build` + `run_rh.py`).  **The wave
adopts the reviewer's binding adjudication**: *(1)* the proof
graph has **three cleanly separated levels**, and the two true
holes are named **WINDOW-LOCAL** on every surface from now on —
**Level A** (local window mathematics, PROVED:
`L*(w) + T(w) ⟹ MasterPositive(w)`, cleanly typed per r305/r320;
"for every canonical finite window, L\* and terminal positivity
imply augmented prefix positivity — this is a window-local
theorem"), **Level B** (the universal window family — THE two
arithmetic *window-local* holes: `∀ w ∈ W_canonical: L*(w) ∧
T(w)`), **Level C** (extraction to Weil/RH, OPEN: source-exact
realization of the actual source + admissible transport +
convergence/stabilization in the correct directed order +
density/continuity for the Weil criterion).  *(2)* **The false
cofinality direction is documented, not solved** — and the sorry
census does **not** fully measure the global distance: the missing
Level-C extraction architecture appears in *no* `sorry`, because
the correct theorem statements do not yet exist.  *(3)*
`SourceExact` is a good firewall but must **not** stand as a free
assumption in a final proof — the named open Lean target is the
target architecture **`CanonicalPrimeWindow` + exactly one
construction theorem `sourceExact_buildPrimeWindow`**
(load-bearing for the global connection, not formalization
polish).  *(4)* the typed per-layer **reviewer assessment**
(supersedes the wave-11 4 → 5 scalar; an assessment, not a
machine-checked fact): finite dictionaries 9.5/10, formal
window-local architecture 9/10, terminal mechanism 6/10, L\*
mechanism 2/10, source/extraction architecture 2–3/10, complete
RH proof path ~3/10, audit quality 9/10 — *"stronger as a
window-local research program; as a complete RH proof chain
currently less far than previously assumed — the removal of a
shortcut that did not exist."*  NO RH CLAIM.

**Lane status after wave 12.**  **L\*: internally PAUSED**
(allocation 0 %, reviewer decision; the selection geometry is not
pursued further — the r322 `ALGORITHM_ARTIFACT` closed the
antiphase fingerprint as a Dykstra-basin property; reopening only
under the six documented conditions of the r322 round notes; as an
*external* object the two-weights sampling formulation
`‖J_N‖ < 1` stays the memo candidate, `rh/problem/`).  **The
terminal: sharpened to the QMAX × M₂ shape route** (the sliding
bound's two-part demand: is F_A bounded / what bounds the
qmax-share; **round R324 in flight — not consumed here**).  **The
extraction layer: the repair fork R325 measured** (primary verdict
`ELEMENTWISE_STABILIZATION_GO` — the repair is a quantifier set,
not a new mesh theory; see the r325 paragraph below; the Level-C
target `CanonicalPrimeWindow` + `sourceExact_buildPrimeWindow`
stays the named open Lean goal).  **The fiber: the sliding cubic bound stands as
the successor contract** (`PRIME.L2.RENYI3.SLIDING_BOUND.01` [O];
`Σq³ ≤ 1.3056·F_A²·(log m)²/m²` on the measured rungs, all four
named violators inside; C_impl 7.97 disclosed 7.5× looser — form,
not sharpness).  **The base fork: closed honestly** (r318 P1
banked as language — a total restatement, no lever; classical sign
regularity dead on MAIN).  **Lean: 5 sorries, the two true
WINDOW-LOCAL holes** (`lstar_subordination`,
`terminal_positive_main` — Level B; plus the `crossing_budget`
mathlib gap, the definitional source bridge and the
pair-coordinate refinement; census 5 → 5 across the r320 repair
with two typed retypes).  The mincut (base 4 / refined 5) is
unchanged.  NO RH CLAIM.

r324 ran the **QMAX × M₂ origin round — the reviewer-sharpened
terminal main lane after the r321 sliding-bound GO**
(`qmax_m2_origin_probe.py`, 36/36, SPEC_SHA `dc36cacb`; contract
`PRIME.L2.RENYI3.QMAX_M2_ORIGIN.01`, experiments-side, exploration
only, no ledger row; coexisting with the parallel wave-12 worker
and R325).  **Binding course correction mid-round, disclosed**: the
round's original F_A-provenance contract was executed and sealed
first (`fa_provenance_probe.py`, 36/36, SPEC_SHA `9a6696f8`,
r324-pre in the inventory — letters `FA_BOUNDED_DISTRIBUTIONAL` +
`QMAX_SHARE_OPEN`), then the reviewer adjudication demoted the
F_A-boundedness hunt as a potentially unnecessarily strong
statement — for the needed power N₂ ≳ m^0.888 it suffices that
F_A·log m = O(m^{0.112−ε}); the cleaner decomposition is
**QMAX × M₂**; the banked pre-work is reproduced bit-near as
anchors.  **S0 (exact)**: the interpolation `M₃ ≤ q_max·M₂` and
`M₂ ≤ q_max` (Fractions-proved, live-warded at slack 0.0), and the
identity `F_A·B·log m == m·q_max` (2.7e-16) — F_A is the rank-local
normalization of m·q_max, no black correction factor; the r321
bracket upper is improved by exactly the factor M₂/q_max ≤ 1.
**S1 (M₂ export)**: the pointwise export `m·M₂ ≤ C₂ = 2.2557`
(mid-ladder freeze) **fails 7/39** exactly on the banked violator
set (kz53/kz67/kz83 at 3.05–3.19) — the reviewer's "export should
carry" expectation is refuted pointwise (trend +0.014 flat,
envelope 3.1938; the distributional r302 export holds, KS
0.0230/0.0158).  **S2 (the multiscale pileup hypothesis, four
measurements)**: scale recomposition exact (2.2e-16); fold
multiplicity ≤ 2 on all 69 live worlds; per-scale mass C_PIL
9.3583 **fails 11/39** (clause refuted); active scales
nsc_rel/log m ≤ 2.026 **holds 0/39** (the O(log m) scale count
certifies); direct C_INF 1.7481 fails 5/39.  **Hardness census
answered**: the named near-critical spikes are *not* multiscale
convergence — one near-maximal source scale dominates the argmax
block (kz67 s0 = 22.4 of G 13.6; nsc at the ladder median): the
pileup tip is a **single heavy scale** — the reviewer hypothesis
refuted in clause (3), confirmed in clause (4).  **Verdict
`PILEUP_GROWS_SUBCRITICAL(+0.172)`** on the sealed six-letter
tree: e(G/log m) +0.158 + e(m·M₂) +0.014 = +0.172 < CRIT 0.224,
stability halves +0.141/−0.160 both subcritical — decided; by the
reviewer contract a sufficient first-class outcome.  **S3
composition (typed MEASURED)**: Σq³ ≤ 8.941·(log m)·m^{+0.172}/m²
⟹ N₂ ≥ N₃ ≥ m^{0.914}/√(8.941·log m) ⟹ N₂ ≥ m^0.888 for all
m ≥ m₀\* = 10^59.6; the measured ladder ends at m = 274 — the gap
is the disclosed extrapolation hypothesis, no cofinal claim.
Protocol: course correction documented in the docstring; the
banked pre-work outcomes disclosed pre-spec; one calibration-pass
solver fix before the record freeze (m₀\* search moved to log
space, sealed rule unchanged); deterministic run1/run2.  NO RH
CLAIM.

r325 ran the **extraction-repair fork — the new main strand after
the R319 red team** (`extraction_order_probe.py`, 18/18, SPEC_SHA
`31277f91` final with record / freeze-run `57e50d36`; contract
`PRIME.EXTRACTION.ORDER_REPAIR_FORK.01`, experiments-side,
exploration only, no ledger row; coexisting with the parallel
wave-12 worker and R324).  The round is mixed analytic + machine:
the three reviewer variants for repairing the RH-connection seam
(elementwise stabilization vs. canonical mesh rebuild vs. explicit
quadrature correction) were adjudicated on their machine-checkable
halves under a sealed tree.  **Verdict: primary
`ELEMENTWISE_STABILIZATION_GO`, co-letters
`CANONICAL_MESH_REBUILD_GO` +
`SIGNED_DEFECT_CORRECTION_GO(SIGN)` +
`ANCHOR_ONLY_INSUFFICIENT`.**  *(A)* On the native dense class
(dyadic step-function autocorrelations — the v749 "Weil form of
step functions" class, seed-sealed) **all three channels**
(comb/arch/pole) of the canonical tower windows stabilize
**exactly** at the predicted finite anchor onset α\* =
(n_g+1)·D₀/2 (comb rel dev 0.0/3.3e-18 above onset, 2.7e-2–5.3e-2
below — a real onset) and the value is **constant under dyadic
mesh refinement** (comb 3.3e-18, arch 1.5e-15, pole 2.0e-17): the
elementwise architecture consumes **no mesh-cofinal ladder and no
transport** — the R310 `finite_forms_converge_to_weil` shape
generalizes from the Lean comb channel to every channel on the
native class (derived tent-read identity `L_cat(F) = −Σ μₙ
(I_D F)(uₙ)` gated at 1e-12).  *(B)* v749 T1.4b reproduced: every
deployed frame-A window **is** a directly rebuilt tower member at
rel dev 0.0; the rebuild ladder is PD per stage with honestly
falling margins (λ_min/c₀ 1.33e-5 → 6.92e-7, reported not
consumed); the dyadic transport identity holds at 8.0e-17 but
variant B never uses it.  *(C)* On the sealed non-native C² class
the per-channel quadrature defect is **one-signed** (comb
negative, arch positive, pole positive, all stages), falls at the
D² class rate per channel (ratios 3.6–9.0; the *total* ratio
2.2–3.3 is reported only — comb/arch sign opposition cancels
there), and sits inside the **rigorous closed-form interpolation
envelope** |E_ch| ≤ (D²/8)‖F″‖·K_ch — controlled on its own, no
wall margin anywhere; the anchor-only floor is real (E_tot
constant +1.77e-4/+1.59e-4 across the anchor sweep at fixed mesh —
the R319 false-floor mechanism on this probe's own elements).
Fork consequence (analytic deliverable of the round): the repair
statement for the next Lean round is the **elementwise quantifier
set** — `CanonicalPrimeWindow = {w | ∃ a h, w = build(a, h)}` with
exactly one construction theorem, per-element predefined finite
stabilization on the native class, and the window-local positivity
premise quantified over the canonical family; **no new mesh theory
and no H_cof ladder**.  The mesh limit survives only as the
classical density/continuity step (the v912 topology caveat),
quantified by leg C.  Nothing here proves any window positivity or
any local lemma.  NO RH CLAIM.

r327 ran the **fold-group mass-cap round — the r324 follow-up on
the terminal main lane** (`group_mass_cap_probe.py`, 34/34,
SPEC_SHA `11e4fd40` final with record / freeze `71f8b7b4`;
contract `PRIME.L2.RENYI3.GROUP_MASS_CAP.01`, experiments-side,
exploration only, no ledger row; coexisting with the parallel
r326 Lean round).  The r324 hardness census had located the
near-critical spikes in a *single heavy source scale*; this round
descends one level: **what caps the mass of the dominant fold
group?**  New exact machinery (warded live on 69 worlds at slack
0.0): the per-group mass ledger on the identical r314 fold
segmentation, the partition `Σ_g gabs == A1_j`, the L1
recomposition, the **two-ancestor bound** `gabs ≤ mult·gmax ≤
2·vmax` (the proven mult-2 cap as a mass statement) and the exact
cap chain `m·q_max ≤ ng·hgn`.  **Q1 anatomy answered — sealed
letter `SOURCE`** (med named ratio 1.057 ≥ 0.5): the r324 heavy
scale *is* a single fold group at the two sharpest spikes, and the
structural finding of the round: on **all 65 rungs** the heaviest
group of the argmax block is **exactly one β/ω fold pair** (one
bulk + one window atom at one position; window-atom histogram
[0, 65, 0]) with median alignment 1.000 — the two ancestors
reinforce; at kz53 the single pair carries 88.8 % of the argmax
block's atom mass (second group 13× lighter): the spike is
literally **one bulk/window coincidence** (kz67 is the exception
shape: 6 groups, bshare 0.416).  **Q2 caps — verdict
`CAP_PARTIAL`** on the sealed six-letter tree: the λ-pair route
`hga ≤ 2·vmaxb` is exact but the source-size coordinate
`2·m·vmaxb/L` **fails the mid-ladder freeze 19/39** (C_LV 1.1838;
kz55 6.85 / kz53 5.77, named 0/4) — the von-Mangoldt size analogy
does not close as a polylog cap; the direct group cap fails at
A = 1 (18/39) and A = 2 (8/39); **the group count certifies**:
`ng(argmax block) ≤ 2.6351·log m` at 0/39 + named 4/4 — the third
certified O(log m) count of the lane.  Exponents: e(hgn/log m)
**−0.200 falling**, e(lvb/log m) −0.261 falling (the heavy-group
coordinate decays with depth — the spike set is a shallow/mid
phenomenon, not growth), e_route +0.106 < 0.224 but halves
+0.411/−1.183 straddle — `HEAVY_GROUP_UNBOUNDED` does **not**
fire.  **Q3**: no certifying cap ⟹ the r324 MEASURED composition
stands unchanged (m₀\* = 10^59.6); slim anchors bit-near (r306 /
r316 / r324 C_NSC 2.0258 · C_PIL 9.3583 · C_INF 1.7481 with
violation counts 0/11/5 exact / r324-pre C_M2 + seven-violator
set exact).  Must-fails e1 group-drop LOUD / e2 wrong-window-max
exact / e3 cap-posthoc AST+toy / e4 qmax-readback AST all CAUGHT.
Deterministic run1/run2; no amendment before or after the freeze
except the record-table insertion.  NO RH CLAIM.

**Wave 13 (2026-08-27)** froze the **extraction repair and the
terminal composition** (rounds r323–r327, incl. the r326 Lean round
and the r323 clean abort) as
**`v970_extraction_and_composition.py`** (ledger
`PRIME.EXTRACTION.ELEMENTWISE.01` [E] +
`PRIME.L2.RENYI3.MEASURED_COMPOSITION.01` [O] + updates to
`PRIME.L2.RENYI3.SLIDING_BOUND.01` [O],
`PRIME.L2.RENYI3.PROVENANCE.01` [O],
`PRIME.REDTEAM.EXTRACTION_AUDIT.01` [E, Lean note] and
`PRIME.LSTAR.SUBORDINATION.01` [O, Lean note]; the FOUR probes
embedded byte-exact in the sealed `--smoke` stage — SPEC SHAs
`9a6696f8` (r324-pre `fa_provenance`, 36/36 — banked) /
`dc36cacb` (r324 `qmax_m2_origin`, 36/36) /
`31277f91` (r325 `extraction_order`, 18/18) /
`11e4fd40` (r327 `group_mass_cap`, 34/34); r324 imports the
embedded r324-pre, r327 the embedded r324-pre AND r324 — plus a
module-own exact section S0: the r324 interpolation
`M₃ ≤ q_max·M₂` with slack exactly 1/36 and double equality at the
flat profile, the F_A identity `F_A·B·log m == m·q_max` as exact
constructional algebra (gauge mutant caught), the r327 group
ledger (partition, cap chain 12/7 ≤ 20/7, two-ancestor bound with
both mutants breaking by exactly 1), the r325 onset formula
`α* = (n_g+1)·D₀/2` as an **exact toy stabilization** (read ==
atom sum 35/16 for every α ≥ 5/8, constant under mesh refinement,
nearest-neighbour mutant caught), the m₀\* log-space solver
rebuilt verbatim (landing 10^59.6; the disclosed pre-freeze
linear-space bug rebuilt and caught), the r328B chain-audit gate
S0-T5d (block-level need 454/501, chain bar 0.188, verdict
survival 0.172 < 0.188 exact, chain-honest m₀\* = 10^238, the
false-fire window mutant caught), six sealed verdict bars with
tipping mutants, and the wave-13 composition gate; NO Lean call —
the r326 Lean round is consumed as a REPORT and the r323 abort as
a NOTE).  **The honest formulation of the new connection (binding,
verbatim)**: *two window-local lemmata (over the canonical family)
+ the elementwise architecture + cited classics (density, v912
continuity, the Weil criterion) ⟹ RH; open links: the arch/pole
transcription (classical), the compression bridge, the source
completion (definitional).*  NO RH CLAIM.

**The adversarial audits (r328A/B/C, read-only).**  Three
independent adversarial audits ran on the wave-13 state on the
same day; reports are with the coordinator (the full booking goes
into `next.txt` by the coordinator).  **A (methodology)**: no
verdict-tipping bug; a verifiability gap and a sequence-distortion
risk named.  **B (chains)**: found the one MODERATE transcription
error of the wave — the r324 composition bar CRIT = 0.224 was
derived from the ATOM-level need 0.888, while the composition
functional N₂ = n_eff is BLOCK-level (need 0.908 in N units,
908/1002 ≈ 0.9062 in m units), so **the chain-honest bar is
CRIT = 0.188 and the chain-honest m₀\* is 10^238** (the record
10^59.6 belonged to the atom-level need).  **The r324 verdict
SURVIVES** (+0.172 < 0.188); the sealed record is NOT edited — the
correction is carried as a chain annotation (ledger, the v970 gate
S0-T5d, the rh-paper S3′ item) and the going-forward bar is 0.188.
Otherwise clean, incl. two independent recomputations.  **C (blind
spots)**: the ladder selection effect, the twin resolution limit
and the missing second living world named.

**The EXT3 fresh anchors (r329, audit repairs 2+3 of r328A/C).**
`ext3_fresh_anchors_probe.py` (33/33, SPEC_SHA `bbfaf1990bd37061`,
under the two-commit freeze protocol: pre-freeze `604b8aaf`, record
`8cbd95f9`) repaired the sequence bias of the r321/r324 lane: 12
genuinely NEW windows (machine-checked disjoint from the 80-kz used
ledger of every prior round, all deeper than every prior window,
N_w 1721..2577), built with the IDENTICAL frame-A construction and
only the anchor CHOICE changed — 6 deliberately SMALL-GAP anchors
inside the core zone range (grel 0.28..0.53 vs ladder median 1.32;
the ladder is large-gap selected by construction, audit C1) + the 6
deepest fresh windows — tested pointwise against the four FROZEN
constants, nothing recalibrated.  **Verdict `EXT3_VIOLATIONS` +
`SMALL_GAP_DIFFERENT`(ρ₂ med ×10.29, F_A med ×3.21)**: the r306
uniform constant C = 1.069 BREAKS on the four sharpest small-gap
anchors (kz51 by 7.1×) — it is a large-gap-family statement; the
r321 sliding bound 1.3056·F_A² covers exactly those spikes with
reserves 5.0–6.1 at the frozen constant (the coordinate carries
genuine out-of-sample evidence) but misses four QUIET anchors
(F_A 0.27–0.39, excess ≤ 1.40×, convention-sensitive: inside under
the insertion rule) — the named successor question is a floor term
below the ladder's sampled F_A range; the two O(log m) counting
constants C_NSC 2.0258 and C_NG 2.6351 HOLD 12/12 (min reserves
1.6/1.8) — the counting side is the lane's most choice-robust
asset.  Four stratum-A candidates were the first POSITIVE_PREFIX
failures ever seen on this construction.  Experiments-side, NO
ledger row, NO RH CLAIM.

**The second living world (r330, audit repair 4 of r328C).**
`dirichlet_secondworld_probe.py` (24/24, SPEC_SHA `665260180de9`,
the first round under the two-commit freeze protocol: pre-freeze
`2b1b69ce`, record `4f2e4f52`) built the missing POSITIVE control —
a Dirichlet-L comb (primitive real χ mod 3; weights `χ(n)·MU_ALL`
bitwise in the v563 gauge, plus the milder magnitude variant
DIRICHLET_ABS without the p = 3 atoms) through the IDENTICAL frozen
machinery on the identical 42-rung ladder, nothing recalibrated.
**Verdict `SECOND_WORLD_SPLITS`** along a clean wall/cascade seam:
SHARED (genuine living-world properties, now n = 2 living
arithmetics vs 2 dead controls) = the r300 FILL class (med 0.436
FILL_LOW vs EPST/SCR FILL_HIGH), the n_eff growth direction
(+0.527) and the mult-2 one-pair fold geometry (42/42, mixed share
1.00); SPLIT (retyped as zeta-side idiosyncrasies INSIDE THIS
FRAME) = the half-filling wall itself (DIR nf None 0/42, depth
ratio med 0.080 — and sharpest: ABS, the mere REMOVAL of the three
p = 3 atoms, already kills the wall at control depth 37 at w9),
the r306 pointwise C₂ bound (32/42 violations, max 27.6×), the
σ decay (the diagonal GROWS, +0.793) and the r299 O-sign class
(O_NEG 27/42).  HONEST LIMIT (binding): the frame keeps zeta's
archimedean term/border/depth — the round cannot distinguish
"wall is zeta-specific" from "wall needs the matched
Γ/conductor frame"; the GRH-faithful Dirichlet window is the named
follow-up.  Experiments-side, NO ledger row, NO RH CLAIM.

**The twin resolution (r331, audit repair 5 of r328C).**
`twin_resolution_probe.py` (24/24, SPEC_SHA `f871fe84acdc0526`
final with record, freeze `8d771fae0a374fa8`, two-commit protocol:
pre-freeze `b8319121`, record `a45d4560`; seed-free deterministic)
measured the METRIC_ONLY reach the audit flagged as possibly
untested: is the r289 rational twin resolution-capable at the
r286 tight anchors (margins down to 1.8e-8), or is its own
construction disturbance larger than the margin there?  **Verdict
`METRIC_ONLY_CONFIRMED_AT_RESOLUTION` (margins 1.675e-4 ..
1.806e-8 over 11 anchors)**: at the r289 precision (1e-8 of the
local gap) the twin footprint `|Δ margin|` lands **3.4+ orders
below EVERY margin** (11/11 RESOLVE, worst ratio 3.7e-4 at kz98,
mp dps-30/45 pair-verified on every sub-1e-5 margin) — the feared
blindness does not occur, because the margin responds to the
rationalization far more weakly than 1:1 in gap units; the
precision sweep (1e-3/1e-6/1e-8/1e-10 on kz9/64/98/119) shows the
twin is BLIND at 1e-3 everywhere (signature broken) and
twin-stable from **t\* = 1e-6** down on all four anchors (kz98
the closest call at ratio 9.9e-2); the genuine diophantine
counter-test (golden-ratio detune, weights bitwise, gap structure
preserved ≤ 1e-4, the 28 exact w9 2-power resonances demonstrably
destroyed 28 → 0) **keeps the wall at every resolving dose** on
w9 AND on the tightest anchor kz119 (Δ +4.7e-11 at margin
1.8e-8) — METRIC_ONLY doubly confirmed, now WITH a resolution
certificate on the tight range instead of an untested
extrapolation.  Honest limits: the f64 construction floor caps
twin precision at ~1e-10 gap (floor censuses disclosed);
near-resonance structure untested; 11 anchors × 4 tolerances is
a finite census, not a family theorem.  Experiments-side, NO
ledger row, NO RH CLAIM.

**The β/ω companion-orbit packing (r333, reviewer rank 1 after the
r327 carrier).**  `companion_orbit_packing_probe.py` (34/34,
SPEC_SHA `02b22672e1c3f91d` final with record, freeze
`5567af41c0d4bc82`, two-commit protocol: pre-freeze `279ba0dd`,
record `55d4a832`) implemented the reviewer's packing theorem
candidate: to every dominant β/ω pair g\* a companion set O(g\*)
constructed CANONICALLY from its source data (same dyadic gmax bin
+ same fold parity + the maximal inner interval, BOUNDARY_F 0.10 —
NOT a similarity search, the e4 mutant proves it), with clauses
(T1) |O| ≥ c₀·m/(log m)^a, (T2) weight comparability, (T3)
overlap ≤ 2.  **Verdict `ORBIT_TOO_SHORT`**, with TWO
theorem-grade positives and ONE structural discovery: (T2) the
BAND DISCIPLINE certifies EXACTLY — inside one dyadic gmax bin the
two-ancestor bound forces `1/4 < gabs(g)/gabs(g') < 4` (derived
algebra, sharp near 1/4; live rmin worst 0.2819, cband worst 3.55
on 69 worlds): b = 0, c₁ = 1/4, the r327 λ-pair failure repaired
at band level; (T3) overlap is partition-exact (membership 1 ≤ 2);
the composed chain `q_max ≤ ng·κ/(|O|·rmin)` is exact live (κ =
L1_atom/L a stable ~3–4 factor, NOT the blocker).  But T1 fails at
all a ∈ (0,1,2) and the composed bound fails pointwise (named
coverage 0/4), because THE BOUNDARY CASE IS THE RULE: on 44/65
rungs the dominant pair sits in the outer 10 % margin of the
masked support (d\* med 0.072, min 0.002); the inner-interval
orbit is EMPTY at kz53/kz83/kz55 (|O| = 0; kz67: 2) while the
one-sided census shows same-key companions DO exist on the edge
side (1 each at the spikes, up to 25 elsewhere) — the spike
coincidence is an EDGE phenomenon; all four named spikes are
SMALL_ORBIT + HEAVY_BAND simultaneously.  Subcritical fallback
UNDECIDED (δ = e(hgn) = −0.004 < 0.112 but halves +0.148/+0.674
straddle).  Worlds: twin identical protocol; EPSTEIN different
band structure (10 vs 13 bins) with live orbit; SCRAMBLE orbit
DECAYS (|O| = 0) as predicted.  Named successor: the ONE-SIDED
edge companion construction (the boundary case as the main case)
or the martingale second route; the r324 MEASURED composition
stays the honest end state.  EXT3 anchors not adopted (r329
record uncommitted at spec time, disclosed).  Experiments-side,
NO ledger row, NO RH CLAIM.

**The L\* capacity fork (r334, reviewer rank 2 — the new L\*
language).**  `fold_capacity_probe.py` (26/26, SPEC_SHA
`d5698a1d2d9c96c3` final with record, two-commit protocol:
pre-freeze `67833ecf`, record `9bf3a0aa`) reopened L\* under the
Fold-Capacity-Packing reading: the polynomial capacity
`Cap_{μ,N}(E) = inf{∫p²dμ : deg p < N, |p| ≥ 1 on E}`, the exact
level-set formula, and the chain (K1) `ν(E) ≤ κ·Cap(E)` + (K2)
layer-cake `≤ C_lev·∫p²dμ` + (K3) `κ·C_lev < 1 ⟹ L*`.  The
machinery is certified (dual active-set solver with KKT
certificate on all 138 723 solves, exact-Fraction enumeration
route to 6.8e-16, point capacity == inverse Christoffel against
the r284 route to 9.7e-16, mp dps-40 wards).  **Verdict
`ALLSET_NEEDED` + `WORLD_SEEN` + `CAPACITY_CONTENT`** (the
restatement rule does not fire, on two clauses): the small exact
all-subset census (4 instances, S ≤ 16, `Cap_abs` by full sign
enumeration) shows the supremum is NOT uniformly
interval-carried (I4: an outermost symmetric two-component pair
wins, ratio 0.7013) and not uniformly the λ shadow (κ_all/λ down
to 0.8504); on w9 the interval clause is carried by the r284
SHALLOW-EDGE PAIR (folds 2/4) at κ_int = 0.999567 — a new
source-pure near-wall coordinate whose margin 1 − κ_int (2.6×
the spectral margin on w9) shrinks along the ladder at the same
locus; the K3 products 0.90..1.06 straddle 1 across the live
rungs (no reserve), and the gate-side structural ceiling is
explicit: the extremal direction's own chain gives
κ·C_lev = 1.0085 ≥ λ_max, so NO (K1, K2) pair valid for all
polynomials can put the product below λ_max — the capacity route
cannot out-margin the spectrum.  Honest content: the interval
clause is a sharp world DISCRIMINATOR (EPST κ_int 1794, SCR
8.5e6 — both dead worlds shout κ > 1 through intervals alone,
which certifies λ > 1 exactly; MAIN/TWIN < 1 everywhere,
twin-stable to 5.5e-9).  Experiments-side, NO ledger row, NO RH
CLAIM.

**The one-sided edge packing (r335, the edge form after the r333
structural finding).**  `edge_packing_dichotomy_probe.py` (36/36,
SPEC_SHA `950c9b9ee96298e0` final with record, freeze
`a1f4c0dd97edcbf4`, two-commit protocol: pre-freeze `ec0fd5d5`,
record `3fa88d02`) tested whether the one-sided edge-companion
construction carries: for every dominant β/ω pair a ONE-SIDED
orbit O_e (same key, window from ps toward the interior, reach =
the masked support radius W/2, parameter-free a-priori) plus the
margin mass census, combined into the exact dichotomy
`q_max ≤ min(Q_e, M_b)` with `Q_e = ng·κ/(|O_e|·rmin_e)` and
`M_b = ng·κ·S_marg/L1` — EITHER many companions OR small edge
mass; EXT3 = the 12 r329 record anchors adopted as pure test rows
(committed `8cbd95f9`).  **Verdict `EDGE_ORBIT_ALSO_SHORT`** — no
clause certifies (T1-edge viol 43/40/38 of 51 at a = 0/1/2, named
0/4; dichotomy C_D 141.84/33.06/7.71 viol 14/11/6) — but with two
structural findings: (1) **the spikes flip to COVERED** — the MASS
arm covers all four named rungs at every a (kz53 B^D_1 12.51 /
kz83 20.47 / kz67 2.62 / kz55 11.86 vs C_D(1) 33.06, reserves
1.6×–12.6×): the r333 spike blocker is resolved pointwise at the
named rungs by the edge-mass clause; what violates now is the
mid-band family (kz73/76/61 + deep kz95/98/109, margin share
0.18–0.28, both arms 6–10); (2) **the inward window misses the
edge companions** — at the spikes |O_e| = 1 (g\* alone; 6
inner-empty rungs gained 0): the r333 one-sided same-key
companions sit OUTWARD of g\* (between the support edge and ps),
so T1 in the inward one-sided form fails harder than the
two-sided T1 and the spike coverage comes entirely from the mass
arm.  T2/T3 stay theorem grade on the edge orbit (rmin_e worst
0.3162 > 1/4, cband_e worst 3.1630 < 4 on 81 live worlds;
membership 1 ≤ 2); both arms and the arm partition are exact live
(arm census O 15 / M 62 of 77 — the mass arm is the rule); δ =
e(hgn) = −0.004 with halves +0.148/+0.674 straddling (exact r333
continuity).  Worlds: twin identical protocol; EPSTEIN operates;
SCRAMBLE |O_e| = 2 vs inner 0 — the one-sided orbit does NOT
decay on SCRAMBLE (honest negative: the support-radius reach is
permissive).  Named successor: the packing language is exhausted
at the edge — the reviewer's martingale second route (named only,
not executed); the r324 MEASURED composition stays the honest end
state.  Experiments-side, NO ledger row, NO RH CLAIM.

**The Chebyshev T+H parity section (r336, reviewer rank 3 — the
last pre-authorized L\* arm).**  `lstar_parity_section_probe.py`
(25/25, SPEC_SHA `9e4bf4f129405a64` final with record, two-commit
protocol: pre-freeze `d38b43c1`) put the parity-section language to
its three kill tests.  The exact coordinate stands: with the signed
cosine moments `c_r = Σ(μ_j − ν_j)cos(rθ_j)` on the affinely
mapped hull, `∫p²d(μ−ν) = Σ a_k a_l (c_{|k−l|} + c_{k+l})/2`
EXACTLY (bit-equal in Fractions on the rational toy; 2.2e-15 on
w9; ≤ 5.9e-15 on the ladder sample) — **the L\* matrix is a finite
Toeplitz+Hankel parity section**, congruent by the exact rational
congruence `M = S·V₊ᵀT_M(c)·V₊·S` to the EVEN parity compression
of the full symmetric Toeplitz section (the mirror of the
v549/v550 odd-sector step-over machinery); section flip at f64
signs clean (+5.7e-8 at 184 / −1.2e-8 at 185, negative inertia
exactly 1).  **Verdict `PARITY_FINITESECTION_CARRIER` +
`NEGATIVE_WINDOW_STEP_OVER(ii)` + `PARITY_WORLD_BLIND` +
`PREDICTOR_ORDERING`** — kill test 1 CARRIES (the symbol attains
−5.07 in 127 negative zones, total measure 0.2766 π, while the
184-section is PD: the positivity is a genuine finite-section
effect, and the Fejér mean is positive on all 183 active parity
modes of MAIN and TWIN while 34 Dirichlet readings sit in Gibbs
zones — clause (ii) steps over, exactly the T151 template); kill
test 2 FAILS (all four sealed statistics world-blind under the
r281 rule; the near-separation — Fejér floor positive on
MAIN/TWIN, negative on EPST/SCR/HL2 — is broken by SMOOTH, a dead
world with positive Fejér floor: its death is not a
local-negative-mass event); kill test 3 PARTIAL (the blind
Fejér predictor lands in the 1-octave band on 42/42 rungs, max
dev 0.224, spearman +0.999, but only 1/4 controls — no law).
The honest central negative: the w9 Fejér floor does NOT dip at
the crossing (+2.16e-3 at 184, +1.99e-3 at 185) — the L\* flip is
not a local-mass event at Fejér resolution; it lives in the
sharper two-atom extremal of r284.  Banked: the Fejér floor as a
SECOND source-pure near-wall coordinate (positive with empty
Schur block on 42/42 rungs, min +4.44e-4; falls ~5× along the
ladder where 1 − κ_int falls ~50× — not the same number in
disguise).  With kill test 2 failed and kill test 3 partial, the
reviewer's own rule reads: the parity-section lane does NOT stay
open as an L\* proof arm; the exact T+H coordinate (C1) is
banked.  Experiments-side, NO ledger row, NO RH CLAIM.

**The fold-martingale second route (r337, the reviewer's
pre-named terminal route after the exhaustion of the packing
language).**  `fold_martingale_probe.py` (37/37, SPEC_SHA
`617a83a6b3105acb` final with record, freeze `5460c2fa54983956`,
two-commit protocol: pre-freeze `30ac0cb5`, record `d0d272e3`)
read the signed fold-group masses of the argmax block as an
adapted process along the SEALED dyadic source-scale filtration
(adjudicated a-priori before the freeze: parameter-free /
theorem-adjacent / canonical coarse-to-fine; a protocol mutant
proves the audit bites), with three exact wards live — the
telescope `X_K == x_{j*}`, Cauchy–Schwarz `zr ≤ √n_act`, and the
identity `q_max = zr·κ·√V/L1`.  **Verdict
`MARTINGALE_STRUCTURE_ONLY`** — the structure is there: the
Z-maximal certification carries at the Azuma form `zr ≤
C_Z(1)·√(log m)` with C_Z(1) = 0.7364 and **0 test violations of
51 (incl. all 12 EXT3), named 4/4** — but the derived maximal
arm `A_z = C_Z·√(log m)·κ·√V/L1` fails pointwise at every a
(viol 18/13/11), so no new m₀\*.  Three structural findings:
(1) the certification is tight but nearly structureless — c\*
med 0.949 (ONE increment dominates the quadratic variation: the
r324 single-heavy-scale finding recurs in signed form), zr med
1.1012 essentially flat, contraction median exactly 0.500;
(2) **the r335 mid-band blocker dissolves** under the martingale
arm — head-to-head at a = 1: kz73 B^D_1 54.40 → B^C_1 0.399,
kz76 53.59 → 0.483, kz61 48.91 → 0.755, kz95 41.96 → 0.352,
kz98 41.68 → 0.400, kz109 40.90 → 0.299 vs C_C(1) 0.5293
(coverage 5/6, improvement 70×–180×; pointwise q_max ≤ A_z on
77/77 rungs); (3) the violation margin **migrates to the deep
anchors and the spikes** — A_z undercuts the mass arm everywhere
(letters A 77 / M 0), the freeze collapses 33.06 → 0.5293, and
the new a = 1 violators are EXT3 kz51/62/54/123/125 (up to
3.95×) plus named kz53/kz67/kz83 (≤ 2×): the r335 mass-arm spike
coverage is lost inside the min (a two-constant certification,
each arm against its own freeze, was not this round's sealed
contract — named residual direction only).  Measured route
growth e(B^C_1) = +0.059 vs the r324 chain's +0.172 (3× flatter,
but the pointwise certification fails).  Worlds: SCRAMBLE zr
1.1159 / contr 0.500 vs ladder med 1.1012 / 0.500 —
DISTINGUISHED by a 1.3 % zr margin only (honest negative: the zr
statistic is essentially world-blind in magnitude; the
arithmetic content lives in the calibrated constants).  The r324
MEASURED composition stays the honest end state.
Experiments-side, NO ledger row, NO RH CLAIM.

**The density martingale and the moment dictionary (r339, the
reviewer's reformulation of the terminal — transport of mass
density through the fold genealogy instead of a maximum).**
`fold_density_dictionary_probe.py` (36/36, SPEC_SHA
`822d6bff564d88ba` final with record, freeze `a845041cdf866543`,
two-commit protocol: pre-freeze `907ebfd1`, record `7863495f`;
INDEPENDENT of r337 by the reviewer's verdict-reading rule — r337
tested the RAW signed process, the density per remaining
descendant is a martingale by algebra) built the canonical fold
genealogy over the m final atoms `a_i = |PΔ_i|` (Leg 0 sealed
BEFORE the freeze: the iterated r270 adjacent pairing —
parameter-free, construction-canonical, position-contiguous;
balanced bisection as sensitivity; a protocol mutant proves the
audit bites).  **The exact layer stands at theorem grade**:
`X_k = d(V_k)/d(root)` is a nonnegative martingale from mass
conservation alone (Fractions BIT-EQUAL on w9+w13, 76 nodes;
per-node identity 2.2e-16 on 81 live worlds; cubic recursion and
envelope exact; Jensen drift ≤ 0 exact), and the moment
dictionary `E[X_∞²] == m·M₂`, `E[X_∞³] == m²·M₃` is exact
against the r324 chain — the terminal target **is**
`E[X_∞³] ≤ C(log m)^A`, and the banked r306 C₂ = 1.069 (0/57,
re-gated) already certifies A = 2 on the 57-rung set.  **Verdict
`LOCAL_INFLATION_SUPERCRITICAL`** — neither Γ budget certifies:
the full worst-case budget `W_F = Π_k max_v Γ(v)` violates
44/43/40 of 51 at a = 1/2/3 (C_F 39.91/8.98/2.02), the good
tree (heavy jumps max R_c > 3/2 removed) violates 33/20/11;
e(W_F) = +0.956 DECIDED supercritical.  **The gap is the
finding**: the true expectation e(m²M₃) = +0.112 grows 8×
slower than the worst-case budget — Γ_max is dominated by
near-leaf degenerate pairs (per-level med profile 1.05 → 3.99
against the balanced-pair ceiling 4; the deep half carries med
0.712 of log W_F) where hardly any path mass sits: the concrete
R341 instruction is that the Bellman argument must weight Γ by
path probability, not take per-level maxima.  The Γ_max census
is world-blind (SCRAMBLE 51.76 below ladder med 223.86 — honest
negative) but the dictionary itself separates (SCRAMBLE m²M₃
20.59 vs ladder med 6.22); the worst W_F violators are the EXT3
deep anchors (kz42 at 18.29×) — the r337 deep-anchor margin
recurs; the drift is strictly ≤ 0 but NOT period-4 (r4c 0.957
vs bar 0.5 — the clause-5 coboundary subtraction does not
apply).  kz53 is the one spike with its Γ argmax mid-tree
(level 2, 4.303 — the r324 single-heavy-scale event as an early
heavy fold); every mid-band rung inflates at the last fold
level.  The r324 MEASURED composition stays the honest end
state.  Experiments-side, NO ledger row, NO RH CLAIM.

**The weighted Cauchy–Binet transport (r340, the new L\* base
architecture after the capping of capacity and parity: transport
instead of norm bound).**  `cauchybinet_hall_probe.py` (26/26,
SPEC_SHA `86587b4895e95689` final with record, freeze
`82c4d1839386a853`, two-commit protocol: pre-freeze `3736af89`,
record in the same change) re-read L\* through the exact signed
Cauchy–Binet split `D_n = det H_n(μ̃) = Σ_E W − Σ_O W` with
`W(I) = (Π|w_i|)·Δ(x_I)²` (Sylvester dictionary gated: `D_n > 0
∀ n ≤ N_w` IS L\*), built the sealed canonical fold-exchange
graph (companions = nearest free μ atoms per side, step cap 2;
exact edge ratio `W(J)/W(I) = (|w_x|/|w_y|)·Π((x−z)/(y−z))²`
bit-equal on 192 small-model edges, w9 log ward 2.3e-14), and
measured the transport by exact integer max-flow plus FULL
`2^|O|−1` Hall-cut enumeration on exact small models.  **Verdict
`CUT_REDUCTION_FAIL` + `SOURCE_EXCHANGE_INSUFFICIENT` +
`TRANSPORT_WORLD_BLIND` + `NEW_COORDINATE`** — the round's
central structural hypothesis (minimal weighted Hall cuts =
Gale-order interval ideals, which would have re-read 1 − κ_int
as a dual Hall slack) is **refuted at the smallest exact scale**:
5 of 9 fully enumerated loci carry NON-ideal argmin cuts
(antichain sections and mid-order singletons pinned to the
shallow negative atoms; the principal-ideal family misses the
exact minimum by 1.17×–10.2×, growing toward half-filling).
Coverage itself is exact almost everywhere (TOY 2/2, F1 3/3,
F2 4/4 with strict reserve `D_n > 0`), but on the PD instance I2
the elementary edge set fails exactly at half-filling depth
`n = N = 5` (deficiency 0.023).  The extremal test fails
honestly: the W9EDGE binding cut carries only 3.86e-3 of its
W-mass on configurations containing both r284 extremal atoms
(folds 2/4) — the extremal pair is a spectral object, not the
Hall-cut carrier.  **The structural finding: the shallow-edge
prefix of EVERY world, MAIN included, is locally indefinite from
n = 2 on its own configuration space** (W9EDGE D signs
+,−,−,−,+,+ at a resolved −9.7e-3 relative deficit; controls
break at n = 1) — the full-window positivity is a GLOBAL rescue,
not a local edge property, the flow-side mirror of r336
`PARITY_WORLD_BLIND`; the ladder edge Hall slack is NEGATIVE on
42/42 rungs (−0.13…−0.18) against uniformly positive margins — a
new sign-opposite near-edge coordinate, not a restatement.  m4
documents the r282 Kasteleyn boundary exactly (the sign gauge
shifts `D_2` by exactly `2Σ_O = 1/3`).  Twin: identical cut form
at 4.8e-8.  No R341 GO: the interval-Hall route to a quotable
theorem does not open; banked are the exact Cauchy–Binet/Hall
machinery, the local-edge indefiniteness of every world, and the
sign-opposite edge coordinate.  Experiments-side, NO ledger row,
NO RH CLAIM.

**The two-atom extremal solved as its own problem (r342, the Ü1
main find of the r338 door-2 full revision).**
`pair_extremal_probe.py` (38/38, SPEC_SHA `b09f8ccdcfdedc16`
final with record, freeze `4ca125c7906149bf`, two-commit
protocol: pre-freeze `bdc0b439`, record `5c80e561`) solved the
r284 two-atom extremal as its own problem: at the binding point
L\* is EXACTLY the determinant condition `c² < (1−d₁)(1−d₂)`
with `d_k = v_k K_N(y_k, y_k) = E_kk`, `c = √(v₁v₂) K_N(y₁, y₂)
= E_12` at the shallow-edge pair (folds 2/4, below log 2 — the
pair on 69/69 windows incl. the 12 r329 EXT3 fresh anchors as
pure test rows, pair mass 0.9352…0.9994).  **Verdict
`PAIR_LAW_FOUND` + `PAIR_WORLD_COMPLETE`** — (i) the component
laws are halves-stable on the 57 (p −0.754 / q −0.645 / c
−0.697 vs the disclosed r338 priors −0.71/−0.66) and extend to
the EXT3 pure test (11/12 per column); (ii) the pair WEIGHTS
are closed-form source objects: the digamma dictionary
`F_A(ξ) = −log π + Re ψ(1/4 + iξ/2)` for the archimedean layer
plus the EXACT two-cell tent closed form of the prime layer
predict v at median rel dev < 1e-4 over all 69 rungs; (iii) the
α composition closes bookkeeping-exactly (α_pair 3.341 vs
α_full 3.332 vs composed 3.333; 42-only 3.059 == r286 3.06)
with the OPEN analytic remainder ρ_r = 2.624 — the cancellation
exponent of the determinant reserve `r_det = 1 − c²/pq`, the
concretized specialist question; (iv) the sealed world census
closes Ü2: the union `{PR ≥ 3} ∪ {κ_int ≥ 1}` is the first
measured world-complete criterion of the L\* lane — with the
honest structure that the PR clause alone FAILS (the dead
worlds are also concentrated at their own N_w, PR 1.80…2.92)
while the κ_int clause alone is world-complete: EPST 1794 /
SCR 8.51e6 (== r334) / **SMOOTH 2.193 / HL2 1964 (first
evaluation)** ≥ 1, live MAIN/TWIN 0.999567 < 1.  **The R343
material, measured:** the exact Schur dressing makes L\*
equivalent to `{λ_rest < 1}` AND `{dressed pair determinant}`;
the rest margin decays parallel to the full margin (slope
−3.28, offset ~20×) while the DRESSED pair reserve `r'_det` is
FLAT (slope +0.018, w9 value 0.303) — the O(1) certificate
candidate a full L\* theorem needs.  Honest negatives: `r_det`
is halves-curved (−0.767); kz56 breaks the R2 band
(`EXT3_FAMILY_BREAK`, the 12–40 % 2×2-uniformity does not
extend unmodified to the small-gap family); the kernel side
(why `v_k K_N → 1` with deficit ~N^−0.7 at the hard edge; r285
sub-classical 0.38 confirmed as the deep-rung late-window
growth) stays without closed form.  Experiments-side, NO ledger
row, NO L\* claim, NO RH CLAIM.

**The dressed pair reserve under fire (r343, the r342 O(1)
certificate candidate: analytic lower bound or refutation).**
`pair_coupling_probe.py` (38/38, SPEC_SHA `9ffc2705f326d625`
final with record, freeze `a24718fcb8c72f85`, two-commit
protocol: pre-freeze `e652071f`, record `95fcdcae`) put the
r342 Schur coordinates under fire.  **The exact spine:** the
resolvent-correlation identity `r'_det = 1 − M₁₂²/(M₁₁M₂₂)`
with M = the pair block of `(I−E)^{-1}` (gated per rung against
the independent full-inverse route + EXACT-Fractions wards on
the r334 rational instances, incl. the exact leading-minors
criterion for `λ_rest < 1`), and the exact spectral
decomposition of the dressing over the rest eigenmodes — the
contract's formula `r'_det = f(bare scalars, coupling vectors)`
measured per rung.  **Verdict `DRESSED_RESERVE_DECAYS(soft:
curvature clause)` — read precisely:** every HARD kill clause
misses (slope +0.018 dead flat; the 12 EXT3 rows — genuinely
BLIND, their dressed columns were never printed in r342 —
12/12 in the 0.5-decade band, incl. kz56 whose dressed reserve
0.1605 is flat while its bare R2 breaks the band: the dressing
absorbs exactly the Klein-gap family effect; 6 fresh EXT4
windows at N_w 2656…3181 — the deepest L\*-lane windows ever
measured, sealed source-pure selection under the family domain
cap z² ≤ 400000 — 6/6 in band, margins positive), but the
sealed FLAT clause fails on the halves-curvature bar (−0.426 >
0.35, a bar inherited from the decay-law convention applied to
a flat O(1) column with ~×2 scatter — sealed before sight,
applied honestly): the candidate is neither refuted on any
family nor certified flat.  **The structural findings:** (i)
the dressing eats bare reserve and coupling in near-exact
proportion (Δd₁/p and Δc/c → 1.000/−1.000; ρ' pinned ~0.846;
`m2'/margin` 1.0003…1.031 — the dressed 2×2 is a near-exact
λ_max carrier); (ii) **the bound autopsy closes Leg C
honestly**: the theorem-grade triangle/Cauchy–Schwarz/operator
chain certifies ZERO of 75 rows — the dressing works by the
SIGNED cancellation Δc ≈ −c (overshoot 73× at w9, to 2.4e11 at
depth): an O(1) bound must control the signed cancellation —
the concretized specialist question; (iii) `r'_det` is a gated
**TOP-2 spectral object** (median top-2 dev 0.063, K_res ≤ 3
on 75/75): its flat content IS the flat gap ratio
`g21 = (1−λ₂)/(1−λ₁)` (median 8.2; kz56 at 16.0 with the
reserve unmoved) — the R344 coordinate; (iv) `λ_rest` alone
separates the six instrumented worlds (dead ≥ 1, live 0.9963;
Schur equivalence 6/6 exact; slope −3.276 == r342 bit-near,
offset ~21×); (v) **`PEELING_STRUCTURE_FOUND`**: the rest block
is again a shallow-edge two-atom extremal (median rest-pair
mass 0.9902; level-2 reserve ≥ 0.23 on the 57), the w9 cascade
peels (2,4)→(6,8)→(10,12) with monotone margins — concentration
breaks at level 4: a level-2 recursion statement, not an
induction; (vi) RESTATEMENT does not fire (median shadow dev
1.0): the dressing is real content.  Honest negatives: the
sealed flatness protocol needs a curvature-aware successor for
O(1) columns (R344 protocol material, not retrofitted); the
EXT4 rows are disclosed-seen, only EXT3 carried blind weight.
Experiments-side, NO ledger row, NO L\* claim, NO RH CLAIM.

**The R\* mass balance and the two-scale freeze (r344, the named
r341 follow-up: each arm on its own sealed freeze).**
`fold_two_scale_balance_probe.py` (41/41, SPEC_SHA
`06375c3ade846178` final with record, freeze
`bf16912825178c44`, two-commit protocol: pre-freeze `d875647c`,
record `dfb8d1c1`; amendment a1 disclosed — two anchor-GATE
reading fixes only, no sealed surface touched) ran the stop
threshold R over the sealed grid (3/2, 8/5, 5/3, 12/7, 7/4, 2)
on the r341 exact path layer (imported verbatim; every exact
ward live at all six grid points, Fractions bit-equal at the
two sealed pins).  **Verdict `TWO_SCALE_PARTIAL`: the
three-arm min-coverage CERTIFIES at R̂ = 7/4 at a = 1 — 0
violations on the 51 test rows, named 4/4, midband 6/6 — but
the sealed choice formula (argmin max(e_H, e_G)) is not
halves-stable (the dyadic halves pick 8/5 vs 7/4): R\* remains
a tuning surface in the exponent language, said honestly.**
The findings: (i) **the mass-balance curve is measured** (four
interior grid points new): hsh med 0.872 → 0.266 → 0.000
across the grid, E3h share med 0.944 → 0.386 → 0.000, W_B med
1.489 → 3.796 → 6.096; at R = 2 the algebraic pair limit fires
as derived (no pair is ever heavy, 71/77 rungs lose the heavy
arm, the sealed eligibility guard fires); **the arm exponents
cross exactly once**, in the gap (12/7, 7/4) — the
equilibration point is real, its location not yet
halves-stable; (ii) **the disjunctive (complementarity)
coverage certifies at EVERY grid point already at a = 1**
(C_H on the hand-off column `(m·q_max)²·msh`, C_G on the
ε-chain column W_B, third arm = the banked r321 sliding bound
`1.3056·F_A²` at the frozen constant), while the r341
component-AND form still fails (13/10/10) — the r337 directive
is exactly what certifies; **the third arm is load-bearing**:
8 of 51 test rows are covered only by the r321 bound,
including the entire EXT3-B deep-anchor family 6/6 (kz51
F_A(ins) 5.54, reserve 5.3× — the r329 out-of-sample reading
reproduces under the insertion convention); (iii) **the
partition is NOT hsh-identifiable** (honest negative: the
sealed rule predicts the covering arm on 7/51 covered rows —
the partition is posthoc relative to hsh, a source-side
predictor remains open); (iv) kz55 is an EARLY-level ε-chain
anomaly (argmax level 1, the r339 kz53 class), NOT absorbed by
rebalancing (W_B 2.888 → 8.459) but covered comfortably by the
third arm (ratio 0.14): the single r341 Form-1 violator is a
sliding-bound row; (v) **the honest composition with a
certifying cover**: `E[X_∞³] ≤ [3.9859 + 1.1409](log m) +
1.3056·5.54²·(log m)²` on the covered set ⟹ m₀\* = 10^22.6 —
beats the r324 route 10^59.6 and the r341 envelope typing
10^24.0 (which had no cover), still above the r306 census
10^13.5; F_Amax = 5.54 is the kz51 insertion value, disclosed.
Worlds: the balance CURVE separates where the r341 point
budget was world-blind (SCRAMBLE stays heavy-flat, twin
collapses, EPSTEIN good-dominated).  Experiments-side, NO
ledger row, NO RH CLAIM.

**The path-weighted Bellman / Reverse-Hölder round (r341, the
terminal main round after the r339 dictionary — the two-arm
theorem candidate with path-probability weights instead of
per-level maxima).**  `fold_bellman_reverse_holder_probe.py`
(38/38, SPEC_SHA `7a9cf6b3fe24d0d7` final with record, freeze
`f0d2c744c46942bc`, two-commit protocol: pre-freeze `86f523e8`,
record `70d01ef8`) executed the r339 build instruction and the
r337 two-constant directive on the sealed r270 pairing tree.
**The exact path layer stands at theorem grade**: the
cubic-tilted tower `E[X_∞³] == E_Q[Π_k Γ(V_k)]` with steps
`p̃_c = p_c R_c³/Γ(v)` (the naive untilted form is FALSE — toy
break 7/18 EXACT, must-fail e1), the Φ–Γ identity
`Σ_c Φ(c)/Φ(v) == Γ(v)` with `Φ = A³/n²`, the stopped leaf
partition at τ = first heavy level with the exact hand-off
`E3h ≤ (m·q_max)²·msh` (the r335 q_max currency entering the
cubic layer exactly), and the good ε-chain envelope
`E3g ≤ W_B = Π(1+ε_k)` — all Fractions BIT-EQUAL on w9 + w13 +
three sealed toys.  **Verdict `PATHWEIGHT_ALSO_SUPERCRITICAL`
(per the sealed tree — the honest letter-reading: what fails is
the freeze-pointwise certification, NOT the growth: e(W_B) =
−0.214 and e(E3h) = +0.313 where r339 e(W_F) was +0.956)** —
the heavy arm fails its own freeze with NAMED 0/4 and the EXT3
deep anchors as worst violators (kz51 7.61×), Form 1 (the
integrated ε-chain) misses by the SINGLE violator kz55 (1.41×
at a = 2; EXT3 clean), Form 2 (Ψ-weight Bellman, prefactor 1)
certifies nowhere (best combo (A, C0) = (3, 1.0): 23 composed
node violations — the disclosed mid-tree risk is real).  **The
findings**: (1) the r339 thesis is CONFIRMED QUANTITATIVELY —
path weighting deflates the worst case 21× at the median
(E_P[ΠΓ] med 12.53 vs W_F med 265.54; path-weighted per-level
profile 1.05…1.80 vs Γ_max 1.05 → 3.99; near-ceiling nodes
carry pm3 med 0.161 of path mass); (2) **the stop at R* = 3/2
is too early on the path-mass scale** — the heavy arm inherits
med 94.4 % of the cubic mass with modal stop level 1 (the
mid-band stops its whole leaf mass at the FIRST fold level),
and the split is threshold-unstable (R_ALT = 7/4 moves the E3h
share to 0.386): R* is a tuning surface, not a canonical
constant — where between 3/2 and 7/4 the arms equilibrate is
the named open question; (3) the good arm is ONE RUNG from
certifying (kz55 only).  Worlds: the path-weighted budget is
also world-blind (honest negative; the dictionary stays the
only sharp separator).  Composition typed ENVELOPE ONLY (the
a = 3 envelopes would give the polylog m₀\* = 10^24.0 vs r324
10^59.6 and the r306 census reading 10^13.5) — **the honest
state of the route remains the r324 measured composition until
an arm certifies**.  Experiments-side, NO ledger row, NO RH
CLAIM.

**Lane status after wave 13.**  **The terminal: the MEASURED
composition is the end state of the route** (r324
`PILEUP_GROWS_SUBCRITICAL` +0.172 < 0.224 record / 0.188
chain-honest per the r328B audit, `Σq³ ≤
8.941·(log m)·m^{+0.172}/m²` ⟹ the N₂ density need for m ≥ m₀\*
(record 10^59.6 at the atom-level need, chain-honest 10^238) — the
gap from 274 to m₀\* the disclosed extrapolation
hypothesis; the certified pieces: the O(log m) scale count C_NSC
2.0258 at 0/39 and the O(log m) group count C_NG 2.6351 at 0/39;
r327 anatomy SOURCE: the heaviest group is exactly ONE β/ω fold
pair on all 65 rungs, kz53 = one bulk/window coincidence at
88.8 %; **the named direction is the coincidence geometry** —
WHERE can the single heavy coincidence sit, a source-geometry
question, not group count or multiplicity; r337 probed the
reviewer's martingale second route — structure present, arm too
weak; r339 then re-formulated the terminal as the density
martingale + moment dictionary: the exact layer is theorem-grade
and the target reads `E[X_∞³] ≤ C(log m)^A`, but the worst-case
Γ budget is supercritical — `LOCAL_INFLATION_SUPERCRITICAL` —
so the named R341 form is Bellman with path-probability
weights, not per-level maxima; r341 then executed exactly that
form — the exact path layer (tilted tower, stopped
decomposition, ε-chain envelope, Φ–Γ identity, heavy hand-off)
is theorem-grade and path weighting deflates the worst case
21×, but `PATHWEIGHT_ALSO_SUPERCRITICAL` per the sealed tree:
no arm certifies on its own freeze — the R* = 3/2 stop sends
94 % of the cubic mass to the heavy arm (threshold-unstable vs
7/4), the good ε-chain misses by the single violator kz55, and
the Ψ-Bellman with prefactor 1 is denied mid-tree; the named
material is the R* mass balance and the kz55/deep-anchor freeze
structure; r344 then measured exactly that balance — the
mass-balance curve over the sealed grid (3/2 … 2) with the
arm-exponent crossing located in (12/7, 7/4), and the r337
complementarity form executed with the path-weighted machinery:
the three-arm min-coverage (hand-off freeze + ε-chain freeze +
the banked r321 sliding bound as third arm) certifies with 0
violations on the 51 test rows at EVERY grid point at a = 1,
composing to m₀\* = 10^22.6 < 10^59.6 — the first certifying
cover of the lane — but `TWO_SCALE_PARTIAL`: the sealed
argmin-max stop rule is not halves-stable (8/5 vs 7/4), the
partition is not hsh-identifiable, and the third arm is
load-bearing for the deep anchors (kz51 F_Amax 5.54 enters the
uniform constant); the named rests are the canonical stop
criterion and a source-side partition predictor).  **The extraction:
the architecture stands** (r325 `ELEMENTWISE_STABILIZATION_GO` +
r326 `RH/Elementwise.lean`: `sourceExact_buildPrimeWindow` PROVED,
the comb stabilization elementwise PROVED,
`weil_nonneg_of_windowlocal` PROVED — extraction without the
ladder, H_cof replaced as the target route; **3 typed Level-C
rests**: the arch/pole kernel transcription (classical), the
bordered compression bridge (a named Prop), the source completion
(definitional)).  **L\*: reopened as the capacity fork** (the r323 fork was cleanly
aborted at the end of its design phase, before any write access —
nothing counted as a round; r334 reopened the lane under the
capacity language and measured its ceiling: `ALLSET_NEEDED` +
`CAPACITY_CONTENT` — the route cannot out-margin the spectrum,
its honest yields are the interval world-discriminator and the
new near-wall coordinate 1 − κ_int; the L\* contract itself
stays [O], nothing moved; r336 then executed the last
pre-authorized arm, the Chebyshev T+H parity section: the exact
parity-section coordinate is certified and banked, but kill test
2 fails — `PARITY_WORLD_BLIND` — and the blind predictor is
ordering-grade only, so the parity lane does not stay open as a
proof arm either; r340 then probed the reviewer's post-cap
transport architecture — the weighted Cauchy–Binet/Hall
transport, no norm bound — and the interval-Hall reduction is
refuted on exact small models (`CUT_REDUCTION_FAIL`), with the
structural finding that the shallow-edge prefix of every world,
MAIN included, is locally indefinite: the wall is a global
rescue, not a local edge property; r342 then solved the r338
two-atom extremal as its own problem — `PAIR_LAW_FOUND`: the
three binding-point scalars (d₁, d₂, c) carry stable laws incl.
the EXT3 pure test, their weights are closed-form source
objects (digamma dictionary + exact tent closed form), the α ≈
3.3 decay is their bookkeeping-exact composition with the open
remainder ρ_r = 2.624, and the union {PR ≥ 3} ∪ {κ_int ≥ 1} is
the first measured world-complete criterion (κ_int alone
suffices; SMOOTH 2.193 / HL2 1964 first evaluated); the Schur
dressing names the R343 coupling-control coordinates (rest
margin ∥ full margin, dressed reserve FLAT ~0.3); r343 then put
the dressed reserve under fire — no family refutes it (EXT3
blind 12/12 in band incl. the kz56 breaker, six fresh EXT4
windows to N_w 3181 in band), but the sealed flatness clause
fails on the halves-curvature bar, so the letter of record is
`DRESSED_RESERVE_DECAYS(soft)`; the triangle/operator bound
chain is measured DEAD on all 75 rows (the dressing is a signed
cancellation Δc ≈ −c — the typed specialist question), `r'_det`
is gated a TOP-2 spectral object (its flat content = the flat
gap ratio g21 ~ 8), λ_rest alone separates the six worlds, and
the rest block peels as a two-atom extremal again
(`PEELING_STRUCTURE_FOUND`, level-2 grade); the L\*
contract still [O]).  **Lean: 8 sorries, all typed** (the
five window-local legacy sorries byte-identical — the two true
WINDOW-LOCAL holes `lstar_subordination` +
`terminal_positive_main` unchanged — plus the three named Level-C
statements; census 5 → 8, the wave-12 census reservation partially
discharged).  The window-local positivity premise is proved
NOWHERE.  The mincut (base 4 / refined 5) is unchanged.  NO RH
CLAIM.

## Folder guide

```
rh/
├── README.md            — this file
├── INVENTORY.json       — machine-readable register of every RH-relevant
│                          file (path, SHA-256, round, ledger IDs, status);
│                          pinned entries drive the drift detector
├── lean/                — Lean 4 formalization pilot (r267 recommendation)
│   ├── README.md        — setup, build, proved-vs-sorry status
│   ├── lean-toolchain, lakefile.toml
│   └── RH/
│       ├── Basic.lean     — abstract chain data (h/F/c/rho/D/tau/tau^aug)
│       ├── Recursion.lean — PROVED core identities (recursion ⇒ drain,
│       │                    bilinear identity, telescope, terminal
│       │                    equivalence, two-branch trivial direction)
│       ├── Inertia.lean   — matrix-theorem layer (v956/v958); since the
│       │                    r305 Lean round PROVED except the Jacobi
│       │                    inertia theorem crossing_budget: psd_base,
│       │                    positive_prefix_firewall (Sylvester
│       │                    criterion, both directions),
│       │                    sylvester_pullback, half_filling_boundary
│       │                    + a general PosDef layer (bordered Schur
│       │                    step, injective principal submatrix)
│       ├── Window.lean    — r273: concrete VonMangoldtWindow structure,
│       │                    derived Hankel/chain/Schur objects, honest
│       │                    opaque MainWindow predicate; wave 5: T4
│       │                    proved; wave 6: the canonical form
│       │                    lstar_subordination (lemma L*, the ONE
│       │                    remaining sorry of this file) + L* ⇒
│       │                    free-window PROVED (free_window_positivity
│       │                    is a corollary of L* since r305)
│       ├── Source.lean    — r310 SOURCE INTERFACE (reviewer plan §6.3):
│       │                    explicit PrimeWindowSpec (prime powers,
│       │                    anchor, mesh level), buildPrimeWindow with
│       │                    Real.log nodes + mathlib Λ weights,
│       │                    MainWindowExplicit; PROVED: nodes_injective,
│       │                    predefined_family, mesh_refinement_compatible
│       │                    (+ mesh_refinement_shrinks),
│       │                    cofinal_prime_windows,
│       │                    finite_forms_converge_to_weil (prime side,
│       │                    stabilization form),
│       │                    rational_window_approximates; r310b
│       │                    (reviewer plan §8): four-stage support
│       │                    chain (primePow_index_injective →
│       │                    nodes_injective → foldedWindow
│       │                    aggregation → support_nodup) + the four
│       │                    source theorems (source_exact,
│       │                    weights_nonnegative, support_canonical,
│       │                    refinement_compatible — all PROVED);
│       │                    ONE documented sorry: the opacity bridge,
│       │                    since r310b in reviewer target form
│       │                    mainWindow_iff_builtFromPrimeSource
│       │                    (mainWindow_explicit_bridge proved from it)
│       ├── Closure.lean   — r305 RECONSTRUCTION: LStar + TerminalPositive
│       │                    predicates, lstar_terminal_implies_master
│       │                    (PROVED: L* + terminal q<1 ⇒ master
│       │                    positivity via Schur bordering + principal
│       │                    submatrices); THE MASTER THEOREM
│       │                    augmented_prefix_positive now a COROLLARY
│       │                    (no own sorry); the terminal hole named
│       │                    terminal_positive_main (sorry)
│       ├── Open.lean      — ladder bookkeeping + measured kill lists
│       │                    (r273: the two false-universal open theorems
│       │                    removed, no sorry anymore)
│       ├── PairBound.lean — PROVED finite pair algebra of the r269/r271
│       │                    fixed form c2PAIR (block decomposition, pair
│       │                    triangle, level-2 refinement); the H5 margin
│       │                    law retyped as pair_margin_main (sorry, r273)
│       └── Counterexamples.lean — r273 reviewer guards: machine-checked
│                            refutations of the three pre-r273 universal
│                            forms (proved, no sorry)
├── paper/
│   └── rh_program.tex   — the focused RH-program paper (English,
│                          self-contained, builds with pdflatex)
├── problem/
│   ├── lstar_problem.tex(+pdf) — the missing lemma L* as a fully
│   │                     standalone problem statement for external
│   │                     mathematicians (no project vocabulary)
│   └── verify_lstar_instance.py — machine check that the standalone
│                         definition IS the campaign object
└── verification/
    ├── make_inventory.py — regenerates INVENTORY.json
    └── run_rh.py         — the RH suite (see below)
```

## The standalone problem (`rh/problem/`)

The r283 attack reduced the open center to ONE scalar inequality, the
missing lemma L\*: for the window pair `(mu, nu)`, every real
polynomial `p != 0` with `deg p < (S+1)/2` satisfies
`int p^2 dnu < int p^2 dmu` (equivalently `lambda_max(E_{N_w}) < 1`).
`rh/problem/lstar_problem.pdf` states exactly this as a
**self-contained open problem** — the two measures defined from first
principles (von Mangoldt comb, tent sampling, archimedean lag
function, circulant spectral density, folding), the concrete
`S = 367` instance tabulated, equivalent Hankel / CD-kernel /
sampling formulations, the neutral numerical record, and the
moment-counting meaning of the degree bound — with **no project
vocabulary**, suitable for handing to external mathematicians
(approximation theory, orthogonal polynomials, sampling theory,
analytic number theory, operator theory).

`rh/problem/verify_lstar_instance.py` (numpy/mpmath only) rebuilds
the instance from the document's formulas alone and gates it against
the repository builders: atom positions/weights agree to `4e-16`
(relative), `lambda_max(E_184) = 0.99983248` and
`lambda_max(E_185) = 1.00003660` reproduced incl. a dps-60 mpmath
ward, the inequality checked on monomials and 500 random polynomials
— 12/12 gates, final line `LSTAR INSTANCE VERIFIED`.  The optional
`--ladder` sweep measures the margin `1 - lambda_max(E_{N_w})` on all
42 admissible windows: **all positive, but shrinking with S** — from
`1.68e-4` (S = 367, the *largest* margin of the family) down to
`1.42e-7` (kz64, S = 1717).  Honest note: the family margin loses
~3 orders of magnitude over the tested range, so the unbounded-family
version of the problem is genuinely uncertain in both directions.
The document claims nothing; L\* stays **[O]**, the open center —
since wave 6 registered in canonical form as its own ledger contract
**`PRIME.LSTAR.SUBORDINATION.01`** (established as the canonical
reduction by `v963_lstar_reduction_dictionary.py` [E]).
NO RH CLAIM.

## The RH suite

```bash
# from the repo root
python rh/verification/run_rh.py                 # full run (~8-9 min measured)
python rh/verification/run_rh.py --fast          # integrity + probes + lean (~36 s)
```

The suite runs, in order:

1. **Integrity** — SHA-256 of every pinned `INVENTORY.json` entry
   (drift in a pinned file = FAIL; unpinned living documents = INFO),
2. **Sealed probes** — the campaign probes r250–r327 from
   `experiments/tfpt-discovery/` in `--smoke` mode (fast, seconds each),
3. **The fifteen v9xx RH modules** — `v955`, `v956`, `v958`, `v959`,
   `v960`, `v961`, `v962`, `v963`, `v964`, `v965`, `v966`, `v967`,
   `v968`, `v969`, `v970` executed by path from
   `verification/` (never copied; `v959` ~3.5 min, `v968` ~50 s, the
   rest seconds — the wave-4/5/6/7/8/9/10/11/12/13 modules embed
   their probes in the sealed `--smoke` stage; skipped under
   `--fast`),
4. **Lean gate** — `lake build` in `rh/lean/` if a Lean toolchain is present.

Output is house-style: `[PASS/FAIL]` per item, final line
`RH SUITE: ALL CHECKS PASSED` or `RH SUITE: FAILURES n`.
It uses the repo venv `experiments/tfpt-discovery/.venv` if present.
Reference runs 2026-08-25: full **25/25 items PASS, wall 493 s, ALL CHECKS
PASSED** (pre-r268/r269 list); after the r268+r269 sync, fast mode
**23/23 items PASS, wall 37 s, ALL CHECKS PASSED**; after the r271 sync,
fast mode **24/24 items PASS, wall 36 s, ALL CHECKS PASSED** (Lean gate
included, 7 intentional `sorry` warnings); after the r270 sync (r270+r271
merged additively), fast mode **25/25 items PASS, ALL CHECKS PASSED**;
after the r272 sync, fast mode **26/26 items PASS, wall 43 s, ALL CHECKS
PASSED** (44/44 pinned entries byte-identical, Lean gate included); after
the r273 sync, fast mode **27/27 items PASS, ALL CHECKS PASSED** (45/45
pinned entries byte-identical); after the r273 Lean-pilot repair, fast
mode **27/27 items PASS, wall 43 s, ALL CHECKS PASSED** (Lean gate
included, **6 intentional `sorry` warnings** — down from 7: the two
Open-edge sorries consolidated into the one master-theorem sorry); after
the r274 sync, fast mode **28/28 items PASS, wall 45 s, ALL CHECKS
PASSED** (47/47 pinned entries byte-identical, Lean gate included,
6 intentional `sorry` warnings); after the r275+r276 sync (merged
additively), fast mode **30/30 items PASS, wall 47 s, ALL CHECKS
PASSED** (48/48 pinned entries byte-identical, Lean gate included,
6 intentional `sorry` warnings); after the r278 sync, fast mode
**31/31 items PASS, wall 50 s, ALL CHECKS PASSED** (49/49 pinned
entries byte-identical, Lean gate included, 6 intentional `sorry`
warnings); after the r277 sync (merged additively next to the
parallel r278 sync), fast mode **32/32 items PASS, wall 48 s, ALL
CHECKS PASSED** (50/50 pinned entries byte-identical, Lean gate
included, 6 intentional `sorry` warnings); after the wave-4 promotion
(v960/v961 pinned into the inventory, MODULES list extended), fast mode
**32/32 items PASS, wall 48 s, ALL CHECKS PASSED** (**52/52** pinned
entries byte-identical, Lean gate included, 6 intentional `sorry`
warnings); after the r279 sync (merged additively next to the parallel
wave-4 promotion), fast mode **33/33 items PASS, wall 49 s, ALL CHECKS
PASSED** (53/53 pinned entries byte-identical, Lean gate included,
6 intentional `sorry` warnings); after the r280 sync, fast mode
**34/34 items PASS, wall 49 s, ALL CHECKS PASSED** (54/54 pinned
entries byte-identical, Lean gate included, 6 intentional `sorry`
warnings); after the r281 sync, fast mode **35/35 items PASS, wall
50 s, ALL CHECKS PASSED** (55/55 pinned entries byte-identical, Lean
gate included, 6 intentional `sorry` warnings); after the r282 sync,
fast mode **36/36 items PASS, wall 53 s, ALL CHECKS PASSED** (56/56
pinned entries byte-identical, Lean gate included, 8 intentional
`sorry` warnings); after the r283 sync, fast mode **37/37 items PASS,
wall 49 s, ALL CHECKS PASSED** (57/57 pinned entries byte-identical,
Lean gate included, 8 intentional `sorry` warnings); after the wave-5
promotion (v962 pinned into the inventory, MODULES list extended to
seven, consumption notes on the r279–r281 probes), see the current
reference run recorded with the wave-5 exit gate; after the r284 sync,
fast mode **38/38 items PASS, wall 53 s, ALL CHECKS PASSED** (59/59
pinned entries byte-identical, Lean gate included, 8 intentional
`sorry` warnings); after the r285 sync (merged additively next to the
parallel `rh/problem/` registration), fast mode **39/39 items PASS,
wall 51 s, ALL CHECKS PASSED** (62/62 pinned entries byte-identical,
Lean gate included, 8 intentional `sorry` warnings); after the r287
sync (the L2 lane, merged additively next to the parallel r286 L\*
round), fast mode **40/40 items PASS, wall 49 s, ALL CHECKS PASSED**
(63/63 pinned entries byte-identical, Lean gate included, 8
intentional `sorry` warnings); after the wave-6 promotion (v963
pinned into the inventory, MODULES list extended to eight,
consumption notes on the r282–r285 probes, the r285
`make_inventory.py` sync repair included, merged additively next to
the parallel r287 sync), **full mode 48/48 items PASS, wall 519 s,
ALL CHECKS PASSED** (64/64 pinned entries byte-identical, all eight
modules green incl. v963 at 1.1 s, Lean gate included, **9
intentional `sorry` warnings** — the wave-6 canonical form
`lstar_subordination` added, the direction L\* ⇒ free-window
positivity proved); after the r286 sync (the L\* margin-scaling
round, merged additively), fast mode **41/41 items PASS, wall 48 s,
ALL CHECKS PASSED** (65/65 pinned entries byte-identical, Lean gate
included, 9 intentional `sorry` warnings); after the r289 sync (the
arch-kernel diophantine relevance test, merged additively next to the
parallel r288 sync), fast mode **43/43 items PASS, wall 52 s, ALL
CHECKS PASSED** (67/67 pinned entries byte-identical, Lean gate
included, 9 intentional `sorry` warnings); after the r290 sync (the
profile-space basin geometry round, merged additively next to the
parallel wave-7 `lstar_problem.tex` update whose refreshed pin the
parallel sync registered), fast mode **44/44 items PASS, wall 57 s,
ALL CHECKS PASSED** (69/69 pinned entries byte-identical, Lean gate
included, 9 intentional `sorry` warnings); after the r291 sync (the
ridge-anatomy round, merged additively), fast mode **45/45 items
PASS, wall 57 s, ALL CHECKS PASSED** (70/70 pinned entries
byte-identical, Lean gate included, 9 intentional `sorry`
warnings); after the r292 sync (the curvature two-form round,
merged additively), fast mode **46/46 items PASS, wall 51 s, ALL
CHECKS PASSED** (71/71 pinned entries byte-identical, Lean gate
included, 9 intentional `sorry` warnings); after the r293 sync
(the metric-reconciliation round, merged additively), fast mode
**47/47 items PASS, wall 53 s, ALL CHECKS PASSED** (72/72 pinned
entries byte-identical, Lean gate included, 9 intentional `sorry`
warnings); after the r294 sync (the F10 stability round, merged
additively), fast mode **48/48 items PASS, wall 53 s, ALL CHECKS
PASSED** (73/73 pinned entries byte-identical, Lean gate included,
9 intentional `sorry` warnings); after the r295 sync (the F10
sp-hardening round, merged additively), fast mode **49/49 items
PASS, ALL CHECKS PASSED** (74/74 pinned entries byte-identical,
Lean gate included, 9 intentional `sorry` warnings); after the
wave-9 promotion (v966 pinned into the inventory, MODULES list
extended to eleven, consumption notes on the r296–r300 probes —
strictly additive next to the in-flight r301 worker), **full mode
65/65 items PASS, wall 489 s, ALL CHECKS PASSED** (81/81 pinned
entries byte-identical, all eleven modules green incl. v966 at
1.4 s, Lean gate included, 9 intentional `sorry` warnings); after
the r302 sync (the unif-target round, merged strictly additively
next to the landed wave-9 state — all pre-existing inventory
entries byte-identical), fast mode **56/56 items PASS, wall 59 s,
ALL CHECKS PASSED** (83/83 pinned entries byte-identical, Lean
gate included, 9 intentional `sorry` warnings); after the r303 sync
(the atom-target regress-audit round, merged strictly additively —
all 83 pre-existing pinned entries byte-identical, the two unpinned
living-document rows refreshed their informational SHAs as
expected), fast mode **57/57 items PASS, wall 62 s, ALL CHECKS
PASSED** (84/84 pinned entries byte-identical, Lean gate included,
9 intentional `sorry` warnings); after the r304 sync (the
short-range-law round, merged strictly additively — all 84
pre-existing pinned entries byte-identical), fast mode **58/58
items PASS, wall 62 s, ALL CHECKS PASSED** (85/85 pinned entries
byte-identical, Lean gate included, 9 intentional `sorry`
warnings); after the wave-10 promotion (v967 pinned into the
inventory, MODULES list extended to twelve, consumption notes on
the r301–r304 probes — strictly additive, all 85 pre-existing
pinned entries byte-identical), **full mode 70/70 items PASS,
wall 515 s, ALL CHECKS PASSED** (86/86 pinned entries
byte-identical, all twelve modules green incl. v967 at 1.3 s,
Lean gate included); after the r305 Lean reconstruction round
(RH/Closure.lean added, Inertia layer proved — Lean-only, no
pinned entry touched; run coexisting with the parallel wave-10
consolidation), fast mode **59/59 items PASS, wall 62 s,
ALL CHECKS PASSED** (Lean gate included, **4 intentional `sorry`
warnings** — down from 9: the master theorem, the fog-free hole
and four of the five Inertia statements are now PROVED; see the
Lean status below); after the r310 source-interface round
(RH/Source.lean added — Lean-only, no pinned entry touched; run
coexisting with the parallel r309 probe round), fast mode
**61/61 items PASS, wall 77 s, ALL CHECKS PASSED** (Lean gate
included, **5 intentional `sorry` warnings** — up from 4 by exactly
the documented opacity bridge `mainWindow_explicit_bridge`; the
four pre-existing holes are byte-identical; see the r310 Lean
status paragraph below); after the r310b source-theorem refinement
round (reviewer plan §8 adjudication — RH/Source.lean restructured
in place, Lean-only, no pinned entry touched), fast mode
**62/62 items PASS, wall 69 s, ALL CHECKS PASSED** (Lean gate
included, **5 intentional `sorry` warnings** — count unchanged: the
one Source.lean sorry MOVED into the reviewer target form
`mainWindow_iff_builtFromPrimeSource`, and the r310 form
`mainWindow_explicit_bridge` is now PROVED from it; the four
pre-existing holes are byte-identical; see the r310b Lean status
paragraph below); after the r313 sync (the Rényi-3 proof-form
fork round, merged strictly additively — all pre-existing pinned
entries byte-identical, run coexisting with the parallel
r309/r311 rounds), fast mode **63/63 items PASS, wall 68 s, ALL
CHECKS PASSED** (95 inventory entries with the r313 probe pinned,
Lean gate included, 5 intentional `sorry` warnings — unchanged);
after the r311 sync (the block-Green nontriviality round, merged
strictly additively — all pre-existing pinned entries
byte-identical, run coexisting with the parallel r309/r313
rounds), fast mode **64/64 items PASS, wall 77 s, ALL CHECKS
PASSED** (96 inventory entries with the r311 probe pinned, Lean
gate included, 5 intentional `sorry` warnings — unchanged);
after the r314 sync (the signed-cubic-flux round part 1, merged
strictly additively — all pre-existing pinned entries
byte-identical, run coexisting with the parallel r311 round),
fast mode **65/65 items PASS, wall 75 s, ALL CHECKS PASSED**
(97 inventory entries with the r314 probe pinned, Lean gate
included, 5 intentional `sorry` warnings — unchanged);
after the r315 sync (the Φ₃-functional round part 2, merged
strictly additively — all pre-existing pinned entries
byte-identical, run coexisting with the parallel r312 round),
fast mode **66/66 items PASS, wall 75 s, ALL CHECKS PASSED**
(98 inventory entries with the r315 probe pinned, Lean gate
included, 5 intentional `sorry` warnings — unchanged);
after the r316 sync (the two-regime-bound round part 3, merged
strictly additively — all pre-existing pinned entries
byte-identical, run coexisting with the parallel r312 round),
fast mode **67/67 items PASS, wall 76 s, ALL CHECKS PASSED**
(99 inventory entries with the r316 probe pinned, Lean gate
included, 5 intentional `sorry` warnings — unchanged);
after the r312 sync (the block-Green membership round, merged
strictly additively — all pre-existing pinned entries
byte-identical, run coexisting with the parallel r314/r315/r316
fiber rounds), fast mode **68/68 items PASS, wall 109 s, ALL
CHECKS PASSED** (100 inventory entries with the r312 probe
pinned, Lean gate included, 5 intentional `sorry` warnings —
unchanged); after the wave-11 promotion (v968 pinned into the
inventory via head-sync, MODULES list extended to thirteen,
consumption notes on the r306–r316 probes — strictly additive,
all 100 pre-existing entries byte-identical), **full mode 81/81
items PASS, wall 579 s, ALL CHECKS PASSED** (97/97 pinned entries
byte-identical, all thirteen modules green incl. v968 at 49.6 s,
Lean gate included, 5 intentional `sorry` warnings — unchanged);
after the r318 sync (the indefinite-fork round, merged strictly
additively — all pre-existing pinned entries byte-identical, run
coexisting with the parallel r317/r319 rounds), fast mode
**70/70 items PASS, wall 115 s, ALL CHECKS PASSED** (103
inventory entries with the r318 probe pinned, Lean gate included,
5 intentional `sorry` warnings — unchanged); after the r321 sync
(the continuous-coordinate round, merged strictly additively —
all pre-existing pinned entries byte-identical, run coexisting
with the parallel r318/r320 rounds), fast mode **71/71 items
PASS, ALL CHECKS PASSED** (104 inventory entries with the r321
probe pinned, Lean gate included, 5 intentional `sorry` warnings
— unchanged); after the r322 sync (the antiphase-sign-law dig
round, merged strictly additively — all pre-existing pinned
entries byte-identical, run coexisting with the parallel
r320/r321 rounds), fast mode **72/72 items PASS, wall 110 s, ALL
CHECKS PASSED** (105 inventory entries with the r322 probe
pinned, Lean gate included, 5 intentional `sorry` warnings —
unchanged); after the r324 sync (the QMAX × M₂ origin round plus
its superseded r324-pre pre-work probe, merged strictly
additively — all pre-existing pinned entries byte-identical, run
coexisting with the parallel wave-12 worker and R325, whose v969
module and r325 probe land in the same window), fast mode
**75/75 items PASS, wall 120 s, ALL CHECKS PASSED** (109
inventory entries with the r324-pre and r324 probes pinned, Lean
gate included, 5 intentional `sorry` warnings — unchanged); after
the r327 sync (the fold-group mass-cap round, merged strictly
additively — all pre-existing pinned entries byte-identical, run
coexisting with the parallel r326 Lean round), fast mode **76/76
items PASS, wall 125 s, ALL CHECKS PASSED** (110 inventory
entries with the r327 probe pinned, Lean gate included; the Lean
gate measured 8 intentional `sorry` warnings in this window — the
parallel r326 Lean round in flight, not a change of this sync).

## Lean status

See `rh/lean/README.md` for the authoritative list. Summary: the Lean 4
project builds green (since the r305 Lean reconstruction round
**4 intentional `sorry`s**, down from 9; since the r310 source-interface
round **5** — the fifth is the documented opacity bridge in
`RH/Source.lean`, since the r310b refinement in the reviewer target
form `mainWindow_iff_builtFromPrimeSource` (the r310 form
`mainWindow_explicit_bridge` is proved from it), type
definitional/technical; **since r320 still 5, two of them RETYPED** —
the R319 audit found the r310b statement types jointly inconsistent
(U1–U3) and r320 repaired them, see the r320 paragraphs below: the wave-6 canonical form
`lstar_subordination` [the base/wall hole], the terminal statement
`terminal_positive_main` [the border/fiber hole, new in
`RH/Closure.lean`], its pair-coordinate refinement `pair_margin_main`
[same hole, r263 dictionary not yet formalized], and the Jacobi inertia
theorem `crossing_budget` [a mathlib gap, exact-rationally certified in
v962 T2]); the **provable
core identities of stage 4 are actually proved** (no `sorry`, no mathlib —
abstract field/ordered-field algebra); since r271 the **finite pair
algebra of the fixed drive bound is proved too** (`RH/PairBound.lean`);
since r273 the pilot is **truth-capable**: the three pre-r273 universal
open statements (refutable over bare bookkeeping — machine-checked guards
in `RH/Counterexamples.lean`) are replaced by the concrete window
structure `VonMangoldtWindow` and ONE master theorem
`augmented_prefix_positive` (`RH/Window.lean`, augmented-matrix positivity
through the half-filling cap, conditioned on the honest opaque
`MainWindow` predicate). Both former open edges are PROVED from the master
theorem by finite matrix algebra (block extraction ⇒ Sylvester prefix
chain; Schur complement ⇒ terminal `q < 1`; the Schur identity
`D_n = B − u^T H_n^{-1} u` proved against mathlib). **Since wave 5** the
hole is additionally stated fog-free in base coordinates:
`free_window_positivity` (`∀ n < cap : 0 < h_n` on `MainWindow` — the
r279 b3 gap statement / v962 T4 reinstform; `master_implies_free_window`
PROVED shows it is the base half of the master `sorry`), the v962 T1
counting layer is **proved for real** (`moment_counting_free_pivots`,
`first_forced_pivot`), T4 is proved (`main_window_reduction`), T2 is
stated (`crossing_budget`, `sorry` with v962 reference), and the wave-5
guard `upper_pinning_not_universal` (the exact one-negative instance with
machine-checked Hankel minors) plus the proved escape `o1_pinning_escape`
extend the counterexample layer; the matrix-level theorems (PSD base,
half-filling, Sylvester pull-back) remain roadmap `sorry`s, each with
round, ledger ID and probe reference in the header. **Since wave 6** the
hole is additionally stated in its canonical two-measure form:
`lstar_subordination` (lemma L\* on `VonMangoldtWindow` — for every real
polynomial `p ≠ 0` with `deg p < cap`: `∫p²dν < ∫p²dμ`; the wave-6
`sorry`, ledger `PRIME.LSTAR.SUBORDINATION.01` [O]), with the direction
L\* ⇒ free-window positivity **proved for real**
(`lstar_implies_free_window` via `lstar_implies_hankel_posDef` and the
quadratic-form dictionary `hankel_quadform` — finite algebra, no
`sorry`) and the composed corollary `lstar_free_window_main`.

**Since the r305 Lean reconstruction round (reviewer plan §3 + §6.1)**
the pilot carries the RECONSTRUCTION THEOREM
(`RH/Closure.lean`, `lstar_terminal_implies_master`, PROVED):
`LStar w → TerminalPositive w → ∀ n ≤ cap, A_{w,n} ≻ 0` — L\* plus the
r258/r260 terminal statement (`B > 0 ∧ q_N < 1`, faithful to
TERMINAL_Q_LAW and the v959 budget telescope) imply the full master
positivity by finite matrix algebra (L\* ⇒ `H_cap ≻ 0`; terminal ⇒
`D_cap = B(1−q) > 0`; Schur bordering ⇒ `A_cap ≻ 0`; principal-submatrix
restriction ⇒ all `A_n ≻ 0` — the matrix form of the backward Riccati
drain `D_n = D_N + Σρ_k`). Consequently `augmented_prefix_positive` is
now a **corollary** (no own `sorry`, NOT a third arithmetic input) and
`free_window_positivity` is a **corollary of L\*** alone. The two true
**window-local** holes (Level B of the wave-12 three-level proof
graph) stand formally alone: `lstar_subordination` (`RH/Window.lean`,
base/wall) and the terminal statement — `terminal_positive_main`
(`RH/Closure.lean`, border/fiber), whose pair-coordinate form is
`pair_margin_main` (`RH/PairBound.lean`; the r263 dictionary
`Z² / m = q_N` connecting the two coordinates is measured on 42/42
rungs, not yet formalized, so both keep their own honest `sorry`).
In the same round the Inertia matrix layer was proved for real
(`psd_base`, `positive_prefix_firewall` — Sylvester's criterion built
from scratch against mathlib (`posDef_of_leading_det_pos`) —,
`sylvester_pullback`, `half_filling_boundary`); only the Jacobi inertia
theorem `crossing_budget` remains a `sorry` (mathlib v4.29.1 carries
neither Jacobi's minor-sign rule nor Sylvester's law of inertia in
matrix-counting form — a formalization gap, not a mathematical one; the
exact-rational v962 certificate stands). NO RH CLAIM.

**Since the r310 source-interface round (reviewer plan §6.3)** the
formerly fully opaque source boundary is an explicit proof interface
(`RH/Source.lean`, conservative route: `MainWindow` stays opaque and
byte-identical, all four existing holes untouched). The exact real
window source is now CONSTRUCTED: `PrimeWindowSpec` (strictly
increasing prime-power atoms, anchor `p^k`, mesh level, arch/border
fields) with the node positions and comb weights DERIVED — nodes
`= Real.log(p^k)`, weights `= Λ(p^k)` via mathlib's exact
`ArithmeticFunction.vonMangoldt` — faithful to the Python builder
chain (`v563_paper2_readouts.von_mangoldt_table`/`build_window` →
`port_integrable_kernel_probe.folded_measure` →
`principal_bessel_probe.window_pack`; the CONSTRUCTION is translated,
never a measured value). PROVED for real: `nodes_injective` (the
`log p^k` positions are pairwise distinct — log monotonicity on the
ordered prime powers), `node_pos`, `combWeight_pos`, and the four
reviewer structure theorems — `predefined_family` (the family atom set
is decided by arithmetic alone: `n` is an atom iff `IsPrimePow n ∧
n ≤ anchor²`; constructive), `mesh_refinement_compatible` (refinement
changes only the mesh level, never the atom data; `rfl`-level by
design, plus the real content `mesh_refinement_shrinks`),
`cofinal_prime_windows` (Euclid via `Nat.exists_infinite_primes`), and
`finite_forms_converge_to_weil` (PROVED in the strong stabilization
form: for compactly supported test data the finite comb forms EQUAL
the Weil prime-side tsum once the window covers the support — comb
channel only; the archimedean-kernel transcription is the documented
classical TODO, the spectral/zero side stays the open content of the
program). The rational certificate layer is connected by
`rational_window_approximates` (PROVED: every real prime window admits
rational certificate data within every positive bound — density of ℚ)
and the mesh-level predicate `RepresentsSpec`. The ONE new `sorry` is
the opacity bridge `mainWindow_explicit_bridge`
(`MainWindow w ↔ ∃ s, RepresentsSpec w s`; type
DEFINITIONAL/TECHNICAL, not arithmetic: nothing about an `opaque`
predicate is provable by design — the bridge IS the honest new
interface, and it becomes `Iff.rfl` if a future wave takes the
invasive route of redefining `MainWindow`). Census 4 → 5, with the
increase fully accounted: the new interface states what was not
formalizable before r310, and everything in it except the one
opacity-forced bridge is PROVED. NO RH CLAIM.

**Since the r310b source-theorem refinement round (reviewer plan §8,
accepted r310 + adjudicated refinements)** the interface carries the
honest folding layer and the four source theorems. (1) THE FOUR-STAGE
SUPPORT CHAIN — the reviewer's warning that `nodes_injective` must not
be overclaimed is heeded structurally: stage 1
`primePow_index_injective` (PROVED: `p^k = q^ℓ ⟹ p = q ∧ k = ℓ`,
unique factorization via `minFac`), stage 2 `nodes_injective`
(unchanged — valid for the UNFOLDED `log p^k` nodes, which is exactly
what `buildPrimeWindow` produces), stage 3 `foldedWindow` (NEW: the
explicit quotient/aggregation construction for the folded source —
equal folded nodes merged via `Finset.image`, weights ADDED over each
fiber, transcribing `folded_measure`/`np.add.at`; with the
mass-conservation structure theorem `foldedWindow_mass`, PROVED: the
aggregation is a quotient, not a projection), stage 4 `support_nodup`
(PROVED: distinctness returns AFTER the aggregation, on the merged
support). FOLDING STATUS, honest: `buildPrimeWindow` contains NO
hidden fold — it is the unfolded chain; the fold is the explicit
stage-3 step, stated for an arbitrary fold map `φ : ℝ → ℝ` of which
the corpus map (`min(j, L−j)` index fold, `cos(2π j/L)` nodes) is an
instance whose exact transcription belongs to the documented classical
arch/border TODO. (2) THE FOUR SOURCE THEOREMS, all PROVED:
`buildPrimeWindow_source_exact` (the built window carries
DEFINITIONALLY the `Λ`/`log` source data — `rfl` by derivation),
`buildPrimeWindow_weights_nonnegative` (all weight channels, comb and
arch, fold-stable), `buildPrimeWindow_support_canonical` (the full
stage-1/2/3/4 chain), `buildPrimeWindow_refinement_compatible` (the
built window is mesh-independent, also under any FIXED fold map; the
corpus fold map itself is mesh-dependent — the mesh enters ONLY
through the fold map, which is why the fold is a separate stage).
(3) THE BRIDGE IN TARGET FORM: the one Source.lean `sorry` moved into
`mainWindow_iff_builtFromPrimeSource`
(`MainWindow w ↔ ∃ s, RepresentsWindow w (buildPrimeWindow s) s.mesh`
— the reviewer's target semantics, via the new window-level predicate
`RepresentsWindow` with `representsSpec_iff` as the definitional
link), and `mainWindow_explicit_bridge` is now PROVED from it; same
type (DEFINITIONAL/TECHNICAL, opacity-forced), census unchanged at 5.
(4) APPROXIMATION WARNING (reviewer verbatim): a warning block at
`rational_window_approximates` records that the rational windows are
certificate objects, not the definition, and that an approximation
proof chain would need error ≪ the shrinking L* margin (1.68e-4 →
1.806e-8, v964) — NOT established; no statement of the file uses
approximation as a proof path; the intended route is the direct real
construction. NO RH CLAIM.

**Since the r320 bridge-repair round (the R319 red-team audit U1–U3,
all three kernel-reproduced before repair).** The r310b statement
TYPES of the source interface were jointly INCONSISTENT — not the
proofs (there are none; the sorries are the holes), the types: (U1)
the bridge + `terminal_positive_main` ⟹ `False` (`RepresentsWindow`
never bound `u`/`B`; a `B = −1` window represented the empty spec),
(U2) the bridge + `lstar_subordination` ⟹ `False` (mesh-level-0
tolerance `log anchor` admitted the TOTAL node collision on the
{2,3,4} spec; `p = X − 1` vanishes on the collided window at
`cap = 2`), (U3) `pair_margin_main` + any main window ⟹ `False`
(free `(Zloc, runs)`: `Zloc = |F| + 1`, `runs = []`). Repairs (r319
recommendations 1–5): `RepresentsSpec`/`RepresentsWindow` retyped
(exact u/B-fidelity; separation discipline — tolerance below half the
minimal node gap and below every comb weight; the spec carries the
transcribable `budget_pos`, the r243 positivity half); the free spec
channels `archWeight`/`border` — found by the r320 verification to
carry the same disease beyond the audit — are guarded by the NEW
opaque `SourceExact` (the r273 convention applied to the spec side;
elimination = the arch/border/fold transcription TODO); the old types
are conserved and machine-refuted as three PERMANENT sorry-free guards
(`RH/Counterexamples.lean`: `old_bridge_terminal_inconsistent`,
`old_bridge_lstar_inconsistent`, `old_pair_margin_forces_empty`, all
quantified over an arbitrary predicate `P` in r273-G3 style); the
`(Zloc, runs)` extraction is now a DEFINITION
(`RH/PairBound.lean`: `signRuns`/`terminalDrive`/`bulkRuns`/
`edgeLocal`, split inequality PROVED as `canonical_split`) and
`pair_margin_main` is stated for the canonical extraction only; a
satisfiability WITNESS (anchor 2, atoms {2,3,4}, mesh level 4, exact
rational certificate window; `witness_represents` sorry-free) proves
the retyped predicate nonempty — `∃ w, MainWindow w` itself stays
deliberately unprovable (two opaque predicates; exactly what blocks
the adversaries). The bridge docstring's false "becomes `Iff.rfl`"
promise is corrected. Census UNCHANGED at 5; `lstar_subordination`,
`terminal_positive_main`, `crossing_budget` byte-identical. Two
wave-12 reviewer notes are binding here: **the sorry census does not
fully measure the global distance** — the missing Level-C extraction
architecture appears in *no* `sorry`, because the correct theorem
statements do not yet exist (the census of five counts the
window-local layer only); and `SourceExact` is a good firewall but
must **not** stand as a free assumption in a final proof — the named
open Lean target is **`CanonicalPrimeWindow` + exactly one
construction theorem `sourceExact_buildPrimeWindow`** (load-bearing
for the global connection, not formalization polish; the extraction
repair fork R325 is in flight on exactly this, not consumed). NO RH
CLAIM.

**The mesh-vs-anchor cofinality seam (r320 documentation).** The word
"cofinal" names two different directions, and the seam is now mapped
(docstring block in `RH/Source.lean` after `cofinal_prime_windows`):
`cofinal_prime_windows` (PROVED, Euclid) is the ANCHOR direction —
for every `N` an anchor whose window carries ≥ `N` atoms; orthogonally
`mesh_refinement_shrinks` (PROVED) drives the mesh below every
positive bound at a fixed anchor. Hypothesis **(H_cof)** (the carrier
lane, `TfptCarrier/CofinalWeil.lean`, v849) demands something else:
PSD certificates along a PRE-FIXED mesh-refinement tower.
`cofinal_prime_windows` does NOT supply this — it produces anchors,
never PSD certificates, and its direction (atom count) is not the
(H_cof) direction (refinement along one pre-fixed tower). The named
OPEN Lean goal: identify the `specFamily` mesh-refinement tower with
the v749 canonical tower (`verification/v749_simpler_tower.py`, stages
`n_k = 2^{2k+1}`, the 2-adic nesting of the v755 Schur recursion) and
transport window positivity along it into the (H_cof) shape. Until
that identification exists, no theorem of `RH/Source.lean` feeds
(H_cof), and none claims to. NO RH CLAIM.

**Since the r326 elementwise-architecture round (the R325 repair set
in Lean).** The Level-C extraction layer exists as named statements
(`RH/Elementwise.lean`), implementing the R325 fork adjudication
(`extraction_order_probe.py`, sealed, primary
`ELEMENTWISE_STABILIZATION_GO`): **(i)** the canonical window
predicate `CanonicalPrimeWindow` (`∃ a m ha, w = buildPrimeWindow
(specFamily a m ha)`) with the PROVED construction theorem
`sourceExact_buildPrimeWindow` — the wave-12 reviewer target: the new
REAL predicate `SourceExactSpec` (node/comb derivation + atom-set
completeness + budget positivity; arch/border deliberately not bound)
is provably satisfied by every canonical family member, so
`SourceExact` is **eliminated as a free assumption from the
extraction route** (the opaque guard itself stays untouched in
`RH/Source.lean`; the relation is the new typed completion sorry —
every family member admits a source-exact completion, genuine
arch/border/budget data on the same atom set). **(ii)** the native
dense class `GridElement` is BUILT (dyadic step-function
autocorrelations, the v749 class: derived autocorrelation, derived
even piecewise-linear interpolant, PROVED compact support), and the
elementwise finite stabilization is PROVED for the comb channel in
the corpus gauge `2Λ(n)/√n` with the anchor onset PREDEFINED from
the element's support alone (`elementAnchor f = max(1,
⌈exp(steps·D0)⌉)`) and NO mesh quantifier — the elementwise form of
`finite_forms_converge_to_weil`; the arch/pole kernel channels enter
as OPAQUE reads (the classical TODO made visible) with their
stabilizations as typed sorrys (classical, S2 — R325 S1 measured
them exact), and the full-form statement
`elementwise_finite_stabilization` is proved from the three
channels. **(iii)** the extraction WITHOUT the ladder is PROVED:
`weil_nonneg_of_windowlocal` — window-local positivity of the
canonical family (typed honestly on the plain full form with the
grid-compatibility guard) implies Weil-form nonnegativity on every
grid element by ONE finite instantiation per element (Euclid anchor
+ the element's own native mesh); **no mesh-cofinal PSD tower and no
transport is consumed** — (H_cof) is REPLACED as the target route
(the seam block in `RH/Source.lean` carries the r326 update; the
carrier-lane documentation stays historically correct). The bordered
⟹ plain compression step is the NAMED Prop
`BorderedCompressionBridge` (parametrized, no truth commitment) with
the composed extraction `weil_nonneg_of_bordered` proved. Census
5 → **8**: the three new sorrys are statements that were NOT
formalizable before this round — the wave-12 reviewer reservation
("the Level-C distance appears in no `sorry`") is thereby partially
discharged; the five pre-existing sorrys are byte-identical, and the
window-local premise (Level B, the two true holes) is untouched and
unclaimed. NO RH CLAIM.

## Canonical sources (never duplicated here)

- `verification/v955_…​.py`, `v956_…​.py`, `v958_…​.py`, `v959_…​.py`,
  `v960_…​.py`, `v961_…​.py`, `v962_…​.py`, `v963_…​.py`,
  `v964_…​.py`, `v965_…​.py`, `v966_…​.py`, `v967_…​.py`,
  `v968_…​.py`, `v969_…​.py` — the fourteen
  certified RH theorem modules (docstrings carry the exact theorem
  stock)
- `verification/status_ledger.csv` — all `PRIME.*` rows: the formal claim
  state including fine types; **the ledger wins**
- `tfpt_prime_front.tex` — the program prose (root)
- `experiments/tfpt-discovery/*_probe.py` — the sealed probes (SHA-256 pinned
  in `INVENTORY.json`)
- `experiments/next.txt` — campaign notes (German, newest at bottom,
  r250–r273)

## TODO (later consolidation steps)

- [x] **Pipeline integration (audit gate, wave 4)**: `bash build.sh audit`
      now runs `rh/verification/run_rh.py --fast` as its own
      "RH workspace" section (must end `RH SUITE: ALL CHECKS PASSED`).
- [ ] **Pipeline integration (rest)**: register `rh/` in `docs_map` /
      `website_map` / manifest.
- [ ] **Lean full build-out**: prove the Inertia layer (needs mathlib —
      declared option, see `rh/lean/README.md`), replace the Open-edge
      `sorry`s only if/when the mathematics exists.
- [x] **Wave 4**: rounds r260–r278 promoted as `v960`/`v961`;
      `INVENTORY.json` and the suite MODULES list extended in the same
      change (probes stay sealed experiments-side, consumed byte-exact in
      the sealed smoke stage).
- [x] **Wave 5**: rounds r279–r281 frozen as the small mathematical
      theory `v962` (four theorems + four refutations); `INVENTORY.json`,
      `make_inventory.py` and the suite MODULES list extended in the same
      change; the Lean layer carries the fog-free central hole
      `free_window_positivity` (T1/T4 proved, T2 stated, N1 guard).
- [x] **Wave 6**: rounds r282–r285 + the standalone problem frozen as
      the L\* reduction dictionary `v963` (canonical form + four-language
      elimination + decomposition bookkeeping + honest margin decay);
      `INVENTORY.json`, `make_inventory.py` (r285 sync repair included)
      and the suite MODULES list extended in the same change; the Lean
      layer carries the canonical form `lstar_subordination` with the
      proved direction L\* ⇒ free-window positivity.
- [x] **Wave 7**: rounds r286–r289 frozen as the L\* coherence census
      `v964` (margin warning resolved harmless + the exact van der
      Corput theorem for L2 + the destructive-coherence carrier map +
      the METRIC_ONLY adjudication); `INVENTORY.json`,
      `make_inventory.py` (head-sync repair included) and the suite
      MODULES list extended in the same change; the Lean layer needs no
      new statement — the hole stays `lstar_subordination`.
- [x] **Wave 8**: rounds r290–r295 frozen as the L\* curvature arc
      `v965` (basin geometry + ridge anatomy + curvature spectroscopy +
      metric reconciliation + the honest F10 rejections: F10_FRAGILE +
      F10_SP_MAJORITY — F10 unpromoted, the rejections are the frozen
      content; new open contract `PRIME.LSTAR.CLOSED_FUNCTIONAL.01`);
      `INVENTORY.json`, `make_inventory.py` and the suite MODULES list
      extended in the same change; the Lean layer needs no new
      statement — the wave carries no new formal target form, the hole
      stays `lstar_subordination`.
- [x] **Wave 9**: rounds r296–r300 frozen as the L2 reduction chain
      `v966` (the DENS fork closed honestly: r296 DENS_WORLD_BLIND,
      the lane routed to L2; the four-round chain r297→r300 reducing
      the delta' > 0.21 target of the generic half to the one
      inequality NEFF_TARGET, with the magnitude no-go catalog and
      the two world-separating classes frozen verbatim; new open
      contract `PRIME.L2.NEFF_TARGET.01`); `INVENTORY.json`,
      `make_inventory.py` and the suite MODULES list extended in the
      same change (r301 in flight, not consumed — the sync is strictly
      additive); the Lean layer needs no new statement — NEFF_TARGET
      is a measured ladder-slope aggregate, not an exact finite
      identity, and a universal ladder form would be refutable (the
      r273 guard pattern); the hole stays `lstar_subordination`.
- [x] **Wave 10**: rounds r301–r304 frozen as the cascade closure
      `v967` (the wave-9 rest targets executed to their end: r301
      NEFF_SPLIT + r302 UNIF_DERIVED, the first DERIVED of the lane;
      the r303 regress audit RETYPING the cascade as an exact
      dictionary around one core S = +0.0547 with the hard reviewer
      rule; the r304 stop case LAW_LONGRANGE — the global-profile
      mixing route documented CLOSED, the period-4 comb the named
      open object; the rejections and retypings ARE the frozen
      content); `INVENTORY.json`, `make_inventory.py` and the suite
      MODULES list extended in the same change (consumption notes on
      the r301–r304 probes; no round in flight at this cut); the
      Lean layer needs no new statement — the stop carries no new
      formal target form (measured ladder aggregates, the r273
      guard pattern); the hole stays `lstar_subordination`.
- [x] **Wave 11**: rounds r305–r316 (incl. r310b) frozen as the
      architecture adjudication `v968` in the reviewer's binding
      four-level structure (Level 1 formal theorems — Lean already
      green, sorries 9 → 5; Level 2 certified finite statements;
      Level 3 negative architecture decisions verbatim — Lane A
      closed as cone language with the mechanism named, B not
      overtaking with the rank-one identities banked as
      `PRIME.SOURCE.RANKONE.UPDATE.IDENTITY.01` [E], C parked;
      Level 4 open mechanisms typed); the claim split is binding
      (`TENSOR.MECHANISM` [O] separate, NOT one TRANSFER title);
      the fiber GO stands with the provenance open
      (`PRIME.L2.RENYI3.PROVENANCE.01` [O], R317 fork at the
      reviewer); `INVENTORY.json` (head-sync: v968 entry +
      consumption notes r306–r316, all pre-existing entries
      byte-identical) and the suite MODULES list extended in the
      same change; the Lean layer ships NO change (r305/r310/r310b
      already merged and green); the holes stay
      `lstar_subordination` + `terminal_positive_main`.
- [x] **Wave 12**: rounds r317–r322 (the red-team / fork morning)
      frozen as `v969_forks_and_redteam.py` (the r319 extraction-chain
      audit U1–U3 + the r320 repair frozen as
      `PRIME.REDTEAM.EXTRACTION_AUDIT.01` [E]; the fiber fork frozen
      as `PRIME.L2.RENYI3.SLIDING_BOUND.01` [O]; the base fork closed
      honestly — P1 banked as language, the antiphase law
      `ALGORITHM_ARTIFACT`); the reviewer adjudication adopted (the
      two holes named WINDOW-LOCAL, the three-level proof graph A/B/C,
      the false cofinality direction documented NOT solved, the typed
      per-layer assessment table); `INVENTORY.json`,
      `make_inventory.py` and the suite MODULES list extended in the
      same change; the Lean layer ships NO change (r319/r320 already
      merged and green); rounds R324/R325 in flight, not consumed.
- [x] **Wave 13**: rounds r323–r327 (the extraction repair and the
      terminal composition) frozen as
      `v970_extraction_and_composition.py` (the r325 elementwise GO
      + the r326 Lean implementation frozen as
      `PRIME.EXTRACTION.ELEMENTWISE.01` [E] with the honest
      new-connection formulation verbatim; the r324 measured chain
      + the r327 SOURCE anatomy frozen as
      `PRIME.L2.RENYI3.MEASURED_COMPOSITION.01` [O] — subcritical
      +0.172, m₀\* record 10^59.6 / chain-honest 10^238 with bar
      0.188 per the r328B chain audit, the extrapolation gap typed;
      the r324-pre pre-work banked; the r323 abort consumed as a
      clean note); `INVENTORY.json`, `make_inventory.py` and the
      suite MODULES list extended in the same change; the Lean
      layer ships NO change (r326 already merged and green, census
      5 → 8); the audit rounds r328A/B/C read-only, reports with
      the coordinator; no measurement round in flight at this cut.
- [x] **The Level-C construction theorem** (wave-12 reviewer target,
      load-bearing): the target architecture `CanonicalPrimeWindow`
      plus exactly one construction theorem
      `sourceExact_buildPrimeWindow`, eliminating `SourceExact` as a
      free assumption (the extraction-repair fork r325 measured its
      machine-checkable halves — primary
      `ELEMENTWISE_STABILIZATION_GO`: the repair is the elementwise
      quantifier set, not a new mesh theory).
      **DONE r326** (`RH/Elementwise.lean`, see the r326 Lean-status
      paragraph): construction theorem PROVED, comb-channel
      elementwise stabilization PROVED, extraction
      `weil_nonneg_of_windowlocal` PROVED; census 5 → 8 with the
      three new typed Level-C sorrys (arch/pole kernel stabilization
      — classical; the source-exact completion — opacity-forced).
      REMAINING open there: the two kernel transcriptions, the
      bordered compression bridge as a theorem, and the Level-B
      window-local premise itself.
- [ ] Wire the paper PDF into the docset once the pipeline integration
      fully lands.
