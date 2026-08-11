# The end-form of a localized Weil-positivity program, and an exact inertia bridge to the positive-index method

**Stefan Hamann — TFPT project, prime/zeta line ("The Prime Front")**
2026-08-11 · research documentation · [fixpoint-theory.com/prime-front](https://fixpoint-theory.com/prime-front)

---

## Abstract

This note summarizes, in plain language, the current state of the RH-side of the TFPT prime/zeta program, for readers of the 2026 Anthropic result that two thirds of the nontrivial zeta zeros are simple and lie on the critical line. **We state the frame first and without hedging: we have no proof of the Riemann Hypothesis, no partial zero-density result, and we claim no progress toward RH.** What we do have, and what this note documents, is: (1) an exact reduction of a Weil-positivity criterion, through a Galerkin discretization of Suzuki's localized Weil operator, to the positive-definiteness of one small matrix update per rung of a window ladder, split exactly into a 7×7 co-block half and a two-scalar half; (2) a rigorous, interval-arithmetic finite head — positivity certified on all 42 reachable rungs of the deployed ladder; (3) a certified floor for the co-block half over the ideal (not merely computed) source objects on the entire reachable surface; (4) a registered, frozen, falsifiable target inequality for the open half, which passed its first genuine blind holdout; (5) a complete measured map of the certificate languages that *fail* on the open half — twelve-plus dead routes, each with a named mechanism; and (6) — new, and the reason this note exists — an exact algebraic dictionary showing that our registered target is, verbatim, a *positive-index statement in the Anthropic currency*: an effective-rank/rank-trace inequality for a fixed 8×8 core, where integrality at dimension 8 replaces the asymptotic positive-proportion problem. The dictionary is machine-verified algebra; whether the target itself is provable is exactly as open as before. Everything here is reproducible: every load-bearing claim is carried by a module in a 901-script, audit-gated verification suite (module identifiers `vN` throughout), and everything exploratory is labelled as such.

---

## 1. What the object is

### 1.1 The window form

The program's central object is a family of finite symmetric matrices ("window forms") indexed by a window depth *h*: Galerkin discretizations, on a ladder of prime-power-adapted windows, of Suzuki's localized Weil operator. This identification is not an analogy; it is theorem-closed at the measure level in the suite (`v643_w1_theorem.py`): the atom layer of the window form *is* the prime measure of Suzuki's screw function, literally — positions at logarithms of prime powers, weights Λ(n)/√n (`v630_suzuki_contact.py`) — and the archimedean layer matches exactly at every lag (identity verified to ~10⁻⁵² at the kernel level, full form equality at 1.28×10⁻¹⁰). The smallest eigenvalue of the window form at depth *h* is the "wall margin" τ_h; empirically τ_h > 0 on every reachable window, decaying like a measured power law (currently p ≈ −1.93 over 95 rungs).

### 1.2 The extraction chain, and what "cofinal positivity ⇒ RH" means

The continuum extraction is the implication chain

> **cofinal finite positivity ⇒ Weil positivity ⇒ RH,**

which is *theorem-grade modulo named classical citations* in the program's bookkeeping (`PRIME.EXTRACTION.CHAIN.01`; the second arrow is classical — Weil 1952, Bombieri 2000, Suzuki). "Cofinal" is the minimal hypothesis, kernel-checked in Lean 4 (`TfptCarrier/CofinalWeil.lean`, 11 declarations): positivity is needed not uniformly and not at every rung, but on a cofinal family, with the hierarchy *uniform ⊊ pointwise ⊊ cofinal* proved strictly in both directions inside the kernel. The honest typing, stated in the source documents and repeated here: this chain is an exact *reformulation* at RH strength (`PRIME.ERRORTERM.SCALE.01`: no unconditional slack anywhere in the reduction). Proving cofinal positivity for all *h* is the Riemann Hypothesis, full stop. Every finite statement below is a finite statement.

### 1.3 The posture

The parent document (`tfpt_prime_front.tex`, ~470 pages) is a *research documentation*, not a results paper: it separates machine-verified, ledger-carried modules from exploratory sandbox probes, and it carries the fence "no claim of progress on the Riemann Hypothesis" as a governing statement (deliberately repeated thirteen times). This note inherits that posture verbatim.

---

## 2. What is rigorously proven

### 2.1 The finite head: 42 interval-rigorous rungs (`v897_certified_interval_ladder.py`)

The positivity σ_h > 0 of the deployed window ladder is **proven on all 42 reachable rungs h = 142…878**, with the earlier informal evaluation-error model *retired*:

- **15 exact-rational rungs (h ≤ 300):** integer Bareiss/Sylvester certificates on a rational grid, certified floors 4.53×10⁻⁵ … 2.14×10⁻⁴. Every transcendental evaluation (Gauss–Legendre nodes, archimedean kernel, atom positions) is enclosed in rigorous `mpmath.iv` intervals — the node-enclosure lemma runs a definite sign change of P₄₈ on every node interval at radius 10⁻⁹⁰ — and the decision path contains *no working-precision rounding anywhere*.
- **27 validated-precision rungs (h > 300):** interval-shifted Cholesky certificates at 120/200 digits, certified floors 5.16×10⁻⁶ … 8.35×10⁻⁵, honestly labelled "validated-precision" and never upgraded to exact-rational.
- **The negative control fires:** the identical machinery *refuses* the Epstein x²+5y² comb (no Euler product) at pivot index 10 — the certificates are reading the arithmetic, not the assembly.

What this is: a rigorous finite computation. What it is not: any statement about h → ∞. The asymptotic tail is untouched and carried as open contracts (`PRIME.PORT.TAIL.01`, `PRIME.PORT.LEADING.SIGN.01`).

### 2.2 The exact end-form: one matrix, two halves (`v901_tangent_schur_bfloor.py`)

Along the soft (minimal-eigenvalue) direction v_h of the rung-h form — backward spectral data only — the normalized update M_h splits *exactly* in tangent-Schur coordinates. Write

```
M = [ n   b* ]      n = v*Mv (scalar),  b ∈ R^7,  B = 7×7 co-block,
    [ b   B  ]      q = b* B⁻¹ b  (the Schur load).
```

Then, **with B ≻ 0: M ≻ 0 ⟺ n − q > 0.** All identities are warded at machine precision (det M = (n−q)·det B at 2.8×10⁻¹⁴; the equivalence boolean-exact on 39/39 ladder steps), and both scalars are *source-only*: no forward margin, no forward sign enters — the sign of the gap is an output, never an input. This is the end-form: the whole positivity question per transition is (i) the co-block half B ≻ 0 and (ii) the scalar half n > q.

### 2.3 The certified B-half: B ⪰ 0.5523·I over ideal objects (`v905_bfloor_ideal_certificate.py`)

The co-block half is now a **certified surface theorem**. The mechanism is worth stating because it is the round's lesson:

- Every classical certified bound applied to B directly is *dead*: Gershgorin (raw and scaled), Brauer–Cassini, Weyl are negative on all 39 steps, best bound-maximum −88.2 against the true measured floor +0.679 (`v901`, verdict CERTFLOOR-DEAD — registered as an honest negative).
- The missing piece was not a better inequality but the *certificate class*. The chain **B ⪰ ½·P_G + c_dom·I, P_G ⪰ c_G·I ⟹ B ⪰ c_B·I**, where P_G is the source-only Christoffel–Darboux Gram co-block of the wall's own positive chain, is decided inequality-by-inequality by *exact-rational LDL* (Fraction arithmetic on the float64 entries, which are exact rationals — the `v897` certificate class) on all 39/39 reachable steps: c_B min/med 0.5914/5.92, efficiency c_B/λ_min(B) ≥ 0.871. The structural surprise: **P_G is Gershgorin-positive on its own** (min +0.33 on 39/39) — the classical bounds that are dead on B are alive on the positive Gram substructure; the single non-classical link is the dominance decision.
- The interval rollout then lifts the statement from the computed matrices to the **ideal source objects**: the full pipeline (lags → FFT density → Lanczos chain → frame → Schur complements → B) carried as rigorous outward enclosures, *fail-closed* at every step (an interval too wide to decide refuses rather than passes). Six former width-refusals were closed by an Ogita–Rump high-accuracy tier (error-free TwoSum block accumulation, radius improvement ~150–300×). Result: **min c_B = 0.5523 on 39/39 ideal-object steps** (93.4% of the float reference), 0 refused / 0 failed / 0 skipped, and the smooth control world *refused* by the identical exact machine (0/35).

Honest scope, stated everywhere in the module: this is a *surface* statement — certified on the reachable 39-step ladder with a frozen frame convention; h-uniformity beyond h ~ 900 is open, and the n > q half is untouched by it.

### 2.4 Grade of these results

Sections 2.1–2.3 are machine-checked, ledger-carried, and audit-gated (the suite's exit gate requires byte-consistency across papers, ledger, and public mirror). They are finite linear algebra and interval arithmetic — rigorous, reproducible, and *not* number theory beyond the finite prime tables they consume. No RH input is used anywhere (firewalled: no zeros, no prime-table oracles; AST-checked inside the probes).

---

## 3. The registered target: n − q ≥ ½·μ₁(h) (`v907_halfgap_registered_target.py`)

For the open half, the program maintains one **registered falsification instrument** — explicitly *not evidence*:

> **Target.** n_h − q_h ≥ ½·μ₁(h), where μ₁(h) = 4·sin²(π/(2h+1)) is the deployed Dirichlet/KMS parity normalization.

The registration discipline is the content:

- The constant is **exactly ½, frozen**, under an explicit **no-adjust clause** that is part of the registered object: a later miss may not be repaired by 0.49, by reweighting, or by rung exclusion — a fail is a first-class result.
- **Surface: 67/67 pass**, margins (ŝ − ½) min/med/max +0.0025/+0.5273/+1.6845; tightest rung kz 98/h 997. The second-tightest rung is *shallow* (h = 199): the tightness is arithmetic scatter, not a depth trend.
- **Registry hash** sha256 `ae292e55…` (67 lines, frozen); every future reachable rung is scored verbatim by the frozen pipeline, pass iff ŝ ≥ ½, no refits.
- **First genuine blind holdout: 28/28 pass.** A 10× deeper prime table (byte-exact against the deployed one on the overlap) opened 28 new faithful rungs h = 1219…2854 outside the registry; margins min/med +0.2232/+0.6567 — the holdout minimum sits ~90× *above* the registered surface minimum. No depth erosion. Controls discriminate 3/3 (smooth, Epstein, scramble all sit far below ½).
- The four candidate algebraic origins of the ½ are typed **open**; the registration *demands* a derivation of the constant before any claim upgrade. The inequality is declared in every print as the wall margin *reparametrized* — its value is falsifiability, not novelty.

A pass is survival of a falsification attempt, never evidence. We repeat this because it is the program's own rule.

---

## 4. The resistance map: what provably fails on the open half

This is the honest core of the document. Over rounds 57–63 the program ran every certificate language it could formulate against the n > q half, under frozen preregistered specs with mandatory negative controls. The results are *measured negatives with named mechanisms* — a map of where the difficulty actually lives. (Modules `v901`, `v902`, `v906`; sandbox probes cited by round note. "Sandbox" = frozen preregistered probe, not promoted, no claim.)

| # | Route | Verdict | Named mechanism | Source |
|---|-------|---------|-----------------|--------|
| 1 | Classical eigenvalue bounds on B | dead, 39/39 | Gershgorin/Brauer/Weyl best max −88.2 vs true floor +0.679; certifiability needed a *class* change (§2.3), not a better inequality | `v901` |
| 2 | Pivot rewrite (edge-Christoffel) | relocation, not progress | exact positive-weight identity holds, but the decisive ratio's distance from 1 *tracks* the margin itself, slope +1.008: the target is the wall's PD premise quantified | `v902` |
| 3 | Gram/Lukács moment completion | dead on every rung | the wall *is* a moment-matrix statement for a signed comb measure, but PSD completion fails at the x = +1 band where comb nodes accumulate; violation outruns the margin by up to 152× | `v902` |
| 4 | Directional Löwner (floor model D for B) | dead 38/39 | Löwner step q ≤ q̄ holds *exactly*, but the loss q̄/q is med **91.3**, max 408: D is an excellent floor model and a catastrophic *inverse* model along b — D⁻¹ inflates exactly the components B⁻¹ suppresses | `v906` |
| 5 | Pointwise tail sign | dead 67/67 | net negativity is a cancellation fact; the pointwise object is positive beyond every cut (MECH-NET-ONLY) | `v906` |
| 6 | Pointwise repair / exception set | dead | no small exception set exists: repair support is diffuse (k₉₀ med 646) and grows like X^0.89 | `v906` |
| 7 | Cumulative/Abel transport | dead 67/67 | transport identities exact, but no order r ≤ 3 fixes the dual-kernel sign, and the honest all-integer envelope misses by a gap law **e^{+1.744α}**; the true remainder itself is ~10⁴ μ₁ units — the hole is oscillation, not constants | `v906` |
| 8 | Segment split at classical breakpoints | dead 67/67 | the **anchor price**: every Chebyshev-class envelope pays the deep anchor at every point of a segment; segment route is *costlier* than global (gain < 1 on 67/67); named missing object: a Brun–Titchmarsh/Selberg-class *difference* envelope | `v906` |
| 9 | Difference envelope (that named object) | dead, measured | deployed BT/monotonicity difference forms beat the absolute form by ~1% (gain med 1.010); on window lengths ~ segment span every classical one-increment envelope is too wide — what is missing is *cancellation across the window* | sandbox CLVII |
| 10 | Oscillation carriers | structural finding | the tail pairing lives on ~3 low-frequency window harmonics (j* ∈ {1,2,3} on 67/67), measurably *below* γ₁ = 14.13 (med distance 13.17 vs uniform-null ~1) — and the 12-bin reduction carries the *full* margin decay (−1.042): dimension reduced, difficulty conserved | `v906` / CLV |
| 11 | Finite harmonic certificate on the 12 bins | dead 67/67 | low-frequency smoothness loses to the envelope curse: the harmonic route is ~7× *costlier* than direct at identical envelope; the entire exponent gain (+1.74 → +1.45) comes from the envelope class and accrues to the direct route equally | sandbox CLX |
| 12 | Krylov/moment defect correction (small frame) | 37/39 but not a theorem shape | Gauss–Radau brackets are exact and optimal; med k* = 3 — but the two misses carry the largest \|b\|, and: | sandbox CLVIII/CLIX |
| 13 | — same, full frame | **unbounded** | k* grows like **h^0.80** (R² 0.77); the organizer is the b-scale: k* ~ 37·log(\|b\|/τ), R² 0.859, and \|b\|/τ itself grows ~ h. The med-k* ≈ 3 of the small frame was the frame size, not a theorem signal | sandbox CLXI |
| 14 | Transitionwise Riccati barrier (induction) | fails broadly, 21/37 | every fail is a down-flow: the μ₁-increment allowance sits 3–6 orders below the pivot flow — and the barrier is **wall-blind**: adversarial worlds that break the wall at rung level still *certify* transitions (cosh 15/36, rescale 5/17). Its entire wall content sits in the base case + the B-half premise | sandbox CLXII |

Two structural positives frame the map: **HEADFLOOR-O1** — the entire h-decay of the wall margin lives in the sign-measured atom tail; the head is O(1) and the open statement collapses to a *one-sided* bound between two O(1) objects whose gap is the margin itself (67/67) — and the 12-bin reduction of row 10, which localizes the open statement to ~12 windowed ψ-fluctuation Fourier coefficients F̂(ω), ω ≲ 5.

The summary the map supports: the open half is not blocked by a missing trick in any tested certificate language; every tested language fails by measurement with a named seat, and the failures are consistent with the object being genuinely RH-hard (as the reduction chain says it is).

---

## 5. The new bridge: the registered target is a positive-index statement

The 2026 Anthropic zeta result ([announcement](https://www.anthropic.com/research/riemann-zeta); [Lean formalization](https://github.com/anthropics/zeta-23-lean)) proves that at least 2/3 of the nontrivial zeros are simple and on the critical line, by (as we read it) moment/inertia control: trace and Frobenius-norm information about a family of operators is converted into a lower bound on the number of eigenvalues of a given sign — a **rank-trace / positive-index argument** — with an information-theoretic ceiling of ≈ 0.68 for what two moments at bandwidth 1 can yield. This section states an exact dictionary from our end-form into that currency. The algebra below is machine-verified (symbolic check in the writing session; the fixed-core identities are elementary).

### 5.1 The dictionary

Fix a ladder step and its 8×8 core M = [[n, b*], [b, B]] from §2.2, with B ≻ 0. Let c = B^{−1/2} b, so q = |c|² = b*B⁻¹b, and for a scale u > 0 define the congruence-scaled core

```
A(u) = [ u·n      √u·c* ]        A(u) = D* M' D,   D = diag(√u, I₇),
       [ √u·c     I₇    ]        M' = diag(1, B^{-1/2}) M diag(1, B^{-1/2}).
```

A(u) is congruent to M (for u > 0, B ≻ 0), so it is positive definite exactly when M is. Its two moments are closed forms: tr A = u·n + 7 and ‖A‖²_F = u²n² + 2uq + 7. Then, **exactly**:

> **(tr A)² − 7·‖A‖²_F = u·(14(n − q) − 6u·n²).**

Choose the specific scale **u_h = 7μ₁(h)/(6n²)**. At that scale the right-hand side becomes **14·u_h·(n − q − ½μ₁(h))**, hence:

> **Effective-rank ratio (tr A)²/‖A‖²_F > 7 ⟺ n − q > ½·μ₁(h).**

The left side is verbatim the Anthropic-style statistic. By the rank-trace inequality — for a symmetric matrix with tr A ≥ 0, the number n₊ of positive eigenvalues satisfies n₊ ≥ (tr A)²/‖A‖²_F — a ratio above 7 forces n₊ ≥ 8. And at dimension 8, **integrality finishes**: n₊ = 8, i.e. A(u_h) ≻ 0, i.e. M ≻ 0. (Premises: B ≻ 0 — the certified half, §2.3 — and n > 0, which holds on the measured surface; tr A > 0 is then automatic.)

So the registered half-gap target of §3 **is** a positive-index-8 statement: *the effective rank of one explicitly scaled 8×8 source-only core exceeds 7.* Not analogous — equal, by the identity above, with the frozen constant ½ appearing exactly as the choice of scale u_h.

### 5.2 Why the fixed-core version differs favorably

In the asymptotic-population setting, two-moment information yields a *proportion* of eigenvalues of the right sign, and the ceiling (≈ 0.68 at bandwidth 1) is the honest limit of that information. For a **fixed core of dimension 8**, the same two moments plus integrality yield an *all-or-nothing* conclusion: n₊ > 7 has no room to be anything but 8. The population problem is replaced by a threshold problem for one small matrix per rung — which is exactly the shape of statement our certificate machinery (exact-rational LDL, interval enclosures, Ogita–Rump tiers) is built to decide, and has decided on every reachable rung so far (§2, §3).

What this does *not* do: it does not make the target provable. The difficulty is conserved (as it was in every relocation of §4): proving ratio > 7 for all h is proving n − q > ½μ₁ for all h. What it does do is put the open statement in the *same currency* as the Anthropic method, so that their machinery — moment bounds, optimized test windows, inertia counting — applies to our fixed core directly.

### 5.3 The source-only congruence program, and the adopted no-go gate

Two concrete lines follow from the dictionary, both consuming only source-side objects:

1. **P_G-preconditioned trace/Frobenius.** The certified chain B ⪰ ½P_G + c_dom·I (§2.3) suggests computing the two moments of A(u_h) after congruence by the *positive Gram substructure* P_G — the object on which classical bounds are alive. The question becomes whether P_G-preconditioning damps the Frobenius mass of the off-diagonal block enough for the ratio to clear 7 with certified enclosures.
2. **Higher moments at dimension 8.** Two moments are far from the information limit of an 8×8 matrix: the Chebyshev–Markov–Stieltjes moment problem at dimension 8 is *finite*. The Gauss/Radau bracket machinery of round CLIX (exact lower/upper brackets on q from k moments b*Bʲb plus the certified floor c_B; brackets warded to 0 violation) is precisely the higher-moment refinement of the two-moment statistic — the measured obstruction being the growing b-scale (row 13 of §4), which any moment-side argument must now confront explicitly.

We have also adopted the Anthropic **no-go as a route gate**: their ≈ 0.68 ceiling for two moments at bandwidth 1 is registered in the program's route bookkeeping as a hard filter — any proposed argument on our side that, unwound, amounts to two-moment/bandwidth-1 information on an asymptotic population is rejected at the gate rather than measured again. The resistance map (§4) is exactly the kind of bookkeeping this gate slots into.

### 5.4 Status of the in-flight probes

At the time of writing, two preregistered probes executing the dictionary program — the **rank-trace census** (the ratio statistic across the full ladder and holdout, with controls) and the **moment-inertia bracket** (higher-moment inertia bounds on the fixed core) — are frozen and running (planned diary notes CLXIII/CLXIV) but have **not landed**; the newest recorded round-63 entry remains the Riccati barrier (CLXII). Their verdicts, whatever they are, will be recorded under the same fail-first discipline as everything above. Nothing in this note depends on their outcome.

---

## 6. What we can offer, and what we ask

**Offer.**

- **Data and frozen registries:** the 67-rung surface and the 28-rung blind holdout with the frozen registry (sha256 `ae292e55…`), the full per-rung tangent-Schur scalars (n, q, b, B), the certified floors, and the complete resistance-map measurements — all reproducible from the public suite.
- **Reproducible modules:** the verification suite (901 scripts, deterministic, audit-gated; the modules cited here are `v897`, `v901`, `v902`, `v905`, `v906`, `v907`), with negative controls (Epstein, scramble, smooth) wired into every load-bearing module.
- **A Lean 4 tree** (`TfptCarrier/`, building against Mathlib): kernel-checked modules for the extraction-chain hypotheses (`CofinalWeil.lean`, `CofinalCurrent.lean`), the grade no-go and Krein defect theorems, and per-rung wall-ladder certificates (`WallLadder/RungKz*.lean`); the composition of the finite-head certificate chain through `krein_cofinal_weil` is in progress.

**Ask.**

1. Their **moment/inertia machinery applied to our fixed core** — in particular whether their optimized bounds sharpen the two-moment threshold of §5.1, or whether higher-moment inertia at dimension 8 has a form their framework already covers.
2. Their **optimized window family tested inside our Galerkin family** — the deployed windows were never optimized in their sense; a better window changes n, q, and μ₁ together, and the dictionary says exactly what "better" must mean (ratio clearance).
3. A joint look at the **origin of the ½** — our frozen constant has four typed open candidate derivations; their scaling analysis may recognize it.

---

## 7. Honest limitations

1. **No RH claim.** Nothing here is progress on the Riemann Hypothesis, and we say so as the governing statement, not as a formality. The reduction is exact *at RH strength*: the open half is RH-hard by the program's own typing.
2. **Surface-only certificates.** Every certified statement (§2) lives on the reachable finite ladder (rungs h ≤ 878 interval-rigorous; the B-half chain exact-rational/interval on the 39 reachable steps, h ≲ 900). The 67-rung registered surface extends to h = 1433 and the blind holdout to h = 2854 at the declared deployed-pipeline (float) level only. h-uniformity is open everywhere.
3. **The target is an instrument.** The registered inequality (§3) has passed one blind holdout; that is survival of a falsification attempt, not evidence. Its constant ½ has no derivation.
4. **The dictionary conserves difficulty.** §5 is a coordinate change with a favorable integrality structure, verified as algebra; it creates a shared language, not a proof path with any certified step beyond what §2 already carries.
5. **External result as cited input.** The Anthropic 2/3-result and its ≈ 0.68 ceiling are used here as cited external inputs, pending independent verification on our side; the route gate of §5.3 is our adoption of their stated no-go, at their stated scope.

---

## References and traceability

- TFPT, *The Prime Front* (research documentation), and the verification suite: [fixpoint-theory.com/prime-front](https://fixpoint-theory.com/prime-front). Modules: `v897_certified_interval_ladder.py`, `v901_tangent_schur_bfloor.py`, `v902_wall_relocation_map.py`, `v905_bfloor_ideal_certificate.py`, `v906_tail_cartography.py`, `v907_halfgap_registered_target.py`. Sandbox rounds cited by diary note: CLV, CLVII, CLVIII, CLIX, CLX, CLXI, CLXII (2026-08-11).
- Anthropic research announcement: <https://www.anthropic.com/research/riemann-zeta> (paper PDF linked there).
- Anthropic Lean formalization: <https://github.com/anthropics/zeta-23-lean>.
- M. Suzuki, on the localized Weil operator and screw functions: J. London Math. Soc. 107 (2023); JFA 281 (2021), 109116. A. Weil (1952); E. Bombieri (2000) for Weil positivity ⟺ RH.
