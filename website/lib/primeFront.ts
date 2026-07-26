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
  | "BAND-PARTIAL";

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
    date: "2026-07-26",
    part: 91,
    title:
      "The band's law — the first prime rescues positivity; named target inequality (T)",
    verdict: "BAND-PARTIAL",
    summary:
      "T91 (BAND-PARTIAL, 19/19, 3/4 closed): the band is the one-atom zone log 2 < a ≤ 0.9253 with inner edge a_neg = 0.7410 — and the geography's cleanest finding: beyond a_neg the prime-free margin changes sign and THE PRIME ATOM BECOMES LOAD-BEARING — the first prime rescues positivity where the archimedean margin is exhausted (the T89 balance point is exactly this sign change, 0.03%). Atom turn-on law exact (k = 2m+1, Beta integrals); uncertainty constants decided: a·t_rms → π (Wirtinger) and a·t_cent → 2Si(π) − 4/π = 2.4306 (new closed constant); band and tight curves are two orbit regions of the same functional with shared exact scale ∫k_ζ = 2θ_RS (Riemann–Siegel). Named target inequality (T): (P_pole + A_arch)(f) ≥ √2·log2·h_f(log2) on the band — provable-shaped as a ZONE EXTENSION beyond Bombieri's log 2 (a self-standing classical target!); RH ⇒ (T), (T) ⇏ RH. Honest open: the super-exponential λ_pf rate remains empirical. Classical: Wirtinger/Rayleigh, Beta integrals, Lambert-W, Si integral, Bombieri, θ_RS. Sandbox; not RH evidence.",
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
  { id: "meaning", label: "What it would mean" },
  { id: "updates", label: "Live updates" },
] as const;
