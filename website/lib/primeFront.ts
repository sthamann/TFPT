/**
 * Live update feed for /prime-front ("The Prime Front").
 *
 * HOW TO POST A NEW ENTRY after a completed agent run:
 *   1. Open this file.
 *   2. Prepend a new object to `PRIME_FRONT_UPDATES` (newest first).
 *   3. Fill every field below; keep `summary` to 1–2 sentences.
 *      `headline` is ONE short sentence (the feed card title); `keyFacts`
 *      are 3–5 short bullet strings distilled from the run record. The
 *      long `title` stays the full diary text (shown under "Read the
 *      full entry").
 *   4. Rebuild / redeploy the website — no other files required for a feed post.
 *
 * Do not invent numbers. Copy verdicts and anchors from experiments/next.txt
 * (or the promoted verification/vN_*.py). Firewall: no "RH evidence" language.
 *
 * PERFORMANCE SPLIT: this file keeps only the newest ~16 entries inline
 * (they ship with the route's first-load JavaScript). Everything older
 * lives in `primeFrontArchive.ts`, which the feed loads lazily via a
 * dynamic import. When this inline window grows past ~20 entries, move
 * the oldest ones verbatim to the top of `PRIME_FRONT_ARCHIVE` and bump
 * `PRIME_FRONT_ARCHIVE_COUNT` below.
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
  | "TWO-OF-THREE"
  | "ARITHMETIC-DONE"
  | "HARNACK-EXPLAINED"
  | "WIDE-RESTRUCTURED"
  | "NET-IMPROVED"
  | "CBS-RESISTS"
  | "TELESCOPE-CARRIES"
  | "ASSEMBLY-GREEN"
  | "SEAMS-CERTIFIED"
  | "BOTH-SHAPED"
  | "THREE-OF-FOUR"
  | "KAPPA-WILD"
  | "SUPPLY-PARTIAL"
  | "SPECTRUM-ONLY"
  | "BOUNDED-STATE"
  | "ONE-CARRIES"
  | "BOTH-RESIST"
  | "PAIR-EXACT"
  | "DENSE-RESISTS"
  | "FINITE-CORE"
  | "HARDY-RESISTS"
  | "PROFILE-RESISTS"
  | "SHARP-CARRIES"
  | "INTERVAL-CARRIES"
  | "ONE-STEP-MISSING"
  | "LEVEL-CARRIES"
  | "ASYMPTOTIC-SHAPED"
  | "ONE-INPUT-MISSING"
  | "PARTIAL-SMOOTHING"
  | "ONE-TERM-MISSING"
  | "ODD-CARRIES"
  | "REBUILD-RESISTS"
  | "ALIGN-RESISTS"
  | "FLOORS-RESIST"
  | "TERMS-RESIST"
  | "ANGLES-RESIST"
  | "ENTRY-RESISTS"
  | "FORM-RESISTS"
  | "PAIRING-RESISTS"
  | "CLOSURE-RESISTS"
  | "DELTA-REDUCED"
  | "FRONT-RESISTS"
  | "TOLERANCE-CARRIES / SECTOR-RESISTS"
  | "ETA-RESISTS"
  | "CASCADE-RESISTS"
  | "VECTOR-RESISTS"
  | "MINORS-RESIST"
  | "RATIO-RESISTS"
  | "SIEVE-RESISTS"
  | "MAP-COMPLETE"
  | "PARTIALLY-PORTABLE"
  | "DEFICIT-VARIES"
  | "PARTIAL-CANCEL"
  | "PHASES-RESIST"
  | "DEGENERATE"
  | "SITS-AT-ZERO"
  | "FIBER-PINNED"
  | "ONE-MODE"
  | "FLOOR-CERTIFIED"
  | "LONG-RANGE"
  | "FORMULAS-EXACT"
  | "NULL-RAY-CORRELATED"
  | "SIGNS-MATCH"
  | "OCCUPATION-FOLLOWS-KERNEL"
  | "TRANSPORT-KILLED"
  | "DENSITY-DRIVEN"
  | "PRIME-FREE-CLOSED-FORM"
  | "UNIT-BOX-SEPARATED"
  | "TWO-LAYER-SPLIT"
  | "DIRECTION-DENSITY-FIXED"
  | "CLOSED-FORM-EXACT"
  | "CLOSED-DELTA-EXACT"
  | "ZERO-COMB-IDENTIFIED"
  | "EXISTENCE-UNCONDITIONAL"
  | "POLE-RANK-ONE"
  | "DET-LAW-DERIVED"
  | "CUTOFF-CLOSES"
  | "UNCONDITIONAL-ENTRY-CERT"
  | "MAPPING-CLOSES"
  | "NEAR-IDENTITY-1D"
  | "FINGERPRINTS-HIT"
  | "DICT-THREE-ENTRIES"
  | "REAL-STRUCTURE-LANDED"
  | "JOINT-EMBEDDING-LANDED"
  | "EQUIVARIANT-DUAL-LANDED"
  | "DUALITY-FORMS-LANDED"
  | "SEAM-REAL-LANDED"
  | "EQUIVARIANT-ORDER-LANDED"
  | "TRANSLATOR-LANDED"
  | "MODES-SEPARATED"
  | "SELECTION-RULE-LANDED"
  | "GAMMA-FREE-SLICE-LANDED"
  | "C-FACTORIZATION-LANDED"
  | "LANDSCAPE-MAPPED"
  | "PERIODS-CLASSICAL"
  | "POLARIZATION-FRAMED"
  | "PERIOD-NORMALIZED"
  | "TRANSPORT-CANONICAL"
  | "GAMMA-TOY-LANDED"
  | "SU2-KINEMATIC-SEPARATED"
  | "SEAM-FORCES-COVER"
  | "UNIFORM-C1"
  | "FLIPS-ARE-TRUNCATION"
  | "N3-PINNED"
  | "CLOCK-SURVIVES-INTERACTING"
  | "SEAM-IDENTIFIED"
  | "COVERED-SEAM-SKELETON"
  | "LATTICE-AUDIT-PASSED"
  | "PRIME-SHADOW-EXACT"
  | "SIX-FRONT-BATCH"
  | "W1-DICTIONARY-DERIVED"
  | "TEN-FRONT-BATCH"
  | "ERRATUM-PLUS-SIX-BATCH";

export type PrimeFrontUpdate = {
  /** ISO date (YYYY-MM-DD) of the agent run. */
  date: string;
  /**
   * Diary part number (Teil N).
   * Use `0` for meta / promotion announcements that span multiple diary parts
   * (string anchors are not supported by this field).
   */
  part: number;
  /** Full diary text of the run (long; rendered as the expandable body). */
  title: string;
  /** One short sentence — the compact feed-card title. */
  headline?: string;
  /** 3–5 short fact bullets distilled from the run record. */
  keyFacts?: readonly string[];
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
    date: "2026-08-02",
    part: 0,
    title:
      "The second round of the day \u2014 the Lerch erratum and six promotions (v643\u2013v648). THE ERRATUM, honest and prominent: the W1 chain (v631, v640\u2013v642) read Suzuki's eq. (1.3) with Lerch coefficient \u22121; the paper's own \u00a72.2 data lock +1/4 (machine evidence: v643 check C0.1 \u2014 the printed origin constant A = \u00bd(log 2\u03c0 \u2212 \u03c8(2)) = 0.70754637 reproduces only with +1/4, and g(0) = 0). Every identity of that chain is a correct identity of its kernel g\u0303 = g \u2212 (5/4)\u00b7Lerch, and every measured number transfers VERBATIM via the exact transfer identity cgal_sm(g\u0303) = \u22124\u00b7cgal_sm(g) \u2014 no check flips, only the labels change: Suzuki's own smooth layer is +\u03c1 (not \u22124\u03c1), the true dictionary is the SINGLE scalar +1/D on both layers (sign-compatible with positivity), and the origin constant vanishes, \u03ba = 0 exactly \u2014 the v563 near-cell scheme IS Suzuki's origin bookkeeping. The erratum is documented as a dated correction across the modules, the ledger, the papers and this site. ON THE CORRECTED READING, SIX PROMOTIONS: (1) v643 proves the W1 MEASURE-LEVEL THEOREM on Suzuki's true screw function \u2014 the projection lemma closes the last named remainder (Suzuki's L\u00b2\u2080 mean-zero condition is automatic on the u-side, D: H\u00b9\u2080 \u2192 L\u00b2\u2080 a bijection; the window parity sector is the odd H\u00b9\u2080 subspace), A_arch = \u2212g\u2033_smooth exactly at every lag (3.4e\u221252), the moment law is derived to D\u2078, and the full form equality holds at 1.28e\u221210 on the common odd sector (11/11, W1-THEOREM). (2) v644 starts W2 honestly \u2014 the classical FEM-density slice verified at rate (H\u00b9 1.000, L\u00b2 2.000), Rayleigh\u2013Ritz monotone from above on all six nested sequences, \u03bb_a = 0\u207a within ~1e\u22129 across five windows (consistent with Suzuki Thm 1.3, no sign statement, ground state even reported), the Mosco/norm-resolvent remainder typed \u2014 started, not closed (7/7, W2-SLICES-PASS). (3) v645 resolves the named RP open item of the \u2124\u2083 orbifold: the twist sector IS reflection positive under the parafermionic Klein-twist pairing \u2014 mirror composed with charge conjugation, \u03c9^{kq} zero-mode phases, \u03b7 = (1,1) unique among 36 \u2014 N-stable 48\u2013384, both axes, including the full mixed OS Gram; the naive pairing keeps its \u22124.5e\u22123 violation as must-fail (17/17, KLEIN-RP-POSITIVE). (4) v646 answers the reverse code question two ways: the \u03c3-fixed syndromes carry the anchor decomposition (0,1,1,2)/norm 6 exactly with the anchor = family sum (the selection equivariance-specific, the weights generic), while the bytecode p-sequence p_n = 2 + 2\u207f has NO code reading (0/11 families; 18/27/81 without a hit \u2014 the bingo is buried) (20/20). (5) v647: the \u03bc\u2084 center is a POWER of both regular clocks (w\u2076 = \u00b1J with \u03c7_w = x\u2074 + ix\u00b2 \u2212 1 exactly primitive \u03b6\u2082\u2084; u\u2075 = \u00b1J central), and c = w\u00b2 is KILLED by a parity theorem \u2014 the compiler clock's census 19\u00d712 + 3\u00d74 has odd counts, impossible even in S\u2082\u2084\u2080; u\u2074 acts freely with census {5:48}, 48 = 240/5 (25/25). (6) v648 types the W3 tool diagnosis: the sign-uncertainty/Mellin-strip toolbox has a real 25-digit dictionary to the critical strip (phase slope = L = log \u03c0 \u2212 \u03c8(\u00bc), the W1 \u03b4\u2080 weight) but its mass lever dies at d = 1 (e^{\u2212E*} = 0.967 where < \u00bd is needed; d \u2265 21 required; 5.2\u00d7 scale mismatch) \u2014 side findings: \u03bb_min > 0 on 67/67 complete windows (min +8.26e\u22124), the incomplete combs are exactly the v618 flip set {1219, 1292, 1445}, and corr(log \u03bb_min, log \u03b5h) = +0.704 (12/12, LP-CEILING-ANALOG-TYPED). Contract state after the round: W1 theorem-closed, W2 started, W3 open with a tool diagnosis, W4 unchanged \u2014 no RH statement.",
    headline:
      "The second round of the day \u2014 the Lerch erratum (honest, same-day) and six promotions: W1 becomes a theorem, W2 starts, the Klein twist carries positivity, and the W3 toolbox is typed.",
    keyFacts: [
      "ERRATUM (same day): the W1 chain read eq. (1.3) with Lerch \u22121, correct is +1/4 (v643 C0.1) \u2014 every number transfers verbatim via cgal(g\u0303) = \u22124\u00b7cgal(g); true dictionary: ONE scalar +1/D, \u03ba = 0 exactly",
      "v643\u2013v648 promoted: W1 theorem (11), W2 start (7), Klein RP (17), RM(1,3) reverse (20), ST31 degree 24/20 (25), sign-uncertainty W3 diagnosis (12) \u2014 92 new checks, all green",
      "v643: A_arch = \u2212g\u2033_smooth exactly (3.4e\u221252), projection lemma closes the L\u00b2\u2080 remainder, form equality 1.28e\u221210 \u2014 W1-THEOREM",
      "v645: \u03b7 = (1,1) unique among 36 \u2014 the \u2124\u2083 twist sector is reflection positive under the Klein pairing (N-stable, full mixed OS Gram); v647: c = w\u00b2 killed by a parity theorem",
      "v648: the W3 mass lever dies at d = 1 (e^{\u2212E*} = 0.967, needs < \u00bd); surface positive 67/67 (min +8.26e\u22124). Suite 636 \u2192 642 scripts",
    ],
    verdict: "ERRATUM-PLUS-SIX-BATCH",
    summary:
      "v643\u2013v648 promoted: v643_w1_theorem.py (PRIME.WEIL.THEOREM.01, 11 checks \u2014 incl. the erratum evidence C0.1), v644_w2_form_density.py (PRIME.WEIL.W2.01, 7), v645_klein_rp.py (QGEO.KLEINRP.01, 17), v646_rm13_reverse.py (E8.RM13REV.01, 20), v647_st31_degree24.py (E8.ST31DEG.01, 25), v648_sign_uncertainty.py (PRIME.SIGNUNC.01, 12). Dated convention erratum on v631/v640\u2013v642 (labels corrected, numbers unchanged \u2014 the transfer identity is exact). Suite 636 \u2192 642 scripts.",
    badge: "sandbox",
    script: "w1_theorem_probe.py",
  },
  {
    date: "2026-08-02",
    part: 0,
    title:
      "The ten-front batch \u2014 eleven modules land at once, with three Lean certificates and two honest negatives promoted as typed results. (1) The PGL\u2082 contract executed: the four external transfer solvers share ONE M\u00f6bius action on the disc \u22127 norm line where they act at all \u2014 the Koide reading is the parabolic step +9 = N_fam\u00b2, exactly intertwined with the base translation; g_QCD is exactly conjugate to it (a translation by +7/2\u03c0 in the 1/\u03b1 coordinate); Boltzmann = pole identically by the single-flow theorem; the preregistered kill criterion 'four incompatible arithmetic actions' is NOT met. (2\u20133) The \u2124[i]-E8 quotient and its stabilizer: the \u03bc\u2084-quotient of the 240 roots IS the hermitian-unimodular \u2124[i]-E8 with the 60-hyperplane system of ST31; G31 is the FULL unitary stabilizer (order 46080, exactly counted) with Molien-unique degrees (8,12,20,24) and the exact clock identities \u03c3 = c\u2074, J = c\u2079 \u2014 while the numerology |G31| = |W(D5)|\u00d7|W(A3)| is KILLED three ways (the real structure is (Z4\u22182^{1+4}).Sp\u2084(2), and W(D5) does not even embed in a rank-4 group); honest sharpening: the compiler clock is NOT the Springer-regular 12-element. The 60\u21926 cascade route on the quotient is killed too. (4\u20136) P canonical, constructed \u2014 and an honest negative: the Lorentz congruence matrix P is the unique minimal-Frobenius operator-compatible class in the full [\u22124,4] census AND, mod sign, the Frobenius-minimal integer congruence transporting the C_V null frame (= the two minimal isotropic rays of J_fix) onto the canonical rank-one rays (a finite 16-member family); but after h-control the fine Hodge invariants carry NO independent information about the C = 1 margin (\u03c1 = \u22120.04, p = 0.75) \u2014 the direct window-level geometry\u2192arithmetic bridge route closes. (7) The code dictionary: the unique equivariant Hamming placement is Reed\u2013Muller RM(1,3) on AG(3,2) \u2014 information bits one per \u03bc\u2084 pair (3 families + 1 anchor), syndrome flag with the \u03c3 3-cycle, decode = projection verbatim (3840/3840) \u2014 and all of it dies for every non-equivariant placement (2/30). (8) The twist OPE: the interacting \u2124\u2083 orbifold slice stands at the abelian vertex level \u2014 h_\u03c3 = 1/36 from three exact routes, two-point exponent 2\u0394 = 2/9 at 0.06%, crossing exactly symmetry-protected on the lattice, closed four-point form, OPE c\u2081 = 2/9 model-bound at 0.19%; reflection positivity is validated in the real \u2124\u2082 class and honestly OPEN for the complex \u2124\u2083 pairing (parafermionic Klein twist). (9\u201311) The W1 closure trio: the boundary equation A_arch = (1/4)g\u2033_smooth \u2212 (5/4)(log \u03c0 \u2212 \u03c8(\u00bc))\u03b4\u2080 is EXACT (residuals < 1e\u221222; \u03c8(\u00bc) = \u2212\u03b3 \u2212 3log2 \u2212 \u03c0/2 sympy-exact) with the derived window-independent 1 + 1/(6d\u00b2) moment law; the frozen dictionary transports UNCHANGED to three fresh windows (atom constant D\u00b2 at machine precision \u2014 the preregistered kill criterion is met nowhere); and the full quadratic form closes at the operator level (block norm 2.7e\u22123 derived, 6.3e\u22126 after measure-level re-binning; pole block rank 1; the literal per-vector 1% bar fails only on two TYPED lattice residuals, certified as typed-residual checks with the probe's numbers unchanged). W1 is theorem-capable at the window level; the L\u00b2\u2080 projection is the last remainder; W2/W3 stay open and unmoved. Plus: Lean certificates for the W1 dictionary identities, the Lorentz congruence, and the G31 order identities (lake build green, no sorry).",
    headline:
      "The ten-front batch — eleven modules land at once, with three Lean certificates and two honest negatives promoted as typed results.",
    keyFacts: [
      "v632–v642 promoted — eleven modules in one batch, incl. two honest negatives as typed results",
      "W1 closure trio: boundary equation exact (< 1e−22), frozen dictionary portable to three fresh windows, full quadratic form closed at operator level (6.3e−6 after re-binning)",
      "ST31: G31 is the full unitary stabilizer of the μ₄ quotient (order 46080), σ = c⁴ and J = c⁹ exact — the |W(D5)|×|W(A3)| numerology killed three ways",
      "Code dictionary: the unique equivariant Hamming placement is RM(1,3) on AG(3,2) — one bit per μ₄ pair, 3 families + 1 anchor, decode = projection verbatim (3840/3840)",
      "Three Lean certificates (WeilDictionary, LorentzCongruence, G31Orders; no sorry). Suite 625 → 636 scripts",
      "Erratum (same day, v643 C0.1): the W1 trio read eq. (1.3) with Lerch −1 (correct: +1/4) — every number transfers verbatim; true dictionary: one scalar +1/D, κ = 0",
    ],
    verdict: "TEN-FRONT-BATCH",
    summary:
      "v632\u2013v642 promoted: v632_ftransfer_pgl2.py (FTR.PGL2.01, 25 checks), v633_orbit60_quotient.py (E8.ZIQUOTIENT.01, 17), v634_st31_structure.py (E8.ST31.01, 56), v635_p_canonicity.py (PRIME.PCANON.01, 18), v636_p_construction.py (PRIME.PCONSTRUCT.01, 25), v637_fine_c1_bridge.py (PRIME.FINEC1.01, 8 \u2014 honest negative), v638_code_semantics.py (E8.CODESEM.01, 27), v639_twist_ope.py (QGEO.TWISTOPE.01, 25), v640_w1_boundary.py (PRIME.WEIL.BOUNDARY.01, 11), v641_w1_portability.py (PRIME.WEIL.PORTABLE.01, 10), v642_w1_matrix.py (PRIME.WEIL.MATRIX.01, 8). Three Lean certificate files (WeilDictionary, LorentzCongruence, G31Orders). Suite 625 \u2192 636 scripts.",
    badge: "sandbox",
    script: "w1_matrix_identity_probe.py",
  },
  {
    date: "2026-08-02",
    part: 0,
    title:
      "The W1 dictionary \u2014 the mystery residual was the pole term all along, and the conversion is now DERIVED. Hours after v630 measured a non-scalar conversion between the TFPT window symbol and Suzuki's screw-function Galerkin matrix, the residual is resolved exactly. The key is a small exact identity with large consequences: term-wise (2n+\u00bd)\u00b2/(n+\u00bc)\u00b2 = 4, so the second derivative of the Hurwitz\u2013Lerch block of the screw function collapses to a geometric series \u2014 [e^{\u2212t/2}\u03a6(e^{\u22122t},2,\u00bc)]\u2033 = 4e^{\u2212t/2}/(1\u2212e^{\u22122t}). That yields a structure theorem: off the prime atoms, g\u2033(t) = \u22122cosh(t/2) \u2212 4e^{\u2212t/2}/(1\u2212e^{\u22122t}). Read it term by term: the second piece is exactly \u22124 times the TFPT archimedean density (the Weil 1952 kernel v563 has always used), and the first piece \u2014 2cosh(t/2) = e^{t/2} + e^{\u2212t/2} \u2014 is the zeta POLE block, the s = 0, 1 weights of the explicit formula. Suzuki's screw function bundles the pole into g; TFPT has always tracked it separately as an exact rank-one piece (v591). So the v630 drift was never a mystery normalization: it was the pole term sitting inside g. The verification: the closed per-lag ratio law \u2212D[4 + (e^t+1)(1\u2212e^{\u22122t})] matches the measured profile, and after subtracting the pole block from g the conversion becomes the SCALAR \u22124D, converging monotonically (1.0006 at lag 16, residual = the declared d\u207b\u00b2 discretization). The atom layer converts with the constant D\u00b2 exactly \u2014 identical for the atoms at log 2 and log 3. Bottom line: W1 of the PRIME.WEIL.OPERATOR contract \u2014 identify the TFPT window form with the Galerkin matrix of Suzuki's localized Weil operator \u2014 is now closed at the measured level on BOTH halves: atoms literal (v630), smooth layer scalar after pole separation (v631). The RH-hard steps W2/W3 are untouched, and the implication map stays explicit: closing W1 does not move W3.",
    headline:
      "The W1 dictionary — the mystery residual was the pole term all along, and the conversion is now DERIVED.",
    keyFacts: [
      "v631_w1_dictionary.py promoted (PRIME.WEIL.DICT.01, 7 checks; probe w1_dictionary_probe.py 7/7, W1-DICTIONARY-DERIVED)",
      "Lerch collapse exact",
      "Structure theorem g″ = −2cosh(t/2) − 4×(TFPT arch density)",
      "Closed ratio law verified",
      "Suite 624 → 625 scripts",
      "Erratum (same day, v643 C0.1): eq. (1.3) carries Lerch +1/4, not −1 — the identities hold for g̃ = g − (5/4)·Lerch, numbers transfer verbatim; Suzuki's smooth layer is +ρ, the true dictionary the single scalar +1/D",
    ],
    verdict: "W1-DICTIONARY-DERIVED",
    summary:
      "v631_w1_dictionary.py promoted (PRIME.WEIL.DICT.01, 7 checks; probe w1_dictionary_probe.py 7/7, W1-DICTIONARY-DERIVED): Lerch collapse exact; structure theorem g\u2033 = \u22122cosh(t/2) \u2212 4\u00d7(TFPT arch density); closed ratio law verified; pole-subtracted conversion scalar \u22124D (\u2192 1.0006); atom constant D\u00b2. W1 closed at measured level; W2/W3 open. Suite 624 \u2192 625 scripts.",
    badge: "sandbox",
    script: "w1_dictionary_probe.py",
  },
  {
    date: "2026-08-02",
    part: 0,
    title:
      "The six-front batch \u2014 E8 is literally a code, the windows sit in one Hodge chamber, and the TFPT atoms ARE Suzuki's prime measure. Five modules land at once. (1) The E8 code: Construction A on the self-dual extended Hamming code [8,4,4] yields E8 exactly \u2014 even, unimodular, shells 240/2160 \u2014 and every single-bit error is exhaustively correctable (16\u00d78 corrupted words, unique nearest codeword every time). 'E8 is an error-correcting code' is now a THEOREM in the suite, giving the robustness narrative its exact anchor next to v625's commuting Hecke channels. (2) The Hodge chamber: transporting every complete window through the exact Lorentz congruence P\u1d40J_det P = J_fix, all 67 windows land in the POSITIVE cone of the cover polarization lattice \u2014 on ONE sheet. The geometric positivity route has its measured half; honest typing: scrambles do NOT leave the chamber (membership is density-layer, per v582), so the fine C = 1 arithmetic lives INSIDE the chamber. (3) The orbifold Casimir: the covered seam's deck classes {1/6, 1/2, 5/6} carry the \u2124\u2083-twist Casimir data EXACTLY \u2014 closed form for the vacuum energy, sector table (1/72, \u22121/24, 1/72), twist gap 1/18. (4) The incidence census: the review's 48\u00d75 \u2192 240 kill test fires \u2014 the compiler-canonical order-12 element gives orbits 19\u00d712 + 3\u00d74, not 20 free \u03b6\u2081\u2082 orbits \u2014 but the positive residue is sharp: the \u03bc\u2084 clock alone acts FREELY with exactly 60 orbits, and 60 = D_start. The cascade start appears as a free orbit count on the E8 roots. (5) The Suzuki first contact: the atom layer of the TFPT window form IS the prime measure of Suzuki's screw function \u2014 positions log n, weights \u039b(n)/\u221an, literally, exactly \u2014 so W1's atomic half of the PRIME.WEIL.OPERATOR contract is CLOSED; the smooth-layer conversion is measured non-scalar, the preregistered residual, now with data.",
    headline:
      "The six-front batch — E8 is literally a code, the windows sit in one Hodge chamber, and the TFPT atoms ARE Suzuki's prime measure.",
    keyFacts: [
      "v626–v630 promoted: v626_e8_code.py (E8.CODE.01, 7 checks — Hamming [8,4,4] → E8, exhaustive syndrome decoding)",
      "v627_hodge_chamber.py (PRIME.HODGECONE.01, 5 checks — 67/67 windows in the positive cone, one sheet)",
      "v628_orbifold_casimir.py (QGEO.ORBCAS.01, 5 checks — exact twist Casimir table, gap 1/18)",
      "v629_root_incidence.py (E8.INCIDENCE.01, 6 checks — naive 48×5 killed, 60 = D_start free clock orbits)",
      "Suite 619 → 624 scripts",
    ],
    verdict: "SIX-FRONT-BATCH",
    summary:
      "v626\u2013v630 promoted: v626_e8_code.py (E8.CODE.01, 7 checks \u2014 Hamming [8,4,4] \u2192 E8, exhaustive syndrome decoding); v627_hodge_chamber.py (PRIME.HODGECONE.01, 5 checks \u2014 67/67 windows in the positive cone, one sheet); v628_orbifold_casimir.py (QGEO.ORBCAS.01, 5 checks \u2014 exact twist Casimir table, gap 1/18); v629_root_incidence.py (E8.INCIDENCE.01, 6 checks \u2014 naive 48\u00d75 killed, 60 = D_start free clock orbits); v630_suzuki_contact.py (PRIME.WEIL.CONTACT.01, 4 checks \u2014 atom identity literal, Galerkin profile measured). Suite 619 \u2192 624 scripts.",
    badge: "sandbox",
    script: "suzuki_contact_probe.py",
  },
  {
    date: "2026-08-02",
    part: 0,
    title:
      "The prime shadow — primes enter AFTER the geometry, exactly. An external note asked: what if primes are not the origin but the readout — the shadow of the finished geometry in discrete arithmetic? The checkable core of that reading turns out to be EXACT, on the compiler's own objects. The E8 theta function, computed from the glue decomposition (θ₂⁸ + θ₃⁸ + θ₄⁸)/2, has shell counts r(2n) = 240·σ₃(n) — the first shell is literally 240, the root count — so Θ_E8 = E₄: the finished lattice's counting function is a modular form. The 'address space' framing is unique factorization, exactly: shell counts factor over coprime addresses (and the must-fail control shows non-coprime does NOT factor: σ₃(4) = 73 ≠ 81). The 'independent check channels' framing is a theorem: the Hecke operators T_p act with eigenvalue 1 + p³ for every prime, they COMMUTE, and the compiler's theta is a simultaneous eigenvector of all of them — each prime reads the same finished object through its own independent channel. And the zeta shadow: L(E₄, s) = ζ(s)·ζ(s−3) — the Riemann zeta function appears as the FACTORIZED SHADOW of the E8 counting function. So the chain μ₄ → D5⊕A3 → E8 → theta → Hecke → ζ → primes is exact at every arrow: within TFPT's narrative, the direction of explanation is fixed — geometry first, primes as readout. Honest scope: these are classical facts (Jacobi, Hecke); the content is that the compiler's own objects realize them verbatim. The bolder framings — primes as compiler eigenfrequencies, RH as maximal coherence — stay typed hypotheses, not adopted.",
    headline:
      "The prime shadow — primes enter AFTER the geometry, exactly.",
    keyFacts: [
      "v625_prime_shadow.py promoted (PRIME.SHADOW.01, 6 checks; probe prime_shadow_probe.py 6/6, PRIME-SHADOW-EXACT)",
      "Θ_E8 = E₄ (shells 240·σ₃(n), n = 1..12 exact)",
      "Multiplicativity over coprime addresses (must-fail σ₃(4) ≠ σ₃(2)²)",
      "Hecke channels T_p with eigenvalue 1+p³ commute, theta simultaneous eigenvector",
      "Suite 618 → 619 scripts",
    ],
    verdict: "PRIME-SHADOW-EXACT",
    summary:
      "v625_prime_shadow.py promoted (PRIME.SHADOW.01, 6 checks; probe prime_shadow_probe.py 6/6, PRIME-SHADOW-EXACT): Θ_E8 = E₄ (shells 240·σ₃(n), n = 1..12 exact); multiplicativity over coprime addresses (must-fail σ₃(4) ≠ σ₃(2)²); Hecke channels T_p with eigenvalue 1+p³ commute, theta simultaneous eigenvector; L(E₄,6) = ζ(6)ζ(3) to 5e−8. Speculative framings typed hypotheses. Suite 618 → 619 scripts.",
    badge: "sandbox",
    script: "prime_shadow_probe.py",
  },
  {
    date: "2026-08-02",
    part: 0,
    title:
      "The external-review lattice audit — d = 4 forced three ways, and the prime front and the cover share one Lorentz lattice. A third external review proposed three derivations; all three verify exactly, and one follow-up gets an honest retype. The four-mark closure: generalize the compiler's glue to D_{d+1} ⊕ A_{d-1} with d marks. Three INDEPENDENT selectors each force d = 4: the discriminant-group census (a cyclic glue needs the D-side cyclic, i.e. d even, and ℤ_d ≅ ℤ_4 pins d = 4 — negative controls at 3, 5, 8), the even-glue norm equation (d+1)/4 + (d−1)/d = 2 (unique positive solution 4), and a Pythagorean gem: (d−1)² + d² = (d+1)² holds exactly at d = 4 — the (3, 4, 5) triple of cycles, marks and carrier is the unique Pythagorean closure. With (1,1,2) the only 3-part partition of 4 and c₃ = 1/(2e₁π) = 1/(8π), the two-input presentation reduces to ONE boundary degree d = 4 plus π — conditional on the glue architecture, honestly typed. The anchor bytecode: one cubic χ = (t−1)²(t−2) and a two-state recursion p_{n+2} = 3p_{n+1} − 2p_n carry the whole discrete source (240 + 8 = 248). And the surprise: the Lorentz congruence. The prime-front determinant form J_det (det 2) and the cover's integer polarization lattice J_fix (det 72) are the SAME rational quadratic form — an explicit integer matrix P with det −6 gives PᵀJ_det P = J_fix exactly, index 6 = p₂(a). A genuine new bridge between prime analysis and Hodge geometry. Honest negative: the review's follow-up test toward B = M†J is ill-posed as stated (B is complex symmetric, det 8 verified). And the review's RH architecture — identify the window forms with Suzuki's localized Weil operator (arXiv 2606.09096, citations verified), then prove the uniform C = 1 contraction — is preregistered as the contract PRIME.WEIL.OPERATOR.01 with kill tests. No RH statement.",
    headline:
      "The external-review lattice audit — d = 4 forced three ways, and the prime front and the cover share one Lorentz lattice.",
    keyFacts: [
      "v624_external_lattice_audit.py promoted (EXTREV.LATTICE.01, 14 checks; probe external_lattice_audit_probe.py 14/14, LATTICE-AUDIT-PASSED)",
      "Four-mark closure forces d = 4 three ways (conditional on the glue architecture)",
      "Anchor bytecode exact (240+8 = 248)",
      "Lorentz congruence PᵀJ_det P = J_fix exact (index 6)",
      "Suite 617 → 618 scripts",
    ],
    verdict: "LATTICE-AUDIT-PASSED",
    summary:
      "v624_external_lattice_audit.py promoted (EXTREV.LATTICE.01, 14 checks; probe external_lattice_audit_probe.py 14/14, LATTICE-AUDIT-PASSED): four-mark closure forces d = 4 three ways (conditional on the glue architecture); anchor bytecode exact (240+8 = 248); Lorentz congruence PᵀJ_det P = J_fix exact (index 6); B-test retyped ill-posed; Suzuki arXiv 2606.09096 + 2607.24830 verified; new contract PRIME.WEIL.OPERATOR.01 (W1–W4 chain, kill tests). Suite 617 → 618 scripts.",
    badge: "sandbox",
    script: "external_lattice_audit_probe.py",
  },
  {
    date: "2026-08-01",
    part: 0,
    title:
      "The covered seam — r⁴ = ω becomes lattice combinatorics, and the 12-grid is the mode grid. The last bedrock residue (the covering-level identification) gets its lattice skeleton. Lift the 16-site seam circle to the μ₃-cover: crossing a mark-bond shifts the sheet. With the uniform weights the covered seam is ONE 48-site NS circle (the walk closes after exactly 48 steps, covering the base three times); with the alternating weights (1,2,1,2) it DISCONNECTS into three circles — a new, fully independent reason the alternating pair dies (v617 killed it by mark equivalence, v620 by the census; connectivity kills it structurally). The centerpiece: the lift L of the base quarter turn satisfies L⁴ = deck and L¹² = id on all 48 points, with proper order 12. The relation r⁴ = ω — established algebraically in v597, derived by Gauss–Manin transport in v614 — is now LATTICE COMBINATORICS: going around the base once (four lifted quarter turns) is exactly one deck step. And the ζ₁₂ spectrum of the Burau rotation is simply the lifted clock's order. The kernel identity of v622 generalizes verbatim to N = 48, and the NS modes split under the deck into three classes of 16 with offsets {1/6, 1/2, 5/6} — the 12-grid, with the offset-1/2 class equal to the BASE seam frequencies exactly (the untwisted sector is the base, verbatim) and the twisted pair carrying ω on bilinears (the t = ω/ω̄ sheets). A bonus structure: the fermionic deck lift satisfies deck³ = −1 — the NS full rotation, the spin double, resonating with the (−1)^F dichotomy of v510/v519. The 12 lifted mark crossings form a single 12-cycle under the clock. What remains: the interacting orbifold statement (twist-field OPE) — the last named piece of the bedrock.",
    headline:
      "The covered seam — r⁴ = ω becomes lattice combinatorics, and the 12-grid is the mode grid.",
    keyFacts: [
      "v623_covered_seam.py promoted (QGEO.COVERSEAM.01, 9 checks; probe covered_seam_probe.py 9/9, COVERED-SEAM-SKELETON)",
      "One 48-site covered seam (uniform weights",
      "Alternating pair disconnects — independent kill)",
      "Lifted clock exact order 12 with L⁴ = deck (r⁴ = ω as lattice combinatorics)",
      "Suite 616 → 617 scripts",
    ],
    verdict: "COVERED-SEAM-SKELETON",
    summary:
      "v623_covered_seam.py promoted (QGEO.COVERSEAM.01, 9 checks; probe covered_seam_probe.py 9/9, COVERED-SEAM-SKELETON): one 48-site covered seam (uniform weights; alternating pair disconnects — independent kill); lifted clock exact order 12 with L⁴ = deck (r⁴ = ω as lattice combinatorics); N = 48 kernel identity; 12-grid mode split, untwisted = base verbatim, twisted pair = ω on bilinears; deck³ = −1 (spin double); 12 marks in a single 12-cycle. GATE.QGEO does not move. Suite 616 → 617 scripts.",
    badge: "sandbox",
    script: "covered_seam_probe.py",
  },
  {
    date: "2026-08-01",
    part: 0,
    title:
      "The seam identification — the physical seam IS the conformal seam. The bedrock premise of the cover program had two named residues after v617/v620; the second one — identify the PHYSICAL seam (the 16-Majorana circle of the RP program) with the CONFORMAL seam (the unit circle with μ₄ marks that forces the curve) — closes now at the kernel and dictionary level. The centerpiece is an exact identity: the v519 chiral seam kernel c(d) = (2/N)/sin(πd/N) for odd distances, zero for even, is EXACTLY the antiperiodic (NS) chiral mode sum on the 16-site circle — for all fifteen distances, including the even-distance zeros. The flat-band structure that the whole seam program is built on IS the NS spectrum of a chiral Majorana on the circle (closed form: Σ sin((2j+1)x) = sin²(8x)/sin(x)). The must-fail control: the Ramond (periodic) mode sum does NOT reproduce the kernel — NS is forced, not a convention. The geometry follows: placing the sites at half-integer angles θ_a = (2a+1)π/16, the four mark-bond midpoints land EXACTLY on μ₄ = {1, i, −1, −i} and no site touches a mark — v519's 'marks at bond midpoints' becomes literal geometry. The group dictionary is exact (the clock is z → iz, the RP reflections are diameter reflections), and the v534 straddle law turns out to be geometry: the straddled cuts {7, 15} are exactly the axes that pass THROUGH marks, the avoiding cuts {3, 11} the axes between marks. What remains of the bedrock: the covering-level identification — the μ₃-cover of the seam double carrying the full interacting theory. The gate does not move, but both v617 residues are now processed: one conditional (N = 3 pinned by the compiler's constants), one closed at kernel level.",
    headline:
      "The seam identification — the physical seam IS the conformal seam.",
    keyFacts: [
      "v622_seam_identification.py promoted (QGEO.SEAMID.01, 10 checks; probe seam_identification_probe.py 10/10, SEAM-IDENTIFIED)",
      "v519 kernel = discrete NS mode sum exactly (all d, even zeros included",
      "Ramond must-fail)",
      "Mark-bond midpoints = μ₄ exactly",
      "Suite 615 → 616 scripts",
    ],
    verdict: "SEAM-IDENTIFIED",
    summary:
      "v622_seam_identification.py promoted (QGEO.SEAMID.01, 10 checks; probe seam_identification_probe.py 10/10, SEAM-IDENTIFIED): v519 kernel = discrete NS mode sum exactly (all d, even zeros included; Ramond must-fail); mark-bond midpoints = μ₄ exactly; group dictionary exact (clock = z→iz, cuts = diameters {π/4, π/2, 3π/4, π}); straddle law = axis-through-marks geometry; residue: covering-level identification. GATE.QGEO does not move. Suite 615 → 616 scripts.",
    badge: "sandbox",
    script: "seam_identification_probe.py",
  },
  {
    date: "2026-08-01",
    part: 0,
    title:
      "The flip mechanics — the two violators were truncation artifacts, and C = 1 is exception-free. The v618 dichotomy (violator set = sign-flip set) demanded a mechanism, and it turned out to be embarrassingly concrete: the atom demand of a window runs to u ≤ 2α, the prime-power data cap sits at U_max = 12.899 — and EXACTLY the two flip windows are the two whose demand exceeds the cap (h = 1219 missing Δu = 0.677 of comb, h = 1445 missing 0.089). The flip set equals the truncation set exactly. Three controls make it airtight: the zero-comb coupling spectrum at the flips is regular-collective (share 5.8%, N50 = 57 — no resonance, so the flip is NOT an arithmetic event); truncation INJECTION into complete windows reproduces the flips with matched sign and magnitude at both scales (removing Δu = 0.089 from healthy windows gives q_real ≈ −5e-5..−1e-4, matching 1445's observed −1.0e-4; removing 0.677 gives −0.18..−0.22, matching 1219's −0.29); and the overshoot is monotone in the truncation depth. Conclusion: on the complete-comb surface (2α ≤ U_max, 67 windows) the C = 1 bound holds with ZERO exceptions (max eps·h = 0.982). The 'sign-flip windows' are retired as data-boundary artifacts — extending the surface requires extending the prime-power data, not the theory.",
    headline:
      "The flip mechanics — the two violators were truncation artifacts, and C = 1 is exception-free.",
    keyFacts: [
      "v619_flip_mechanics.py promoted (PRIME.FLIPMECH.01, 6 checks; probe flip_mechanics_probe.py 6/6, FLIPS-ARE-TRUNCATION)",
      "Flip set = truncation set exactly (2α > U_max = 12.899)",
      "No zero resonance",
      "Injection reproduces both scales",
      "Suite 614 → 615 scripts",
    ],
    verdict: "FLIPS-ARE-TRUNCATION",
    summary:
      "v619_flip_mechanics.py promoted (PRIME.FLIPMECH.01, 6 checks; probe flip_mechanics_probe.py 6/6, FLIPS-ARE-TRUNCATION): flip set = truncation set exactly (2α > U_max = 12.899); no zero resonance; injection reproduces both scales; C = 1 exception-free on the 67 complete windows (max 0.982). Suite 614 → 615 scripts.",
    badge: "sandbox",
    script: "flip_mechanics_probe.py",
  },
  {
    date: "2026-08-01",
    part: 0,
    title:
      "The cyclic-N census — the seam admits {3, 5}, and the compiler's constants pin N = 3. v617 forced the base (μ₄ marks, Möbius-rigid) and the uniform weights, leaving 'why cyclic-3' as the first bedrock residue. The census over N = 2..6 makes it conditional. Full ramification everywhere (the structure the whole cover program uses) holds iff gcd(N,4) = 1 — that is N ∈ {3, 5}. The genus table is 1, 3, 3, 6, 7, and the primitive character sheets have the compiler's rank 3 exactly for N ∈ {3, 5, 6}. Then Chevalley–Weil decides: N = 3 has BOTH primitive sheets Lorentz ((1,2) and (2,1) — precisely the v599/v613 signature); N = 5 has a Lorentz pair among four sheets; and N = 6 has NO primitive Lorentz sheet — its Lorentz content sits entirely on the order-3 characters, i.e. N = 6 reduces to N = 3. So the seam-admissible orders with a primitive Lorentz sheet are exactly {3, 5}. The pinning is then integral: the deck determinant on the rank-3 module is det(1 − deck) = N³ — 27 for N = 3, which is EXACTLY the machine-checked v597 constant, versus 125 for N = 5; and the saturation analogue is 81 = 3⁴ (the v566/v600 compiler self-code index) versus 625. Among {3, 5}, the compiler's integral constants select N = 3 uniquely. The residue moves one level down: why the compiler is 3-adic is N_fam = 3, anchored in the E8 glue (240 = 16·5·3).",
    headline:
      "The cyclic-N census — the seam admits {3, 5}, and the compiler's constants pin N = 3.",
    keyFacts: [
      "v620_cyclic_n_census.py promoted (QGEO.NCENSUS.01, 9 checks; probe cyclic_n_census_probe.py 9/9, N3-PINNED)",
      "Full ∞-ramification iff gcd(N,4)=1",
      "Genus table 1,3,3,6,7",
      "Primitive Lorentz sheets exactly for N ∈ {3,5} (N=6 reduces to N=3)",
      "Suite 613 → 614 scripts",
    ],
    verdict: "N3-PINNED",
    summary:
      "v620_cyclic_n_census.py promoted (QGEO.NCENSUS.01, 9 checks; probe cyclic_n_census_probe.py 9/9, N3-PINNED): full ∞-ramification iff gcd(N,4)=1; genus table 1,3,3,6,7; primitive Lorentz sheets exactly for N ∈ {3,5} (N=6 reduces to N=3); pinning det(1−deck) = N³: 27 (v597) vs 125; saturation 81 vs 625 — N = 3 unique. GATE.QGEO does not move. Suite 613 → 614 scripts.",
    badge: "sandbox",
    script: "cyclic_n_census_probe.py",
  },
  {
    date: "2026-08-01",
    part: 0,
    title:
      "The interacting semigroup — the clock survives, the circle does not. v524 built the Klein–Landau local symmetric semigroup on the free net: seven transfer steps τ_k, all Hermitian, even steps positive, odd steps indefinite. What happens on the interacting survivor? The free baseline reproduces exactly at N = 16 (even steps PD; odd steps k = 1, 3, 5 indefinite with EXACTLY ZERO one-particle diagonal — the free chirality datum; k = 7 parity-trivial). Then the interaction speaks: at g > 0 exact Hermiticity survives PRECISELY on the steps {4, 6, 7} and is lost on {1, 2, 3, 5} — the semigroup contracts to the μ₄-commensurate window. The mechanism is exact Clifford combinatorics: the survivor interaction is α₄-invariant EXACTLY and not α₁/α₂-invariant — translation invariance breaks to the quartet stabilizer {0, 4, 8, 12}, and the surviving steps are the symmetry-protected ones (k = 5 shows quartet containment alone does not protect: even parity AND alignment are both needed). The gamma-relevant centerpiece: the clock step k = 4 stays Hermitian and positive definite over the whole coupling grid — the μ₄ clock is a positive transfer operator on the interacting net. The reading: the reconstructed rotation group of the interacting toy is the CLOCK TOWER, not the full circle — exactly the structure the TFPT compiler needs, emerging from RP + the alignment bit alone.",
    headline:
      "The interacting semigroup — the clock survives, the circle does not.",
    keyFacts: [
      "v621_interacting_semigroup.py promoted (WOIT.GAMMA.SEMI.01, 6 checks; probe interacting_semigroup_probe.py 6/6, CLOCK-SURVIVES-INTERACTING)",
      "Free v524 pattern at N = 16 reproduced",
      "At g > 0 Hermiticity survives exactly on {4, 6, 7}",
      "Clock step PD (11,0,0) over g ∈ {1/32..8}",
      "Suite 612 → 613 scripts",
    ],
    verdict: "CLOCK-SURVIVES-INTERACTING",
    summary:
      "v621_interacting_semigroup.py promoted (WOIT.GAMMA.SEMI.01, 6 checks; probe interacting_semigroup_probe.py 6/6, CLOCK-SURVIVES-INTERACTING): free v524 pattern at N = 16 reproduced; at g > 0 Hermiticity survives exactly on {4, 6, 7}; clock step PD (11,0,0) over g ∈ {1/32..8}; mechanism: α₄-invariance exact (quartet stabilizer). Gamma proper stays open on A_hol. Suite 612 → 613 scripts.",
    badge: "sandbox",
    script: "interacting_semigroup_probe.py",
  },
  {
    date: "2026-08-01",
    part: 0,
    title:
      "The uniform constant, frozen — C = 1 on the declared surface, and the violators are exactly the sign flips. The equidistribution conjecture (sec:theory-open) asked for |q_real/q_model| ≤ C·h⁻¹ uniformly; the measured constant is now a frozen, reproducible module. On the declared frame-A surface (69 floor-passed windows, h = 142..1445) the model value q_model keeps ONE sign on the whole ladder (0.039..0.112 — no model zero crossing anywhere), and on every lock-sign window (67 of 69) the deviation times h stays below 1: C = 1, measured max 0.982 at h = 184, tertile medians 0.61/0.45/0.39 non-increasing with depth. The sharpened typing is the real news: exactly TWO windows carry a q_real SIGN FLIP (h = 1219 blowing up to 9.2e3, and the edge window h = 1445 at 3.5) — and the violator set of the C = 1 bound equals the flip set EXACTLY. The sign predicts the bound violation window-sharp. This comes with an honest correction: the earlier diary note attributed the h = 1219 blow-up to a q_model zero crossing (wrong — q_model never crosses zero on the ladder) and missed the edge flip at h = 1445 entirely. Scrambled combs break the bound by more than four orders of magnitude: the constant is genuine arithmetic placement. No uniformity proof, no RH statement — the conjecture now has a frozen measured anchor.",
    headline:
      "The uniform constant, frozen — C = 1 on the declared surface, and the violators are exactly the sign flips.",
    keyFacts: [
      "v618_uniform_constant.py promoted (PRIME.UNIFC.01, 6 checks; probe uniform_constant_probe.py 6/6, UNIFORM-C1)",
      "69 windows, q_model one-signed (0.039..0.112)",
      "Eps·h ≤ 0.982 on all 67 lock-sign windows (C = 1)",
      "Tertiles 0.61/0.45/0.39",
      "Suite 611 → 612 scripts",
    ],
    verdict: "UNIFORM-C1",
    summary:
      "v618_uniform_constant.py promoted (PRIME.UNIFC.01, 6 checks; probe uniform_constant_probe.py 6/6, UNIFORM-C1): 69 windows, q_model one-signed (0.039..0.112); eps·h ≤ 0.982 on all 67 lock-sign windows (C = 1); tertiles 0.61/0.45/0.39; exactly two q_real sign-flip windows (1219, 1445) = the violator set exactly (XXIII note corrected); scramble breaks by > 1e4. Suite 611 → 612 scripts.",
    badge: "sandbox",
    script: "uniform_constant_probe.py",
  },
  {
    date: "2026-08-01",
    part: 0,
    title:
      "The seam forcing round — the conformal seam axioms force the cover. The cover program spent eighteen rounds showing that the curve y³ = x⁴ − 1 CARRIES the compiler; the bedrock question (QGEO.SYM.01) runs the other way: why THIS geometry? Now the direction reverses. Z4-Möbius rigidity: any four points on the sphere that a Möbius transformation of order 4 permutes cyclically are Möbius-equivalent to μ₄ = {1, i, −1, −i} — the multiplier is forced to a primitive fourth root, and the cross-ratio census is the harmonic orbit {2, 1/2, −1} exactly (must-fail control: 4/3, the disc-7 jet datum, is NOT in the orbit — two distinct conformal data). The weight census: among all 81 ways to put μ₃-monodromy on the four marks, exactly four are clock-equivariant, and the demand that the four marks be EQUIVALENT kills the alternating pair — the weights are FORCED uniform, j ∈ {1, 2}, which is precisely the conjugate sheet pair (t = ω vs ω̄) of the period dictionary. Uniform weight gives y³ = (x−1)(x−i)(x+1)(x+i) = x⁴ − 1 exactly, full ramification at infinity (the separated infinity-cusp of v603), and Riemann–Hurwitz genus 3. And the seam geometry matches piece by piece: the marks lie ON the unit circle, the doubling reflection fixes seam and marks, and the OS-cut (bond) reflection z ↦ −i·z̄ permutes the marks by exactly the (k, 5−k) pattern of the v599 real structure, meeting the circle at the two v519 cut points. What remains of the bedrock: why cyclic-3 (N_fam = 3, anchored in the E8 glue 240 = 16·5·3), and the physical-seam ↔ conformal-seam identification. The gate does not move — but the residue just got much smaller.",
    headline:
      "The seam forcing round — the conformal seam axioms force the cover.",
    keyFacts: [
      "v617_seam_cover_forcing.py promoted (QGEO.SEAMFORCE.01, 14 checks; probe seam_cover_forcing_probe.py 14/14, SEAM-FORCES-COVER)",
      "Z4-Möbius rigidity (multiplier forced ±i",
      "Harmonic CR orbit {2, 1/2, −1}, 4/3 excluded)",
      "Weight census 81 → 4 → 2 (uniform j ∈ {1,2} forced)",
      "Suite 610 → 611 scripts",
    ],
    verdict: "SEAM-FORCES-COVER",
    summary:
      "v617_seam_cover_forcing.py promoted (QGEO.SEAMFORCE.01, 14 checks; probe seam_cover_forcing_probe.py 14/14, SEAM-FORCES-COVER): Z4-Möbius rigidity (multiplier forced ±i; harmonic CR orbit {2, 1/2, −1}, 4/3 excluded); weight census 81 → 4 → 2 (uniform j ∈ {1,2} forced); y³ = x⁴−1 exact; genus 3 (Riemann–Hurwitz); bond reflection = v599 (k,5−k) pattern with the two v519 cut points. GATE.QGEO does not move. Suite 610 → 611 scripts.",
    badge: "sandbox",
    script: "seam_cover_forcing_probe.py",
  },
  {
    date: "2026-08-01",
    part: 0,
    title:
      "The internal SU(2) separates — kill test (6) gets its first slice. The last untouched gamma kill test asks whether the internal SU(2) 'remains a spacetime factor'. At the kinematic twistor level the answer is now exact linear algebra. Twistor space C⁴ is a two-dimensional quaternionic module, and the split is clean: LEFT quaternion multiplications are the spacetime side (they move the fibers of PT → S⁴), RIGHT quaternion multiplications are the internal side (they preserve every fiber — they act trivially on spacetime). The surprise of the conventions: Woit's Euclidean structure ρ_tw IS the internal j-direction, verbatim. The internal su(2) commutes exactly with the full spacetime action and intersects it in {0}: the internal factor SEPARATES — the kinematic branch of kill test (6) does not fire. Then the structure of the clock itself: RHO4 = left-e^{iπ/4} ∘ right-e^{iπ/4} EXACTLY — the μ₄ clock entangles a spacetime half-quarter turn with an internal half-quarter phase, neither factor alone is the clock, the ALE deck is purely spacetime-side (DECK4 = left-i), and RHO4² = deck × internal quarter phase. Finally the breaking: the complex-linear part of the internal algebra is exactly the U(1) generated by right-i — choosing the Euclidean section breaks the internal SU(2) to U(1) kinematically, and the v519 mark torsor is exactly the μ₄ subgroup of that internal U(1). The dynamical half (a genuine gauge action on the reconstructed net) stays open on A_hol — but the kinematic geometry is now settled, and it is consistent with the electroweak reading.",
    headline:
      "The internal SU(2) separates — kill test (6) gets its first slice.",
    keyFacts: [
      "v616_su2_internal_kinematic.py promoted (WOIT.SU2.KIN.01, 17 checks; probe su2_internal_kinematic_probe.py 17/17, SU2-KINEMATIC-SEPARATED)",
      "Ρ_tw = right-j verbatim",
      "Internal su(2) fiber-preserving, commutes with spacetime, intersection {0} — kill test (6) kinematic branch does not fire",
      "RHO4 = left-u ∘ right-u exact (u = e^{iπ/4})",
      "Suite 609 → 610 scripts",
    ],
    verdict: "SU2-KINEMATIC-SEPARATED",
    summary:
      "v616_su2_internal_kinematic.py promoted (WOIT.SU2.KIN.01, 17 checks; probe su2_internal_kinematic_probe.py 17/17, SU2-KINEMATIC-SEPARATED): ρ_tw = right-j verbatim; internal su(2) fiber-preserving, commutes with spacetime, intersection {0} — kill test (6) kinematic branch does not fire; RHO4 = left-u ∘ right-u exact (u = e^{iπ/4}); DECK4 = left-i; internal SU(2) → U(1) by the complex choice, mark torsor = internal μ₄. Dynamical half stays live on A_hol. Suite 609 → 610 scripts.",
    badge: "sandbox",
    script: "su2_internal_kinematic_probe.py",
  },
  {
    date: "2026-08-01",
    part: 0,
    title:
      "The interacting gamma slice — chirality survives the interaction, exactly where RP survives. The WOIT gamma milestone asks: does the one-generation-without-mirrors structure survive on the INTERACTING algebra? v608 answered at the free level; this round answers on the only interacting model the corpus has — the 16-Majorana Fidkowski–Kitaev seam toy, at the ONE interaction that reflection positivity allows (the alignment bit δ = π/2 with positive coupling, the unique survivor of v534's selection law). The results, over the whole coupling grid g ∈ {1/32 … 8}: the chiral OS Grams stay positive definite on all four cuts (the survivor anchor), and at the SAME forced twist the MIRROR state splits sector-exactly — even sector positive, odd sector STRICTLY NEGATIVE DEFINITE, bounded away from zero (max odd eigenvalue ≤ −1.5e−3). Every mirror fermion vector has strictly negative OS norm at the interacting level: no mirror mode survives reconstruction — kill test (3) cannot fire on the surviving interaction class. The mark transport also survives: the quarter turn maps the mark-A quartet algebra exactly onto the mark-B sites, and both the mark Gram and the transport form are Hermitian positive definite at every coupling — kill test (7) does not fire. The odd transfer spectrum stays multiplicity-free at every g (honest finding: at strong coupling g = 8 the eigenvalues pair up to gaps ~2e−5 without closing — a near-doubling tendency worth watching). And the control that makes the story sharp: the RP-dead member (m = 2) has NO gamma home — its straddled cut is indefinite. The dynamical selection of the alignment bit and the chirality data CO-LOCATE on the same interaction: reflection positivity does not just pick the bit, it picks the member where chirality is protected. Fence: one toy, one interaction class; gamma proper stays open on A_hol; kill test (6) untouched.",
    headline:
      "The interacting gamma slice — chirality survives the interaction, exactly where RP survives.",
    keyFacts: [
      "v615_gamma_toy_interacting.py promoted (WOIT.GAMMA.TOY.01, 11 checks; probe gamma_toy_interacting_probe.py 11/11, GAMMA-TOY-LANDED)",
      "Chiral/mirror parents I ± iC",
      "Chiral Grams PD (37,0,0) on all four cuts, g ∈ {1/32..8} (v534 regression)",
      "Mirror odd sector (0,8,0) strictly ND, bounded ≤ −1.5e−3 — kill test (3) cannot fire on the survivor",
      "Suite 608 → 609 scripts",
    ],
    verdict: "GAMMA-TOY-LANDED",
    summary:
      "v615_gamma_toy_interacting.py promoted (WOIT.GAMMA.TOY.01, 11 checks; probe gamma_toy_interacting_probe.py 11/11, GAMMA-TOY-LANDED): chiral/mirror parents I ± iC; chiral Grams PD (37,0,0) on all four cuts, g ∈ {1/32..8} (v534 regression); mirror odd sector (0,8,0) strictly ND, bounded ≤ −1.5e−3 — kill test (3) cannot fire on the survivor; mark transport α₄ exact 16/16, Gram + transport form PD at every g — kill test (7) does not fire; odd spectrum multiplicity-free (min gap 2.3e−5); dead member m = 2 indefinite. Kill tests (3)/(6)/(7) stay live on A_hol. Suite 608 → 609 scripts.",
    badge: "sandbox",
    script: "gamma_toy_interacting_probe.py",
  },
  {
    date: "2026-08-01",
    part: 0,
    title:
      "The Gauss–Manin transport — the deck step IS the boundary monodromy. v613 left one residue: the dictionary r = deck∘rotation was established at the conjugacy level (equal characteristic polynomials, an explicit conjugator) — but WHICH deck step, and why, was carried by a matrix, not by geometry. Now it is geometry. Transport the segment cycles along the rigid quarter-rotation loop of the four branch points (p_k(τ) = e^{iπτ/2}·p_k — the rotation braid): the discriminant λ(τ) = e^{2πiτ} winds ONCE around zero, and the transported periods factorize EXACTLY — I_m(τ) = e^{iπτ(m+1)/2}·e^{−2πiτj/3}·I_m(0). At τ = 1 the second factor is the deck step, and it is ω = t⁴ on EVERY sheet — the BOUNDARY MONODROMY of the local system (four punctures of weight t each). The must-fail controls make the counting sharp: t³ = 1 ≠ ω and t⁵ ≠ ω — the exponent 4 = |μ₄| is forced by the puncture count. Consequence: the canonical transport matrix in the segment basis is EXACTLY ω·M — no conjugator freedom left, the v613 dictionary is CANONICAL, and (ωM)⁴ = ω·1 re-derives v597's r⁴ = ω one more level down. The picture is simple: the braid fixes the boundary of the disk, the rigid rotation does not — their homological difference is exactly one boundary loop's worth of local-system monodromy. Certificates: direct branch-tracked integration at four intermediate τ matches the closed factorization to 5e−28 (dps 80), and the transported segment at τ = 1 lands on ω × its static neighbor to 2.7e−81.",
    headline:
      "The Gauss–Manin transport — the deck step IS the boundary monodromy.",
    keyFacts: [
      "v614_transport_deck.py promoted (QGEO.TRANSPORT.01, 12 checks; probe transport_deck_probe.py 12/12, TRANSPORT-CANONICAL)",
      "Transported periods factorize exactly (rotation substitution identity)",
      "Deck factor = ω = t⁴ on every sheet, uniform across all three rows",
      "Must-fail t³/t⁵ controls (exponent 4 = |μ₄| forced)",
      "Suite 607 → 608 scripts",
    ],
    verdict: "TRANSPORT-CANONICAL",
    summary:
      "v614_transport_deck.py promoted (QGEO.TRANSPORT.01, 12 checks; probe transport_deck_probe.py 12/12, TRANSPORT-CANONICAL): transported periods factorize exactly (rotation substitution identity); deck factor = ω = t⁴ on every sheet, uniform across all three rows; must-fail t³/t⁵ controls (exponent 4 = |μ₄| forced); canonical transport matrix = ωM (char = char(r) exact); 25-digit certificates (max 5e−28, τ=1 at 2.7e−81). v613 residue F1 closed. GATE.QGEO does not move. Suite 607 → 608 scripts.",
    badge: "sandbox",
    script: "transport_deck_probe.py",
  },
];

/**
 * Number of entries in `primeFrontArchive.ts` (kept here so the feed can
 * show honest totals before the archive chunk has loaded).
 */
export const PRIME_FRONT_ARCHIVE_COUNT = 225;

/**
 * Display order for the narrative sections (anchor ids), grouped for the
 * on-page navigation. Order = reading order: the results first, then the
 * diary arc as it happened, the live feed last.
 */
export const PRIME_FRONT_SECTION_GROUPS = [
  {
    label: "Start here",
    sections: [
      { id: "big-picture", label: "Big picture" },
      { id: "suzuki-w1", label: "The W1 theorem" },
      { id: "uniform-constant", label: "C = 1" },
    ],
  },
  {
    label: "The diary, as it happened",
    sections: [
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
      { id: "compression", label: "The compression" },
      { id: "certification-sprint", label: "Certification sprint" },
      { id: "harnack-telescope", label: "Harnack & telescope" },
      { id: "meaning", label: "What it would mean" },
      { id: "prime-shadow", label: "Prime shadow" },
      { id: "e8-code", label: "E8 code" },
      { id: "hodge-chamber", label: "Hodge chamber" },
      { id: "sixty-lines", label: "Sixty lines" },
    ],
  },
  {
    label: "Archive",
    sections: [{ id: "updates", label: "Live updates" }],
  },
] as const;
