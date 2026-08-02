import type { Metadata } from "next";
import Link from "next/link";
import { ArrowLeft } from "lucide-react";
import { SITE_URL } from "@/lib/utils";
import { PRIME_FRONT_SECTION_GROUPS } from "@/lib/primeFront";
import { HonestyBanner } from "@/components/primefront/HonestyBanner";
import { StorySixty } from "@/components/primefront/StorySixty";
import { WhereWeAre } from "@/components/primefront/WhereWeAre";
import { DiarySection } from "@/components/primefront/DiarySection";
import { CompilerSchematic } from "@/components/primefront/CompilerSchematic";
import { ShellCensus } from "@/components/primefront/ShellCensus";
import { AperyBridge } from "@/components/primefront/AperyBridge";
import { KillChainTimeline } from "@/components/primefront/KillChainTimeline";
import { PrimeDetector } from "@/components/primefront/PrimeDetector";
import { NeighborStepping } from "@/components/primefront/NeighborStepping";
import { TwoMelodies } from "@/components/primefront/TwoMelodies";
import { CenterAtlas } from "@/components/primefront/CenterAtlas";
import { SpectrumLadder } from "@/components/primefront/SpectrumLadder";
import { WeilArcMap } from "@/components/primefront/WeilArcMap";
import { ModuleLadder } from "@/components/primefront/ModuleLadder";
import { ConeCoverage } from "@/components/primefront/ConeCoverage";
import { DoorsPanel } from "@/components/primefront/DoorsPanel";
import { HeatFamilyBalance } from "@/components/primefront/HeatFamilyBalance";
import { CrossingMap } from "@/components/primefront/CrossingMap";
import { HandoverCrossing } from "@/components/primefront/HandoverCrossing";
import { InstrumentRace } from "@/components/primefront/InstrumentRace";
import { ClosureMapGrid } from "@/components/primefront/ClosureMapGrid";
import { TwoDoorsConvergence } from "@/components/primefront/TwoDoorsConvergence";
import { ReductionCascade } from "@/components/primefront/ReductionCascade";
import { TelescopeRungs } from "@/components/primefront/TelescopeRungs";
import { DiaryFeed } from "@/components/primefront/DiaryFeed";
import { W1DictionaryMap } from "@/components/primefront/W1DictionaryMap";
import { HodgeConeMap } from "@/components/primefront/HodgeConeMap";
import { StatusBadge } from "@/components/primefront/StatusBadge";

export const metadata: Metadata = {
  title: "The Prime Front — Research Diary",
  description:
    "A plain-language research diary of TFPT's prime / zeta line. The story in short: primes → the classical explicit formula → a window matrix built from primes inside the theory's E8 bookkeeping → proved identical to Suzuki's Weil operator (the W1 theorem, v643, machine-verified after a same-day erratum) → one open positivity statement (W3) where everything RH-hard still sits. Also: the uniform constant C = 1, the Lorentz congruence, E8 as a literal error-correcting code, and the G31 reflection group. Machine-verified modules v535–v673 inside a 667-script suite. No claim of progress toward the Riemann Hypothesis.",
  keywords: [
    "TFPT prime front",
    "E8 census",
    "Hecke from geometry",
    "Weil structure",
    "relative trace formula",
    "uniform constant C = 1",
    "Lorentz congruence",
    "Hodge chamber",
    "E8 error-correcting code",
    "Shephard-Todd G31",
    "Suzuki screw function",
    "Weil operator W1",
    "zeta research diary",
    "Stefan Hamann",
  ],
  alternates: { canonical: `${SITE_URL}/prime-front` },
  openGraph: {
    type: "article",
    title: "The Prime Front — TFPT research diary",
    description:
      "Primes, the E8 census, and the load-bearing modules v535–v673: Hecke from geometry, the phase-2 certified map, C = 1 exception-free, the Lorentz congruence, the E8 code dictionary, and the Suzuki W1 identification theorem-closed (after a same-day convention erratum). Honest fence: closing W1 does not move W3. No claim of progress toward the Riemann Hypothesis.",
    url: `${SITE_URL}/prime-front`,
    siteName: "TFPT — Topological Fixed-Point Theory",
    locale: "en_US",
  },
  twitter: {
    card: "summary_large_image",
    title: "The Prime Front — TFPT",
    description:
      "Research diary of the prime / zeta line. Machine-verified v535–v673: the induction sprint, the phase-2 map, C = 1, the Lorentz bridge, the E8 code, and the Suzuki W1 theorem (erratum-corrected one-scalar dictionary) — with the RH-hard step W3 explicitly open. Not RH evidence.",
  },
};

export default function PrimeFrontPage() {
  return (
    <>
      <section className="relative isolate overflow-hidden pt-12 pb-8 sm:pt-16">
        <div aria-hidden className="absolute inset-0 grid-bg pointer-events-none" />
        <div
          aria-hidden
          className="absolute -top-40 left-1/2 -z-10 h-[520px] w-[1000px] -translate-x-1/2 rounded-full opacity-35 blur-3xl"
          style={{
            background:
              "radial-gradient(closest-side, rgba(56,189,248,0.28), rgba(16,185,129,0.12), transparent)",
          }}
        />
        <div className="relative mx-auto max-w-5xl px-4 sm:px-6 lg:px-8">
          <nav aria-label="Breadcrumb" className="mb-6">
            <Link
              href="/"
              className="inline-flex items-center gap-1.5 text-sm text-slate-400 transition-colors hover:text-slate-200"
            >
              <ArrowLeft size={14} aria-hidden />
              Back to overview
            </Link>
          </nav>

          <p className="font-mono text-[11px] uppercase tracking-[0.2em] text-sky-300/90">
            Research diary · 243 agent runs · 5818 sandbox checks · suite 667
            scripts
          </p>
          <h1 className="mt-3 font-serif text-4xl font-semibold leading-tight text-slate-50 sm:text-5xl md:text-6xl">
            The Prime Front
          </h1>
          <p className="mt-4 max-w-2xl text-base leading-relaxed text-slate-300 sm:text-lg">
            The Prime Front is TFPT&apos;s number-theory line: a research
            diary of how the compiler&apos;s E₈ bookkeeping kept producing
            the classical machinery of primes. One identification is now a
            machine-verified theorem; the RH-hard step is honestly open.
            This page tells the whole story in ordinary language.
          </p>

          <StorySixty />

          <WhereWeAre />

          <div className="mt-6">
            <HonestyBanner />
          </div>

          <BlindPrimeCallout />

          <nav aria-label="On this page" className="mt-8 space-y-3">
            {PRIME_FRONT_SECTION_GROUPS.map((g) => (
              <div
                key={g.label}
                className="flex flex-wrap items-baseline gap-2"
              >
                <span className="w-full font-mono text-[10px] uppercase tracking-[0.18em] text-slate-500 sm:w-44 sm:shrink-0">
                  {g.label}
                </span>
                <span className="flex min-w-0 flex-wrap gap-2">
                  {g.sections.map((s) => (
                    <a
                      key={s.id}
                      href={`#${s.id}`}
                      className="rounded-full border border-slate-700/50 bg-slate-900/40 px-3 py-1 font-mono text-[11px] text-slate-400 transition hover:border-slate-500 hover:text-slate-200"
                    >
                      {s.label}
                    </a>
                  ))}
                </span>
              </div>
            ))}
          </nav>
        </div>
      </section>

      <BigPictureSection />

      {/* 01 — Suzuki W1 */}
      <DiarySection
        id="suzuki-w1"
        plain="The matrix this diary builds from primes turned out to be, exactly, the matrix of an operator the number theorist Suzuki had defined independently — and that identification is now a proved, machine-checked theorem."
        eyebrow="01 · The Suzuki identification · v630 / v631 / v640–v644 / v648"
        title="W1 closes as a theorem — after an honest erratum"
        badge="machine-verified"
        visual={<W1DictionaryMap />}
      >
        <p>
          The RH architecture preregistered in v624 (contract
          PRIME.WEIL.OPERATOR.01, citations web-verified: Suzuki
          arXiv:2606.09096 and 2607.24830) starts with W1: identify the TFPT
          window form with the Galerkin matrix of Suzuki&apos;s localized
          Weil operator. First contact (v630):{" "}
          <strong className="font-medium text-emerald-200">
            the atom layers are the same object, literally
          </strong>{" "}
          — positions log n, weights Λ(n)/√n, exact on all 40 atoms — while
          the smooth-layer comparison measured a non-scalar conversion: the
          preregistered residual, with data.
        </p>
        <p>
          Hours later, v631 resolved it: the residual{" "}
          <em>is the zeta pole term</em>, and the follow-up rounds made the
          dictionary sturdy — v640 closed the boundary cells symbolically,
          v641 froze the dictionary and ran it unchanged on three fresh
          windows (a preregistered kill test:{" "}
          <strong className="font-medium text-slate-200">portable</strong>),
          v642 lifted it to the full quadratic form at operator level.
        </p>
        <p>
          <strong className="font-medium text-amber-200">
            Erratum (2026-08-02, corrected the same day):
          </strong>{" "}
          that chain read Suzuki&apos;s eq. (1.3) with Lerch coefficient −1;
          the paper&apos;s own §2.2 data lock +1/4 (v643, check C0.1). All
          the chain&apos;s identities are correct identities of its kernel
          g̃ = g − (5/4)·Lerch, and every measured number transfers verbatim
          via the exact identity cgal(g̃) = −4·cgal(g) — only the labels
          change: Suzuki&apos;s own smooth layer is{" "}
          <span className="font-mono text-slate-200">+ρ</span> (not −4ρ),
          the dictionary is the{" "}
          <strong className="font-medium text-slate-200">
            single scalar +1/D on both layers
          </strong>{" "}
          (sign-compatible with positivity), and the origin constant
          vanishes, κ = 0 exactly.
        </p>
        <p>
          On the corrected reading, v643 proves the{" "}
          <strong className="font-medium text-emerald-200">
            measure-level W1 theorem
          </strong>
          : Suzuki&apos;s L²₀ mean-zero condition is automatic on the u-side
          (the projection lemma — the last named remainder closes), A_arch =
          −g″_smooth exactly at every lag (3.4e−52), and the full form
          equality holds at 1.28e−10 on the common odd sector. v644 starts
          W2 honestly (classical FEM density at rate; Rayleigh–Ritz monotone
          from above on nested spaces; λ_a = 0⁺ within ~1e−9, no sign
          statement; the Mosco remainder named). And v648 types the W3 tool
          diagnosis: the sign-uncertainty toolbox has a real 25-digit
          dictionary to the critical strip, but its mass lever dies at
          d = 1 — while the W3 surface itself is empirically positive on all
          67 complete windows (min λ_min = +8.26e−4).
        </p>
        <p className="text-slate-400">
          The honest map: W1 (theorem-closed) → W2 (started, not closed) →
          W3 (uniform positivity — the RH-hard step, open; the toolbox
          diagnosis closes one candidate route) → W4 (classical given
          W2+W3). Closing W1 does not move W3. No RH claim.
        </p>
      </DiarySection>

      {/* 02 — C = 1 */}
      <DiarySection
        id="uniform-constant"
        plain="One measured constant controls every window with complete data — and the only two exceptions turned out to be missing data, not broken mathematics."
        eyebrow="02 · The uniform constant · v618 / v619"
        title="C = 1, exception-free — and the two violators were the data's edge"
        badge="machine-verified"
      >
        <p>
          The equidistribution conjecture of the theory-open section asks for{" "}
          <span className="font-mono text-slate-200">
            |q_real/q_model| ≤ C·h⁻¹
          </span>{" "}
          uniformly. The measured constant is now frozen:{" "}
          <strong className="font-medium text-emerald-200">C = 1</strong>. On
          the declared surface (69 floor-passed windows, h = 142…1445) the
          model value keeps <em>one</em> sign on the whole ladder — no model
          zero crossing anywhere — and on every lock-sign window{" "}
          <span className="font-mono text-slate-200">eps·h ≤ 0.982</span>,
          with tertile medians 0.61 / 0.45 / 0.39 falling with depth (v618).
        </p>
        <p>
          Exactly two windows violated the bound, and both carried a{" "}
          <span className="font-mono text-slate-200">q_real</span> sign flip.
          v619 found the mechanism, and it is disarmingly concrete: a
          window&apos;s atom demand runs to{" "}
          <span className="font-mono text-slate-200">u ≤ 2α</span>, the
          prime-power data cap sits at{" "}
          <span className="font-mono text-slate-200">U_max = 12.899</span> —
          and the two flip windows are exactly the two whose demand exceeds
          the cap. Injecting the same truncation into healthy windows
          reproduces the flips in sign <em>and</em> magnitude at both scales.
        </p>
        <p>
          On the complete-comb surface — 67 windows — the C = 1 bound holds
          with{" "}
          <strong className="font-medium text-slate-200">
            zero exceptions
          </strong>
          . The “sign-flip windows” are retired as data-boundary artifacts:
          extending the surface needs more prime-power data, not new theory.
          Scrambled combs break the bound by four orders of magnitude — the
          constant is genuine arithmetic placement. No uniformity proof, no
          RH statement.
        </p>
      </DiarySection>

      {/* 03 — Hook */}
      <DiarySection
        id="hook"
        plain="While checking its own bookkeeping, the project kept finding classical prime-number objects; this page asks which of those finds are mechanism and which are coincidence."
        eyebrow="03 · Hook"
        title="What if the bookkeeping secretly speaks number theory?"
        badge="sandbox"
      >
        <p>
          Prime numbers look elementary: indivisible integers. Their global
          pattern is anything but. The Riemann Hypothesis asks for a precise
          spectral order behind that pattern — and this diary does{" "}
          <em>not</em> claim to approach that hypothesis.
        </p>
        <p>
          TFPT is a discrete compiler: two axioms build an E₈ lattice and read
          off Standard-Model structure. While exploring that lattice&apos;s
          shell census, the suite found classical modular objects — thetas,
          Hecke eigenvalues, Apéry congruences — sitting inside compiler-native
          counts. The question became: which of those links are mechanism, and
          which are beautiful coincidence?
        </p>
        <p className="text-slate-400">
          What follows is the arc from first surprise (Teil 11) through a
          four-stage kill of the “archimedean-from-seam” slogan, through Hecke
          from geometry, to the July&nbsp;25 reframe: a{" "}
          <Link
            href="/verification"
            className="text-emerald-300 underline decoration-emerald-400/30 underline-offset-2 hover:text-emerald-200"
          >
            relative-trace identity (v538)
          </Link>
          , a{" "}
          <Link
            href="/verification"
            className="text-emerald-300 underline decoration-emerald-400/30 underline-offset-2 hover:text-emerald-200"
          >
            Weil structure with two named obstructions (v539)
          </Link>
          , an{" "}
          <Link
            href="/verification"
            className="text-emerald-300 underline decoration-emerald-400/30 underline-offset-2 hover:text-emerald-200"
          >
            amplitude route with a positive linear carrier (v540)
          </Link>
          , a{" "}
          <Link
            href="/verification"
            className="text-emerald-300 underline decoration-emerald-400/30 underline-offset-2 hover:text-emerald-200"
          >
            matching-lemma and transport-ledger package with two named limits
            (v541)
          </Link>
          , and the consolidated stand below —{" "}
          <a
            href="#program-status"
            className="text-sky-300 underline decoration-sky-400/30 underline-offset-2 hover:text-sky-200"
          >
            what remains TFPT-specific is exactly one object: I5 — now
            geographically framed
          </a>
          .


        </p>
      </DiarySection>

      {/* 04 — Compiler */}
      <DiarySection
        id="compiler"
        plain="The whole theory starts from two fixed numbers, which force one specific eight-dimensional lattice — everything on this page is read off that lattice."
        eyebrow="04 · The compiler in one picture"
        title="Two axioms, one lattice completion"
        badge="sandbox"
        visual={<CompilerSchematic />}
      >
        <p>
          The discrete compiler starts with two numbers only: the seam constant{" "}
          <span className="font-mono text-slate-200">c₃ = 1/(8π)</span> and the
          carrier rank{" "}
          <span className="font-mono text-slate-200">g_car = 5</span>. From
          those, the theory forces a split{" "}
          <span className="font-mono text-slate-200">D₅ ⊕ A₃</span> completed by
          a four-element glue group{" "}
          <span className="font-mono text-slate-200">μ₄</span> to the unique even
          unimodular lattice in eight dimensions —{" "}
          <span className="font-mono text-slate-200">E₈</span>.
        </p>
        <p>
          Everything on this page is about what that lattice&apos;s point counts
          know — and what they do not know — about primes and L-functions.
          Classical theorems stay classical; the TFPT contribution is the
          in-suite mechanics that make those objects appear from frozen
          geometry.
        </p>
      </DiarySection>

      {/* 05 — Signed census */}
      <DiarySection
        id="census"
        plain="Counting lattice points with signs produces, unexpectedly, a classical formula tied to the Gaussian integers and to pi."
        eyebrow="05 · First discovery · Teil 11"
        title="The signed census: θ₃² · θ₄⁶ as a tensor factor"
        badge="sandbox"
        visual={<ShellCensus />}
      >
        <p>
          Colour every E₈ shell point by its μ₄ “glue class” (four colours).
          Ordinary counting recovers the classical Eisenstein series{" "}
          <span className="font-mono text-slate-200">
            1 + 240 Σ σ₃(n) qⁿ
          </span>
          . The surprise is the <em>signed</em> difference between opposite
          colours:
        </p>
        <p className="rounded-xl border border-slate-700/40 bg-slate-900/50 px-4 py-3 font-mono text-sm text-sky-200">
          Θ₀ − Θ₂ = θ₃(q)² · θ₄(q)⁶
        </p>
        <p>
          Here θ₃² is the theta series of the Gaussian integers ℤ[i] — the
          classical object whose L-value{" "}
          <span className="font-mono text-slate-200">L(1, χ₄) = π/4</span>{" "}
          produces π. It appears as a literal tensor factor of the
          compiler&apos;s signed glue census (classical Jacobi theta algebra;
          the probe content is the in-suite correlation).
        </p>
        <p>
          Three character channels sit on the same shells:{" "}
          <strong className="font-medium text-slate-200">total</strong> (all
          primes, ζ(s)ζ(s−3)),{" "}
          <strong className="font-medium text-slate-200">signed</strong> (entire
          L-series — the glue character kills the pole), and{" "}
          <strong className="font-medium text-slate-200">spinor</strong>{" "}
          (2-adic).
        </p>
      </DiarySection>

      {/* 06 — Bridges */}
      <DiarySection
        id="bridges"
        plain="The same counting reproduces a famous congruence from the proof that zeta(3) is irrational — and look-alike controls fail it."
        eyebrow="06 · Surprise bridges · Teil 12"
        title="The census “knows” the Apéry numbers"
        badge="sandbox"
        visual={<AperyBridge />}
      >
        <p>
          The cuspidal piece of the signed count is the weight-4 form{" "}
          <span className="font-mono text-slate-200">
            f₈ = η(2τ)⁴ η(4τ)⁴
          </span>{" "}
          — classically the Beukers / Ahlgren–Ono form tied to Apéry&apos;s
          proof that ζ(3) is irrational. For every odd prime p ≤ 97 the probe
          checks{" "}
          <span className="font-mono text-slate-200">
            A((p−1)/2) ≡ a_p mod p²
          </span>
          . Via Teil 11, the signed E₈ count at odd prime shells satisfies the
          same congruence. Placebos on nearby eta products fail.
        </p>
        <p className="text-slate-400">
          Beautiful, form-specific, and still sandbox: a correlation inside the
          suite, not a new proof of irrationality.
        </p>
      </DiarySection>

      {/* 07 — Kill chain */}
      <DiarySection
        id="kill-chain"
        plain="The project deliberately tests its own favourite explanations to destruction — and publishes the failures as first-class results."
        eyebrow="07 · Honesty as a method · Teile 14, 19–25"
        title="The kill chain — presented as a feature"
        badge="sandbox"
        visual={<KillChainTimeline />}
      >
        <p>
          An early slogan said the seam&apos;s measured angle 2π was “the
          self-dual temperature.” Teil 14{" "}
          <strong className="font-medium text-rose-200">deflated</strong> that:
          the steps parameter is compiler-specific; the angle is universal
          Bisognano–Wichmann / Unruh conversion.
        </p>
        <p>
          Then the whole “archimedean term from the seam” route was killed in
          four preregistered stages: mode density → interval cut → dictionary →
          scattering phase. Lesson, typed and kept: the seam is a discrete μ₄
          clock, not a hidden Gamma factor. The archimedean piece of the
          explicit formula is treated as a classical externum for recovery
          work.
        </p>
        <p className="text-slate-400">
          Killing your own favourite story on purpose is the method. Null
          results are first-class outcomes.
        </p>
      </DiarySection>

      {/* 08 — Predict */}
      <DiarySection
        id="predict"
        plain="The lattice can say exactly whether a number is prime and what kind of prime it is — but not where the next prime lies; that limit is stated, not hidden."
        eyebrow="08 · Can it predict primes? · Teil 21"
        title="Three honest channels — and one missing operator"
        badge="sandbox"
        visual={<PrimeDetector />}
      >
        <p>
          <strong className="font-medium text-slate-200">(a) Exact geometric
          primality.</strong>{" "}
          n&gt;1 is prime if and only if the E₈ shell at norm 2n has exactly{" "}
          <span className="font-mono text-slate-200">240(1+n³)</span> vectors —
          the classical σ₃ criterion, checked with zero false
          positives/negatives to 10⁴.
        </p>
        <p>
          <strong className="font-medium text-slate-200">(b) Per-prime
          properties.</strong>{" "}
          Glue characters predict arithmetic type: χ₄-fibre ⟺ p = a²+b²
          (100% for p&lt;1000 in the probe). The compiler says what a prime
          does, not where the next one sits.
        </p>
        <p>
          <strong className="font-medium text-slate-200">(c) Positional
          prediction</strong>{" "}
          needs the zero spectrum. Measured budget:{" "}
          <span className="font-mono text-slate-200">x_max ≈ 0.31 · T</span>.
          That Hilbert–Pólya operator does not exist in the suite — stated
          honestly, not as a near miss.
        </p>
      </DiarySection>

      {/* 09 — Hecke */}
      <DiarySection
        id="hecke"
        plain="The standard machinery of modular forms emerges from stepping between lattice neighbours — the first machine-verified module of this arc."
        eyebrow="09 · The mechanism · Teile 27–32 · v535"
        title="Hecke from geometry — first load-bearing result"
        badge="machine-verified"
        visual={<NeighborStepping />}
      >
        <p>
          Kneser p-neighbours — isotropic lines in E₈/pE₈ — carry the Hecke
          structure of the census. The count of lines is{" "}
          <span className="font-mono text-slate-200">σ₃(p) · #P³(𝔽_p)</span>{" "}
          (enumerated at p = 2, 3, 5, 7: 135 / 1120 / 19656 / 137600).
        </p>
        <p>
          The frozen marked neighbour-sum operator is an affine Hecke element{" "}
          <span className="font-mono text-slate-200">
            ν_p = a · Id + b · T_p
          </span>{" "}
          with <span className="font-mono text-slate-200">b = σ₃(p) + a_p</span>
          . Prime fingerprints fall out of geometry:{" "}
          <span className="font-mono text-emerald-300">
            a₃ = −4, a₅ = −2, a₇ = 24
          </span>
          . Census redundancy is purely 2-adic oldform structure (dim 7 = 5+2);
          recovery is newform projection.
        </p>
        <p>
          Promoted as{" "}
          <span className="font-mono text-emerald-300">
            verification/v535_hecke_from_geometry.py
          </span>{" "}
          (HECKE.GEOM.01, 25/25, AUDIT OK). Classical theorems (Kneser, Hecke,
          Atkin–Lehner, multiplicity one) are classical; the claim is the
          in-suite mechanics. No RH statement. Later joined by v536–v539 — see
          the July&nbsp;25 arc below.
        </p>
      </DiarySection>

      {/* 10 — Eichler */}
      <DiarySection
        id="eichler"
        plain="The lattice count splits into a smooth background plus an interference term, and the interference is exactly the square of a modular coefficient."
        eyebrow="10 · The Eichler layer · Teile 33, 36 · v536"
        title="Smooth background + coherent interference"
        badge="machine-verified"
        visual={<TwoMelodies />}
      >
        <p>
          Once the neighbour operator is frozen, the geometric count splits as
          an elementary Witt piece plus{" "}
          <span className="font-mono text-slate-200">exactly a_p²</span> — like
          a smooth melody with a coherent flicker on top. Two-sided
          confirmation (mod-p geometry on one side, eta-product on the other)
          holds at p ≤ 5; closed forms extend the identity to p ≤ 100.
        </p>
        <p className="text-slate-400">
          Promoted as v536 (Eichler trace layer). Together with v535 and the
          half-integral bridge v537, it becomes one finite relative-trace
          identity — v538.
        </p>
      </DiarySection>

      {/* 11 — Weight drop */}
      <DiarySection
        id="weight-drop"
        plain="The easy half of the connection to zeta closes with classical tools; the hard half — the actual zeros — is untouched, and the page says so."
        eyebrow="11 · Two-channel weight drop · Teile 35, 39"
        title="Abelian channel closed; cuspidal channel remains"
        badge="sandbox"
        visual={<CenterAtlas />}
      >
        <p>
          Rankin–Selberg translates{" "}
          <em>only</em> the abelian shadow into GL(1) products of{" "}
          <span className="font-mono text-slate-200">{"{1, χ₄}"}</span>. The
          centre atlas shows the ξ-line (centre 1/2) is reached exactly by
          weight ≤ 1 theta factors: Mellin(θ₃) → ζ(2s) and{" "}
          <span className="font-mono text-slate-200">
            ζ_ℚ(i) = ζ(s) L(s, χ₄)
          </span>
          . That abelian weight drop is factorisation + Mellin — typed closed.
        </p>
        <p>
          The cuspidal channel — where the{" "}
          <span className="font-mono text-slate-200">a_p</span> live — still
          sits at centre 2 and needs its own bridge. Possessing ζ as a{" "}
          <em>function</em> is not possessing its zeros as a{" "}
          <em>spectrum</em>.
        </p>
      </DiarySection>

      {/* 12 — Stage 4 */}
      <DiarySection
        id="stage-4"
        plain="A map of the terrain shows the finite machine has only a two-point spectrum — far too small for the operator RH would need."
        eyebrow="12 · Stage-4 map · Teil 40"
        title="Two-point spectrum; infinitely many still missing"
        badge="sandbox"
        visual={<SpectrumLadder />}
      >
        <p>
          The operator algebra on the census forms is commutative with a{" "}
          <em>two-point</em> Gelfand spectrum (the σ₃-system and the a_p-system,
          with oldform copies). A Hilbert–Pólya carrier would need an unbounded
          / non-commutative operator with infinitely many eigenvalues.
        </p>
        <p>
          Only two candidate classes remained inside the suite&apos;s early
          vocabulary: seam modular flow, and adelic Bost–Connes-style
          completion. Each has preregistered kills. The July&nbsp;25 reframe
          then changed the game: instead of hunting one infinite operator, the
          diary switched to a <em>family</em> plus a relative trace formula —
          see the next section.
        </p>
        <p className="rounded-xl border border-slate-700/40 bg-slate-900/40 px-4 py-3 text-sm text-slate-400">
          Verdict of the terrain map:{" "}
          <span className="font-mono text-slate-200">TERRAIN-MAPPED</span>.
          Cartography, not a proof attempt. Distance to RH — stated without
          theatre — remains large.
        </p>
      </DiarySection>

      {/* 13 — July 25 arc */}
      <DiarySection
        id="july-25-arc"
        plain="Instead of hunting one infinite operator, the diary switched to a family of values inside a classical trace-formula frame — and isolated exactly two named obstructions."
        eyebrow="13 · The July 25 arc · Teile 51–64 · v538 / v539"
        title="From finite machine to a Weil structure with two named obstructions"
        badge="sandbox"
        visual={<WeilArcMap />}
      >
        <p>
          The day&apos;s reframe: stop asking for one infinite operator. Ask for
          a <strong className="font-medium text-slate-200">family</strong> of
          central values plus a{" "}
          <strong className="font-medium text-slate-200">
            relative trace formula
          </strong>
          . Sandbox Teile&nbsp;51–64 built that map; two pieces graduated into
          the load-bearing suite.
        </p>

        <p>
          <strong className="font-medium text-slate-200">(i) The reframe.</strong>{" "}
          v535, v536 and v537 are three projections of{" "}
          <em>one</em> finite relative-trace identity — promoted as{" "}
          <span className="font-mono text-emerald-300">v538</span>. The
          half-integral object g is literally the quaternary lattice form{" "}
          <span className="font-mono text-slate-200">
            n = (x²+y²)/2 + 2z² + u² + 2w²
          </span>
          .
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            (ii) First non-collapsing infinite carrier.
          </strong>{" "}
          The Waldspurger family kernel K_D (classical Waldspurger periods) has
          rank that grows 8→192 as the discriminant window opens — the
          preregistered “rank saturates” kill did{" "}
          <em>not</em> fire. Its spectrum <em>is</em> the central-value family
          at the GL(2) centre&nbsp;2 — not the ξ-line. This is not RH evidence.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            (iii) Self-generating Hilbert space.
          </strong>{" "}
          The weight |d|<sup>−5/2</sup> is only the critical line of the family
          measure, not a canonical measure. The RTF forces |d|<sup>−1</sup>. At
          that weight, a positive pairing is cutoff-independent and builds its
          own space by GNS:{" "}
          <span className="font-mono text-slate-200">ℓ²(d, b²/|d|)</span>.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            (iv) Both infinities in one object.
          </strong>{" "}
          The double series Z(s,w) packs p-towers (closed Euler factors) and the
          d-family into one bilaterally verified identity (errors down to
          ~8×10<sup>−8</sup>; classical Goldfeld–Hoffstein multiple Dirichlet
          series named classical). Residual probes find only classical GL(1)
          shadows — no new ξ-sector.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            (v) Exact GL(1) core.
          </strong>{" "}
          Inside the positive form, the trivial Sato–Tate isotype is exactly{" "}
          <span className="font-mono text-slate-200">
            G₀ = ζ_p(w−3)² / ζ_p(2w−6)
          </span>
          . Fibre decomposition after character patterns kills twist-mix
          (~10<sup>−16</sup>).
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            (vi) Final linear relation — two obstructions.
          </strong>{" "}
          Machine-checked (and promoted as{" "}
          <span className="font-mono text-emerald-300">v539</span>):
        </p>
        <p className="rounded-xl border border-amber-400/35 bg-amber-500/10 px-4 py-3 font-mono text-sm text-amber-100">
          Q_fam = 2Q_ζ − 2Q_ζ(♭) + Arch + Corr
        </p>
        <p>
          The claim is{" "}
          <em>
            “Weil structure fully identified up to two explicitly isolated
            obstructions”
          </em>{" "}
          — the obstructions are the verified content, not a footnote. Not
          “almost RH.”
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            (vii)–(viii) Progression — see status.
          </strong>{" "}
          Sandbox Teil&nbsp;64 resolves Corr as the det/det₂ Jacobian; the
          categorical minus remains. Teile&nbsp;65–66 close the π-digit front;
          Teil&nbsp;67 takes a first amplitude-Dirac step. Full consolidated
          stand in the next section.
        </p>

        <p className="rounded-xl border border-slate-700/40 bg-slate-900/40 px-4 py-3 text-sm text-slate-400">
          Display markers: sandbox probes stay{" "}
          <span className="font-mono text-amber-200">[sandbox]</span>;
          v538–v541 are{" "}
          <span className="font-mono text-emerald-200">[machine-verified]</span>.
          Ledger fine types live in the verification suite only. No silent
          marker upgrades.
        </p>
      </DiarySection>

      {/* 14 — Program status (consolidated) */}
      <DiarySection
        id="program-status"
        plain="Everything TFPT-specific has been compressed into one remaining object, I5 — provably equivalent to what RH needs, and honestly not proved."
        eyebrow="14 · Where the program stands"
        title="Prime and Riemann Front: From Finite Hecke Structure to an Infinite Stabilized Trace Space"
        badge="sandbox"
        visual={
          <div className="space-y-3">
            <ProgramStatusCallout />
            <ModuleLadder />
          </div>
        }
      >
        <p>
          <strong className="font-medium text-slate-200">
            A single arithmetic machine.
          </strong>{" "}
          Three load-bearing modules are projections of{" "}
          <em>one</em> finite relative-trace identity{" "}
          <span className="font-mono text-emerald-300">[machine-verified · v538]</span>
          . Geometric side (E₈ lattice / glue census) and spectral side (modular
          coefficients / twisted L-values) are computed independently and agree
          exactly on the verified finite modules.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            A canonical infinite state space.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox]</span> A family
          over fundamental discriminants: rank grows without saturation
          (&gt;&nbsp;6000 active discriminants; Hecke integral and exact). The
          Hilbert space is not chosen by hand — a positive pairing produces the
          canonical GNS construction. Independent cutoffs converge
          sub-percent; the prediction chain closes at ~0.09%.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Packing the two infinities.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox]</span> The double
          series Z(s,w) packs two infinite directions — Euler towers over p, and
          the discriminant family — in the classical frame of multiple Dirichlet
          series / relative trace formulas (Goldfeld–Hoffstein named classical).
          Geometric = lattice counting; spectral = L-values × Euler factors.
          Fence: this does <em>not</em> transport to ξ; the system stays at
          weight&nbsp;4. RH concerns GL(1) at Re&nbsp;=&nbsp;1/2.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The Weil structure and its two obstructions.
          </strong>{" "}
          Compared with the Weil structure of the explicit formula (Weil 1952):{" "}
          <span className="font-mono text-emerald-300">
            [machine-verified · v539]
          </span>{" "}
          isolates <em>two</em> discrepancies — an exponential correction factor,
          and a categorical minus tied to the square nature of the family. Both
          isolations are the verified content of v539; neither is a footnote.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Exact determinant stabilization.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T64]</span>{" "}
          exp(−Σ p<sup>−u</sup>) = det(1−K)/det₂(1−K) (Hilbert–Carleman),
          verified ~1.5×10<sup>−16</sup>. In a consistent det₂ formulation the
          correction vanishes identically — it is the regularisation Jacobian of
          the GL(1) transition, not a new object. At u&nbsp;=&nbsp;1/2 the
          ordinary determinant does not exist while det₂ does, so stabilization
          is <em>necessary</em>, not optional: the stabilized transition is the
          only currently identified determinant convention that can reach the
          critical line. This resolves obstruction&nbsp;2 in sandbox; it does
          not rewrite the v539 isolation claim.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The remaining categorical minus.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox]</span> Survives
          every tested canonical projection. Source identified: the channel
          collects even prime powers, bound to the square form of the twist
          family (Waldspurger: L<sub>central</sub> ~ b(d)²) — the construction
          sees the <em>square</em>, not the linear amplitude. Krein signature
          confirms genuine indefiniteness; 32/32 sign characters do not remove
          it selectively. Not a normalisation error, incomplete character sum, or
          regularisation artefact — a categorical feature of the square level.
        </p>

        <p className="rounded-xl border border-violet-400/35 bg-violet-500/10 px-4 py-3 text-sm leading-relaxed text-violet-50 sm:text-base">
          At T86 — 91 probes in — and with seven promoted modules (v535–v541),
          the matching lemma is closed on ALL atom classes (window-certificate
          format,
          modulo proven classics only). What remains TFPT-specific is exactly
          ONE object: I5 in one-family form —{" "}
          <span className="font-mono text-violet-100">
            {"Q_cert + Δ₂ + A_fam − A_shift ≥ 0"}
          </span>{" "}
          — provably equivalent to Weil positivity ⟺ RH.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The remaining mathematical problem — absolute compression.
          </strong>{" "}
          Rest list:{" "}
          <em>(1)</em> I5 in one-family form, typed ⟺ Weil positivity ⟺ RH.
          Matching Lemma: closed on all atom classes (window 10⁶ proved ·
          signed pairing on non-coherent · λ-channel on coherent · 2-line
          exact){" "}
          <span className="font-mono text-amber-200">[sandbox · T86]</span>
          {" / "}
          <span className="font-mono text-emerald-200">
            [machine-verified · v541]
          </span>
          . I5 core = explicit Gaussian-mode family in the one-atom band
          (log&nbsp;2,&nbsp;0.925]; a<sub>neg</sub>&nbsp;=&nbsp;0.7486
          (recalibrated by T93); named target (T); finite-block zone extension
          certified only on an 8-dim window — structurally the wrong instrument
          for the full claim{" "}
          <span className="font-mono text-amber-200">[sandbox · T87–T93]</span>.
          Blind demo: 753/753 primes, zero errors, no division{" "}
          <span className="font-mono text-amber-200">[sandbox · T94]</span>.
          Induction skeleton: identities proved, zones 2–4 closed, zone-5 tip =
          equality problem, asymptotics = one named arithmetic bound (A){" "}
          <span className="font-mono text-amber-200">[sandbox · T99–T101]</span>.
          Milestone at T101: 2428/2428 sandbox checks (3139/3139 at T125).
          Fence: fits/extrapolations marked; geography locates; it does not
          attack; I5 remains ⟺ RH. This is not RH evidence.
        </p>

        <p>
          <strong className="font-medium text-slate-200">Current status.</strong>
        </p>
        <ul className="list-disc space-y-1.5 pl-5 text-sm text-slate-400">
          <li>
            <span className="text-rose-200">Not present:</span> RH proof,
            almost-RH, a zeros operator, a proof of I5.
          </li>
          <li>
            <span className="text-emerald-200">Present:</span> mechanism through
            v541{" "}
            <span className="font-mono text-emerald-200">
              [machine-verified]
            </span>
            ; Matching Lemma closed on all atom classes{" "}
            <span className="font-mono text-amber-200">[sandbox · T86]</span>;
            I5 geography complete{" "}
            <span className="font-mono text-amber-200">[sandbox · T87–T93]</span>
            ; a<sub>neg</sub> recalibrated to 0.7486; blind lattice demo{" "}
            <span className="font-mono text-amber-200">[sandbox · T94]</span>;
            induction skeleton: zones 2–4 closed; asymptotics = arithmetic bound
            (A){" "}
            <span className="font-mono text-amber-200">[sandbox · T99–T101]</span>
            ; I5 in one-family form — the single remaining TFPT-specific object.
          </li>
          <li>
            <span className="text-sky-200">Since then:</span> the induction
            sprint{" "}
            <span className="font-mono text-amber-200">
              [sandbox · T102–T125]
            </span>{" "}
            compressed that named bound from one matrix inequality to one sign
            plus one accounting convention, with certified steps to zone 155,921
            and 400/400 certified rungs, and the finale assembled the whole
            chain end to end on 52 zones — told in{" "}
            <a
              href="#compression"
              className="text-sky-300 underline decoration-sky-400/30 underline-offset-2 hover:text-sky-200"
            >
              sections 22–24
            </a>
            . Series complete: 3139/3139 sandbox checks at T125.
          </li>
        </ul>

        <p>
          <strong className="font-medium text-slate-200">On π — closed.</strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T65–T66]</span>{" "}
          Not the decimal digits of π are the key. Digit probes are closed:
          π digits at prime places are uniform; density-detrended
          cross-correlations null; full placebo battery (e, √2, log&nbsp;2,
          crypto, p±1) plus blind windows — no replication (PI-NULL 16/16;
          FOUR-LEVEL-NULL 23/23). The π spike is classical continued-fraction
          structure (355/113); compiler constant 1/(8π) is not an outlier
          (z&nbsp;=&nbsp;−0.90). Contrast kept: π-driven Cramér randomness
          reproduces prime <em>density</em> perfectly, but only true primes
          track Hardy–Littlewood pair correlation (corr&nbsp;+0.81 vs flat) —
          what makes primes primes is arithmetic, not randomness. The key now
          sits in one-family I5 — not in π digits.
        </p>
      </DiarySection>

      {/* 15 — Amplitude / transport route */}
      <DiarySection
        id="amplitude-route"
        plain="A promising route through squared quantities hits a measured wall, and the wall's exact size is named inside the verified claim."
        eyebrow="15 · The amplitude route · Teile 67–72 · v540"
        title="From the Dirac square root to a positive linear carrier — and the measured wall"
        badge="machine-verified"
        visual={
          <div className="space-y-3">
            <AmplitudeRouteCallout />
            <ConeCoverage />
          </div>
        }
      >
        <p>
          <strong className="font-medium text-slate-200">
            Square level closed — every squaring deletes.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T67–T69]</span>{" "}
          The amplitude Dirac D = [[0,V],[Vᵀ,0]] with D² = family kernel exists
          exactly and is Hecke-equivariant; signs of b(d) are a genuine
          metaplectic residue (52% mixed fibres — Kohnen depth, classical).
          Geometric polarisation b = N₊ − N₋ is exact; Θ = N₊+N₋ is a pure
          Siegel–Weil Eisenstein eigenform (σ₃ eigenvalues, null cusp); the
          family is the difference of two positive counting families
          (b² = Θ² − 4N₊N₋) — yet the minus is polarisation-invariant. By a
          Cauchy–Littlewood lemma, <em>every</em> coefficient bilinear form
          inherits even-k deletion (theorem-like across five channels); at the
          same time the minus is exactly the square-class double-counting of the
          towers (inclusion–exclusion bookkeeping). No full-weight carrier can
          live on the square level.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            A positive linear carrier — plus-only ζ-balance.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T70–T71]</span>{" "}
          The linear positive measure stands: Θ(d) = −48·L(−1,χ_d) exact (Cohen
          1975); full weights [1,1,1,1]; Weil balance{" "}
          <span className="font-mono text-slate-200">
            Q = Q_ζ(g₋) + Q_ζ(g₊)
          </span>{" "}
          is plus-only (~1e-15). The ζ(2s) factor appears only as the squarefree
          sieve of the carrier — bookkeeping, not weight deletion.           Functional
          equation exact:{" "}
          <span className="font-mono text-slate-200">
            {"Λ_Θ(s) = 8^{1−s}Λ_Θ†(5/2−s)"}
          </span>{" "}
          (rel ~1e-40; Fricke closed via Jacobi inversions). Plus survives
          reflection; the mirror family has a rigid sign law. The guaranteed
          cone is FE-selfdual; the Weil cone is not. Overlap 5/24; violations sit
          exactly at the first spectral node. Fence: this is Euler-region
          positivity (absolute convergence), not a central-line statement.
          Classical named: Cohen, Shintani, Siegel–Weil, Jacobi/Fricke, Weil
          1952.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The cone library saturates — one measurable gap.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T72]</span>{" "}
          Twisting absorbs the sign class n&nbsp;≡&nbsp;0,1&nbsp;mod&nbsp;4 (gap
          −26% mean, −90% max), but coverage saturates at 5/24. The pin h(0)&nbsp;&gt;&nbsp;0
          blocks every Weil element against every twist. Final compression: the
          entire residual distance to the Weil cone is the FE-covariant gap
          functional λ* on the atoms n&nbsp;≡&nbsp;6&nbsp;mod&nbsp;8 — no finite
          theta library can erase it (Farkas/LP certificates). Promoted as{" "}
          <span className="font-mono text-emerald-300">v540</span> with λ*
          named inside the claim. The doors that furnish this wall are the next
          section. This is not RH evidence.
        </p>
      </DiarySection>

      {/* 16 — Doors furnished */}
      <DiarySection
        id="doors-furnished"
        plain="The remaining gap gets doors: no-go theorems for what cannot work, a calculus for the gap functional, and one named inequality that would close it."
        eyebrow="16 · The doors get furnished · Teile 73–81"
        title="Two no-go theorems, a λ* calculus, a window proof — and one named inequality"
        badge="sandbox"
        visual={
          <div className="space-y-3">
            <DoorsFurnishedCallout />
            <DoorsPanel />
          </div>
        }
      >
        <p>
          <strong className="font-medium text-slate-200">
            Uniform cone route closed — hybrids live per direction.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T73]</span>{" "}
          Even in the continuum, the direction-uniform cone route is closed
          (sign-constancy lemma; window pin on [0,&nbsp;δ<sub>h</sub>)). But for
          every one of the 19 uncovered Weil directions an explicit,
          verified per-direction hybrid cone exists — the constructive residue
          is a per-test-function certificate machine (HYBRID-GAINS).
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Door A — vacuum structure, spectral no-go.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T74]</span> The
          spectral Dirac phase carries the metaplectic sign datum (~93%) but
          has a Hecke defect. L2 no-go: every ε-equivariant spectral functional
          is exactly re-signing invariant — the spectral world is provably
          sign-blind. The Dirac vacuum sits on the atom of minimal
          Waldspurger-normalised mass. Final Door&nbsp;A requirement list
          R1–R5: the sought polarisation must act coefficient-wise, be
          non-multiplicative in tower depth, non-bilinear, not a support
          reshuffle, and sign-seeing while Euler-compatible — a Krein /
          Gupta–Bleuler <em>quotient</em> in which the ♭-piece is a null/gauge
          sector, not a summand. Classical named: Dirac sea, Krein, Kohnen,
          McKean–Singer.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Door B — λ* gets its own calculus.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T75]</span> The
          gap functional now has closed form, its own functional equation (orbit
          invariant tanh(σ²ω²/2) under positive multipliers ⋊ dilations),
          critical width σ<sub>c</sub>&nbsp;=&nbsp;√2, and convexity (averaging
          lowers; no scale wedge — an earlier fp artefact honestly corrected).
          The target inequality is measured both sides: safety factor 273 at
          ω&nbsp;=&nbsp;1 falling to &lt;&nbsp;1 at ω&nbsp;=&nbsp;4.2. Two named
          open inequalities remain (hull positivity = transport wall;
          universal λ*-vs-A). Classical: Fejér / support functionals, Mellin /
          dilation semigroup.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Door C — a universal recipe, and its named lemma.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T76]</span> The
          hybrid recipe certifies 91/91 nontrivial adversarial Weil directions
          (100%); cost is polynomial / window-extensive (λ<sub>m</sub> ~
          m<sup>5/2</sup> = Eisenstein law). Lattice discreteness is the floor —
          δ<sub>h</sub>&nbsp;→&nbsp;0 does not break S1. Conjecture form with
          named core lemma: the{" "}
          <em>Matching Lemma</em> on the log lattice (a Diophantine divisor-sum
          problem). Implication architecture typed,{" "}
          <em>not claimed</em>: Matching Lemma ⇒ value-side representability ⇒
          [if the value→spectral transport held — the open wall] ⇒ Weil
          positivity. Any RH content would relocate into a universality proof
          plus the transport wall. This is not RH evidence.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Matching Lemma — classically shaped, then window-proved.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T77–T78]</span>{" "}
          T77: classically shaped (Gronwall/Robin; lemma not yet proven). T78
          WINDOW-PROVED (25/25): machine-proved on [4,&nbsp;10⁶] — exact-integer
          inequality chain, full enumeration over 939 870 clash atoms, 0
          violations, exact margin 0.082159; four structure laws at 0
          tolerance. Tail honestly open: Robin 1983 + constants miss by factor
          6.16 — missing residual ingredient named: a correlation lemma
          (thinning × cancellation). T80 RESERVE-PARTIAL: signed envelope
          character-exact; tail closed to ~10²³; last gap confined to the
          χ₋₄-coherent class. Classical: Gronwall 1913, Robin 1983,
          Alaoglu–Erdős, Pólya–Vinogradov, Landau.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Avoidance fails — theorem-shaped.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T81]</span>{" "}
          AVOIDANCE-FAILS (31/31): the same multiplicativity that makes the
          lever exact proves the counter-lever: coherent targets are reachable
          only by coherent rescalings — the freedom does not exist on coherent
          demand. Salvage: T76 recipe already minimally coherent (0 unforced
          keys, 100/100); forced coherent clash closed per certificate (83/90,
          worst ratio 0.18); on avoidant demand the class vanishes identically.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The transport ledger closes — the wall is I5.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T79]</span>{" "}
          LEDGER-CLOSES (22/22):{" "}
          <span className="font-mono text-slate-200">
            {"Q_Weil = Q_cert + Δ_arch + Δ₂"}
          </span>{" "}
          (identity ~7e-16 on 100/100); Δ<sub>pole</sub> ≡ Δ<sub>conv</sub> ≡ 0
          proved; the odd-prime side equals the certified combination exactly.
          The wall is one named inequality I5 (prime↔archimedean coupling),
          typed ⟺ Weil positivity ⟺ RH. Fence: equivalence typing, not progress
          toward proving it. This is not RH evidence. Classical: Weil 1952,
          Guinand, digamma terms.
        </p>

        <p className="rounded-xl border border-slate-700/40 bg-slate-900/40 px-4 py-3 text-sm text-slate-400">
          After T81 the constructive recipe chain saturated. Three new
          perspectives (next section) retype I5 and reopen the Z[i] channel —
          without claiming RH progress.
        </p>
      </DiarySection>

      {/* 17 — Three perspectives */}
      <DiarySection
        id="three-perspectives"
        plain="Three reframes land at once — a term thought external is internal, the wall cuts across the symmetry rather than along it, and the hardest class lives exactly where the compiler is most at home."
        eyebrow="17 · Three new perspectives · Teile 82–84"
        title="The arch term was internal, the wall is transversal, and the last class is the compiler's home"
        badge="sandbox"
        visual={
          <div className="space-y-3">
            <ThreePerspectivesCallout />
            <HeatFamilyBalance />
          </div>
        }
      >
        <p>
          <strong className="font-medium text-slate-200">
            The archimedean term was never external.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T82]</span>{" "}
          ARCH-INTERNAL (22/22): Δ<sub>arch</sub> is exactly the internal
          Γ-difference via Legendre duplication{" "}
          <span className="font-mono text-slate-200">
            {"((2π)^{−s}Γ(s) = ½Γ_R(s)Γ_R(s+1))"}
          </span>
          ; battery 18/18 rel 6.6e-15. The family carries its Γ factor as the
          Mellin signature of its heat-sum nature (verified also outside the
          convergence region). Bonus: the raw heat sum reproduces pole 5/2 and
          residue 8<sup>−3/2</sup> exactly from pure counting. Consequence: I5
          changes type. New form — for all autocorrelations h:
        </p>
        <p className="rounded-xl border border-amber-400/35 bg-amber-500/10 px-4 py-3 font-mono text-sm text-amber-100">
          {"Q_cert(h) + Δ₂(h) + A_fam(h) − A_shift(h) ≥ 0"}
        </p>
        <p>
          Self-consistency of <em>one</em> heat family: atom expansion vs Mellin
          signature of the same theta objects. Nearest classical relative:
          Connes 1999 semi-local trace-formula positivity (named context, not
          used). Fence: type change ≠ proof; I5 remains ⟺ RH. This is not RH
          evidence.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The wall is FE-transversal, not FE-positional.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T83]</span>{" "}
          INVARIANT-NULL (27/27): FE symmetrization is fully absorbed — 5/24 and
          λ* already are the numbers of the symmetric sector (the test region
          was right). Depth find: the product of the two FE reflections
          J₁∘J<sub>1/2</sub> = e<sup>±u</sup> is the unit-line shift — the
          centre delta is literally the transport operator. The value side is
          invariant under the whole infinite-dihedral ladder; the spectral cone
          only under its own reflection. Bonus: explicit-formula null test
          ≤&nbsp;2e-12 validates the whole convention. Classical: Mellin
          involutions, infinite dihedral group.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The last class is the compiler&apos;s home.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T84]</span>{" "}
          LIFT-WORKS-UNANCHORED (29/29): the coherent class equals primitive
          ℤ[i]-norms exactly; Grossencharacter phases replace Mertens divergence
          by L(1,λ)-convergence — the lifted chain never crosses; frontier jumps
          from ~10²³ to ~10<sup>5.9·10¹²</sup>. Circle closes: the last gap sat
          in the ℤ[i] sector — the origin object of the series (μ₄-glue, χ₄,
          θ₃²) — and exactly there the compiler&apos;s own character structure
          supplies the control that is provably impossible over ℚ. T85{" "}
          <span className="font-mono text-amber-200">
            [sandbox · LEMMA-CLOSES-LAMBDA]
          </span>{" "}
          closes the coherent class via the λ-channel (90/90 certificates;
          3.6× window margin). Fence: provably-shaped, not a formal proof; I5
          untouched. Classical: Hecke 1918/1920, Grossencharacter L-functions,
          Landau. This is not RH evidence.
        </p>
      </DiarySection>

      {/* 18 — I5 geography */}
      <DiarySection
        id="i5-geography"
        plain="The remaining object I5 sits precisely between two decades-old research programmes, which frame the same gap from opposite sides."
        eyebrow="18 · I5 geography · Teile 87–101"
        title="The I5 geography: two decades-old programs frame the same gap"
        badge="sandbox"
        visual={
          <div className="space-y-3">
            <I5GeographyCallout />
            <CrossingMap />
          </div>
        }
      >
        <p>
          <strong className="font-medium text-slate-200">
            The Connes dictionary is exact at the core.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T87]</span>{" "}
          DICTIONARY-EXACT-CORE (22/22): Q<sub>cert</sub> atoms equal Connes&apos;
          finite orbit terms (rel&nbsp;6e-16); A<sub>fam</sub> − A<sub>shift</sub>{" "}
          equals Connes&apos; archimedean W<sub>∞</sub> term including the
          principal-value constant (rel&nbsp;3e-23); the internal kernel is
          exactly the Riemann–Siegel phase derivative{" "}
          <span className="font-mono text-slate-200">{"k_ζ = 2θ'_RS"}</span>{" "}
          (rel&nbsp;0.0); the dihedral shift is the scaling group R<sup>*</sup>
          <sub>+</sub> (Weyl commutation exact). Classical: Connes 1999,
          Connes–Consani 2021, Weil, Guinand, Bombieri, Burnol, Yoshida.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Two positivity programs are complementary.
          </strong>{" "}
          Connes–Consani prove positivity <em>before</em> the primes (a
          prime-free Sonin window whose boundary sits exactly at the first
          prime atom u&nbsp;=&nbsp;log&nbsp;2). The compiler proves positivity{" "}
          <em>after</em> the primes (atom certificates). The sectors touch; they
          do not overlap. The I5 coupling lives in the crossing region outside
          both — around t*&nbsp;≈&nbsp;2π. Transferable-shaped: Sonin
          compression (same kernel, same group — verified). Not transferable:
          positivity itself (fence).
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The tight set is eight parametrized curves.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T88]</span>{" "}
          TIGHT-SET-PARAMETRIZED (27/27): band-limited zero plateau (support
          below log&nbsp;2, consistent with the T87 boundary), safe zones, and
          eight nearly vertical tight curves whose spacing follows the Γ-side
          density law (ratio 1.01&nbsp;±&nbsp;0.08 — smooth Γ-density
          description, zero-free; no spectral identification). Zero negatives
          on true autocorrelations (earlier T76 negatives were missing
          arch/p&nbsp;=&nbsp;2 bookkeeping, as the ledger predicted). Validation
          at machine precision (null test 4.3e-13, ledger 1.1e-14). The RH
          content of I5 concentrates on these low-dimensional curves.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The crossing, measured: no sharp wall at log&nbsp;2 — a thin band
            where atom and arch balance, and a 1–4-dimensional residual
            subspace that carries the minimum.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T89]</span>{" "}
          CROSSING-MAPPED (19/19): window compression of the Weil form is
          exactly Bombieri&apos;s classical object (Bombieri 2000 —
          unconditional positivity classically known precisely up to support
          width log&nbsp;2: Yoshida 1992, Bombieri 2000, Connes–Consani 2021).
          Proven zone reproduced (prime side identically zero below log&nbsp;2,
          λ<sub>min</sub>&nbsp;≥&nbsp;0 everywhere). The boundary is{" "}
          <em>not</em> sharp — margin falls smoothly (classical
          band-limitation), no collapse at log&nbsp;2; atom turn-on is
          (a−log&nbsp;2)³-soft. The residual subspace controlled by neither
          program grows softly 0&nbsp;→&nbsp;4 dimensions (of 32) across the
          crossing and carries the global margin minimum — the I5 core, for the
          first time, as a concrete finite-dimensional object per window width.
          The real attackable content is the atom↔arch balance in the thin band
          log&nbsp;2&nbsp;&lt;&nbsp;a&nbsp;≲&nbsp;1.0; the full-form minimum
          lives in the pole-coupled DC direction that Connes–Consani exclude by
          vanishing conditions; the pole-free CC sector keeps measured comfort
          to a&nbsp;≈&nbsp;0.75. Honest self-correction: BOUNDARY-SHARP retyped
          (“small at the boundary” ≠ “collapses at the boundary”). Fence:
          margins beyond the proven zone are measurements; no spectral
          identification. Classical: Weil 1952, Guinand, Yoshida 1992,
          Bombieri 2000/2003, Connes–Consani 2021/2023,
          Connes–Consani–Moscovici, Suzuki.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The residual, dissected: explicit Gaussian modes — and three
            directions under no control.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T90]</span>{" "}
          CORE-DISSECTED (17/17): the 1–4 residual vectors are n-stable across
          the discretisation ladder (angles ≤ 2.1°) and explicit — closed
          Gauss×cos / Gauss×sin fits with 99+% capture, the I5 core as a small
          family of concrete even/odd Gaussian-modulated modes. The coverage
          matrix (certificate extension / Connes–Consani pole vanishing / Sonin
          projection, against ten vectors along the a-ladder) leaves{" "}
          <em>three of ten</em> vectors controlled by no structure at all. That
          decides the core question: the residual is <em>not</em> merely pole
          coupling — pole projection clears only the a&nbsp;=&nbsp;0.75 window;
          from a&nbsp;≈&nbsp;0.85 genuine pole-free atom↔arch content remains,
          and at a&nbsp;=&nbsp;1.2 an odd atom-coupled mode appears that the
          even analyses could not see. Requirement line: an I5 idea must deliver
          positivity for an explicit family of Gaussian-modulated modes in the
          thin band — and it is not reducible to pole cleanup. Classical cited:
          Bombieri, Connes–Consani vanishing, Slepian.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The core, dissected: explicit Gaussian modes; and the band&apos;s
            law: the first prime rescues positivity where the archimedean
            margin ends — with a named, provable-shaped target inequality
            extending the classical zone.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T91]</span>{" "}
          BAND-PARTIAL (19/19, 3/4 closed): the band is the one-atom zone
          log&nbsp;2&nbsp;&lt;&nbsp;a&nbsp;≤&nbsp;0.9253 with inner edge
          a<sub>neg</sub>&nbsp;=&nbsp;0.7486{" "}
          <span className="font-mono text-sky-300">
            (recalibrated by T93; was 0.7410)
          </span>
          . Beyond a<sub>neg</sub> the prime-free margin changes sign and the
          prime atom becomes load-bearing — the first prime rescues positivity
          where the archimedean margin is exhausted (exactly the T89 balance
          point). Atom turn-on law exact (k&nbsp;=&nbsp;2m+1, Beta integrals);
          uncertainty constants decided: a·t<sub>rms</sub>&nbsp;→&nbsp;π
          (Wirtinger) and a·t<sub>cent</sub>&nbsp;→&nbsp;2Si(π)&nbsp;−&nbsp;4/π
          =&nbsp;2.4306. Band and tight curves are two orbit regions of the same
          functional with shared exact scale{" "}
          <span className="font-mono text-slate-200">{"∫k_ζ = 2θ_RS"}</span>.
          Named target (T), for a in the band and ‖f‖&nbsp;=&nbsp;1:{" "}
          <span className="font-mono text-slate-200">
            {"(P_pole + A_arch)(f) ≥ √2·log2·h_f(log2)"}
          </span>{" "}
          — provable-shaped as a zone extension beyond Bombieri&apos;s
          log&nbsp;2 (a self-standing classical target!); RH&nbsp;⇒&nbsp;(T),
          (T)&nbsp;⇏&nbsp;RH. Honest open: the super-exponential λ<sub>pf</sub>{" "}
          rate remains empirical. Classical: Wirtinger/Rayleigh, Beta integrals,
          Lambert-W, Si integral, Bombieri, θ<sub>RS</sub>.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Zone-extension attempt — and the self-check that found the
            checker&apos;s bug.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T92]</span>{" "}
          T-SKELETON (36/36): certified is Q&nbsp;≥&nbsp;6.7e-12&nbsp;&gt;&nbsp;0
          on the 8-dim window subspace over the whole scanned region (error
          certificate, 1057-point covering). The full extension does{" "}
          <em>not</em> stand — λ<sub>min</sub> collapses geometrically (factor
          ~11 per mode); the complement would need ~5000 modes: the finite-block
          route is structurally the wrong instrument. Pearl:{" "}
          <span className="font-mono text-slate-200">
            {"k(0) = −γ − 3log2 − π/2 − log π"}
          </span>{" "}
          exact.{" "}
          <span className="font-mono text-amber-200">[sandbox · T93]</span>{" "}
          MIXED (41/41): T92&apos;s calibration flags against the band map were
          resolved by a third independent implementation — the T89/T91 map
          survives; the only real bug was in the checker (constants imported
          untranslated, accidentally certifying a harder four-atom region); one
          constant precision-improved by ~1% (a<sub>neg</sub>).{" "}
          <em>
            The self-check found the checker&apos;s bug — and
            precision-improved one constant by 1%. The anchor discipline works
            in both directions.
          </em>
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The relay, measured: each prime rescues the direction the previous
            zone exhausts — and arrives before it is needed. The proof target
            shifts from fragile minima to robust counterfactuals.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T95]</span>{" "}
          T-CONTINUUM-NUMERIC (28/28): C1 fully proved — |h<sub>f</sub>
          (log&nbsp;2)|&nbsp;≤&nbsp;1/2 via disjoint support intervals;
          ‖S‖&nbsp;=&nbsp;1/2 exact with characterised eigenspace; atom-extremal
          directions satisfy the target with margin (“the directions that
          maximize the atom cost are provably safe”); continuum margin curve
          positive everywhere; extremizer is not the two-bump — binding
          mechanism is atom rescue. Lower bound open; missing instrument named.{" "}
          <span className="font-mono text-amber-200">[sandbox · T96]</span>{" "}
          EDGE-ARTIFACT + RELAY-CONFIRMED (21/21): the T95 “edge” at α* was a
          map artefact — value exactly reproduced, edge reading withdrawn
          (λ<sub>min</sub>&nbsp;&gt;&nbsp;0 on all [0.38,&nbsp;0.86]; margin
          collapses exponentially, λ&nbsp;∼&nbsp;exp(−49α)).{" "}
          <em>
            Second self-correction of the weekend, same anchor discipline.
          </em>{" "}
          Without the log&nbsp;3-atom, λ<sub>min</sub> crashes to −0.445; the
          loser is the anti-double-bump at distance log&nbsp;3 (alignment
          −0.99); rescue identity to 5e-15. Handover windows all positive:
          +0.025/+0.009/+0.011/+0.007 for the first four atoms. Strategy shift:
          margin problem, not edge problem; numerics exhausted past
          α&nbsp;≈&nbsp;0.55; the counterfactual is the proof target
          (O(0.1)-sizes). Classical: Paley–Wiener, Prolate, Galerkin/Richardson.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The induction takes shape: self-similarity puts the hypothesis
            inside the decomposition; the target is now one scalar inequality
            per zone — and the third self-correction of the weekend replaced a
            circular lemma by an exact identity.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T97]</span>{" "}
          ALIGNMENT-ONLY with certified half-step (105/105): alignment is sharp
          (sign alignment ⟺ coupling window nonempty, without exception); the
          t=0 killer loss on the anti-bump space is proved (
          <span className="font-mono text-slate-200">
            {"k_eff = (1−cos(tu))k(t)"}
          </span>
          , gain ×2–4.8); structure pearl: the E₀ block is literally the same
          form on the smaller window — “the induction hypothesis appears inside
          its own decomposition as self-similarity” (7e-14).{" "}
          <span className="font-mono text-amber-200">[sandbox · T98]</span>{" "}
          LAW-CONFIRMED-MECHANISM-OPEN (44/44): the conjectured one-vector lemma
          was circular (Douglas range inclusion — the law is forced by
          positivity itself); three T97 premises honestly refuted; replacement
          target:{" "}
          <span className="font-mono text-slate-200">
            {"D_k(α) ≤ μ_k/2"}
          </span>{" "}
          — exact scalar inequality, no constant, no vector; holds in all four
          zones, saturates at zone tips. Certificate upgrades: E₋ 43%→93% mean
          (whole zone in 3 of 4) via the probability-measure identity on the
          archimedean wings; E₊ certified for the first time. Skeleton: 8
          pieces proved, 2 certificates, 3 refuted, 3 open. Classical: Douglas
          1966, Schur, Slepian–Pollak–Landau.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Induction closeout (T99–T101).
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T99]</span>{" "}
          DECAY-LAW-FOUND (23/23): exact parity selection rule (
          <span className="font-mono text-slate-200">
            {"J₋Q₋₀J₀ = −Q₋₀"}
          </span>
          ) — the fragile near-null mode is excluded from the binding channel by
          symmetry; recursive inequality with only 1.01–1.20× loss; termination
          is arithmetic (240/240 in ≤&nbsp;4 steps to the classical zone).{" "}
          <span className="font-mono text-amber-200">[sandbox · T100]</span>{" "}
          REMAINDER-CLOSES-ZONES (27/27) — the 100th probe: closure 11/24 →
          24/24 (6/6 in every zone); the drift was a lattice artefact; one lever
          gained 1.7–69× (“the Bessel step threw the induction data away
          twice”); zones 2–4 fully closed; zone-5 tip typed as an equality
          problem (Fredholm shape, simple degeneration). Classical:
          Bessel/Parseval, Slepian, Schur test, Fredholm alternative.{" "}
          <span className="font-mono text-amber-200">[sandbox · T101]</span>{" "}
          CROWDING-TRENDS (31/31) — the fork across 16 zones (n&nbsp;=&nbsp;2..29):
          collapsed law{" "}
          <span className="font-mono text-slate-200">
            {"w_k = 0.0838·(atom gap)/μ_k"}
          </span>{" "}
          (fit) — “the handoff window is the atom spacing divided by the atom
          strength”; primitives flat; D<sub>k</sub>&nbsp;≤&nbsp;μ<sub>k</sub>/2
          holds 64/64 and never fails; only the closing instrument loses (
          <span className="font-mono text-slate-200">{"r ~ exp(−0.16k)"}</span>,
          fit/extrapolation). Core: “The crowding sits in the proof family, not
          in the mathematics it is trying to prove — the most hopeful version of
          the verdict.” Asymptotics would need (A) the arithmetic lower bound of
          the collapsed law [the localized hardness], (B) uniform relative
          margin, (C) a better bulk instrument, (D) a finite check.
        </p>

        <p className="rounded-xl border border-slate-700/40 bg-slate-900/40 px-4 py-3 text-sm text-slate-400">
          The relay induction, audited across 16 zones: the mathematics trends
          self-sustaining — flat primitives, a collapsed one-parameter law for
          the handoff, and a target inequality that never fails; what loses the
          race is the current instrument, and the hardness is localized in one
          arithmetic lower bound. Status: identities proved, zones 2–4 closed,
          zone-5 tip = equality problem, asymptotics = one named arithmetic
          bound (A). I5 remains ⟺ RH; the geography locates where any attack
          must work, it does not perform one. Milestone at T101: 2428/2428
          sandbox checks. All laws marked as fits/extrapolations. This is not RH
          evidence. What happens to that named arithmetic bound over the next
          twenty-three parts is{" "}
          <a
            href="#compression"
            className="text-sky-300 underline decoration-sky-400/30 underline-offset-2 hover:text-sky-200"
          >
            the compression story below
          </a>
          .
        </p>
      </DiarySection>

      {/* 19 — The mechanism (T102) */}
      <DiarySection
        id="mechanism"
        plain="The moment where positivity gets hard is not a mystery: it is an exact crossing of two finite quantities, manufactured by a known coupling."
        eyebrow="19 · The mechanism · Teil 102"
        title="The onset is manufactured — an anchored crossing, not a singularity"
        badge="sandbox"
        visual={<HandoverCrossing />}
      >
        <p>
          <strong className="font-medium text-slate-200">
            The mechanism is exact.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T102]</span>{" "}
          MECHANISM-IDENTIFIED (42/42): the k-th atom acts on the three
          induction blocks E₋/E₀/E₊ as exactly{" "}
          <span className="font-mono text-slate-200">
            diag(−1/2, 0, +1/2)
          </span>
          , so the atom strength μ_k enters the handoff exactly once,
          linearly. A two-sided sandwich over the Schur profile σ_k(δ)
          brackets the handover: 2w_k lands between the two crossings in
          16/16 zones (ratio 0.749…0.940), and the onset is{" "}
          <em>anchored</em> at δ_c = 2w_k (R² 0.968). Honest correction to
          T96: the essential-singularity reading is compatible but no longer
          singled out — the onset is a crossing of two finite quantities.
        </p>
        <p>
          <strong className="font-medium text-slate-200">
            The binding constraint flips.
          </strong>{" "}
          No concentration condition binds: the bare E₋ form is strongly
          positive (2.65…3.52 — that is 4–14× the atom line μ_k/2), and the
          classical ceilings (Cauchy–Schwarz, Landau–Pollak/prolate) are
          saturated near 97%. The onset is manufactured entirely by the{" "}
          <em>Schur dressing against E₀⊕E₊</em> — the coupling to the
          induction hypothesis — which takes 35.7%…97.3% of the bare
          eigenvalue. The window law is a pure μ-power:{" "}
          <span className="font-mono text-slate-200">
            w ~ μ^(−0.563 ± 0.098)
          </span>
          , q = −1.84 ± 0.37 (fits); no log, no g_k.
        </p>
        <p>
          <strong className="font-medium text-slate-200">
            C/g was a proxy — refuted three ways.
          </strong>{" "}
          The T101 decomposition with the atom gap g_k is triply negative:
          causally impossible (Q_(k−1) is blind to the next atom&apos;s
          position), statistically dispensable (a causal law fits better),
          and arithmetically only a ceiling (C_k ≤ g_k·μ_k exactly; the
          16-zone extrapolation violates the ceiling from k = 69 — checked
          over 18 120 prime-power atoms to n = 200 000). The hard core
          localizes to one scalar per zone: a lower bound on σ_k just above
          atom entry. That is a probe-level typing of where the hardness
          sits, not progress on it. Sandbox; not RH evidence. Classical:
          Schur complement, Cauchy–Schwarz, Landau–Pollak/prolate.
        </p>
      </DiarySection>

      {/* 20 — The instrument rebuilt (T101 → T103) */}
      <DiarySection
        id="instrument-race"
        plain="What looked like a mathematical failure was a blunt tool — rebuilding the tool closed most of the map without changing the mathematics."
        eyebrow="20 · The instrument rebuilt · Teile 101 → 103"
        title="The race, rerun: the slope halves and the map jumps to 44/64"
        badge="sandbox"
        visual={
          <div className="space-y-3">
            <InstrumentRace />
            <ClosureMapGrid />
          </div>
        }
      >
        <p>
          <strong className="font-medium text-slate-200">
            T101 lost the race by the instrument, not the math.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T101]</span>{" "}
          Across 16 zones the primitives are flat and the target inequality
          never fails — but the closing instrument&apos;s race quantity r_k
          decays as r ~ exp(−0.1622k) (fit, ± 0.0562): only zone 2 closes
          throughout, 7/64 cells on the zone × wing-fraction map.
        </p>
        <p>
          <strong className="font-medium text-slate-200">
            T103 rebuilds the tool.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T103]</span>{" "}
          INSTRUMENT-IMPROVED (29/29), pure tool-building at door C: the T101
          race curve is first reproduced exactly, then re-run with a
          θ-weighted band sum (certified weights, chain ρ ≤ b_band ≤ b_tail ≤
          b_t99 at 64/64 samples) and full m(Λ) exploitation. The demand
          Λ_ok stays bounded across all 16 zones —{" "}
          <span className="font-mono text-slate-200">0.771…3.640</span>{" "}
          instead of 2.3…376, a reduction of 3.0×…103.4× — at the honest
          price of explicit modes growing 2 → 232. The new race slope is{" "}
          <span className="font-mono text-slate-200">−0.0748 ± 0.0116</span>{" "}
          (fit, 2.2× flatter): r_k falls only 9.33 → 2.70 and{" "}
          <em>never leaves the spectrum</em>. The closure map jumps 7/64 →
          44/64 with one fixed, k-uniform instrument (Λ₀ = 3, r = 2 — no
          zone tuning).
        </p>
        <p>
          <strong className="font-medium text-slate-200">
            Measured verdict: the loss is in the wing, not the bulk.
          </strong>{" "}
          θ-weighting and finite rank are exhausted — the bulk is not
          low-rank (effective rank up to 0.579·dim E₋). What remains is the
          wing slack S = 1 − ρ: the pencil is nearly saturated and S falls
          0.2091 → 0.0392. Named next levers: a wing-adapted
          prolate/Slepian basis, or a Fredholm shape of the equality
          argument. All laws are fits; the tool-problem is half solved, the
          mathematics unchanged. Sandbox; not RH evidence.
        </p>
      </DiarySection>

      {/* 21 — Two doors converge (T102 + T103 → T104) */}
      <DiarySection
        id="two-doors"
        plain="Two independent attacks ended up pointing at the same single object, so the remaining hardness is one scalar per zone, not many."
        eyebrow="21 · Convergence · Teile 102–104"
        title="Two doors, one object: the wing near-null direction"
        badge="sandbox"
        visual={<TwoDoorsConvergence />}
      >
        <p>
          Two independent attacks — door A on the handoff law&apos;s lower
          bound (T102) and door C on the closing instrument (T103) — end
          their day pointing at the <em>same</em> object. T102 finds the
          onset manufactured by the Schur dressing against E₀⊕E₊; T103 finds
          the remaining instrument loss in the wing slack S = 1 − ρ. The
          dressing and the slack are one object seen from two sides: the
          wing near-null direction.
        </p>
        <p>
          The hard core is thereby localized to one scalar per zone — a
          lower bound on the Schur profile σ_k(δ_ref) just above atom entry,
          i.e. quantitative Weil positivity at the atom edge. Provable-shaped
          next to it: resolvent edge-regularity{" "}
          <span className="font-mono text-slate-200">
            σ_k(δ) ≥ σ_k(δ_ref)·(δ/δ_ref)^q_k
          </span>
          .
        </p>
        <p className="rounded-xl border border-slate-700/40 bg-slate-900/40 px-4 py-3 text-sm text-slate-400">
          Result — T104 (SCHUR.PROFILE.BOUND, two independent arms:{" "}
          <span className="font-mono text-slate-300">
            schur_profile_bound_probe.py
          </span>{" "}
          +{" "}
          <span className="font-mono text-slate-300">
            schur_profile_chain_probe.py
          </span>
          ) is CHAIN-PARTIAL: the naive margin route is dead, exact
          spectral-split chains close 16/16 with finite data, and the hard
          core moves to a bare_k lower bound plus the soft dressing scalar L.
          From here the diary runs twenty parts of pure compression on exactly
          that core —{" "}
          <a
            href="#compression"
            className="text-sky-300 underline decoration-sky-400/30 underline-offset-2 hover:text-sky-200"
          >
            the next three sections
          </a>
          . Sandbox; not RH evidence.
        </p>
      </DiarySection>

      {/* 22 — The compression (T105–T112) */}
      <DiarySection
        id="compression"
        plain="Twenty diary parts squeeze one big matrix inequality down to a single boundary value — and a supposed wall turns out to be the measuring grid itself."
        eyebrow="22 · The compression · Teile 105–112"
        title="From one matrix inequality down to one boundary value — then a wall that turned out to be a ruler"
        badge="sandbox"
      >
        <p>
          <strong className="font-medium text-slate-200">
            The bare bound gets a closed form; the avoidance law becomes a
            theorem.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T105]</span>{" "}
          ONE-OF-TWO (28/28): the T104 arm discrepancy dissolves — one currency,
          exact split bare = μ_k/2 + b0 — and three classical steps (Bessel;
          Legendre; Cauchy–Schwarz at the pole pair) collapse into a closed lower
          bound, positive 16/16 at 81.1–92.7% of the measured value, with no
          eigenvalue and no induction data as input. The avoidance law is
          upgraded to a theorem, and an exact parity superselection appears: two
          channels that never mix.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Parity halves the object.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T106]</span>{" "}
          DENSITY-MAPPED (32/32): the Weil pole splits exactly into a positive
          rank-1 lift in the even channel and a negative rank-1 pressure in the
          odd one. The even channel closes 16/16; the odd channel — the one with
          the <em>better</em> density — carries all remaining hardness. Two
          routes are honestly killed on the way (the density chain, invariant
          amplification). What is left is one Loewner statement on half the
          dimensions.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            One scalar, then one number.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T107]</span>{" "}
          SCALAR-TRACTABLE (30/30): the matrix statement becomes exactly one
          scalar ratio{" "}
          <span className="font-mono text-slate-200">r = κ/ε ≤ 1</span>, measured
          r = 0.005…0.18 — two orders of magnitude of room instead of five
          decimal places. The symbol route is structurally dead: the certified
          symbol bound is negative (−2.54…−0.81) where the truth is positive, so
          Grenander–Szegő cannot deliver here in principle.{" "}
          <span className="font-mono text-amber-200">[sandbox · T108]</span>{" "}
          EPSILON-IDENTITY (44/44): ε turns out to be exactly the square of the
          last Cholesky pivot — the classical Szegő–Levinson prediction error —
          so its positivity <em>coincides with</em> the induction&apos;s own
          positivity rather than being an extra demand. (R) drops to two scalars,
          and on eight zones what remains is literally one boundary value of one
          explicit vector.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Both scalars certified; the circle closes on the measured zones.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T109]</span>{" "}
          BOUNDARY-CERTIFIED (29/29): the mechanism is not decay but
          cancellation (the source sits <em>on</em> the boundary), so
          Combes–Thomas is refuted with the exact Green row and replaced by a
          residual certificate that carries the cancellation instead of bounding
          it away; ω is cracked unconditionally by a graded matrix cap. The chain
          closes 16/16 on exactly one strict-margin input, 10²–10⁶ weaker than
          the conclusion.{" "}
          <span className="font-mono text-amber-200">[sandbox · T110]</span>{" "}
          MARGIN-PROPAGATES (28/28): that input propagates — certified base case,
          15 certified handover steps, and the atom entry is structurally free
          (the new atom <em>raises</em> the floor on 15/15). Three sharp gaps
          stay: no reserve, no scalar step law, no k-uniformity.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Then the ladder is driven deep — and the wall is measured, not
            extrapolated.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T111]</span>{" "}
          CROSSING-CONFIRMED (23/23): 199 zones, 117 handovers to n = 521. The
          crossing sits at n* ≈ 462 — an upper bound; T110&apos;s n ≈ 170 was a
          fit artefact — and it splits into three separate walls: the margin wall
          n ≈ 462, a twin-prime ladder wall at 521→523 (purely arithmetic), and a
          requirement wall at n = 727. Decisive detail: the handover mechanism
          itself never fails, 117/117 at retention 1.000000. The chain tears at
          a ratio, never at a step.{" "}
          <span className="font-mono text-amber-200">[sandbox · T112]</span>{" "}
          SCALING-PARTIAL (20/20): rebuilt in a frame whose cell width follows
          the local prime gap, two of the three walls fall structurally — but the
          margin wall is frame-invariant at exponent −0.974. It is the substance
          of the requirement, not the geometry.
        </p>

        <p className="rounded-xl border border-slate-700/40 bg-slate-900/40 px-4 py-3 text-sm text-slate-400">
          Eight parts, one direction: every stage removed something and named
          what was left. The honest state at T112 is a wall that survives every
          reframing — and the next part asks the question that decides its
          status: what currency is it measured in? Sandbox; not RH evidence.
        </p>
      </DiarySection>

      {/* 23 — The certification sprint (T113–T119) */}
      <DiarySection
        id="certification-sprint"
        plain="The compressed claim is then certified step by step, down to extreme depths, with the computer checking every rung."
        eyebrow="23 · The certification sprint · Teile 113–119"
        title="The wall dissolves, the depth explodes, and the last inequality gets a textbook address"
        badge="sandbox"
        visual={<DepthReachCallout />}
      >
        <p>
          <strong className="font-medium text-slate-200">
            The currency question: the wall measures the ruler.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T113]</span>{" "}
          SUBSTANCE-CONFIRMED (27/27): the falling ratio carries the same
          exponent −1.168 ± 0.259 in <em>all five</em> currencies (raw, /λ_max,
          /trace density, /D, /D²), so it cannot be normalised away. But the
          substance is not the expected one. Under refinement there is{" "}
          <em>no plateau</em>: both eigenvalues carry the same power of the cell
          width (D^1.83 and D^1.76). The continuum window form has no gap at all
          — the quantity that was falling measures the{" "}
          <strong className="font-medium text-slate-200">discretization</strong>
          , not the spectrum. And the positive floor survives only as a
          cancellation of relative size ~1e-7, so norm perturbation theory is
          five orders of magnitude too coarse. Consequence: the T109 requirement
          chain was dividing by an artifact.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Remove the division, and the wall dissolves.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T114]</span>{" "}
          WALL-DISSOLVES (22/22): rebuilt without the margin division, the exact
          Schur complement (Albert 1969) certifies 27/27 ladder steps —{" "}
          <strong className="font-medium text-emerald-200">
            eleven of them beyond the old wall, up to n = 1331
          </strong>{" "}
          — and all seven zones where the T109 chain tore, including the wall
          zone n = 449 itself. The exact object is O(0.1) with no cancellation
          (λ_min(S) = 0.068–0.154, i.e. 42–67% of the block scale), while the
          same quantity via the norm bound is negative by a factor 2.4e5–9.6e7:
          an O(1) numerator divided by a 1e-6 artifact floor. Every norm route
          <em>had</em> to fail. From here on chains stop at the compute cap
          (h ≤ 1500), never at a step.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Compression breaks the cap: a certified step at zone 155,921.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T115]</span>{" "}
          TRANSPORT-BLOCKED (26/26), with a large compression gain. Transport
          between the ladder&apos;s non-nested grids certifies only mild
          refinement — and for a principled reason, not a weak bound: on nested
          ladders, where the transport error is exactly zero, the Schur floor
          itself falls like ρ^−1.7, so no bound can undo it. The two-scale
          compression, however, keeps the step bit-exactly margin-free and{" "}
          <strong className="font-medium text-emerald-200">
            certifies a step at n = 155,921
          </strong>{" "}
          (117× deeper than T114), on a fine lattice of 93,470 cells compressed
          to 1490; the longest chain runs 10 steps, 33 certified steps over four
          chains, every certificate 10⁵–10¹¹× above the numerical noise floor.
          The stopper is cost on 3 of 4 chains and a failing step on none. The
          remaining list is the shortest of the series so far: three points, only
          one of them an inequality.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The induction step is a boundary process — and the prime comb refuses
            to be compressed.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T116]</span>{" "}
          RICCATI-PARTIAL (33/33): the global pole rides exactly in a
          12×12+12+1 state (bordered elimination, no truncation anywhere), and
          the Riccati march runs 169,236 prepends to 1,354,088 cells — 903× the
          old cap — at flat cost, 76 µs per step. Declared honestly as a
          cost-geometry demonstration, not a Weil certificate. What breaks is
          unexpected: the full symbol does not decay at all, because every
          incoming cell couples back to every prime power in the window. The
          prime comb is the object that refuses compression. Bonus, and the hinge
          of the next three parts: the one remaining inequality acquires textbook
          shape — ε is the error of a piecewise-constant Galerkin method, and
          classical Aubin–Nitsche duality hits the measured exponent.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Theorem-shaped, and no rate lost.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T117]</span>{" "}
          THEOREM-SHAPED (23/23): the family is exactly nested, so ε is a
          Galerkin best-approximation error of <em>one</em> bilinear form and its
          monotonicity is a theorem rather than a fit; positivity becomes a
          non-membership statement. The direction trap is handled in the open —
          Céa and Aubin–Nitsche give <em>upper</em> bounds — and the two-level
          chain delivers a certified <em>lower</em> bound on 19/19 pairs
          (bound/ε ∈ [0.111, 0.185]) at rate θ&apos; = 1.74 against θ = 1.76: no
          power of D is lost. Self-correction: T116&apos;s factor-120 jumps were
          a sweep artifact — all 23 prime-power entries actually{" "}
          <em>raise</em> ε. What remains is three named analytic lemmas about one
          symbol, each with a classical address, two of them constants rather
          than rates.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Two of three lemmas stand — and the arithmetic half closes as a
            theorem.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T118]</span>{" "}
          TWO-OF-THREE (36/36): the first route is refuted with a reason (the
          exact two-grid symbol is a <em>harmonic</em> mean, which the comb dips
          make vacuous on 14/14 windows), and the rescue is a classical shift
          onto the oscillation Gram, whose symbol is the{" "}
          <em>arithmetic</em> mean of the same aliasing pair — the low-frequency
          negativity is suppressed quadratically instead of poisoning the
          statement. On a 15,680-point FFT lever the certified floor rises
          logarithmically and crosses zero on 3 of 5 zones: the failing windows
          were under-resolved, not obstructed. Saturation turns out to be an
          identity here, so its constant is computed (0.252–0.336) rather than
          assumed.{" "}
          <span className="font-mono text-amber-200">[sandbox · T119]</span>{" "}
          ARITHMETIC-DONE (27/27): the arithmetic half closes as a theorem —
          positivity of the symbol for every cell width below an{" "}
          <em>explicit</em>{" "}
          <span className="font-mono text-slate-200">
            D₀(α) = exp(−(Ξ(α) + B))
          </span>
          , with Ξ the prime-power atom count and B = −1.0474 universal (drift
          ≤ 3.1e-4). The energy route to the remaining statement is proven{" "}
          <em>empty</em> — a tautology — which is what makes that statement
          genuinely new content. And the sharpest identity of the run,{" "}
          <span className="font-mono text-slate-200">κ_end = 1/(1+R)</span>{" "}
          exactly (1.1e-16), reduces everything to one discrete Harnack
          inequality with a classical address.
        </p>

        <p className="rounded-xl border border-slate-700/40 bg-slate-900/40 px-4 py-3 text-sm text-slate-400">
          Worth stating plainly, because it is the load-bearing fence of the
          whole sprint: this chain is classical numerical analysis from end to
          end and{" "}
          <strong className="font-medium text-slate-200">
            contains no zeta input anywhere
          </strong>
          . As a conditional lemma the material is essentially complete — what is
          missing is a proof of the Harnack statement, not a discovery. Sandbox;
          not RH evidence.
        </p>
      </DiarySection>

      {/* 24 — The Harnack pair, the telescope, and the assembly (T120–T125) */}
      <DiarySection
        id="harnack-telescope"
        plain="The final assembly runs the whole chain end to end on 52 zones; what is still missing is uniformity — one sign plus one declared convention."
        eyebrow="24 · The Harnack pair, the telescope, the assembly · Teile 120–125"
        title="Why the last constant is one, why the coupling resists, the direction flip that made the recursion carry — and the finale that composed the whole chain"
        badge="sandbox"
        visual={<TelescopeRungs />}
      >
        <p>
          <strong className="font-medium text-slate-200">
            The Harnack core becomes proof-shaped — at the price of one more
            defect.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T120]</span>{" "}
          HARNACK-EXPLAINED (21/21): the mysterious constant R ≈ 1 is not a
          constant at all but a <em>symmetry</em>. The two increment families are
          the odd and the even half of <em>one</em> sequence, offset by a single
          fine cell — so given one sign on the corner cells their difference is a
          sum over disjoint neighbour pairs, and{" "}
          <span className="font-mono text-slate-200">|R − 1| ≤ 0.04745</span>{" "}
          follows <em>unconditionally</em>, dominating the measurement on
          3112/3112 sign-pure rows. Two negative results count as much: the
          per-cell version of the inequality is provably <em>false</em> (ratios
          up to 4.8e3 — only the summed form holds, which is exactly the form
          proved), and the discrete maximum principle is false too, so the sign
          route is closed and needs a corner-localized decay estimate instead.
          The honest cost: the defect count rises 3 → 4, because in the real
          frame the window grows with the zone, and the unconditional D₀
          criterion then covers exactly 3 of 1492 zones — below the first
          handover the ladder actually uses.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Against the real ladder: a section statement survives where a symbol
            statement dies.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T121]</span>{" "}
          WIDE-RESTRUCTURED (21/21): the chain never needed the symbol to be
          positive — it needs the finite <em>section</em> to be, at the
          frame&apos;s own resolution. Over 16 windows of the real ladder
          (n = 7…283,303) the section is positive 16/16 while the symbol infimum
          is negative on 8 of them; where it is negative the section eigenvalue
          sits on the positive side at 0.64–3.23× its size. The mechanism has a
          classical address (Christoffel functions): a polynomial of that degree
          cannot fit inside one comb dip. The measurement discipline is kept
          explicit — dense rows carry a Cholesky certificate, large rows are
          labelled measurements, because a Ritz value can refute positivity but
          never prove it. The net balance is then computed instead of feared:
          only a poly-log deficit (α^−1.57, uniform in the resolution) rather
          than collapse, decomposed exactly into two repairable steps — and one
          link of the chain is honestly refuted on 42/42 rows.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            Both repairs land: the refuted link was an identity all along.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T122]</span>{" "}
          NET-IMPROVED (20/20): the term that broke T121&apos;s link is not an
          error term but the reflection half of an{" "}
          <strong className="font-medium text-slate-200">isometry</strong> — the
          oscillation block is exactly a compression of one window Toeplitz form.
          From that plus a certified cell envelope and Parseval follows a
          certified band floor, non-vacuous 36/36 up to α = 6.31, where
          T121&apos;s budget tore at 3.45. The structural version of the
          Rayleigh step is <em>sharp</em> (slack 1.00–1.03, drift α^−0.002): no
          α-dependence passes through it any more. Certified deficit halves to
          α^−0.729 and is exactly uniform in the resolution; with the measured
          coupling the balance reads α^−0.113 against the chain&apos;s own
          ceiling α^−0.116 — statistically indistinguishable. The chain closes up
          to its own ceiling as soon as the coupling step stops being a worst
          case.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The coupling resists — structurally, and the reason is the result.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T123]</span>{" "}
          CBS-RESISTS (20/20): the band margin closes almost for free (that
          reduction costs under 1%) and the oscillation block gains{" "}
          <em>uniform</em> certified positivity on every row of the ladder. But
          every certified route to the coupling needs the coarse form from below,
          and there the numbers are brutal: λ_min of the coarse block is
          1.7e-5–2.9e-4, three to five orders below the certified envelope. The
          diagnosis is exact — the near-null direction of the window form is{" "}
          <em>smooth</em>, so it lives in the coarse space, and its eigenvalue
          comes from a cancellation <em>inside</em> the form that no pointwise
          symbol minorant can see. So the worst-case coupling step must stay, and
          the theorem now says why. The entire remaining gap (α^0.5) is that one
          coupling — and it is identified, not estimated: it is the same object
          as the lid the chain throws away, four to eight recursion levels deep.
          Verdict: the two-level argument cannot be tightened. It has to become a
          recursion.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The telescope carries — because the rung is a maximum.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T124]</span>{" "}
          TELESCOPE-CARRIES (28/28): the nested level chain is <em>one</em>{" "}
          window form on nested spaces (nesting exact to 2.4e-14), so the whole
          two-level system of the previous part is literally one rung of a
          ladder, every rung is the saturation identity, and the rungs telescope:
          their sum is the whole quantity (9.2e-10). Then the key of the part —
          the <strong className="font-medium text-slate-200">
            direction flip
          </strong>
          . Read as a residual in the inverse norm, the rung is a{" "}
          <em>maximum</em>, not a minimum: every test vector gives a{" "}
          <em>lower</em> bound, and the denominator wants the form from{" "}
          <em>above</em> — exactly where the certified machinery works. The
          resulting certified rung bound holds on 400/400 rungs with a drift
          seven times weaker than T123&apos;s. T123&apos;s obstruction is still
          there and is now irrelevant: the new bound never forms the object that
          was obstructed. The recursion closes in the right direction, coarse to
          fine, with base case pure semidefiniteness — and the balance moves
          +0.444 of the α^0.5 gap, about nine tenths of the way to the ceiling.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            What that leaves — stated as the diary states it.
          </strong>{" "}
          The coupling defect{" "}
          <em>collapses</em>: it is no longer a quantitative estimate but a sign
          that the coarse-to-fine induction already carries, and three quantities
          (the coarse minimum, its condition number, and the coupling constant)
          leave the chain entirely. One new defect is booked, honestly negative:
          the solution-free version of the rung fails, because the coupling term
          cancels the data instead of perturbing it — so the certified bound
          carries the solution along, on the same bookkeeping standard as the rest
          of the chain. The Harnack pair from T120/T121 is untouched.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            The finale: the chain composes, and the Harnack pair leaves the
            spine.
          </strong>{" "}
          <span className="font-mono text-amber-200">[sandbox · T125]</span>{" "}
          ASSEMBLY-GREEN (34/34): every certificate of the arc mounted on one
          ladder, end to end. All five stages — base Cholesky, telescope rungs, ε
          lower bound, the κ chain, the margin-free handover — complete on 52 of
          52 zones, and 30 of those zones form a run that{" "}
          <em>literally composes</em>: consecutive prime-power atoms on{" "}
          <em>one</em> common resolution, so the new window of a zone{" "}
          <em>is</em> the old window of the next one (an integer identity on
          29/29 pairs), the incoming atom&apos;s block is bit-exactly zero, and
          the output of every step is the input of the next (residual 2.9e-14).
          That is a way <em>around</em> the frame seam rather than through it: no
          two consecutive gap ratios are dyadic, so per-zone frames can never be
          refined onto one another — the run picks one frame for all of it and
          pays with its length. The load-bearing finding is a change of shape:
          the weakest stage of the whole mounting is the κ chain, but that chain
          is the second, independent route — the{" "}
          <strong className="font-medium text-slate-200">
            load-bearing spine is Cholesky-certified on all 52 zones
          </strong>
          , so the Harnack pair no longer carries anything. 430 completed
          Cholesky certificates, 440 identities, the certified margin a factor
          32–8.7e10 above the <em>declared</em> floating-point bound, and the
          base case an <em>equivalence</em> rather than a bridge.
        </p>

        <p>
          <strong className="font-medium text-slate-200">
            And then the accounting.
          </strong>{" "}
          Theorem V-final is printed with line-by-line attribution — 25 lines: 10
          identities, 9 Cholesky certificates, 3 window certificates (all three
          inside the route the spine does not use), and exactly{" "}
          <em>one</em> hypothesis, which is an accounting convention — plus a
          five-point statement of what is <em>not</em> claimed. Of the 31 links
          of the chain, 90.3% are an identity or a completed Cholesky; on the
          spine it is 96.2%, with zero window certificates. The series ledger:
          24 refuted routes across 18 parts, four of them killed by the{" "}
          <em>same</em> cause — a bound that needs the coarse form from below,
          which is why the exit was a direction flip and not a sharper estimate —
          and a seven-station cascade whose last station is new. What remains is
          named as what it is: <em>uniformity</em> in the zone index, not size.
        </p>

        <p className="rounded-xl border border-slate-700/40 bg-slate-900/40 px-4 py-3 text-sm text-slate-400">
          Where the series ends: a finite, machine-certified chain whose
          load-bearing spine is 96.2% identity or Cholesky certificate and whose
          only hypothesis is a declared accounting convention; a relay mechanism
          that never failed a step; and certified single steps as deep as zone
          155,921. What does not stand is uniformity in the zone index — and
          that, not any missing estimate, is the honest distance to any infinite
          statement. The distance to RH remains large. Series complete at 125
          parts and 3139/3139 sandbox checks. A second phase — the full proof
          (T126+) — has since opened on exactly that uniformity gap: the seam
          architecture is finished (T126 · SEAMS-CERTIFIED, 31/31), the two
          genuinely new inequalities it left are both proof-shaped and both
          changed under dissection (T127 · BOTH-SHAPED, 28/28), and the
          resulting four-point list is worked cheapest first (T128 ·
          THREE-OF-FOUR, 27/27): three of the four points stand at their
          preregistered bars — the exception list derived and closed, the
          retention bound exact bookkeeping, the boundary-layer exclusion now
          a proof — while the kappa bar was missed honestly, by 3.6%,
          systematically in the ratio. T129 (KAPPA-WILD, 28/28) then tested
          the preregistered kappa law and broke it productively: the fitted
          law falls once on 331 fresh transports, but the structure
          underneath is solved with identities — flat is exactly 1, linear is
          exactly 2, everything above is curvature, and the curvature chain
          is a per-transport theorem on all 436 — while the two affordable
          deep seams carry complete certificates on the graded space with a
          measured 8% false-positive rate declared before the results. T130
          (ONE-OF-TWO, 30/30) then attacked the two named pieces and exactly
          one stands: the graded-to-uniform bridge stands as an identity —
          the matrix-form Céa/Strang defect reproduces the uniform floor on
          84 pairs with zero overshoot, explains the 8% false positives
          completely (they sat in a bracket the bridge makes unnecessary),
          and carries both deep seams to positive fine floors at up to 3.8×
          the factorization cap — while the curvature bound honestly broke
          its frozen shape band on 13/545 disjoint test transports and is
          reduced to a uniform bound on one exponent; one new small gap is
          named (M25, a certified floor for the old window block) with an
          elegant candidate: the chain&apos;s own telescope supplies it. T131
          (SUPPLY-PARTIAL, 25/25) then tested exactly that candidate and
          built the self-supply loop, one number short of closed: the
          epsilon-to-floor secular sandwich is a theorem (sharp to ~1.3,
          sign half an equivalence via Albert) that replaces the Lanczos
          estimate on all 84 bridge pairs — exposing that the old Ritz
          value overestimated the floor by up to 7.9× — sign constancy is
          proved via Perron–Frobenius on the inverse (575/575), the
          one-hump honestly broke at depth, and M25 reduces to positivity
          of the pole-free section with nine decades of slack. Two
          irreducibles remain (the word &ldquo;for all&rdquo;, the RH
          address). The reverse-flow probes T132/T133 carried the certified
          toolkit back to the theory side; T134 (PARTIAL, 21/21) then closed
          the existence half of the pole-free floor — every Cholesky pivot
          positive on 79/79 windows, yet all six cheap lower-bound routes
          fail by sign, not size, leaving one named opening: an M-matrix
          question; and T135 (BOUNDED-STATE, 13/13) found a bounded faithful
          seam state (m_cert = 12) where the Weil window provably has none.
          T136 (ONE-CARRIES, 30/30) then closed one of that question&apos;s
          three items outright — Varga&apos;s regular splitting makes the
          Jacobi radius an identity and Collatz–Wielandt bounds it a priori to
          1.00–1.03 of the measured gap on 900/900 blocks — while the exact
          three-way bookkeeping puts the whole D-degradation in the margin and
          M17 closes negatively (the bad subspace is delocalised); T137
          (BOTH-RESIST, 22/22) made the perturbation explicit (comb-generated
          anti-diagonal stripes at the prime-power atoms, each a perfect
          matching, amplitude certified) and killed a whole family by
          certificate: ρ(|E|) ≥ 1.32 from below on 35/35 windows, so no
          absolute-value envelope can ever supply the floor — the residue is
          one named object, a sign-preserving bound. Both identity blocks are
          now load-bearing as v543 and v544. T138 (PAIR-EXACT, 26/26) kept
          the signs and found the mechanism: the sign of every inter-stripe
          coupling follows the interval geometry of its two edges (nested
          positive, disjoint mostly negative), the couplings cancel globally
          to 0.29–0.70, and the m-paired Neumann certificate — majorise only
          the outer factor of an exact identity — removes the arithmetic wall
          on all 77 dead blocks (pool 563 → 875/900), while the best
          certified signed bound cuts the overshoot by three orders of
          magnitude without reaching a floor. T139 (DENSE-RESISTS, 30/30)
          then asked the classics for the missing decay lemma and got an
          arithmetic refusal: the Jaffard hypothesis fails as certified (the
          archimedean lags sit on the borderline, the prime-power atoms lay
          O(1) spikes), the never-truncating layer series is killed from
          below (every layer series ≥ 1.03 &gt; target on 15/15), yet the
          sign law of T138 now FOLLOWS from one exact telescoping identity —
          derived, not measured — and the m-ladder certifies 619/620 border
          blocks. The measured core has shrunk from &ldquo;a decay law is
          missing&rdquo; to &ldquo;one signed inequality at stripe distance
          b ≤ 16 is missing&rdquo;. T140 (FINITE-CORE, 31/31) then attacked
          exactly that inequality and gave it an exact finite core per zone:
          the telescope identity lifts to the form level (Gram = CHCᵀ exactly,
          rank ≤ h−1), the spectral radius is an exact eigenvalue
          ρ(W) = λ_max(K^½HK^½) of a closed-geometry coverage kernel K times
          the mass-plus-Dirichlet form H, the checkerboard split replaces the
          O(nb) Weyl steps by three D-independent ones, and all the
          D-dependence sits in the geometry (blocks zone-uniform at D^0.13,
          λ_max(K) ~ D^−2.99) — while the additive band reading is closed
          negative from below at every b ≤ 16, with rank inflation as its
          mechanism. T141 (HARDY-RESISTS, 22/22) then attacked that Hardy
          ingredient itself: four exact identities put it in classical
          two-weight shape (the K⁺ Rayleigh form, the signed crossing kernel
          L_N = DᵀBD, Q = B₊1 identifying T140&apos;s Cauchy–Schwarz step as a
          Gershgorin step, and DKDᵀ = L_Δ with the endpoint edge mass as
          diagonal), but it does not close — and the way it fails is the
          finding: the certified constant is NOT zone-uniform (D^−0.366 ±
          0.036 against a bar of 0.25, spread 14.21× over a D range of
          11.21×) while the object it bounds IS (D^−0.229 ± 0.007), so the
          growth is manufactured by the diagonal profile; the additive shape
          is dead as a shape at its own exact Weyl floor (1.694–3.855× the
          target) and the joint shape fails at the normalisation alone
          (Ω = 20.71–2723.99, Gantmacher–Krein returning). The rest list
          collapses to R1b: one closed conductance profile with Y ⪯ K⁺ and
          Ω ≈ 1. The identity blocks of both parts are now load-bearing as
          v545. T142 (PROFILE-RESISTS, 24/24) then constructed that profile
          instead of guessing it: the capacity decomposition
          K⁻¹ = DᵀJ⁻¹D + xxᵀ/cap (verified to 5.5e-12 on 26/26 windows)
          exhibits the optimal Hardy weight exactly — its Dirichlet half is
          an orthogonal projection, so Ω = 1 EXACTLY where T141 had guessed
          20.7–2724 — and the weight was never a free choice but geometry
          (Miclo, Maz&apos;ya). The certified chain now reads
          Λ·Ω = 2.2671–2.4536 against a target of ~1 (0/26, flat in D at
          D^0.006 ± 0.004): the failure is a constant factor ~2.3, not a
          power — and the rank ladder closes the whole comparison path (the
          chain stalls at 1.39–1.41 for every rank up to 128 and first
          crosses 1 at r*/m = 0.995–0.998, the tautology Y = K⁻¹ itself).
          Structural consequence: since ρ(W) = 1 − Θ(D³) with equality only
          at Y ∝ K⁻¹, no comparison argument can deliver D-uniformity — the
          next move is the sharp capacity-Rayleigh route (R1c). T143
          (SHARP-CARRIES, 24/24) then ran exactly that route, and it
          carries: the exact capacity-Rayleigh form
          1 − ρ(W) = inf_v [dᵀ(J⁻¹−B)d + (xᵀv)²/cap − Σₖ sₖvₖ²] /
          [dᵀJ⁻¹d + (xᵀv)²/cap] is an identity on all 26 windows, with a
          structurally new bookkeeping — the minimiser is orthogonal to the
          equilibrium charge (the capacity rank one carries nothing), the
          mass share is NEGATIVE (the site masses help), and the gap is the
          cancellation between the Green share (exactly 1) and the crossing
          share (1.0029–1.0698), so the naive split ρ ≤ ρ_mass + ρ_long is
          certified dead at 6–67× ρ. Maz&apos;ya&apos;s capacity criterion,
          applied to the GAP form, lands inside its window [1/4, 1] on all
          26 windows (Φ_sup·λ = 0.5438–0.6457) — for a form that is not
          Markovian — and the route&apos;s loss factor is zone-uniform
          (D^−0.048 ± 0.010 against the 0.25 bar). The supremum lives on
          closed families (0.9606–1.0000 of the best; in node coordinates
          INTERVALS dominate the extremal sets by 8.3–129.5 — the
          one-dimensional structure Muckenhoupt&apos;s 1972 closed form
          needs), and on 6 border blocks it is certified by enumerating
          every subset. Miclo&apos;s constructive chain loses a factor
          46–2561: the criterion&apos;s conclusion holds while its classical
          proof mechanism does not. D-uniformity is thereby reduced to ONE
          named inequality — cap_E(A) ≥ |A|·λ₀/c₀ with an absolute c₀, a
          non-Markovian Maz&apos;ya capacity bound. T144 (INTERVAL-CARRIES,
          31/31) then ran the interval route at exactly that inequality, and
          it carries: the interval class is exhausted exactly on every
          window (11,390,676 intervals via a Cholesky prefix-sum identity,
          verified to 1e-11), the closed two-weight sup B_res·λ̂ =
          0.6694–0.7813 lands inside the Maz&apos;ya window on the whole
          surface and flat in D (D^0.013 ± 0.005), the family restriction
          falls entirely — a max-density-subgraph bound (Charikar 2000;
          Goldberg 1984) covers all 2^m sets, Ψ_all·λ̂ = 0.7399–0.8515, the
          flattest number of the probe — and the Markov perturbation route
          is certified dead (the positive couplings carry 36.7–64.5% of the
          off-diagonal mass). The certified chain
          λ ≥ 1/(c₀·κ_up·c_glob·B_res) delivers 0.1002–0.2653 of the exact
          gap per window with exactly ONE unproven input — the absolute
          Maz&apos;ya constant c₀, whose sharpest shape S1′ is a
          Muckenhoupt-type hypothesis. T145 (ONE-STEP-MISSING, 33/33) then
          ran the proof attempt itself, transcribing Maz&apos;ya&apos;s
          classical proof step by step onto the non-Markovian form: M1, M2,
          M3 are theorems (verified in the direction used), the Markov
          property sits in exactly ONE line — M4, truncation subadditivity —
          and that line is void here (3.2–14.9% positive off-diagonals
          carrying 18.3–58.7% of the mass; the conductance share of the
          minimiser&apos;s energy is negative on 61/64 windows, while the
          same code on the Stieltjes surrogate E−P reproduces the classical
          step). But M4 splits, and the mass half is a theorem that
          dominates: σ_tot = 0.2145–0.4425 &lt; 1 on every window, so the
          proof survives with the Markov theorem replaced by one exact
          number — c₀ = 12σ_tot = 2.573–5.310 on the energy route,
          c₀ = G_dy = 2.156–6.394 on the sign-free Green route, best
          c₀ = 2.248–4.227, flat in D (D^0.028 ± 0.017), and S1′ is
          CERTIFIED per window as cert_λ_max(R) ≤ c₀·Ψ on 64/64. An
          explicit no-go (R = aaᵀ+εI, entrywise nonnegative, proved density
          ceiling 4+ε, λ_max ~ log m) proves that a bound on the
          minimiser&apos;s LEVEL PROFILE is necessary — no weaker hypothesis
          can carry the inequality. S2′ is retired as a requirement, R5
          downgraded to a per-window certificate, all 24 border blocks
          closed; the mathematical rest is ONE named lemma — L1, an
          a-priori delocalization bound for the layer-cake constant. The
          identity spine of the capacity cycle T142–T145 is load-bearing as
          v546. T146 (LEVEL-CARRIES, 30/30) then ran the proof attempt for
          exactly L1, and the level lemma STANDS: closed on the measurement
          surface as a chain of theorems and certified window inequalities
          with no step reading the minimiser. Four hand-written closed
          profiles lose against the true minimiser by 5–9 orders of
          magnitude (T143&apos;s cancellation anatomy from the other side);
          what survives is the Green column itself, and the proof lever is
          the resolvent identity ψ = λRψ — the Θ(D³) smallness, everywhere
          else the obstacle, turns into the tool: ‖ψ‖_∞ ≤ λ_up·max_j‖Re_j‖
          on every window, where the trivial ceiling would give only
          Γ ≤ √m = 5.1–34.2 and the measured Γ = 1.8356–2.2797 sits a
          factor 2.7–16.0 below it — that factor IS the delocalization,
          read off the form. Davis–Kahan is instrumented and DISCARDED
          (informative on 0/64 windows: the spectrum bottom is a
          near-degenerate block, λ₂/λ̂ = 1.25–3.45, with the separation in
          the denominator), the cake base proves FREE (Maz&apos;ya&apos;s
          classical dyadic 8 falls to 2β² ≈ 2, gain 3.76 at zero cost), and
          the composite constant c₀^ap = 2β²·Γ·min(1,Γ₁) + ε =
          3.9042–4.8488 on 64/64 windows — against T145&apos;s measured
          2.156–6.394, and remarkably at the size of Maz&apos;ya&apos;s
          classical Dirichlet value 4 — closes the chain
          λ ≥ 1/(κ_up·c₀^ap·Ψ) at loss factor 0.0422–0.1586. The stress
          tests discriminate correctly (the no-go grows unboundedly,
          m^0.420, while the Dirichlet control stays flat), and the one
          genuine remainder is D-uniformity for ALL D — the asymptotic
          delocalization of the Green columns. The a-priori core of T146
          is load-bearing as v547. T147 (ASYMPTOTIC-SHAPED, 39/39) then
          attacked exactly that asymptotics, and the last open object
          became an exact identity: Γ = √Q★ · Sw, verified to 2.3e-16 on
          all 85 windows (m = 26–1491) — all-D uniformity is EQUIVALENT to
          the boundedness of two named scalars, the purely spectral
          Sw = λ_up·‖R‖_F (effective bottom multiplicity n_eff =
          1.251–2.831, certified ≤ 4.6438 by LDLᵀ inertia at a certified
          cut) and the purely geometric, scale-free Q★ =
          m·max_j(R²)_jj/tr(R²) (certified ≤ 2.8634 against the localised
          extreme m = 1491). The classical decay route is computationally
          dead (Demko–Moss–Smith q = 0.994815–0.999999, the envelope above
          the bar on every window; the columns delocalised, half-mass
          radius 0.048–0.105 of the window) — delocalization itself IS the
          bound, and T146&apos;s unexplained factor 2.7–16.0 is identified
          as √(m/Q_B) up to λ_up/λ_lo. The mechanism is named: the form is
          a Toeplitz-minus-Hankel section whose bottom eigenvectors are
          Fourier modes (Szegő/Widom; Kac–Murdock–Szegő 1953, Widom 1958,
          Basor–Ehrhardt 2009), and the sharp prediction Q_B ≤ 2|B| HOLDS
          (measured 1.375–1.839 × |B|), with the chain down to it pure
          arithmetic (sparsity ≤ 8.773 per bottom mode). The 3 open R4
          border blocks CLOSE (8/8 rebuilt border windows, c₀^ap =
          2.9135–3.0383, same mechanism), σ_tot does not fall but acquires
          an address (the truncation leaves the invariant subspace), and
          the mandatory stress separates the forms: the no-go diverges and
          hits the closed prediction Q_B = m/H_m to 1.0000 while the
          Dirichlet control stays flat. The honest limit: the
          mechanism&apos;s hypothesis is translation covariance, which the
          Jacobi whitening breaks (Toeplitz defect 0.21–0.54), and the
          spectral factor needs its own a-priori input — one statement, a
          Szegő theorem for the diagonally reweighted section, lifts the
          whole chain to all D. T148 (szego_bottom_probe.py,
          ONE-INPUT-MISSING) then dissected exactly that statement — Sw
          closes on the surface, the lifting hypothesis is isolated to one
          named scalar, TV(log Λ) — and T149 (weight_smoothing_probe.py,
          PARTIAL-SMOOTHING) eliminated that scalar exactly through the
          gauge freedom, refuting the hypothesis and relocating the missing
          input to the deep modes of the pure Toeplitz-minus-Hankel
          section; T150 (mode_ladder_probe.py, ONE-TERM-MISSING) named the
          mechanism — parity — with one growing gate left, and T151
          (odd_ladder_probe.py, ODD-CARRIES) closed that bound by a Sobolev
          reroute — linear per-mode, non-growing constant — leaving one
          scalar (the bottom pencil ratio) as the fit; T152
          (pencil_ratio_probe.py, ONE-TERM-MISSING) then refuted both
          archimedean gifts — positivity is a cancellation between geometry
          and arithmetic, not a property of either half — while certifying
          the floor via a Schur two-block criterion with a fixed 16-mode
          block, leaving one term (an m-free ceiling on K_bot); T153
          (psi_ladder_probe.py, REBUILD-RESISTS) then carried out the mapped
          rebuild and refuted it by exactly its own factor — Psi is pinned
          to within 1.20–1.30× by the first parity sine (a₁ = 8/π²), so the
          reserve never existed — while Psi closes as a term, the
          archimedean part turns out positive on the bulk parity block, and
          the two remaining terms lose to one Green/alignment estimate; T154
          (green_align_probe.py, ALIGN-RESISTS) closed the ceiling exactly at
          fixed size and left two uniformity floors, and T155
          (bottom_floor_probe.py, FLOORS-RESIST) made the complement floor an
          EXACT 12×12 certificate — recovering the collapse price 78.8–100%
          at fixed size — while proving that no symbol argument can ever
          deliver the block floor (the full symbol infimum is negative on
          every window); T156 (twelve_by_eight_probe.py, TERMS-RESIST) then
          reduced both remaining terms to single scalars — the t₁-loss is an
          exact closed function of two dimensionless numbers whose ceiling
          makes the first term one m-free lower bound on an angle (p₁ =
          0.2010–0.3282), and the block positivity is an alignment fact, not
          Fejér mass damping (arch and atom parts of nearly equal size cancel
          on the minimiser); and T157 (angle_floor_probe.py, ANGLES-RESIST)
          attacked both angles — neither falls, but both change shape: the
          sine-block confinement theorem proves the bottom eigenvector lives
          97.1–98.4% inside the first sixteen parity sines from the certified
          floor and ladder alone, the resolvent route pins 97.9% of the
          measured angle with one measured fixed-size scalar left, the first
          term becomes a bound on ONE diagonal entry of a 16×16 Schur
          complement the chain already builds (Cauchy–Schwarz misses it by
          5.4e2–5.2e5 — the cancellation is nearly complete), the arch half
          is uniformly certified via an executed adaptive Lipschitz ceiling,
          and the four instrument candidates of T155/T157 are promoted as
          v552; and T158 (schur_entry_probe.py, ENTRY-RESISTS) found the
          cancellation-seeing bound and proved it is a THEOREM — the Thomson
          dual form turns the entry into a Dirichlet maximum, so every trial
          vector bounds it from the right side (which is exactly why
          Cauchy–Schwarz missed by up to five orders: it evaluates a maximum
          at one direction), a Cholesky ladder of strictly positive terms
          pins the entry to 1/g₁₆ = 2.9670–7.9664, tight to 1.13–1.27 and
          flat, and the entry is exactly two-dimensional once one Green
          column is granted — while T157&apos;s growth pointer is refuted
          (the atom mass grows faster than the arch floor in every band;
          what carries the inequality is the off-band coupling a dyadic
          argument discards) and the measured-step count falls from three
          to one; T159 (exact_form_probe.py, FORM-RESISTS) then executed
          the cancellation algebraically — seven machine-checked
          identities including the gauge identity (the form is exactly
          blind to the lag mass, which is where the h²-sized halves live)
          and closed Dirichlet-kernel weights, the honest negative that
          the h² does not telescope away, and the sign law that the
          archimedean bulk block is a raw symmetric Z-matrix with a
          closed Collatz–Wielandt floor above target (theorem count
          6 → 13; cores of T158+T159 promoted as v553); and T160
          (pairing_probe.py, PAIRING-RESISTS) then attacked the pairing
          through its correlation structure and located the hardness
          exactly: W¹ does not oscillate at all (3 sign blocks, closed
          bound, the Leibniz device refuted by being worse than the
          pricing it was meant to beat), three closed sign-definite
          moment laws evaluate the smooth arch half to the
          double-precision floor (slack 1.95e-8, floor-limited — a
          second declared numerical horizon), and a machine-checked
          SAMPLING IDENTITY shows the atom half IS a Λ-weighted prime
          sum at 32 explicit frequencies, needed to depth h⁻² and
          measured to cancel only to 0.00–0.37 of the trivial bound —
          the h² cancellation is the intrinsic arithmetic hardness of
          the problem; and T161 (classical_closure_probe.py,
          CLOSURE-RESISTS) then ran the circularity triage, and the
          answer is the most consequential of the phase: the chain is
          NOT circular — the 32 frequencies are exactly the Fourier
          harmonics of the log-window (a theorem), the measured
          cancellation is fully the PNT main term, and the required
          depth has δ = 1.148–1.881 against RH strength 1/2, below the
          boundary term of every partial summation — the chain needs
          MORE than RH-strength input would supply, so the h² sits in
          the SPLIT, not in the primes; the analytic arch half closes
          m-free (scale-free Bernstein rate (3+√5)/2, closed head split
          Ψ, degree schedule K(h) = O(log h)), R-B is refuted in all
          four readings and replaced by a certified fraction bound
          ≤ 1/4, a re-split moves the demand to δ_eff = 0.98–1.38 —
          still above 1/2 — and the closed cores of T160+T161 are
          promoted as v554; and T162 (third_split_probe.py,
          DELTA-REDUCED) then delivered the third split and measured its
          end: a Mellin ladder of closed cell moments lowers the demand
          1.88 → 1.38 → 0.93 but saturates at its second term (an
          asymptotic series), one Abel step makes the demand prime-free
          with every arithmetic input in the Chebyshev constant κ, and
          the Fejér split pushes the proof demand BELOW the RH threshold
          1/2 on all 18 windows — at a price in the 1/s ceiling growing
          h^2.86, so the hardness is relocated into the price, not
          closed; R-A′ closes (a Lerch/Frullani integral at machine
          precision), R-B′ is refuted (the Gram form is indefinite), and
          the remaining question is a Pareto front; and T163
          (pareto_front_probe.py, FRONT-RESISTS) then surveyed exactly
          that front and answered it with a THEOREM rather than a failed
          search: the exchange law δ_bnd = 1/2 + log(2κg₁₆TV/P)/log X is
          an identity at all 27×50 grid points of the Fejér knob — price
          and demand are two coordinates of ONE object, and the crossing
          condition is exactly P &gt; 2κg₁₆TV — while a new four-line
          theorem (w₀ = ‖a‖²; TV ≥ |w₀| by telescoping; x₁ = 1 forces
          ‖a‖² ≥ 1/μᴾ₁) gives EVERY admissible trial vector a
          total-variation floor ~ h², verified on all 1674 trial vectors
          built (slack 5.0–23.2): a flat-price sub-1/2 demand is provably
          impossible in the parity sector — 0/27 windows cross at any
          flat cap, only the full chain-derived price P_aff ~ h^3.05
          reaches 27/27, and the crossing price h^1.91 sits 5–152×
          strictly inside the certificate the chain already accepts, with
          the margin widening. The hardness was never in the primes:
          T161&apos;s granularity, T162&apos;s h^2.86 and the crossing&apos;s
          h^1.91 are all the reciprocal smallest parity eigenvalue — the
          spectral gap of the parity Laplacian meeting the entry
          normalisation. R-C‴ is closed negatively by the inequality, the
          successor R-E is named with two prime-free arms (a growing 1/s
          ceiling for the downstream chain, or a sector whose gap does
          not vanish like h⁻²), and the theorem cores of T162+T163 are
          promoted as v555; and T164 (sector_change_probe.py,
          TOLERANCE-CARRIES / SECTOR-RESISTS) then decided both arms of
          exactly that successor. Arm A: the chain spends its entry
          ceiling at power exactly one — d log(1−F)/d log U =
          −1.0005…−1.0000 with d log r/d log U = +1.0000 exactly, so an
          h^ε ceiling costs h^−ε of margin for every ε &gt; 0, and the
          declared tolerance rule&apos;s ε* = 0.50 is disowned on the spot
          as the ε = 0 tailwind of this surface divided by power one (no
          free tolerance) — yet the O(1) gate turns out to be discharged
          window by window by a Cholesky identity: the single constant
          U_ref = 4.90 = max 1/g₁₆ carries every sub-surface window out
          of sample, 1/g₁₆ = 1.75–5.33 is flat (h^+0.061 over a 9× lever
          arm) and every ladder rung is strictly increasing (Schur 1917),
          so what T163 called R2″ collapses onto ONE quantifier:
          sup_m 1/g₁₆(m) &lt; ∞ — the m-freedom of a certified flat list,
          the same grammatical form as every other open uniformity item
          and no longer an existence question about trial vectors. Arm B,
          negative by a THEOREM: x₁ = 1 fixes only the scale, Q and TV
          are homogeneous of degree two, so Q/TV — and with it δ_bnd —
          is the same number in every sector (a gauge, checked to
          8.8e-10 on 5670 sector-by-vector combinations); the full space
          is strictly worse (μ₀ = 0 exactly, no finite floor); and every
          floor-flattening shift pays the identical exponent back as a
          transfer factor — floor × transfer = 1/μᴾ₁ identically, the
          exponents summing to +1.997 at every shift: the h² lives in
          the ratio, never in the floor. The one genuine surprise,
          immediately priced: an unconstrained ascent on Q/TV over the
          full window overshoots the crossing bar 2κ by 31–960× on 9/9
          windows (δ_bnd = −0.124…−0.090, alignment efficiency
          η = 7.9–11.1% of the Abel ceiling ‖C‖_∞ ~ h^+1.19) — and
          T163&apos;s TV floor HOLDS on the unconstrained optimiser
          (TV·μᴾ₁ = 8.30–11.72 ≥ 1, a family T163 never searched) at a
          price h^+3.30, worse than T163&apos;s crossing price h^+1.91:
          the difficulty is relocated, not reduced, and the binding axis
          is the gauge-invariant ALIGNMENT between Δw and the partial
          sums C. R-B‴ is narrowed (the −0.172 belongs to the T162
          quarter-bar object, under the 0.25 bar on 27/27, and mostly
          regresses on the zone scale, not on h), R-D is settled in type
          (the whole T156 kernel is gauge-invariant — no fifth device
          from a sector change), and the new successor R-F is named,
          gauge-invariant and prime-free. And T165 (alignment_eta_probe.py,
          ETA-RESISTS in the strong form) then closed exactly that
          successor — by certificate, not by failed search: the
          gauge-invariant identity P_pr = g₁₆·R·(TV/(t₁v)²)/μᴾ₁
          (machine-checked to 3.3e-16) makes demand and price ONE
          equation with four named factors, so with T163&apos;s floor
          every crossing vector pays ≥ 2κg₁₆/μᴾ₁ — the old KMS h² — and
          the predicates &apos;R &gt; 2κ&apos; and T163&apos;s crossing
          criterion agree on 105/105 vectors: R-F is strictly STRONGER
          than the R2″ demand, was never independent, and its m-uniform
          form IS R-E-A. The price exponent decomposes exactly (h^+3.261
          = +1.997 KMS + 1.279 overshoot + 0.093 floor − 0.108 g₁₆ — 61%
          is the same h² the series has met since T152), the η(Cap)
          ceiling is a theorem and falls short of the requirement by
          5.5–392× at the preregistered tolerance (0/27 ladder points
          reach the bar), and a genuinely decoupled ν-surface settles the
          old confound: the quarter-bar drift is ZONE DEPTH, not window
          length (factor 19) — at the honest cost of retiring
          U_ref = 4.90 (sup 1/g₁₆ = 5.7327 off-recipe). Exactly ONE
          genuine open object remains: inf_m g₁₆(m) &gt; 0, a lower bound
          on the 16-step Schur-cascade gain g₁₆/g₁ ≥ c·h^(3−ε) uniform in
          m — a quantifier over a certified flat list, a cancellation and
          hence provably beyond absolute-value budgets. The theorem cores
          of T164+T165 are promoted as v556. And T166
          (schur_cascade_probe.py, CASCADE-RESISTS, 30/30) then dissected
          exactly that lower bound, completely: four identities (all
          theorems, re-verified on a 63-window union of the frame-A and
          decoupled-ν surfaces) turn the ladder into one readable object —
          the whole gain is a NEAR-COLLINEARITY, g_K/g₁ = 1/(1 − R_K²);
          rung 2 alone carries a median 59% of it (K_half between 2 and 5
          on 63/63); and the gain is invariant under the spectral
          normalisation B → DBD: the h³ is a property of the arithmetic
          Gram block alone. The cancellation matches none of the prepared
          stories — no entry of the 2×2 block cancels and neither half is
          collinear alone, yet the sum degenerates to seven digits: it
          lives in the 2×2 GRAM DETERMINANT (pieces 2.6–630 against a
          full determinant 1.0e-3–1.7e-2), the same arch-against-atom
          mechanism as T159/T160, one level up, in a quadratic functional
          of the lag sums — and an anti-fitting scramble destroys the
          effect by a factor 4569, pinning it to where the prime powers
          sit. The best closed route reaches gain h^+1.319 against the
          target h^+3.110 — an exponent gap of h^+1.791, so the verdict
          is CASCADE-RESISTS: the one missing inequality IS the
          cancellation, now stated in three equivalent dresses (an m-free
          bound on one Gram-minor ratio, a closed β at rung 2 to 3.6e-4,
          or a closed near-null vector of A in the low modes). Two honest
          contract corrections are recorded — F2 is strictly weaker than
          inf g₁₆ &gt; 0 (B₁₁ grows h^+3.110, not h³), and the
          confinement route was a sub-surface artefact (h^+2.714 on a
          reduced eigen-horizon vs h^+1.963 on the full union) — and
          U_ref moves to 7.45 on the union (declared off-recipe,
          CERT-WINDOW). And T167 (null_vector_probe.py, VECTOR-RESISTS,
          39/39) then attacked the most constructive dress and delivered
          the strongest reduction of the series: the pivot identity
          uᵀQu = 1/g_K + δᵀQδ (for δ₁ = 0) makes every candidate&apos;s
          excess an exact number, and at K = 2 the closed vector
          (1, −Q₂₁/Q₂₂) is EXACT — a theorem, ρ = 0 on 63/63 windows,
          with B₁₁g₂ = 1/(1 − r₁₂²) as an identity and a certified gain
          exponent h^+2.921 on frame A (= h^(3−0.079)) — so the third
          dress collapses onto the second. Perturbation theory dies for
          a remarkable reason: the Jacobi series genuinely diverges, but
          the Kato series CONVERGES fast (radius 0.067) — to the WRONG
          object: the exact bottom eigenvector is itself a useless trial
          (ρ = 6.04, overlap with mode 1 only 0.083, relative spectral
          gap 5.4e-07). Against the working hypothesis the threshold is
          MILDEST at K = 2 (by 1.90× and h^+0.353) — exactly where the
          vector is free — while at K = 6 no closed family comes within
          264× (the curves diverge at h^+1.138). The scramble control
          sharpens to a change of TYPE: with scrambled prime positions
          the 2×2 diagonal itself loses positivity on 8/8 windows
          (Q₁₁ = −1.35e5…−1.02e4), so g₂ does not even exist there. And
          the unification is an IDENTITY, not a narrative:
          eps_ent(K) = ρ·(1/g_K)/S_K holds to machine precision —
          (1/g_K)/S_K is T166&apos;s determinant ratio (dress a), the
          1 − r₁₂² (dress b) and the entry threshold (dress c): there is
          ONE inequality, not three. What remains is R1, a single
          scalar: an m-free upper bound on 1 − Q₁₂²/(Q₁₁Q₂₂) ≤ C·h^(−3+ε)
          — three closed lag sums (required relative accuracy 1.25e-05
          median, 3.19e-08 at the worst window). The theorem cores of
          T166+T167 are promoted as v557. And T168
          (lagrange_minors_probe.py, MINORS-RESIST, 39/39) then ran the
          Lagrange identity on that scalar — and the sum of squares is
          REAL: the h×h arithmetic kernel A is positive definite on all
          63 windows (λ_min = 5.0e-6–8.5e-4, 8–12 orders above the
          solver floor), so the Lagrange identity is a genuine sum of
          strictly positive squares (ρ_LAG = 1.000 exactly) —
          hard-fenced as a statement about one finite matrix per window,
          the Weil criterion never tested, assumed or reverse-inferred.
          The anatomy is exact and surprising: the mode vectors are
          EUCLIDEAN-ORTHOGONAL (r₁₂ = 0 to 6.2e-17 — the
          near-parallelism is created entirely by the arithmetic
          metric), the Wronskian minors telescope in closed form and are
          MAXIMAL in norm (½ΣM² = 1/(μ₁μ₂) exactly ~ h⁴ — the minors
          carry zero arithmetic, all smallness sits in the PSD kernel at
          one closed vector), and in the eigenbasis the sum is THIN (90%
          from 4–101 pairs, top pair always (λ_min, λ_max), raw
          κ = 0.30–0.94). The exponent ledger closes with every factor
          but ONE at the right power, and that factor is identified
          exactly: ν₂ with minimiser t* = Q₁₂/Q₁₁, needed to 6.05e-3
          median / 2.82e-4 worst — the same accuracy T167 measured from
          the other end: by T168-TH7 that factor IS the target
          rewritten, the hardness is SELF-SIMILAR under exact
          reformulation. And T169 (tstar_ratio_probe.py, RATIO-RESISTS,
          41/41) then PROVED that self-similarity: every genuinely
          closed candidate for t* misses by O(1) and DIVERGES from the
          threshold at h^+1.99 (the block is atom-dominated,
          â₁₁^arch &lt; 0 on 63/63 — the archimedean half is not even
          PSD), while the only candidate that meets the threshold,
          √(â₂₂/â₁₁), does so by a new theorem (T169-TH7) because it
          reintroduces the determinant itself — R collapses in closed
          form to 2√(â₁₁â₂₂)·det Â/[(â₁₁+â₂₂)(√(â₁₁â₂₂)+â₁₂)]: the
          t*-language is provably exactly as hard as T167&apos;s
          determinant, THE LOOP CLOSES AS AN IDENTITY. The real gains:
          the first CERT-UNIF in weeks (ν₁ ≤ max(â₁₁,â₂₂)+|â₁₂|,
          Gershgorin, unconditional, losing only 1.10–1.19), the chain
          rebuilt R4-free so the Weil-shaped positivity of A_h never
          enters it (frame-A trend h^−2.948, ε = 0.052 inside the 0.5
          carry window), and the one open object now in STANDARD
          analytic shape for the first time — R1 is a BILINEAR VON
          MANGOLDT SUM against closed Dirichlet weights. That matters:
          T161&apos;s beyond-RH triage applied to LINEAR sums; for
          BILINEAR forms the large sieve gives unconditional square-root
          cancellation. And T170 (bilinear_sieve_probe.py,
          SIEVE-RESISTS, 40/40) then ran that toolbox — and the
          classification completed as a THEOREM: the form is written
          down exactly (det Â = det B − D(B,S) + det S, on the
          reference window a FIVE-ORDER cancellation between three
          ~200-sized pieces: 6.18, −198.8, 192.6 summing to 4.22e-3)
          and it collapses back onto the linear hardness for structural
          reasons — a RANK-3 THEOREM (the kernel is the polarisation of
          the determinant: rank K ≤ 3 for every window, every h, every
          X, so the bilinear form IS the rank-3 polynomial
          S₁₁S₂₂ − S₁₂² in THREE linear Λ-sums) and a TYPE-II COLLAPSE
          (the closed weights see n only through log n, so every
          Vaughan Type II block has effective rank O(1)). Gained: 32
          frequencies reduce to three linear functionals, by theorem.
          Not gained: the precision, which is what binds (5.3e-5 per
          sum at h = 285, sharpening h^−3 — an RH yardstick 1.1e3 too
          coarse). No unconditional route exceeds δ = +0.996 against
          the target 3.0 — the shortfall is a theorem, not a
          measurement — and the scramble moves the truth by 2.76 in the
          exponent while rank and bounds stay put: R1 is finally
          CLASSIFIED as a NEAR-DEGENERACY, not a size — two explicit
          finite Λ-sum vectors becoming collinear at rate h^−3, beyond
          any size-bounding tool. The theorem cores of T168+T169+T170
          are promoted as v558. And T171 (final_map_probe.py,
          MAP-COMPLETE, 43/43) then assembled the capstone — and the
          map is complete: all sixteen links of the reduction chain,
          from the I5 floor to the final R1 shape, reproduce in ONE
          connected run on 12 windows over two frames (thirteen as
          theorems, three as per-window certificates — with the honest
          exclusion said out loud: one deep window has an indefinite
          low block, so the ladder links carry 11/12 while the identity
          links carry 12/12); all eight classified no-go routes fail on
          instances exactly as classified (the size budget overshoots,
          the symbol infimum is negative, the sector change is an
          invariance, perturbation theory addresses the wrong object,
          the atom band outruns the archimedean band, the sieve
          certifies δ = +1.044 against a target of 3.0, the scramble
          flips the type, and the L_P gain is identically one); and the
          precision ledger closes — the needed joint precision
          (2.284e-7 at h = 1444, sharpening h^−3) stands five orders
          beyond the RH yardstick (1.2e5×) and three beyond the best
          unconditional bound (3.1e3×). The file adds ZERO new uniform
          statements — which is precisely why it is promotable: the
          capstone is load-bearing as v559 (PRIME.PHASE2.CAPSTONE.01).
          Phase 2 now stands as a certified map with exactly one open
          object: R1, classified as a near-degeneracy, not a size. T172
          (frame_beyond_probe.py, PARTIALLY-PORTABLE, 44/44) then tested
          how far that map carries beyond frame A — and the map is
          portable while the numbers are not: 13 of the 16 links transfer
          UNCHANGED to a gap-blind frame, to ν = 3 and ν = 8, to
          non-prime-power anchors and to both congruence classes mod 4
          (none breaks; the three that shift are exactly the three that
          carry a number), the indefiniteness is localised at the SIEVE
          HORIZON rather than at any frame or zone, R1&apos;s
          near-degeneracy persists on all 54 windows while its rate
          spans h^−1.83 to h^−2.87 across the 9 legs, and the sharpest
          new result is a placement finding: scrambling atom positions
          at a fixed Λ-value multiset raises the angle by a median 6.6e4
          and removes the decay — the collapse belongs to the actual
          prime-power placement, not the value multiset. T173
          (frame_rate_probe.py, DEFICIT-VARIES, 40/40) then showed the
          demanded rate is itself a frame datum: q = 1 − s is an exact
          identity and q &lt; 1 on every frame of an eleven-member
          family, so the h^−3 target was the q = 1 idealisation; with
          the calibrated dimensionless gap functional the deficit
          between demand and delivery becomes THE number of the phase —
          +0.155 ± 0.102, 1.8× flatter than either side, invariant under
          the anchor-to-grid rule and the lever split, not yet constant
          in ν (2.4σ) — and no frame closes it (best upper edge +0.170):
          frame shopping is over, by numbers. R1 now stands frame-free:
          the delivered near-degeneracy closes slower than the relative
          gap it must pay, by a tenth and a half in the exponent. T174
          (cancellation_identity_probe.py, PARTIAL-CANCEL, 37/37) then
          exhausted the gauge route by theorems: everything
          multiplicative cancels exactly — including the entire
          μᴾ₁ ≈ h⁻² channel, which is exactly absent from the deficit —
          the one shared channel (the KMS ladder shape against k²) is
          certified unconditionally at under five percent of the
          deficit, and the additive arch/comb mixture provably admits
          NO factorisation: neither term alone has a positive Schur
          floor on a single one of 84 cells. The bigger result is a
          direct measurement: on a frame-rule-free rectangle the
          deficit stands at +0.1111 ± 0.0222 — five sigma from zero,
          0.4σ from T173&apos;s +0.155 — and its driver is identified as
          comb density per lag cell (the deficit falls monotonically
          +0.281 → +0.062 as the window sees more prime powers per
          cell; on legs log ν and log density are collinear at
          r = −0.921, so T173&apos;s ν-driver was the density channel
          under another name). P6 does not become a theorem — only its
          multiplicative half — but the deficit no longer needs an
          invariance argument: it is measured on a D-rule-free surface.
          The theorem cores of T172–T174 are load-bearing as v560. T175
          (phase_placement_probe.py, PHASES-RESIST, 38/38) then measured
          the placement phases directly, and they are real and — for the
          first time in the series — CAUSAL: displacing the six heavy
          low atoms inside their own lag cells moves the ratio linearly
          over four declared decades with an exact rebuild control
          (dlog R/dδ = 879), but no phase formula is certifiable — the
          response is provably non-smooth, since the Schur floor is a
          near-degeneracy (GAP down to 1.6e-9) and the fitted first
          harmonic falls a factor ~300 short. The per-anchor scatter the
          phases were meant to explain largely dissolved into the error
          bar (jackknife excess scatter 0.1224 → 0.0000, robust pull
          1.41 → 1.22 — an OLS independence assumption, not physics; the
          five-sigma deficit untouched), and the dense limit sharpened
          into the endgame question: the deficit falls monotonically
          with comb density over three decades and the densest reachable
          bin is consistent with ZERO — but true zero, power-law
          approach and a low plateau are indistinguishable under this
          sieve. T176 (dense_limit_probe.py, SITS-AT-ZERO, 24/24) then
          ran that larger sieve — the last decidable measurement of
          phase 2: ATOM_MAX raised twentyfold to 2.5e7 (the sieve turned
          out never to have been the binding constraint), the density
          ceiling from 361 to 6120, two new preregistered ratio-4 bins
          both consistent with zero (pooled +0.0496 ± 0.0372) — no sign
          change, no stabilisation above zero, the plateau window
          narrowed by a factor 3.6 (0.264 → 0.074); the derivative
          identity closed (Feynman–Hellmann, retiring T175&apos;s
          measured slope as a separate claim), the 300× harmonic anomaly
          explained as a POLE (the intervention runs along a
          near-degeneracy that is already there at δ = 0), and one
          anti-promotion: the bin estimator is ladder-dependent, so
          every quoted deficit must cite its rung ladder. True zero,
          power-law approach and low plateau remain indistinguishable —
          the next decade of window would cost a factor one hundred —
          and the phase-2 measurement programme closes here as planned;
          the exact cores are load-bearing as v562, and the work
          continues in the classification papers and the backflow
          lines. Sandbox; not RH evidence.
        </p>
      </DiarySection>

      {/* 25 — Meaning */}
      <DiarySection
        id="meaning"
        plain="Four honestly calibrated levels of what this could mean — from done-and-verified to dream-not-claimed."
        eyebrow="25 · What would it mean"
        title="Four calibrated levels — with caveats"
        badge="sandbox"
      >
        <ol className="space-y-3">
          <MeaningLevel
            level="Done"
            badge="machine-verified"
            title="Mechanism + finite RTF + Weil structure + amplitude route + proof package"
            body="v535–v537: Hecke, Eichler, half-integral bridge. v538: one finite relative-trace identity. v539: Weil structure fully identified up to two explicitly isolated obstructions. v540: amplitude Dirac + geometric polarisation with the Cohen seed Θ(d) = −48·L(−1,χ_d), the universal even-k deletion as square-class double counting, the positive linear carrier with plus-only ζ-balance, and the exact FE — with the open boundary λ* named inside the claim. v541: the T78–T85 proof package — matching lemma proved exact-integer on [4, 10⁶], transport ledger closes exactly (Δ_pole ≡ Δ_conv ≡ 0 proven), character-exact signed envelope, arch internal via Legendre duplication, coherent class closed by the λ-equivariant CM channel — with the two named limits inside the claim."
          />
          <MeaningLevel
            level="Near"
            badge="sandbox"
            title="I5 geography complete; the induction compressed to one sign plus one convention, and assembled end to end"
            body="After 230 probes (5360/5360 sandbox checks) and the promoted modules v535–v564 plus v569–v570, v573 and v576–v596 of this arc: matching lemma closed; I5 geography complete; and the induction that would carry I5 is compressed from one matrix inequality (T104) to a sign the coarse-to-fine recursion already carries plus one declared accounting convention (T124/T125) — with the relay mechanism certified step by step, 400/400 rungs, a single certified step at zone 155,921 (T115), and the finale assembling the whole chain on 52 zones, its load-bearing spine 96.2% identity or Cholesky certificate with the Harnack pair no longer in it (T125). What remains TFPT-specific is exactly ONE object: I5 in one-family form ⟺ Weil positivity ⟺ RH; what is missing for any infinite statement is uniformity in the zone index — now phase 2 of the diary (T126+), where the seam architecture is finished, both remaining inequalities are proof-shaped, three of the four resulting points stand at their preregistered bars (T128), the kappa law falls once while the curvature theorem underneath it stands on all 436 transports (T129), the graded-to-uniform bridge stands as an identity carrying both deep seams while the curvature bound is reduced to one exponent (T130), and the self-supply loop is built one number short of closed — the epsilon-to-floor sandwich and Perron sign constancy are theorems, and M25 reduces to positivity of the pole-free section with nine decades of slack (T131) — two irreducibles remain; the identity block underneath that map is now load-bearing v542, the two reverse-flow parts landed (T132: the Beurling–Deny triad as an operator discriminator, spectrum-only; T133: the certificate audit that hardened v379 to an exact positive-mixture Gram), T134 closed the existence half of the pole-free floor while every cheap route fails by sign (the surviving opening is an M-matrix question), T135 showed the seam DtN admits a bounded faithful state where the Weil window provably has none, T136 closed the a-priori radius item of that M-matrix question by Varga's identity while the exact bookkeeping put the whole degradation in the margin and M17 closed negatively, and T137 made the long-lag support an explicit arithmetic stripe set and certified the whole absolute-value envelope family DEAD (ρ(|E|) ≥ 1.32 from below on 35/35) — the thirteen mature statements of both parts are load-bearing as v543 and v544, T138 found the mechanism of the compensation (the sign law is interval geometry, and the m-paired certificate removes the arithmetic wall on all 77 dead blocks), T139 refuted the classical decay lemma at its hypothesis for an arithmetic reason while DERIVING that sign law from one exact telescoping identity — the residue is one named object: a signed inequality at stripe distance b ≤ 16 — and T140 gave that inequality an exact finite core per zone (ρ(W) = λ_max(K^½HK^½), a closed-geometry coverage kernel times a mass-plus-Dirichlet form) with all the D-dependence in the geometry; the residue is now a zone-uniform discrete Hardy inequality — which T141 then attacked directly and which RESISTS, with its resistance located: four exact identities put it in classical two-weight shape, but the certified constant is not zone-uniform (D^−0.366 ± 0.036) while the exact object it bounds is (D^−0.229 ± 0.007), the additive shape is dead as a shape at its own exact Weyl floor and the joint shape fails at the normalisation alone, so the residue collapses to one closed conductance profile with Y ⪯ K⁺ and Ω ≈ 1 (the identity blocks of T140 and T141 are load-bearing as v545) — and T142 then CONSTRUCTED that profile: the capacity decomposition exhibits the optimal Hardy weight exactly (Ω = 1 by a projection identity, against T141's guessed 20.7–2724), the certified chain misses by a constant factor 2.27–2.45 flat in D, and the rank ladder closes the whole comparison path — no comparison argument can deliver D-uniformity — and T143 then ran the sharp capacity-Rayleigh route itself, which CARRIES: the exact form is an identity on all 26 windows, Maz'ya's capacity criterion applied to the gap form lands inside its window [1/4, 1] with a zone-uniform loss factor (D^−0.048 ± 0.010), the supremum lives on closed families (intervals in node coordinates, certified by full enumeration on the small border blocks), so the residue is now ONE named inequality — a non-Markovian Maz'ya capacity bound cap_E(A) ≥ |A|·λ₀/c₀, whose interval structure points at Muckenhoupt's 1972 two-weight calculus — and T144 then ran that interval route, which CARRIES: the interval class is exhausted exactly (11.4 million intervals via a Cholesky prefix-sum identity), the closed two-weight sup lands inside the Maz'ya window flat in D (B_res·λ̂ = 0.6694–0.7813, D^0.013 ± 0.005), the family restriction falls entirely (a max-density-subgraph bound covers all 2^m sets — the flattest number of the probe), the Markov perturbation route is certified dead, and the certified chain λ ≥ 1/(c₀·κ_up·c_glob·B_res) has exactly ONE unproven input: the absolute Maz'ya constant c₀, whose sharpest shape S1′ is a Muckenhoupt-type hypothesis — and T145 then ran the proof attempt itself: the Maz'ya proof transcribes step by step, the Markov property sits in exactly one line (M4), that line splits with the mass half a theorem that dominates (σ_tot = 0.2145–0.4425 < 1 everywhere), c₀ becomes explicit (best 2.248–4.227, flat at D^0.028 ± 0.017) and S1′ is CERTIFIED per window on 64/64 windows, while an explicit no-go proves that an a-priori bound on the minimiser's level profile — the level lemma L1 — is necessary and cannot be replaced by any weaker hypothesis (the identity spine of T142–T145 is load-bearing as v546) — and T146 then ran the proof attempt for exactly L1, and the level lemma STANDS: closed on the measurement surface as a chain of theorems and certified window inequalities with no step reading the minimiser — the proof lever is the resolvent identity ψ = λRψ itself (the Θ(D³) smallness turns from obstacle into tool), Davis–Kahan is instrumented and discarded (the spectrum bottom is a near-degenerate block), the cake base is free (Maz'ya's classical dyadic 8 falls to 2) and c₀^ap = 3.9042–4.8488 on 64/64 windows lands at the size of Maz'ya's classical Dirichlet value 4, with the chain closing at loss factor 0.0422–0.1586 (the a-priori core is load-bearing as v547); the one genuine remainder is D-uniformity for ALL D — the asymptotic delocalization of the Green columns — which T147 then reduced to an exact identity: Γ = √Q★ · Sw splits all-D uniformity into a purely spectral and a purely geometric factor, both certified on the surface (Sw ≤ 4.6438 by LDLᵀ inertia, Q★ ≤ 2.8634), the classical decay route is computationally dead — delocalization itself IS the bound — the mechanism is named (a Toeplitz-minus-Hankel section whose bottom eigenvectors are Fourier modes, Szegő/Widom) with the sharp prediction Q_B ≤ 2|B| holding at 1.375–1.839 × |B|, the 3 open R4 border blocks close by the same mechanism, and one statement — a Szegő theorem for the diagonally reweighted section — lifts the chain to all D — which T148 then dissected: the second factor Sw closes on the surface via an LDLᵀ layer-cake certificate (Sw ≤ 1.9587, flat in m, certified per D-stratum), at the price of an honest negative — the arithmetic Toeplitz symbol has f(0) < 0 on all 48 windows, so positive definiteness is a finite-section effect of the minus-Hankel part and the lumping, and the KMS order-2 hypothesis survives only as a measurement (α = 1.64–1.99) — while the lifting statement RESISTS with its hypothesis isolated to a single named scalar: the total variation of the log-whitening weight, proved to be the roughness and not the conditioning by a controlled weight-class experiment (two BV classes flat, a TV ~ m class diverging as x^1.994 at identical κ_Λ = 4); the first end-to-end number arrives — the certified chain delivers 8.36–15.86% of the true gap on all 48 windows (median 11.74%), and a-priori-shaped factors still give a valid lower bound at a factor 10.5 — verdict ONE-INPUT-MISSING: exactly one scalar side (Q★, the lifting statement) lacks an m-free certified statement, target ν_L ≲ 34 against 282 measured (the identity/certificate core of T147+T148 is load-bearing as v548) — which T149 then attacked through the gauge freedom (PARTIAL-SMOOTHING, 30/30, 44 windows m = 277–1393, 9 preregistered gauges): the whitening diagonal is a free gauge (an identity plus one Rayleigh step), and the constant gauge — the geometric mean of the Jacobi diagonal — ELIMINATES the blocking scalar exactly (TV(log Λ̃) = 0.0000 against T148's 11.93 at x^0.444) at a certified sandwich price (σ ≤ 5.5789, κ̃_up ≤ 2.3146 against 1.2647); yet the two-regime decomposition refutes the hypothesis itself — gauges that kill the flutter and keep the macro profile move ν_L̃ by at most 0.9%, only macro-removers improve it at all (factor 1.27, target met on 21/44 windows) — so the roughness ν_L̃ responds to was never in the multiplier: the missing input is RELOCATED, not closed, and with the constant gauge it is now exactly the smoothness of the deep modes of the pure Toeplitz-minus-Hankel section, in the ladder form ν_k ≤ C·k² with m-free C that the Dirichlet control singles out (classically ν_k = πk² exactly; measured C = 18.66–44.61, x^0.272), with the flutter amplitude measured flat (0.064–0.182) as the second lever — the family-maximum chain delivers 9.52–18.45% of the true gap certified (median 12.68%) and improves the a-priori side by up to 1.28× — which T150 (MODE.LADDER, ONE-TERM-MISSING, 36/36, 72 prime-power windows m = 50–1393) then ran, and the mechanism acquired a name — parity: the form is exactly the compression of the full symmetric Toeplitz section onto its antisymmetric parity sector (U₋ᵀT_M U₋ = T − H, certified, cross block ≤ 4.0e-16), and the sign-changing symbol's entire negative inertia sits in the EVEN sector (72/72, LDLᵀ) — T148's honest negative is explained rather than worked around; the flutter amplitude became a certified form functional (0.0606–0.1993, falling, atom budget 4B√N closed), the whitening diagonal an explicit zone functional, and the arithmetic atoms turn out to be co-responsible for the positivity — the archimedean section alone has 2–7 negative eigenvalues, so the additive perturbation route is structurally dead while the multiplicative gauge step closes; one gate remains: the ladder constant still grows, C ≤ 43.391 = 13.81π certified per stratum at x^0.258 against the flatness bar 0.25, the gap from π exactly the arithmetic leakage of the bottom mode; the identity/certificate core of T149+T150 is promoted as v549 — which T151 (ODD.LADDER, ODD-CARRIES, 27/27, 72 prime-power windows m = 96–1491) then closed by REROUTING OFF THE LADDER: the odd grid steps over the symbol's negative window (θ_c/θ_1 = 0.328–0.407 on 72/72 — positivity is a grid fact), the bottom spectrum is certified against the parity Laplacian (λ_k ≤ S·μᴾ_k, S = 1.1019–2.3870, LDLᵀ, with f(0) < 0 absorbed by the Rayleigh floor), and a discrete Sobolev step at the odd sector's virtual node yields a per-mode bound LINEAR in k with a NON-GROWING constant (C_S = 11.5137–19.5731, trend x^0.020±0.007) where the quadratic ladder grew; the symbol route is a computed dead end (the section form averages the symbol against a Fejér kernel instead of sampling it, margin vacuous on 72/72 — the local model is a matrix statement, not a symbol one; T150's rest 3 answered negatively), the archimedean minimum is closed AND attained exactly (min Λ^arch = c^arch₀ − c^arch_{M−1}, deviation 0.0; rest 2 done), end to end moves to 2.01e-2–3.52e-2 of the true gap with the bottleneck relocated from Q★ to Ψ, and ONE scalar remains a fit: the bottom pencil ratio R = K_bot/κ = 3.3634–9.7108 (flat, x^0.037±0.015) — the T145 no-go breaks exactly there (x^1.986); the identity/certificate core of T151 is promoted as v550 — and T152 (PENCIL.RATIO, ONE-TERM-MISSING, 37/37, 60 prime-power windows m = 149–1445) then attacked exactly that scalar and refuted both hoped-for gifts from the smooth kernel: the archimedean part is itself NEGATIVE in the odd sector (λ_min = −2.81…−1.84, O(−m²) in pencil normalisation) and the atom part too (−1408…−25.6), so positivity is a CANCELLATION between geometry and arithmetic, not a property of either half — yet the floor closes structurally: a Schur two-block criterion (Schur 1917, both Cholesky floors subtracted) with a FIXED 16-mode low block certifies κ ≥ 0.225–0.250 on every window, flat (m^0.006, quartile medians identical), consuming one unproven block inequality (B_HH ⪰ t·I; entrywise/Gershgorin provably insufficient), so R ≤ K_bot/t = 4.41–8.47 (m^0.099) is certified on both ends with ONE term missing — an m-free ceiling on K_bot (every structural route loses 4–6 orders of magnitude; the no-go breaks exactly there at m^1.994) — and the Ψ map finds the real prize: Ψ is 1/λ_min(E) to within 1.19, and 89–92% of it sits on the eight modes the certified ladder already controls, a 3.46–5.29× end-to-end reserve — and T153 (PSI.LADDER, REBUILD-RESISTS, 33/33, 18 prime-power windows m = 96–1430) then carried out exactly that rebuild and refuted it by its own factor: Ψ is pinned between a₁/λ₁ and 1/λ₁ with a₁ = 0.7695–0.8306 = exactly 8/π² (the first parity sine), so it is determined to within 1.204–1.300 and the reserve never existed, the head/tail replacement falling 2.31–6.12× BELOW the hard lower bound (the maximising subset is aligned with v₁ — the head IS the value of Ψ); the theorem-grade collapse Ψ ≤ const/(t·μᴾ₁) nevertheless CLOSES Ψ as a term (retiring Charikar's licence and every per-window diagonalisation, needing nothing beyond the certified floor t = 0.2400–0.2500, flat), the end-to-end fraction moves to 1.01e-2–3.92e-2 (net gain 0.50–1.11, median 0.62 — the retired level constant and the collapse cost cancel), a T152 sign REVERSES (on the bulk parity block the archimedean part is POSITIVE, 1.0362–1.4143 > t = 0.25, after an m-free peeling of at most 8 modes; the atoms are the negative part there), four block candidates die (entrywise and block Gershgorin, mode-index Toeplitz+Hankel, recursive Schur — scale-invariant) while the live Kato-type route certifies positivity per window at a distance shrinking x^-1.778, and two inverse-iteration steps (A⁻¹L_P)²t_k give a FLAT ceiling K^F = 1.432–2.369 on a fixed-size-8 certificate against the true K_bot = 1.102–1.896: two terms remain, and both lose to the SAME missing object — a Green/alignment estimate of where the bottom eigenvector of the section sits — and T154 (GREEN.ALIGN, ALIGN-RESISTS, 29/29, 12 prime-power windows h = 50–1077) then attacked exactly that estimate, and the ceiling CLOSES, exactly and at fixed size: the sixteen-column certificate span{t₁..t₈} + A⁻¹L_P·span{t₁..t₈} needs NO residual argument at all (Ritz values are upper bounds for the eigenvalues of the same index, Courant–Fischer 1920 / Cauchy 1829 — the direction correction of the part; Temple/Kato is a floor device) and carries K^F = 1.1019–1.9964, flat (x^0.094±0.058), agreeing with the inertia-certified K_bot to 5.17e-7 on every window — the eight size-m LDLᵀ counts are retired from the ceiling step, and the T145 no-go stress confirms the instrument explodes (x^2.29) exactly where flatness is false; the floor half is refuted at fixed size (the Temple/Kato correction is 6.3e3–5.6e9 too large, x^4.50 — residual O(1) against a target O(m⁻²)) and the obstruction is NAMED: seven of the eight bottom directions of A agree with the bottom of L_P to 0.15–1.35° (median), ONE sits at 82.93–89.79°, and that single misalignment IS the collapse price — one Cholesky of A − γI per window recovers it in full (4.408–7.985), moving per-window end-to-end to 4.45e-2–3.13e-1, inside the target band 3e-2–3e-1, while the m-free-in-shape number stays 1.01e-2–3.92e-2 (the closure buys uniformity and cost, not size — two numbers, never conflated); the missing fixed-size ingredient reduces to ONE arithmetic-free number (the L_P floor on the complement of the eight bottom Ritz directions — only the tridiagonal parity Laplacian and one 8-dim subspace — measured flat at 5.93–8.25 μᴾ₁, worth 91–100% of the price), the Kato route to R1 fails by the SAME geometry (loss 1.97–121 attained 24.9–89.8° apart, the minimiser sitting 96.5–100% on the eight modes above whatever cut is chosen; the Hankel-reflection culprit refuted twice over), and both remaining terms are now UNIFORMITY terms with certified per-window numbers (R1′ the m-free block floor, R2′ the m-free bottom-mode floor); the fixed-size ceiling certificate is promoted as v551 — and T155 (BOTTOM.FLOOR, bottom_floor_probe.py, FLOORS-RESIST, 27/27, 16 prime-power windows h = 142–1293) then attacked exactly those two floors: the complement floor — the one m-sized object left in the bottom-mode chain — becomes an EXACT fixed-size certificate, min over v ⊥ W of vᵀL_Pv ≥ μᴾ_{K+1} − λ_max(M^{1/2}(I − GGᵀ)M^{1/2}) with G = T_K Q_W, an identity plus one Rayleigh bound valid for every m and every subspace, reproducing the size-m eigenproblem to 0.999999783–0.999999958 on 16/16 windows with K = 12 sufficient everywhere (both closed-form controls hit to 1e-9, including the configuration in which the certificate is worthless); the defect is localised at mode 1 itself (uncovered fraction 0.539–0.611 — the 83–90° direction of T154 is demystified: it lives on modes 9–12, which is exactly why K = 12), the collapse price 3.250–8.471 is recovered 78.8–100.0% at fixed size (in full on 12/16 windows; end to end 3.28e-2–2.83e-1 CERTIFIED AT FIXED SIZE, the declared 4e-2 bar missed by 1.22 at the bottom and reported), and two repairs are refuted with their mechanism (the pencil-ceiling chain dies on κ = λ_max(B) growing x^2.82; wider subspaces raise the floor but explode the residual, Temple 0/32); on the block side λ_min(B_HH) = 0.2430–0.4249 is certified positive while the direction-aware 2×2 split dies on the coupling (‖B_wd‖ too large by 16.8–908), the local atom norm falls short by 1.50–5.12, the deeper cut does not terminate — and the strongest negative statement of the series arrives: the arch reserve IS the symbol infimum (to 5e-4, a theorem candidate, Szegő 1915/Widom 1958) but the FULL symbol infimum is −714.2…−7.6, negative on every window against a positive section floor, so NO symbol argument can ever produce the block floor: the mechanism must be Fejér cancellation in the finite section. Two open terms remain, deliberately not merged: R2″, an m-free upper bound on ONE 12×8 object; R1″, the m-free atom part on the bulk with every symbol argument excluded — and T156 (TWELVE.EIGHT, twelve_by_eight_probe.py, TERMS-RESIST, 37/37, 16 prime-power windows h = 142–1293) then attacked exactly those two objects, and both are now single scalars: on span{t₁, A⁻¹L_P t₁} the coupling t₁ᵀAy₁ = μᴾ₁ is an identity with no A in it, so the t₁-loss is an EXACT closed function F(P, r) of the Kantorovich product P = 5.6e2–1.0e6 and the inverse moment ratio r = 2.7158–3.4089 (flat, verified to 2e-16), the two-line theorem r ≤ 1/(Ls) ≤ 1/p₁ is tight to a factor 1.03–1.16, and R2″ is therefore ONE m-free lower bound on the angle p₁ = cos²∠(t₁, e₁(A)) = 0.2010–0.3282 (55.1–63.4°, flat) — everything above it is a theorem, with one separate measured debt (the 2×2 model dominates the 8-dimensional defect on 16/16 real windows but FAILS on 8/8 no-go sizes, so it stays MEASURED); the interlacing route is refuted in both forms (a rank hole at K = 12; empty at K = 8); and on the block side the mechanism is identified against expectation — the expected arch inequality ≥ 1 is FALSE (inf = 0.8226–1.3973, below 1 on 3/12; the weaker inf ≥ t survives at factor 3.29–5.59), the Fejér damping is worth only a factor 5–31 where the split needs 3–1981 growing, the additive split is DIVERGENT (atom norm h^2.31), and the positivity is an ALIGNMENT fact: arch 1.47–91.71 and atom −91.27…−1.12 cancel to 0.2661–0.4436 on the minimiser, which sees only 6.2e-3–0.50 of the atom operator at 52.7–90.0° from the atom-extremal vector; the balance is 9 THEOREM rungs (3 new), 6 CERTIFIED, 3 MEASURED (the zero-MEASURED target missed), the end-to-end number is unchanged at 3.28e-2–2.83e-1, and the no-go breaks on THREE axes including a collapse of p₁ itself to 7e-15–7e-10. and T157 (ANGLE.FLOOR, angle_floor_probe.py, ANGLES-RESIST, 32/32, 16 prime-power windows h = 142–1293) then attacked exactly those two angles, and neither falls, but both change shape: the tail of the resolvent route becomes a THEOREM — the sine-block confinement ‖γ_H‖² ≤ λ₁/(t·μᴾ₁₇) ≤ (S/t)/ρ₁₇ = 0.0165–0.0293, from the certified pencil floor and ladder ceiling alone, so the bottom eigenvector lives 97.1–98.4% inside the first sixteen parity sines (T146's measured '98 percent' replaced by one line; with the pencil ceiling κ in place of the flat ladder S the bound would be vacuous, 43.5–38632) — and the floor p₁ ≥ ĝ₁²(1 − tail) = 0.1968–0.3228 is 97.9% of the measured angle, with ĝ₁² = 0.2010–0.3282 the ONE measured fixed-size scalar left (the classical Rayleigh angle bound is empty on 16/16, its block version is empty at every J although the block mass p_blk(2) = 0.9369–0.9993 shows t₁ lives in the bottom pair, and the angle-free Cauchy–Schwarz route loses exactly the Kantorovich product P = 5.6e2–1.0e6); the structural gem: 1/s = (S_L)₁₁ = 2.3359–6.2049 flat — the whole first term is a bound on ONE diagonal entry of the 16×16 Schur complement the chain already forms, and the inversion-free Cauchy–Schwarz ceiling misses it by 5.4e2–5.2e5 because the cancellation is nearly complete; the arch half of the second term is uniformly certified for the first time (an executed adaptive Lipschitz ceiling, 12/12 windows, 255–3139 evaluations, cost h^0.85), the two extremal vectors sit at OPPOSITE ends of the band (atom θ/π = 0.0158–0.3762 against arch 0.9901–0.9995 — a proof must use the θ-growth of the arch ratio, not its infimum at π), and the alignment term stays a per-window domination with quotient 1.0003–1.0907, a 7.3e-4 margin and a shrinking trend; the balance is 9 THEOREM, 2 CERT-UNIF, 2 CERT-WINDOW, 3 MEASURED — all three measured steps now fixed-size in their statement — and the no-go breaks on FIVE axes (p₁ x^−4.818; S x^2.289; the confinement to vacuum; the resolvent identically zero; (S_L)₁₁ x^2.268); the four instrument candidates of T155/T157 are promoted as v552 — and T158 (SCHUR.ENTRY, schur_entry_probe.py, ENTRY-RESISTS, 36/36, 28 prime-power windows h = 254–1393) then found exactly that cancellation-seeing bound, and it is a THEOREM: the Thomson dual form s = max_x(2x₁ − xᵀBx) turns the entry into a Dirichlet maximum, so every trial vector bounds it from the right side — which is exactly why Cauchy–Schwarz missed by 3.13e3–5.18e5: it evaluates a maximum at a single direction, the wrong variational structure — and the Cholesky ladder g_K = Σ y_j² of strictly positive terms (monotone partial sums on 28/28, starting at T157's route 1/g₁ = â) pins the entry at K = 16 to 1/g₁₆ = 2.9670–7.9664 against the true 2.3359–6.3868, tight to 1.1323–1.2738 and flat, while the Green route is the identity itself — span{t₁, A⁻¹L_P t₁} attains s exactly (L_P t₁ = μᴾ₁ t₁), so the entry is two-dimensional the moment one Green column is granted, and the fixed sine truncation pays exactly the factor 1.13–1.27 for not needing it; the honest negative: T157's growth pointer is REFUTED — the growth is real (θ^−1.259…−1.438) and the binding vector does sit in the lowest dyadic band (0.68–0.999 of its mass), but the atom negative mass grows FASTER (θ^−1.546…−1.744), band-local domination fails in every band on 21/21, and the off-band coupling that actually carries the inequality exceeds the margin by 660–7.7e5 — R1″ is now a question about the sign structure of the off-band arch entries; the T156 debt is relabelled MEASURED → CERT-UNIF (both sides fixed-size certified, margin 1.1969–1.3568 on 28/28), the end to end survives the substitution at a cost of 1.000–1.724, the no-go breaks five-fold, and the measured-step count falls from three to one (the balance: 6 THEOREM / 3 CERT-UNIF / 4 CERT-WINDOW / 1 MEASURED); the five T158 candidates (P1 the dual form plus the positive ladder, P2 the two-dimensionality, P3 the 1/g₁₆ sharpness, P4 the negative result on the growth pointer, P5 the relabelling) stay PENDING, to be bundled with T159 — and T159 (EXACT.FORM, exact_form_probe.py, FORM-RESISTS, 41/41, 24 prime-power zones, 23 M1 windows h = 142–1293) then executed the cancellation algebraically, and the algebra was delivered while the bound was not: SEVEN machine-checked identities, all to 1e-12 of the absolute scale (the honest bar: relative to the cancelled total they hold only to 1.0e-12–1.2e-8 — the cancellation eats 4.0–8.1 digits in double precision, which is the quantitative reason an m-free bound cannot be read off numerically) — the y-reduction, the exact lag sum xᵀB_LL·x = Σ c_d·w_d, the CLOSED 256-term Dirichlet-kernel weights (1.42e-15–2.85e-15), the GAUGE IDENTITY Σ w_d = 0 exactly (Toeplitz-minus-Hankel annihilates constant lag vectors, so the form is BLIND to the lag mass — exactly where the h²-sized halves live), the two closed scalars w₀ = Σx_k²/μᴾ_k and 2w₀ − w₁ = ‖x‖², the p-fold Abel identity exact to p = 5, and the TmH signature; the honest answer to the kernel question is NEGATIVE with an exponent: the h² does NOT telescope away — both halves grow h^+2.104, i.e. like the weight w₀ itself, and it cancels only in the sum — and four routes are closed with exponents (T158's own fixed sixteen-vector h^+3.510 because cond(B_LL) ~ h^+2.887, the one-sine rung h^+3.066, the preregistered ansatz family h^+3.001, ℓ¹×sup Abel pricing h^+3.604); the genuinely new structural win is a SIGN LAW: the archimedean bulk block is exactly a symmetric Z-matrix, raw, on every window, with a closed sign-based Collatz–Wielandt floor 0.82–1.16 comfortably above the 0.25 target — but the atom block obeys no sign law (0.31–0.49, noise), so the criterion is vacuous on the full block; the balance jumps to 13 THEOREM / 6 CERT-UNIF / 4 CERT-WINDOW / 1 MEASURED, a numerical horizon is declared (cond(B_LL) > 1e12 past h = 1292), and the theorem cores of T158+T159 are load-bearing as v553 — and T160 (PAIRING, pairing_probe.py, PAIRING-RESISTS, 46/46, 19 of 20 prime-power windows h = 199–1256, one dropped at the declared horizon) then attacked exactly that pairing through its correlation structure, and the watershed of phase 2 arrived: W¹ does NOT oscillate (J = 3 sign blocks on every window, flat, with a closed bound J ≤ 130 from the trigonometric zero count — so the Leibniz/block device is refuted by being 5.2–16.5× WORSE than the ℓ¹×sup pricing it was meant to beat, and the head peel is refuted at K* = M on 19/19); what carries the smooth half is three closed sign-definite MOMENT LAWS (m₀ = 0; m₁ = −[S0² + 2ΣP_j²] ≤ 0; m₂ = −[2S1 − (M−1)S0]² ≤ 0, a perfect square; closed for every even p), against which a fixed polynomial witness reproduces the arch half to slack 1.95e-8 of the O(1) total — FLOOR-LIMITED at 4.12–428× the double-precision floor, a second declared numerical horizon: the arch half is neither established nor refuted on any double-precision surface — while the machine-checked SAMPLING IDENTITY shows the atom half IS, identically, a finite combination of Λ-weighted prime sums Σ Λ(n)·n^(−1/2)·cos(t·log n) at 32 explicit frequencies t = π(k±l)/α, needed to relative depth 2.2e-6–1.1e-4 = h^(−2) and measured to cancel only to 0.00–0.37 of the trivial bound: the h² cancellation is the INTRINSIC ARITHMETIC HARDNESS of the problem, not an assembly artefact — the geometric half is evaluated to the arithmetic floor and nothing geometric is left to trade; the total-variation bound U3 became a genuine THEOREM (measured 3.9101–5.4153 against the closed constant 3.9105–5.4153 — four digits), the Z-law is reduced to ONE prime-free trigonometric inequality (the pointwise strengthening refuted — the monotone weight is essential), ρ = 1.0036–1.0140 > 1 flat on every destructive direction is the sole surviving R1″ fact (its inequality form refuted: the margin loses to the arch cross-coupling by 1e4), the composite chain recovers 3.79e-10–2.10e-7 of the exact value and would recover 1.0000 if the atom half were evaluated as exactly as the arch half now is, the balance moves to 14 THEOREM / 4 CERT-UNIF / 3 CERT-WINDOW / 3 MEASURED with 8 refuted families, and the T145 no-go breaks on 6/6 axes including the two new ones — and T161 (CLASSICAL.CLOSURE, classical_closure_probe.py, CLOSURE-RESISTS, 35/35, 19 log-spaced zones h = 50–1445) then ran the two classical rests plus the circularity triage, and the triage returned its most consequential answer: THE CHAIN IS NOT CIRCULAR — the 32 frequencies satisfy t·(2α) = 2πj exactly (a THEOREM: they are the Fourier harmonics of the log-window), the measured cancellation is FULLY the PNT main term ((√X−1)/(¼+t²) matches S(t) to 0.014 on the largest window, a Mellin factor and not an arithmetic saving), and the required depth is δ = 1.1482–1.8809 against RH strength 1/2 — in absolute terms 0.012–0.31 of the LAST TERM of the sum, below the boundary term every partial-summation bound carries, so no strengthening of the ψ(x)−x input (zero-free region, RH, or beyond) reaches it: the chain does not secretly need RH, it needs MORE than RH-strength input would supply, which localises the h² in the SPLIT rather than in the primes (a re-split against the smooth prime term moves the demand to δ_eff = 0.9839–1.3767 — it MOVES, but stays above 1/2 on 18/18); meanwhile R-A's analytic half CLOSES m-free — A = D·Âhat exactly, only the 1/s head binds, the scale-free Bernstein rate ρ* = (3+√5)/2 = 2.618034, the closed head split c^arch = Ψ + D·Ĝ with Ψ D-free and m-free (no peeling), and an explicit degree schedule K(h) = O(log h) (the 'fixed degree' hope refuted with a number; the one residual is the prime-free log-moment ΣΨ_d·w_d, outside the polynomial ladder) — R-B is REFUTED in all four readings (3168/4320 pairs, the aggregate by SIGN) and replaced by a certified off-diagonal FRACTION bound 0.1035–0.1713 ≤ 1/4, flat; the balance moves to 20 THEOREM / 10 CERT-UNIF / 4 CERT-WINDOW / 9 MEASURED / 3 REFUTED, and the closed cores of T160+T161 are load-bearing as v554 (PRIME.SAMPLING.HARM.01). — and T162 (THIRD.SPLIT, third_split_probe.py, DELTA-REDUCED, 30/30) then ran exactly that search, and the third split EXISTS while the exhaustion SATURATES: the archimedean Mellin ladder of the explicit formula lowers the demand 1.88 → 1.38 → 0.93 in closed cell moments, forced and not fitted, but it is an asymptotic series turning around at K* = 2 (past it the residual rises ×15–25); one Abel step makes the demand PRIME-FREE — δ_bnd = 1/2 + log(2κ·‖Δw‖₁/|Q|)/log X exactly, every arithmetic input in the single Chebyshev constant κ = 0.038821 — and its optimal level is exactly one, for the closed reason 32π/α > 1; the FEJÉR split (tapering the trial vector, which by self-adjointness IS pairing against Fejér-averaged Λ-mass) pushes the proof demand BELOW the RH threshold 1/2 on all 18 windows (δ_bnd = 0.133–0.417) at a price in the 1/s ceiling growing h^2.86 — the hardness is RELOCATED into the price, not closed; alongside, R-A′ CLOSES (the log-moment agrees on three independent routes to machine precision, via a Lerch/Frullani integral whose d = 1 term peels as 2·log 2 and reproduces Ψ₁ = −log 2) and R-B′ is REFUTED (the 16×16 Gram form is indefinite on every window — the a-weighted quarter bar survives as a contribution bound); the remaining R2″ question is a PARETO FRONT — does an operating point exist where demand and price are simultaneously affordable? — and T163 (PARETO.FRONT, pareto_front_probe.py, FRONT-RESISTS, 32/32) then answered exactly that question with a THEOREM rather than a failed search: the exchange law δ_bnd(x) = 1/2 + log(2κ·g₁₆·TV(x)/P(x))/log X is an IDENTITY at all 27×50 grid points of the Fejér knob (price and demand are two coordinates of one object; the knob is measured monotone in both, so its curve IS the front, and its endpoints are the chain's own ladder rungs σ = 1 ↔ K = 1, σ = ∞ ↔ K = 16), the front does NOT cross at any flat price (0/27 at caps 1.25/2/10 and at every flat rung tier; only the full chain-derived P_aff = g₁₆B₁₁ ~ h^3.05 reaches 27/27, with the crossing price P_cross ~ h^1.91 sitting 5–152× strictly INSIDE that already-accepted certificate, margin widening), and a new four-line theorem makes the resistance structural: w₀ = ‖a‖², TV ≥ |w₀| by telescoping, and the entry normalisation x₁ = 1 forces TV(x) ≥ 1/μᴾ₁ = 1/(4sin²(π/N)) ~ h² for EVERY admissible trial vector (verified on all 1674 built, slack 5.0–23.2) — so every sub-1/2 demand pays P > 2κg₁₆/μᴾ₁ ~ h², already 5.5–590× above the largest flat cap: R-C‴ ('bounded total variation at bounded price') is CLOSED NEGATIVELY by an inequality, the h² of T162, the crossing and T161's granularity are ALL the reciprocal smallest parity eigenvalue — the spectral gap of the parity Laplacian (KMS 1953) meeting the entry normalisation, never the primes — and the mode sweep K = 16 → 64 confirms it from the other side (a price dividend, a WORSE demand, the TV exponent unmoved at h^1.81); the closed cores of T162+T163 are load-bearing as v555 (PRIME.PARETO.TV.01), and the successor R-E is named with two prime-free arms (arm A: the downstream chain tolerates a growing 1/s ceiling; arm B: the entry functional in a sector whose gap does not vanish like h⁻²), with R-B‴ (h-uniform positivity margin) and R-D (fifth device) open beside it — and T164 (SECTOR.CHANGE, sector_change_probe.py, TOLERANCE-CARRIES (arm A) / SECTOR-RESISTS (arm B), 28/28) then decided BOTH arms of exactly that successor: arm A — the T156 kernel spends the entry ceiling at power exactly one (d log(1−F)/d log U = −1.0005…−1.0000, d log r/d log U = +1.0000 exactly; the declared rule's ε* = 0.50 is disowned as the surface's own tailwind h^+0.391 divided by power one), yet the O(1) gate is DISCHARGED window by window by a Cholesky identity — the single constant U_ref = 4.9008 = max 1/g₁₆ carries all 9 sub-surface windows out of sample, 1/g₁₆ = 1.7527–5.3286 is flat (h^+0.061 over a 9× lever arm, split halves −0.010/+0.103) and every g_K is strictly increasing (Schur 1917) — so R2″ collapses onto ONE quantifier, sup_m 1/g₁₆(m) < ∞, the m-freedom of a certified flat list; arm B — NEGATIVE BY A THEOREM: the entry normalisation is a gauge (Q and TV homogeneous of degree two, x₁ = 1 fixes only the scale, so Q/TV and δ_bnd are the same numbers in every sector, to 8.8e-10 over 5670 sector-by-vector combinations), the full space is strictly worse (μ₀ = 0 exactly), every floor-flattening shift pays the identical exponent back as a transfer factor (floor × transfer = 1/μᴾ₁ identically, exponent sum +1.997 at every shift to 1e-9), and the whole T156 kernel is gauge-invariant (3.5e-16 on 54 combinations — no fifth device from a sector change, R-D settled in type); the surprise, immediately priced: an unconstrained ascent on Q/TV overshoots the crossing bar 2κ = 0.0776 by 31–960× on 9/9 windows (η = 7.9–11.1% of the Abel ceiling ‖C‖_∞ = 30.3–930.6 ~ h^+1.185) while T163's TV floor HOLDS on the unconstrained optimiser (TV·μᴾ₁ = 8.30–11.72 ≥ 1) at a price h^+3.299 — worse than the crossing price h^+1.91, so the binding axis is the gauge-invariant, prime-free ALIGNMENT between the weight increments and their partial sums, named R-F; R-B‴ is narrowed (the −0.172 belongs to the T162 quarter-bar object, 0.1001–0.1623 < 0.25 on 27/27, regressing −0.046 on log h against −0.477 on log α — an independent α/h surface is needed) — and T165 (ALIGNMENT.ETA, alignment_eta_probe.py, ETA-RESISTS in the strong form, 30/30) then closed exactly that alignment successor BY CERTIFICATE: the gauge-invariant exchange identity P_pr = g₁₆·R·(TV/(t₁v)²)/μᴾ₁ (machine-checked to 3.3e-16) makes demand and price ONE equation — four named factors: the quantifier, the demand, T163's floor ≥ 1, the KMS h² scale — so with the floor every crossing vector pays ≥ 2κg₁₆/μᴾ₁ and the two R-F clauses cannot be chosen independently; the predicates 'R > 2κ' and T163's crossing criterion agree on 105/105 vectors (R-F is strictly stronger than the R2″ demand and its m-uniform form IS R-E-A — the alignment question was never independent); the price exponent decomposes exactly (h^+3.261 = +1.997 KMS + 1.279 overshoot + 0.093 floor − 0.108 g₁₆; the bar-tight counterfactual stays h^+1.982), the η(Cap) ceiling is a theorem falling short of the requirement by 5.5–392× at Cap = 10 (0/27 ladder points reach the bar), the free optimum's anatomy is measured (low-harmonic: 93.7–96.8% of its energy in the first 32 parity modes vs a 2.4–22.5% baseline, essentially orthogonal to the PNT main term, heavy edges on prime-power cells, a cancelling sum with participation 8–18%), and a genuinely decoupled ν-surface (25 windows over 6 zones, ν ∈ {4,5,6,8,11,16}) settles R-B‴'s confound: the quarter-bar drift is ZONE DEPTH, not window length (factor 19) — at the honest cost of retiring U_ref = 4.90 (sup 1/g₁₆ = 5.7327 off-recipe, +17%); end-to-end on the union (52 windows) every crossing vector has P_pr ≥ 54.6–5865 > Θ_TOL, closed for every h > 84.9; balance 11 THEOREM / 1 CERT-UNIF / 4 CERT-WINDOW / 6 MEASURED, and the honest bottom line is that exactly ONE genuine open object remains in this line: inf_m g₁₆(m) > 0, equivalently a lower bound on the sixteen-step Schur-cascade gain g₁₆/g₁ ≥ c·h^(3−ε) uniform in m — a quantifier over a certified flat list, a cancellation and hence provably beyond absolute-value budgets; the theorem cores of T164+T165 are promoted as v556 (PRIME.GAUGE.PPR.01) — and T166 (SCHUR.CASCADE, schur_cascade_probe.py, CASCADE-RESISTS, 30/30) then dissected exactly that cascade lower bound on a 63-window union of the frame-A and decoupled-ν surfaces: four identities (all theorems) turn the fifteen ladder increments into one readable object — the gain is the near-collinearity 1/(1 − R_K²), rung 2 alone carries a median 59% (K_half ∈ {2..5} on 63/63), and the gain is invariant under B → DBD, so the h³ is a property of the arithmetic Gram block alone; the cancellation matches none of the prepared stories (no entry of the 2×2 block cancels, neither half is collinear alone) and lives in the 2×2 GRAM DETERMINANT — pieces 2.6–630 against a full determinant 1.0e-3–1.7e-2, the same arch-against-atom mechanism as T159/T160 one level up — with an anti-fitting scramble destroying the effect by a factor 4569; the best closed route reaches gain h^+1.319 against the target h^+3.110, an exponent gap of h^+1.791, so the one missing inequality IS the cancellation, restated in its sharpest form as one Gram-minor ratio (equivalently a closed β at rung 2 to 3.6e-4, or a closed near-null vector of A in the low modes); two honest contract corrections are recorded (F2 is strictly weaker than inf g₁₆ > 0 since B₁₁ ~ h^+3.110, and the confinement route was a sub-surface artefact exposed by raising the eigen-horizon), U_ref moves to 7.45 on the union (declared off-recipe), and the nine candidates P166.1–P166.9 stay PENDING — and T167 (NULL.VECTOR, null_vector_probe.py, VECTOR-RESISTS, 39/39) then attacked the most constructive dress and closed it as a construction: the pivot identity uᵀQu = 1/g_K + δᵀQδ (δ₁ = 0) makes every candidate's excess an exact number, at K = 2 the closed vector (1, −Q₂₁/Q₂₂) is EXACT (a theorem, ρ = 0 on 63/63, B₁₁g₂ = 1/(1 − r₁₂²) an identity, certified gain exponent h^+2.921 = h^(3−0.079) on frame A) so the third dress collapses onto the second; perturbation theory is closed off structurally (the Kato series converges fast, radius 0.067 — to the WRONG object: the exact bottom eigenvector is itself a useless trial with overlap 0.083 and relative spectral gap 5.4e-07); the threshold is mildest at K = 2 (1.90×, h^+0.353), against the hypothesis and exactly where the vector is free, while at K = 6 the accuracy curves diverge (separation h^+1.138, ratio already 264×); the scramble sharpens to a change of TYPE (the 2×2 diagonal itself loses positivity on 8/8 windows — g₂ does not exist there); and the unification is an IDENTITY: eps_ent(K) = ρ·(1/g_K)/S_K to machine precision — the determinant ratio, the scalar 1 − r₁₂² and the entry threshold are ONE inequality, not three. The rest is R1, a single scalar (an m-free upper bound on 1 − Q₁₂²/(Q₁₁Q₂₂) ≤ C·h^(−3+ε), three closed lag sums; required relative accuracy 1.25e-05 median / 3.19e-08 worst); the theorem cores of T166+T167 are promoted as v557 (PRIME.CASCADE.VECT.01) — and T168 (LAGRANGE.MINORS, lagrange_minors_probe.py, MINORS-RESIST, 39/39) then ran the Lagrange identity on that scalar, and the sum of squares is REAL: the h×h arithmetic kernel is positive definite on all 63 windows (hard-fenced as a per-window statement, the Weil criterion never tested, assumed or reverse-inferred), the mode vectors are Euclidean-orthogonal (the near-parallelism is created entirely by the arithmetic metric), the Wronskian minors telescope in closed form and are MAXIMAL in norm (all smallness sits in the PSD kernel at one closed vector; in the eigenbasis the sum is thin and raw), the exponent ledger closes with every factor but ONE at the right power, and that factor is the single ratio t* = Q₁₂/Q₁₁ of two closed lag sums, needed to the same accuracy T167 measured from the other end — by T168-TH7 the target itself rewritten: the hardness is SELF-SIMILAR under exact reformulation — and T169 (TSTAR.RATIO, tstar_ratio_probe.py, RATIO-RESISTS, 41/41) then PROVED that self-similarity: every genuinely closed candidate for t* misses by O(1) and diverges from the threshold at h^+1.99 (the block is atom-dominated, its archimedean diagonal entry negative on 63/63 — the closed family was structurally hopeless), while the only candidate meeting the threshold, √(â₂₂/â₁₁), does so by a new identity (T169-TH7) that reintroduces det Â — the loop closes as an identity, T167's scalar, T168's factor and T169's candidate are the same object; the real gains: the first CERT-UNIF in weeks (ν₁ ≤ max(â₁₁,â₂₂)+|â₁₂|, Gershgorin, unconditional, uniform in h), the chain rebuilt R4-free so the Weil-shaped positivity of A_h never enters it (frame-A trend h^−2.948, ε = 0.052 inside the 0.5 carry window), and the one open object in STANDARD analytic shape for the first time — R1 is a BILINEAR VON MANGOLDT SUM against closed Dirichlet weights; T161's beyond-RH triage applied to LINEAR sums, and for BILINEAR forms the large sieve gives unconditional square-root cancellation — and T170 (BILINEAR.SIEVE, bilinear_sieve_probe.py, SIEVE-RESISTS, 40/40) then ran exactly that toolbox, and the classification completed as a THEOREM: the bilinear form is written down exactly (Â = B − S; det Â = det B − D(B,S) + det S with det S a genuine double von Mangoldt sum against a closed antisymmetric-quadratic wedge kernel; on the reference window the three pieces are 6.18 / −198.8 / 192.6 summing to 4.22e-3 — a five-order cancellation no tool bounding det S alone can see) and it collapses back onto the linear hardness for structural reasons, both theorems — the kernel is the polarisation of the determinant on 2×2 symmetric matrices, rank 3 for every window, every h, every X, so the form IS the rank-3 polynomial S₁₁S₂₂ − S₁₂² in three linear Λ-sums, and the Vaughan Type II blocks have effective rank O(1) because the closed weights see n only through log n; gained: 32 frequencies reduce to THREE linear functionals, by theorem; not gained: the precision, which binds (5.3e-5 per sum at h = 285, sharpening h^−3, an RH yardstick 1.1e3 too coarse); no unconditional route exceeds δ = +0.996 against the target 3.0 — the shortfall is a theorem, not a measurement — and the scramble control localises the arithmetic entirely in the JOINT VALUES of (S₁₁, S₂₂, S₁₂) against the archimedean block (rank-3 and the unconditional bounds unchanged, the truth moves 2.76 in the exponent); as a free by-product the chain becomes R4-free (1 − r₁₂² = det Â/(â₁₁â₂₂) is an identity — the Weil fence is never approached), and R1 is finally CLASSIFIED as a NEAR-DEGENERACY, not a size: an unconditional certificate that two explicit finite Λ-sum vectors become collinear at rate h^−3, beyond the reach of any size-bounding tool. The theorem cores of T168+T169+T170 are promoted as v558 (PRIME.BILINEAR.RANK.01) — and T171 (FINAL.MAP, final_map_probe.py, MAP-COMPLETE, 43/43) then assembled the capstone, and the map is complete: all sixteen links of the reduction chain, from the I5 floor to the final R1 shape, reproduce in ONE connected run on 12 windows over two frames (13 theorems, 3 per-window certificates; the one indefinite deep window excluded out loud — the ladder links carry 11/12, the identity links 12/12), all eight classified no-go routes fail on instances exactly as classified, and the precision ledger closes — the needed joint precision (2.284e-7 at h = 1444, sharpening h^−3) stands 1.2e5× beyond the RH yardstick and 3.1e3× beyond the best unconditional exponent; the file adds ZERO new uniform-in-m statements, which is precisely why it is promotable, and the capstone is load-bearing as v559 (PRIME.PHASE2.CAPSTONE.01): phase 2 is a certified map with one open object — R1, classified as a near-degeneracy, not a size — and T172 (FRAME.BEYOND, frame_beyond_probe.py, PARTIALLY-PORTABLE, 44/44) then tested how far the map carries beyond frame A: 13 of the 16 links transfer UNCHANGED to a gap-blind frame, to ν = 3 and ν = 8, to non-prime-power anchors and to both congruence classes mod 4 (0 broken; the three that shift are exactly the number-carrying ones), the indefiniteness is localised at the sieve horizon rather than at any frame or zone, R1's near-degeneracy persists on all 54 windows while its rate spans h^−1.83 to h^−2.87, and the scramble at fixed Λ-value multiset removes the decay entirely — the collapse belongs to the actual prime-power placement — and T173 (FRAME.RATE, frame_rate_probe.py, DEFICIT-VARIES, 40/40) then showed the demanded rate is itself a frame datum: q = 1 − s is an exact identity with q < 1 on every frame of an eleven-member preregistered family (the h^−3 target was the q = 1 idealisation), the calibrated dimensionless gap functional makes the deficit between demand and delivery THE number of the phase — +0.155 ± 0.102, 1.8× flatter than either side, invariant under the anchor-to-grid rule and the lever split, not yet constant in ν — and no frame closes it (best 2-s.e. upper edge +0.170): frame shopping is over, by numbers, and R1 stands frame-free — the delivered near-degeneracy closes slower than the relative gap it must pay, by a tenth and a half in the exponent. and T174 (CANCEL.IDENTITY, cancellation_identity_probe.py, PARTIAL-CANCEL, 37/37) then exhausted the gauge route by theorems: everything multiplicative cancels exactly — including the entire μᴾ₁ ≈ h⁻² channel, which is exactly absent from the deficit — the one shared channel is certified unconditionally at under five percent of the deficit, and the additive arch/comb mixture provably admits no factorisation (neither term alone has a positive Schur floor on a single one of 84 cells); the bigger result is a direct measurement: on a frame-rule-free rectangle the deficit stands at +0.1111 ± 0.0222 — five sigma from zero, 0.4σ from T173's +0.155 — with its driver identified as comb density per lag cell (log ν and log density collinear at r = −0.921 on legs: T173's ν-driver was the density channel under another name); P6 does not become a theorem (only its multiplicative half) but the deficit no longer needs an invariance argument, and the theorem cores of T172+T173+T174 are promoted as v560 (PRIME.FRAME.DEFICIT.01). and T175 (PHASE.PLACEMENT, phase_placement_probe.py, PHASES-RESIST, 38/38) then measured the placement phases directly: they are real in log R (F = 10.80 against a within-rung-scrambled null at 1.69, +10.3% held out, a composite placebo at 0.81/−0.0275) and — for the first time in the series — CAUSAL under an exact-rebuild intervention (dlog R/dδ = 879, linear over four declared decades), but no phase formula is certifiable: the response is provably non-smooth, because the Schur floor is a near-degeneracy (GAP = 1.6e-9–4.8e-6) and the fitted first harmonic falls a factor ~300 short; the per-anchor heterogeneity dissolved substantially into the error bar (jackknife excess scatter 0.1224 → 0.0000, robust pull 1.41 → 1.22 — an OLS independence assumption, not physics; the 5σ deficit untouched), and the deficit(dens) curve falls monotonically over three decades to a densest reachable bin CONSISTENT WITH ZERO — true zero, power-law approach and low plateau undecidable under this sieve, whose ceiling dens ≤ 361 is a theorem about the caps. And the reverse flow resumed at the far end (T177, CP.INVARIANT, DEGENERATE): the gauge-degree toolkit turned on the four CP ledger rows of the physics side — Tier-1 invariance of the π/3 structure is empty by algebra, two inequivalent Jarlskog-class invariants carry π/3 (neither coupled to the frame), one new exact identity strengthens the E8 channel split (Σρ^d = 4ρ against Σρ^m = 4, the phase channel sheet-blind), and four ledger rows were narrowed with no marker moves — the exact cores promoted as v561. And T176 (DENSE.LIMIT, dense_limit_probe.py, SITS-AT-ZERO, 24/24, landed after T177) then ran the larger sieve — the last decidable measurement of phase 2: ATOM_MAX 1.2e6 → 2.5e7 (factor 20.8; the sieve was never the binding constraint), density ceiling 361 → 6120, two new ratio-4 bins BOTH consistent with zero (+0.0376 ± 0.0410 and +0.1046 ± 0.0881; pooled +0.0496 ± 0.0372) — no sign change, no stabilisation above zero, the plateau window narrowed 3.6× (0.264 → 0.074 at 2σ), all six old bins reproducing at 0.0σ; R1 closed as the Feynman–Hellmann identity (3.4e-6 at the U-minimum; dDEL/DEL carries 84%), the 300× harmonic anomaly explained as a pole (PD survives 2% of the phase period; the two lowest modes already 13.4% apart at δ = 0 — the intervention runs along a crossing that is already there), and one anti-promotion: the bin estimator is ladder-dependent at the 0.20 level, so every quoted deficit must cite its rung ladder. True zero / power-law approach / low plateau stay indistinguishable (the next decade of window would cost a factor one hundred in the sieve) — the phase-2 measurement programme closes as planned, the exact cores are promoted as v562, and the work continues in the classification papers and the backflow lines. Not RH evidence."
          />
          <MeaningLevel
            level="Big if"
            badge="sandbox"
            title="A genuinely new functor"
            body="A compiler functor that is Hecke-translating, Euler-preserving, and ξ-carrying without smuggling ζ from outside. Kills are preregistered; classical named pieces stay classical (Weil 1952, Waldspurger, Cohen 1975, Shintani, Cauchy–Littlewood, Hilbert–Carleman det₂, Connes 1999, Connes–Consani 2021 as context)."
          />
          <MeaningLevel
            level="Dream"
            badge="sandbox"
            title="Riemann Hypothesis"
            body="Not claimed. Not evidenced. I5 ⟺ RH is an equivalence typing of the irreducible core — not a proof claim. Crypto unaffected."
          />
        </ol>
      </DiarySection>

      {/* 26 — Prime shadow */}
      <DiarySection
        id="prime-shadow"
        plain="Within TFPT's narrative the direction of explanation is fixed: the geometry comes first, and primes read it out — even zeta appears as the shadow of a counting function."
        eyebrow="26 · The prime shadow · v625"
        title="Primes enter after the geometry — exactly"
        badge="machine-verified"
      >
        <p>
          An external note asked: what if primes are not the origin but the{" "}
          <em>readout</em> — the shadow of the finished geometry in discrete
          arithmetic? The checkable core of that reading is exact, on the
          compiler&apos;s own objects. The E₈ theta function, computed from
          the glue decomposition, <em>is</em> the Eisenstein series{" "}
          <span className="font-mono text-slate-200">Θ_E8 = E₄</span>: shell
          counts{" "}
          <span className="font-mono text-slate-200">r(2n) = 240·σ₃(n)</span>,
          and the first shell is literally 240 — the root count.
        </p>
        <p>
          The “address space” reading is unique factorization, exactly (shell
          counts factor over coprime addresses; the must-fail control shows
          non-coprime does not). The “independent check channels” reading is
          a theorem: the Hecke operators{" "}
          <span className="font-mono text-slate-200">T_p</span> act with
          eigenvalue{" "}
          <span className="font-mono text-slate-200">1 + p³</span> for every
          prime, they <em>commute</em>, and the compiler&apos;s theta is a
          simultaneous eigenvector of all of them. And{" "}
          <span className="font-mono text-slate-200">
            L(E₄, s) = ζ(s)·ζ(s−3)
          </span>
          : the Riemann zeta function appears as the factorized shadow of the
          E₈ counting function.
        </p>
        <p className="text-slate-400">
          Honest scope: these are classical facts (Jacobi, Hecke) — the
          content is that the compiler&apos;s own objects realize them
          verbatim, fixing the direction of explanation inside TFPT&apos;s
          narrative: geometry first, primes as readout. The bolder framings
          (“primes as compiler eigenfrequencies”, “RH as maximal coherence”)
          stay typed hypotheses, not adopted.
        </p>
      </DiarySection>

      {/* 27 — E8 code */}
      <DiarySection
        id="e8-code"
        plain="The E8 lattice literally is an error-correcting code, and the compiler's own symmetries pick out its bits."
        eyebrow="27 · The error-correcting code · v626 / v638"
        title="E8 is literally a code — and the compiler reads its bits"
        badge="machine-verified"
      >
        <p>
          “E₈ is an error-correcting code” is now a theorem in the suite:
          Construction A on the self-dual extended Hamming code{" "}
          <span className="font-mono text-slate-200">[8,4,4]</span> — the
          Reed–Muller code RM(1,3) — yields E₈ exactly (even, unimodular,
          shells 240/2160), and every single-bit error is exhaustively
          correctable: 16×8 corrupted words, a unique nearest codeword every
          time (v626).
        </p>
        <p>
          v638 then made the dictionary compiler-native instead of decorative.
          Stage 1 killed the naive coordinate placement; among all 30
          placements of the code exactly <em>one</em> (up to the anchor
          orientation) carries both compiler symmetries. On that placement
          the code&apos;s coordinates organise as{" "}
          <strong className="font-medium text-slate-200">
            four μ₄ pairs — one bit per pair
          </strong>
          : the family 3-cycle rotates three of the pairs and fixes the
          fourth, the anchor —{" "}
          <span className="font-mono text-slate-200">
            3 families + 1 anchor
          </span>
          , read off the code. The placement reproduces the v629 root
          censuses verbatim, and the syndrome space factors along the same
          structure.
        </p>
        <p className="text-slate-400">
          Compiler ties, typed: 8 = rank, 4 = d = |μ₄|, 16 = the carrier
          half-spinor. Robustness language on this page now has an exact
          anchor — a code with a decoder, not a metaphor.
        </p>
      </DiarySection>

      {/* 28 — Lorentz congruence / Hodge chamber */}
      <DiarySection
        id="hodge-chamber"
        plain="Two programmes that seemed unrelated were computing on the same lattice all along — and every data window lands in one chamber of it."
        eyebrow="28 · One Lorentz lattice · v624 / v627 / v635–v637"
        title="The prime form and the cover lattice are the same geometry"
        badge="machine-verified"
        visual={<HodgeConeMap />}
      >
        <p>
          The surprise of the third external review (v624): an explicit
          integer matrix P with det −6 satisfies{" "}
          <span className="font-mono text-slate-200">
            Pᵀ J_det P = J_fix
          </span>{" "}
          exactly — the prime-front determinant form (det 2) and the cover
          polarization lattice of the geometry program (det 72) are the{" "}
          <em>same</em> rational quadratic form, the cover an index-6
          sublattice. A genuine new bridge between prime analysis and Hodge
          geometry: the two programs this diary has been running were
          computing on one lattice.
        </p>
        <p>
          v627 measured what that buys: transported through the congruence,{" "}
          <strong className="font-medium text-emerald-200">
            all 67 complete windows land in the positive cone
          </strong>{" "}
          of the cover polarization lattice — on one sheet, det S &gt; 0
          everywhere. v635/v636 then closed the “found matrix” worry: P is
          the unique minimal, operator-compatible congruence in its census
          class, and it is <em>constructed</em> from canonical operator data
          (the null-cone rays are exactly ker C_V and fix C_V), not guessed.
        </p>
        <p className="text-slate-400">
          Honesty, twice: chamber membership is a density-layer statement —
          scrambled combs do not leave the chamber (v582) — so the fine
          C&nbsp;=&nbsp;1 arithmetic lives <em>inside</em> the chamber; and
          v637 closed the tempting follow-up as a preregistered negative: the
          fine Hodge-ray invariants do <em>not</em> predict the C = 1 margin
          window-wise beyond the trivial h-trend.
        </p>
      </DiarySection>

      {/* 29 — ST31 */}
      <DiarySection
        id="sixty-lines"
        plain="The quotient of the lattice is a classical sixty-line reflection group — and a tempting numerical coincidence about its size is killed, not celebrated."
        eyebrow="29 · The sixty lines · v633 / v634"
        title="The μ₄ quotient is a classical reflection group — and a numerology is buried"
        badge="machine-verified"
      >
        <p>
          v629 left a sharp positive residue: the μ₄ clock acts freely on the
          240 E₈ roots with exactly 60 orbits — and 60 = D_start, the
          cascade&apos;s starting value. v633 built the quotient properly:
          the 60 orbits are 60 lines in{" "}
          <span className="font-mono text-slate-200">ℤ[i]⁴</span> carrying a
          hermitian form (the “complex E₈”), and the 60 order-2 unitary
          reflections they define generate a group of order{" "}
          <span className="font-mono text-slate-200">46080</span>.
        </p>
        <p>
          v634 identified it: exactly 60 reflections, invariant degrees{" "}
          <span className="font-mono text-slate-200">(8, 12, 20, 24)</span> —
          the Shephard–Todd group{" "}
          <strong className="font-medium text-slate-200">G31</strong>, pinned
          by fingerprints computed, not cited. The compiler sits inside it
          canonically:{" "}
          <span className="font-mono text-emerald-300">
            σ = c⁴, J = c⁹
          </span>{" "}
          — the family 3-cycle and the μ₄ clock are <em>powers of one
          order-12 element</em>, and the μ₄ center is a power of the clock.
        </p>
        <p className="text-slate-400">
          The kill that keeps it honest: |G31| = 46080 = |W(D₅)|·|W(A₃)| —
          the compiler&apos;s glue factors — looks like destiny and is{" "}
          <strong className="font-medium text-rose-200">numerology</strong>.
          The abstract isomorphism dies on cheap invariants, and stronger:
          W(D₅) embeds in <em>no</em> rank-4 group at all, so the order
          coincidence carries no subgroup structure. Also reproduced one
          level deeper: the compiler clock is <em>not</em> the ζ₁₂-regular
          element of G31 (census 19×12 + 3×4) — the v629 kill stands.
        </p>
      </DiarySection>

      {/* 30 — Live updates */}
      <section
        id="updates"
        aria-labelledby="updates-heading"
        className="scroll-mt-24 border-t border-slate-800/60 py-12 sm:py-16"
      >
        <div className="mx-auto max-w-5xl px-4 sm:px-6 lg:px-8">
          <div className="flex flex-wrap items-center gap-2">
            <span className="font-mono text-[10px] uppercase tracking-[0.18em] text-slate-500">
              30 · Live updates
            </span>
            <StatusBadge badge="sandbox" />
          </div>
          <h2
            id="updates-heading"
            className="mt-3 font-serif text-2xl font-semibold text-slate-50 sm:text-3xl md:text-4xl"
          >
            Feed — one entry per completed agent run
          </h2>
          <div className="mt-8">
            <DiaryFeed />
          </div>
        </div>
      </section>

      <footer className="border-t border-slate-800/60 py-10">
        <div className="mx-auto max-w-5xl px-4 text-sm text-slate-500 sm:px-6 lg:px-8">
          Numbers and verdicts are taken from{" "}
          <code className="font-mono text-slate-400">experiments/next.txt</code>{" "}
          (2026-07-23 … 2026-08-02 diary and backflow rounds) and the promoted
          modules{" "}
          <code className="font-mono text-slate-400">
            verification/v535_*.py
          </code>{" "}
          …{" "}
          <code className="font-mono text-slate-400">v648_sign_uncertainty.py</code>.
          Exploratory probes live under{" "}
          <code className="font-mono text-slate-400">
            experiments/tfpt-discovery/
          </code>
          .
        </div>
      </footer>
    </>
  );
}

function BigPictureSection() {
  return (
    <section
      id="big-picture"
      aria-labelledby="big-picture-heading"
      className="scroll-mt-24 border-t border-slate-800/60 py-12 sm:py-16"
    >
      <div className="mx-auto max-w-5xl px-4 sm:px-6 lg:px-8">
        <div className="flex flex-wrap items-center gap-2">
          <span className="font-mono text-[10px] uppercase tracking-[0.18em] text-slate-500">
            00 · Start here
          </span>
          <StatusBadge badge="sandbox" />
        </div>
        <h2
          id="big-picture-heading"
          className="mt-3 font-serif text-2xl font-semibold leading-tight text-slate-50 sm:text-3xl md:text-4xl"
        >
          The big picture, in plain words
        </h2>
        <p className="mt-4 max-w-3xl text-base leading-relaxed text-slate-300">
          The rest of this page is a diary, written as the work happened, in the
          language of the work. This section says the same thing in ordinary
          words. Every number below is copied from the diary. Every limit is
          named.
        </p>

        <div className="mt-8 grid gap-4 lg:grid-cols-2">
          <article className="rounded-2xl border border-sky-400/30 bg-sky-500/[0.07] p-5 sm:p-6">
            <p className="font-mono text-[10px] uppercase tracking-[0.18em] text-sky-300/90">
              What we are trying to do
            </p>
            <h3 className="mt-2 font-serif text-xl text-sky-50">
              A relay that has to work forever
            </h3>
            <div className="mt-3 space-y-3 text-sm leading-relaxed text-slate-300">
              <p>
                There is one quantity in this story that must never go negative.
                It is assembled out of the prime numbers, one prime at a time.
              </p>
              <p>
                Picture a bridge built span by span. Each new span has to carry
                the load the previous one can no longer hold. Each prime is a
                span. The handover from one span to the next is the whole
                question: does it keep working, prime after prime, without end?
              </p>
              <p>
                If it does, a classical statement called{" "}
                <strong className="font-medium text-slate-100">
                  Weil positivity
                </strong>{" "}
                holds — and that statement is equivalent to the Riemann
                Hypothesis. So this is where the difficulty of the Riemann
                question actually sits, in the form the diary can touch. We do
                not prove it. We test the relay span by span, and we write down
                precisely what a proof would still need.
              </p>
            </div>
          </article>

          <article className="rounded-2xl border border-emerald-400/30 bg-emerald-500/[0.07] p-5 sm:p-6">
            <p className="font-mono text-[10px] uppercase tracking-[0.18em] text-emerald-300/90">
              What actually works today
            </p>
            <h3 className="mt-2 font-serif text-xl text-emerald-50">
              The mechanism has never failed a step
            </h3>
            <div className="mt-3 space-y-3 text-sm leading-relaxed text-slate-300">
              <p>
                Every handover that has been checked, checks out. 117 of 117 in
                the deep ladder{" "}
                <span className="font-mono text-amber-200">[T111]</span>, 400 of
                400 rungs of the nested ladder{" "}
                <span className="font-mono text-amber-200">[T124]</span>, and a
                single certified step as deep as{" "}
                <strong className="font-medium text-emerald-100">
                  zone 155,921
                </strong>{" "}
                <span className="font-mono text-amber-200">[T115]</span>.
              </p>
              <p>
                When a chain of steps stops, it stops because the computer runs
                out of budget — never because a step failed.
              </p>
              <p>
                It did not always look like that. For a while there seemed to be
                a wall near zone 462{" "}
                <span className="font-mono text-amber-200">[T111]</span>. Then
                T113 showed the wall was measuring our own grid rather than the
                mathematics, and T114 removed the division that produced it. The
                wall dissolved, and eleven steps opened up beyond it at once.
              </p>
            </div>
          </article>
        </div>

        <div className="mt-8 grid gap-4 lg:grid-cols-2 lg:items-start">
          <div className="space-y-3 text-base leading-relaxed text-slate-300">
            <p className="font-mono text-[10px] uppercase tracking-[0.18em] text-slate-500">
              The great compression
            </p>
            <h3 className="font-serif text-xl text-slate-50">
              One huge inequality, squeezed down to two small ones
            </h3>
            <p className="text-sm">
              At T104 the thing left to prove was a{" "}
              <strong className="font-medium text-slate-100">
                matrix inequality
              </strong>
              : a statement about every direction of a large window at once.
              Twenty-one parts later the load-bearing part of it is one sign and
              one accounting convention.
            </p>
            <p className="text-sm">
              Nothing was assumed away on the route. At each stage the older,
              bigger statement was shown to <em>follow</em> from the smaller one
              — and several times a step of our own was refuted and replaced,
              which is why the list is shorter and not longer.
            </p>
            <p className="text-sm text-slate-400">
              What is left today: on the load-bearing spine, one <em>sign</em>{" "}
              that the induction already carries on its own plus one declared{" "}
              <strong className="font-medium text-slate-200">
                accounting convention
              </strong>{" "}
              <span className="font-mono text-amber-200">[T124, T125]</span>. The{" "}
              <strong className="font-medium text-slate-200">
                Harnack pair
              </strong>{" "}
              — two statements that hold on every window measured so far —
              survives as a second, independent route and no longer holds the
              spine up. Neither is proved. Both now have classical addresses.
            </p>
          </div>
          <ReductionCascade />
        </div>

        <aside className="mt-8 rounded-2xl border border-rose-400/35 bg-rose-500/[0.08] p-5 sm:p-6">
          <p className="font-mono text-[10px] uppercase tracking-[0.18em] text-rose-300">
            The honest box — what this is not
          </p>
          <ul className="mt-3 space-y-2.5 text-sm leading-relaxed text-slate-300">
            <li>
              <strong className="font-medium text-rose-100">
                No zeta input anywhere.
              </strong>{" "}
              The chain is classical numerical analysis from end to end. The
              zeta function is never fed in, so nothing here can be a
              circular re-derivation of it.
            </li>
            <li>
              <strong className="font-medium text-rose-100">
                No infinite statement.
              </strong>{" "}
              Every certificate is a finite computation on a finite window. An
              infinite claim would need the two open statements — that is the
              whole point of naming them.
            </li>
            <li>
              <strong className="font-medium text-rose-100">
                A measured ladder is not all zones.
              </strong>{" "}
              Where this page says “certified”, it means the zones that were
              actually run. Trends beyond them are labelled as fits or
              extrapolations, never as results.
            </li>
            <li>
              <strong className="font-medium text-rose-100">
                Sandbox, not verified physics.
              </strong>{" "}
              This entire diary is exploratory. Only the promoted modules{" "}
              <Link
                href="/verification"
                className="font-mono text-emerald-300 underline decoration-emerald-400/40 underline-offset-2 hover:text-emerald-200"
              >
                v535–v673 (this front)
              </Link>{" "}
              are machine-verified and cited in the papers. Sandbox probes never
              move a status marker.
            </li>
            <li>
              <strong className="font-medium text-rose-100">
                Not RH, and not almost-RH.
              </strong>{" "}
              “I5 ⟺ RH” is a <em>typing</em> of what the remaining core is
              equivalent to — not a proof claim and not a measure of distance.
              No progress toward the Riemann Hypothesis is claimed. Cryptography
              is unaffected.
            </li>
          </ul>
        </aside>

        <div className="mt-8">
          <p className="font-mono text-[10px] uppercase tracking-[0.18em] text-slate-500">
            Where it stands now
          </p>
          <dl className="mt-3 grid grid-cols-2 gap-3 sm:grid-cols-4">
            <BigPictureStat
              term="240"
              desc="agent runs in the diary — series complete at 125 parts, phase 2's measurement programme closed, backflow rounds ongoing"
              tone="sky"
            />
            <BigPictureStat
              term="5818"
              desc="sandbox checks across 261 probes, all passing"
              tone="amber"
            />
            <BigPictureStat
              term="v535–v673"
              desc="machine-verified modules of this front, inside the 667-script suite (all green)"
              tone="emerald"
            />
            <BigPictureStat
              term="96.2%"
              desc="of the load-bearing spine is an identity or a Cholesky certificate (T125 finale)"
              tone="violet"
            />
          </dl>
          <p className="mt-4 max-w-3xl text-sm leading-relaxed text-slate-400">
            The chain is written as a{" "}
            <strong className="font-medium text-slate-200">
              conditional theorem
            </strong>
            : everything after the word “suppose” is proved or certified. The
            finale{" "}
            <span className="font-mono text-slate-300">
              grand_assembly_probe.py
            </span>{" "}
            assembled it end to end on 52 zones and changed what has to be
            supposed:{" "}
            <strong className="font-medium text-slate-200">
              the series is complete, and the chain&apos;s spine needs no Harnack
              pair
            </strong>{" "}
            — that pair survives only as a second, independent route, and the one
            hypothesis the spine still carries is a declared accounting
            convention. What is missing for any infinite statement is uniformity
            in the zone index, not size — and that is now the program: the
            series is complete, and a second phase (T126+) attacks the two
            remaining genuinely new inequalities, the direction lemma and the
            zone-uniform seam floor. Sandbox; not RH evidence.
          </p>
        </div>
      </div>
    </section>
  );
}

const BIG_PICTURE_TONE = {
  sky: "border-sky-400/30 bg-sky-500/[0.08] text-sky-100",
  amber: "border-amber-400/30 bg-amber-500/[0.08] text-amber-100",
  emerald: "border-emerald-400/30 bg-emerald-500/[0.08] text-emerald-100",
  violet: "border-violet-400/30 bg-violet-500/[0.08] text-violet-100",
} as const;

function BigPictureStat({
  term,
  desc,
  tone,
}: {
  term: string;
  desc: string;
  tone: keyof typeof BIG_PICTURE_TONE;
}) {
  return (
    <div className={`rounded-xl border px-3 py-3 ${BIG_PICTURE_TONE[tone]}`}>
      <dt className="font-mono text-xl leading-none">{term}</dt>
      <dd className="mt-1.5 text-[11px] leading-tight text-slate-400">
        {desc}
      </dd>
    </div>
  );
}

function ProgramStatusCallout() {
  return (
    <aside className="rounded-2xl border border-violet-400/35 bg-violet-500/10 p-5 sm:p-6">
      <p className="font-mono text-[10px] uppercase tracking-[0.18em] text-violet-300/90">
        Strongest current sentence
      </p>
      <p className="mt-3 font-serif text-lg leading-relaxed text-violet-50 sm:text-xl">
        After the T102–T176 sprint (its identity cores promoted module by
        module), what remains TFPT-specific is exactly one object: I5 — and the
        induction that would carry it is now
        a conditional theorem, assembled end to end on 52 zones (T125), whose
        load-bearing spine is 96.2% identity or Cholesky certificate, whose only
        hypothesis is a declared accounting convention, and whose seam
        architecture is finished (T126): both remaining inequalities are
        proof-shaped (T127), three of the four resulting points stand at their
        preregistered bars (T128 · THREE-OF-FOUR), the kappa law that T128
        preregistered falls once on 331 fresh transports while the theorem
        underneath it stands (T129 · KAPPA-WILD), the graded-to-uniform
        bridge stands as an identity — 84 pairs, zero overshoot, both deep
        seams carried to positive fine floors — while the curvature bound is
        reduced to a uniform bound on one exponent (T130 · ONE-OF-TWO), and
        the self-supply loop is built one number short of closed: the
        epsilon-to-floor sandwich and sign constancy via Perron–Frobenius
        are theorems, and M25 reduces to positivity of the pole-free
        section with nine decades of slack (T131 · SUPPLY-PARTIAL). The phase
        then reversed direction: that identity block is now load-bearing v542
        (nine per-instance identities and theorems, nothing uniform), the
        Beurling–Deny triad became an operator discriminator for the seam DtN
        (T132 · SPECTRUM-ONLY), and a certificate audit hardened v379 to an
        exact positive-mixture Gram certificate (T133 · MIXED). T134 (PARTIAL)
        then closed the existence half of the pole-free floor — all Cholesky
        pivots positive on 79/79 windows, every cheap lower-bound route
        failing by sign, not size, with the one surviving opening named (an
        M-matrix question) — and T135 (BOUNDED-STATE) found the bounded
        faithful seam state the Weil window provably lacks (m_cert = 12, from
        the pre-declared set, with the honest caveats stated). T136
        (ONE-CARRIES) then closed one item of that M-matrix question outright —
        Varga&apos;s regular splitting makes the Jacobi radius an identity and
        Collatz–Wielandt bounds it a priori, sharp to 1.00–1.03 on 900/900 — with
        the exact bookkeeping putting the whole degradation in the margin and M17
        closing negatively; and T137 (BOTH-RESIST) made the long-lag support an
        explicit arithmetic stripe set and certified the entire absolute-value
        envelope family DEAD (ρ(|E|) ≥ 1.32 from below on 35/35), leaving one
        named residue: a sign-preserving bound. Thirteen of those statements are
        now load-bearing as v543 and v544.
      </p>
      <ul className="mt-4 space-y-2 text-xs leading-relaxed text-slate-400">
        <li>
          <strong className="font-medium text-slate-200">
            The mechanism never failed.
          </strong>{" "}
          Every certified handover step holds — 117/117 on the deep ladder at
          retention 1.000000 (T111), 400/400 rungs of the nested ladder (T124).
          Chains stop at the compute budget, never at a step (T114 onward).
        </li>
        <li>
          <strong className="font-medium text-slate-200">
            Depth reached, certified.
          </strong>{" "}
          Margin-free steps up to n = 1331 (T114) and a single certified step at
          n = 155,921 (T115); a boundary march to 1.35 million cells at flat cost
          (T116, declared a cost-geometry demonstration, not a Weil certificate).
        </li>
        <li>
          <strong className="font-medium text-slate-200">
            The wall was a ruler.
          </strong>{" "}
          The apparent margin wall near zone 462 (T111) is currency-invariant
          (T113) but measures the discretization, not the spectrum — and it
          dissolves once the artifact division is removed (T114).
        </li>
        <li>
          <strong className="font-medium text-slate-200">
            What is left, named.
          </strong>{" "}
          On the spine: one sign the coarse-to-fine induction already carries
          plus one declared accounting convention (T124/T125). The Harnack pair
          (T120/T121, unconditional |R − 1| ≤ 0.047 plus a section statement) is
          a second, independent route and no longer load-bearing. Honestly booked
          against it: the solution-free version of the rung fails, the
          unconditional D₀ criterion covers only the first three zones (T120),
          and uniformity in the zone index is open (T125).
        </li>
        <li>
          <strong className="font-medium text-slate-200">
            Series balance — complete.
          </strong>{" "}
          125 parts, 3139/3139 sandbox checks, seven promoted modules; the finale{" "}
          <span className="font-mono text-slate-300">
            grand_assembly_probe.py
          </span>{" "}
          (ASSEMBLY-GREEN, 34/34) assembled the chain end to end — 52 zones, 430
          completed Cholesky certificates, 24 refuted routes and a seven-station
          cascade booked. The mandate T ≤ 125 is fulfilled; the next step is a
          decision about extraction, not another probe. The full sprint is told
          in{" "}
          <a
            href="#compression"
            className="text-sky-300 underline decoration-sky-400/30 underline-offset-2 hover:text-sky-200"
          >
            sections 22–24
          </a>
          . Not almost-RH. This is not RH evidence.
        </li>
        <li>
          <strong className="font-medium text-slate-200">
            Phase 2: the full proof.
          </strong>{" "}
          The first post-series probe{" "}
          <span className="font-mono text-slate-300">
            uniformity_seams_probe.py
          </span>{" "}
          (T126 · UNIFORMITY.SEAMS — SEAMS-CERTIFIED, 31/31) finishes the seam
          architecture: the partition question dissolves into a monotone
          free-resolution construction (re-gridding only at record-small gaps,
          all 12 real seams certified), the continuum step closes via an exact
          fractional-Dirichlet identity, and the full-proof map contains
          exactly two genuinely new inequalities — a direction lemma (U3) and a
          zone-uniform seam floor (U5). T127 ·{" "}
          <span className="font-mono text-slate-300">
            two_inequalities_probe.py
          </span>{" "}
          (TWO.INEQUALITIES — BOTH-SHAPED, 28/28) dissects both, and both
          changed: the U5 separating variable is pure interval geometry (not
          arithmetic in the atoms), U5-as-stated is refuted and replaced by a
          band covering 99.26% of real seams with the 8 exceptions enumerated
          (Mersenne/Fermat pairs and one twin), and U3 collapses to a coarse
          floor via the harmonic-mean identity. T128 (TEML —{" "}
          <span className="font-mono text-slate-300">teml_probe.py</span>,
          THREE-OF-FOUR, 27/27) works the resulting four-point list cheapest
          first, and three of the four points stand at their preregistered
          bars: the exception list is derived and closed (the affordable seams
          certified, the open five each lacking exactly one Cholesky beyond
          the cap), the retention bound is exact bookkeeping, and the
          boundary-layer exclusion is now a proof with an 11.6× floor margin —
          while the fourth point, the protrusion concentration kappa, missed
          its own bar by 3.6%, systematically in the ratio, which is what
          distinguishes a bound from a fit. T129 (KAPPA.DEEP.SEAMS,{" "}
          <span className="font-mono text-slate-300">
            kappa_deep_seams_probe.py
          </span>
          , KAPPA-WILD, 28/28) is the most productive break of the phase: the
          fitted kappa law C* = 0.3506 — frozen from T128&apos;s calibration
          groups, labelled a fit — falls once on 331 fresh transports (the bar
          not moved; the bare bar 2 would have fallen 78 times), but the
          structure underneath is solved with identities: the flat border
          vector sits at kappa = 1 exactly, a linear density profile at
          kappa = 2 exactly, everything above is curvature, and the curvature
          chain is a per-transport theorem with an explicit constant holding
          on all 436 transports — with a matrix-free depth block to 17,669
          fine cells (12× the cap) showing zero violations. The two affordable
          deep seams (n = 127, 256) carry complete certificates on the graded
          two-scale space, honestly downgraded to GRADED-CERTIFIED: a measured
          8% false-positive rate of the compression, declared before the
          results. A narrow verification module on the five identity/theorem
          pieces is ready in principle; two irreducibles remain (the word
          &ldquo;for all&rdquo;, the RH address). T130 (CURVATURE.BRIDGE,{" "}
          <span className="font-mono text-slate-300">
            curvature_bridge_probe.py
          </span>
          , ONE-OF-TWO, 30/30) then attacked the two named pieces and exactly
          one stands: the graded-to-uniform bridge stands as an identity —
          the matrix-form Céa/Strang defect S_graded − S_uniform = RᵀX⁻¹R
          (verified to 1.3e-13, matrix-free) reproduces the directly computed
          uniform floor on 84 pairs at 21 zones to 2.5e-11 relative with zero
          overshoot, explains the 8% false positives completely (16/84 at
          bracket level, 0/84 at floor level — they sat in a bracket the
          bridge makes unnecessary), and carries both deep seams (n = 127,
          256) to positive fine floors 0.2255 and 0.2614 at 1.7× and 3.8×
          the factorization cap, conditional on one measured number per grid
          (M25, new and named — with the elegant candidate that the
          chain&apos;s own telescope floor supplies it) — while the curvature
          bound honestly broke its frozen shape band on 13/545 disjoint test
          transports and is reduced to a zone-uniform bound on one exponent
          (M22). T131 (SELF.SUPPLY,{" "}
          <span className="font-mono text-slate-300">
            self_supply_probe.py
          </span>
          , SUPPLY-PARTIAL, 25/25) then tested exactly that candidate and
          built the self-supply loop, one number short of closed: the
          epsilon-to-floor secular sandwich is a theorem — sharp to a
          factor ~1.3, with the sign half an equivalence (Albert 1969) —
          that replaces the Lanczos estimate on all 84 bridge pairs with
          zero brackets lost and exposes that the old Ritz value
          overestimated λ_min(X) by up to 7.9×; sign constancy became a
          theorem via Perron–Frobenius on S⁻¹ (entrywise positive
          575/575); the one-hump honestly broke at depth (530/575, the
          certified sag majorant replaces it) and S* rose to 1.8472 over
          its frozen 1.1926; M17 assembled but vacuous, with the loss
          attributed to the U-metric mismatch; M25 is REDUCED-TO-POLE-FREE
          with nine decades of slack. The identity block is since promoted as
          load-bearing v542, and the reverse-flow probes T132 (BD.SEAM) and
          T133 (CERT.FLOOR) carried the certified toolkit back to the theory
          side. T134 (POLE.FREE.FLOOR,{" "}
          <span className="font-mono text-slate-300">
            pole_free_floor_probe.py
          </span>
          , PARTIAL, 21/21) then attacked that floor and closed its
          existence half: every Cholesky pivot of the pole-free form is
          positive on 79/79 windows, but all six cheap lower-bound routes
          fail by sign, not size — the anatomy names the one surviving
          opening (an M-matrix question: the lumped Stieltjes comparison
          S_B = S + L_Δ is certified on 900/900 blocks), and the whitening
          honestly corrects T131&apos;s diagnosis. T135 (COMB.COMPRESS,{" "}
          <span className="font-mono text-slate-300">
            comb_compressibility_probe.py
          </span>
          , BOUNDED-STATE, 13/13) ported the T116 Riccati machinery to the
          seam DtN and found a bounded faithful state (m_cert = 12, error
          falling out to h = 1e5) where the Weil window provably has none —
          driver weight summability, value partly circular, QEC.SEAM.01 not
          advanced. T136 (M.MATRIX.PAIR,{" "}
          <span className="font-mono text-slate-300">
            m_matrix_pair_probe.py
          </span>
          , ONE-CARRIES, 30/30) took that M-matrix question and closed one of
          its three items outright: the classical row-sum criterion is vacuous
          at the border (900/900 negative row sums), but Varga&apos;s regular
          splitting makes ρ(J) = τ/(1+τ) an identity and the Collatz–Wielandt
          bound at the anchor vector reproduces the measured radius to
          1.00–1.03, flat in D and in the zone on 900/900 — while the exact
          three-way split λ_min ~ D^−0.56 × D^2.72 × D^0.12 puts the whole
          degradation in the margin and M17 closes negatively (the bad subspace
          is delocalised, its dimension tracking the atom count). T137
          (LONG.LAG,{" "}
          <span className="font-mono text-slate-300">long_lag_probe.py</span>,
          BOTH-RESIST, 22/22) then made the perturbation explicit — the support
          is a comb stripe set at the prime-power atoms, each stripe a perfect
          matching, the amplitude certified — and killed a whole family by
          certificate: ρ(|E|) ≥ 1.32 from below on 35/35 windows, so no
          Gershgorin bound, no row-sum norm and no rescaling can ever supply the
          floor; the Green function is bracketed two-sidedly, and the residue is
          one named object, a sign-preserving bound. Thirteen mature statements
          of both parts are promoted as v543 and v544. T138 (SIGN.COMPENSATION,{" "}
          <span className="font-mono text-slate-300">
            sign_compensation_probe.py
          </span>
          , PAIR-EXACT, 26/26) kept the signs: the coupling sign follows the
          interval geometry (nested positive, disjoint mostly negative), the
          m-paired Neumann certificate removes the arithmetic wall on all 77
          dead blocks (pool 563 → 875/900), and the margin question returns
          one level down as ρ(W_S). T139 (GREEN.DECAY,{" "}
          <span className="font-mono text-slate-300">green_decay_probe.py</span>
          , DENSE-RESISTS, 30/30) refuted the classical decay lemma at its
          hypothesis — arithmetically: the certified kernel envelope decays
          like d^−0.15..−0.04 where Jaffard needs p &gt; 1 — while deriving
          T138&apos;s sign law from one exact telescoping identity and killing
          the never-truncating layer series from below; the core shrinks to one
          signed inequality at stripe distance b ≤ 16. T140 (SIGNED.BAND,{" "}
          <span className="font-mono text-slate-300">signed_band_probe.py</span>
          , FINITE-CORE, 31/31) attacked exactly that inequality: the telescope
          identity lifts to the form level (Gram = CHCᵀ exactly, rank ≤ h−1),
          the spectral radius is the exact finite-core eigenvalue
          ρ(W) = λ_max(K^½HK^½) of a closed-geometry coverage kernel times a
          mass-plus-Dirichlet form, the checkerboard split replaces the O(nb)
          Weyl steps by three D-independent ones, the additive band reading is
          closed negative from below at every b ≤ 16 (rank inflation
          1.2–10.1× is the mechanism), and all the D-dependence sits in the
          geometry — the remaining ingredient is a zone-uniform discrete Hardy
          inequality. T141 (DISCRETE.HARDY,{" "}
          <span className="font-mono text-slate-300">
            discrete_hardy_probe.py
          </span>
          , HARDY-RESISTS, 22/22) attacked exactly that ingredient: four exact
          identities put it in classical two-weight shape, but the certified
          constant is not zone-uniform (D^−0.366 ± 0.036 against a bar of
          0.25) while the exact object it bounds is (D^−0.229 ± 0.007) — the
          growth sits in the diagonal profile — the additive shape is dead as
          a shape at its own exact Weyl floor (1.694–3.855× the target), and
          the joint shape fails at the normalisation alone (Ω = 20.71–2723.99).
          The rest collapses to one closed conductance profile (R1b); the
          identity blocks of T140 and T141 are promoted as v545. T142
          (CONDUCTANCE.PROFILE,{" "}
          <span className="font-mono text-slate-300">
            conductance_profile_probe.py
          </span>
          , PROFILE-RESISTS, 24/24) then constructed that profile instead of
          guessing it: the capacity decomposition K⁻¹ = DᵀJ⁻¹D + xxᵀ/cap
          exhibits the optimal Hardy weight exactly — Ω = 1 EXACTLY by a
          projection identity, against T141&apos;s guessed 20.7–2724 — so the
          weight was never a free choice but geometry; the certified chain
          misses the target by a constant factor 2.27–2.45, flat in D, and
          the rank ladder closes the entire comparison path (no truncation
          below r*/m = 0.995 ever crosses 1) — no comparison argument can
          deliver D-uniformity, so the next move is the sharp
          capacity-Rayleigh route (R1c). T143 (SHARP.CAPACITY,{" "}
          <span className="font-mono text-slate-300">
            sharp_capacity_probe.py
          </span>
          , SHARP-CARRIES, 24/24) then ran exactly that route, and it
          carries: the exact capacity-Rayleigh form is an identity on all 26
          windows, Maz&apos;ya&apos;s capacity criterion applied to the gap
          form lands inside its window [1/4, 1] (Φ_sup·λ = 0.5438–0.6457)
          with a zone-uniform loss factor (D^−0.048 ± 0.010), the supremum
          lives on closed families — intervals in node coordinates,
          certified by full enumeration on the small border blocks — and
          D-uniformity is reduced to ONE named inequality, cap_E(A) ≥
          |A|·λ₀/c₀ with an absolute c₀ (a non-Markovian Maz&apos;ya
          capacity bound, pointing at Muckenhoupt&apos;s 1972 two-weight
          calculus). T144 (CAPACITY.INEQUALITY,{" "}
          <span className="font-mono text-slate-300">
            capacity_inequality_probe.py
          </span>
          , INTERVAL-CARRIES, 31/31) then ran the interval route at exactly
          that inequality, and it carries: the interval class is exhausted
          exactly (11,390,676 intervals via a Cholesky prefix-sum identity),
          the closed two-weight sup B_res·λ̂ = 0.6694–0.7813 lands inside
          the Maz&apos;ya window flat in D (D^0.013 ± 0.005), the family
          restriction falls entirely (a max-density-subgraph bound covers
          all 2^m sets — Ψ_all·λ̂ = 0.7399–0.8515, the flattest number of
          the probe), the Markov perturbation route is certified dead, and
          the certified chain λ ≥ 1/(c₀·κ_up·c_glob·B_res) has exactly ONE
          unproven input: the absolute Maz&apos;ya constant c₀, whose
          sharpest shape S1′ is a Muckenhoupt-type hypothesis. T145
          (MAZYA.PROOF,{" "}
          <span className="font-mono text-slate-300">
            mazya_proof_probe.py
          </span>
          , ONE-STEP-MISSING, 33/33) then ran the proof attempt itself, and
          the transcription works: M1–M3 are theorems, the Markov property
          sits in exactly one line (M4), that line splits — its mass half
          is a theorem and dominates, σ_tot = 0.2145–0.4425 &lt; 1
          everywhere — so c₀ becomes explicit (best 2.248–4.227, flat at
          D^0.028 ± 0.017) and S1′ is certified per window as
          cert_λ_max(R) ≤ c₀·Ψ on 64/64 windows; an explicit no-go proves
          a level-profile hypothesis is necessary, S2′ retires, R5 is
          downgraded, all 24 border blocks close, and the mathematical rest
          is ONE named lemma — L1, an a-priori delocalization bound for the
          minimiser&apos;s layer-cake constant; the identity spine of
          T142–T145 is promoted as v546. T146 (LEVEL.LEMMA,{" "}
          <span className="font-mono text-slate-300">
            level_lemma_probe.py
          </span>
          , LEVEL-CARRIES, 30/30) then ran the proof attempt for exactly
          L1, and the level lemma STANDS: a chain of theorems and certified
          window inequalities with no step reading the minimiser — the
          resolvent identity ψ = λRψ is the lever, Davis–Kahan is
          instrumented and discarded, the cake base is free
          (Maz&apos;ya&apos;s dyadic 8 falls to 2), and
          c₀^ap = 2β²·Γ·min(1,Γ₁) + ε = 3.9042–4.8488 on 64/64 windows —
          at the size of Maz&apos;ya&apos;s classical Dirichlet value 4 —
          closes the chain at loss factor 0.0422–0.1586; the a-priori core
          is promoted as v547. T147 (GREEN.ASYMPTOTIC,{" "}
          <span className="font-mono text-slate-300">
            green_asymptotic_probe.py
          </span>
          , ASYMPTOTIC-SHAPED, 39/39) then made the remaining all-D
          asymptotics an exact identity — Γ = √Q★ · Sw, a purely spectral
          factor times a purely geometric one, both certified on the
          surface, the mechanism named (Szegő/Widom, Q_B ≤ 2|B| holding) —
          and T148 (SZEGO.BOTTOM,{" "}
          <span className="font-mono text-slate-300">
            szego_bottom_probe.py
          </span>
          , ONE-INPUT-MISSING, 31/31) dissected the lifting statement: Sw
          closes on the surface via an LDLᵀ layer-cake certificate
          (≤ 1.9587, flat), the Toeplitz symbol is refuted at 0 (an honest
          negative: positive definiteness is a finite-section effect), the
          hypothesis of the lift is isolated to one named scalar —
          TV(log Λ), the roughness and not the conditioning — and the
          certified chain delivers 8.4–15.9% of the true gap end to end;
          the identity/certificate core of T147+T148 is promoted as v548.
          T149 (WEIGHT.SMOOTHING,{" "}
          <span className="font-mono text-slate-300">
            weight_smoothing_probe.py
          </span>
          , PARTIAL-SMOOTHING, 30/30) then attacked exactly that input
          through the gauge freedom: the whitening diagonal is a free
          gauge, and the constant gauge eliminates the blocking scalar
          exactly — TV(log Λ̃) = 0.0000 at a certified sandwich price
          (σ ≤ 5.5789, κ̃_up ≤ 2.3146) — which refutes the hypothesis
          itself: gauges that kill the flutter move ν_L̃ by at most 0.9%,
          so the roughness sits in the form, not in the multiplier; the
          missing input relocates to the deep modes of the pure
          Toeplitz-minus-Hankel section (ladder form ν_k ≤ C·k² with
          m-free C). T150 (MODE.LADDER,{" "}
          <span className="font-mono text-slate-300">
            mode_ladder_probe.py
          </span>
          , ONE-TERM-MISSING, 36/36) then ran exactly that ladder and named
          the mechanism — parity: the form is exactly the compression of
          the full symmetric Toeplitz section onto its antisymmetric parity
          sector, and the sign-changing symbol&apos;s entire negative
          inertia sits in the EVEN sector (72/72, LDLᵀ) — T148&apos;s
          honest negative explained rather than worked around; the flutter
          amplitude became a certified form functional with a closed atom
          budget, and the atoms are co-responsible for the positivity (the
          archimedean section alone is indefinite). One gate remains — the
          ladder constant still grows (C ≤ 43.391 = 13.81π, x^0.258
          against the bar 0.25), the gap from π exactly the arithmetic
          leakage of the bottom mode; the identity/certificate core of
          T149+T150 is promoted as v549. T151 (ODD.LADDER,{" "}
          <span className="font-mono text-slate-300">
            odd_ladder_probe.py
          </span>
          , ODD-CARRIES, 27/27) then closed exactly that inequality by
          rerouting off the ladder: the odd grid steps over the
          symbol&apos;s negative window (θ_c/θ_1 = 0.328–0.407 on 72/72 —
          positivity is a grid fact), the bottom spectrum is certified
          against the parity Laplacian (λ_k ≤ S·μᴾ_k, S = 1.10–2.39,
          LDLᵀ), and a discrete Sobolev step at the odd sector&apos;s
          virtual node yields a per-mode bound LINEAR in k with a
          non-growing constant (C_S = 11.51–19.57, x^0.020) where the
          quadratic ladder grew; two of T150&apos;s three rest items are
          settled (the archimedean minimum closed and attained exactly;
          the symbol route a computed dead end — the local model is a
          matrix statement, not a symbol one), and one scalar remains a
          fit: the bottom pencil ratio R = K_bot/κ = 3.36–9.71 (flat,
          x^0.037), exactly where the T145 no-go breaks (x^1.986); the
          identity/certificate core of T151 is promoted as v550. T152
          (PENCIL.RATIO,{" "}
          <span className="font-mono text-slate-300">
            pencil_ratio_probe.py
          </span>
          , ONE-TERM-MISSING, 37/37) then attacked exactly that scalar:
          both hoped-for archimedean gifts are refuted — the arch section
          is itself negative in the odd sector (O(−m²) in pencil
          normalisation) and the atom part too, so positivity is a
          cancellation between geometry and arithmetic — while the Schur
          two-block criterion with a fixed 16-mode block certifies the
          floor κ ≥ 0.225–0.250 flat on every window, leaving R ≤ K_bot/t
          = 4.41–8.47 certified on both ends and one term missing (an
          m-free ceiling on K_bot); the Ψ map locates 89–92% of the
          dominant loss on the eight modes the certified ladder already
          controls (a 3.5–5.3× reserve). T153 (PSI.LADDER,{" "}
          <span className="font-mono text-slate-300">
            psi_ladder_probe.py
          </span>
          , REBUILD-RESISTS, 33/33) then carried out exactly that rebuild
          and refuted it by its own factor: Psi is pinned to within
          1.20–1.30× by the first parity sine (a₁ = 8/π²) — the reserve
          never existed — while the theorem-grade collapse closes Psi as a
          term, the archimedean part turns out positive on the bulk parity
          block (after an m-free peeling of at most 8 modes), and the two
          remaining terms — the block inequality and a Green/alignment
          estimate of two inverse-iteration columns, which already yields a
          flat ceiling on a fixed-size certificate — lose to the same
          missing object. T154 (GREEN.ALIGN,{" "}
          <span className="font-mono text-slate-300">
            green_align_probe.py
          </span>
          , ALIGN-RESISTS, 29/29) then attacked exactly that estimate, and
          the ceiling closes exactly at fixed size: the sixteen-column
          certificate (the parity sines plus one Green step A⁻¹L_P·t_k)
          needs no residual argument at all — Ritz values are upper bounds
          by Courant–Fischer — and agrees with the inertia-certified true
          constant to 5.2e-7 on every window, retiring the size-m
          factorisations from the ceiling step; the obstruction acquires a
          geometry — seven of eight bottom directions agree with the parity
          Laplacian to a degree, ONE sits at 83–90°, and that single
          misalignment IS the collapse price, recovered in full per window
          by one Cholesky (per-window end-to-end 4.45e-2–3.13e-1, inside
          the target band; the m-free-in-shape number unchanged at
          1.01e-2–3.92e-2) — both remaining terms are now uniformity terms
          with certified per-window numbers, the most accessible
          arithmetic-free (only the tridiagonal parity Laplacian and one
          8-dim subspace, worth 91–100% of the price); the certificate core
          is promoted as v551. T155 (BOTTOM.FLOOR,{" "}
          <span className="font-mono text-slate-300">
            bottom_floor_probe.py
          </span>
          , FLOORS-RESIST, 27/27) then attacked exactly those two floors:
          the complement floor becomes an EXACT fixed-size certificate — a
          12×12 problem in the exact Kac–Murdock–Szegő numbers plus one
          12×8 overlap matrix reproduces the size-m eigenproblem to
          0.999999783–0.999999958 on 16/16 windows, the defect is localised
          at the bottom mode itself (uncovered 0.539–0.611; the 83–90°
          direction lives on modes 9–12, which is exactly why K = 12), and
          the collapse price is recovered 78.8–100% at fixed size (end to
          end 3.28e-2–2.83e-1 CERTIFIED AT FIXED SIZE; the declared 4e-2
          bar missed by 1.22 and reported) — while on the block side
          λ_min(B_HH) = 0.2430–0.4249 is certified positive but the full
          symbol infimum is negative on every window (−714.2…−7.6): NO
          symbol argument can ever deliver that floor, the mechanism must
          be Fejér cancellation in the finite section (the arch half alone
          IS its symbol infimum to 5e-4 — a theorem candidate). T156
          (TWELVE.EIGHT,{" "}
          <span className="font-mono text-slate-300">
            twelve_by_eight_probe.py
          </span>
          , TERMS-RESIST, 37/37) then attacked exactly those two objects,
          and both remaining terms are now single scalars: on
          span{"{"}t₁, A⁻¹L_P t₁{"}"} the coupling is the identity t₁ᵀAy₁
          = μᴾ₁ (no A in it), so the t₁-loss is an EXACT closed function
          F(P, r) of two dimensionless numbers (verified to 2e-16; P =
          5.6e2–1.0e6, r = 2.7158–3.4089 flat), the two-line ceiling
          r ≤ 1/(Ls) ≤ 1/p₁ is tight to 1.03–1.16, and R2″ is ONE m-free
          lower bound on the angle p₁ = cos²∠(t₁, e₁(A)) = 0.2010–0.3282
          (55.1–63.4°) — with the separate measured debt that the 2×2
          model&apos;s domination of the 8-dimensional defect fails on the
          no-go family; on the block side the mechanism is identified
          against expectation — the expected arch inequality ≥ 1 is
          refuted (0.8226–1.3973, below 1 on 3/12; the weaker ≥ t form
          survives at factor 3.29–5.59), the Fejér damping is worth only
          5–31 where the split needs 3–1981 growing, the additive split is
          divergent (atom norm h^2.31), and the positivity is an ALIGNMENT
          fact: arch 1.47–91.71 and atom −91.27…−1.12 cancel to
          0.2661–0.4436 on the minimiser, at 52.7–90.0° from the
          atom-extremal vector; the no-go breaks on three axes, including
          a collapse of p₁ itself (7e-15–7e-10). T157 (ANGLE.FLOOR,{" "}
          <span className="font-mono text-slate-300">
            angle_floor_probe.py
          </span>
          , ANGLES-RESIST, 32/32) then attacked exactly those two angles,
          and neither falls, but both change shape: the sine-block
          confinement theorem ‖γ_H‖² ≤ (S/t)/ρ₁₇ = 0.0165–0.0293 proves,
          from the certified floor and ladder alone, that the bottom
          eigenvector lives 97.1–98.4% inside the first sixteen parity
          sines (a measured statement replaced by one line), the resolvent
          route pins p₁ ≥ 0.1968–0.3228 — 97.9% of the measurement — with
          ĝ₁² the one measured fixed-size scalar left, the whole first
          term becomes a bound on ONE diagonal entry of the 16×16 Schur
          complement the chain already forms (1/s = (S_L)₁₁ =
          2.3359–6.2049 flat; Cauchy–Schwarz misses it by 5.4e2–5.2e5 —
          the cancellation is nearly complete), the arch half is uniformly
          certified via an executed adaptive Lipschitz ceiling (12/12,
          cost h^0.85) with the two extremals at opposite ends of the band
          (a proof must use the θ-growth, not the infimum at π), and the
          alignment term stays a per-window domination with a 7.3e-4
          margin and a shrinking trend; the no-go family breaks on five
          axes, and the four instrument candidates of T155/T157 are
          promoted as v552. T158 (SCHUR.ENTRY,{" "}
          <span className="font-mono text-slate-300">
            schur_entry_probe.py
          </span>
          , ENTRY-RESISTS, 36/36) then found exactly that bound, and it is
          a THEOREM: the Thomson dual form turns the entry into a Dirichlet
          maximum — every trial vector bounds it from the right side, which
          is why Cauchy–Schwarz missed by 3.13e3–5.18e5 (a maximum evaluated
          at one direction — the wrong variational structure) — and the
          Cholesky ladder of strictly positive terms pins the entry to
          1/g₁₆ = 2.9670–7.9664, tight to 1.1323–1.2738, flat; the entry is
          exactly two-dimensional (span{"{"}t₁, A⁻¹L_P t₁{"}"} attains s to
          1e-9, because L_P t₁ = μᴾ₁ t₁), the T157 growth pointer is
          REFUTED (the atom mass grows faster than the arch floor in every
          band; band-local domination fails 21/21; what carries the
          inequality is the off-band coupling a dyadic argument discards),
          the T156 debt moves MEASURED → CERT-UNIF, and the measured-step
          count falls from three to one; what remains is one m-free upper
          bound on ONE 16×16 quadratic form whose halves cancel to depth
          h², and the sign structure of the off-band arch entries. T159
          (EXACT.FORM,{" "}
          <span className="font-mono text-slate-300">
            exact_form_probe.py
          </span>
          , FORM-RESISTS, 41/41) then executed the cancellation
          algebraically — seven machine-checked identities, among them
          the gauge identity (the form is exactly blind to the lag mass,
          which is where the h²-sized halves live) and the closed
          Dirichlet-kernel weights; the h² does not telescope away, and
          the archimedean bulk block turns out to be a raw symmetric
          Z-matrix with a closed sign-based floor above target (theorem
          count 6 → 13, cores promoted as v553). T160 (PAIRING,{" "}
          <span className="font-mono text-slate-300">
            pairing_probe.py
          </span>
          , PAIRING-RESISTS, 46/46) then attacked the pairing through
          its correlation structure, and the hardness is located: the
          machine-checked sampling identity shows the atom half IS a
          Λ-weighted prime sum at 32 explicit frequencies, needed to
          depth h⁻² and measured to cancel only to 0.00–0.37 of the
          trivial bound, while closed moment laws evaluate the smooth
          arch half down to the double-precision floor — the h²
          cancellation is the intrinsic arithmetic hardness of the
          problem. T161 (CLASSICAL.CLOSURE,{" "}
          <span className="font-mono text-slate-300">
            classical_closure_probe.py
          </span>
          , CLOSURE-RESISTS, 35/35) then ran the circularity triage,
          and the chain is NOT circular: the 32 frequencies are exactly
          the Fourier harmonics of the log-window (a theorem), the
          measured cancellation is fully the PNT main term, and the
          required depth has δ = 1.148–1.881 against RH strength 1/2 —
          below the boundary term of every partial summation, so the
          chain needs MORE than RH-strength input would supply: the h²
          sits in the SPLIT, not in the primes. The analytic arch half
          closes m-free (rate (3+√5)/2, closed head split, degree
          schedule O(log h)); R-B is refuted and replaced by a
          certified fraction bound ≤ 1/4; the closed cores are promoted
          as v554. T162 (THIRD.SPLIT,{" "}
          <span className="font-mono text-slate-300">
            third_split_probe.py
          </span>
          , DELTA-REDUCED, 30/30) then delivered the third split and
          measured its end: the archimedean Mellin ladder lowers the
          demand 1.88 → 1.38 → 0.93 in closed cell moments but saturates
          at K* = 2 (an asymptotic series — past the optimum the residual
          rises ×15–25), one Abel step makes the demand prime-free
          (δ_bnd = 1/2 + log(2κ‖Δw‖₁/|Q|)/log X exactly, optimal level
          one for the closed reason 32π/α &gt; 1), and the Fejér split
          pushes the proof demand BELOW the RH threshold 1/2 on all 18
          windows (δ_bnd = 0.133–0.417) — at a price in the 1/s ceiling
          growing h^2.86: relocated, not closed. R-A′ closes (the
          log-moment as a Lerch/Frullani integral, machine precision on
          three routes); R-B′ is refuted (the 16×16 Gram form is
          indefinite on every window; the a-weighted quarter bar
          survives). T163 (PARETO.FRONT,{" "}
          <span className="font-mono text-slate-300">
            pareto_front_probe.py
          </span>
          , FRONT-RESISTS, 32/32) then surveyed exactly that front, and
          it resists AS A THEOREM: the exchange law makes price and
          demand two coordinates of one exact identity (δ_bnd &lt; 1/2 ⟺
          P &gt; 2κg₁₆TV), the knob curve IS the front (monotone, both
          endpoints the chain&apos;s own ladder rungs), no flat price
          crosses (0/27), and the four-line TV-floor theorem — the entry
          normalisation x₁ = 1 meeting the smallest parity eigenvalue —
          forces every admissible trial vector to pay ~ h²: the hardness
          was never in the primes, it is the spectral gap of the parity
          Laplacian; R-C‴ closes negatively, the cores of T162+T163 are
          promoted as v555, and the successor R-E has two prime-free
          arms. T164 (SECTOR.CHANGE,{" "}
          <span className="font-mono text-slate-300">
            sector_change_probe.py
          </span>
          , TOLERANCE-CARRIES / SECTOR-RESISTS, 28/28) then decided both
          arms: the chain spends its ceiling at power exactly one (no
          free tolerance), yet the O(1) gate is discharged window by
          window by a Cholesky identity — the last analytic term
          collapses onto one quantifier — and no sector can help, as a
          theorem: the entry normalisation is a gauge, floor × transfer
          reproduces the h² identically. T165 (ALIGNMENT.ETA,{" "}
          <span className="font-mono text-slate-300">
            alignment_eta_probe.py
          </span>
          , ETA-RESISTS, 30/30) then closed the alignment successor by
          certificate: the P_pr identity makes demand and price one
          equation, every crossing vector pays the KMS h², and exactly
          ONE genuine open object remains — the cascade lower bound
          inf_m g₁₆(m) &gt; 0; the theorem cores of T164+T165 are
          promoted as v556. T166 (SCHUR.CASCADE,{" "}
          <span className="font-mono text-slate-300">
            schur_cascade_probe.py
          </span>
          , CASCADE-RESISTS, 30/30) then dissected exactly that lower
          bound: the whole gain is a near-collinearity
          (g_K/g₁ = 1/(1 − R_K²)), rung 2 alone carries a median 59%,
          and the gain is invariant under the spectral normalisation —
          the h³ belongs to the arithmetic Gram block alone; the
          cancellation lives in a 2×2 Gram determinant (no entry
          cancels, neither half is collinear alone — the same
          arch-against-atom mechanism as T159/T160, one level up), an
          anti-fitting scramble destroys it by a factor 4569, and the
          best closed route reaches h^+1.32 against the target h^+3.11
          — an exponent gap: the one missing inequality IS the
          cancellation, in three equivalent dresses. T167 (NULL.VECTOR,{" "}
          <span className="font-mono text-slate-300">
            null_vector_probe.py
          </span>
          , VECTOR-RESISTS, 39/39) then closed the most constructive
          dress as a construction: at K = 2 the closed vector
          (1, −Q₂₁/Q₂₂) is EXACT — a theorem — so the third dress
          collapses onto the second, and the machine-checked identity
          eps_ent(K) = ρ·(1/g_K)/S_K unifies all three dresses into ONE
          inequality; the threshold is mildest exactly where the vector
          is free (K = 2), the Kato series converges to the wrong
          object, the scrambled control loses positivity of the 2×2
          block itself, and the rest is one scalar — an m-free upper
          bound on 1 − r₁₂²; the theorem cores of T166+T167 are promoted
          as v557. T168 (LAGRANGE.MINORS,{" "}
          <span className="font-mono text-slate-300">
            lagrange_minors_probe.py
          </span>
          , MINORS-RESIST, 39/39) then ran the Lagrange identity on that
          scalar: the sum of squares is real (the arithmetic kernel is
          positive definite on all 63 windows, hard-fenced per window),
          the Wronskian minors are closed and maximal in norm — all
          smallness sits in the PSD kernel at one closed vector — and
          the one open factor is a single ratio t* = Q₁₂/Q₁₁, by
          T168-TH7 the target itself rewritten: the hardness is
          self-similar under exact reformulation. T169 (TSTAR.RATIO,{" "}
          <span className="font-mono text-slate-300">
            tstar_ratio_probe.py
          </span>
          , RATIO-RESISTS, 41/41) then PROVED the self-similarity: the
          only candidate meeting the threshold, √(â₂₂/â₁₁), collapses by
          a new identity back onto det Â — the loop closes as an
          identity; the gains are the first CERT-UNIF in weeks
          (Gershgorin on ν₁, unconditional), an R4-free chain that never
          touches the Weil-shaped positivity, and R1 in standard shape
          for the first time — a bilinear von Mangoldt sum against
          closed Dirichlet weights. T170 (BILINEAR.SIEVE,{" "}
          <span className="font-mono text-slate-300">
            bilinear_sieve_probe.py
          </span>
          , SIEVE-RESISTS, 40/40) then ran that toolbox, and the
          classification completed as a theorem: the form is exact (a
          five-order cancellation between three ~200-sized pieces on
          the reference window) and collapses back onto the linear
          hardness — rank K ≤ 3 for every window (the kernel is the
          polarisation of the determinant, so the form IS a rank-3
          polynomial in three linear Λ-sums) and Vaughan Type II blocks
          of effective rank O(1); no unconditional route exceeds
          δ = +0.996 against the target 3.0, and R1 is finally
          classified as a NEAR-DEGENERACY — two explicit Λ-sum vectors
          becoming collinear at h^−3 — beyond any size-bounding tool;
          the theorem cores of T168+T169+T170 are promoted as v558.
          T171 (FINAL.MAP,{" "}
          <span className="font-mono text-slate-300">
            final_map_probe.py
          </span>
          , MAP-COMPLETE, 43/43) then assembled the capstone, and the
          map is complete: all sixteen links of the reduction chain
          reproduce in one connected run (13 theorems, 3 certificates),
          all eight classified no-go routes fail as classified, and the
          precision ledger closes 1.2e5× beyond the RH yardstick — zero
          new uniform statements, the capstone promoted as v559: phase
          2 is a certified map with one open object (R1). T172
          (FRAME.BEYOND,{" "}
          <span className="font-mono text-slate-300">
            frame_beyond_probe.py
          </span>
          , PARTIALLY-PORTABLE, 44/44) then tested how far the map
          carries beyond frame A: 13 of the 16 links transfer UNCHANGED
          to a gap-blind frame, to ν = 3 and ν = 8, to non-prime-power
          anchors and to both congruence classes mod 4 (0 broken — the
          three that shift are exactly the number-carrying ones), the
          indefiniteness is localised at the sieve horizon, R1&apos;s
          near-degeneracy persists on all 54 windows while its rate is
          frame-bound (h^−1.83 to h^−2.87), and the scramble at fixed
          Λ-value multiset removes the decay — the collapse belongs to
          the actual prime-power placement. T173 (FRAME.RATE,{" "}
          <span className="font-mono text-slate-300">
            frame_rate_probe.py
          </span>
          , DEFICIT-VARIES, 40/40) then showed the demanded rate is
          itself a frame datum (q = 1 − s exact, q &lt; 1 on all eleven
          frame data — the h^−3 target was the q = 1 idealisation) and
          made the deficit between demand and delivery the number of
          the phase: +0.155 ± 0.102, invariant under the anchor-to-grid
          rule and the lever split, not yet constant in ν — no frame
          closes it, frame shopping is over by numbers, and R1 stands
          frame-free. T174 (CANCEL.IDENTITY,{" "}
          <span className="font-mono text-slate-300">
            cancellation_identity_probe.py
          </span>
          , PARTIAL-CANCEL) then exhausted the gauge route by theorems —
          everything multiplicative cancels exactly, the additive
          arch/comb mixture admits no factorisation — and measured the
          deficit frame-free at +0.1111 ± 0.0222 (5.0σ), its driver the
          comb density per lag cell; the theorem cores of T172–T174 are
          load-bearing as v560. T175 (PHASE.PLACEMENT,{" "}
          <span className="font-mono text-slate-300">
            phase_placement_probe.py
          </span>
          , PHASES-RESIST) then measured the placement phases directly:
          they are real and causal for R (an intervention with exact
          rebuild control, dlog R/dδ = 879, linear over four decades),
          but no phase formula is certifiable — the response is
          non-smooth, the Schur floor a near-degeneracy — the
          heterogeneity dissolved substantially into the error bar (the
          5σ deficit untouched), and the densest reachable density bin
          is consistent with zero, undecidable under this sieve. T176
          (DENSE.LIMIT,{" "}
          <span className="font-mono text-slate-300">
            dense_limit_probe.py
          </span>
          , SITS-AT-ZERO) then ran the larger sieve — the last
          decidable measurement of phase 2: the density ceiling rose
          from 361 to 6120, both new bins are consistent with zero
          (pooled +0.0496 ± 0.0372), the plateau window narrowed by a
          factor 3.6, and the measurement programme closed as planned
          (the exact cores are load-bearing as v562; the work continues
          in the classification papers and the backflow lines).
        </li>
      </ul>
    </aside>
  );
}

/** Depth reach of the certified chain across the sprint. Numbers from T110–T116. */
const DEPTH_REACH = [
  { part: "T110", label: "certified handover steps", value: "15" },
  { part: "T111", label: "handovers, retention 1.000000", value: "117/117" },
  { part: "T114", label: "deepest certified zone", value: "n = 1331" },
  { part: "T115", label: "deepest certified zone", value: "n = 155,921" },
  { part: "T116", label: "boundary march, flat cost", value: "1.35M cells" },
] as const;

function DepthReachCallout() {
  return (
    <aside className="rounded-2xl border border-emerald-400/35 bg-emerald-500/10 p-5 sm:p-6">
      <p className="font-mono text-[10px] uppercase tracking-[0.18em] text-emerald-300/90">
        Depth reach · what the chain certified, part by part
      </p>
      <dl className="mt-4 space-y-2.5">
        {DEPTH_REACH.map((d) => (
          <div
            key={`${d.part}-${d.value}`}
            className="flex flex-wrap items-baseline gap-x-2 border-b border-emerald-400/15 pb-2 last:border-0 last:pb-0"
          >
            <dt className="font-mono text-[11px] text-emerald-200">{d.part}</dt>
            <dd className="font-mono text-sm text-emerald-50">{d.value}</dd>
            <dd className="min-w-0 flex-1 text-right text-[11px] text-slate-400">
              {d.label}
            </dd>
          </div>
        ))}
      </dl>
      <p className="mt-4 text-xs leading-relaxed text-slate-400">
        Chains stop at the compute budget, never at a step. T116&apos;s march is
        a cost-geometry demonstration, not a Weil certificate — declared as such
        in the probe. Sandbox; not RH evidence.
      </p>
    </aside>
  );
}

function AmplitudeRouteCallout() {
  return (
    <aside className="rounded-2xl border border-sky-400/35 bg-sky-500/10 p-5 sm:p-6">
      <p className="font-mono text-[10px] uppercase tracking-[0.18em] text-sky-300/90">
        Transport compressed · v540
      </p>
      <p className="mt-3 font-serif text-lg leading-relaxed text-sky-50 sm:text-xl">
        The transport problem, compressed to a single measurable object: the
        FE-covariant gap functional λ* on the atoms n&nbsp;≡&nbsp;6&nbsp;mod&nbsp;8.
      </p>
      <p className="mt-4 text-xs leading-relaxed text-slate-400">
        Square closed → linear plus-carrier stands → cone library saturates at
        5/24. Euler-region positivity ≠ central-line claim. This is not RH
        evidence.
      </p>
    </aside>
  );
}

function DoorsFurnishedCallout() {
  return (
    <aside className="rounded-2xl border border-amber-400/35 bg-amber-500/10 p-5 sm:p-6">
      <p className="font-mono text-[10px] uppercase tracking-[0.18em] text-amber-300/90">
        The wall, named
      </p>
      <p className="mt-3 font-serif text-lg leading-relaxed text-amber-50 sm:text-xl">
        The wall is no longer a metaphor: it is one named inequality —
        provably equivalent to RH. Everything else on the route is either
        proved or classically provable-shaped.
      </p>
      <p className="mt-4 text-xs leading-relaxed text-slate-400">
        After T81 the recipe chain saturated. Three new perspectives (next)
        retype I5 without claiming RH progress.
      </p>
    </aside>
  );
}

function ThreePerspectivesCallout() {
  return (
    <aside className="rounded-2xl border border-sky-400/35 bg-sky-500/10 p-5 sm:p-6">
      <p className="font-mono text-[10px] uppercase tracking-[0.18em] text-sky-300/90">
        New I5 form — fence
      </p>
      <p className="mt-3 font-mono text-sm leading-relaxed text-sky-100">
        {"Q_cert(h) + Δ₂(h) + A_fam(h) − A_shift(h) ≥ 0"}
      </p>
      <p className="mt-3 font-serif text-base leading-relaxed text-sky-50">
        Self-consistency of one heat family. Type change ≠ proof; I5 remains ⟺
        RH. Coherent class closed via λ-channel.
      </p>
      <p className="mt-4 text-xs leading-relaxed text-slate-400">
        T86{" "}
        <span className="font-mono text-amber-200">
          [sandbox · LEMMA-FULLY-CLOSED]
        </span>
        : all atom classes covered. Geography next. This is not RH evidence.
      </p>
    </aside>
  );
}

function BlindPrimeCallout() {
  return (
    <aside className="mt-8 rounded-2xl border border-emerald-400/40 bg-emerald-500/10 p-5 sm:p-6">
      <p className="font-mono text-[10px] uppercase tracking-[0.18em] text-emerald-300/90">
        Blind demo · T94 · sandbox
      </p>
      <p className="mt-3 font-serif text-xl leading-relaxed text-emerald-50 sm:text-2xl">
        753/753 primes predicted — zero errors, no division.
      </p>
      <p className="mt-3 text-sm leading-relaxed text-slate-300">
        All primes in the never-seen window [1,000,001–1,010,000], from pure
        lattice counting (Jacobi 1834: odd n prime ⟺ r₄(n)&nbsp;=&nbsp;8(n+1)),
        with an AST-enforced prediction path containing no divisibility tests.
        Predictions committed by MD5 before any truth was computed. Honesty: a
        sieve is ~820× faster — structure completeness, not algorithmic
        progress. Spectral prediction (zeros&nbsp;→&nbsp;primes) remains bound
        to I5/RH.
      </p>
    </aside>
  );
}

function I5GeographyCallout() {
  return (
    <aside className="rounded-2xl border border-violet-400/35 bg-violet-500/10 p-5 sm:p-6">
      <p className="font-mono text-[10px] uppercase tracking-[0.18em] text-violet-300/90">
        Geography — not an attack
      </p>
      <p className="mt-3 font-serif text-lg leading-relaxed text-violet-50 sm:text-xl">
        Identities proved; zones 2–4 closed; zone-5 tip = equality; asymptotics
        = one named arithmetic bound (A).
      </p>
      <p className="mt-4 text-xs leading-relaxed text-slate-400">
        I5 remains ⟺ RH; the geography locates where any attack must work, it
        does not perform one. 2428/2428 sandbox checks at T101 (5360/5360
        today). This is not RH evidence.
      </p>
    </aside>
  );
}

function MeaningLevel({
  level,
  badge,
  title,
  body,
}: {
  level: string;
  badge: "sandbox" | "machine-verified";
  title: string;
  body: string;
}) {
  return (
    <li className="rounded-2xl border border-slate-700/45 bg-slate-950/40 p-4">
      <div className="flex flex-wrap items-center gap-2">
        <span className="font-mono text-[10px] uppercase tracking-widest text-slate-500">
          {level}
        </span>
        <StatusBadge badge={badge} />
      </div>
      <h3 className="mt-1 font-serif text-lg text-slate-100">{title}</h3>
      <p className="mt-1 text-sm leading-relaxed text-slate-400">{body}</p>
    </li>
  );
}
