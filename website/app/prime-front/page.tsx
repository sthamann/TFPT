import type { Metadata } from "next";
import Link from "next/link";
import { ArrowLeft } from "lucide-react";
import { SITE_URL } from "@/lib/utils";
import { PRIME_FRONT_SECTIONS } from "@/lib/primeFront";
import { HonestyBanner } from "@/components/primefront/HonestyBanner";
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
import { UpdateFeed } from "@/components/primefront/UpdateFeed";
import { StatusBadge } from "@/components/primefront/StatusBadge";

export const metadata: Metadata = {
  title: "The Prime Front — Research Diary",
  description:
    "A plain-language research diary of TFPT's prime / zeta investigation: signed E8 census, Hecke from geometry (v535–v537), the relative-trace identity (v538), Weil structure with two isolated obstructions (v539), the amplitude route with positive linear carrier whose open boundary is the FE-covariant gap functional λ* on n ≡ 6 mod 8 (v540), the matching-lemma/transport-ledger package with two named limits (v541), the I5 geography through T94, and a blind lattice prime demo (753/753). Not RH evidence.",
  keywords: [
    "TFPT prime front",
    "E8 census",
    "Hecke from geometry",
    "v535",
    "v538",
    "v539",
    "v540",
    "v541",
    "Weil structure",
    "relative trace formula",
    "zeta research diary",
    "Stefan Hamann",
  ],
  alternates: { canonical: `${SITE_URL}/prime-front` },
  openGraph: {
    type: "article",
    title: "The Prime Front — TFPT research diary",
    description:
      "Primes, the E8 census, and load-bearing results v535–v541. Program status: I5 geography through T93 (a_neg = 0.7486); blind lattice demo 753/753; residual object I5. No claim of progress toward the Riemann Hypothesis.",
    url: `${SITE_URL}/prime-front`,
    siteName: "TFPT — Topological Fixed-Point Theory",
    locale: "en_US",
  },
  twitter: {
    card: "summary_large_image",
    title: "The Prime Front — TFPT",
    description:
      "Research diary of the prime / zeta line. Machine-verified: v535–v541. Induction closeout T99–T101: zones 2–4 closed; asymptotics = one named arithmetic bound (2428 sandbox checks). Residual: I5 — not almost-RH.",
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

          <HonestyBanner />

          <p className="mt-8 font-mono text-[11px] uppercase tracking-[0.2em] text-sky-300/90">
            Research diary · Teile 11–121 · 3037 sandbox checks
          </p>
          <h1 className="mt-3 font-serif text-4xl font-semibold leading-tight text-slate-50 sm:text-5xl md:text-6xl">
            The Prime Front
          </h1>
          <p className="mt-4 max-w-2xl text-base leading-relaxed text-slate-300 sm:text-lg">
            Primes are the simplest objects with the most complex behaviour.
            This page tells — in ordinary language — how a discrete physics
            compiler&apos;s bookkeeping kept speaking number theory, what
            survived preregistered kills, and what is still honestly open.
          </p>

          <BlindPrimeCallout />

          <nav
            aria-label="On this page"
            className="mt-8 flex flex-wrap gap-2"
          >
            {PRIME_FRONT_SECTIONS.map((s) => (
              <a
                key={s.id}
                href={`#${s.id}`}
                className="rounded-full border border-slate-700/50 bg-slate-900/40 px-3 py-1 font-mono text-[11px] text-slate-400 transition hover:border-slate-500 hover:text-slate-200"
              >
                {s.label}
              </a>
            ))}
          </nav>
        </div>
      </section>

      {/* 1 — Hook */}
      <DiarySection
        id="hook"
        eyebrow="01 · Hook"
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

      {/* 2 — Compiler */}
      <DiarySection
        id="compiler"
        eyebrow="02 · The compiler in one picture"
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

      {/* 3 — Signed census */}
      <DiarySection
        id="census"
        eyebrow="03 · First discovery · Teil 11"
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

      {/* 4 — Bridges */}
      <DiarySection
        id="bridges"
        eyebrow="04 · Surprise bridges · Teil 12"
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

      {/* 5 — Kill chain */}
      <DiarySection
        id="kill-chain"
        eyebrow="05 · Honesty as a method · Teile 14, 19–25"
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

      {/* 6 — Predict */}
      <DiarySection
        id="predict"
        eyebrow="06 · Can it predict primes? · Teil 21"
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

      {/* 7 — Hecke */}
      <DiarySection
        id="hecke"
        eyebrow="07 · The mechanism · Teile 27–32 · v535"
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

      {/* 8 — Eichler */}
      <DiarySection
        id="eichler"
        eyebrow="08 · The Eichler layer · Teile 33, 36 · v536"
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

      {/* 9 — Weight drop */}
      <DiarySection
        id="weight-drop"
        eyebrow="09 · Two-channel weight drop · Teile 35, 39"
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

      {/* 10 — Stage 4 */}
      <DiarySection
        id="stage-4"
        eyebrow="10 · Stage-4 map · Teil 40"
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

      {/* 11 — July 25 arc */}
      <DiarySection
        id="july-25-arc"
        eyebrow="11 · The July 25 arc · Teile 51–64 · v538 / v539"
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

      {/* 12 — Program status (consolidated) */}
      <DiarySection
        id="program-status"
        eyebrow="12 · Where the program stands"
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
          After 91 probes and seven promoted modules (v535–v541), the matching
          lemma is closed on ALL atom classes (window-certificate format,
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
          Milestone: 2428/2428 sandbox checks. Fence: fits/extrapolations marked;
          geography locates; it does not attack; I5 remains ⟺ RH. This is not
          RH evidence.
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

      {/* 13 — Amplitude / transport route */}
      <DiarySection
        id="amplitude-route"
        eyebrow="13 · The amplitude route · Teile 67–72 · v540"
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

      {/* 14 — Doors furnished */}
      <DiarySection
        id="doors-furnished"
        eyebrow="14 · The doors get furnished · Teile 73–81"
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

      {/* 15 — Three perspectives */}
      <DiarySection
        id="three-perspectives"
        eyebrow="15 · Three new perspectives · Teile 82–84"
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

      {/* 16 — I5 geography */}
      <DiarySection
        id="i5-geography"
        eyebrow="16 · I5 geography · Teile 87–101"
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
          must work, it does not perform one. Milestone: 2428/2428 sandbox
          checks. All laws marked as fits/extrapolations. This is not RH
          evidence.
        </p>
      </DiarySection>

      {/* 17 — The mechanism (T102) */}
      <DiarySection
        id="mechanism"
        eyebrow="17 · The mechanism · Teil 102"
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

      {/* 18 — The instrument rebuilt (T101 → T103) */}
      <DiarySection
        id="instrument-race"
        eyebrow="18 · The instrument rebuilt · Teile 101 → 103"
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

      {/* 19 — Two doors converge (T102 + T103 → T104) */}
      <DiarySection
        id="two-doors"
        eyebrow="19 · Convergence · Teile 102–104"
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
          T105 (BARE.AVOIDANCE.CORE) is running. Sandbox; not RH evidence.
        </p>
      </DiarySection>

      {/* 20 — Meaning */}
      <DiarySection
        id="meaning"
        eyebrow="20 · What would it mean"
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
            title="I5 geography complete — one object remains"
            body="After 101 probes (2428/2428 sandbox checks) and v535–v541: matching lemma closed; I5 geography complete; induction skeleton closed through T101 — identities proved, zones 2–4 closed, zone-5 tip = equality problem, asymptotics = one named arithmetic bound (A). What remains TFPT-specific is exactly ONE object: I5 in one-family form ⟺ Weil positivity ⟺ RH. Geography locates; it does not attack. Not RH evidence."
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

      {/* 21 — Live updates */}
      <section
        id="updates"
        aria-labelledby="updates-heading"
        className="scroll-mt-24 border-t border-slate-800/60 py-12 sm:py-16"
      >
        <div className="mx-auto max-w-5xl px-4 sm:px-6 lg:px-8">
          <div className="flex flex-wrap items-center gap-2">
            <span className="font-mono text-[10px] uppercase tracking-[0.18em] text-slate-500">
              21 · Live updates
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
            <UpdateFeed />
          </div>
        </div>
      </section>

      <footer className="border-t border-slate-800/60 py-10">
        <div className="mx-auto max-w-5xl px-4 text-sm text-slate-500 sm:px-6 lg:px-8">
          Numbers and verdicts are taken from{" "}
          <code className="font-mono text-slate-400">experiments/next.txt</code>{" "}
          (2026-07-23…26 diary) and the promoted modules{" "}
          <code className="font-mono text-slate-400">
            verification/v535_*.py
          </code>{" "}
          …{" "}
          <code className="font-mono text-slate-400">
            v541_matching_lemma_ledger.py
          </code>
          . Exploratory probes live under{" "}
          <code className="font-mono text-slate-400">
            experiments/tfpt-discovery/
          </code>
          .
        </div>
      </footer>
    </>
  );
}

function ProgramStatusCallout() {
  return (
    <aside className="rounded-2xl border border-violet-400/35 bg-violet-500/10 p-5 sm:p-6">
      <p className="font-mono text-[10px] uppercase tracking-[0.18em] text-violet-300/90">
        Strongest current sentence
      </p>
      <p className="mt-3 font-serif text-lg leading-relaxed text-violet-50 sm:text-xl">
        After 105 probes and seven promoted modules, what remains TFPT-specific
        is exactly one object: I5 — with an induction skeleton whose identities
        are proved, zones 2–4 closed, and asymptotics localized to one named
        arithmetic bound (A), whose mechanism T102 identifies (Schur dressing
        against E₀⊕E₊ — one scalar per zone).
      </p>
      <p className="mt-4 text-xs leading-relaxed text-slate-400">
        Zone-5 tip = equality problem; crowding sits in the proof family, not
        the mathematics — T103&apos;s improved instrument closes 44/64
        (probe-level), and T104 (CHAIN-PARTIAL, two independent arms) closes
        the exact spectral-split chain 16/16 with finite data, moving the hard
        core to a bare_k lower bound plus one soft dressing scalar L. T105
        (ONE-OF-TWO) certifies the bare bound in closed form and upgrades the
        avoidance law to a theorem — the whole remaining hardness is one
        Friedrichs-angle statement. T106 (DENSITY-MAPPED) splits that
        statement by parity: the even channel closes 16/16 and all remaining
        hardness localizes in the odd channel — one Loewner statement on half
        the dimensions. T107 (SCALAR-TRACTABLE) reduces that statement to
        exactly one scalar ratio r = κ/ε ≤ 1 with two orders of magnitude of
        measured room (r = 0.005…0.18); the symbol route is structurally
        dead, and exactly one object remains to certify. T108
        (EPSILON-IDENTITY) makes ε&apos;s positivity an exact identity — the
        square of the last Cholesky pivot (Szegő), coinciding with the
        induction positivity itself — and reduces (R) to two scalars; what
        remains is literally one boundary value of an explicit vector. T109
        (BOUNDARY-CERTIFIED) certifies both scalars — ω unconditionally via a
        graded matrix cap, the boundary value via a residual certificate that
        carries the cancellation — so the whole chain for (R) now closes 16/16
        conditional on exactly one strict-margin input that is 10²–10⁶ weaker
        than the conclusion. T110 (MARGIN-PROPAGATES) closes the induction
        circle end-to-end on the measured zones — certified base case, 15
        certified handover steps, atom entry structurally free — with three
        sharp gaps remaining (no reserve, no scalar step law, no
        k-uniformity). T111 (CROSSING-CONFIRMED) drives the ladder deep and
        measures the wall instead of extrapolating it: the crossing sits at
        n* ≈ 462 (an upper bound; n ≈ 170 was a fit artefact) and splits into
        three separate walls — margin, twin-prime ladder death, vacuous
        requirement — while the handover mechanism itself never fails
        (117/117 at retention 1.000000); depth, not n, is the operating
        variable. T112 (SCALING-PARTIAL) then rebuilds the construction in
        the gap-coupled scaled frame: two of the three walls fall
        structurally, but the margin wall is frame-invariant at exponent
        −0.974 — the hardness is now one limit operator plus one convergence
        rate, with the prime-gap dependence disclosed. T113
        (SUBSTANCE-CONFIRMED) settles the currency question with a
        fundamental reinterpretation: the margin wall is real in every
        currency, but it measures the discretization, not the spectrum — the
        continuum window form has no gap, the positive floor is a
        cancellation (~1e-7 relative), and the T109 requirement chain
        divides by an artifact margin. T114 (WALL-DISSOLVES) then builds
        the repair and the wall dissolves: the exact Schur complement
        (Albert 1969) certifies every step margin-free — eleven beyond the
        old wall, up to n = 1331, and all seven torn zones — the wall was
        an O(1) numerator divided by an artifact floor; chains now stop
        only at compute cost, never at a step, and the remaining core is
        two objects (ε relative to κ via the exact Szegő identity, and
        transporting an O(0.1) Schur complement between grids). T115
        (TRANSPORT-BLOCKED) then splits that core: transport certifies
        only mild refinement (the Schur floor itself falls like ρ^−1.7 on
        nested ladders — no bound can fix that), but the two-scale
        compression breaks the cap — a certified margin-free step at
        n = 155,921 (117× deeper), chains of 10 steps, the stopper always
        cost and never a step — leaving the shortest list ever: three
        points, only one of which is an inequality (the classical
        Szegő–Levinson prediction-error bound for one symbol). T116
        (RICCATI-PARTIAL) then runs the boundary reformulation: the
        induction step IS a boundary process — the global pole rides
        exactly in a 12×12+12+1 state and the Riccati march runs 169,236
        steps to 1.35 million cells at flat cost (a cost-geometry demo,
        not a Weil certificate) — but the prime comb itself refuses
        compression, and the one remaining inequality acquires textbook
        shape (ε is a Galerkin error; Aubin–Nitsche hits the measured
        exponent). T117 (THEOREM-SHAPED) then makes that inequality
        theorem-shaped: ε is an identity plus a certified lower-bound
        chain that loses no power of D (θ&apos; = 1.74 against θ = 1.76);
        the jumps have closed forms — T116&apos;s factor-120 drops were a
        sweep artifact, prime-power entries actually raise ε — and what
        remains are three named analytic lemmas about one symbol, each
        with a classical address, two of them constants, not rates. T118
        (TWO-OF-THREE) then measures the three against their addresses
        and two of them stand: saturation is an identity here, and the
        Schur floor is rescued by the CBS shift onto the oscillation
        Gram, whose arithmetic aliasing symbol suppresses the prime comb
        quadratically (certified on 7/14 windows; on a 15,680-point FFT
        lever the floor rises logarithmically and crosses zero — the
        failing windows were under-resolved, not obstructed); the corner
        lemma corrects itself (log is the boundary of all powers) and
        the chain is repaired on ‖y‖²; exactly one genuinely new
        analytic statement remains — a D^1.75 lower bound for the
        oscillation mass. T119 (ARITHMETIC-DONE) then closes the
        arithmetic half as a theorem — an explicit D₀(α) with a universal
        constant B = −1.0474 — proves the energy route to (H2) empty
        (genuinely new content), and reduces everything to one discrete
        Harnack inequality via the exact identity κ_end = 1/(1+R). T120
        (HARNACK-EXPLAINED) then explains that inequality: R ≈ 1 is
        parity symmetry — the two increment families are one sequence
        shifted by a single fine cell, giving an unconditional
        |R−1| ≤ 0.047 certificate, while the per-cell inequality is
        provably false — but the honest defect count rises 3 → 4: the
        frame&apos;s window grows with the zone, so D₀ covers almost
        nothing of the real ladder, and the CBS constant decays with the
        width. T121 (WIDE-RESTRUCTURED) then runs the σ_z route against
        the real ladder: it survives as a SECTION statement — positive
        16/16 up to α = 6.28, half of them where the symbol infimum is
        negative — and dies as a symbol statement; the net α-balance
        shows only a poly-log deficit (α^−1.57, uniform in D) instead of
        collapse, exactly decomposed into two repairable steps (Rayleigh
        + CBS, nothing in γ), and one link of the chain was honestly
        refuted (the odd-sector Hankel term); both repairs are under
        test in T122 (net_balance_repair_probe.py, running). Milestone:
        3037/3037 checks. Not almost-RH. This is not RH evidence.
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
        does not perform one. 2428/2428 sandbox checks. This is not RH
        evidence.
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
