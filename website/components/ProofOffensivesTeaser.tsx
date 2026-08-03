import Link from "next/link";
import { ArrowRight } from "lucide-react";

/**
 * Compact home teaser for the newest Prime-Front content: the August 3
 * proof offensives (v682–v700) and their four visual schematics on
 * /prime-front#august-offensives. Server component — no client JS; the
 * subtle motion is pure CSS.
 */
export function ProofOffensivesTeaser() {
  return (
    <section
      id="august-offensives-teaser"
      className="relative scroll-mt-20 py-10 sm:py-12"
      aria-labelledby="august-offensives-teaser-heading"
    >
      <div className="mx-auto max-w-6xl px-4 sm:px-6 lg:px-8">
        <Link
          href="/prime-front#august-offensives"
          className="group block rounded-2xl border border-sky-400/25 bg-gradient-to-r from-sky-500/[0.07] via-slate-900/40 to-emerald-500/[0.06] p-5 transition-colors hover:border-sky-300/40 sm:p-6"
        >
          <div className="flex flex-wrap items-center gap-2">
            <span
              className="h-1.5 w-1.5 animate-pulse rounded-full bg-sky-300"
              aria-hidden
            />
            <span className="font-mono text-[10px] font-semibold uppercase tracking-widest text-sky-300/90">
              New · August 3 · the Prime Front
            </span>
            <span className="font-mono text-[10px] text-slate-500">
              v682–v700 · 19 modules in one day
            </span>
          </div>

          <h2
            id="august-offensives-teaser-heading"
            className="mt-3 font-serif text-xl font-semibold leading-snug text-slate-50 sm:text-2xl"
          >
            The proof offensives of August 3 — now with four pictures
          </h2>

          <p className="mt-2 max-w-3xl text-sm leading-relaxed text-slate-400">
            The Hecke-SOS blueprint with its one named gap (Z1), the T-B
            closure map (60 of 70 windows closed
            unconditionally-modulo-citations), and two sandbox explorations —
            the positivity corridor of the primes and the deck-sector split of
            the arch density. Promoted results and exploration clearly
            badged; no RH claim.
          </p>

          <span className="mt-3 inline-flex items-center gap-1.5 text-sm font-semibold text-sky-300 transition-colors group-hover:text-sky-200">
            See the four schematics
            <ArrowRight
              size={14}
              className="transition-transform group-hover:translate-x-0.5"
              aria-hidden
            />
          </span>
        </Link>
      </div>
    </section>
  );
}
