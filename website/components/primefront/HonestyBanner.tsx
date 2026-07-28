import Link from "next/link";

export function HonestyBanner() {
  return (
    <aside
      role="note"
      aria-label="Research honesty notice"
      className="relative overflow-hidden rounded-2xl border border-amber-400/40 bg-gradient-to-br from-amber-500/15 via-slate-950/80 to-slate-950/90 p-5 sm:p-6"
    >
      <div
        aria-hidden
        className="pointer-events-none absolute -right-8 -top-8 h-32 w-32 rounded-full bg-amber-400/10 blur-2xl"
      />
      <p className="font-mono text-[10px] font-semibold uppercase tracking-[0.2em] text-amber-300">
        Honesty first
      </p>
      <p className="mt-2 max-w-3xl text-sm leading-relaxed text-amber-50/95 sm:text-base">
        Research diary. Everything here is exploratory sandbox work unless
        explicitly marked as machine-verified (currently:{" "}
        <Link
          href="/verification"
          className="font-mono text-emerald-300 underline decoration-emerald-400/40 underline-offset-2 hover:text-emerald-200"
        >
          v535–v541
        </Link>
        ). Load-bearing v539 isolates <em>two</em> obstructions; load-bearing
        v540 consolidates the amplitude/linear route with the open boundary{" "}
        <em>inside</em> the claim; load-bearing v541 promotes the
        matching-lemma/transport-ledger package. Sandbox T86–T93 complete the
        I5 geography (with honest self-corrections through T98); T94 is a blind
        lattice prime demo; T99–T101 close the induction arc; and the sprint
        T102–T125 compresses the remaining arithmetic bound step by step — from
        one matrix inequality down to one sign plus one declared accounting
        convention, assembled end to end in the T125 finale, with certified steps
        as deep as zone 155,921. The series is complete at 125 parts; phase 2
        (“the full proof”, T126+) is now open with 3170/3170 sandbox checks and
        exactly two genuinely new inequalities remaining. What remains
        TFPT-specific is exactly one object: I5 in one-family form — by the
        closed ledger equivalent to Weil positivity ⟺ RH, an equivalence typing
        only, not “almost RH.” No claim of progress toward the Riemann
        Hypothesis is made.
      </p>
    </aside>
  );
}
