"use client";

import { motion } from "motion/react";
import { StatusBadge } from "./StatusBadge";

export function TwoDoorsConvergence() {
  return (
    <div className="rounded-2xl border border-slate-700/50 bg-slate-950/60 p-4 sm:p-5">
      <p className="mb-3 font-mono text-[10px] uppercase tracking-widest text-emerald-300/90">
        Two doors, one object · T102 + T103
      </p>

      <div className="grid grid-cols-2 gap-2">
        <motion.div
          initial={{ opacity: 0, y: 8 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true, amount: 0.4 }}
          transition={{ duration: 0.4 }}
          className="rounded-xl border border-violet-400/40 bg-violet-500/10 p-3"
        >
          <p className="font-mono text-[10px] uppercase tracking-wider text-violet-300/90">
            Door A · handover mechanism
          </p>
          <p className="mt-1.5 text-xs leading-relaxed text-violet-100">
            The onset is manufactured by the Schur dressing against E₀⊕E₊ —
            it takes 35.7%…97.3% of the bare eigenvalue.
          </p>
          <p className="mt-1.5 font-mono text-[10px] text-slate-500">
            T102 · arithmetic_bound_probe.py
          </p>
        </motion.div>

        <motion.div
          initial={{ opacity: 0, y: 8 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true, amount: 0.4 }}
          transition={{ duration: 0.4, delay: 0.15 }}
          className="rounded-xl border border-sky-400/40 bg-sky-500/10 p-3"
        >
          <p className="font-mono text-[10px] uppercase tracking-wider text-sky-300/90">
            Door C · instrument race
          </p>
          <p className="mt-1.5 text-xs leading-relaxed text-sky-100">
            The remaining loss is the wing slack S = 1 − ρ (falls
            0.2091 → 0.0392); the pencil is nearly saturated, the bulk is not
            low-rank.
          </p>
          <p className="mt-1.5 font-mono text-[10px] text-slate-500">
            T103 · instrument_probe.py
          </p>
        </motion.div>
      </div>

      {/* converging arrows */}
      <svg
        viewBox="0 0 320 52"
        className="w-full"
        aria-hidden
      >
        <motion.path
          d="M 80 4 C 80 28, 120 40, 152 46"
          fill="none"
          stroke="rgba(196,181,253,0.7)"
          strokeWidth={1.6}
          initial={{ pathLength: 0 }}
          whileInView={{ pathLength: 1 }}
          viewport={{ once: true, amount: 0.5 }}
          transition={{ duration: 0.7, delay: 0.5, ease: "easeInOut" }}
        />
        <motion.path
          d="M 240 4 C 240 28, 200 40, 168 46"
          fill="none"
          stroke="rgba(125,211,252,0.7)"
          strokeWidth={1.6}
          initial={{ pathLength: 0 }}
          whileInView={{ pathLength: 1 }}
          viewport={{ once: true, amount: 0.5 }}
          transition={{ duration: 0.7, delay: 0.65, ease: "easeInOut" }}
        />
        <motion.circle
          cx={160}
          cy={47}
          r={3.5}
          fill="#34d399"
          initial={{ opacity: 0, scale: 0 }}
          whileInView={{ opacity: 1, scale: 1 }}
          viewport={{ once: true, amount: 0.5 }}
          transition={{ duration: 0.3, delay: 1.3 }}
          style={{ transformBox: "fill-box", transformOrigin: "center" }}
        />
        <motion.text
          x={64}
          y={30}
          fontSize="8"
          className="fill-slate-500 font-mono"
          initial={{ opacity: 0 }}
          whileInView={{ opacity: 1 }}
          viewport={{ once: true, amount: 0.5 }}
          transition={{ duration: 0.3, delay: 0.9 }}
        >
          same object,
        </motion.text>
        <motion.text
          x={222}
          y={30}
          fontSize="8"
          className="fill-slate-500 font-mono"
          initial={{ opacity: 0 }}
          whileInView={{ opacity: 1 }}
          viewport={{ once: true, amount: 0.5 }}
          transition={{ duration: 0.3, delay: 0.9 }}
        >
          two sides
        </motion.text>
      </svg>

      <motion.div
        initial={{ opacity: 0, y: 8 }}
        whileInView={{ opacity: 1, y: 0 }}
        viewport={{ once: true, amount: 0.4 }}
        transition={{ duration: 0.45, delay: 1.4 }}
        className="rounded-xl border border-emerald-400/40 bg-emerald-500/10 p-3 text-center"
      >
        <p className="font-mono text-[10px] uppercase tracking-wider text-emerald-300/90">
          The wing near-null direction
        </p>
        <p className="mt-1.5 text-sm leading-relaxed text-emerald-50">
          One scalar per zone: a lower bound on σ_k just above atom entry.
        </p>
      </motion.div>

      <div className="mt-2 flex flex-wrap items-start gap-2 rounded-xl border border-dashed border-sky-400/35 bg-slate-900/50 p-3">
        <StatusBadge badge="sandbox" />
        <p className="min-w-0 flex-1 text-xs leading-relaxed text-slate-400">
          T104 · SCHUR.PROFILE.BOUND — CHAIN-PARTIAL (two independent arms,
          21/21 + 47/47): the naive margin route is dead, exact spectral-split
          chains close 16/16 with finite data, and the hard core moves to a
          bare_k lower bound plus the soft dressing scalar L. T105 ·
          BARE.AVOIDANCE.CORE — ONE-OF-TWO (28/28): bare is certified in closed
          form and the avoidance law becomes a theorem, leaving one
          Friedrichs-angle statement. Everything after that — the twenty-one
          parts that compress this one statement down to one sign plus one
          declared accounting convention, drive the certified ladder to zone
          155,921, and finally assemble the whole chain end to end (T125 ·
          GRAND.ASSEMBLY — ASSEMBLY-GREEN, 34/34: all five stages on 52 of 52
          zones, 430 completed Cholesky certificates, the load-bearing spine
          96.2% identity-or-Cholesky with the Harnack pair no longer in it) — is
          told in{" "}
          <a
            href="#compression"
            className="text-sky-300 underline decoration-sky-400/30 underline-offset-2 hover:text-sky-200"
          >
            sections 20–22
          </a>
          . Series complete at 125 parts / 3139 sandbox checks; the mandate
          T ≤ 125 is fulfilled. Phase 2 — the full proof — is now open: T126 ·
          UNIFORMITY.SEAMS (SEAMS-CERTIFIED, 31/31) finishes the seam
          architecture, T127 · TWO.INEQUALITIES (BOTH-SHAPED, 28/28) dissects
          the two genuinely new inequalities it left — U5-as-stated is refuted
          and replaced by a band plus an enumeration, U3 collapses to a coarse
          floor —           and T128 · TEML (THREE-OF-FOUR, 27/27) works the resulting
          four-point list cheapest first: three of the four points stand at
          their preregistered bars (the exception list derived and closed, the
          retention bound exact bookkeeping, the boundary-layer exclusion now
          a proof with an 11.6× floor margin), while the kappa bar was missed
          honestly — by 3.6%, systematically in the ratio. T129 ·
          KAPPA.DEEP.SEAMS (KAPPA-WILD, 28/28) is the most productive break
          of the phase: the fitted kappa law falls once on 331 fresh
          transports — bar frozen, violation counted — but the theorem
          underneath stands: flat is exactly 1, linear is exactly 2,
          everything above is curvature, and the curvature chain is a
          per-transport theorem on all 436; the two affordable deep seams
          carry complete certificates on the graded space, honestly
          downgraded with a measured 8% false-positive rate declared before
          the results. Two irreducibles remain (the word &ldquo;for
          all&rdquo;, the RH address); T130 (curvature_bridge_probe.py) is
          running.
        </p>
      </div>

      <p className="mt-3 text-xs leading-relaxed text-slate-500">
        The convergence is a measured typing of where the hardness sits — not
        progress on it. Sandbox; not RH evidence.
      </p>
    </div>
  );
}
