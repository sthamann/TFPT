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
          T104 · SCHUR.PROFILE.BOUND — CHAIN-PARTIAL: the hard core moved to
          a bare_k lower bound + the soft dressing scalar L. T105 ·
          BARE.AVOIDANCE.CORE — ONE-OF-TWO (28/28): bare is certified in
          closed form and the avoidance law is a theorem. T106 ·
          FRIEDRICHS.ANGLE — DENSITY-MAPPED (32/32): the parity split closes
          the even channel 16/16 and localizes all remaining hardness in the
          odd channel — one Loewner statement on half the dimensions. T107 ·
          ODD.CHANNEL.CLOSURE — SCALAR-TRACTABLE (30/30): (R) ⟺ one scalar
          ratio r = κ/ε ≤ 1, measured r = 0.005…0.18 — two orders of magnitude
          of room; the symbol route is structurally dead. T108 ·
          RATIO.CERTIFICATE — EPSILON-IDENTITY (44/44): ε&apos;s positivity is
          an exact identity (Szegő pivot / Q-energy) coinciding with the
          induction positivity itself; (R) is down to two scalars — what
          remains is literally one boundary value of an explicit vector. T109 ·
          BOUNDARY.DECAY — BOUNDARY-CERTIFIED (29/29): both scalars certified —
          ω unconditionally via a graded matrix cap, the boundary value via a
          residual certificate that carries the cancellation; the chain closes
          16/16 on exactly one strict-margin input, 10²–10⁶ weaker than the
          conclusion. T110 · MARGIN.PROPAGATION — MARGIN-PROPAGATES (28/28):
          the induction circle closes end-to-end on the measured zones —
          certified base case, 15 certified handover steps (graded Loewner
          minorant, retention 1.0000), atom entry structurally free; three
          sharp gaps remain (no reserve, no scalar step law, no k-uniformity).
          T111 · DEEP.ZONE.STRESS — CROSSING-CONFIRMED (23/23): the deep
          ladder (199 zones, 117 handovers to n = 521) measures the crossing
          at n* ≈ 462 (upper bound; the n ≈ 170 extrapolation was a fit
          artefact) and splits the wall into three — margin wall n ≈ 462,
          twin-prime ladder wall 521→523, requirement wall n = 727 — while
          the mechanism never fails (117/117 at retention 1.000000): depth,
          not n, is the operating variable. T112 · ADAPTIVE.SCALING —
          SCALING-PARTIAL (20/20): in the gap-coupled frame D_k = g_k/(2ν)
          two walls fall structurally (ladder death and the depth dependence
          of the requirement wall; 461→463 certifies with the reserve
          opening ~1000×), but the margin wall is frame-invariant at
          exponent −0.974 — the substance of the requirement, not geometry;
          the hardness is now one limit operator plus one convergence rate,
          with the prime-gap dependence disclosed. T113 (
          <span className="font-mono text-slate-300">
            limit_operator_probe.py
          </span>
          ) running.
        </p>
      </div>

      <p className="mt-3 text-xs leading-relaxed text-slate-500">
        The convergence is a measured typing of where the hardness sits — not
        progress on it. Sandbox; not RH evidence.
      </p>
    </div>
  );
}
