"use client";

import { motion } from "motion/react";

/** Twenty-three promoted modules; check counts from verification/script_registry.csv. */
const MODULES = [
  { id: "v535", claim: "HECKE.GEOM.01", title: "Hecke from geometry", checks: 25 },
  {
    id: "v536",
    claim: "HECKE.GEOM.EICHLER.01",
    title: "Eichler trace layer",
    checks: 23,
  },
  {
    id: "v537",
    claim: "HECKE.GEOM.HALFINT.01",
    title: "Half-integral bridge",
    checks: 20,
  },
  {
    id: "v538",
    claim: "HECKE.GEOM.RTF.01",
    title: "Relative-trace identity",
    checks: 18,
  },
  {
    id: "v539",
    claim: "RTF.GNS.WEIL.01",
    title: "Weil structure · two obstructions",
    checks: 25,
  },
  {
    id: "v540",
    claim: "RTF.GNS.AMP.01",
    title: "Amplitude route · linear carrier",
    checks: 34,
  },
  {
    id: "v541",
    claim: "RTF.GNS.LEDGER.01",
    title: "Matching lemma · transport ledger",
    checks: 33,
  },
  {
    id: "v542",
    claim: "PRIME.MARGIN.IDENT.01",
    title: "Margin-chain identities · phase 2",
    checks: 44,
  },
  {
    id: "v543",
    claim: "PRIME.MMATRIX.IDENT.01",
    title: "Lumped M-matrix pair · phase 2",
    checks: 35,
  },
  {
    id: "v544",
    claim: "PRIME.LONGLAG.SUPP.01",
    title: "Long-lag support structure · phase 2",
    checks: 24,
  },
  {
    id: "v545",
    claim: "PRIME.HARDY.IDENT.01",
    title: "Hardy-core identities · phase 2",
    checks: 36,
  },
  {
    id: "v546",
    claim: "PRIME.CAPCHAIN.IDENT.01",
    title: "Capacity-chain identities · phase 2",
    checks: 26,
  },
  {
    id: "v547",
    claim: "PRIME.LEVEL.LEMMA.01",
    title: "Level-lemma identities · phase 2",
    checks: 20,
  },
  {
    id: "v548",
    claim: "PRIME.GREEN.SZEGO.IDENT.01",
    title: "Green/Szegő identities · phase 2",
    checks: 21,
  },
  {
    id: "v549",
    claim: "PRIME.GAUGE.PARITY.IDENT.01",
    title: "Gauge/parity identities · phase 2",
    checks: 21,
  },
  {
    id: "v550",
    claim: "PRIME.ODD.SECTOR.IDENT.01",
    title: "Odd-sector identities · phase 2",
    checks: 19,
  },
  {
    id: "v551",
    claim: "PRIME.RITZ.CEIL.01",
    title: "Fixed-size Ritz ceiling certificate · phase 2",
    checks: 16,
  },
  {
    id: "v552",
    claim: "PRIME.ANGLE.INSTR.01",
    title: "Four fixed-size angle instruments · phase 2",
    checks: 18,
  },
  {
    id: "v553",
    claim: "PRIME.EXACT.FORM.IDENT.01",
    title: "Exact-form identities · phase 2",
    checks: 25,
  },
  {
    id: "v554",
    claim: "PRIME.SAMPLING.HARM.01",
    title: "Sampling/harmonics identities · phase 2",
    checks: 21,
  },
  {
    id: "v555",
    claim: "PRIME.PARETO.TV.01",
    title: "Pareto/total-variation identities · phase 2",
    checks: 23,
  },
  {
    id: "v556",
    claim: "PRIME.GAUGE.PPR.01",
    title: "Gauge/P_pr identities · phase 2",
    checks: 23,
  },
  {
    id: "v557",
    claim: "PRIME.CASCADE.VECT.01",
    title: "Cascade/vector identities · phase 2",
    checks: 23,
  },
  {
    id: "v558",
    claim: "PRIME.BILINEAR.RANK.01",
    title: "Bilinear/rank identities · phase 2",
    checks: 27,
  },
] as const;

const TOTAL_CHECKS = MODULES.reduce((s, m) => s + m.checks, 0);
const MAX_CHECKS = Math.max(...MODULES.map((m) => m.checks));
const SANDBOX_PROBES = 170;
const SANDBOX_CHECKS = 4514;

export function ModuleLadder() {
  return (
    <div className="rounded-2xl border border-slate-700/50 bg-slate-950/60 p-4 sm:p-5">
      <div className="mb-4 flex flex-wrap items-center justify-between gap-2">
        <p className="font-mono text-[10px] uppercase tracking-widest text-emerald-300/90">
          Load-bearing modules · checks per module
        </p>
        <span className="font-mono text-[10px] text-slate-500">
          {TOTAL_CHECKS} checks · {MODULES.length} modules
        </span>
      </div>

      <ul className="space-y-2">
        {MODULES.map((m, i) => (
          <li key={m.id}>
            <div className="flex items-baseline justify-between gap-2">
              <span className="font-mono text-[11px] text-emerald-200">
                {m.id}
              </span>
              <span className="min-w-0 flex-1 truncate text-xs text-slate-300">
                {m.title}
              </span>
              <span className="font-mono text-[10px] text-slate-500">
                {m.checks}
              </span>
            </div>
            <div className="mt-1 h-1.5 overflow-hidden rounded-full bg-slate-800">
              <motion.div
                className="h-full rounded-full bg-gradient-to-r from-emerald-500/70 to-emerald-300/70"
                initial={{ width: 0 }}
                whileInView={{ width: `${(m.checks / MAX_CHECKS) * 100}%` }}
                viewport={{ once: true, amount: 0.5 }}
                transition={{ duration: 0.6, delay: i * 0.07, ease: "easeOut" }}
              />
            </div>
            <p className="mt-0.5 font-mono text-[9px] uppercase tracking-wider text-slate-600">
              {m.claim}
            </p>
          </li>
        ))}
      </ul>

      <div className="mt-4 grid grid-cols-3 gap-2 text-center">
        <Stat value={String(TOTAL_CHECKS)} label="load-bearing checks" tone="emerald" />
        <Stat
          value={`${SANDBOX_CHECKS}`}
          label={`sandbox checks · ${SANDBOX_PROBES} probes`}
          tone="amber"
        />
        <Stat value="1" label="object remains · I5" tone="violet" />
      </div>

      <p className="mt-3 text-xs leading-relaxed text-slate-500">
        Sandbox probes never move a marker; only the twenty-four modules above
        are cited in the papers and the ledger. I5 is an equivalence typing (⟺ Weil
        positivity ⟺ RH), not a proof claim.
      </p>
    </div>
  );
}

const TONE = {
  emerald: "border-emerald-400/30 bg-emerald-500/10 text-emerald-100",
  amber: "border-amber-400/30 bg-amber-500/10 text-amber-100",
  violet: "border-violet-400/30 bg-violet-500/10 text-violet-100",
} as const;

function Stat({
  value,
  label,
  tone,
}: {
  value: string;
  label: string;
  tone: keyof typeof TONE;
}) {
  return (
    <div className={`rounded-xl border px-2 py-2 ${TONE[tone]}`}>
      <p className="font-mono text-lg leading-none">{value}</p>
      <p className="mt-1 text-[10px] leading-tight text-slate-400">{label}</p>
    </div>
  );
}
