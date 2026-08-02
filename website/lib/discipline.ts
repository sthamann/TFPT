// GENERATED FILE -- DO NOT EDIT BY HAND.
// Written by verification/make_discipline_stats.py from the suite's own
// files (status_ledger.csv, v*.py sources, tfpt_research_contracts.tex,
// the committed Lean tree). Single source of truth: the verification suite.
// Regenerate with:  python3 verification/make_discipline_stats.py
// (run automatically by `bash build.sh gen`; freshness enforced by
//  verification/audit_sync.py, so /method can never drift from the suite).
// The replay and anchor entries are certified live by
// verification/v649_discipline_audit.py (META.DISCIPLINE.01) on every
// suite run: python3 verification/run_all.py.

export interface DisciplineStats {
  /** D1: status ledger census. */
  ledger: {
    rows: number;
    active: number;
    retired: number;
    /** Public display classes ([E]/[C]/[O]/[X] + axioms + unclassified). */
    dist: { E: number; C: number; O: number; X: number; AXIOM: number; OTHER: number };
    /** Rows documenting kills / honest negatives / retypes (keyword census). */
    killRows: number;
  };
  /** D2: static census over all verification/v*.py sources. */
  suite: {
    modules: number;
    checkSites: number;
    mustfailOccurrences: number;
    mustfailModules: number;
    eMarks: number;
    cMarks: number;
    seededModules: number;
    sympyModules: number;
  };
  /** D3: preregistration census. */
  contracts: {
    ledgerRows: number;
    withKillCriteria: number;
    texMentions: number;
  };
  /** D4: deterministic replay (certified by v649 on every suite run). */
  replay: {
    sample: string[];
    passed: number;
    total: number;
  };
  /** D5: classical anchors recomputed by v649 + citation census. */
  anchors: {
    recomputed: number;
    citations: Record<string, number>;
    arxivModules: number;
  };
  /** D6: committed Lean proof modules. */
  lean: {
    modules: number;
  };
  /** D7: open gates (display class O, active). */
  openGates: {
    count: number;
  };
}

export const DISCIPLINE: DisciplineStats =
{
  "ledger": {
    "rows": 715,
    "active": 710,
    "retired": 5,
    "dist": {
      "E": 169,
      "C": 222,
      "O": 282,
      "X": 37,
      "AXIOM": 2,
      "OTHER": 3
    },
    "killRows": 214
  },
  "suite": {
    "modules": 643,
    "checkSites": 6040,
    "mustfailOccurrences": 1469,
    "mustfailModules": 233,
    "eMarks": 3579,
    "cMarks": 1323,
    "seededModules": 77,
    "sympyModules": 396
  },
  "contracts": {
    "ledgerRows": 47,
    "withKillCriteria": 29,
    "texMentions": 20
  },
  "replay": {
    "sample": [
      "v55_coxeter_cycle",
      "v91_spine_tetrahedron",
      "v108_pascal_ladder",
      "v205_xi_threequarter",
      "v305_witness_independence",
      "v405_seam_equiv_omega",
      "v505_celestial_wp5e_beta_equivariant_ledger",
      "v555_pareto_tv_identities",
      "v600_joint_embedding",
      "v634_st31_structure",
      "v643_w1_theorem",
      "v648_sign_uncertainty"
    ],
    "passed": 12,
    "total": 12
  },
  "anchors": {
    "recomputed": 5,
    "citations": {
      "Jacobi": 17,
      "Hecke": 8,
      "Construction A": 2,
      "Suzuki": 16,
      "Weil": 34,
      "Eisenstein": 23
    },
    "arxivModules": 40
  },
  "lean": {
    "modules": 54
  },
  "openGates": {
    "count": 278
  }
};

/** Replay verdict enum of v649 (frozen). */
export const DISCIPLINE_VERDICT = "DISCIPLINE-CERTIFIED";
