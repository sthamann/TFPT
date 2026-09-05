/**
 * Release metadata for the published PDFs. The version, release date, byte
 * size, and SHA-256 hash are written here so reviewers can verify integrity
 * and tell two snapshots apart without opening the file.
 *
 * To regenerate these values after a re-build:
 *
 *   npm run release:write   (rewrites bytes + sha256 in place)
 *
 * which invokes scripts/release-hashes.mjs against website/public.
 */

export interface ReleaseAsset {
  /** Path under /public, including the leading slash. */
  href: string;
  /** Series version this snapshot was built against. */
  version: string;
  /** ISO-8601 release date for the artifact. */
  releaseDate: string;
  /** Byte size of the PDF (raw bytes). */
  bytes: number;
  /** SHA-256 hash of the PDF, lowercase hex. */
  sha256: string;
  /** Optional changelog entry — short, one-line. */
  changelog?: string;
}

const COMMON = {
  version: "TFPT 5.4",
  releaseDate: "2026-09-05",
};

export const RELEASE_ASSETS: Record<string, ReleaseAsset> = {
  "/papers/introduction.pdf": {
    href: "/papers/introduction.pdf",
    ...COMMON,
    bytes: 5471443,
    sha256:
      "0d943193a17eda9ebad97687e7f2c39929e12551fb0c4663c763b96658decab4",
    changelog:
      "Compiler-closure reading guide: two axioms, the dependency DAG, the predictions and the proof ledger.",
  },
  "/papers/tfpt_1_architecture_e8.pdf": {
    href: "/papers/tfpt_1_architecture_e8.pdf",
    ...COMMON,
    bytes: 1495359,
    sha256:
      "2b2eb1b667ad38219b56c017b8391ba999e34bdd4431327a447f567c95c7e387",
    changelog:
      "Architecture: the two axioms, the D₅ × A₃ → E₈ construction, and the EM fixed point with existence + uniqueness.",
  },
  "/papers/tfpt_2_standard_model.pdf": {
    href: "/papers/tfpt_2_standard_model.pdf",
    ...COMMON,
    bytes: 1343289,
    sha256:
      "fbe1edd14e615d9266f1a1c7936bcef370d9bc6b31c03c03694124345bdc7ca9",
    changelog:
      "The Standard Model in one φ₀-ladder, the flavor residue matrix, and the derived solar angle θ₁₂.",
  },
  "/papers/tfpt_3_e8_audit_bootstrap.pdf": {
    href: "/papers/tfpt_3_e8_audit_bootstrap.pdf",
    ...COMMON,
    bytes: 2049112,
    sha256:
      "49fd9aeef1104a44fc9a760c6b039039a8c806dbbb9fb45b269596cb75e83f8f",
    changelog:
      "The seven E₈ slices as an audit raster, the cascade bridge, and the Möbius bootstrap.",
  },
  "/papers/tfpt_4_frontier.pdf": {
    href: "/papers/tfpt_4_frontier.pdf",
    ...COMMON,
    bytes: 851025,
    sha256:
      "c153fc9239c62b5136730140f31f65f4a130c90e88b81a94710d5e7884d698c2",
    changelog:
      "Honest status of η_B, m_p/m_e, Koide, dark matter and full quantum gravity — not forced onto the ladder.",
  },
  "/papers/tfpt_horizon_readouts.pdf": {
    href: "/papers/tfpt_horizon_readouts.pdf",
    ...COMMON,
    bytes: 934667,
    sha256:
      "0121387c4f5297562acbcb384dc82307b4ab352ae658a4188eca0b7639e03cd1",
    changelog:
      "Appendix H — the horizon unit system: c₃ = 1/(8π) as the universal horizon thermal code.",
  },
  "/papers/origin_theory.pdf": {
    href: "/papers/origin_theory.pdf",
    ...COMMON,
    bytes: 1146770,
    sha256:
      "e7b2db2904c442d8c8ce807421290eeae079e371dc9439bc24a5ff41b12bc864",
    changelog:
      "Origin Theory: the (5,3) skeleton, the triply-forced 8, the order-30 Coxeter cycle, and the gapped unique attractor.",
  },
  "/papers/tfpt_research_contracts.pdf": {
    href: "/papers/tfpt_research_contracts.pdf",
    ...COMMON,
    bytes: 2677318,
    sha256:
      "22b794fad8551f017412feb14298a071f53274c2b20462c9bb21eef62533d390",
    changelog:
      "Research contracts separate compiler Rest = v_geo ⊕ G_net ⊕ F_transfer from Rest_TOE. Round 3 (v1022–v1025) adds finite certificates and obstruction results for T2/T4/T6–T8 but closes no T-gate: T1–T8 and TFPT.TOE.COMPLETE.01 remain [O], and no shared complete 3+1D parent is constructed. v1022/v1025 require native python-flint==0.9.0. Wave 4: v998–v1001 + lattice-fundamental quasilocal-family amendment + ALPHA relative-det note.",
  },
  "/papers/tfpt_safeguards.pdf": {
    href: "/papers/tfpt_safeguards.pdf",
    ...COMMON,
    bytes: 619366,
    sha256:
      "b5a71660e24e98c26a8384ba69b6fcd81988e173ad7b16393bf5552c85bb8561",
    changelog:
      "Safeguards: the verification discipline — the status calculus, no-free-pattern + reverse audit (and v431: the 5/8 'overhead' degrees are the forced two-family ladder 6·spine ⊔ det-ladder, not diffuse slack), the over-determination map (v427) with its honest self-correction (v428: the seven arithmetic witnesses compress one (2,3,5)/E₈ object; the genuine multiplication is the input forced four ways + the foreign α⁻¹) and its unconditional floor (v432: ~10⁻¹⁰ from disjoint pieces only, ~20 orders above the v100 conditional; hardened by v436 to an assumption-minimal 1/94,500 ≈ 4.40σ counting floor with a monotone concession ladder), the firewall + No-Unit theorem, frozen predictions + null model, the independent Wolfram and Lean paths, and the red team.",
  },
  "/papers/tfpt_5_redteam.pdf": {
    href: "/papers/tfpt_5_redteam.pdf",
    ...COMMON,
    bytes: 870920,
    sha256:
      "0955c94ceb2ede8b636035564c9f1c4bfb33c64e6060c81737454788257ace08",
    changelog:
      "The adversarial audit: Targets A–E, what survives, what each target reduces to, and the kill tests.",
  },
  "/papers/note_e8_gaussian_code.pdf": {
    href: "/papers/note_e8_gaussian_code.pdf",
    ...COMMON,
    bytes: 370071,
    sha256:
      "04f4e876bc9ad29bf4794fdce5e04d1e985b35372c1811a36397a71c4c08cdb9",
    changelog:
      "Working note N1 — the Gaussian code bridge: E₈ over ℤ[i] via Construction A over the extended Hamming code, the canonical four-bit quotient L/(1+i)L ≅ F₂⁴, and the G₃₁ quartic companion (v689/v690 + 65 Lean theorems).",
  },
  "/papers/note_hilbert_polya_truncations.pdf": {
    href: "/papers/note_hilbert_polya_truncations.pdf",
    ...COMMON,
    bytes: 470497,
    sha256:
      "dfa91acbfc0ebb63e188ca95f1d16ce108368346cd8889b702e195c61712c64e",
    changelog:
      "Working note N2 — a computable, zeta-free truncation family for the Weil measure: measurements on a Hilbert–Pólya candidate (v714, v716–v721, v727–v734; prediction freeze; no RH claim).",
  },
  "/papers/changelog.pdf": {
    href: "/papers/changelog.pdf",
    ...COMMON,
    bytes: 3318678,
    sha256:
      "0529060b15d43678ffc2f9b97a2d8435f7422afc4106035373dd3e6f6bd2bcd2",
    changelog:
      "The canonical dated changelog of every change to the theory, the suite, the papers and the website.",
  },
};

export function getReleaseAsset(href: string): ReleaseAsset | undefined {
  return RELEASE_ASSETS[href];
}

export function formatBytes(bytes: number): string {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(0)} KB`;
  return `${(bytes / (1024 * 1024)).toFixed(2)} MB`;
}

export function formatHashShort(sha256: string): string {
  return `${sha256.slice(0, 8)}…${sha256.slice(-4)}`;
}
