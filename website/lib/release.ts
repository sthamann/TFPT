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
  releaseDate: "2026-08-06",
};

export const RELEASE_ASSETS: Record<string, ReleaseAsset> = {
  "/papers/introduction.pdf": {
    href: "/papers/introduction.pdf",
    ...COMMON,
    bytes: 4144902,
    sha256:
      "6d375a31c56d20d0d580017d0b0e98e9643582b60305f94381d9c4d1c822ee43",
    changelog:
      "Compiler-closure reading guide: two axioms, the dependency DAG, the predictions and the proof ledger.",
  },
  "/papers/tfpt_1_architecture_e8.pdf": {
    href: "/papers/tfpt_1_architecture_e8.pdf",
    ...COMMON,
    bytes: 1372480,
    sha256:
      "b9f489b88cf6aec422e47302fdb360fa715c45a5ccec2e84a96f85c9a3ffa844",
    changelog:
      "Architecture: the two axioms, the D₅ × A₃ → E₈ construction, and the EM fixed point with existence + uniqueness.",
  },
  "/papers/tfpt_2_standard_model.pdf": {
    href: "/papers/tfpt_2_standard_model.pdf",
    ...COMMON,
    bytes: 1312734,
    sha256:
      "f5bcee8d3b28cc0d8ccea5b2bb94333d6cf857cb3a363dc4690f02c55713027f",
    changelog:
      "The Standard Model in one φ₀-ladder, the flavor residue matrix, and the derived solar angle θ₁₂.",
  },
  "/papers/tfpt_3_e8_audit_bootstrap.pdf": {
    href: "/papers/tfpt_3_e8_audit_bootstrap.pdf",
    ...COMMON,
    bytes: 1630531,
    sha256:
      "c95b97228f62b6544e6386c6820a3f269591456522fbf28c0ecafc51a0957482",
    changelog:
      "The seven E₈ slices as an audit raster, the cascade bridge, and the Möbius bootstrap.",
  },
  "/papers/tfpt_4_frontier.pdf": {
    href: "/papers/tfpt_4_frontier.pdf",
    ...COMMON,
    bytes: 825879,
    sha256:
      "13c041354c25da62144cda3c6ee245c8d62a172d8a5a94a96c42d85b569d1a86",
    changelog:
      "Honest status of η_B, m_p/m_e, Koide, dark matter and full quantum gravity — not forced onto the ladder.",
  },
  "/papers/tfpt_horizon_readouts.pdf": {
    href: "/papers/tfpt_horizon_readouts.pdf",
    ...COMMON,
    bytes: 932723,
    sha256:
      "7849b10243e085735b73bd1d7ae8c72ea0456097ada4b8b9ed753186f0073cd9",
    changelog:
      "Appendix H — the horizon unit system: c₃ = 1/(8π) as the universal horizon thermal code.",
  },
  "/papers/origin_theory.pdf": {
    href: "/papers/origin_theory.pdf",
    ...COMMON,
    bytes: 1136242,
    sha256:
      "a8a149763ac138f2f585cd0ac46f256bc6a14b5ebe188b9e32e56b34dca24f4e",
    changelog:
      "Origin Theory: the (5,3) skeleton, the triply-forced 8, the order-30 Coxeter cycle, and the gapped unique attractor.",
  },
  "/papers/tfpt_research_contracts.pdf": {
    href: "/papers/tfpt_research_contracts.pdf",
    ...COMMON,
    bytes: 2048662,
    sha256:
      "75a175b17a823f3204ad8f85315de296f7a147a80a88ccd3946519d3d9b3fae6",
    changelog:
      "Research contracts for the remaining interfaces: v_geo (scale anchor, [O]), G_net (metric inclusion, [C] closed modulo cited theorems), F_transfer (downstream transport, typed runnable solvers [C]).",
  },
  "/papers/tfpt_safeguards.pdf": {
    href: "/papers/tfpt_safeguards.pdf",
    ...COMMON,
    bytes: 618897,
    sha256:
      "0df8df1d06a3882c6d07c60fc818ce2188f6bc89285e614a74e018c1b3805760",
    changelog:
      "Safeguards: the verification discipline — the status calculus, no-free-pattern + reverse audit (and v431: the 5/8 'overhead' degrees are the forced two-family ladder 6·spine ⊔ det-ladder, not diffuse slack), the over-determination map (v427) with its honest self-correction (v428: the seven arithmetic witnesses compress one (2,3,5)/E₈ object; the genuine multiplication is the input forced four ways + the foreign α⁻¹) and its unconditional floor (v432: ~10⁻¹⁰ from disjoint pieces only, ~20 orders above the v100 conditional; hardened by v436 to an assumption-minimal 1/94,500 ≈ 4.40σ counting floor with a monotone concession ladder), the firewall + No-Unit theorem, frozen predictions + null model, the independent Wolfram and Lean paths, and the red team.",
  },
  "/papers/tfpt_5_redteam.pdf": {
    href: "/papers/tfpt_5_redteam.pdf",
    ...COMMON,
    bytes: 864575,
    sha256:
      "9fa75fd434d1ed5c87d875f711cbfc9c5dabcd3cda516abfbcabe1c5f04cf083",
    changelog:
      "The adversarial audit: Targets A–E, what survives, what each target reduces to, and the kill tests.",
  },
  "/papers/note_e8_gaussian_code.pdf": {
    href: "/papers/note_e8_gaussian_code.pdf",
    ...COMMON,
    bytes: 345119,
    sha256:
      "43e9620e4b6dc52f3c8e7f3e73d8f82ac5d89ab86bb2311a33f6c8a763c68428",
    changelog:
      "Working note N1 — the Gaussian code bridge: E₈ over ℤ[i] via Construction A over the extended Hamming code, the canonical four-bit quotient L/(1+i)L ≅ F₂⁴, and the G₃₁ quartic companion (v689/v690 + 65 Lean theorems).",
  },
  "/papers/note_hilbert_polya_truncations.pdf": {
    href: "/papers/note_hilbert_polya_truncations.pdf",
    ...COMMON,
    bytes: 470497,
    sha256:
      "2fc9f3144b2f7c7581cdcdb0b1ae762828a1176fec054a09bcbeeff0a5652d4e",
    changelog:
      "Working note N2 — a computable, zeta-free truncation family for the Weil measure: measurements on a Hilbert–Pólya candidate (v714, v716–v721, v727–v734; prediction freeze; no RH claim).",
  },
  "/papers/changelog.pdf": {
    href: "/papers/changelog.pdf",
    ...COMMON,
    bytes: 2415259,
    sha256:
      "9f95968b7ea060cf8762e20351553c04f9f39cbd95d040ff99410bfcdd45ce13",
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
