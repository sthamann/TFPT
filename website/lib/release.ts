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
  releaseDate: "2026-08-11",
};

export const RELEASE_ASSETS: Record<string, ReleaseAsset> = {
  "/papers/introduction.pdf": {
    href: "/papers/introduction.pdf",
    ...COMMON,
    bytes: 4594235,
    sha256:
      "25dd310070ec177eb9d99d696f50c6e595da8913350f7e2f6b1c9cc1415d8f07",
    changelog:
      "Compiler-closure reading guide: two axioms, the dependency DAG, the predictions and the proof ledger.",
  },
  "/papers/tfpt_1_architecture_e8.pdf": {
    href: "/papers/tfpt_1_architecture_e8.pdf",
    ...COMMON,
    bytes: 1486305,
    sha256:
      "5e20361c27878ed54ea3d7d067bc040c1c29f66690072e0f3e96f87fe723c9a8",
    changelog:
      "Architecture: the two axioms, the D₅ × A₃ → E₈ construction, and the EM fixed point with existence + uniqueness.",
  },
  "/papers/tfpt_2_standard_model.pdf": {
    href: "/papers/tfpt_2_standard_model.pdf",
    ...COMMON,
    bytes: 1328735,
    sha256:
      "34bb8740ed9fd7c08e0af136125214fe585d716c20950bac975a1ac45eda231d",
    changelog:
      "The Standard Model in one φ₀-ladder, the flavor residue matrix, and the derived solar angle θ₁₂.",
  },
  "/papers/tfpt_3_e8_audit_bootstrap.pdf": {
    href: "/papers/tfpt_3_e8_audit_bootstrap.pdf",
    ...COMMON,
    bytes: 1878927,
    sha256:
      "7c9e05432c4885104a8fc12ddaec4d2a35d09da608ea485caa000c56be8bc0c7",
    changelog:
      "The seven E₈ slices as an audit raster, the cascade bridge, and the Möbius bootstrap.",
  },
  "/papers/tfpt_4_frontier.pdf": {
    href: "/papers/tfpt_4_frontier.pdf",
    ...COMMON,
    bytes: 825795,
    sha256:
      "ffa420cfb8274bef89ac609eb4d087f4602d765f4381599f5f3f5df2099e2052",
    changelog:
      "Honest status of η_B, m_p/m_e, Koide, dark matter and full quantum gravity — not forced onto the ladder.",
  },
  "/papers/tfpt_horizon_readouts.pdf": {
    href: "/papers/tfpt_horizon_readouts.pdf",
    ...COMMON,
    bytes: 932362,
    sha256:
      "ec0e44b97a1e77626c4b88787e86442430e6d7c345ff1f145805a28174930bde",
    changelog:
      "Appendix H — the horizon unit system: c₃ = 1/(8π) as the universal horizon thermal code.",
  },
  "/papers/origin_theory.pdf": {
    href: "/papers/origin_theory.pdf",
    ...COMMON,
    bytes: 1136168,
    sha256:
      "4e900c14f29778f663e57bdaa68de119ecc419bcb4397f0231febc3cd2d7f633",
    changelog:
      "Origin Theory: the (5,3) skeleton, the triply-forced 8, the order-30 Coxeter cycle, and the gapped unique attractor.",
  },
  "/papers/tfpt_research_contracts.pdf": {
    href: "/papers/tfpt_research_contracts.pdf",
    ...COMMON,
    bytes: 2290282,
    sha256:
      "0d958c2e964f791727b60befb0dcd274172030f88d257d084b5712a7fd7d2554",
    changelog:
      "Research contracts for the remaining interfaces: v_geo (scale anchor, [O]), G_net (metric inclusion, [C] closed modulo cited theorems), F_transfer (downstream transport, typed runnable solvers [C]).",
  },
  "/papers/tfpt_safeguards.pdf": {
    href: "/papers/tfpt_safeguards.pdf",
    ...COMMON,
    bytes: 618281,
    sha256:
      "a44075bd8a9bb65c544e39326306267bc0c58c93c8653abd51a02f979219d421",
    changelog:
      "Safeguards: the verification discipline — the status calculus, no-free-pattern + reverse audit (and v431: the 5/8 'overhead' degrees are the forced two-family ladder 6·spine ⊔ det-ladder, not diffuse slack), the over-determination map (v427) with its honest self-correction (v428: the seven arithmetic witnesses compress one (2,3,5)/E₈ object; the genuine multiplication is the input forced four ways + the foreign α⁻¹) and its unconditional floor (v432: ~10⁻¹⁰ from disjoint pieces only, ~20 orders above the v100 conditional; hardened by v436 to an assumption-minimal 1/94,500 ≈ 4.40σ counting floor with a monotone concession ladder), the firewall + No-Unit theorem, frozen predictions + null model, the independent Wolfram and Lean paths, and the red team.",
  },
  "/papers/tfpt_5_redteam.pdf": {
    href: "/papers/tfpt_5_redteam.pdf",
    ...COMMON,
    bytes: 864013,
    sha256:
      "d193955b0a9a0ec6c6c3132bb536f2e1bbdcc40e759845535469f065a115498f",
    changelog:
      "The adversarial audit: Targets A–E, what survives, what each target reduces to, and the kill tests.",
  },
  "/papers/note_e8_gaussian_code.pdf": {
    href: "/papers/note_e8_gaussian_code.pdf",
    ...COMMON,
    bytes: 347652,
    sha256:
      "16bb889275621a7642d42e3746b43de77703d7e7bda0515c495c81d5f65ea36c",
    changelog:
      "Working note N1 — the Gaussian code bridge: E₈ over ℤ[i] via Construction A over the extended Hamming code, the canonical four-bit quotient L/(1+i)L ≅ F₂⁴, and the G₃₁ quartic companion (v689/v690 + 65 Lean theorems).",
  },
  "/papers/note_hilbert_polya_truncations.pdf": {
    href: "/papers/note_hilbert_polya_truncations.pdf",
    ...COMMON,
    bytes: 470497,
    sha256:
      "6ff5c5fd7da1e14d73c1e8441cca7288640bca473946b84effa014988b7dc3dd",
    changelog:
      "Working note N2 — a computable, zeta-free truncation family for the Weil measure: measurements on a Hilbert–Pólya candidate (v714, v716–v721, v727–v734; prediction freeze; no RH claim).",
  },
  "/papers/changelog.pdf": {
    href: "/papers/changelog.pdf",
    ...COMMON,
    bytes: 2733788,
    sha256:
      "271fe6be233b38a9d786bf603001ec5fd073b49bc5c1453d97fd184d31bcb7f9",
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
