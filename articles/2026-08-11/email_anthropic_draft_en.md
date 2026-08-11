# Draft email — to the Anthropic zeta research team

**Subject:** An exact positive-index-8 form of a registered Weil-positivity target — possible overlap with your moment/inertia method

Dear Anthropic zeta team,

I lead the prime/zeta line of the TFPT project ([fixpoint-theory.com/prime-front](https://fixpoint-theory.com/prime-front)), a machine-verified research program around a Galerkin discretization of Suzuki's localized Weil operator. I'm writing because a recent exact identity on our side lands, verbatim, in the currency of your two-thirds result — and I'd like to ask whether a technical exchange makes sense.

The finding: our program's one open inequality is carried by a fixed 8×8 core M = [[n, b*],[b, B]] per window rung, where (with B ≻ 0) M ≻ 0 ⟺ n − b*B⁻¹b > 0, and we maintain a registered, frozen target n − q ≥ ½μ₁(h) (no-adjust clause, registry-hashed; 67/67 on the surface, first blind holdout passed 28/28 with minimum margin +0.223). The identity: for an explicit congruence scaling A(u_h) of that core, **(tr A)²/‖A‖²_F > 7 ⟺ n − q > ½μ₁(h)** — exactly your rank-trace/positive-index statistic. And because the core is fixed at dimension 8, integrality replaces the asymptotic-population problem: a ratio above 7 forces n₊ = 8, i.e. full positive-definiteness. The B-half premise is certified on our whole reachable surface (B ⪰ 0.5523·I over ideal source objects, exact-rational/interval certificates, 39/39).

To be plain about status: we have no RH proof and claim no progress toward RH. Our reduction chain (cofinal finite positivity ⇒ Weil positivity ⇒ RH) is an exact reformulation at RH strength; the finite head is rigorous (42 interval-certified rungs), the scalar half is open, and our measured map of failed certificate routes says it is hard exactly where your no-go (the ≈0.68 two-moment/bandwidth-1 ceiling, which we've adopted as a route gate) predicts. The registered target is a falsification instrument, not evidence.

Three concrete things we'd propose:

1. Your moment/inertia bounds applied to our fixed core — does the two-moment threshold at dimension 8, or a higher-moment refinement, fall inside machinery you already have?
2. Testing your optimized cosine window inside our Galerkin window family — the dictionary makes "better window" a precise ratio-clearance statement.
3. Comparing Lean formalization approaches — we maintain a Lean 4 tree (against Mathlib) with the extraction-chain hypotheses kernel-checked and the finite-head certificate composition in progress.

I attach a self-contained summary (proven / registered / open, with module identifiers and the full dictionary algebra); everything is reproducible from the public suite. If any of this is of interest, I'd be glad to share data, registries, and code in whatever form is useful.

With appreciation for your work — the formalization especially,

Stefan Hamann
TFPT project — [fixpoint-theory.com/prime-front](https://fixpoint-theory.com/prime-front)

*Attachment: paper_endform_and_inertia_bridge_en.md (research summary, 2026-08-11)*
