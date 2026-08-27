#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_inventory.py -- regenerate rh/INVENTORY.json.

The inventory is the machine-readable register of every RH-relevant
file in the repository that the rh/ workspace REFERENCES (it never
copies them).  Pinned entries carry the SHA-256 of the file content
at inventory time; run_rh.py uses them as a drift detector.  Living
documents (ledger, program paper, notes) are registered unpinned
(sha256 recorded for information, drift reported as INFO only).

Usage:  python rh/verification/make_inventory.py   (from repo root)
"""

import hashlib
import json
import os
import sys
from datetime import date

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(REPO, "rh", "INVENTORY.json")

VER = "verification"
EXP = "experiments/tfpt-discovery"

_R285_STATUS = "L* decomposition round (33/33): L* split candidate-wise into (D) the per-atom diagonal condition v_k K_{N_w}(y_k) < 1 - delta and (C) the coherence condition assist = lambda_max/maxdiag - 1 < delta', with EXACT bookkeeping (multiplicative identity lambda = maxdiag x (1+assist), budget equivalence margin > 0 <=> rho < 1 gated at every degree; additive Weyl form measured coarse, slack 0.287).  VERDICT DECOMPOSITION_EXACT(ladder budgets 42/42: maxdiag_Nw 0.9343..0.9941 with sp(N, maxdiag) +0.65 -- the per-atom budget TIGHTENS toward 1 with N; assist_Nw 0.0059..0.0702 with sp -0.65 -- the coherent assist SHRINKS: no fixed delta/delta' budget supported; controls rip at (C): at their crossing maxdiag 0.26..0.38, assist 1.69..2.99, offdiag share 0.63..0.75) + PERATOM_CLASSICAL(binding atom = the shallow-u ARCH edge below the first prime on 42/42 rungs, w9 top-5 = folds 2/4/6/8/10 near-degenerate; classical arcsine yardstick typed COMPARISON_ONLY: measured growth at the binding atom p = 0.38 is SUB-classical (bulk law 1), n*_bind = 58 vs n_DIAG = 187 (1.68 octaves), coverage 0/42 -- (D) survives BECAUSE the discrete kernel grows sub-classically: the provability route for (D) must be a discrete bound) + the b3 SURPRISE: (D) is NOT world-blind at window scale (controls' n_DIAG = 50..91 << 184 vs MAIN 187 > 185: (D) itself separates; the clean split 'D blind + C carries all arithmetic' is REFUTED as stated, while at the DEATH degree the split is clean) => sealed COHERENCE_UNSTRUCTURED by the blind_D clause (the assist itself separates: MAIN 0.0195 at crossing vs controls 1.69..2.99, K_D2 MAIN_SEPARATING; z-scores: MAIN z_v = -3.15 destructive vs controls +4.95..+12.44 constructive -- sign-clean, band-MIXED) + ENSEMBLE_POSITION(MAIN's wall assist is the extreme LOW outlier of both conservation-gated ensembles (sign 16 reps, scramble 12 reps, pct 0.00): the near-single-atom wall is MAIN-specific; early coherence scramble-generic pct 0.75).  mp dps-60 ward incl. maxdiag_184; NO L* claim, NO RH CLAIM"

_R286_STATUS = "L* margin-scaling round (29/29): the decay law of the spectral L* margin 1 - lambda_max(E_{N_w}) and the large-S counterexample hunt, run on the STANDALONE document pipeline (rh/problem/verify_lstar_instance.py imported verbatim; flagship AND every new anchor builder-cross-checked <= 2.9e-15).  VERDICT MARGIN_LAW(POWER preferred, MAD 0.608 vs exp 0.850; alpha = 3.05 vs S / 3.06 vs N_w / 2.39 vs z; halves 3.67/2.07, curvature -1.60 -- the law FLATTENS at large S, the honest anti-extrapolation flag; driver: the decay is the cancellation sharpening c_w -> 1 (exact identity margin = (1-maxdiag)(1-c_w), c_w range 0.9944..0.999991; slope of log(1-c_w) -2.377 vs the diagonal's -0.693); extremal shallow-edge band persists 42/42, PR median 1.81) + NO_CROSSING_EXTENDED(the sealed extension rule -- the document admissibility rule with the window cap lifted, z^2 <= table cap -- admits 83 candidates; the first 15 by (N_w, kz): N_w 942..1218, S 1883..2435, z 101..547, ALL margins sign-safe positive under the sealed tier protocol (f64 two-route equivalence + mp dps-30/45 staggered on every sub-1e-5 rung: 32/42 + 15/15 mp-verified, worst f64-vs-mp rel dev 6.6e-6): min +1.806e-8 at kz=119 (z=529), max +3.314e-7 -- NO counterexample, the triple-verification protocol never fired) + QN_RECONCILED(the r258 terminal bordered ratio q_N measured on the same 57 rungs: q_N < 1 on 57/57 with min(1-q_N) = 0.019 and spearman(1-q_N, N_w) = +0.35 = FLAT (r258 reproduced) while the spectral margin falls 4.0 decades (spearman -0.87) -- TWO DIFFERENT OBJECTS; exact Rayleigh lambda_max >= ||B[:, N_w-1]||^2 on 42/42 (median terminal Rayleigh 0.448, far below lambda ~0.999998); degree sector: the extremal coefficient vector is FULL-SECTOR (terminal-degree share median 1.9e-8, top-decile 0.002, degree PR median 188) -- the terminal pivot cannot see the sector the spectral margin lives in; HARMLESS adjudication (calibration amendments a1/a2 disclosed in the spec, the original symmetric rule's outcome printed permanently): no crossing + extended census O(1) + margin slope -3.33 >= pre-wall loading-speed slope -4.44 - 1.0 (margin/v8 in 19..652, median 144: the margin falls SLOWER than the local loading speed, the benign direction) + q_N flat-positive => the margin decay is the spectral SIGNATURE of the known O(1) wall offset (kinematic eigenvalue-density shrinkage), not an approach to failure) + CENSUS_EXTENDED(minC - N_w on the 15 new anchors {0:1, 1:10, 2:2, 3:1, 4:1}, max +4 <= 8 -- the O(1) wall offset SURVIVES beyond the family cap; equivalence margin > 0 <=> minC >= N_w two-routes 57/57) + DETECTOR_LEDGER(EXTRAP_CALIBRATED: the 42-rung power fit predicts all 15 extension margins inside the 1-decade band, max dev 0.78 dec; no-mp/short-fit/admissibility-cap/posthoc-window mutants all flagged; scopes + fragment audit clean).  NO L* claim, NO RH claim."

_L2_STATUS = "L2 standard-framework census (25/25, the 30-pct lane): WHICH classical deterministic cancellation framework captures the level-2 block sums P_j?  Four sealed frames, each an exact per-rung bound on |sum P_j| = |R|, tested against the r272 flip condition delta' > 0.21 and the 7-exception certification.  VERDICT FRAMEWORK_PARTIAL: F1 DISCREPANCY/Erdos-Turan (Abel form D* x V) dies loudly (both factors grow, delta' -0.21, cert 0/7, paircorr FIRE); F2 VAN DER CORPUT at the frozen H = ceil(sqrt(m)) -- an EXACT world-blind THEOREM, no constant, no fit -- delivers delta' +0.31 > 0.21 and certifies 6/7 exceptions (all but kz15; ladder 38/42, demand max +0.04); F3 Salem-Zygmund sub-Gauss (measured constant C med 0.90 max 1.68, spread 1.86 STABLE, controls covered) delta' +0.38, 6/7; F4 MDS/Azuma against the sign-bin chain filtration (md_ratio 0.29, C4 max 0.68 stable, covered) delta' +0.41, 6/7 -- sealed key adjudicates F4, the honest edge is F2 (exact theorem at 6/7).  kz15 razor misses 0.07-0.12 dec on every surviving frame.  AUTOCORRELATION_PROFILE surprise: only 5/44 rungs (0 exceptions) have a positive root-scale cancellation deficit -- the block sequence REINFORCES on 39/44 worlds ((sum P)^2 > sum P^2): the c3 cancellation is already absorbed INTO the P_j magnitudes; the classical frames win by sqrt(m)-count economy, not by lag structure. L2 lemma candidate: the vdC form as a provable chain statement + the r269 kz15 exact-finite certificate.  mp dps-60 wards; r272 regressions exact.  NO RH claim."

_R288_STATUS = "L* destructive-coherence anatomy (31/31): WHAT about the exact prime-comb source makes the nu interference in the mu frame DESTRUCTIVE (the r285 z_v = -3.15 address)?  All E signs come from the mu-CD kernel K_n between nu positions; phases measured against the zeros of the mu-orthonormal P_{N_w} (Jacobi eigenvalues of the chain; Sturm counts, midpoint alternation and mp dps-60 edge wards exact).  VERDICT SIGN_ANATOMY(the CARRIER MAP: neg fraction 0.503 at N_w but band-structured -- fold distance 1-2 is 87 pct POSITIVE (in-phase), distance 3-4 is 80 pct NEGATIVE (antiphase); the top mode at N_w rides the in-phase adjacent pairs (band 1-2 = +0.68 of |T|, ARCH-ARCH +0.74) while the wall destructivity in the source frame at the crossing (X_v = -0.0517, C_off -0.1046) is carried by the ANTIPHASE next-nearest 3-4 fold ARCH-ARCH pairs (-0.105); controls at their crossing ALL constructive, C_off +0.30..+1.00 -- total sign reversal) + SAMPLING_BLIND(the pure sinc/CD sign formula reaches 0.728 plain / 0.878 |E|-weighted on MAIN, controls 0.81..0.83 -- below the 0.90 bar: only a PARTIAL sign mechanism) + SOURCE_SEPARATOR_NOT_FOUND(the sealed phase-dispersion candidate K_S1 = R = 0.352 is WORLD_BLIND (dead worlds bracket MAIN from both sides) and ensemble-GENERIC (pct 0.44/0.25); chain correlations sp(K_S1, assist_cross) only +0.12/+0.24 -- the plain dispersion does not predict the wall; honest negative) + DIFFERENT_OBJECTS(the r287 block object REINFORCES on w9 (canc = -4.5e-2 < 0, pipeline gates exact) while the kernel-side wall interference is destructive: the r287 magnitude absorption and the kernel destructivity are TWO coordinates) + the DOSE-CURVE FINDING: under conservation-faithful union jitter the phase field barely rotates (turn rate 0.24/dose == static prediction 0.23 -- the chain zeros CO-MOVE; quarter turn would need dose ~1.0) yet z_v flips destructive -> constructive already at dose 0.005 and the depth collapses to 0.71 N_w: the destructive coherence sits in hyper-fine alignments far below the median phase resolution -- the r276 metric firewall is sharper than any phase-rotation account, MEASURED (r276 P2_JIT quoted COMPARISON_ONLY, different channel).  Ensembles = the 28 r285 replicates with identical seeds (records re-gated: assist_cross med 1.89/3.99, MAIN pct 0.00/0.00); equidistribution statistics typed MEASUREMENT_ONLY (forbidden as proof premises).  NO L* claim, NO RH claim."

_R289_STATUS = "L* arch-kernel diophantine relevance test (42/42; strict reviewer contract: no proof attack, no Baker/Matveev application -- need documentation only): the four questions on the r288 hyper-fine coherence, decided by a sealed THREE-WAY twin adjudication.  VERDICT METRIC_ONLY: the RATIONAL twin (every tent center log n replaced by a small-denominator rational multiple of Delta via CF convergents, denominators <= 5.7e4, position cost <= 2.1e-9 = 1e-8 x local gap, weights/cells/on-node family conserved exactly) KEEPS the FULL signature identically -- minC 184, crossing 185, z_v -3.149, carrier band 3-4 -0.105, AA -0.056, mp dps-60 ward at the 1.7e-4 margin -- every exact log-linear-form relation destroyed with ZERO effect; the SHUFFLE twin (tent-split fractions permuted among the atoms, fraction distribution/cells/weights preserved exactly, assignment destroyed) loses totally (depth 0.21, z_v +4.4) but at effective metric dose 0.12 gap where the plain jitter ladder already loses (collapse between theta 1e-3 and 3e-3): fully metric-explained, Baker unnecessary on all tested scales.  Q1 (kernel-degree anatomy): X_v decomposes exactly over CD degrees (sum identity 8e-16); the destructivity is carried by the MID degree band 50-75 pct (-3.17 z-units all pairs / -2.48 carrier set), NOT the terminal degrees (-0.99); at dose 0.005 with degree and frame FIXED no kernel band moves O(1) (all <= 0.01 z) -- the O(1) mover is the CROSSING RELOCATION itself (pivot cascade 185 -> ~131), while the phase field stands still (dphi 0.0012; r288 anchor reproduced with identical seeds).  Q2 (sub-gap attribution): the tent-split fractions {log n / Delta} = {46 log2 n} are the ONLY sub-grid entry of the source into the weights (two-node completeness 4.2e-14; fraction DIFFERENCES underdetermine, global-shift mutant 0.80 LOUD; folded weights typed DOWNSTREAM_COMPLETE); gradient chain d/dfrac = Delta d/du exact (r278 Hellmann-Feynman re-gate 3.1e-07); comb channel at 0.005 (positions grid-frozen): terminal |dlg| 0.47 = O(1) while dphi 0.0067 -- 29x below the sealed 100x bar, honest SUBGAP_CARRIER NOT_SEALED (the r288 200x was a union-channel statement).  Q3 (linear forms, DOCUMENTATION_ONLY): Delta = log2/46 EXACT; 8/70 atoms on-node => 28 EXACT 2-power resonances {46 log2 2^(a-b)} = 0 inside the window; min nonzero |{46 log2(n/m)}| = 1.34e-6 at (181, 241); attribution tops at the small primes n = 3, 7, 25, 8, 29, 17; need exponent ~8.9 vs Baker/Matveev-class ~7e4 (~7900x weaker in the exponent) -- diophantine input only matters for the UNBOUNDED family; fraction detectors WORLD_BLIND.  The measured METRIC THRESHOLD: the coherence needs the fraction PROFILE at ~1e-3 of the local gap and is indifferent to its number-theoretic exactness.  NO L* claim, NO RH claim."

_R290_STATUS = "L* profile-space geometry of the working set (32/32; SPEC_SHA f953dd71): after r289 METRIC_ONLY the round maps WHERE in weight-profile space the destructive coherence survives.  Profile = the signed lag density on the fixed grid (linear in the tent-split fractions at fixed cells, r289 completeness); distances = theta_eq in the exact LAG coordinate with the ANALYTIC reference 0.5 sum m g / Delta = 125.75 (amendment a1, demanded by the sealed self-consistency gate: the density-L1 realization was ill-conditioned, ~37 pct seed noise).  VERDICT BASIN_GEOMETRY(onset map over five sealed interpolation paths MAIN -> SMOOTH/SCR/EPST/ENS-replicate/RIDGE at the NEAR = 0.90 bar: world-directed axes die at theta_eq 4e-5..2.3e-4 -- 5..50x BELOW the 1e-3..3e-3 jitter threshold: the rim is strongly ANISOTROPIC; all onsets GRADUAL (soft shoulder s 0.95 -> 0.77..0.89, no cliff; |kappa| <= 0.08 quasi-linear decay); ISOTROPY = TUBE (killfrac 1.00 at 3e-3 over 16 conservation-gated random directions, 0.38/0.62/1.00 at 5e-4/1e-3/2e-3 -- NEAR-radius ~5e-4..2e-3, consistent with the r285 ensemble pct 0.00); SMOOTH direction = a privileged killer axis (onset 1.36e-4) yet ORTHOGONAL to the first-order wall gradient (cos -3e-5): collective, not gradient-linear) + RIDGE_MAP(the r280 OPT axis rebuilt (theta_up 3.873e-5 anchor exact, endpoint minC 185): the crest LIFTS the wall (minC 185) at extension factors 1..8 = theta_eq 1.5e-4..1.2e-3 -- inside a tube where half the random directions already kill -- and only dies at 2.4e-3) + ALL_FUNCTIONALS_BLIND(the four sealed source-pure profile scalars -- antiphase 3-4 autocorrelation (r288 carrier), total variation, perturbative r278 gradient alignment, mid-band 50-75 pct deviation energy (r289 flip carrier) -- ALL fail to predict the survival depth over the 187-point test set (best sp +0.263, bars 0.5/0.8) and none separates the worlds non-trivially: the coherence functional, if it exists as a closed profile scalar, is NOT among these classes; the working set stays implicitly characterized).  Channel identity gated bitwise on MAIN (union, minC, crossing, z_v, rho vs crossing_from_B); r288/r289 records + r289 jitter-ladder anchor (exact seeds) + r285 control flips 25/21/27/25 + 28 replicates (ENS_SIGN as exact density sign-flips, rep-0 identity gated) reproduced.  NO L* claim, NO RH claim."

_R291_STATUS = 'L* ridge anatomy (30/30; SPEC_SHA bb512c17): the r280/r290 raising axis as an analysis object.  VERDICT RIDGE_SECTION(k_min = 9 at the matched dose; MATCHED-DOSE BUDGET_THRESHOLD PERFECT: one scalar m_star in (1.280, 1.291] separates lift from no-lift over all 18 factor-1 subset cases -- the side-selected first-order wall budget -sum c_j dose with a ~1.3 second-order resistance factor is the near-threshold story; carriers at factor 1: PRIME + HEAD + XIPOS; top-9 atoms n = 2, 3, 5, 13, 11, 4, 29, 7, 89 (small-prime head-heavy, NOT one-atom); GLOBAL rule violated by exactly one overdrive retraction, TOP6 at factor 8 (margin 9.09, minC back to 184, fully alive)) + RIDGE_NO_FIXPOINT(PLATEAU(185): the step-controlled r280 recipe saturates at minC 185 from iterate 1 and NEVER reaches 186; the axis decoheres after iterate 3 (cos +0.92/+1.00/+0.43 then negative) and every full step is rejected thereafter -- a one-degree crest plateau, no stationary better world) + LIFT_MAIN_SPECIFIC(the same iteration from EPSTEIN terminates at step 0 (theta_up 5.03e10 = first-order flat wall c_25 = -2.0e-11); the matched-dose ladder 7.75e-5..0.1 never lifts EPSTEIN) + NONLOCAL_BLIND(best F5S ridge projection sp +0.471 < 0.6 and below the trivial size baseline |sp(F0)| 0.881; F6 rank-2 lethality Gram (train/test split-sealed) +0.137, F7 antiphase deviation pair correlation -0.463 -- the third functional class family stays blind, the working set remains implicitly characterized) + SMOOTH_ANATOMY(SMOOTH_COLLECTIVE_2ND_ORDER: direction PR 112 of dim 368 delocalized; margin curvature d2 -23.3 Richardson-stable vs same-length random -0.98, ratio 23.7 >= 10 -- a genuine quadratic wall valley, NOT an r259-style higher-order resummation gap).  w9/metric/control/ridge anchors gated (REF 125.75, flips 25/21/27/25, theta_up 3.873e-5, r290 map facs 1/4/8 + death at 16 reproduced); amendments a1-a3 disclosed (r290 pinned metric quadruple, 1e-9 conservation-gate headroom, matched-dose/global reporting split).  NO L* claim, NO RH claim.'

_R292_STATUS = "L* curvature two-form round (36/36; SPEC_SHA 050821ff): after the r291 budget-threshold anatomy the round measures whether the working set is characterized by a curvature TWO-FORM.  The margin Hesse form at MAIN over a sealed 29-direction set (3 world axes + ridge + the 9 ridge atoms + 16 fresh conservation-gated random directions, theta_eq-normalized; 3-step diagonals + the full 406-pair polarization matrix, expansion-identity crosscheck 15/15).  VERDICT HESSIAN_SPECTRUM(HESSIAN_LOW_RANK(1) in the sealed L2 spectral metric: ONE eigendirection carries 92.5 pct of sum|lam| = 0.452, a pure on-support DENS combination -- NOT the SMOOTH axis (|cos| 0.07) and NOT the ridge (RIDGE_CURV_FLAT, L2-rank 28/29: the lift axis sits in the flattest curvature sector, exactly the r291 expectation); all 29 diagonals negative -- a genuine multi-directional quadratic valley) + EPSTEIN_CONTRAST(EPST_CURV_FLAT ratio 5.4e-15: the first-order flat EPSTEIN wall is ALSO second-order structureless in the log-h channel -- MAIN's wall has curvature structure, EPSTEIN's has none) + THRESHOLD_NOT_EXPLAINED(the measured log|h_184| quadratic term along the 18 r291 matched-dose sub-directions is NEGATIVE on all 18 (q_R median -0.911; one-sided cubic-triple solve, linear identity vs the analytic budget max rel dev 0.0033) -- a second-order ASSIST, not resistance: the corrected budgets b - q/2 still separate PERFECTLY but the bracket midpoint moves to 1.761, AWAY from the naive flip level 1 -- the r291 ~1.3x factor is a near-flip nonlinearity invisible to low-order jets; m* stays an emergent threshold object) + RETRACTION_HIGHER_ORDER(B2(8) = 27.6 and B3(8) = 150.9 both over-predict lift vs the measured TOP6 factor-8 retraction: not explained at second or third order, honest) + CURVATURE_BLIND(the split-sealed two-form functionals F10 curvature energy / F11 budget+curvature / F12 Rayleigh reach sp +0.884/+0.672/+0.186 on the disjoint 94-point test corpus vs the size baseline |sp(F0)| = 0.907 -- the FOURTH honest negative under the sealed rule, but the closest miss of the whole functional program: 0.023 below the baseline where the r290/r291 classes lost by 0.4-0.6; AUC(dead) 0.097).  w9/metric/control/ridge anchors + the full r291 threshold regression (18 margins max abs dev 5e-4, bracket (1.280, 1.291], TOP6 retraction, top-9 atoms) + SMOOTH d2 -23.31 (ratio 23.7) + EPSTEIN first-order flatness re-gated; amendments a1-a4 disclosed (Hessian step-ladder halved at the measured margin-valid scale, one-sided estimator for the side-selected budget kink, fractional-dose conservation-gate pad 1e-7, L2-unit-consistent ridge display).  NO L* claim, NO RH claim."

_R293_STATUS = "L* metric reconciliation round (43/43; SPEC_SHA 33c44cc6): the r292 near-miss adjudicated.  Leg A seals THREE metrics BEFORE any functional evaluation (M1 theta_eq, M2 L2 on the signed density, M3 = the curvature metric sqrt(x |H_tr| x) with the spectral-absolute-value rule; distortion factors max/min M2-M1 11.0 / M3-M1 59.5 / M3-M2 123.5, largest-deviation sector DENS = the r292 top-eigenaxis sector; positive-mass share of H_tr 0.0052).  VERDICT METRIC_RECONCILED(F10 in M2): THE HAMMER LANDED -- evaluated against the size baseline OF ITS OWN METRIC, the split-sealed curvature energy F10 beats a baseline for the first time in four rounds (sp +0.884 vs sp(F0_M2) -0.860, home margin +0.024) AND holds FUNCTIONAL_BEYOND_BASELINE (partial spearman sp(F10, s | F0_M2) = +0.423 >= 0.4 on the 94-point test split, sign-stable +0.826 on the train-side replica) -- the FIRST predictive closed functional of the working set; PROMOTION CANDIDATE flagged for the next consolidation wave (nothing promoted here) + MIX_CONTROL(MIX_IS_CAUSE: the r292 mixed margin -0.023 vs the reconciled home margin +0.024 -- the theta_eq-vs-L2 metric mix WAS the cause of the r292 near-miss) + the m* resolution: all 8 selected flips (4 closest lifters + 4 closest non-lifters by |b - 1.2853|) are SIMPLE ZEROS of h_184 (alpha = 1.000 +- 0.003 on both sides; 225 gated wall evaluations, conservation exact, death watch clean, every flip lands at minC 185 with h_185 < 0) but MSTAR_NO_LAW: the per-direction thresholds m*_dir = b x f* range 1.13..1.39 (XIPOS 1.1295 .. TOP6 1.3931) -- NO fixed safety distance to the singularity in budget, theta_eq, L2 or curvature units (gap spreads 3.0/3.0/3.8/10.2 vs bar 1.25): the r291 'one scalar m*' bracket was a fixed-dose artifact, the ~1.3 factor is TYPICAL, not universal + RETRACTION_SIDE(the TOP6 factor-8 retraction IS a second simple h_184 zero at f_ret = 7.107, lift window (1.226, 7.107) -- on the side the monotone budget prognosis (8.07 and rising) is structurally blind to; h_185 sign +1 there) + BASELINE_ANATOMY(F0_DIRECTIONAL: median within-band |sp| 0.370 > 0.3, dead-only sp -0.890 / alive-only -0.005 -- the baseline is NOT a pure radius artifact, the r290-r292 absolute-sp bars stand as drawn).  w9/metric/control/ridge anchors + the full r291 threshold regression + the COMPLETE r292 spectroscopy (trace share 0.925, lam_top -0.418, SMOOTH |cos| 0.07, ridge rank 28/29) and contest (F10 +0.884 / F0_M1 -0.907 / AUC 0.097 on the bit-identical 94-point split) re-gated; amendment a1 disclosed (the Leg-C death watch redrawn from the raw h_185 sign to the real criterion minC < 184 -- diagnostic display fix, no bar moved).  NO L* claim, NO RH claim."

_R294_STATUS = "L* F10 stability round (42/42; SPEC_SHA 88c6fd1e): the r293 promotion candidate stress-tested BEFORE consolidation wave 8, exactly as the r293 worker demanded.  Leg A (the sealed promotion bar, fixed before any corpus was built): FIVE fresh seed-disjoint 147-point corpora (seeds 294000+1000k; same conservation discipline and r293 test-family mixture; tags pairwise disjoint and disjoint from the r292 training tags, gated) against the UNCHANGED hash-sealed r293 H_tr (sha256 gate 3447ed198a56, re-training mutant CAUGHT): the |sp| WIN over the home L2 baseline replicates on EVERY corpus (5/5: F10 +0.675..+0.787 vs F0_M2 -0.660..-0.714, margins +0.003..+0.115) BUT the partial-correlation channel does NOT replicate at its r293 magnitude (partials +0.248..+0.555, median 0.299; only 2/5 clear the sealed +0.3 bar, two miss by <= 0.003) => sealed verdict F10_FRAGILE (win 5/5, part 2/5, full 2/5) -- promotion precondition NOT met, recommendation NO for wave 8; the r293 partial +0.423 was the TOP of the fresh distribution, not the center.  Leg B robustness: TRAIN_ROBUST (leave-3-out jackknife sigma 0.0101 <= 0.012) + RANK2_CARRIES (the top DENS eigenaxis alone reaches +0.855 and misses the 0.860 baseline by 0.005; two whitened eigendirections suffice; top-axis DENS coefficient-mass share 0.989) + STEP_STABLE (half FD ladder sp(F10_mid, F10_half) = 0.9995; the double ladder MARGIN_INVALID with 4 non-finite margins = the r292-a1 NaN boundary reproduced and excluded by the sealed rule).  Leg C window transport (w7/w11, complete window-own chains: walls NREF_w 59/63, extraction dose caps 0.157/0.068 after the measured w11 full-ridge collapse, ladders unhalved, 67-point window corpora, window-local survival readout, 0 censored): KNIFE-EDGE -- w11 beats its own L2 baseline by +0.002, w7 loses by 0.009 => WINDOW_PARTIAL(w11); the solid transport finding is MECHANISM CONSTANCY: the top curvature eigenaxis is a DENS combination on w7/w9/w11 alike (coefficient-mass shares 0.795/0.989/0.834).  Leg D: L2_NOT_DENS -- the DENS-neutralized baseline (six DENS eigenaxes rescaled to the typical curvature per unit L2, scales c_j 2.2..75.6) is WEAKER than plain L2 (sp -0.828 vs -0.860) and the partial vs it RISES to +0.586: the F10 information is not the DENS-sector length distortion; why L2 is the right coordinate stays open.  Anchors complete (w9; the FULL r293 contest bit-identical incl. tag SHA a76cc578 and partials +0.423/+0.826; r292 spectroscopy 0.925/-0.418); must-fails m1 hash-seal / m2 seed-reuse (146 tags) / m3 conservation / m4 window category error (dimension-blocked, disclosed control) all CAUGHT.  NO promotion, NO L* claim, NO RH claim."

_R295_STATUS = "L* F10 sp-hardening round (44/44; SPEC_SHA 4d7d8095): the two r294 fronts executed -- harden the partial-free |sp| statement on a LARGE ensemble, and dissect the partial channel by direction family.  Leg A (the sealed THREE-CLAUSE hardening bar, fixed before any corpus was built and untouchable): TWENTY fresh seed-disjoint 147-point corpora (seeds 300000+1000k; the FULL forbidden-seed set of r292/r293/r294 -- 192 seeds -- enumerated and gated, zero overlap; r294 corpus protocol verbatim) against the UNCHANGED H_tr (hash gated == the published r294 seal 3447ed198a56): F10 wins the |sp| contest on 14/20 corpora (margins -0.079..+0.125, median +0.028, six losses incl. -0.041/-0.057/-0.079) => sealed verdict F10_SP_MAJORITY -- the HARDENED bar (>= 18/20 AND margin median >= +0.02 AND no loss beyond 0.02) fails on TWO clauses: the r294 5/5 was itself TOP-OF-DISTRIBUTION, exactly the way the r293 partial was; the partial-free ranking statement is a DOCUMENTED REGULARITY (19/25 combined fresh-corpus census), NOT a theorem candidate; promotion recommendation NO for wave 8; Leg C (candidate statement) VOID by the sealed rule.  PARTIAL20 median +0.279 IQR [+0.182, +0.401] range [-0.133, +0.566] (the r293 +0.423 confirmed top-of-distribution; PARTIAL_STD = the full 147-point-mix partial sealed as the standardized statistic of future rounds).  Leg B partial anatomy (pooled 2921 points, MAIN once): NO family reaches the STRONG 0.4 bar -- PATH +0.156 / WORLD +0.245 / FRAC +0.067 NULL / DENS +0.104 (regime split disclosed; ATOM/RIDGE structurally absent by the unchanged split seal): the beyond-distance information is diffuse, strongest on structured deviations toward dead worlds, near-null on random FRAC rays; the r293-composition-matched subsample (87 points per corpus, the PATH-heavy r293 mix 0.46/0.29/0.17) lifts the median partial to +0.346 -- the sealed gain clause (+0.067 >= +0.05) PASSES but the level clause fails by 0.004 (0.346 < 0.35) => R293_LUCK under the sealed adjudication: composition explains PART of the r293/r294 partial discrepancy, the rest was n=94 sampling noise.  Leg D bycatch: L2_VIA_CONSERVATION formally (the eta_0 conservation projection leaves the contest exactly invariant, |sp| shifts 0.0000) but near-tautological (projected-out eta_0 share median 2.9e-30, max 3.2e-3): 'why L2' stays open in substance, weakness disclosed in the verdict.  Anchors complete (w9; the FULL r293 contest bit-identical incl. tag SHA a76cc578, partial +0.423/+0.826 and the family census 1/5/40/25/15/8; r292 spectroscopy 0.925/-0.418; the COMPLETE r294 core rebuilt bitnear: 5-corpus table, jackknife sigma 0.0101, rank sp +0.855/+0.863/+0.884, DENS share 0.989); must-fails m1 hash / m2 seed collision (26 overlaps flagged on the r294 base) / m3 bar sharpness (every clause has teeth, both directions) / m4 conservation all CAUGHT.  NO promotion, NO L* claim, NO RH claim."

_R296_STATUS = "L* DENS-identity round (39/39; SPEC_SHA ffb413c8): the reviewer hard fork on the r292 curvature axis -- is the 92.5-pct DENS top eigendirection e_top mathematically identifiable and coupled to the L* extremizer?  ONE sealed candidate library, hash-frozen BEFORE any cosine: T1 = grad_delta lambda_max(E_184) via Hellmann-Feynman from the extremal eigenvector phi* (envelope theorem dlam = int p*^2 d(dnu) - lam int p*^2 d(dmu), fold-map pullback with the s_l envelope; FD ward worst rel dev 3.8e-7 on 6 sealed coordinates at the re-measured margin-valid step; lam2/lam3 gradients as excluded controls), the von Mangoldt family (all/prime-only/prime-power-only), the local-metric family (fold-gap deviations / 1-over-log / r289 fraction profile), the kernel diagonals at N_w (Christoffel / E-diag / CD-diag), the 4-dim moment-gradient subspace T5, 8 nulls + a decoy; PRIMARY metric = the in-span cosine cos(e_top, P_D T) on the 29-direction span (e_top exists only there), raw/lag/capture honesty columns per candidate; bars 0.80/0.40 sealed.  VERDICT DENS_WORLD_BLIND(T5_MOMENTS 0.825 downgraded; WORLD-DOWNGRADE) + T1_COUPLING(cos_span +0.394 -- THE COUPLING NUMBER: below the sealed 0.40 bar, above the 0.331 noise max, below noise+pad 0.431; the lam3 control couples HARDER (-0.574): what coupling exists is a near-crossing eigenvalue-BAND property, not top-specific) + COS_TABLE(every arithmetic T2/T3/T4 member <= 0.38, all within reach of the noise reference; T3_LOCALGAP degenerate -- the gap profile is constant on the 367/368-occupied fold grid, disclosed) + WORLD_CONTRAST(the T5 moment-subspace root share HOLDS at EPST 0.571 / SCR 0.603 >= 0.40 => construction-adjacent low-moment smoothness, the sealed downgrade fires PART -> BLIND; DISCLOSED TENSION: both control values sit BELOW their own null-quadruple maxima 0.758/0.691 -- the absolute 0.40 collapse bar is near-unreachable for a 4-dim subspace in a 12-dim span, the downgrade is conservative and the sealed rule stands unmoved; EPST's margin channel at its OWN wall degree 25 carries Richardson-STABLE curvature structure -- r292's flatness was the log-h channel, a different object, disclosed) + WINDOW_TRANSPORT(w7 0.930 / w9 0.825 / w11 0.942 for T5 -- but the w7/w11 values sit AT their window null maxima 0.923/0.734 (12-dim spans); T1_w 0.179/0.314 below window noise; window DENS shares 0.795/0.834 == r294, class rule 3/3 hold, no window downgrade).  THE FORK ANSWER: the DENS top axis is NOT identifiable with any sealed arithmetic density, and at the sealed bar NOT the L*-gradient direction (T1 misses PART by 0.006, honest); the only above-chance structure is its low-moment smoothness (0.825 vs ~0.4 chance at w9); raw-grid overlaps <= 0.08 with every candidate at captures <= 0.36 -- the axis is dominated by structure no on-support arithmetic profile in the library represents.  Anchors complete (w9 + metric + flips + ridge/top9 + lambda_max(E_184) = 0.99983248 dev 4.1e-9 + r283 top-5 + EPSTEIN first-order flatness + the FULL r292 spectroscopy rebuilt verbatim incl. the DENS05/04/03/01 loadings); amendments a1-a4 disclosed pre-freeze (evenness gate conditioning at libm precision, FD coordinate/step re-anchoring at the measured margin-valid scale, the r294 0.989 share correctly re-typed as the 18-direction training-form record, hardening repetition pinned to the ACTUAL best member incl. the subspace rule); must-fails m1 wrong-normalization / m2 unprojected / m3 seal-break / m4 decoy all CAUGHT/HELD.  NO L* claim, NO RH claim."

_R297_STATUS = "L2 vdC chain-provenance round (27/27; SPEC_SHA e42a76eb; reviewer routing after the r296 DENS close: resources to the L2 front): can the vdC INPUT (the P-variance / autocorrelation scaling) be derived source-purely from the chain, turning the v964-S0 van der Corput theorem plus input into delta' > 0.21 as a theorem on the generic half?  VERDICT TARGET_INEQUALITY_FROZEN(the weakest sufficient input bound, all constants live: sigma := slope(S_F/M^2) <= sigma* = 2 x (slope(eps_c2/M) - 0.21) - slope(pref) = -0.516; measured sigma -0.714 -- the truth clears the target with 0.198 sigma margin = 0.099 in delta': what is missing is not room, it is provenance; exact slope additivity 2.6e-16; m2 pad-sharpness: sigma* composes to delta' = 0.21 EXACTLY, the pad-dropped mutant to 0) + ROUTE_TABLE(three source-pure provenance routes, every chain step EXACT on 47 worlds incl. broken-arithmetic controls -- world-blind by the same algebra: B1 PAIR-IDENTITY (A(0) <= gapmax x eps_c2 via the r271 alternation identity) delta' -0.112, break: gapmax/M slope -0.535 vs needed <= -1.178 -- the MAX pair gap falls far slower than the mean amplitude (B_amp -0.81): the Jacobi near-balance is a MEAN statement, its extreme is not chain-controlled; B2 NODE-DENSITY (A(0) <= 2 maxM sumM; r273 root-scale trace) delta' -0.097, the boundary POSITIONS are near-equidistributed with FALLING discrepancy (D_rank med 0.024, slope -0.42 -- provable terrain) but the MASS imbalance n x maxM/sumM GROWS (slope +0.244): node density supplies the count, not the decay; B3 FEJER/PARSEVAL (per-run Cauchy-Schwarz + the symbolic sum rule sum w x^2 pi_k^2 = h_{k+1} + a_k^2 h_k + b_k^2 h_{k-1} PROVED exact in Fractions) delta' -0.097 -- with THE SURPRISE OF THE ROUND: the chain is EXACTLY ORTHOGONAL w.r.t. the WINDOW measure (cross/diag devs 0.000 across kz15/kz20/mains/EPST/SCR) while BORDER and WDIFF break (cross med ~1.1): the Parseval attachment point EXISTS; what remains is the measure transfer (the vdC input lives on the BORDER union) and the signed-to-absolute gap (med 4.3, slope +0.217, mixed-sign weights)) + CHAIN_PROVENANCE_PARTIAL(B3_PARSEVAL, gap: slope miss 0.307) + KZ15_NOTE(structurally outside for all three routes, misses 0.57/0.86/0.91 dec -- the r270 exact-finite certificate remains the permanent closure).  Leg C void; the common obstruction named: every universal majorization of A(0) through MAGNITUDES pays a growing max/mean imbalance factor; the SIGNED structure carrying sigma = -0.71 is what magnitude chains cannot see.  The two assets left for the lane: the exact window-measure orthogonality (a bordered-form transfer statement window -> border union is the concrete named object) and the falling rank-coordinate equidistribution.  Anchor regression exact (r287 delta'(F2) +0.309, cert 6/7 + 38/42, kz15 miss 0.096 dec; r272 sp +0.67 / gamma_true +0.453 / sl_c2 +0.196; r270 reserve 0.0268; v964-S0 vdC set re-computed green in Fractions); amendments a1 (index bugfix after smoke crash) / a2 (reporting precision) disclosed, no bar moved; mp dps-60 deep wards.  NO L2 lemma claim, NO RH claim."

_R298_STATUS = "L2 window-border transfer round (28/28; SPEC_SHA 05e831be; the r297 asset-1 attachment point executed): does an exact, sign-preserving transfer statement exist from the WINDOW measure (where the chain is exactly orthogonal and the Parseval sum rule is proved) to the BORDER form (where the vdC input S_F lives)?  VERDICT DECOMPOSITION_EXACT(S_F = B(omega,omega) + B(Delta,omega+beta) with Delta = border (-) window and B the frozen positional Fejer block kernel -- S_F re-derived as the Fejer quadratic form B(beta,beta) of the border measure, kernel reproduction worst 7.5e-16 against the sealed r287 pipeline; the identity exact to 8.8e-16 of the decomposition scale on 47 worlds incl. broken-arithmetic controls; per-block attribution t_j = PDelta_j phi_j sums exactly; the LINEAR corollary MEASURED: the window image of the drive functional VANISHES at float precision (2.2e-12 main / 6.4e-12 deep) -- the drive is 100 pct transfer at the linear level) + TRANSFER_DOMINANT(the sealed three-way adjudication: NOT subdominant, NOT the r259/r261 near-cancellation -- the window main term is EMPTY: MAIN/M^2 med -3.94 dec, two decades BELOW S_F med -2.00 dec, slope -1.386 collapsing; T > 0 on 42/42, share med 0.99; the omega-Delta cross term negligible (med -1.4e-4) => S_F ~ B(PDelta,PDelta): THE VDC INPUT IS THE FEJER ENERGY OF THE DIFFERENCE MEASURE ITSELF; sigma = -0.714 is the decay exponent of this ONE named source-pure quantity; the in-T sign cancellation sum|t_j|/|T| med 1.69 grows at +0.207 -- the r297 imbalance factor localized INSIDE T; m2 mustfail: the magnitude version |Delta| breaks the identity by 1.7e+0 of scale LOUD -- the r297 magnitude error certified; controls EPST/SCR same class, WORLD-BLIND) + IMBALANCE_NOTE(BROAD_BASED: dropping the top-1/top-3 mass runs lowers the imbalance slope +0.244 -> +0.188/+0.197, not below the 0.5x bar -- no outlier carrier; the connection stat sp(block mass, |t_j|) med +0.69: the mass-imbalance carriers ARE the T carriers -- the two r297 gap halves are ONE structure) + Leg C VOID (no candidate theorem; the precisely narrowed wave-9 gap object: prove the Fejer-energy decay of Delta).  Anchors bit-near (sigma -0.714 / sigma* -0.516 / sl_c2 +0.196 / sl_pref +0.489; D_rank med 0.024 sl -0.42; imb sl +0.244; WINDOW orth HOLDS 5.1e-7, BORDER cross med 1.101; sum rule re-proved exact in Fractions, wrong-b^2 mutant CAUGHT); amendments a1 (toy-constant hand-arithmetic fix after smoke) / a2 (reporting text) disclosed, no bar moved; mp dps-60 deep wards.  NO L2 lemma claim, NO RH claim."

_R300_STATUS = "L2 diag-target round (31/31; SPEC_SHA 55218b5d; the r299 frozen rest pair executed: can the diagonal decay sl_D <= sigma* = -0.516 -- margin 0.055, THIN -- be derived source-purely, and the ratio flatness with it?): the exact participation anatomy of the ONE positive quantity D = sum PDelta_j^2 plus three derivation routes and the exact kernel-envelope ratio majorant.  VERDICT DIAG_ANATOMY(n_eff = L1^2/D med 37.41 slope +0.963; the EXACT decomposition sl_D = 2 x sl_L1 - sl_neff = 2 x (+0.196) - (+0.963) = -0.571 (dev 6.7e-16): the entire diagonal decay is the participation growth net of a mildly growing L1 mass -- 'many small instead of few large' is now a measured identity; second coordinates sl_mx -0.542, sl_fill -0.225; class shares A/Bv/cross 0.64/0.16/+0.20; D band shares 0.27/0.26/0.46 -- D itself is BROADBAND: the r299 LOWPASS character of B is pure kernel weighting) + ROUTE_TABLE(B1 chain-norm CS on the |dw| measure valid 47/47 but composed -0.384 FAILS sigma* by 0.133 AND the |dw| identity census BREAKS (cross med 0.932) -- the proven sum rule does NOT attach to the difference measure, the chain-norm shortcut is closed honestly; B2 max x L1 valid 47/47 composed -0.346 FAILS by 0.170 -- the max IS atom-controlled (a single c-value distance, mx/maxatom med 1.07): the irreducible loss is the fill decay -0.225, invisible to every max x mass bound; B3 NEFF: needed slope(n_eff) >= 2 sl_L1 - sigma* = +0.908 (disclosed: an exact reparametrization of DIAG_TARGET), measured +0.963 MET margin 0.055; D_rank(Delta) slope -0.117 falling, sp(D_rank, n_eff) -0.81 -- the equidistribution bridge is real but correlational) + DIAG_SPLIT(B3(NEFF_TARGET)) + RATIO_BOUNDED_STRUCTURAL(the exact kernel envelope F_H <= min(H, 1/(H sin^2)) gives B/D <= R_env med 1.61 with slope -0.122 <= +0.05 on 47/47 -- the ratio half of the r299 rest pair is settled at the structural level; the lobe-width heuristic REFUTED, sp(B/D, q50) -0.47) + REST_FROZEN(ONE inequality remains: NEFF_TARGET slope(n_eff) >= +0.908, measured +0.963, margin 0.055) + WORLD_SENSITIVE(the FILL class separates MAIN (FILL_LOW) from both dead controls (FILL_HIGH) -- the second world-separating class of the lane after the r299 O-sign class).  Anchors bit-near (r297 sigma/sigma*, r298 full set, r299 full set incl. sl_D -0.571 / ratio 1.29/-0.168 / full overlap 42/42); m1 participation mutant 8/3 exact, m2 swapped max/L1 toy 2 exact + real 0.84/0.72 of D LOUD, m3 wrong-h sum rule CAUGHT in Fractions, m4a/m4b scope-flagged; amendment a1 (reconstruction scale floor on the degenerate SMOOTH Delta == 0) disclosed, no bar moved; mp dps-60 deep wards.  NO L2 lemma claim, NO RH claim."

_R301_STATUS = "L2 neff-target round (32/32; SPEC_SHA 6f8cc404; the r300 frozen rest executed: why does n_eff = L1^2/D grow ~N, and can slope(n_eff) >= 2 sl_L1 - sigma* = +0.908 -- measured +0.963, margin 0.055 THIN -- be derived source-purely?): the Renyi order anatomy, the jackknife margin stability and three derivation routes.  VERDICT NEFF_ANATOMY(the order family N_2/N_3/N_4/N_inf med 37.41/27.88/24.04/15.47 grows at nearly ONE exponent +0.963/+0.926/+0.894/+0.738 (tail N_4/N_2 slope only -0.069): ECHT anti-concentration, not a broadening tail; positional IQR flat -0.005 -- the count grows, the footprint stays) + STABILITY(jackknife min/max +0.936/+0.979, 0/42 below NEED -- NO MARGIN_FRAGILE, the margin survives every single-rung removal; half-ladder lo/hi +0.982/+0.802: the deep half flattens below +0.908, the honest anti-extrapolation flag) + ROUTE_TABLE(B1 THE STRUCTURAL FINDING: n_eff = n_act/(1 + CV^2) EXACTLY (module-own Fractions re-proof, ward <= 6.8e-16 on 47 worlds) with the PERFECT count link n_act == m on 42/42 -- the effective-carrier count IS the constructive level-2 block count, growing at +1.002; B2 the exact weighted-discrepancy interval chain (|w_j/W - lambda_j| <= 2 delta_w, a star-discrepancy theorem, wards 0.0 on 47 worlds) FAILS by 0.834: the |dc|-WEIGHTED discrepancy (med 0.167, slope -0.017) is 10x the raw D_rank and near-flat -- the r300 bridge correlation lived on the raw positions, the max-based equidistribution shortcut is closed honestly; B3 the difference profile decorrelates at lag ONE (l_loc med 1.0, slope -0.036), n_eff_atom slope +0.942 coupled sp +0.96 -- but the LOCAL sum-rule census BREAKS 24/24 (cross med 0.93): the |dw| non-orthogonality is scale-free, no local energy identity attaches) + NEFF_SPLIT(B1(UNIF_TARGET) | B3(MIX_TARGET)) + REST_FROZEN(strictly smaller than r300's: the ONE remaining inequality is UNIF_TARGET -- slope(1 + CV^2) <= sl_nact - NEED = +0.094, measured +0.039, margin 0.055: the growth statement became a BOUNDEDNESS statement about the block profile of |PDelta|; NEFF_DERIVED missed only the sealed CV_FLAT clause, sl_cv2p +0.039 > 0) + WORLD_SENSITIVE(the UNIF class does NOT separate: MAIN CV^2 1.03 bracketed by EPST 0.72 and SCR 2.76 -- unlike the O-sign and FILL classes, quasi-uniformity is not where the arithmetic lives; r300 PART/FILL classes reproduced EXACTLY).  Anchors bit-near (r297 sigma/sigma*, r298 full set, r299 full set, r300 full set incl. n_eff 37.41 / sl_neff +0.963 / adds 6.7e-16); m1a/m1b Fractions identity mutants break 9 and 7/24 exact, m2a halved discrepancy prefactor 0.3 exact, m2b transposed aggregation matrix 37/55 exact + colsum flag, m3 overlapping sub-ladders CAUGHT, m4 synthetic one-block collapse to n_eff == 1 exact, m5a/m5b scope-flagged; smoke-stage identifier rename disclosed, no bar moved; mp dps-60 deep wards.  NO L2 lemma claim, NO RH claim.; consumed by v967 (wave 10, embedded byte-exact, smoke stage)"

_R302_STATUS = "L2 unif-target round (30/30; SPEC_SHA 36df9424; the r301 frozen rest executed: why is the block profile of |PDelta| quasi-uniform -- slope(1 + CV^2) <= sl_nact - NEED = +0.094, measured +0.039, margin 0.055 -- and can it be derived source-purely?): the exact CV^2 class decomposition, the distribution-stationarity census, the exact per-half depth attribution and three derivation routes.  VERDICT CV2_ANATOMY(between shares SIGN/BAND/PAT med 0.012/0.216/0.111 -- no sealed carrier rule fires: the mild rise is STRUCTURELESS, within-class) + PROFILE_STATIONARY(normalized pooled KS(G1,G3) 0.043 / KS(G2,G3) 0.021 / KS(G1,G2) 0.026, all far under the sealed 0.125 over a 3x depth range -- the |PDelta| block profile has ONE stable normalized distribution: CV^2 boundedness is its second moment, and the r299 pointwise-cconv negative was the wrong convergence notion, convergence in DISTRIBUTION) + DEPTH(the r301 half-ladder flag resolved: the deep-half flattening of n_eff (+0.982 -> +0.802) is carried ENTIRELY by the CV^2 head (+0.228), NOT the count (+1.030, additivity exact per half); the sealed finite-size model fires TRANSIENT_1_OVER_N (m2 group medians 2.07/2.02/2.00 fall along an exact 1/N law onto A = 1.973, held-out dev 0.002 <= 0.05, approached from ABOVE) -- DEPTH_HALF_MISS stands as an honest flag, DEPTH_CAVEAT does NOT fire; jackknife of sl_cv2p 0/42 above UNIF_NEED) + ROUTE_TABLE(B1 THE STRUCTURAL FINDING: the exact coherence identity 1 + CV^2 = n_act chi/(surv^2 n_eff_atom) with chi = D/Q the in-block coherence factor (module-own Fractions re-proof, ward <= 6.0e-16 on 47 worlds, per-block Cauchy-Schwarz cap chi <= kmax slack 0.0) -- chi med 0.630 DESTRUCTIVE and falling (-0.060), surv flat (-0.020), the lag-1 anti-correlation rho_1 = -0.22 is its atom-level mechanism; composed slope +0.039 == sl_cv2p exact (5.4e-16); B2 the local-pattern route closed honestly on BOTH clauses (PAT12 within-share 0.685 > 0.5 -- the node pattern does NOT determine |PDelta|; k-profile ks 0.181 > 0.125); B3 the localized-perturbation gain g/N med 1.079 vs the sealed damping bar 1.0 -- the recursion responds AT the polynomial-degree rate, miss 0.079 honest) + UNIF_DERIVED(B1 + A2: the sealed clause fired -- chi non-growing + stationary profile + 1/N transient + count link 42/42) + CANDIDATE_THEOREM(wave-10 candidate NOT promoted, every hypothesis typed MEASURED: the ONE remaining measured growth statement is ATOM_TARGET -- slope(n_eff_atom) >= sl_nact + sl_chi - 2 sl_surv - UNIF_NEED = +0.888, measured +0.942, margin 0.055: the block-level uniformity became an atom-level anti-concentration statement about the c-difference profile itself) + WORLD_NOTE(the coherence identity holds on 47 worlds by algebra -- world-blind as the charter demands; UNIF non-separation reproduced, MAIN 1.03 bracketed by EPST 0.72 / SCR 2.76).  Anchors bit-near (r297 sigma/sigma*, r300 participation set, r301 full set incl. orders/count-link/halves/jackknife); m1 wrong-class-weight 3/4 exact, m2 unnormalized-KS mustfail re-anchored to the relative form (a1 disclosed, original 0.335 < 0.5 printed permanently; CAUGHT 8x), m3 lag-0 double counting 5/7 exact, m4 two-point synthetic CV^2 729/196 LOUD on every route, m5a/m5b scope-flagged; smoke-stage m3 fix + a2 label reporting disclosed, no adjudication bar moved; mp dps-60 deep wards.  NO L2 lemma claim, NO RH claim.; consumed by v967 (wave 10, embedded byte-exact, smoke stage)"

_R303_STATUS = "L2 atom-target round with HARD REGRESS AUDIT (26/26; SPEC_SHA 375e9f2b; the two mandatory questions of round 303 adjudicated under sealed rules): (Q1) is the thrice-identical 0.055 margin of the r300/r301/r302 cascade ONE algebraic number?  VERDICT MARGIN_CHAIN + REGRESS_CONFIRMED: the four target margins m_D/m_NEFF/m_UNIF/m_ATOM all equal +0.0547 with invariance devs <= 9e-16 (module-own margin_chain, re-proved exact in Fractions on a rational slope set: all four margins 3/50 EXACT; the frozen halves_slope estimator is exactly log-additive over the r300/r301/r302 product identities) -- the cascade is an exact DICTIONARY around ONE measured core S = sigma* - sl_D = +0.0547 (the r299 DIAG margin), NOT three reductions; the r297 sigma level is NOT the core: sigma* - sigma = +0.1976 = S + ratio surplus (sl_D - sigma) +0.1429 EXACT, and the charter's 1/2-conversion conjecture is REFUTED ((sigma*-sigma)/2 = +0.0988 != +0.0547; the factor 2 lives only in NEED = 2 sl_L1 - sigma* and cancels in every margin difference); further reduction rounds may NOT be counted as progress on the inequality.  (Q2) does the slack follow from the lag-1 mixing structure -- the FIRST non-tautological mechanism test of the lane: sealed synthetic dc rearrangements per rung (exact marginal = bitwise multiset permutation, rho_1 steered by seeded greedy swaps to matched / zero / flipped targets, 8 replicates, 1008 builds, convergence 1008/1008, fixed points <= 0.024).  VERDICT MIXING_INSUFFICIENT with a MONOTONE mechanism ladder: chi real 0.630 -> (a) matched 0.764 -> (b) zero 1.029 -> (c) flipped 1.342; end margin +0.055 -> +0.057 -> +0.032 -> -0.044 -- flipping the mixing sign KILLS the target inequality, removing it shrinks the margin, matching rho_1 reproduces the slopes and the end margin (+0.057 vs +0.055) but NOT the destructive-coherence level (chi miss 0.134 > 0.05): the within-block structure beyond lag 1 carries the rest; BYCATCH THEOREM-GRADE FACT: n_eff_atom = L1a^2/Q is a pure MARGINAL functional (invariance 1.0e-15 on all 1008 builds) -- ATOM_TARGET's growth lives in the value distribution alone, the mixing lives in chi/CV^2/D: the r302 A2 stationarity and the mixing are complementary, not redundant.  RHO_SIGN census honest: S1 < 0 on 41/42 only (one rung positive -- the adjacent-pi orthogonality connection cannot be a per-rung theorem in raw form; Fractions-exact negative certificates on the two smallest rungs).  CONSEQUENCE (ii): the reduction cascade r297->r302 is CLOSED as a dictionary; the honest next object is the joint short-range law of the dc profile (rho_1..rho_k / within-block coherence), not another coordinate change.  Mustfails m1 wrong-1/2-factor breaks by exactly sl_L1 = 1/5 (Fractions) / m2 unmatched marginal LOUD 5.9e-1 / m3 seed collision CAUGHT / m4 swapped D/Q breaks 48/35 EXACT / m5a+m5b scope-flagged; NO amendment after freeze; mp dps-60 deep wards.  NO L2 lemma claim, NO RH claim.; consumed by v967 (wave 10, embedded byte-exact, smoke stage)"

_R304_STATUS = "L2 short-range-law round (37/37; SPEC_SHA 2cc5d23f; the reviewer contract executed under the r303 hard rule -- NO coordinate change, only the named mechanism object): (Q1) is the dc-profile anti-correlation genuinely short-range and summable?  MEASURED AND FROZEN: the centered lag profile rho_1..16 over all 42 rungs is a STABLE PERIOD-4 COMB (med -0.222/-0.140/+0.089/+0.130, then recurring strong lags at k = 4m and 4m+2 up to 16: +0.13..+0.16 at 4m, -0.14..-0.15 at 4m+2; halves-stable 3/3 sign agreement; world-specific -- EPST/SCR differ), NOT a decaying tail; the sealed short-range rule (|med rho_j| <= 0.05 for all j > k0 <= 8) finds NO k0 => VERDICT LONGRANGE_STRUCTURE + LAW_LONGRANGE (the reviewer stop case fires on the LETTER of the rule); the reviewer condition SPLITS: the net covariance NC(16) = 1 + 2 sum med rho_k = 0.712 < 1 HOLDS (net-negative) while summability-with-small-tail FAILS (SUM(16) = 1.563, no decay); the zero-sum tautology 1 + 2 sum_all rho_k = 0 for mean-removed finite profiles disclosed and gated live (F1, 3.5e-15) -- the truncated NC is the content.  THE STRUCTURAL COUNTERPOINT (F2, exact): chi = D/Q = 1 + 2 sum_k T_k/Q with T_k the within-block lag-k cross sums (recomposition ward 4.7e-16 on every non-degenerate world) shows the chi-RELEVANT structure IS short-range: real shares c_1 -0.345 / c_2 -0.139 / c_3 +0.068 / c_4 +0.028, sum(5..16) -0.026, tail(>16) 0.000 -- the long-range comb is INVISIBLE to chi; the r303 chi-level miss 0.134 is attributed to k <= 3 (dominated by lag 2: +0.181).  (Q2) does lag-k matching close the r303 gap?  LAGK_MATCH (rho_1..8 steered jointly per rung, module-own swap_arrange_k, 252 builds, convergence 252/252 at tol 0.02, marginal bitwise, natom invariance 8.5e-16): chi 0.652 vs real 0.630 -- the destructive LEVEL now reproduces (|d| 0.022 <= 0.05; the r303 rho_1-only family missed by 0.134) BUT the slopes break (sl_cv2p miss 0.028, sl_D miss 0.027, both > 0.02; margin +0.082): the exact COMPLEMENT of r303 (a) -- the mechanism is TWO-SCALE (level from the short lags, slopes from the rung-wise lag-1 trend); no lag-matching family reproduces both within the sealed bands.  The graduated reviewer ladder (absolute rho_1 targets +0.2/0.0/-0.1/-0.2/-0.3, 840 builds, 840/840 converged): chi 1.361/1.019/0.889/0.765/0.632 MONOTONE in rho_1; margin +0.036/+0.074/+0.019/+0.045/+0.059 NOT monotone -- HONEST NEGATIVE: absolute-level targets do not kill the inequality; the r303 kill (-0.044) came from per-rung FLIPPED targets: the margin responds to the rung-wise trend of rho_1, not its level.  SIGN_PATTERN: ladder med signs of S_1..4 = -/-/+/+ with Fractions-EXACT certificates on both smallest rungs (kz18, kz23: 4/4 match).  Leg 0 anchors bit-near (r297-r302 full set; margin chain +0.0547 x 4 invariance <= 8.7e-16; the r303 three-family ladder REPLICATED with the same seeds: chi 0.764/1.029/1.342, margins +0.057/+0.032/-0.044; RHO_SIGN 41/42, alt 0.523).  CONSEQUENCE (iii): the lane's global-profile mixing route is CLOSED (documented stop): L2 <=> anti-concentration of an explicit block field with LONG-RANGE (period-4 comb) structure; what survives as honest structure: the exact within-block short-range decomposition of chi, the two-scale mechanism split, NC < 1 and the exact sign pattern -- typed MEASURED, no proof target claimed on the closed route.  Mustfails m1 factor-1 net covariance 1/2 != 0 EXACT + chi decomposition 6/7 vs 5/7 EXACT / m2 unmatched marginal LOUD 5.9e-1 / m3 uncentered field 3/4 EXACT + real 0.028 / m4 family smuggle dev 0.332 > 0.02 CAUGHT / m5a+m5b scope-flagged; pre-record G32 print fix disclosed, NO amendment after freeze; seeds 2100 collision-free; mp dps-60 deep wards.  NO L2 lemma claim, NO RH claim.; consumed by v967 (wave 10, embedded byte-exact, smoke stage)"

_R299_STATUS = "L2 Fejer decay round (32/32; SPEC_SHA f432e944; the r298 gap object executed: can the decay of the Delta Fejer energy B(PDelta,PDelta) -- sigma <= sigma* = -0.516 -- be derived source-purely from position equidistribution plus a controllable signed-mass statement?): the module-own EXACT spectral representation B = (1/L) sum_k F_H(theta_k) |Dhat(k)|^2 (frozen L rule, recomposition + Parseval + magnitude-majorization wards <= 7.1e-16 on 47 worlds) plus three source-pure derivation routes.  VERDICT SPECTRUM_MAP(LOWPASS: band shares LOW/MID/EDGE med 0.93/0.04/0.02, band slopes -0.758/-0.014/-0.979, q50 med 0.19 main-lobe units -- sigma is a LOW-FREQUENCY phenomenon, not broadband) + ROUTE_TABLE(B1 ERDOS-TURAN/ABEL valid 47/47 but composed slope +1.948 FAILS sigma* by 2.464 (position/kernel factor +1.504 GROWS, Abel mass term V^2 +0.444 grows -- the r297 magnitude wall recurs at the frequency level; MASS_TARGET frozen and NOT MET by 2.464); B2 PAIR/DIAGONAL valid 47/47: sl_D -0.571 MEETS sigma* -0.516 (margin 0.055), O < 0 on 13/42 only (the PDelta pair field REINFORCES in the majority, the r287 border pattern carries over), ratio B/D med 1.29 slope -0.168 FALLING; B3 OVERLAP valid 47/47 with THE STRUCTURAL FINDING: FULL-SUPPORT overlap on 42/42 rungs -- the border and window unions share EVERY position, Delta_fresh == 0 identically: Delta IS a pure c-value difference measure on one shared node set; cconv med 0.86 slope +0.045, the relative c-difference does NOT fall -- the decay is aggregate cancellation of a non-converging difference profile, class OV_DOM) + DECAY_SPLIT(B2 DIAG_TARGET: prove sl_D <= sigma* (a magnitude-density bound on ONE positive quantity, measured margin 0.055) + ratio flatness (measured -0.168); B3 CVALUE_TARGET: prove the c-value convergence rate) + CANDIDATE_THEOREM(conditional, wave-9 candidate NOT promoted: Leg-B bound + frozen rest hypotheses => sigma <= sigma* => r297 target => v964-S0 vdC => delta' > 0.21 generic half) + CANC_LOCUS(the in-T cancellation gap G = B_abs - B is ENTIRELY low-band (1.01/-0.01/-0.00) and NOT carried by the mass blocks (sp(Mb, c_j) med -0.01, coincidence 3/42 vs the +0.69 energy coupling)) + WORLD_SENSITIVE(disclosed: the O-sign class separates MAIN (reinforcing, O_POS) from EPST/SCR (cancelling, O_NEG) -- the first world-separating sign class of the L2 lane).  Anchors bit-near (sigma -0.714 / sigma* -0.516 / D_rank 0.024 sl -0.42 / imb +0.244; r298 share 0.99, sl_MAIN -1.386, canc 1.69 sl +0.207, overlap w9 367; identity/kernel/attribution 8.8e-16/7.5e-16/5.4e-16); Fractions ET section exact (Abel resummation dev 0, position-bound equality witness, dropped-terminal 0 < 2 and halved-prefactor 2 > 1 mutants CAUGHT, v964-S0 T4 anchor); m1 wrong-Fejer-weight 2/3 exact, m2 magnitude-pair 1.3e+0/4.6e+0 LOUD, m3 tolerant-overlap-split CAUGHT with silent recomposition; amendments a1 (O-sign class into the world-blindness conjunction per the sealed spec text) / a2 (full-support overlap count print) disclosed, no bar moved; mp dps-60 deep wards.  NO L2 lemma claim, NO RH claim."

_R306_STATUS = "L2 Renyi-3 round (27/27; SPEC_SHA 3bb365e1; reviewer plan section 5 -- the sharpest new fiber attack, explicitly sanctioned as NEW information after the r304 documented stop: POINTWISE instead of slope, cubic moment instead of the closed max/discrepancy routes; contract PRIME.L2.ATOM.RENYI3.01, experiments-side, NO ledger row): prove or refute sum_j q_j^3 <= C (log m)^A / m^2 pointwise on every rung of the atomic |PDelta| profile, A in {0,1,2} preregistered, C frozen on the first 5 rungs.  VERDICT RENYI3_GO(C = 1.069, A = 2): the inequality HOLDS on ALL 57 rungs (42 core + the 15 r286 extension anchors N_w 942..1218) with GROWING reserve (trend -0.322 on both frozen estimators; deepest 15 rungs at reserve 2.6..5.7); A = 0 fails (the third shape moment M_3 = m^2 S_3 rises +0.153 before decaying -- the disclosed algebraic prior), A = 1 fails on exactly the two near-critical rungs kz53/kz67 (1.4-1.9 percent under C_2): polylog-SQUARED is the honest sharp form.  Via the exact Renyi/Hill chain (two Lagrange witnesses, Fractions dev 0; Lean-ready statement printed) the GO fixes the NEW POINTWISE PROOF TARGET of the fiber: n_eff = N_2 >= N_3 >= m/(1.034 log m) -- asymptotically above every demanded power m^0.888.  Anatomy: the cubic mass is a NARROW-BAND functional (97 percent within a factor 4 of q_max) carried ever more BROADLY (top-8 share 0.780 falling -0.268; strict triples 92 percent of the cube) -- no heavy-few obstruction; minimal sufficient statement named SHAPE3_TARGET (M_3 polylog-bounded, the third-moment sibling of the r302-proven M_2 stationarity 1.973, M_2 med here 2.03).  WORLD_SENSITIVE with the right sign: EPSTEIN holds, SCRAMBLE BREAKS the bound 1.67x -- the law is recursion structure, not block combinatorics.  Leg 0 anchors bit-near (Renyi family 37.41/27.88/24.04/15.47 at +0.963/+0.926/+0.894/+0.738; core slack S = +0.0547 dev 6.7e-16; count link 42/42 AND 15/15 on the extension).  Honest negatives: C_2 frozen on the shallowest rung, kz53/kz67 margin 1.4-1.9 percent; the GO fixes a proof TARGET, it proves nothing beyond the 57 rungs; the atom-level (dc) sibling stays unformalized.  Must-fails m1-m4 CAUGHT/LOUD (1/64 exact; prereg breach set ward + C_full/C_cal 2.10; L1^3 factor; one-atom 212x) + scope mutants flagged.  Deterministic run1/run2; NO RH CLAIM in either direction; consumed by v968 (wave 11, embedded byte-exact, smoke stage)"

_R307_STATUS = "L* fixed-head kill round (28/28; SPEC_SHA ec2bb008; reviewer adjudication section 4 solution C -- the cheapest kill test of proof architecture C (fixed head + contractive tail), deciding whether the r284 near-one-atom anatomy is descriptive or proof-relevant; contract PRIME.LSTAR.FIXED_HEAD.01, experiments-side, NO ledger row): in the mu-orthonormal frame E = E_H + E_T with H = the first r FLAT ARCHIMEDEAN EDGE atoms below the first prime (u < log 2); two sealed window-fixed orderings (ORDER_FOLD ascending fold; ORDER_W9X = the frozen w9 extremal rank permutation, regated exactly); lambda_max(E_T(r)) measured for r = 0..16 on 77 windows (42 core + 15 r286 anchors + 20 NEW deeper windows N_w 1218..1650, S up to 3299).  VERDICT HEAD_GROWS: NO fixed r <= 8 exists (best candidate min reserve 6.54e-6 at r = 8 FOLD, four decades below the sealed 1e-2 bar; a fixed r_min <= 16 exists on 15/77 windows only, all N_w <= 434; sp(reserve(8), S) = -0.93 with TS slope -3.15 == the r286 margin-law class: removing a FIXED head changes the prefactor, not the decay class) BUT the FULL edge (13..96 atoms, GROWING with the window, sp(edge, N_w) = +0.98) restores a macroscopic reserve on 77/77 (min 2.19e-3, median amplification 3.1e4): the near-one-atom anatomy is REAL and proof-relevant AS A COORDINATE, refuted AS A FIXED DIMENSION -- architecture C in its fixed-head form is dead; the honest successor object is the window-dependent head of size ~ log2/(2 Delta_w), whose own reserve decays slowly (TS -0.87, disclosed).  SCHUR_LEDGER exact: det(I - E) = det(I - E_T) det(S_H) worst 3.4e-8 on 57/57 at r = 8 (rational-EXACT on the JF9 toy: S_H = 5079569347/6986834571); border-augmented identity det(I - E_T) det(S_aug) = det(I - E) D_N worst 2.1e-6 with D_N = (5/7)(1 - q_N) from the INDEPENDENT r244 bordered-chain route (q_N < 1 on 57/57, min gap 0.0195, r286 reproduced; the Christoffel W-identity t^T (I - E)^{-1} t = sum rho gated 4.4e-13 on w9 and exactly on the toy); lambda_min(S_H(8)) > 0 on 57/57 BUT == the wall margin to 4 digits and anchor-MONOTONE (sp -0.92, TS -3.33) along the r286 margin decay -- the hoped-for simplification is EMPTY: S_H carries the whole problem; the border is inert (S_aug positive 57/57, degradation factor >= 0.873 -- BORDER_BREAKS_HEAD did not fire).  TWIN_METRIC_OK (r289 rational twin verbatim: fold sets == MAIN, r_min 3/3 == 3/3 both orders, max |dreserve| 1.0e-7 over r = 0..16 -- diophantine input irrelevant on the head coordinate).  CONTROL_COLLECTIVE (EPST lambda 2191 -> 2189 at r = 16, SCR 1.4e7 -> 3.2e4: no small head rescues a dead world -- the architecture reads MAIN-specific structure).  Tail anatomy: after edge removal the tail rides its own near-single-atom Christoffel profile on the bulk atoms (maxdiag 0.862..0.992, gain <= 1.089); Gershgorin/Schur-test row sums 1.5..3.6, < 1 on 0/57 -- the classical bounds stay dead on the tail.  Structural exactness: reserve monotone in r + interlacing lambda_max(E_T) <= lambda_max(E) EXACT on 77 x both orders (sign safety inherits from the r286 mp-verified margins); mp spots 3/3 at dps 30 (worst rel dev 2.8e-14).  Must-fails: m1 target-inverse ordering AST-FLAGGED (ordering constructors clean of eigensolvers) / m2 Woodbury wrong sign breaks the det identity 0.163 LOUD / m3 monomial frame hides the toy wall (4.36 -> 0.85, rel 0.805) LOUD / m4 r = 0 degenerate exact.  Four pre-spec sizing scratches disclosed in the spec (they saw the likely outcome; bars sized before freeze, no rule tuned after any evaluation).  Deterministic run1/run2; NO L* claim, NO RH CLAIM in either direction; consumed by v968 (wave 11, embedded byte-exact, smoke stage)"

_R308_STATUS = "L* block-Green DISCOVERY round (30/30; SPEC_SHA d5147850; reviewer adjudication par.4 solution A, plan round '305': automated search for a matrix-valued discrete Green identity on fourfold fold blocks -- the route replacing the failed scalar profile functionals (r290-r295), the cancellation-losing absolute-value majorants (r297/r299/r300) and the irrelevant diophantine relations (r289 METRIC_ONLY); experiments-side, NO ledger row): does the augmented form Q(p,t) = int p^2 d(mu-nu) + 2t u(p) + B t^2 (v958 bordered object, B = S_{N-2} + 5/7) admit an EXACT decomposition Q = sum_r Delta_r^T G_r Delta_r + (5/7)t^2 with sealed fourfold-block coordinates (first/second differences, the r288 antiphase layer_i - layer_{i+2} coordinates, local gross-mass mean, border t) and source-pure blocks G_r?  VERDICT IDENTITY_EXISTS + BLOCK_INDEFINITE_MAIN + FEAS_DIAG world separation + CONTROL_LEDGER + R282_DEMARCATION: the identity EXISTS everywhere (exact Fractions with entrywise reconstruction on TOY4-hand/SM0-SM3/MINI16, dof 15..219; f64 residuals 4.0e-14..8.1e-14 on w9, the r289 rational twin, EPST/SCR/SMOOTH and all 57 rungs at DEG_A = 8; kernel exclusion 100 percent) -- existence is the disclosed EASY half (linear in G, hugely underdetermined); the sealed eigenvalue-free minimum-Frobenius selection is indefinite EVERYWHERE (even on the positive-class SM0: the SELECTION, not the world, produces indefiniteness; w9 -1.3e-02 rel, 199/364 negative blocks) => BLOCK_INDEFINITE_MAIN on the letter of the sealed rule; THE ROUND'S CENTRAL MEASUREMENT (sealed diagnosis clause, amendment a1 calibration-stage: uniform two-stage Dykstra schedule 200/2000 for every world): the SDP-like feasibility census SEPARATES THE WORLDS -- w9@DEG_A and twin@DEG_A converge to GENUINELY ALL-PSD block families (min eig rel +6.6e-16 / +2.0e-17 at affine residual 3.7e-14, 200 steps) while EPST/SCR stall at -0.45/-0.49 after 2000 steps: on the restricted DEG_A subspace, where every world's target is positive definite and NOTHING is theorem-forced, MAIN and its diophantine-trivialized twin admit an all-psd fourfold-block Green decomposition and the two hard controls (at this schedule) do not -- the first block-level world discriminator of the L* lane, ONE-SIDED (non-convergence proves no infeasibility) and diagnosis-grade only; SMOOTH@A marginal (-1.3e-05, its crossing 28 == DEG_B boundary role); w9/twin@DEG_B = 28 stay OPEN (-4.2e-05 after 2000).  DEG_B disclosure: control crossings 26/22/28, so control negativity at DEG_B is FORCED by theorem given existence -- gated live on EPST/SCR (held), never sold as discovery.  Anchors: w9 367/263/104, crossing 185; r288 carrier pin z_v -3.149 destructive; r282 pin N_2(MAINLIKE) exact; B_w9 = 8.368649 with 182 nonnegative rho_k.  R282 demarcation: the Kasteleyn/SOS refutations concern the FULL signed configuration class (every cell, h_n itself); this identity lives on the RESTRICTED subspace deg p <= 8/28 << N_w plus border -- different language, no contradiction.  Must-fails: cholesky-of-target FLAGGED, omitted remainder breaks by EXACTLY 5/7 at (t,t), posthoc family FLAGGED, partial coverage exact defect 4.5e+01 LOUD.  Honest: a feasibility census fixes a TARGET (find a sealed CONSTRUCTIVE psd G rule on MAIN, extend to DEG_B, prove the world split), it proves neither L* nor any cofinal statement.  Deterministic run1/run2; NO RH CLAIM in either direction; consumed by v968 (wave 11, embedded byte-exact, smoke stage)"

_R309_STATUS = "L* paired-cone PILOT round (33/33; SPEC_SHA f8d99877; reviewer adjudication par.4 solution B / plan round 306 -- the third full proof architecture: positive and negative source processed PAIRWISE as rank-one updates of the augmented form, base and border together; contract PRIME.LSTAR.PAIRED_CONE.01, experiments-side, NO ledger row): the r308/v958 augmented form A = [[H, u], [u^T, B]] is rebuilt EXACTLY as an ordered signed rank-one stream (sealed base = 5/7 t-floor + reserved/unmatched mu atoms; one signed budget step |S_{N-2}|; nu atoms greedy-matched to their nearest free mu partner; every border atom split exactly into a +/- pair), and at every negative step the reviewer reserve r = 1 - L_bb + L_ab^2/(1 + L_aa) is measured -- the Sherman-Morrison cancellation as an EXPLICIT POSITIVE SQUARE, with the det dictionary det(M1) = (1 + L_aa) det(M), det(M2) = r det(M1) exact on the toys (hand toy PC4: reserve 7/9, dets 9/14 -> 27/14 -> 3/2) and f64-warded on w9 (recon 4.8e-16, mp ward dps 60 dev 3.0e-13).  VERDICT DECOMP_EXACT + MAIN_RESERVE_POSITIVE + CONE_PARTIAL + WORLD_DISCRIMINATOR(all three sealed stats WORLD_BLIND) + MUSTFAIL_LEDGER: MAIN survives the ordered chain at BOTH caps with macroscopic margin (w9 min r +0.9011 at DEG_A = 8 / +0.3232 at DEG_B = 28; the r289 rational twin bit-near; 57/57 ladder rungs all-positive, min +0.7357@kz97) and w9's two tightest pair reserves sit EXACTLY on the r284 extremal folds {2, 4}; the reviewer's decisive cone question SPLITS: the sealed mass-ratio cone r >= 1 - c_1 v/w_a (c_1 = 8.5436 frozen on the 5 shallowest rungs) is SOUND -- ZERO violations on all 20714 pointwise test steps (52 rungs + w9 + twin) -- but certifies only 56 percent of the steps (bounds > 0 on 11643/20714), while the sharper local shield cone r >= 1 - c_2 v/SM2 overshoots exactly where it is positive (54/54 positive bounds violated, every violated step has true r >= +0.901): the reserve is PREDICTABLE in the sound direction from one source quantity but no sealed local phi is yet a certificate -- the named successor object is a phi with a base-depth term (the DEG_A base carries most of the reserve; not retro-fitted); the info test is SILENT by the sealed letter (pooled tightest decile best |sp| 0.50 at pair distance; on w9 alone the mass-ratio/shield features reach -0.61/-0.58, reported) and the reviewer-B4 O(log N) fold-bit register is measured IRRELEVANT at this cap (max |sp| 0.32).  Honest world census: at DEG_A the ordered PAIR reserve is WORLD-BLIND (controls' pair phases all positive, EPST +0.5712 / SCR +0.7796 / SMOOTH +0.9627); the ENTIRE EPST/SCR indefiniteness at DEG_A is the negative BUDGET step (S_{N-2} = -3.992 / -5.237 past the flip -- the r243 budget coordinate reads the wall, amendment-a1 disclosure: the design-time 'every DEG_A target is PD' sentence was corrected, spec-prose + reporting only, no bar or rule moved); DEG_B control breaks forced by theorem, disclosed (EPST/SCR budget step; SMOOTH boundary -7.08 at pair step 4).  Anti-circularity at r275 sharpness: the o1/o2 + forced-tail-sum pin re-run exact; target-inverse cone mutant (Minv_true) and posthoc-pairing mutant (r_true) FLAGGED by the AST scope audit; wrong-sign Sherman-Morrison caught by the exact det dictionary (dev 8/9 EXACT); order-break mutant changes the reserve by EXACTLY 4/9 while the reconstruction stays exact (the r-sequence, not the sum, is the ordered object).  Honest: a passing cone fixes a TARGET for proof work; this round proves neither L* nor any cofinal statement.  Deterministic run1/run2; NO RH CLAIM in either direction; consumed by v968 (wave 11, embedded byte-exact, smoke stage)"

_R311_STATUS = "L* block-Green NONTRIVIALITY round (34/34; SPEC_SHA fac7a8df; the reviewer's MANDATORY round before any constructive G-search; contract PRIME.LSTAR.BLOCKGREEN.NONTRIVIALITY.01, experiments-side, NO ledger row): is the r308 psd block family stronger and more local than Q >= 0, or a sparse-Gram/chordal restatement whose SDP separates the worlds by force?  VERDICT STRICT_SOURCE_CONE (with the binding SEPARATION MECHANISM clause) + SPAN_LEDGER + SIGN_EQUIVALENCE_LEDGER + R308_DEMARCATION -- BOTH reviewer scenarios are true at once, on different layers: (1) THE r308 WORLD SEPARATION IS FULLY THE BUDGET SIGN -- the sparsity graph of the Delta coordinates (band-3 atom cliques + universal border vertex t) is CHORDAL on every size incl. the full w9 S=367 graph (MCS + Tarjan-Yannakakis, maximal cliques == the 364 sealed blocks; C4/C5 control non-chordal), each block spans its full local R^5, so on EVALUATION space the image cone is exactly the chordal pattern-psd cone (Agler et al./Grone et al. -- restatement there); the control targets A_sys are INDEFINITE at DEG_A through the one negative budget diagonal (EPST -3.992 / SCR -5.237, the r309 amendment quantified: lambda_min rel -0.86/-1.00 vs w9/twin +8.4e-4 PD), so the r308 Dykstra separation was theorem-forced; the trivial common dual Y_t = e_t e_t^T (ONE rational entry -- the r243 budget positivity S_{N-2} >= 0, the KNOWN wall reading) closes the r308 one-sided stalls two-sidedly, predicts the two NEW blind scramble worlds 2/2 (seeds 2/3, S = -6.20/-6.15, both stall), and the DECISIVE ablations seal it: EPST/SCR with (t,t) <- |S| become PD and the Dykstra CONVERGES (+2.5e-16 / -5.1e-11), w9 with the EPST budget transplanted dies (-0.449) -- flip ONE scalar and the separation inverts; the exact SCHUR-TAIL LEMMA (s(A_sys@cap d) = S_{N-2} - sum_{k<=d} F_k^2/h_k, exact on SM0/SM1) shows the wall sign enters the target mechanically as the budget rho-tail, and decides all four r308 small-model OPEN stalls as trivially target-indefinite (even the positive-class SM0: s = -(rho_3..rho_8) = -7.08e-2 EXACT); SMOOTH's marginality is a BOUNDARY target (S = +3.48, lambda_min rel +2e-16).  BUT (2) THE CONE IS GENUINELY STRICT on the compressed subspace -- the round's second discovery: 0/16 sealed random in-span PD samples decompose (w9A 0/6, MINI16 0/6, SM1 0/4; staged Dykstra 200/2000/20000), every MINI16/SM1 stall carries a VALID polished-numeric Farkas dual (worst compression -1e-16..-1.4e-8 at <Y,Q> = -1), and MINI#0 + SM1#0 carry FULLY EXACT in-span rational certificates (exact Chebyshev basis, den 1e6, d = 1e-6, <Y,Q> < 0 exact rational; MINI16 = the real-w9 miniature with EXACT full span rank 55/55) -- C_w is STRICTLY smaller than span cap S_+: block membership is a genuinely stronger property than positivity, which MAIN + twin HAVE at DEG_A and generic positive forms LACK; the ambient span layer is also exact (codims TOY4 0 / SM0 11 / SM1 11 / SM2 2 / SM3 1 / MINI16 0 / w9A 4 / w9B 236, exact annihilators with all compressions zero, SM1 rank-one ambient counter-model with exact rational Farkas).  Protocol audits: r308 feas_diag AST-clean (objective-free, no target spectra), start ablation clean (zero start converges, random start -8.5e-9 at the amended 1e-8 bar; EPST stalls from every start).  w9@DEG_B stays honestly OPEN (PD-thin +2.4e-6, -1.78e-5 after 20000, dual diverges -- one-sided).  Consequence: the planned R312 rank-one-conductivity round premised on the r308 world discriminator is DEAD in that form (the discriminator was the wall sign); the re-scoped GO object is cone MEMBERSHIP (why the bordered-Hankel world targets sit inside the strict cone while generic PD forms do not).  Must-fails: wrong-sign Farkas rejected exact / wall-z FLAGGED / wrong-graph + C4/C5 caught / eigvecs_target sampler FLAGGED.  Three disclosed calibration amendments (dual divergence guard + sealed stage-3 consistency; START_NEAR bar; certificate bookkeeping + exact-Chebyshev-basis certification -- the four-way verdict selection rule never moved).  Deterministic run1/run2; NO L* claim, NO RH CLAIM in either direction; consumed by v968 (wave 11, embedded byte-exact, smoke stage)"

_R312_STATUS = "L* block-Green MEMBERSHIP round (32/32; SPEC_SHA 6c32f749; THE ONE construction round granted by the reviewer tree after r311's STRICT_SOURCE_CONE; contract PRIME.LSTAR.BLOCKGREEN.MEMBERSHIP.01, experiments-side, NO ledger row): does a source-pure constructive rule write down an explicit psd block family G_r = sum_l c_{r,l} v_l v_l^T (c >= 0) for the bordered-Hankel world targets, over the VORAB-sealed 22-generator library in the D basis (fold distances 1/2, antiphase 3/4, the r307 arch coordinates, border-vs-interior pairs D_a +- D6, symmetric/antisymmetric four-block modes), with an explicit coefficient formula c_{r,l}(source quantities)?  VERDICT COEFFICIENT_SIGN_WALL + MEMBERSHIP_ANATOMY + W9B_STATUS(OPEN) + R311_DEMARCATION -- the library SPAN carries the w9 identity exactly (unconstrained rel 2.4e-13; abstract span 15/21 exact) but the CONE does not: deterministic Lawson-Hanson NNLS stalls at rel 4.7e-4 (support 35, not capped) with a VERIFIED I-polished Farkas certificate (delta 1e-4, min normalized column product +9.7e-5, <y,q> = -0.9975), and the certificate mass localizes on the border/budget row (t-row fraction 0.58254 >= WALL_FRAC 0.5; at deg 28 even 0.725) -- the identity forces negative coefficients exactly on the wall coordinates: the sign question of the c IS the wall in the constructive language too; twin type-identical (METRIC_ONLY holds); EXACT GRADE on the miniatures: SM1 and MINI16 carry FULLY EXACT I-polished rational Farkas certificates (den 100, delta 1/10000, <Y,q> = -0.9988 / -1.0000 exact rational < 0, all 154/286 exact column products >= 0), and on SM1 even the library SPAN fails exactly ([M_lib | q] rank 41, EXISTS FALSE -- the sharpest small-model form of the wall); rung census 1/57 NNLS-feasible (worst rel 8e-4 at kz12) -- the sign wall is uniform across the ladder; the four sealed formula candidates PHI_A..PHI_D (masses, tent fractions, budget share, per-type table; calibrate-on-w9 / freeze / verify-on-57+twin) never fire (best calibration rel 7.7e-3 = PHI_D >> FORM_BAR 1e-8).  MEMBERSHIP ANATOMY (the round's discovery): the converged r308 Dykstra family on w9@A is high-rank {rank3: 17, rank4: 79, rank5: 268} and library-ALIGNED (top-eigvec alignment med 0.9973) with library-cone share med 0.976 -- 97.6 percent of the psd family lives in the sealed rank-one cone and the WHOLE membership obstruction sits in the remaining 2.4 percent of SDD-forbidden negative cross mixtures; dual anatomy of the 10 r311 sample duals: the generic-missing directions are low/mid-degree INTERIOR mass (|z_t| med 0.000 max 0.087, dominant degree indices 1..6), not border mass.  w9@deg28 stays honestly OPEN, two-sidedly hardened (staged anchor -1.784e-5 == the r311 record; library NNLS@B infeasible, same wall class; +80000 extension steps improve only -1.78e-5 -> -7.66e-6, a slow tail consistent with asymptotic feasibility, NOT decided; strengthened dual POCS at OMEGA 0.1 / 8000 iters runs but polishes invalid = no certificate).  Leg 0 bit-near on the identical rng stream (0/6 + 0/6 + 0/4 strictness, 10/10 polished duals, both exact certificates -0.9999977 / -0.9999930, Y_t values, chordality 364 cliques == blocks).  CONSEQUENCE (sealed reviewer tree): NOT a GO -- lane A closes as cone language with the membership mechanism NAMED (block-psd membership WITHOUT rank-one-SDD membership, obstructed at the budget/border sign), resources move to the fiber.  Must-fails: eigvecs_target library mutant FLAGGED / negative-c REJECTED loud / log-relation consumer EXPOSED by the MAIN-vs-twin comparison (0.886 >= 0.3) / exact recomposition break CAUGHT.  Two disclosed calibration amendments (a1 f64 I-polish certificate bookkeeping, the twin of the sealed exact I-polish; a2 the m3 harness probes the r289 rationalization variable u/Delta) -- the four-way adjudication tree never moved.  Deterministic run1/run2; NO L* claim, NO RH CLAIM in either direction; consumed by v968 (wave 11, embedded byte-exact, smoke stage)"

_R313_STATUS = "L2 Renyi-3 proof-form FORK round (32/32; SPEC_SHA 6505dd10; reviewer plan: triple incidence vs four-step Floquet -- two sealed proof architectures for the r306 pointwise cubic law sum q^3 <= C (log m)^A / m^2, executed side by side; contract PRIME.L2.RENYI3.PROOF_FORK.01, experiments-side, NO ledger row; reviewer instruction binding: the universal constant may be COARSER than the r306 sharp 1.069, and the theorem must not depend on von-Mangoldt arithmetic -- sought class RecursiveDifferenceProfile with MAIN + EPSTEIN in, SCRAMBLE out).  VERDICT BOTH_PARTIAL(named remnants): (R3A -- triple incidence) the raw atomic presentation is EXACT (P_j = sum of local c-differences; fold groups = beta/omega pairs at one position; block ward 3.9e-16 on 61 live worlds) and the cube splits EXACTLY into T1 diagonal / T2 two-equal / T3 near (shared fold) / T4 far (fully separated); med signed shares +0.38 / +2.15 / -1.14 / -0.42 -- the types CANCEL against each other, they are not independently small; proof hope (i) is EXACT (fold multiplicity UNIFORMLY 2 on 57/57 -- the class bound MULT_CAP = 2 is structural) and hope (iii) is REAL (TC_far med 0.069, the far triples telescope away 93 percent of their abs mass, falling); per-type prereg at A = 2 (first-5-frozen C): T2 HOLDS 0/57 (C 1.910, trend -0.400), T3 HOLDS 0/57 (C 1.053, trend -0.299), T1 FAILS 2/57 (kz55 2.54x, kz53 1.70x, trend -0.122) and T4 FAILS 1/57 (kz55 1.84x, trend -0.186) -- shallow-calibration artifacts of the r306 A <= 1 kind (trends fall), but the sealed rule stands: COMPOSITION REFUTED-IN-FORM; world census: SCRAMBLE violates exactly T2 (2.85x) but EPSTEIN also violates T2/T3 while its TOTAL holds the r306 bound -- MAIN-frozen per-type constants are NOT world-portable.  (R3B -- four-step Floquet) the local transfer step is REAL and EXACT (r_{k+1} = t_k + ap_k r_k + bp_k r_{k-1}, ward 8.0e-13/1.9e-12/6.6e-14, SMOOTH deg-skipped; det dictionary 4.4e-16) but the four-step monodromy on the split cubic mode space EXPANDS: med sr 1.0600, contracting rungs 0/57, product exponent +1.917 vs required -2.0 -- after r304 (single step) and this round (period-4 Floquet) NO transfer-operator contraction form of the cubic law remains open at the chain level.  CLASS finding (the sharpest honest negative): the sealed four-block comb property P4 excludes SCRAMBLE but ALSO EPSTEIN; every other sealed candidate (P1 fold, P2 recursion, P3 multiplicity, P5 mass, P6 assignment, CP4 telescope) is world-BLIND -- NO boolean property separates worlds the way the r306 total bound does: the separation lives in the SIZE of the type terms, the class needs a QUANTITATIVE membership functional; the m4 shuffle instance is rejected by P6 (293/300 atoms outside their interval) while mass conservation still holds on it.  Anchors bit-near r306/r304 (C-table 11.875/3.564/1.069, A2 0/57 + trend -0.322, comb -0.222/-0.140/+0.089/+0.130, NC(16) 0.829).  Must-fails m1 double-diagonal break == 3 S3 EXACT / m2 unsplit conservation mode toy break 1/2 EXACT (live clause vacuous, disclosed) / m3 magnitude telescope toy break 1 EXACT + real dev 9.4e-1 / m4 shuffle rejected / m5a+m5b scope-flagged; two disclosed smoke-stage code fixes (degenerate-guard extension, mass denominator), one reporting edit -- no bar or rule moved.  Deterministic run1/run2; NO RH CLAIM in either direction; consumed by v968 (wave 11, embedded byte-exact, smoke stage)"

_R314_STATUS = "L2 signed-cubic-flux round, part 1 -- EXACT ALGEBRA ONLY (29/29; SPEC_SHA 841b3196; reviewer contract after the r313 adjudication: not a better triple class but the SIGNED CUBE BOOKKEEPING -- with x_j = (PDelta)_j, sigma_j = sign(x_j), expand |x_j|^3 = sigma_j x_j^3 over the r313 fold genealogy into the signed cubic tensor C_{abc} = sum_j sigma_j x_{j,a} x_{j,b} x_{j,c} and decompose C_cube = DeltaF + C_collision + C_boundary; contract PRIME.L2.RENYI3.SIGNED_CUBIC_FLUX.01, experiments-side, NO ledger row; quantification deferred to R315/R316 by contract).  VERDICT SIGNED_CUBE_IDENTITY with a VANISHING boundary: sum_j |x_j|^3 = DeltaF + C_collision + C_boundary holds EXACTLY -- on w9 decided in EXACT Fractions on the real data (35 blocks / 150 fold groups: brute tensor enumeration == Newton aggregates == flux telescope with every dev the rational number 0; Fraction(float) is exact on dyadic f64 atoms) and to 4.5e-17 in f64 on all 57 rungs + mains + EPST/SCR (bar 1e-13).  THE DIVERGENCE FORM IS LOCAL: far triples enter as telescoping edge fluxes dF = 3 G (s1^2 - s2) consuming only the transported prefix state (s1, s2) along the position-ordered fold groups (telescope Newton-vs-path worst 4.1e-16); the opening-flux lemma F_1 = G^3 - 3 G G^2 + 2 G^3 == 0 is cubic algebra, so C_boundary == 0 -- the sealed language class closes WITHOUT remainder (disclosed: the razor acts before the genealogy; an unmasked presentation would re-introduce a boundary term).  THE COLLISION SPACE IS EXACTLY COUNTABLE by the banked r313 multiplicity-2 asset: config counts full/pair/far == n / 3n(n-1) / n(n-1)(n-2) (sum n^3) and atom-triple count 3 p1 p2 - 2 p3 == 8(n + 3n(n-1)) EXACT on 61/61 live worlds (w9: 150 full / 1656 pair / 1770 far configs, 14448 collision atom triples).  MEASURED PREVIEW for R315/R316 (labeled, no bound claim): med signed shares DeltaF -0.4226 (falling -0.422), C_pair +0.5980, C_full +0.8537 (rising +0.213) -- the two collision terms together (+1.45) carry the cube against the negative far flux; flux cancellation FC med 0.629 falling -0.141; world table w9 -0.452/+0.823/+0.629 FC 0.617 / w13 -0.541/+0.430/+1.111 / EPSTEIN -2.695/-2.652/+6.347 FC 0.101 / SCRAMBLE -0.171/+0.856/+0.315 FC 0.693; SCRAMBLE localization census names the FULL-collision column (dev 0.63 vs MAIN med) -- but EPSTEIN sits far from MAIN in every share while holding the r306 total bound: single shares do NOT separate worlds, Phi_3 must combine the divergence-form STRUCTURE (sealing discipline: R315 defines Phi_3 from structure, NOT from these record tables).  Anchors bit-near r313 (type shares +0.3808/+2.1542/-1.1442/-0.4226, TC_far 0.069/-0.050, fold multiplicity == 2 on 57/57, block ward 3.9e-16).  Must-fails: m1 sigma dropped, toy break 2 EXACT + w9 dev 1.8e0 LOUD / m2 double-counted fold group, breaks 2 and 98 EXACT + w9 4.6e-1 CAUGHT / m3 unordered support, per-edge break (0, 108, 108) EXACT with the total permutation-blind disclosed + w9 edgewise 7.7e-1 LOUD / m4 mult-blind count, break 56 EXACT (8 vs 64) + w9 8.75e-1 CAUGHT / m5a+m5b scope-flagged.  ONE disclosed reporting-only edit after smoke (m3 print format); no bar, band, rule or verdict rule moved.  Deterministic run1/run2; nothing here bounds anything -- part 1 is bookkeeping only; NO RH CLAIM in either direction; consumed by v968 (wave 11, embedded byte-exact, smoke stage)"

_R315_STATUS = "L2 signed-cubic-flux round, part 2 -- THE Phi_3 MEMBERSHIP FUNCTIONAL (29/29; SPEC_SHA 92d35a3a; reviewer contract part 2: define Phi_3 IN ADVANCE from the r314 divergence-form STRUCTURE, normalized to the Renyi-3 scale NORM = m^2/((log m)^2 L^3), and run the blind world test; contract PRIME.L2.RENYI3.SIGNED_CUBIC_FLUX.01, experiments-side, NO ledger row).  THREE SEALED CANDIDATES: Phi3a = NORM(|COLL| + |BND| + |DFLUX|) [reviewer raw], Phi3b = NORM|COLL| [combined signed collision sum], Phi3c = NORM|COLL|*FCIX [flux-corrected]; honesty conditions machine-audited (AST identifier scan: no cubic-target read-back; sealed literal scan: no r314 world-table number in any builder -- the e1/e2 mutants prove both audits bite); the literal |cube - DeltaF| form is adjudicated READ-BACK-ADJACENT and demoted to a diagnostic column (== Phi3b to 1.6e-15).  VERDICT PHI3_ALL_BLIND on the sealed letter: no candidate passes all four reviewer criteria -- C0 a/b/c = 2.6261/1.5052/0.9400 (first-5 freeze, r306 protocol); K3 boundedness FAILS 2/57 on all three candidates and K2 world-rank FAILS, and the ENTIRE failure sits on kz55 + kz67 (a 4.61/2.73, b 2.55/1.89, c 2.43/1.73; kz55 alone tops SCRAMBLE 2.39 and blows the b/c bands 66/88 vs bar 30) -- the SAME near-critical family as r306 A<=1 and the r313 T1/T4 constants, with every trend FALLING (-0.40..-0.77): a shallow-calibration artifact of the known kind, but the mid-ladder rule was not sealed and is not applied post hoc.  The LOCAL CAUSE is named by the probe's own coordinates: kz55/kz67 carry FCIX 0.955/0.915 vs med 0.629 -- exactly the rungs where the intra-block flux cancellation DIES (the R316-relevant discovery: the FCIX -> 1 stratum is the obstruction, and FCIX is source-pure and computable in advance).  What HOLDS of the reviewer table: EPSTEIN holds all three frozen bounds (2.3545/1.3615/0.1375 <= C0; Phi3c separates it 7x downward via its world FC 0.101), the rational twin w13 is near-identical to w9 (factor 1.10/1.07/1.17 <= 3.0 -- the strongest twin result of the lane), the SCRAMBLE cause is named COLL (component devs dflux 1.33 / coll 3.69 / fcix 0.10, K4 fires), and candidate a passes its band bar 22.85 <= 30.  Leg C contrasts both fired: the multiplicity-3 control (first atom of every mult-2 fold group split into two half atoms -- G1/blocks/positions/mass preserved BITWISE, dev 0.0) breaks the collision count bound 48762 vs 14448 (uniform toy 216 vs 64, break 152 EXACT) and is rejected by the cap while Phi_3 is provably BLIND to it (the multiplicity cap is a necessary CLASS side condition, disclosed); the assignment shuffle (seed 315001, r313-m4 base) is rejected 289/300 with the flux profile broken edgewise 1.0e0 and mass matched 4.0e-17.  Anchors bit-near r314/r306 (shares -0.4226/+0.5980/+0.8537, FC 0.629/-0.141, mult == 2 on 57/57, identity wards 4.5e-17, r306 C_2 1.069 with 0/57).  Must-fails: e1 cubic read-back AST-CAUGHT / e2 r314-table calibration literal-CAUGHT / e3 mult-2-built control ACCEPTED as normal (sharpness 14448 == 14448) / e4 wrong normalization exponent LOUD (factor m = 35 exact) / m5a+m5b scope-flagged.  NO amendment after freeze except the record-table insertion.  Deterministic run1/run2; honest: ALL_BLIND stands -- the divergence-form structure with first-5-frozen constants is NOT the membership functional; the kz55/kz67 + FCIX-stratum diagnosis is post-hoc census, no rule; nothing here bounds anything; NO RH CLAIM in either direction; consumed by v968 (wave 11, embedded byte-exact, smoke stage)"

_R316_STATUS = "L2 signed-cubic-flux round, part 3 -- THE TWO-REGIME BOUND (35/35; SPEC_SHA 5c28b12b; reviewer contract part 3: formulate the pointwise Renyi-3 bound as a coupled two-regime statement over the source-pure FCIX stratum with theta = 0.85 sealed from the r315 record gap and a MID-LADDER calibration sealed in advance (CAL_START = n//3, N_CAL = 5, small-m ranks certified individually); contract PRIME.L2.RENYI3.SIGNED_CUBIC_FLUX.01, experiments-side, NO ledger row).  FOUR SEALED REGIME MAJORANTS, all >= rho_2 by exact algebra warded live (worst slack 6.5e-16 on 69 worlds): PhiL2 = NORM(|COLL|+|BND|+FE) [flux-index form], PhiL1 = NORM(|COLL|+|BND|+|DFLUX|) [r315 raw], PhiH2 = NORM(CNT3+|BND|+FE) [the multiplicity-2 collision COUNTING majorant, CNT3 = sum_j (n_j+3n_j(n_j-1)) gmax_j^3], PhiH1 = (m QMAX/log m)^2 [concentration]; purity machine-audited (AST + record-literal scans; e1 theta-posthoc, e2 cubic read-back, e2b literal calibration all CAUGHT).  VERDICT TWO_REGIME_DEAD(L fails; H fails) on the sealed letter -- and the anatomy is the round's find: (1) the FCIX stratum is NOT the near-critical family: the mid-ladder freeze is killed on the regime-L side by kz53 (rho_2 1.049 with BULK-NORMAL FCIX 0.654) and deep spikes kz83/kz105 -- the obstruction family cuts ACROSS the flux-cancellation stratum; (2) rho_2 itself is not below its mid-ladder window (kz53 1.049 / kz83 0.779 vs cal-window max 0.458): ANY tight pointwise majorant fails at kz53 under a mid-ladder freeze -- the r306 first-5 protocol was not a shallow-calibration artifact but load-bearing (the profile is non-monotone with recurring near-critical spikes); (3) the REGIME-H MECHANICS OBJECT delivered as census: the FCIX outliers kz55/kz67 are near-ONE-BLOCK worlds -- TOP1 cube share 0.558/0.785 vs 0.18 regime-L med (dev 2.68, the named discriminator; kz55 additionally 10x compensation kappa 0.105); concentration (candidate i), not counting tightness (dev 0.53), distinguishes them -- but kz67 misses the mid-window H constants by 1.4x/2.6x and kz55 fell into the small-m set (rank 20/65), so the H test stratum was ONE rung.  What HOLDS: anchors bit-near (r314 shares/FC/mult, r315 C0 2.6261/1.5052/0.9400 + FCIX outliers 0.955/0.915, r306 C_2 1.069 with 0/57); EXT2 extension clean (8 deep anchors to N_w 1393, all mult 2, NO new H member -- the FCIX stratum did not grow with depth); the 21-rung small-m certificate table (C_small 1.0694 = the r306 shallow maximum); the class machinery intact (EPSTEIN/twin/main admitted, twin band 1.04 = the strongest of the lane, SCRAMBLE rejected via COLL attribution 3.69 + edgewise flux break 290/300 with matched mass, mult-3 control rejected via cap with all four majorants provably blind).  Must-fails: e1 AST+toy CAUGHT / e2 AST-CAUGHT / e2b literal-CAUGHT / e3 split overlap 1 EXACT vs sealed 0 / e4 mult-3 cap + count break 48762 vs 14448 / m5a+m5b flagged.  NO amendment after freeze except the record-table insertion.  Deterministic run1/run2; R317 direction typed census-grade: the regime coordinate must capture the near-critical family itself (TOP1/QMAX separates kz55/kz67; kz53 needs a second coordinate) or the spikes become a certified exception family (r287-F2 pattern); nothing here bounds anything; NO RH CLAIM in either direction; consumed by v968 (wave 11, embedded byte-exact, smoke stage)"

_R317_STATUS = "L2 exception-families round -- reviewer fork (b): THE EXCEPTION-FAMILY CENSUS (38/38; SPEC_SHA 04fbe5c0; reviewer contract after the r314-r316 trilogy: the statement form 'all sufficiently large generic windows + explicitly classified exception FAMILIES', with the HARD GATE sealed in advance -- at most TWO source-pure exception families plus the generic theorem, any uncovered violator / family growth / world leak fires WHAC_A_MOLE by contract; contract PRIME.L2.RENYI3.EXCEPTION_FAMILIES.01, experiments-side, NO ledger row).  TWO SEALED CLASS FUNCTIONALS, both source-pure (AST + record-literal audited; consumed inputs: ONE r316 two-regime-state column + rank order): F_A = rank-local QMAX spike ratio (CONCENTRATION), F_B = rank-local PhiL2 spike ratio (DIVERGENCE-MASS SPIKE), window W = 5; thresholds frozen by the sealed largest-gap rule (k <= 6, gap >= 1.25) BEFORE any bound evaluation, class table printed before any bound table.  VERDICT WHAC_A_MOLE(kz53, kz83) on the sealed letter -- the hard gate fired exactly as the contract demands: the sealed rule recovers B = {kz55, kz67} (the r315/r316 FCIX pair, F_B 7.23/4.96 with a clean 1.78 gap; finite certificates C_B = 1.0536 <= the r306 shallow max, PhiL2-warded, non-growing, 0 EXT2 members) but class A stays EMPTY (best gap 1.233 misses the sealed 1.25 bar by 0.017), and on the 63-rung complement the mid-ladder generic constant C_gen = 0.4579 (frozen at kz34, m_0 = 84) is violated by exactly the two uncovered rungs kz53 (rho 1.0493) and kz83 (0.7791) -- a third exception form would be needed, NO third class added, abort by contract, recommendation back to fork (a).  THE CENSUS BEHIND THE LETTER IS THE ROUND'S FIND: the F_A TOP-3 are EXACTLY kz53 (2.47), kz83 (2.39), kz67 (2.38) -- the rank-local QMAX ratio ranks the complete mid/deep near-critical family on top, REFUTING at census level the r316 conjecture that kz53 needs a second coordinate; but the distribution below is a CONTINUUM (1.93, 1.90, 1.74, ...), not a gap-separated family, and the F_B continuum over-ranks kz12 (2.79, rho 0.38 -- a PhiL2 spike that is NOT a rho spike) above kz53 (2.70): threshold classification on these coordinates cannot cover kz53/kz83 without swallowing harmless rungs -- the exception-family FORM (sealed thresholds) is the wrong statement shape for a spike CONTINUUM.  What HOLDS: anchors bit-near incl. the COMPLETE r316 anatomy (ladder n = 65, H = {55, 67}, split 21|5|39, C_L2 0.7476, the 8-violator set, TOP1 0.558/0.785, rho anchors, C_small 1.0694); generic complement trend FALLS (-0.170, reserve med 2.06 -- with the B pair excepted the rest has growing reserve except at the two named spikes); 21 small-m certificates; world machinery intact (twin band 1.04, SCRAMBLE rejected via COLL 3.69 + edgewise shuffle break 288/300 with matched mass).  Must-fails: e1 third-class-after-sight AST-FLAGGED + seal_family_count REFUSES 3 families + gate maps uncovered violators to WHAC_A_MOLE / e2 threshold-after-sight AST-FLAGGED + toy thr break / e3 class-readback (cm/S3) AST-CAUGHT / e4 circular certificate caught by the declared-set ward EXACT + LOUD toy / m5a+m5b scope-flagged.  Protocol disclosure: prefilled placeholder record tables removed COMPLETELY before the first run (the r316 error class, disclosed in the spec); NO amendment after freeze except the record-table insertion.  Deterministic run1/run2; R318 direction typed census-grade: fork (a) with the QMAX local ratio as a continuous regime COORDINATE (bound rho_2 BY the coordinate, not classify by threshold), or the generic constant at the r306 first-5 scale where no exception family is needed on the measured ladder; nothing here bounds anything; NO RH CLAIM in either direction; consumed by v969 (wave 12, embedded byte-exact, smoke stage)"

_R318_STATUS = "L* INDEFINITE-FORK round -- the NEW idea class after the L*-language stop: PONTRYAGIN/KREIN INDEX vs SIGN REGULARITY (25/25; SPEC_SHA f2d98683; reviewer perspective shift after r312's 97.6/2.4 membership anatomy: 'the signedness itself is the mathematical object that is not yet understood' -- L* reformulated from squeezing out positivity to 'where may the negative signature live', with the core-question inversion: are the 2.4 percent negative cross mixtures the FINGERPRINT of the right indefinite theorem?  ONE theoretical representation test under the binding stop rule (both routes world-blind/restatement => FORK_DEAD, lane stop); contract PRIME.LSTAR.INDEFINITE_FORK.01, experiments-side, NO ledger row).  VERDICT P2_MAIN_SPECIFIC(fingerprint verbatim: the cross-mixture residual of the converged r308 Dykstra family lives LAWFULLY on the ANTIPHASE pair (D3, D4) with fixed sign -1, modal share med 0.699, 12/12 sealed rungs + the r289 rational twin EXACT ((2,3), -1, 0.692 -- METRIC_ONLY holds), while the dead controls break by PATTERN, not merely by bar: EPST/SCR put their residual on the ARCH-MEAN x BORDER pair (4, 5) with shares 0.953/1.000 and d6-class 0.962/1.000 vs MAIN 0.027 -- honest caveat typed with the letter: the control fingerprints are measured on NON-psd-feasible 200-step iterates (feas -0.45/-0.49, the r308 non-convergence), labeled ITERATE; dig site named per the stop rule: WHY does the 2.4 percent negative cross mixture live on the antiphase pair with a fixed negative sign, stably over the ladder and the twin) + INDEX_LANGUAGE(the P1 route closes as LANGUAGE: the index bookkeeping is EXACT -- spectral count #{sigma(B_n) > 1} == mp pivot count #{h_k < 0} at window AND guarded deep depth on all seven worlds (SCR deep 66/65 inside its gray band 6, the single tolerated case) -- but the index statement 'n_+(window) == N_w' is a total RESTATEMENT of L* (equivalence with minC >= N_w on 3 exact toys + 7 real worlds, both truth values realized), the global inequality n_+ >= N_w is VACUOUS (carried by the mu-channel majority S_+/N_w = 1.23..1.96 on every world incl. controls), and the negative-subspace invariants are ALL WORLD_BLIND under the r281 distance rule (NEG_LOW/NEG_MED/NEG_PR; NEG_LOW is NOT a proxy of the crossing location -- the negative directions do not live at minC, but their location carries no world separation either); Krein angular test typed INDEPENDENT (spearman(ANG, rho_win) +0.536 < 0.9, reported only)) + SIGNATURE_TABLE(window index defect w9/w13/twin 0/0/0 EXACT -- n_+ == N_w at the window, the index form of half-filling survival -- vs EPST 55 / SCR 37 / SMOOTH 4 / HL2 31: the control flips ARE negative directions inside the window; deep inertia 43@262 / 36@236 / 43@262 / 80@225 / 66@272 / 6@360 / 58@278; ladder 42/42 window defect 0) + FINGERPRINT(consensus (2,3) sign -1 med share 0.699, per-rung 0.605..0.780 over kz 9/12/13/14/15/18/20/22/23/29/32/33, cone-share med 0.970..0.982 -- the 97.6/2.4 anatomy is ladder-uniform) + R312_DEMARCATION.  P2 SR half honest negative: NO census object is MAIN-specifically sign-regular (E./A.principal score 1.000 but world-blind/not separating; every orientation-sensitive census coin-flip ~0.50 on MAIN) => the variation-diminishing implication chain to L* FAILS AT ITS PREMISE ON MAIN (PREMISE_FAILS_ON_MAIN).  Leg 0 bit-near (w9 367/263/104/184/184, B 8.368649, mp budget 104 == S_-, crossing 185; r312 pins: Dykstra CONV +6.6e-16, cone share 0.9760/0.9520/0.9778, alignment 0.9973).  Must-fails: m1 wrong signature convention 3 != 2 exact CAUGHT / m2 transposed census exact CAUGHT / m3 circular-lever mutant AST-FLAGGED / m4 single-rung consensus REJECTED (>= 10 rungs forced).  Protocol disclosures: prefilled placeholder record tables removed COMPLETELY before the first run (the r316 error class); smoke pass 1 aborted on an fp_consensus tiebreak-type harness bug (fixed at smoke stage); ONE reporting-only amendment a1 (the verdict-letter dig-site label derived from the MEASURED modal pair -- the draft had pre-named the border column, contradicted by d6-class 0.027; no bar, rule or tree moved).  Deterministic run1/run2; R319 direction typed census-grade: the antiphase (D3, D4) negative cross mixture as the object of the indefinite-theorem hunt -- its sign law is world-separating in SHAPE, exactly the precondition the reviewer's inversion needs; nothing here proves anything -- NO L* claim, NO RH CLAIM in either direction; consumed by v969 (wave 12, embedded byte-exact, smoke stage)"

_R321_STATUS = "L2 continuous-coordinate round -- reviewer fork (a) in the r317-sharpened form: THE SLIDING CUBIC BOUND (39/39; SPEC_SHA e68883ad; the r317 abort demanded the statement-shape change: bound rho_2 BY the single source-pure coordinate F_A = rank-local QMAX ratio (r317 import, W = 5) instead of classifying by threshold; contract PRIME.L2.RENYI3.CONTINUOUS_COORDINATE.01, experiments-side, NO ledger row).  SEALED: the exact concentration bracket qmax PhiH1 <= rho_2 <= PhiH1 = (F_A B)^2 with the source-pure baseline B = medloc x m/log m (derived algebra, warded two-sided live on 69 worlds); a three-form MONOTONE g-family with mid-ladder calibration and derivation-strength precedence (G_SQ b F^2 with b = (max cal B)^2 -- consumes NO target value; G_TT c1 + c2 F^3, the r317-named heuristic; G_LIN a F), the fit-free upper-envelope rule on the declared test split, gain bar 1.5 and a four-branch verdict tree sealed in advance.  VERDICT SLIDING_BOUND_GO(G_SQ) on the sealed letter: rho_2 <= 1.3056 x F_A^2 pointwise with 0/39 mid-ladder test violations AND all four named r316/r317 violators INSIDE (kz53/kz83/kz67/kz55 at reserves 7.0..9.6 -- exactly what killed every flat mid-ladder constant since r316 is absorbed by the coordinate); reserve min/med 2.71/5.35, trend -0.341 falling (growing reserve); the test envelope is STRICTLY monotone in F_A (bin Spearman +1.000, top bin 1.0536 at F_med 2.15), bulk Spearman +0.84 -- stronger than the r317 continuum reading; gain 15.95.  CANDIDATE THEOREM (sliding cubic bound) printed: sum q^3 <= 1.3056 F_A(w)^2 (log m)^2/m^2 for m >= 73 + 21 small-m certificates (C_small 1.0694); COROLLARY: F_A <= 2.47 measured => uniform C_impl = 7.97 (disclosed 7.5x looser than the r306 first-5 constant 1.069 -- the round buys FORM: one gliding bound, no regimes, no exceptions, not sharpness).  THE TWO HONEST STRUCTURAL FINDS: (1) the pure-algebra transfer route does NOT close -- B is NOT bounded by its cal max (test max 1.4088 = 1.23x cal max 1.1426, trend +0.122 rising): the G_SQ certificate holds on the measured rho_2 directly because the bracket's qmax slack falls faster than B^2 rises; the proof-side object is the qmax-share rho_2/(F_A^2 B^2); (2) SCRAMBLE is NOT rejected by the coordinate bound alone (its inflated cube comes with an inflated concentration ratio F_ins 2.00, covered at 5.21) -- the rejection is carried by the r317 class side condition (COLL attribution 3.69, seeded shuffle 321001 with 294/300 atoms displaced and matched mass), honestly typed: the sliding bound is not world-separating by itself.  Anchors bit-near incl. the complete r316 anatomy AND the complete r317 census (F_A top-3 (53, 83, 67) = 2.47/2.39/2.38 ordered, gap-B {55, 67} THR_B 3.7157, class A EMPTY, complement C_gen 0.4579 with violators (53, 83)); w9/w13/EPSTEIN admitted at the sliding bound (rho_2 0.458/0.461/0.368 <= g(F_ins) 0.90/1.55/1.29), twin band 1.04.  Must-fails: e1 g-posthoc AST + toy CAUGHT / e2 coordinate cubic read-back AST-CAUGHT / e3 envelope-on-cal declared-set ward EXACT / e4 monotonicity break LOUD 1.0625 / m5a+m5b scope-flagged.  Protocol disclosure: prefilled placeholder record tables removed COMPLETELY before the first run (the r316/r317 error class); NO amendment after freeze except the record-table insertion.  Deterministic run1/run2; R322 direction typed census-grade: the provenance question is now source-pure, local and split in two -- (a) is F_A bounded (max 2.47, the near-critical family is its top), (b) what bounds the qmax-share rho_2/(F_A^2 B^2) -- the median-of-max shape law (the r302 M_2 stationarity 1.973 pins the second moment of the same shape); nothing here bounds anything beyond the measured rungs; NO RH CLAIM in either direction; consumed by v969 (wave 12, embedded byte-exact, smoke stage)"

_R322_STATUS = "L* antiphase-sign-law round -- THE DIG after the r318 fork decision (25/25; SPEC_SHA 761b51d4; contract PRIME.LSTAR.ANTIPHASE_SIGN_LAW.01, experiments-side, NO ledger row; the three dig questions: (1) is the (D3, D4) antiphase sign law of the block-Green cross-mixture residual iteration- and construction-invariant, and what do FAIR control objects carry (the r318 caveat: control fingerprints were ITERATE-grade); (2) is the law identity-forced (exact forced-component analysis of the linear identity constraints on the (D3, D4) sector); (3) does the theorem candidate 'psd cone part + antiphase-supported fixed-sign residual' certify on 57 rungs + twin?).  VERDICT ALGORITHM_ARTIFACT on the sealed tree: the law breaks under the sealed random-start variants -- their accepted near-feasible points (staged to 20000 Dykstra steps, min eig rel -8.5e-9/-1.2e-9, r311 START_NEAR class) carry modal pair (0, 2) with shares 0.310/0.264 and cone share med 0.782: the sign law AND the r312 97.6/2.4 cone anatomy are properties of the LEAST-NORM-PROXIMAL DYKSTRA BASIN, not of the psd solution set as a whole; the canonical protocol variants LSTART/LONG/ZERO/REV all carry the law exactly ((2, 3), -1, share 0.692; projection order irrelevant) AND coincide within rel 2e-7 -- exactly 9/15 variant pairs distinct = precisely the pairs involving the random starts: the sampled psd intersection has (at least) two basins and the law is the fingerprint of the proximal one.  HARDENING of the r318 caveat on FAIR objects: the least-norm census reproduces the r318 shape contrast construction-fairly (w9/twin (2,3) -1 0.698 vs EPST/SCR (4,5) -1 1.000 d6 1.000 -- the r318 control-iterate pattern is the least-norm pattern), BUT the budget-ablated CONVERGED controls split: EPST-abl CARRIES the MAIN law stronger than MAIN itself (share 0.742 vs 0.692, CONV +2.5e-16) while SCR-abl breaks the share bar (0.379, d6 0.371) => FAIR_CONTROL_CARRIES 1/2 -- the r318 pattern dichotomy does NOT survive fair convergence; the separation on fair objects is a mixed share margin, not a dichotomy.  FORCED COMPONENT exact: TOY4 calibrator machine-exact (free fraction 0, value 4/7 == the r308 G10 hand pin); SM1 phi_23 exact free fraction 0.47671 with forced value +5.101e-1 POSITIVE; MINI16 exact rowspace membership FALSE; w9 f64 free 0.635 value +0.12 -- the negative sign law is carried ENTIRELY by the free/selection directions (the disclosed D3 = -(D2 + 2 D1) in-block dependence is the structural reason); exact_forced FALSE.  ANATOMY: sign consistency 1.000 (every (2,3)-dominant block has R[2,3] < 0; 92.3 percent of ALL blocks), position-structured (thirds 0.289/0.860/0.926 -- a mid/tail law of the fold order), MAG scale-stable over the 12 sealed rungs (spearman vs S -0.168).  COUPLING: K1 block-local antiphase mass pairing g_j g_{j+2} + g_{j+1} g_{j+3} spearman +0.8114 (the strongest measured correlate, below the sealed 0.9 bar), K2 5/7-reserve -0.832, K3 r288 carrier |z_v| +0.252 -- all UNCOUPLED.  CANDIDATE: 57/57 rungs + twin certified censally at the sealed bars (CONV + pair + sign + share >= 0.5 + cone med >= 0.95; worst share 0.605@kz23), but demoted by the main letter (basin-selected family) and by L5 (a fair control carries the pattern); missing links L1-L5 all open.  Must-fails: m1 unconverged-iterate census REJECTED as ITERATE (the r318 caveat machine-enforced) / m2 permuted-family census (3,4) != (2,3) CAUGHT / m3 wrong-identity-sign -4/7 != +4/7 CAUGHT exact / m4 read-back mutant AST-FLAGGED.  Protocol disclosures: draft record block removed COMPLETELY before the first run (the r316/r318 error class); ONE amendment a1 (gate typing only: the invariance gate is MEASURED and the sealed tree adjudicates -- the negative branch is a verdict letter, not a probe failure; the letter is identical across passes; no bar, acceptance rule, statistic or tree moved).  Deterministic run1/run2; R323 direction typed census-grade: the selection-geometry question -- WHY does the least-norm-proximal basin organize its indefiniteness on the antiphase pair (best reporting-grade hook: the block-local antiphase mass pairing at +0.81) -- and whether ANY canonical source-pure selection rule reproduces the basin (the r312 COEFFICIENT_SIGN_WALL says no formula so far); nothing here proves anything -- NO L* claim, NO RH CLAIM in either direction; consumed by v969 (wave 12, embedded byte-exact, smoke stage)"

_R324PRE_STATUS = "PRE-WORK PROBE of round r324, SUPERSEDED mid-round by the binding reviewer course correction (the F_A-boundedness hunt demoted as a potentially unnecessarily strong statement -- the cleaner decomposition is QMAX x M_2; the round's binding contract is qmax_m2_origin_probe.py below); stays sealed and citable, its banked record (M_2 stationarity, m2 column + freeze + violator set, share majorant) is consumed as ANCHORS by the r324 probe.  Original record: L2 F_A-provenance round -- the terminal main lane after the r321 SLIDING_BOUND_GO (36/36; SPEC_SHA 9a6696f8; contract PRIME.L2.RENYI3.FA_PROVENANCE.01, experiments-side, NO ledger row; the two provenance halves of the r321 sliding bound executed: (a) is F_A bounded -- distributional route (sealed stationary-null resample: 200 seeded replicates of full synthetic QMAX ladder columns drawn iid from the pooled r302-stationary normalized block profile, donors = ladder minus the named spike ranks, declared-set warded) and structural route (NEW EXACT factorization QMAX = s x a with s = in-block survival <= 1, a = atom-mass share of the argmax block, warded 2.2e-16 live); (b) what bounds the qmax share sigma = rho_2/(F_A B)^2 -- NEW EXACT intermediate inequality sigma <= m2/(m qmax) (from sum q^3 <= qmax sum q^2, Fractions-proved, warded two-sided live on 69 worlds), converting the r302 M_2 stationarity into the pointwise transfer rho_2 <= m2 x F_A x B/log m).  VERDICT FA_BOUNDED_DISTRIBUTIONAL + QMAX_SHARE_OPEN + NOT_COMPOSED on the sealed trees: (1) the r302 stationarity EXPLAINS the coordinate -- KS(F_A real, F_A* resampled) = 0.0982 <= 0.125 and the ladder maximum 2.47 sits INSIDE the stationary max law (p_max 0.090 >= 0.05; envelope C_F 3.357, q95 2.583): F_A <= C_F is a distribution-statement candidate; DISCLOSED tension: per-rank p at the named spikes 0.005/0.000/0.000 -- WHERE the extremes sit (exactly the near-critical rungs) is not explained by the exchangeable null, the rank-identity question is the sharpened residue; (2) the structural census REFUTES a single mechanism: kz53 is a SURVIVAL spike at an edge-adjacent atom-poor block (F_s 1.72, edge distance 0.01, 6 atoms) while kz67 is an ATOM-MASS concentration at an interior atom-rich block (F_a 2.09, 12 atoms), kz83 mixed -- no common carrier, the exact majorant F_A <= a/medloc does not mid-ladder certify (4 violations); (3) share letter OPEN with the mechanism measured: stationarity reproduces on the 65-rung ladder TIGHTER than the r302 core (KS 0.023/0.016 vs 0.043; m2_inf census 2.081 vs rec 1.973), sigma trend -0.584 falling (the r321 'share falls faster than B^2 rises' mechanism measured directly), but C_M2 = 2.2557 is violated by exactly the named spikes (m2 3.05..3.19, +35..42 pct) -- THE F_A SPIKE AND THE m2 SPIKE ARE THE SAME OBJECT seen in the third vs second moment; (4) strongest bycatch: the composed bound rho_2 <= 2.2557 x F_A x B/log m holds 0/39 MEASURED with named coverage 4/4 (reserves 1.23..1.81), LINEAR in F_A where r321 was quadratic -- census, not certificate (the m2 clause fails by seal).  Slim anchor set bit-near (r306 C_2 1.069 0/57, r316 rho quartet + C_small, r317 top-3 EXACT ordered + gap-B, r321 b 1.3056 + 0/39 + named 4/4).  Must-fails: e1 resample-peek AST + toy CAUGHT / e2 spike-donor pool declared-set ward EXACT / e3 share-posthoc rho AST + toy CAUGHT / e4 QQ overlap ward EXACT / m5a+m5b scope-flagged.  Protocol disclosure: prefilled placeholder record tables removed COMPLETELY before the first run (the r316/r317/r321 error class); NO amendment after freeze except the record-table insertion.  Deterministic run1/run2; R325 direction typed census-grade: ONE local question left -- the second-moment excess of the near-critical rungs vs the stationary bulk (m2 tail law via the same resample machinery, or the atom-level genealogy of the heterogeneous argmax blocks); nothing here bounds anything beyond the measured rungs; NO RH CLAIM in either direction; consumed by v970 (wave 13, embedded byte-exact, smoke stage)"

_R324_STATUS = "L2 QMAX x M2 origin round -- the reviewer-sharpened terminal main lane after the r321 SLIDING_BOUND_GO (36/36; SPEC_SHA dc36cacb; contract PRIME.L2.RENYI3.QMAX_M2_ORIGIN.01, experiments-side, NO ledger row; BINDING COURSE CORRECTION mid-round, disclosed: the F_A-boundedness hunt of the sealed pre-work probe (fa_provenance_probe.py, r324-pre) demoted as a potentially unnecessarily strong statement that creates an artificial wall -- for the needed power N_2 >~ m^0.888 asymptotically F_A log m = O(m^{0.112-eps}) suffices; the cleaner decomposition is QMAX x M_2; the banked pre-work reproduced bit-near as anchors and disclosed pre-spec).  S0 EXACT (Fractions + live on 69 worlds): the interpolation M_3 <= qmax x M_2 AND M_2 <= qmax (worst slack 0.0), hence rho_2 <= (m qmax)(m M_2)/(log m)^2 -- the r321 bracket upper improved by exactly the factor M_2/qmax <= 1; the identity F_A x B x log m == m x qmax warded 2.7e-16: F_A is the rank-local normalization of m qmax, NO black correction factor.  S1 M_2-EXPORT: pointwise m M_2 <= C_2 = 2.2557 (mid-ladder freeze) FAILS 7/39 exactly on the banked violator set (kz53 3.1528 / kz28 / kz67 3.1938 / kz109 / kz76 / kz61 / kz83 3.0490) -- the reviewer's 'export should carry' expectation REFUTED pointwise at the mid-ladder freeze (disclosed pre-spec, sealed rule adopted verbatim, not tuned); test trend +0.014 FLAT, envelope C_2env 3.1938; the DISTRIBUTIONAL r302 export GO (stationarity KS 0.0230/0.0158 <= 0.125).  S2 MULTISCALE PILEUP (the reviewer mechanism hypothesis, four measurements on the canonical dyadic source-scale decomposition s = floor(log2(vmax/|atom|)) of the argmax block): (1) scale recomposition EXACT 2.2e-16; (2) fold multiplicity <= 2 admitted on ALL 69 live worlds (the r314-proven cap); (3) per-scale mass C_PIL 9.3583 FAILS 11/39 (worst kz67 22.4) -- clause REFUTED; (4) active scales nsc_rel/log m max 2.026, C_NSC 2.0258 HOLDS 0/39 -- the O(log m) scale count CERTIFIES pointwise; direct m qmax <= C_INF log m (C_INF 1.7481) FAILS 5/39 named 1/4; the exact pileup chain m qmax <= nsc x pil warded 0.0 live.  HARDNESS CENSUS ANSWERED: the named near-critical spikes are NOT multiscale convergence -- ONE near-maximal source scale dominates the argmax block (kz53 s0 = 13.78 of G 13.02 at nsc 5 == ladder med; kz83 s0 12.71; kz67 s0 22.38 at nsc 7): the pileup tip is a SINGLE HEAVY SCALE -- the reviewer hypothesis REFUTED in clause (3) and CONFIRMED in clause (4).  VERDICT PILEUP_GROWS_SUBCRITICAL(+0.172) on the sealed six-letter tree (TARGET_LEAK / FA_RESTATEMENT / QMAX_MULTISCALE_GO / M2_GO_QMAX_OPEN / PILEUP_GROWS_SUBCRITICAL / PILEUP_SUPERCRITICAL): growth exponents fit-free dyadic halves-slope over the 39 test rungs e(G/log m) +0.158 + e(m M_2) +0.014 = e_tot +0.172 < CRIT 0.224 (= 2 x (1 - 0.888); the reviewer's 0.112 is the F_A-side equivalent -- this route is LINEAR in m qmax where r321 was quadratic in F_A, the budget doubles), stability halves +0.141/-0.160 both subcritical (the second half FALLING) -- DECIDED; by the reviewer contract a sufficient first-class outcome, not a failure.  S3 COMPOSITION typed MEASURED: sum q^3 <= 8.941 x (log m) x m^{+0.172}/m^2 => N_3 >= m^{0.914}/sqrt(8.941 log m), N_2 >= N_3 (r306 exact chain) => N_2 >= m^0.888 for all m >= m_0* = 10^59.6 => slope(n_eff) >= +0.908 => sigma <= sigma* = -0.516 => the r297 target => the v964-S0 vdC theorem => delta' > 0.21 => the terminal_positive_main edge; all measured m <= 274 closed by the finite certificates (r306 0/57 + this round's tables); the gap (274, 10^59.6) is the DISCLOSED extrapolation hypothesis -- no cofinal claim.  Slim anchors bit-near (r306 C_2 1.069 0/57, r316 rho quartet + C_small + n 65, r321 b 1.3056 0/39 named 4/4, r324-pre KS 0.0230/0.0158 + m2 med/max 2.051/3.194 + C_M2 2.2557 + the seven-violator set EXACT).  Must-fails: e1 scale-drop recomposition break LOUD exactly 1 / e2 3-fold toy genealogy REFUSED by mult_ward while all live worlds admit at exactly 2 / e3 M_3 read-back AST-CAUGHT / e4 scale-cutoff posthoc AST + toy CAUGHT twice / m5a+m5b scope-flagged.  Protocol disclosures: the course correction documented in the docstring; the banked pre-work outcomes disclosed pre-spec (the genuinely open quantities -- scale profiles, pil/nsc columns, all S2 constants and violation counts, the hardness census, both exponents -- computed only after the freeze); ONE calibration-pass fix BEFORE the record freeze (the m_0* solver searched linear space to 1e18 and missed the honest 10^59.6 -- moved to log space, the sealed rule unchanged); NO amendment after freeze except the record-table insertion.  Deterministic run1/run2 (byte-identical up to the runtime line); R326 direction typed census-grade: the single-heavy-scale anatomy of the near-critical argmax blocks (what caps the TOP source scale's mass at one block -- the mult-2 cap holds uniformly, so the cap must come from the per-group mass law, the local mass balance of the dominant fold group), and whether the subcritical exponent is stable under deeper anchors; nothing here bounds anything beyond the measured rungs; NO RH CLAIM in either direction; consumed by v970 (wave 13, embedded byte-exact, smoke stage)"

_R325_STATUS = "Extraction-repair fork round -- the new main strand after the R319 red team: the three reviewer variants for repairing the RH-connection seam (rh/ cofinality runs in the ANCHOR direction; H_cof demands the MESH-refinement direction, measured inadmissible with false floors) measured in their machine-checkable halves (18/18; SPEC_SHA 31277f91 final with record, freeze-run 57e50d36; contract PRIME.EXTRACTION.ORDER_REPAIR_FORK.01, experiments-side, NO ledger row).  VERDICT: primary ELEMENTWISE_STABILIZATION_GO + co-letters CANONICAL_MESH_REBUILD_GO + SIGNED_DEFECT_CORRECTION_GO(SIGN) + ANCHOR_ONLY_INSUFFICIENT on the sealed tree.  (A) On the NATIVE dense class (dyadic step-function autocorrelations, the v749 'Weil form of step functions' class, elements sealed by seed) ALL THREE channels (comb/arch/pole) of the canonical tower windows stabilize EXACTLY at the predicted finite anchor onset alpha* = (n_g+1)D0/2 (comb rel dev 0.0/3.3e-18 above onset, 2.7e-2..5.3e-2 below -- real onset) and the value is CONSTANT under dyadic mesh refinement (comb 3.3e-18, arch 1.5e-15, pole 2.0e-17): the elementwise architecture consumes NO mesh-cofinal ladder and NO transport -- the R310 finite_forms_converge_to_weil shape generalizes from the Lean comb channel to every channel on the native class, incl. the derived tent-read identity L_cat(F) = -sum mu_n (I_D F)(u_n) gated at 1e-12.  (B) v749 T1.4b reproduced: every deployed frame-A window IS a directly rebuilt tower member at rel dev 0.0 (three picks); the direct rebuild ladder is PD per stage with honestly falling margins (lambda_min/c_0 1.33e-5 -> 3.34e-6 -> 6.92e-7, reported not consumed); the dyadic transport identity holds at 8.0e-17 but variant B never uses it.  (C) On the sealed non-native C^2 class (quartic L=2.3, Hann L=2.9) the per-channel quadrature defect is ONE-SIGNED (comb NEGATIVE, arch POSITIVE, pole POSITIVE, all stages), falls at the D^2 class rate per channel (ratios 3.6..9.0; the TOTAL ratio 2.2..3.3 is reported only -- comb/arch sign opposition cancels there, an honest structural record), and sits inside the RIGOROUS closed-form interpolation envelope |E_ch| <= (D^2/8)||F''|| K_ch -- controlled on its own, no wall margin anywhere; the C^0 tent element's comb defect is EXACTLY 0 (its off-grid kink cells contain no prime-power atom, recorded fluke); the anchor-only floor is real: at fixed mesh the anchor ladder stabilizes at E_tot CONSTANT +1.77e-4/+1.59e-4 (the R319 false-floor mechanism on this probe's own elements).  Protocol disclosures: prefilled record block removed before the first run (r316/r317 class); three smoke-harness fixes BEFORE the freeze (grid-element closing zero knot -- np.interp jump, the kernel reads were wrong by ~1e-6 while the lag identity was machine-exact; decay gate retyped total -> per-channel after the smoke run measured the sign cancellation; smoke index); only post-freeze change the record insertion.  Must-fails m1 anchor-only adversary / m2 nearest-neighbour atom placement / m3 corrupted transport weights / m4 membership corruption all CAUGHT.  Deterministic run1/run2.  The fork consequence (analytic, next.txt note): the repair statement for the next Lean round is the elementwise quantifier set (canonical family + per-element predefined finite stabilization + window-local positivity premise), NOT a new mesh theory; nothing here proves any window positivity, any local lemma, or anything about zeta zeros; NO RH CLAIM in either direction; consumed by v970 (wave 13, embedded byte-exact, smoke stage)"

_R327_STATUS = "L2 fold-group mass-cap round -- the r324 follow-up question on the terminal main lane: WHAT caps the mass of the dominant fold group of the argmax block (34/34; SPEC_SHA 11e4fd401ff77583 final with record, freeze 71f8b7b423a3f847; contract PRIME.L2.RENYI3.GROUP_MASS_CAP.01, experiments-side, NO ledger row; the r324 hardness census had located the near-critical spikes in a SINGLE heavy source scale -- this round descends one level to the fold GROUP).  NEW EXACT MACHINERY (module-own, source-pure, all warded live on 69 worlds): the per-group mass ledger on the identical r314 fold segmentation (G1/mult/gblk == genealogy at 0.0), the partition sum_g gabs == A1_j (3.4e-16), the L1 recomposition sum_j |sum_g G1_g| == L (0.0), the TWO-ANCESTOR bound gabs_g <= mult_g x gmax_g <= 2 x vmax (slack 0.0 -- the proven mult-2 cap as a mass statement) and the exact CAP CHAIN m qmax <= ng x hgn (slack 0.0).  Q1 ANATOMY ANSWERED -- sealed letter SOURCE (med named ratio 1.057 >= RATIO_BAR 0.5): the r324 heavy scale IS a single fold group at the two sharpest spikes, and the STRUCTURAL FINDING of the round: on ALL 65 rungs the heaviest group of the argmax block is EXACTLY ONE beta/omega fold pair (window-atom histogram [0, 65, 0]: one bulk + one window atom at one position) with median alignment 1.000 -- the two ancestors REINFORCE, no internal cancellation; kz53: the single pair carries 88.8 pct of the argmax block's atom mass at gap 0.076 (the second group 13x lighter) -- the spike is literally ONE bulk/window coincidence; kz67 is the exception shape (6 groups, bshare 0.416); position thirds [21, 29, 15] mid-heavy (the r322 mid/tail echo).  Q2 CAP CANDIDATES -- VERDICT CAP_PARTIAL on the sealed six-letter tree (TARGET_LEAK / RECOMP_BREAK / LAMBDA_PAIR_CAP_GO / BALANCE_CAP_GO / HEAVY_GROUP_UNBOUNDED / CAP_PARTIAL): (i) the lambda-pair route hga <= 2 vmaxb is EXACT but the source-size coordinate lvb = 2 m vmaxb/L FAILS the mid-ladder freeze 19/39 (C_LV 1.1838; kz55 6.85 / kz53 5.77 / kz83 4.61 / kz67 4.13, named 0/4) -- the von-Mangoldt size analogy does NOT close as a polylog cap; (ii) the direct group cap hgn <= C_HG (log m)^A fails at A = 1 (18/39) AND A = 2 (8/39, kz53 3.2x), minimal A = None; (iii) THE GROUP COUNT CERTIFIES: ng(argmax block) <= 2.6351 x log m with 0/39 violations and named 4/4 (ng med 5 max 13) -- the third certified O(log m) count of the lane (after the r324 scale count 2.0258).  EXPONENTS (fit-free dyadic halves-slope): e(hgn/log m) = -0.200 FALLING and e(lvb/log m) = -0.261 FALLING (the heavy-group coordinate DECAYS with depth -- the spike set is a shallow/mid phenomenon, not growth), e(ngl) +0.293, e_route +0.106 < CRIT 0.224 but halves +0.411/-1.183 STRADDLE -- UNDECIDED by the sealed stability rule: HEAVY_GROUP_UNBOUNDED does NOT fire, and the honest reading is that the cap route is subcritical in the median but pointwise-blocked by the same named spike family as every mid-ladder constant since r316.  Q3: no certifying cap -> the r324 MEASURED composition stands unchanged (sum q^3 <= 8.941 (log m) m^{+0.172}/m^2, N_2 >= m^0.888 for m >= m_0* = 10^59.6, the gap (274, 10^59.6) the disclosed extrapolation hypothesis); measured cap-route envelope max(ngl x hgL) 2.616.  Slim anchors bit-near (r306 C_2 1.069 0/57, r316 n 65 + rho quartet + C_small, r324 C_NSC 2.0258/C_PIL 9.3583/C_INF 1.7481 with violation counts 0/11/5 EXACT + named scale masses 3/3 + e_tot +0.172, r324-pre C_M2 2.2557 + seven-violator set EXACT).  Must-fails: e1 group-drop recomposition break LOUD exactly 1 / e2 wrong-window-max two-ancestor break exactly 1 CAUGHT / e3 cap-posthoc AST(rho) + toy CAUGHT twice / e4 qmax-readback AST-CAUGHT (the q_max record cannot enter the cap side) / m5a+m5b scope-flagged.  Protocol: smoke pass 1 34/34 (0.4 s) NO amendment; calibration pass 1 = first full evaluation 34/34 (36.4 s) NO amendment; only post-freeze change the record-table insertion; deterministic run1/run2 byte-identical up to the runtime line.  R-next direction typed census-grade: the spike cap, if it exists, is a statement about WHERE the single heavy beta/omega coincidence can sit (the bulk/window pairing at one position -- a source-geometry question), not about group counting (certified) or group multiplicity (proven); alternatively the honest end state of the route is the r324 MEASURED composition.  Nothing here bounds anything beyond the measured rungs; NO RH CLAIM in either direction; consumed by v970 (wave 13, embedded byte-exact, smoke stage)"

# (path, role, round, ledger_ids, status, pin) -- kept in exact
# sync with rh/INVENTORY.json: every rh-sync/promotion wave that
# touches the inventory MUST extend this list in the same change
# (a stale list here would silently DROP rows on regeneration).
ENTRIES = [
    # -- certified RH theorem modules (load-bearing suite) --
    (f"{VER}/v955_tau_iiks_toda_dictionary.py", "verification_module",
     "r224-r226",
     ["PRIME.PORT.TAU.FINITE.IIKS.01",
      "PRIME.PORT.TAU.NOPOLE.COFINAL.01",
      "PRIME.PORT.TAU.01"],
     "[E] finite IIKS/Toda dictionary; [O] no-pole cofinal", True),
    (f"{VER}/v956_signedmoment_halffilling_duality.py", "verification_module",
     "r228-r231",
     ["PRIME.FREEMOMENT.JFRACTION.01",
      "PRIME.JACOBI.DUAL.REVERSAL.01",
      "PRIME.FREEMOMENT.POSITIVEPREFIX.01"],
     "[E] half-filling/duality/L-gauge; [O] positive prefix", True),
    (f"{VER}/v958_bordered_tau_readout_dictionary.py", "verification_module",
     "r244-r253",
     ["PRIME.PORT.RHP.BORDERED.READOUT.01",
      "PRIME.PORT.RHP.FULLSOURCE.BASEFIBER.01"],
     "[E] bordered dictionary + error formulas + PSD base; [O] "
     "base/fiber campaign", True),
    (f"{VER}/v959_coupledtau_terminal_dictionary.py", "verification_module",
     "r256-r259",
     ["PRIME.PORT.RHP.COUPLEDTAU.TERMINAL.01",
      "PRIME.PORT.COUPLEDTAU.TERMINAL_CROSSRATIO.01",
      "PRIME.PORT.FULLSOURCE.PREFIX_RESUMMATION.01"],
     "[E] coupled-tau recursion + bilinear form + telescope; [O] "
     "the two last edges", True),
    (f"{VER}/v960_terminal_surface_closure.py", "verification_module",
     "r260-r275",
     ["PRIME.PORT.COUPLEDTAU.SURFACE_CLOSURE.01",
      "PRIME.PORT.COUPLEDTAU.TERMINAL_CROSSRATIO.01"],
     "[E] terminal-surface closure: 42/42 census certified (41 "
     "mechanism + kz15 exact-finite), two-branch theorem, phase "
     "bounds, universal pair theorem H1-H4 (H5 open); [O] cofinal "
     "front (lemma L2)", True),
    (f"{VER}/v961_midpoint_orientation_dictionary.py", "verification_module",
     "r274-r278",
     ["PRIME.PORT.RHP.MIDPOINT.ORIENTATION_DICTIONARY.01",
      "PRIME.PORT.RHP.MIDPOINT.ORIENTED_THEOREM.01",
      "PRIME.PORT.WALL.METRIC_FIREWALL.01"],
     "[E] midpoint-orientation dictionary (Casoratian = h_n, "
     "augmented telescope) + Maslov census GO (R2, atom-Sturm "
     "refuted) + metric-firewall measurements; [O] oriented "
     "midpoint theorem (rounds 279-281 executed and consumed in "
     "wave 5 via v962) + global metric firewall", True),
    (f"{VER}/v962_halffilling_pinning_theory.py", "verification_module",
     "r279-r281",
     ["PRIME.WALL.HALFFILLING_PINNING_THEORY.01",
      "PRIME.PORT.RHP.MIDPOINT.ORIENTED_THEOREM.01",
      "PRIME.PORT.REPRESENTATION.CONTEST.01",
      "PRIME.PORT.RHP.FULLSOURCE.QUASIDEFINITENESS.01"],
     "[E] the half-filling pinning theory: T1 moment counting "
     "(free pivots exactly h_0..h_{N_w-1}, why-half-filling "
     "answered by counting) + T2 crossing budget (#(h<0) = S_-, "
     "Jacobi/Sylvester, world-blind) + T3 two-sided parity "
     "(h-blind, census bilanz + 87376-case exhaustion) + T4 main "
     "window reduction (the entire open statement is minC >= N_w "
     "<=> forall n < N_w : h_n > 0) plus the four refutations "
     "NO_UNIVERSAL_O1_PINNING / NO_EXTREMALITY / "
     "NO_GENERIC_MASLOV_OBSTRUCTION / NO_SIMPLE_OFFSET_LAW; probes"
     " r279/r280/r281 embedded byte-exact (smoke stage); Lean: "
     "free_window_positivity = the fog-free central sorry, T1/T4 "
     "proved, T2 stated, upper_pinning_not_universal guard; [O] "
     "representation contest + full-source quasi-definiteness (in "
     "flight, not consumed)",
     True),
    (f"{VER}/v963_lstar_reduction_dictionary.py", "verification_module",
     "r282-r285",
     ["PRIME.LSTAR.REDUCTION_DICTIONARY.01",
      "PRIME.LSTAR.SUBORDINATION.01",
      "PRIME.PORT.REPRESENTATION.CONTEST.01",
      "PRIME.PORT.RHP.FULLSOURCE.QUASIDEFINITENESS.01"],
     "[E] the L* reduction dictionary: the r283 A2 chain (mu-frame "
     "congruence minor_k(D_mu - G) == D_k(mutilde) exact; frame "
     "contraction h > 0 through the window <=> lambda_max(E_m) < 1;"
     " crossing exactly at minC + 1; pigeonhole ceiling minC <= "
     "S_+; capacity-as-counting refuted by the rank-one pair) => "
     "L* IS THE CANONICAL REDUCTION (int p^2 dnu < int p^2 dmu for "
     "deg p < N_w, equivalently lambda_max(E_{N_w}) < 1); the r285 "
     "decomposition bookkeeping (lambda = maxdiag x (1 + assist) "
     "exact, budget equivalence sign-exact); the r282 four-language"
     " elimination as named negative gates (CONTEST_ALL_DEAD -- SOS"
     " iff empty negative register, Kasteleyn orientation iff S_- ="
     " 0, Hamiltonian-PSD == h > 0, dual-pair sync by theorem; the "
     "common reason: every classical language forces positivity "
     "exactly for the positive measure class); typed measurements: "
     "one-atom-vs-collective dichotomy, (D) separates, sub-"
     "classical p = 0.38, ensemble LOW_OUTLIER, the first two "
     "MAIN_SEPARATING detectors, the honest DCXX margin decay (min "
     "1.4175e-7 at z = 233; r286 in flight, NOT consumed); probes "
     "r282/283/284/285 embedded byte-exact (smoke stage); Lean: "
     "lstar_subordination = the wave-6 canonical sorry, "
     "lstar_implies_free_window PROVED; [O] the L* contract "
     "PRIME.LSTAR.SUBORDINATION.01 with the standalone document "
     "rh/problem/",
     True),
    (f"{VER}/v964_lstar_coherence_census.py", "verification_module",
     "r286-r289",
     ["PRIME.LSTAR.COHERENCE_CENSUS.01",
      "PRIME.PORT.L2.VDC_LEMMA.01",
      "PRIME.LSTAR.SUBORDINATION.01",
      "PRIME.PORT.COUPLEDTAU.TERMINAL_CROSSRATIO.01"],
     "[E] the L* coherence census: the r286 margin-scaling "
     "clarification (15 new anchors N_w 942-1218 ALL mp-sign-safe "
     "positive, min +1.806e-8 -- NO counterexample; O(1) census "
     "offset survives, max +4; flattening power law alpha ~ 3.05, "
     "driver c_w -> 1; HARMLESS quantified; q_N reconciled; "
     "EXTRAP_CALIBRATED 15/15), the r287 vdC theorem (the exact, "
     "constant-free, arithmetic-free van der Corput inequality at "
     "H = ceil(sqrt(m)) delivers delta' = +0.309 > 0.21 world-blind"
     " and certifies 6/7 + 38/42; F1 discrepancy dead; module-own "
     "exact proof route: Fejer window-sum identity + covering "
     "identity + Cauchy-Schwarz, wrong-prefactor must-fail), the "
     "r288 carrier map (antiphase next-nearest ARCH-ARCH pairs, "
     "z_v = -3.149, total control reversal, finest alignments below"
     " phase resolution; honest negatives SAMPLING_BLIND / "
     "SOURCE_SEPARATOR_NOT_FOUND / DIFFERENT_OBJECTS) and the r289 "
     "METRIC_ONLY adjudication (the rational twin keeps the full "
     "signature identical -- the diophantine route excluded; metric"
     " coherence threshold 1e-3..3e-3 gap; tent-split fractions the"
     " only sub-gap entry, completeness 4.2e-14; Baker ~7900x too "
     "weak AND unnecessary); probes r286/287/288/289 embedded "
     "byte-exact (smoke stage); Lean: no new statement -- the hole "
     "stays lstar_subordination; [O] PRIME.PORT.L2.VDC_LEMMA.01 "
     "(the chain origin of the P-variance scaling) + the L* "
     "contract at the wave-7 state (open front: the profile "
     "functional; r290 in flight, not consumed)",
     True),
    (f"{VER}/v965_lstar_curvature_arc.py", "verification_module",
     "r290-r295",
     ["PRIME.LSTAR.CURVATURE_ARC.01",
      "PRIME.LSTAR.CLOSED_FUNCTIONAL.01",
      "PRIME.LSTAR.SUBORDINATION.01"],
     "[E] the L* curvature arc: the r290 basin geometry (the working "
     "set is a soft-shouldered anisotropic TUBE in the exact theta_eq "
     "LAG coordinate -- killfrac 0.38/0.62/1.00 at 5e-4/1e-3/2e-3; "
     "world axes kill 5-50x earlier; the SMOOTH axis a privileged yet "
     "gradient-orthogonal killer; the r280 ridge REAL, minC 185 at "
     "factors 1..8), the r291 ridge anatomy (the lift is a FIRST-ORDER "
     "BUDGET phenomenon with one threshold in (1.280, 1.291] over all "
     "18 matched-dose cases, k_min 9, one TOP6@8 retraction; NO fixed "
     "point -- plateau 185; LIFT_MAIN_SPECIFIC; SMOOTH collective-"
     "quadratic ratio 23.7), the r292 curvature spectroscopy (all 29 "
     "diagonals negative, RANK-1 DENS valley 92.5 pct, lam_top -0.418;"
     " not SMOOTH, not the ridge; m* and the retraction invisible to "
     "jets; EPSTEIN second-order structureless), the r293 metric "
     "reconciliation (three metrics sealed first; F10 home margin "
     "+0.024, MIX_IS_CAUSE; all 8 flips SIMPLE h_184 zeros alpha = 1; "
     "MSTAR_NO_LAW; the TOP6 retraction a second simple zero at f_ret "
     "7.107) -- plus the module-own exact S0 (polarization identity, "
     "budget telescope, doubling law, sealed decision bars, pure "
     "Fractions with mutant must-fails); probes r290..r295 embedded "
     "byte-exact (smoke stage); Lean: no new statement -- the hole "
     "stays lstar_subordination; HONEST NEGATIVES load-bearing: "
     "ALL_FUNCTIONALS_BLIND over five class families, r294 "
     "F10_FRAGILE (win 5/5, part 2/5) and r295 F10_SP_MAJORITY "
     "(14/20) -- F10 UNPROMOTED per the sealed bars, R293_LUCK, "
     "L2_VIA_CONSERVATION near-tautological; [O] "
     "PRIME.LSTAR.CLOSED_FUNCTIONAL.01 (loss-corpus forensics, the "
     "rank-2 DENS core, why L2) + the L* contract at the wave-8 state "
     "(no round in flight at this cut)",
     True),
    (f"{VER}/v966_l2_reduction_chain.py", "verification_module",
     "r296-r300",
     ["PRIME.L2.REDUCTION_CHAIN.01",
      "PRIME.L2.NEFF_TARGET.01",
      "PRIME.PORT.L2.VDC_LEMMA.01",
      "PRIME.LSTAR.CLOSED_FUNCTIONAL.01",
      "PRIME.PORT.COUPLEDTAU.TERMINAL_CROSSRATIO.01"],
     "[E] the L2 reduction chain: the r296 DENS fork closed honestly "
     "(DENS_WORLD_BLIND -- the coupling number cos(e_top, grad "
     "lambda_max) = +0.394 below the sealed 0.40 bar, lam3 couples "
     "harder -0.574, every arithmetic candidate <= 0.38, the moment "
     "subspace construction-adjacent on both dead controls; the lane "
     "routed to L2 per the pre-adjudicated reviewer fork), the r297 "
     "target inequality (sigma <= sigma* = 2(sl_c2 - 0.21) - sl_pref "
     "= -0.516, measured -0.714, margin 0.198 -- provenance missing, "
     "not room; CHAIN_PROVENANCE_PARTIAL(B3), the chain exactly "
     "orthogonal on the WINDOW measure), the r298 exact sign-"
     "preserving transfer (S_F = B(omega,omega) + B(Delta,omega+beta) "
     "at 8.8e-16 on 47 worlds; TRANSFER_DOMINANT: the window main "
     "term EMPTY, S_F ~ B(PDelta,PDelta) -- the vdC input IS the "
     "Fejer energy of the difference measure), the r299 decay split "
     "(LOWPASS 0.93; FULL-SUPPORT overlap 42/42 -- Delta a pure "
     "c-value difference; ET/Abel dead +1.948; the live edge sl_D "
     "-0.571 <= sigma* margin 0.055 with falling ratio -0.168) and "
     "the r300 diagonal anatomy (sl_D = 2 sl_L1 - sl_neff EXACT; "
     "both magnitude routes closed honestly; RATIO_BOUNDED_"
     "STRUCTURAL: R_env falls -0.122) -- plus the module-own exact "
     "S0 (sigma*-composition, Fejer-block decomposition identity, "
     "participation identity with the machine-checked DIAG<->NEFF "
     "equivalence, kernel envelope with equality witnesses, the five "
     "sealed verdict bars as exact decision logic, pure Fractions "
     "with tipping mutants); probes r296..r300 embedded byte-exact "
     "(smoke stage); Lean: no new statement -- NEFF_TARGET is a "
     "measured ladder-slope aggregate, not an exact finite identity, "
     "and a universal ladder form would be refutable (the r273 "
     "guard pattern); the hole stays lstar_subordination; THE CHAIN: "
     "[NEFF_TARGET open] => DIAG => sigma <= sigma* => r297 target "
     "=> v964 vdC => delta' > 0.21 generic -- ONE inequality "
     "remains; [O] PRIME.L2.NEFF_TARGET.01 (slope(n_eff) >= +0.908, "
     "measured +0.963, margin 0.055; D_rank bridge sp -0.81 "
     "correlational; r301 in flight, not consumed)",
     True),
    (f"{VER}/v967_l2_cascade_closure.py", "verification_module",
     "r301-r304",
     ["PRIME.L2.CASCADE_CLOSURE.01",
      "PRIME.L2.NEFF_TARGET.01",
      "PRIME.L2.REDUCTION_CHAIN.01",
      "PRIME.PORT.L2.VDC_LEMMA.01",
      "PRIME.PORT.COUPLEDTAU.TERMINAL_CROSSRATIO.01"],
     "[E] the cascade closure and the documented lane stop: the "
     "wave-9 rest targets executed to their end -- r301 NEFF_SPLIT "
     "(n_eff = n_act/(1 + CV^2) EXACT with the perfect count link "
     "n_act == m on 42/42; the Renyi family at ONE exponent: echt "
     "anti-concentration; jackknife-stable margin; the rest "
     "relocated onto UNIF_TARGET) and r302 UNIF_DERIVED(B1 + A2), "
     "the FIRST DERIVED of the lane (PROFILE_STATIONARY with an "
     "exact 1/N transient onto 1.973; the coherence identity 1 + "
     "CV^2 = n_act chi/(surv^2 n_eff_atom) EXACT; chi = 0.63 "
     "destructive, falling; the rest relocated onto ATOM_TARGET) "
     "-- then the r303 REGRESS AUDIT: all four margins are ONE "
     "algebraic number S = sigma* - sl_D = +0.0547 (invariance <= "
     "9e-16; the 1/2-conjecture refuted) => REGRESS_CONFIRMED, the "
     "cascade r297 -> r302 RETYPED as an exact reduction DICTIONARY "
     "(coordinate finding, not six proof steps; the hard rule "
     "binds: a round counts only with NEW information); the first "
     "causal coordinate found (rho_1 ladder monotone on 1008 "
     "sealed builds, the flip KILLS the inequality at margin "
     "-0.044) but MIXING_INSUFFICIENT (chi-level miss 0.134); "
     "n_eff_atom a pure MARGINAL functional (1e-15 on all builds); "
     "and the r304 SHORT-RANGE LAW: the global lag profile is a "
     "STABLE, WORLD-SPECIFIC PERIOD-4 COMB, no k0 <= 8 -- "
     "LONGRANGE_STRUCTURE + LAW_LONGRANGE, the sealed reviewer "
     "stop case fires; the reviewer condition SPLITS (NC(16) = "
     "0.712 < 1 HOLDS, summability FAILS at 1.563); the "
     "chi-relevant structure IS short-range (chi = 1 + 2 sum "
     "T_k/Q exact, the r303 gap at k <= 3, lag-2-dominated); "
     "lag-8 matching hits the chi level (0.022) but breaks the "
     "slopes (0.028/0.027): THE MECHANISM IS TWO-SCALE -- plus "
     "the module-own exact S0 (margin-invariance algebra with the "
     "halved-need mutant breaking by exactly sl_L1 = 1/5, "
     "coherence identity in product form, chi-lag decomposition, "
     "zero-sum tautology with the truncated-NC content, period-4 "
     "comb certificate, -/-/+/+ sign-pattern certificate with the "
     "kz18/kz23 Fractions pins, the four sealed verdict bars as "
     "exact decision logic with tipping mutants, the lane-stop "
     "composition gate; pure Fractions, no probe imports); probes "
     "r301..r304 embedded byte-exact (smoke stage); Lean: no new "
     "statement -- the stop carries no new exact target form (the "
     "period-4 comb and the two-scale split are measured "
     "finite-window aggregates, not exact finite identities; a "
     "universal ladder form would be refutable, the r273 guard "
     "pattern); the hole stays lstar_subordination; THE DOCUMENTED "
     "STOP: the global-profile mixing route of the L2 lane is "
     "CLOSED -- L2 generic <=> anti-concentration of an explicit "
     "block field with long-range period-4 structure; return with "
     "new tools; what stands: the exact identity dictionary, the "
     "two-scale split, NC < 1, the sign pattern, the marginal "
     "invariance; [O] PRIME.L2.NEFF_TARGET.01 retyped to the "
     "documented stop state; no round in flight at this cut",
     True),
    (f"{VER}/v968_architecture_adjudication.py", "verification_module",
     "r305-r316",
     ["PRIME.LSTAR.ARCHITECTURE_ADJUDICATION.01",
      "PRIME.SOURCE.RANKONE.UPDATE.IDENTITY.01",
      "PRIME.SOURCE.TENSOR.MECHANISM.01",
      "PRIME.L2.RENYI3.PROVENANCE.01",
      "PRIME.LSTAR.SUBORDINATION.01",
      "PRIME.LSTAR.CLOSED_FUNCTIONAL.01"],
     "[E] the architecture adjudication: the architecture day "
     "r305-r316 frozen in the reviewer's binding FOUR-LEVEL STRUCTURE "
     "-- LEVEL 1 real formal theorems (Lean, already green, NO Lean "
     "change ships with this module: lstar_terminal_implies_master "
     "(r305: L* + terminal q_N < 1 => full master positivity; "
     "augmented_prefix_positive and free_window_positivity "
     "corollaries), the four closed Inertia theorems, the real "
     "PrimeWindow construction (r310), source exactness + mass "
     "conservation of any folding + the strong stabilization form "
     "(r310b), the r309 rank-one update identities, Hill monotonicity "
     "+ the Renyi-3 => N_2 bridge; sorry census 9 -> 5, the two TRUE "
     "holes stand alone: lstar_subordination + "
     "terminal_positive_main); LEVEL 2 certified finite statements "
     "(Renyi-3 on 57 rungs at C = 1.069 / A = 2; fixed-head census "
     "77/77; block-Green identity + strict cone with exact rational "
     "Farkas certificates; paired-cone soundness on 20714 steps; fold "
     "multiplicity exactly 2; 93 percent far-triple cancellation; the "
     "signed cubic identity with vanishing boundary, Fractions-exact; "
     "the 21 small-m certificates); LEVEL 3 negative architecture "
     "decisions (FIXED_HEAD_DEAD; PAIRED_CONE_NO_INDUCTION + "
     "WORLD_BLIND, B does not overtake A; the r308 discriminator = "
     "budget sign / chordal restatement; COEFFICIENT_SIGN_WALL -- Lane "
     "A closed as cone language with the mechanism named: block-psd "
     "membership without rank-one SDD, obstructed at the budget/border "
     "sign; FLOQUET_EXPANDING; TRIPLE_TYPE_MAJORANTS_WRONG; "
     "PHI3_ALL_BLIND; TWO_REGIME_DEAD with the anatomy: the "
     "obstruction family cuts across the FCIX stratum, first-5 "
     "load-bearing, kz55/kz67 near-one-block); LEVEL 4 open mechanisms "
     "(the constructive PSD-G rule as a question, the signed cubic "
     "flux theorem, the quantitative source functional, "
     "crossing_budget, the definitional source bridge, L*, the "
     "terminal) -- plus the module-own exact S0 (Hill/Lagrange chain "
     "with the bridge, the opening-flux lemma F_1 == 0 with the flux "
     "telescope, the collision counting identity with the exact "
     "factor-8 mutant break, an r311-class exact rational Farkas "
     "mini-certificate, the TEN sealed verdict bars as exact decision "
     "logic with tipping mutants, the four-level composition gate; "
     "pure Fractions, no probe imports); probes r306..r316 embedded "
     "byte-exact (smoke stage); THE CLAIM SPLIT (binding, NOT one "
     "TRANSFER title): PRIME.SOURCE.RANKONE.UPDATE.IDENTITY.01 [E] "
     "identities only, PRIME.SOURCE.TENSOR.MECHANISM.01 [O] the open "
     "mechanism; PRIME.L2.RENYI3.PROVENANCE.01 [O]: the GO stands, the "
     "provenance open (trilogy: identity exact, functional blind, "
     "two-regime dead; R317 fork at the reviewer); universality "
     "classes: base fine-metric (METRIC_ONLY, twin invariance), fiber "
     "recursive (EPSTEIN holds, SCRAMBLE breaks) -- shared algebra, "
     "not necessarily shared mechanism; reviewer maturity 4 -> 5 "
     "(formal architecture 9.5) typed as a REVIEWER ASSESSMENT; "
     "nothing is proved cofinally, L* stays THE open center; no round "
     "in flight at this cut (R317 waits on the reviewer fork)",
     True),
    (f"{VER}/v969_forks_and_redteam.py", "verification_module",
     "r317-r322",
     ["PRIME.REDTEAM.EXTRACTION_AUDIT.01",
      "PRIME.L2.RENYI3.SLIDING_BOUND.01",
      "PRIME.L2.RENYI3.PROVENANCE.01",
      "PRIME.LSTAR.SUBORDINATION.01"],
     "[E]/[O] the red-team morning and the two forks: the r319 "
     "extraction-chain audit (three kernel-checked type "
     "inconsistencies three levels above the boxes -- U1 bridge x "
     "terminal via B = -1, U2 bridge x L* via the mesh-level-0 "
     "total node collision and p = X - 1, U3 pair_margin_main "
     "forces its source predicate empty -- plus the mesh-vs-anchor "
     "cofinality seam; the honest chain verdict: two lemmata => "
     "window-local master positivity ONLY) and the r320 repair "
     "(bridge retype with u/B fidelity + separation discipline + "
     "budget_pos; SourceExact for the arch/border finding; "
     "pair_margin_main retyped onto the canonical extraction with "
     "signRuns_sum + canonical_split PROVED; three permanent "
     "sorry-free guards; sorry-free witness at anchor 2, atoms "
     "{2,3,4}, mesh level 4 via 2^7 < 3^5 < 2^8; census 5 -> 5 "
     "with two typed retypes) frozen as "
     "PRIME.REDTEAM.EXTRACTION_AUDIT.01 [E]; the fiber fork frozen "
     "as PRIME.L2.RENYI3.SLIDING_BOUND.01 [O] (the r321 theorem "
     "candidate sum q^3 <= 1.3056 F_A^2 (log m)^2/m^2 for m >= 73 "
     "with all four named violators inside at reserves 7.0..9.6; "
     "honest: C_impl 7.97 = 7.5x looser than r306, form not "
     "sharpness; the provenance question split in two: F_A "
     "bounded? / qmax-share?; the r317 WHAC_A_MOLE abort "
     "documented) and the base fork closed honestly (r318 "
     "P2_MAIN_SPECIFIC with P1 banked as language; r322 "
     "ALGORITHM_ARTIFACT: the antiphase sign law and the 97.6/2.4 "
     "anatomy are least-norm-proximal Dykstra-basin properties, "
     "the selection-geometry question the named residual); the "
     "module-own exact S0 (concentration bracket, the U1-U3 "
     "witnesses in exact arithmetic, the r320 separation/"
     "canonical-split certificates, six sealed verdict bars with "
     "tipping mutants, the wave-12 composition gate; NO Lean "
     "call); probes r317/r318/r321/r322 embedded byte-exact "
     "(smoke stage; r321 imports the embedded r317, r322 the "
     "embedded r318); the Lean rounds r319/r320 consumed as "
     "REPORTS (artifacts in rh/lean/, re-verified by lake build + "
     "run_rh.py); reviewer adjudication addendum: the two true "
     "holes named WINDOW-LOCAL (Level B of the three-level proof "
     "graph A/B/C), the false cofinality direction documented NOT "
     "solved, SourceExact -> CanonicalPrimeWindow + "
     "sourceExact_buildPrimeWindow the named open Lean target; "
     "rounds R324 (terminal QMAX x M_2) and R325 (extraction "
     "repair) in flight, deliberately NOT consumed at this cut",
     True),
    (f"{VER}/v970_extraction_and_composition.py", "verification_module",
     "r323-r327",
     ["PRIME.EXTRACTION.ELEMENTWISE.01",
      "PRIME.L2.RENYI3.MEASURED_COMPOSITION.01",
      "PRIME.L2.RENYI3.SLIDING_BOUND.01",
      "PRIME.L2.RENYI3.PROVENANCE.01"],
     "[E]/[O] the extraction repair and the terminal composition: "
     "the r325 extraction-repair fork (ELEMENTWISE_STABILIZATION_GO "
     "-- the repair statement is a QUANTIFIER statement, not a mesh "
     "theory: on the native v749 class all three channels stabilize "
     "exactly at the predicted finite anchor onset alpha* = "
     "(n_g + 1) D_0/2 and are constant under dyadic mesh "
     "refinement; no mesh-cofinal ladder, no transport, H_cof is "
     "not needed and is replaced as the target route; B = "
     "construction substrate, C = quantified density step, the "
     "anchor-only floor real) and the r326 Lean implementation "
     "(rh/lean/RH/Elementwise.lean, consumed as a REPORT: "
     "sourceExact_buildPrimeWindow PROVED, the comb stabilization "
     "elementwise PROVED, weil_nonneg_of_windowlocal PROVED -- "
     "extraction without the ladder; sorry census 5 -> 8 with three "
     "typed Level-C statements, the wave-12 census reservation "
     "partially discharged) frozen as "
     "PRIME.EXTRACTION.ELEMENTWISE.01 [E] with the honest "
     "new-connection formulation verbatim (two window-local lemmata "
     "over the canonical family + the elementwise architecture + "
     "cited classics => RH; open links: arch/pole transcription, "
     "compression bridge, source completion); the terminal lane "
     "frozen as PRIME.L2.RENYI3.MEASURED_COMPOSITION.01 [O] (r324 "
     "PILEUP_GROWS_SUBCRITICAL +0.172 < 0.224 after the reviewer "
     "course correction with the r324-pre pre-work banked: F_A "
     "de-black-boxed via the exact identity F_A x B x log m == m x "
     "q_max, the O(log m) scale count certified C_NSC 2.0258 0/39, "
     "the per-scale mass refuted -- the near-critical windows are "
     "ONE heavy source scale; the composition typed MEASURED: sum "
     "q^3 <= 8.941 (log m) m^{+0.172}/m^2 => N_2 >= m^0.888 for m "
     ">= m_0* = 10^59.6, the gap (274, 10^59.6) the disclosed "
     "extrapolation hypothesis; r328B independent chain audit "
     "carried as a chain annotation, the record sealed: the bar "
     "0.224 used the ATOM-level need 0.888, the block-level need "
     "is 908/1002 = 0.9062 -- chain-honest bar 0.188, the verdict "
     "survives 0.172 < 0.188, m_0* chain-honest 10^238, exact "
     "gate S0-T5d; r327 CAP_PARTIAL with SOURCE "
     "anatomy: the heaviest group of the argmax block is EXACTLY "
     "ONE beta/omega fold pair on all 65 rungs, kz53 = one "
     "bulk/window coincidence at 88.8 pct; the Lambda-pair cap "
     "refuted, the group count certified C_NG 2.6351 0/39 -- the "
     "third O(log m) count; direction: the coincidence geometry); "
     "the module-own exact S0 (interpolation with slack exactly "
     "1/36, the F_A identity, the r327 group ledger, the r325 "
     "onset formula as an exact toy stabilization, the m_0* "
     "log-space solver rebuilt verbatim with the disclosed "
     "linear-space bug caught, the r328B chain-audit gate S0-T5d, "
     "six sealed verdict bars with "
     "tipping mutants, the wave-13 composition gate; NO Lean "
     "call); probes r324-pre/r324/r325/r327 embedded byte-exact "
     "(smoke stage; r324 imports the embedded r324-pre, r327 the "
     "embedded r324-pre and r324); the r323 abort consumed as a "
     "NOTE (clean, before any write access, nothing counts as a "
     "round) and the r326 Lean round as a REPORT (artifacts in "
     "rh/lean/, re-verified by lake build + run_rh.py); the "
     "window-local positivity premise proved nowhere; the two true "
     "window-local Lean holes unchanged; L* stays paused; the "
     "r328A/B/C adversarial audits read-only, reports with the "
     "coordinator; no measurement round in flight at this cut",
     True),
    # -- read-only deployed core consumed by every probe --
    (f"{VER}/v563_paper2_readouts.py", "readonly_core",
     "phase2",
     ["PRIME.PAPER2.READOUT.01"],
     "[E] deployed builder core", True),
    (f"{VER}/v881_carleson_port_geometry.py", "verification_module",
     "port lane",
     ["PRIME.PORT.SCHUR.01",
      "PRIME.PORT.INTEGRABLE.01",
      "PRIME.PORT.LIMIT.01"],
     "[E]/[O] port geometry (upstream of the tau lane)", True),
    (f"{VER}/v894_diagonal_refinement.py", "verification_module",
     "r26 contract",
     ["PRIME.CASE.PAIRCORR.CONTRACT.01"],
     "[O] the PAIRCORR wall contract (the measured wall)", True),
    # -- the sealed campaign probes (experiments side, frozen) --
    (f"{EXP}/centered_basefiber_probe.py", "sealed_probe",
     "r250",
     ["PRIME.PORT.RHP.CENTERED_BASEFIBER.01"],
     "measurement: OUTER_MODEL_FAILS + FIBER_BEYOND_MODEL", True),
    (f"{EXP}/corner_provenance_probe.py", "sealed_probe",
     "r251a",
     ["PRIME.PORT.CORNER.PROVENANCE.01"],
     "adjudication: CORNER_INTERFACE_ONLY(FALL3)", True),
    (f"{EXP}/targetreadout_error_probe.py", "sealed_probe",
     "r251",
     ["PRIME.PORT.RHP.TARGETREADOUT.ERROR.01"],
     "consumed by v958: READOUT_FORMULAS_EXACT", True),
    (f"{EXP}/base_gauge_constant_probe.py", "sealed_probe",
     "r252",
     ["PRIME.PORT.RHP.BASE.GAUGE_CONSTANT.01"],
     "consumed by v958: NO_GLOBAL_GAUGE + CONTOUR_R1_EXACT", True),
    (f"{EXP}/schlesinger_pairing_probe.py", "sealed_probe",
     "r253",
     ["PRIME.PORT.RHP.FIBER.SCHLESINGER_PAIRING.01"],
     "consumed by v958: TAU_RATIO_EXACT + POLE_REMOVAL_EXACT", True),
    (f"{EXP}/offdiag_gram_probe.py", "sealed_probe",
     "r254",
     ["PRIME.PORT.RHP.FIBER.OFFDIAG_GRAM.01"],
     "measurement: OFFDIAG_EXTENSIVE_CONFIRMED", True),
    (f"{EXP}/nodebands_base_probe.py", "sealed_probe",
     "r255",
     ["PRIME.PORT.RHP.BASE.NODEBANDS.01"],
     "measurement: DRIFT_NOT_BAND_LOCALIZED", True),
    (f"{EXP}/baseborder_factorial_probe.py", "sealed_probe",
     "r256",
     ["PRIME.PORT.FIBER.BASEBORDER.FACTORIAL.01"],
     "consumed by v959: EFFECTS_MIXED_COUPLEDTAU + firewall", True),
    (f"{EXP}/coupledtau_probe.py", "sealed_probe",
     "r257",
     ["PRIME.PORT.RHP.FULLSOURCE.COUPLEDTAU.01"],
     "consumed by v959: COUPLED_RECURSION_EXACT", True),
    (f"{EXP}/budget_anatomy_probe.py", "sealed_probe",
     "r258",
     ["PRIME.PORT.RHP.BUDGET.ANATOMY.01"],
     "consumed by v959: TELESCOPE_EXACT + TERMINAL_Q_LAW", True),
    (f"{EXP}/parametrix_pass_probe.py", "sealed_probe",
     "r259",
     ["PRIME.PORT.RHP.PARAMETRIX.PASS.01"],
     "consumed by v959: PARAMETRIX_FAILED (three honest negatives)", True),
    (f"{EXP}/terminal_crossratio_probe.py", "sealed_probe",
     "r260",
     ["PRIME.PORT.COUPLEDTAU.TERMINAL_CROSSRATIO.01"],
     "measurement: TERMINAL_CROSSRATIO_MEASURED(42/42), no GO; "
     "consumed by v960 (wave 4, embedded byte-exact, smoke stage)", True),
    (f"{EXP}/prefix_resummation_probe.py", "sealed_probe",
     "r261",
     ["PRIME.PORT.FULLSOURCE.PREFIX_RESUMMATION.01"],
     "adjudication: three exact representations, three no-gos", True),
    (f"{EXP}/terminal_triangle_probe.py", "sealed_probe",
     "r262",
     ["PRIME.PORT.COUPLEDTAU.TERMINAL_TRIANGLE.01"],
     "measurement: TRIANGLE_PARTIAL_IMPROVED(38/42); consumed by "
     "v960 (wave 4, embedded byte-exact, smoke stage)", True),
    (f"{EXP}/cancellation_adjudication_probe.py", "sealed_probe",
     "r263",
     ["PRIME.PORT.COUPLEDTAU.CANCELLATION_ADJUDICATION.01"],
     "adjudication: TWO_BRANCH_DECOMPOSITION_EXACT + RHP "
     "reduction; consumed by v960 (wave 4, embedded byte-exact, "
     "smoke stage)", True),
    (f"{EXP}/quenched_opening_probe.py", "sealed_probe",
     "r264",
     ["PRIME.PORT.RHP.QUENCHED.OPENING.01"],
     "campaign charter frozen (SHA 7b9e751d) + RHP dictionary", True),
    (f"{EXP}/s_monotonicity_probe.py", "sealed_probe",
     "r265",
     ["PRIME.PORT.RHP.QUENCHED.S_MONOTONICITY.01"],
     "adjudication: DEFINITENESS_WALL_EQUIVALENT", True),
    (f"{EXP}/border_resolvent_identity_probe.py", "sealed_probe",
     "r266",
     ["PRIME.PORT.RHP.QUENCHED.BORDER_RESOLVENT_IDENTITY.01"],
     "exploration: border-dressed resolvent identity search", True),
    (f"{EXP}/ranktrace_adjudication_probe.py", "sealed_probe",
     "r267",
     ["PRIME.PORT.EXTERNAL.RANKTRACE_CEILING.ADJUDICATION.01"],
     "external adjudication: arXiv:2608.13637 CEILING_IS_OUR_WALL", True),
    (f"{EXP}/drive_local_asymptotics_probe.py", "sealed_probe",
     "r268",
     ["PRIME.PORT.RHP.QUENCHED.DRIVE_LOCAL_ASYMPTOTICS.01"],
     "measurement: DRIVE_LOCAL_PARTIAL(b3F0.20, cert 1/7: kz22); "
     "consumed by v960 (wave 4, embedded byte-exact, smoke stage)", True),
    (f"{EXP}/phase_bulk_bound_probe.py", "sealed_probe",
     "r269",
     ["PRIME.PORT.RHP.QUENCHED.PHASE_BULK_BOUND.01"],
     "measurement: PHASE_BOUND_PARTIAL(c3F0.30, cert 6/7; kz15 "
     "rest 0.12 dec); consumed by v960 (wave 4, embedded "
     "byte-exact, smoke stage)", True),
    (f"{EXP}/kz15_boss_probe.py", "sealed_probe",
     "r270",
     ["PRIME.PORT.RHP.QUENCHED.KZ15_BOSS.01"],
     "adjudication: KZ15_EXACT_ONLY (b3 interval certificate "
     "closes kz15 margin +0.0268 width 1.5e-92 dps 640; b1 level-2"
     " cert 5/7, kz15 miss 0.06 dec; SPLIT_RULE_FAILED) -- surface"
     " 41/42 mechanism + kz15 exact-finite; consumed by v960 (wave"
     " 4, embedded byte-exact, smoke stage)", True),
    (f"{EXP}/universal_pair_theorem_probe.py", "sealed_probe",
     "r271",
     ["PRIME.PORT.COUPLEDTAU.UNIVERSAL_PAIR_THEOREM.01"],
     "theorem form + measurement: UNIVERSAL_STILL_PARTIAL(c2PAIR "
     "5/7; b2LEVEL2 leaves kz39 0.002 / kz15 0.06 dec) + Lean "
     "rh/lean/RH/PairBound.lean; consumed by v960 (wave 4, "
     "embedded byte-exact, smoke stage)", True),
    (f"{EXP}/l2_scaling_anatomy_probe.py", "sealed_probe",
     "r272",
     ["PRIME.PORT.COUPLEDTAU.L2_SCALING_ANATOMY.01"],
     "pure anatomy: TREND_CARRIER(A_pairs, alpha +1.01 beta -0.81 "
     "gamma -0.01) + BOUND_COARSENESS_CONFIRMED (truth rest falls,"
     " margin_true rises; H5 not contradicted) + SLACK_RANKING(c3 "
     "beyond-blind-level-2) + FLIP_CONDITION(delta 0.21 of "
     "available gamma_true 0.45, TRUTH_ALLOWS); consumed by v960 "
     "(wave 4, embedded byte-exact, smoke stage)", True),
    (f"{EXP}/euler_mechanism_probe.py", "sealed_probe",
     "r273",
     ["PRIME.PORT.COUPLEDTAU.EULER_MECHANISM.01"],
     "mechanism round (Euler perturbation ladder, 4 surgeries x 3 "
     "thetas x 3 reps x 42 rungs): PERTURBATION_INSENSITIVE (no "
     "stage collapses gamma_true; ~0.45 is near the generic sqrt "
     "baseline) + FIREWALL_MAP(pp = 0.00 -- every surgery kills "
     "free-prefix positivity at degrees 33-101: the wall, not the "
     "cancellation rate, is the arithmetic fingerprint); consumed "
     "by v960 (wave 4, embedded byte-exact, smoke stage)", True),
    (f"{EXP}/wronskian_dictionary_probe.py", "sealed_probe",
     "r274",
     ["PRIME.PORT.RHP.MIDPOINT.WRONSKIAN_DICTIONARY.01"],
     "exact dictionary D_n <-> Wronskian (32/32): "
     "WRONSKIAN_DICTIONARY_GO (base Casoratian == h_n with c' = 1;"
     " h-free midpoint form == node polynomial L -- the r231 "
     "sign-blindness as an identity; augmented border-paired "
     "Casoratian == h_n S_n, so D_{n+1} = B - W^aug/W^base; h-free"
     " normalization c = 1 through the dual second kind, exact on "
     "toys and at mp on the real comb incl. kz15/kz52 wards) + "
     "ORIENTATION_PREVIEW (Pruefer band 0 < dTheta < pi == h_n > 0"
     " exact; MAIN in-band through the free window; control exits "
     "25/21/27 degree-exact; w9 winding 262/366 measured, first "
     "exit at 184 = N_w) -- the Maslov round is freed; consumed by"
     " v961 (wave 4, embedded byte-exact, smoke stage)", True),
    (f"{EXP}/kyp_memory_probe.py", "sealed_probe",
     "r275",
     ["PRIME.PORT.COUPLEDTAU.KYP_MEMORY.01"],
     "two-sided no-go: KYP_INFEASIBLE(uniform classes -- exact o1 "
     "fiber forcing + o2 descent on 46/46 worlds) + "
     "KYP_IS_TARGET_INVERSE(Riccati memory forced to carry the "
     "exact signed tail sums, budget = Z^2 + sum b^2, fingerprint "
     "sp +1.00; certifies SCRAMBLE -- friendliness signature) + "
     "CURRENCY_READING(QM-AM overhead med +0.15 dec, generic "
     "census 8/42) -- the l^1 -> l^2 currency reform does NOT "
     "carry in state-space/KYP form; consumed by v960 (wave 4, "
     "embedded byte-exact, smoke stage)", True),
    (f"{EXP}/minimal_firewall_probe.py", "sealed_probe",
     "r276",
     ["PRIME.PORT.WALL.MINIMAL_FIREWALL.01"],
     "dose-response round (wall survival depth s = nf/N_w vs "
     "surgery dose; 5 surgeries -- neighbor swap / support jitter "
     "/ family decoupling / grid sign permutation / grid magnitude"
     " permutation -- x SINGLE + theta 0.02..0.25 x windows "
     "w9/kz18/kz55 = 360 worlds, exact h-chain flips "
     "mp-counter-checked at dps 40 incl. the v956 boundary flip AT"
     " N_w = 184): FIREWALL_LAW (P1/P2/B1/B2 GRADED + P3 "
     "INTERMEDIATE; shallow saturating deficit laws D ~ theta^b, b"
     " +0.04..+1.09 -- NOT all-or-nothing) + PROPERTY_RANKING "
     "(support exactness most wall-critical: P2_JIT tol 0.343 < "
     "B2_MAG 0.389 < B1_SIGN 0.536 < P1_SWAP 0.621 < P3_FAM 0.700;"
     " jitter at 2 percent of the local atom gap already costs 3/4"
     " of the depth) + CONTINUUM_VS_JUMP(CONTINUUM, 48/90 MIDDLE "
     "stages: single ops survive at median 0.88-1.00, the v956 "
     "control level 0.11-0.15 is reached only inside the graduated"
     " ladder) + V956_PLACEMENT(min 0.109 at B2_MAG@0.15, "
     "CONTROL_TOUCHED) + N_TRANSPORT(MAP_TRANSPORTS, sp "
     "+0.84/+0.86); consumed by v961 (wave 4, embedded byte-exact,"
     " smoke stage)", True),
    (f"{EXP}/maslov_census_probe.py", "sealed_probe",
     "r277",
     ["PRIME.PORT.RHP.MIDPOINT.MASLOV_CENSUS.01"],
     "blind Pruefer/Maslov census (29/29): MASLOV_CENSUS_GO (rule "
     "R2 = interlacing/reality of the Jacobi zeros, selected by "
     "the sealed training cascade on w9 + 4 smallest-N rungs + "
     "scramble(seed 2) variants; blind 37/37 rungs + w12/w13/w26 "
     "SAFE full depth, controls fire 26/22/28 = flips 25/21/27 + 1"
     " within the sealed +-1; r259 energy-branch separation "
     "17/6/9; NOT h-equivalent -- 79 pattern mismatches vs the 78 "
     "h re-entry pivots) + CENTRAL FINDING: the RAW atom-counted "
     "Sturm census is NOT the winding quantity -- on MAIN it "
     "breaks at n = 56/48 with h positive (zeros escape the atom "
     "hull, then pair inside single gaps: signed Sturm separation "
     "genuinely fails), while the interlacing/reality structure "
     "breaks at exactly nf + 1 everywhere (w9 continuation fire "
     "185 = N_w + 1); STURM_CHAIN_VERIFIED honestly not awarded "
     "(c1 refuted); Sturm rotation #eig < x* == n - #signchanges "
     "exact; kz15 + kz52 mp census wards; the oriented midpoint "
     "theorem is the named next round; consumed by v961 (wave 4, "
     "embedded byte-exact, smoke stage)", True),
    (f"{EXP}/metric_stability_probe.py", "sealed_probe",
     "r278",
     ["PRIME.PORT.WALL.METRIC_STABILITY.01"],
     "metric stability round (from the r276 dose-response "
     "continuum to the local law; exact Hellmann-Feynman gradients"
     " d log h_n/du_j = <wdot_j, pihat_n^2>/h_n and border/CD "
     "gradients dF_n/du_j = -<wdot_j, pihat_n B_n> through the "
     "sealed tent-grid channel, FD/Richardson/mp-gated to 1e-7; "
     "on-node geometry disclosed: the 2-power family sits EXACTLY "
     "on tent nodes (w9 alpha = log 16), one-sided derivative "
     "pairs side-selected): METRIC_STABILITY_LAW(TAME curvature, "
     "validity window theta* med 8.0e-4 [w9] / 1.2e-4 [kz55] -- "
     "25x..170x BELOW the smallest r276 dose 0.02) + "
     "GRADIENT_EXPLAINS_DOSE(PARTIAL, ratio 0.41, flip dev 3.00 --"
     " the 0.02 collapse is already a nonlinear cascade) + "
     "U_PROFILE(PREDICTIVE, sp(A, lethality) -0.82; BOTTOM-loaded:"
     " small primes 2/3/5 carry the gap-weighted sensitivity, hull"
     " correlation -0.80..-0.83, margin map tracks at +0.83) + "
     "N_TREND(SHRINKING, b_L -1.09, decade 17.3 -- no fragility "
     "growth, NO uniform constant) + "
     "MAINWINDOW_PREDICATE(PERTURBATIVE_ONLY: MetricNear(c, MAIN, "
     "theta) <= theta*(w) inherits margin >= margin_MAIN - L_D(w) "
     "theta per window; MAIN positivity stays the open center); "
     "consumed by v961 (wave 4, embedded byte-exact, smoke stage)", True),
    (f"{EXP}/oriented_theorem_probe.py", "sealed_probe",
     "r279",
     ["PRIME.PORT.RHP.MIDPOINT.ORIENTED_THEOREM.01"],
     "oriented-midpoint proof round (32/32): the two-sided index "
     "bilanz (left degree n vs dual degree S-1-n against the S "
     "common nodes of the gauge polynomial L): STEP5_OPEN (no "
     "non-restatement candidate kips at N_w: C1 tight-two-sided "
     "never fires anywhere -- a world-blind THEOREM; C2 union "
     "reality / C3 shared gap / C4 left-in-D die at degrees 1..4 "
     "on every world; C5 = first crossing is the h restatement, 78"
     " pattern mismatches = the r274 re-entry pivots) + "
     "OBSTRUCTION_REFUTED (3 exact rational toy counterexamples: "
     "the two-sided compatibility state does NOT forbid the "
     "crossing) + TWO_SIDED_PARITY_THEOREM (T2 node-identity "
     "signs, 246212 gates; T3 forced-gap parity: union occupancy "
     "ODD in every weight-agreement gap, EVEN in every "
     "disagreement gap, at EVERY degree, h-sign-blind; T3' census "
     "bilanz c + c# = (S-1) - |D| + 2 scD, exhaustion-proved k <= "
     "8) + CROSSING_BUDGET_THEOREM (Jacobi/Sylvester: #(h_n < 0 "
     "over the FULL continuation) == S_- = #negative atoms EXACTLY"
     " on every world: w9 104, w13 98, EPST/SCR/SMOOTH 141/94/6, "
     "kz15 121, kz52 551) + "
     "MAIN_SPECIFICITY(LOCALIZATION_IS_THE_ARITHMETIC: the "
     "machinery is world-blind on all controls pre+post flip; the "
     "ENTIRE arithmetic is WHERE the fixed budget is spent -- "
     "mains/wards at min C = N_w +0/+2/+1/+0, controls at "
     "25/21/27); amendments a1-a3 disclosed (R2 anchor == min C + "
     "1, localization equality NOT universal -- w13 crosses at N_w"
     " + 2); the b3 gap statement (Lean-suitable): forall n < N_w "
     ": h_n > 0, i.e. min C >= N_w -- the localization of the "
     "Jacobi-fixed budget is the named open center; consumed by v962 "
     "(wave 5, embedded byte-exact, smoke stage)", True),
    (f"{EXP}/budget_localization_probe.py", "sealed_probe",
     "r280",
     ["PRIME.PORT.WALL.BUDGET_LOCALIZATION.01"],
     "budget-localization anatomy round (29/29): "
     "LOCALIZATION_CENSUS (ALL 42 rungs: offsets min C - N_w "
     "distribution {+0: 18, +1: 10, +2: 6, +3: 6, +4: 1, +5: 1}, "
     "max +5 at kz43, sp(offset, N) = +0.096 over N in "
     "[142, 878] -- O(1), no N-growth; v956/r279 anchors "
     "0/2/2/3/1 + kz15 +1 + kz52 +0 exact; half-filling N_w == "
     "(S+1)//2 on 42/42; mp dps-40 arbitration exact; controls "
     "25/21/27; dose link: min C == nf definitional, the r276 "
     "P2_JIT depth medians 0.250/0.207 reproduced with the exact "
     "seeds -- dose curve and budget localization are ONE "
     "coordinate) + MAIN_NOT_MAXIMAL (the round's core finding: "
     "the w9 extremal localization min C = N_w is NOT a local "
     "maximum of the localization functional -- the r278-gradient "
     "directions OPT/OPT_SAFE/SMALLPRIME push min C 184 -> 185 "
     "PAST half-filling at theta 7.7e-5..1.6e-4, all mp-"
     "confirmed; the b3 variational hypothesis is REFUTED in its "
     "local-maximum form; kz15/w13 offset windows: first-order "
     "raisers exist, none verifies at the tested doses -- "
     "UNRESOLVED) + CRITICALITY_STRUCTURED (cos_w(grad h_minC, "
     "grad h_minC-1) = -0.97 on all three worlds: the crossing "
     "log-gradient is anti-aligned with the last free pivot -- "
     "the h-sign flip of a raw-gradient lockstep; rho_crit 245 = "
     "the 1/|h_minC| inflation of the shallow crossing) + "
     "DUAL_TRANSLATION_EXACT (toys exact in rationals: "
     "h_{S-m} h#_{m-1} == 1, sign h#_k == sign h_{S-1-k}, D_n == "
     "prod h_k; w9 mp: min C >= N_w <=> ALL S_- = 104 dual "
     "negative pivots confined to the lower dual half, and w9 "
     "SATURATES the bound exactly: max neg dual pivot 182 == "
     "S - 1 - N_w) + DUAL_RESTATEMENT (dual negative-mass "
     "fraction 0.424, no isolated carrier) + "
     "MOMENT_PERTURBATION_RANGE(X_weyl = 2 of N_w = 184: the "
     "natural mu/nu Weyl minor-perturbation argument carries 1.1 "
     "percent of the way; PAIRCORR detector: Weyl X = 1/2/6/2 on "
     "EPST/SCR/SMOOTH/HL2 (SMOOTH exceeds MAIN) and the raw-"
     "moment rest zone are BOTH world-blind -- neither carries "
     "the localization); amendments a1-a3 disclosed (Richardson "
     "FD protocol, dual-mirror slice bug fix, wording); NO "
     "localization claim from the census; consumed by v962 "
     "(wave 5, embedded byte-exact, smoke stage)", True),
    (f"{EXP}/halffilling_pinning_probe.py", "sealed_probe",
     "r281",
     ["PRIME.PORT.WALL.HALFFILLING_PINNING.01"],
     "halffilling-pinning round (28/28): UPPER_PINNING_MEASURED (C"
     " = 5 over the 42-rung census; the v956 boundary statement "
     "'the wall dies immediately past the cap' is carried by TWO "
     "computational paths + an mp dps-40 ward on five windows = "
     "MEASUREMENT, not a symbolic theorem; the only PROVEN upper "
     "side is the pigeonhole ceiling minC <= S - S_- from the r279"
     " budget theorem -- w9: 184 <= 263 = p, 79 degrees above N_w;"
     " world-blind O(1) upper pinning REFUTED exactly: the "
     "single-tiny-negative toy measure has minC = S-1, offset = "
     "N_w - 2, unbounded in S -- any O(1) upper pinning must "
     "consume the comb structure) + OFFSET_UNSTRUCTURED (six "
     "sealed source-pure tail-coupling candidates vs offset over "
     "the 42 rungs: max |spearman| = 0.273 << 0.75 -- no "
     "predictive offset formula in this census; honest negative) +"
     " RELAY_CONDITION_FOUND(B_n = [beta_n > 0], RESTATEMENT: no "
     "h-blind margin proxy predicts the flip location on all five "
     "worlds; the r280 anti-alignment cos_w = -0.97 is the h-sign "
     "flip of a raw-gradient lockstep, cos_raw = +0.96..+0.98 at "
     "every crossing; all four dead-world crossings type CREEPING "
     "-- crossing type does not separate the worlds) + "
     "PINNING_DECOMPOSED (the program finding: COUNTING THEOREM "
     "exact -- the free pivots are exactly h_0..h_{N_w-1}, h_{N_w}"
     " is the first forced pivot (toys in rationals, w9 mp dps-200"
     " L-recurrence residual 7.6e-157, forced-fraction profile "
     "exact), so 'why half-filling' = THE FREE MOMENTS END THERE; "
     "the lower side minC >= N_w == 'all free pivots positive' "
     "stays the open center).  Regressions: r280 census "
     "distribution + anchors + flips + cosines exact; v956 quotes "
     "byte-gated.  NO RH claim.; consumed by v962 (wave 5, "
     "embedded byte-exact, smoke stage)",
     True),
    (f"{EXP}/representation_contest_probe.py", "sealed_probe",
     "r282",
     ["PRIME.PORT.REPRESENTATION.CONTEST.01"],
     "representation-contest round (30/30): the north-star "
     "question 'what IS h_n' run as a four-language contest with "
     "one brutal gate (structural forcing of MAIN positivity + "
     "provable rejection of an early-dying control; calibrated on "
     "the positive-class toy POSLIKE). K1 RHP/MONODROMY typed "
     "REPRESENTATION_RESTATEMENT (Pontryagin split h_n = P_n - N_n"
     " exact at every degree on all four toys; the "
     "|m|^2/sum-of-squares structure exists iff the negative "
     "register is empty -- true exactly on POSLIKE, false on every"
     " signed toy incl. MAINLIKE with N/P = 0.538 at the last free"
     " degree; w9 negfrac 0.498 at the wall; the surviving clause "
     "N < P is h > 0 verbatim) + K2 HODGE/INTERSECTION typed "
     "REPRESENTATION_NOT_CONSTRUCTIBLE (the Kasteleyn/orientation "
     "question answered by FULL EXHAUSTION: the coherent linear "
     "orientation class over the Cauchy-Binet configuration space "
     "is EXACTLY {+-sign w} on all four toys (2^S gauges each), "
     "the quadratic/edge class collapses at the n = 3 cells; the "
     "unique surviving gauge changes the value by exactly 2 x "
     "negative mass (exact defect at n = 1; w9 0.1152, S_- = 104) "
     "-- coherence and value preservation are simultaneously "
     "possible iff S_- = 0: the orientation exists exactly for the"
     " positive class, and MAIN is signed; the r261 no-go upgraded"
     " from pairing failure to class obstruction) + K3 DE "
     "BRANGES/CANONICAL SYSTEM typed REPRESENTATION_RESTATEMENT "
     "(Hamiltonian PSD through degree n <=> h_0..h_n > 0 exact "
     "bookkeeping; the only independent shadow, R2 interlacing, "
     "lags +1 on every flipping world -- toys minC + 1 exact, "
     "controls 26/22/28/26 == flip + 1, w9 SAFE through n = 183: "
     "the shadow detects, it never forbids; the |beta| repair "
     "changes the measure at k = 4 exactly) + K4 OPERATOR/DUAL "
     "PAIR typed REPRESENTATION_WORLD_BLIND (h_n h#_{S-1-n} == 1 "
     "exact on all four toys => sign sync holds IDENTICALLY in "
     "every world through every flip -- the v956 dual pair adds no"
     " second constraint; global positive structures are killed by"
     " the r279 budget theorem #(h<0) = S_-) + CONTEST_ALL_DEAD "
     "(the program finding: the pivot positivity has NO "
     "representation in the four classical languages at this "
     "resolution -- each language forces positivity exactly for "
     "the positive measure class, and MAIN is signed; the lower "
     "side stays the open center).  Detector: D1 negfrac / D2 "
     "defect / D4 sync WORLD_BLIND, D3 fire/N separates but is "
     "DECLARED h-reading (decoration).  Amendment a1 (smoke stage,"
     " gate-side m4 loudness form) disclosed.  NO RH claim.; "
     "consumed by v963 (wave 6, embedded byte-exact, smoke stage)",
     True),
    (f"{EXP}/fullsource_quasidefiniteness_probe.py", "sealed_probe",
     "r283",
     ["PRIME.PORT.RHP.FULLSOURCE.QUASIDEFINITENESS.01"],
     "fullsource-quasidefiniteness round (26/26): the big RHP "
     "attack (full source -> monodromy property -> "
     "quasi-definiteness, three sealed candidates). VERDICT "
     "MECHANISM_PARTIAL(A2 monodromy factorization over the "
     "positive/negative part: in the mu-orthonormal frame the full"
     " Hankel form is I - M_n with the node Gram E_n = B_n B_n^T "
     "EXACTLY the v881/r224 IIKS kernel (mu-CD kernel dressed by "
     "the nu weights); steps s1-s4 EXACT -- split, unit-triangular"
     " congruence minor_k(D_mu - G) == D_k(mutilde) (rationals on "
     "five toys, slogdet dev 3.4e-9 on w9), frame contraction h>0 "
     "window <=> lambda_max(E_m) < 1, monotone rank-one loading =>"
     " lambda_max crosses 1 EXACTLY at minC+1; WORLD TEST exact on"
     " all seven worlds: MAIN w9 185 == minC+1, w13 171 == N+3, "
     "kz15 205 == N+2, controls EPST/SCR/SMOOTH/HL2 26/22/28/26 =="
     " flip+1 == the r277 R2 fire degrees -- the contraction IS "
     "the split-source form of the r277 reality/interlacing "
     "condition; s5 capacity: the pigeonhole ceiling minC <= S_+ "
     "RE-PROVED in the frame (null-polynomial witness, exact on "
     "five toys), capacity-as-COUNTING REFUTED exactly (rank-one "
     "pair eps/big: same rank, same support, opposite window fate;"
     " the deciding value is metric 1.25e-13 vs 1.25e+2); MISSING "
     "LEMMA L* formal: for MAIN every real polynomial p with deg p"
     " < N_w has int p^2 dnu < int p^2 dmu, i.e. "
     "lambda_max(E_{N_w}) < 1 (measured margin 1 - rho_184 = "
     "1.68e-4, top eigs 0.99983/0.99874/0.99597 -- broadly tight, "
     "not single-mode); given s1-s4 this ONE scalar <=> full "
     "free-window quasi-definiteness) + A3_COLLAPSES_TO_A2 "
     "(Fredholm zero-freeness of det(I - s E_{N_w}) on (0,1] IS "
     "the contraction scalar: s* = 1/rho dictionary + "
     "degree-resolved slogdet co-location 185/26) + A1_DEAD "
     "(Schwarz reflection symmetry of the real jump data holds "
     "identically on every world = WORLD_BLIND, cannot carry; the "
     "classical additional condition -- residue "
     "positivity/Herglotz -- fails on MAIN ITSELF (S_- = 104 "
     "negative residues); the windowed eta* anatomy is world-blind"
     " by the sealed distance rule; the additional condition "
     "collapses onto L*) + DETECTOR_LEDGER (K_A2 rho_20 "
     "WORLD_BLIND -- honest: the early-degree contraction VALUE "
     "does not separate; K_A3 1/rho_40 MAIN_SEPARATING; "
     "triangle/Gershgorin bound dead at n = 21 << 185; "
     "target-inverse and posthoc-window mutants FLAGGED; full "
     "source consumed, no truncation).  The RHP language holds "
     "exactly ONE invariant on this question and it is the "
     "contraction scalar; the lower side (L* itself) stays the "
     "open center.  NO RH claim.; consumed by v963 (wave 6, "
     "embedded byte-exact, smoke stage)",
     True),
    (f"{EXP}/lstar_two_measure_probe.py", "sealed_probe",
     "r284",
     ["PRIME.PORT.LSTAR.TWO_MEASURE_ANATOMY.01"],
     "L* two-measure anatomy round (30/30): the r283 missing lemma"
     " L* read as a TWO-MEASURE SUBORDINATION problem (nu 104 "
     "atoms strictly mu-subordinate, 263 atoms, on deg < N_w).  "
     "VERDICT SHIELDING_BLIND (the source-pure pairing geometry is"
     " a builder-channel property: on w9 EVERY nu atom sits one "
     "fold step from a mu atom, interlacing 100/104 all-mu-"
     "neighbored, and the four dead controls mirror the same "
     "near-perfect pairing -- all four sealed geometry statistics "
     "WORLD_BLIND by the r281 distance rule) + NYQUIST_ORDERING("
     "MAX sp +0.998 over 46 worlds, HONEST: carried by the trivial"
     " N_w scaling of the ladder; the sealed band |log2(n_pred/"
     "crossing)| <= 1 with n_pred = L/(4s) FAILS on all four "
     "controls, devs 1.48..3.06 octaves -- the controls die at "
     "22..28 DESPITE near-perfect pairing: the support-geometry "
     "form of the resolution/shielding hypothesis is refuted where"
     " it matters) + CHRISTOFFEL_RANGE(the round's core: exact "
     "sandwich max_k v_k K_n(y_k) <= lambda_max(E_n) <= trace(E_n)"
     " gated everywhere; on w9 the trace/triangle side crosses at "
     "n_CS = 10 while the pure single-atom Christoffel bound "
     "crosses at n_DIAG = 187, only TWO degrees above the true "
     "crossing 185: MAIN's wall crossing is a NEAR-SINGLE-ATOM "
     "Christoffel event (coherent assist gain rho/maxdiag = 1.031 "
     "at N_w; aggregate slack trace/rho = 50.2 -- destructive nu "
     "coherence, the cancellation story in the two-measure "
     "coordinate), while on ALL FOUR controls n_DIAG is never "
     "reached by their crossing: control death is a COLLECTIVE "
     "mode) + EXTREMAL_ANATOMY(the near-breaking direction is "
     "99.7 percent TWO background atoms at the shallow-u hull edge"
     " below the first prime (folds 2 and 4, u = 0.03/0.06, "
     "x -> +1, the two weakest-shielded nu atoms, PR 1.89), NOT "
     "the small primes -- SMALLP_DEPLETED, the r278 chain-"
     "sensitivity profile and the extremal near-null direction are"
     " different coordinates; controls at crossing are MORE "
     "diffuse (PR 4.97..9.72); top-3 eigenvectors ONE_BAND, three "
     "oscillatory modes of one edge band) + DETECTOR_LEDGER(all "
     "seven sealed stats WORLD_BLIND -- the separation lives in "
     "the range structure, not in any sealed scalar).  mp dps-60 "
     "ward: rho_184 = 0.99983248 < 1 < 1.00003660 = rho_185.  "
     "Regressions: r283 spectral records, r281 42-rung census "
     "distribution + anchors, control flips 25/21/27/25.  NO RH "
     "claim.; consumed by v963 (wave 6, embedded byte-exact, "
     "smoke stage)",
     True),
    (f"{EXP}/christoffel_decomposition_probe.py", "sealed_probe",
     "r285",
     ["PRIME.PORT.LSTAR.CHRISTOFFEL_DECOMPOSITION.01"],
     _R285_STATUS + "; consumed by v963 (wave 6, embedded "
     "byte-exact, smoke stage)",
     True),
    (f"{EXP}/lstar_margin_scaling_probe.py", "sealed_probe",
     "r286",
     ["PRIME.PORT.LSTAR.MARGIN_SCALING.01"],
     _R286_STATUS + "; consumed by v964 (wave 7, embedded byte-exact, smoke stage)",
     True),
    (f"{EXP}/l2_deterministic_cancellation_probe.py", "sealed_probe",
     "r287",
     ["PRIME.PORT.L2.DETERMINISTIC_CANCELLATION.01"],
     _L2_STATUS + "; consumed by v964 (wave 7, embedded byte-exact, smoke stage)",
     True),
    (f"{EXP}/destructive_coherence_probe.py", "sealed_probe",
     "r288",
     ["PRIME.PORT.LSTAR.DESTRUCTIVE_COHERENCE.01"],
     _R288_STATUS + "; consumed by v964 (wave 7, embedded byte-exact, smoke stage)",
     True),
    (f"{EXP}/arch_kernel_diophantine_probe.py", "sealed_probe",
     "r289",
     ["PRIME.PORT.LSTAR.ARCH_KERNEL_DIOPHANTINE_ANATOMY.01"],
     _R289_STATUS + "; consumed by v964 (wave 7, embedded byte-exact, smoke stage)",
     True),
    (f"{EXP}/profile_functional_probe.py", "sealed_probe",
     "r290",
     ["PRIME.PORT.LSTAR.PROFILE_FUNCTIONAL.01"],
     _R290_STATUS + "; consumed by v965 (wave 8, embedded byte-exact, smoke stage)",
     True),
    (f"{EXP}/ridge_anatomy_probe.py", "sealed_probe",
     "r291",
     ["PRIME.PORT.LSTAR.RIDGE_ANATOMY.01"],
     _R291_STATUS + "; consumed by v965 (wave 8, embedded byte-exact, smoke stage)",
     True),
    (f"{EXP}/curvature_form_probe.py", "sealed_probe",
     "r292",
     ["PRIME.PORT.LSTAR.CURVATURE_FORM.01"],
     _R292_STATUS + "; consumed by v965 (wave 8, embedded byte-exact, smoke stage)",
     True),
    (f"{EXP}/metric_reconciliation_probe.py", "sealed_probe",
     "r293",
     ["PRIME.PORT.LSTAR.METRIC_RECONCILIATION.01"],
     _R293_STATUS + "; consumed by v965 (wave 8, embedded byte-exact, smoke stage)",
     True),
    (f"{EXP}/f10_stability_probe.py", "sealed_probe",
     "r294",
     ["PRIME.PORT.LSTAR.F10_STABILITY.01"],
     _R294_STATUS + "; consumed by v965 (wave 8, embedded byte-exact, smoke stage)",
     True),
    (f"{EXP}/f10_sp_hardening_probe.py", "sealed_probe",
     "r295",
     ["PRIME.PORT.LSTAR.F10_SP_HARDENING.01"],
     _R295_STATUS + "; consumed by v965 (wave 8, embedded byte-exact, smoke stage)",
     True),
    (f"{EXP}/dens_identity_probe.py", "sealed_probe",
     "r296",
     ["PRIME.PORT.LSTAR.DENS_IDENTITY.01"],
     _R296_STATUS
     + "; consumed by v966 (wave 9, embedded "
     "byte-exact, smoke stage)",
     True),
    (f"{EXP}/vdc_chain_provenance_probe.py", "sealed_probe",
     "r297",
     ["PRIME.PORT.L2.VDC_CHAIN_PROVENANCE.01",
      "PRIME.PORT.L2.VDC_LEMMA.01"],
     _R297_STATUS
     + "; consumed by v966 (wave 9, embedded "
     "byte-exact, smoke stage)",
     True),
    (f"{EXP}/window_border_transfer_probe.py", "sealed_probe",
     "r298",
     ["PRIME.PORT.L2.WINDOW_BORDER_TRANSFER.01",
      "PRIME.PORT.L2.VDC_LEMMA.01"],
     _R298_STATUS
     + "; consumed by v966 (wave 9, embedded "
     "byte-exact, smoke stage)",
     True),
    (f"{EXP}/fejer_decay_probe.py", "sealed_probe",
     "r299",
     ["PRIME.PORT.L2.FEJER_DECAY.01",
      "PRIME.PORT.L2.VDC_LEMMA.01"],
     _R299_STATUS
     + "; consumed by v966 (wave 9, embedded "
     "byte-exact, smoke stage)",
     True),
    (f"{EXP}/diag_target_probe.py", "sealed_probe",
     "r300",
     ["PRIME.PORT.L2.DIAG_TARGET.01",
      "PRIME.PORT.L2.VDC_LEMMA.01"],
     _R300_STATUS
     + "; consumed by v966 (wave 9, embedded "
     "byte-exact, smoke stage)",
     True),
    (f"{EXP}/neff_target_probe.py", "sealed_probe",
     "r301",
     ["PRIME.PORT.L2.NEFF_TARGET.01",
      "PRIME.L2.NEFF_TARGET.01",
      "PRIME.PORT.L2.VDC_LEMMA.01"],
     _R301_STATUS,
     True),
    (f"{EXP}/unif_target_probe.py", "sealed_probe",
     "r302",
     ["PRIME.PORT.L2.UNIF_TARGET.01",
      "PRIME.L2.NEFF_TARGET.01",
      "PRIME.PORT.L2.VDC_LEMMA.01"],
     _R302_STATUS,
     True),
    (f"{EXP}/atom_target_probe.py", "sealed_probe",
     "r303",
     ["PRIME.PORT.L2.ATOM_TARGET.01",
      "PRIME.L2.NEFF_TARGET.01",
      "PRIME.PORT.L2.VDC_LEMMA.01"],
     _R303_STATUS,
     True),
    (f"{EXP}/shortrange_law_probe.py", "sealed_probe",
     "r304",
     ["PRIME.PORT.L2.SHORTRANGE_LAW.01",
      "PRIME.L2.NEFF_TARGET.01",
      "PRIME.PORT.L2.VDC_LEMMA.01"],
     _R304_STATUS,
     True),
    (f"{EXP}/renyi3_probe.py", "sealed_probe",
     "r306",
     [],
     _R306_STATUS,
     True),
    (f"{EXP}/fixed_head_probe.py", "sealed_probe",
     "r307",
     [],
     _R307_STATUS,
     True),
    (f"{EXP}/block_green_probe.py", "sealed_probe",
     "r308",
     [],
     _R308_STATUS,
     True),
    (f"{EXP}/paired_cone_probe.py", "sealed_probe",
     "r309",
     [],
     _R309_STATUS,
     True),
    (f"{EXP}/blockgreen_nontriviality_probe.py", "sealed_probe",
     "r311",
     [],
     _R311_STATUS,
     True),
    (f"{EXP}/blockgreen_membership_probe.py", "sealed_probe",
     "r312",
     [],
     _R312_STATUS,
     True),
    (f"{EXP}/renyi3_proof_fork_probe.py", "sealed_probe",
     "r313",
     [],
     _R313_STATUS,
     True),
    (f"{EXP}/signed_cubic_flux_probe.py", "sealed_probe",
     "r314",
     [],
     _R314_STATUS,
     True),
    (f"{EXP}/phi3_functional_probe.py", "sealed_probe",
     "r315",
     [],
     _R315_STATUS,
     True),
    (f"{EXP}/two_regime_bound_probe.py", "sealed_probe",
     "r316",
     [],
     _R316_STATUS,
     True),
    (f"{EXP}/exception_families_probe.py", "sealed_probe",
     "r317",
     [],
     _R317_STATUS,
     True),
    (f"{EXP}/indefinite_fork_probe.py", "sealed_probe",
     "r318",
     [],
     _R318_STATUS,
     True),
    (f"{EXP}/continuous_coordinate_probe.py", "sealed_probe",
     "r321",
     [],
     _R321_STATUS,
     True),
    (f"{EXP}/antiphase_sign_law_probe.py", "sealed_probe",
     "r322",
     [],
     _R322_STATUS,
     True),
    (f"{EXP}/fa_provenance_probe.py", "sealed_probe",
     "r324-pre",
     [],
     _R324PRE_STATUS,
     True),
    (f"{EXP}/qmax_m2_origin_probe.py", "sealed_probe",
     "r324",
     [],
     _R324_STATUS,
     True),
    (f"{EXP}/extraction_order_probe.py", "sealed_probe",
     "r325",
     [],
     _R325_STATUS,
     True),
    (f"{EXP}/group_mass_cap_probe.py", "sealed_probe",
     "r327",
     [],
     _R327_STATUS,
     True),
    # -- the standalone problem statement (rh/problem/) --
    ("rh/problem/lstar_problem.tex", "problem_statement",
     "r283 companion",
     ["PRIME.PORT.RHP.FULLSOURCE.QUASIDEFINITENESS.01"],
     "the missing lemma L* as a fully standalone open problem for "
     "external mathematicians (no project vocabulary): mu/nu defined "
     "from first principles (von Mangoldt comb, tent sampling, "
     "archimedean lag function, circulant spectral density, "
     "folding), the S = 367 instance tabulated, Hankel / CD-kernel / "
     "sampling equivalents, neutral numerical record (42/42 hold; "
     "margins 1.68e-4 down to 1.42e-7, shrinking with S), "
     "moment-counting meaning of the degree bound (S+1)/2; "
     "the 2026-08-25 update paragraph (57/57 instances positive, "
     "flattening power law, the rational-twin fact: the problem is "
     "NOT diophantine, the metric threshold 1e-3..3e-3 gap); "
     "since 2026-08-27 the freeze memo section 'The frozen state "
     "(August 2026): one object' (reviewer decision: the L* lane "
     "FROZEN as specialist problem): the fully typed chain "
     "(dictionary weights / candidate deficits / theorem-grade "
     "pinning / computability closure), the one object (phi_D = "
     "log(pq) + (17/12) ln N and phi_K = log(c^2) + (17/12) ln N "
     "are one shared wander, corr 0.999998, 0.88 nats each vs "
     "0.0017 nats difference -- factor ~500 -- whose difference IS "
     "pointwise the log reserve with rho_r = 2.624 its decay; "
     "two-block interference: pair geometry +0.887 at 2.4x "
     "amplitude against the anti-correlated weight wander -0.72), "
     "the three specialist questions in final form (q1 backward-CS "
     "= the cancellation ratio, q2 sub-classical Christoffel "
     "growth with candidates unresolvable below 10^5.4+, q3 the "
     "resolution paradox), and the honest census frame (85 windows "
     "to N_w 7942, all margins positive, pool exhausted at "
     "10^3.90); explicitly NO claim -- L* stays [O]", True),
    ("rh/problem/verify_lstar_instance.py", "problem_check",
     "r283 companion",
     ["PRIME.PORT.RHP.FULLSOURCE.QUASIDEFINITENESS.01"],
     "standalone rebuild of the document's mu/nu construction "
     "(numpy/mpmath only) cross-checked against the repository "
     "builders: atoms/weights agree to 4e-16 rel, lambda_max(E_184) "
     "= 0.99983248 / lambda_max(E_185) = 1.00003660 reproduced with "
     "dps-60 mpmath ward, inequality gated on monomials + 500 "
     "random polynomials; 12/12 gates, final line LSTAR INSTANCE "
     "VERIFIED; optional --ladder sweeps all 42 admissible windows "
     "(all margins positive, min 1.42e-7 at kz64)", True),
    ("rh/problem/lstar_problem.pdf", "problem_statement",
     "r283 companion",
     [],
     "compiled PDF of lstar_problem.tex (recompiled artefact, "
     "registered unpinned)", False),
    # -- frozen libraries the modules embed byte-exact --
    (f"{EXP}/tau_symbolic_probe.py", "frozen_library",
     "r224",
     ["PRIME.PORT.TAU.SYMBOLIC.01"],
     "embedded in v955", True),
    (f"{EXP}/lax_conditioned_probe.py", "frozen_library",
     "r225",
     ["PRIME.PORT.LAX2.CONDITIONED.01"],
     "embedded in v955", True),
    (f"{EXP}/hirota_sign_probe.py", "frozen_library",
     "r226",
     ["PRIME.PORT.HIROTA.SIGN.01"],
     "embedded in v955; library for v956/v959", True),
    (f"{EXP}/fermiedge_classify_probe.py", "frozen_library",
     "r227",
     [],
     "builder library for v956 (classification round)", True),
    (f"{EXP}/holedual_probe.py", "frozen_library",
     "r228",
     ["PRIME.PORT.SIGNEDMOMENT.HOLEDUAL.01"],
     "embedded in v956", True),
    (f"{EXP}/pontryagin_maxpos_probe.py", "frozen_library",
     "r229",
     ["PRIME.PORT.PONTRYAGIN.MAXPOS.01"],
     "embedded in v956", True),
    (f"{EXP}/jfraction_probe.py", "frozen_library",
     "r230",
     ["PRIME.PORT.FREEMOMENT.JFRACTION.01"],
     "embedded in v956", True),
    (f"{EXP}/rhp_midpoint_probe.py", "frozen_library",
     "r231",
     ["PRIME.PORT.RHP.FREEMOMENT.MIDPOINT.01"],
     "embedded in v956", True),
    (f"{EXP}/szego_equilibrium_probe.py", "frozen_library",
     "r232a",
     [],
     "frozen library (equilibrium QP) for v959 probes", True),
    (f"{EXP}/principal_bessel_probe.py", "frozen_library",
     "r243",
     [],
     "frozen library (budget definition B = S_{N-2} + 5/7)", True),
    (f"{EXP}/bordered_hankel_probe.py", "frozen_library",
     "r244",
     ["PRIME.PORT.RHP.BORDERED.READOUT.01"],
     "embedded in v958", True),
    (f"{EXP}/bordered_finite_rank_probe.py", "frozen_library",
     "r245",
     [],
     "embedded in v958 (Schlesinger rank-1 insertion)", True),
    (f"{EXP}/border_centering_probe.py", "frozen_library",
     "r248",
     [],
     "embedded in v958 (centering congruence)", True),
    # -- living documents (registered, unpinned) --
    (f"{VER}/status_ledger.csv", "ledger",
     "-",
     ["(all PRIME.* rows; the ledger wins on any status disagreement)"],
     "living document -- drift is expected, reported as INFO", False),
    ("tfpt_prime_front.tex", "program_prose",
     "-",
     [],
     "living document -- the program prose of the prime front", False),
    ("experiments/next.txt", "campaign_notes",
     "notes DLXXIX-DXCVII",
     [],
     "living document -- German campaign notes (newest at bottom)", False),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    entries = []
    missing = []
    for rel, role, rnd, ids, status, pin in ENTRIES:
        p = os.path.join(REPO, rel)
        if not os.path.isfile(p):
            missing.append(rel)
            continue
        entries.append({
            "path": rel,
            "role": role,
            "round": rnd,
            "ledger_ids": ids,
            "status": status,
            "pin": pin,
            "sha256": sha256(p),
            "bytes": os.path.getsize(p),
        })
    doc = {
        "workspace": "rh/",
        "generated": date.today().isoformat(),
        "claim_boundary": ("Research documentation.  NOT evidence for or "
                           "against the Riemann Hypothesis in either "
                           "direction.  NO RH CLAIM.  Mincut base 4 / "
                           "refined 5.  The status ledger wins on any "
                           "disagreement."),
        "note": ("rh/ references these files, it never copies them.  "
                 "pin=true entries are frozen: a SHA-256 mismatch is a "
                 "FAIL in rh/verification/run_rh.py.  pin=false entries "
                 "are living documents: drift is reported as INFO."),
        "entries": entries,
    }
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print("wrote %s  (%d entries, %d pinned)" % (
        OUT, len(entries), sum(1 for e in entries if e["pin"])))
    if missing:
        print("MISSING (not registered):")
        for m in missing:
            print("  -", m)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
