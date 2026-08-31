#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_rh.py -- the RH-specific verification suite of the rh/ workspace.

Runs, in order:
  (1) INTEGRITY  -- SHA-256 of every pinned rh/INVENTORY.json entry
                    (pinned drift = FAIL; unpinned living documents = INFO),
                    (2) PROBES     -- the sealed campaign probes r250-r464+r467+r468 from
                    experiments/tfpt-discovery/ in --smoke mode,
  (3) MODULES    -- the twenty RH verification modules v955/v956/v958/
                    v959/v960/v961/v962/v963/v964/v965/v966/v967/v968/
                    v969/v970 plus the wave-14 set v978/v979/v980/v981/
                    v982, executed BY PATH from verification/ (never
                    copied; v959 ~3.5 min, v968 ~50 s, the rest seconds:
                    wave-4..13 modules embed their probes in the sealed
                    --smoke stage; the wave-14 modules re-derive their
                    exact layers from scratch and gate frozen record
                    aggregates) -- skipped under --fast,
  (4) LEAN       -- `lake build` in rh/lean/ if a Lean toolchain is present.

Output is house-style: [PASS]/[FAIL] per item, final line
  RH SUITE: ALL CHECKS PASSED    or    RH SUITE: FAILURES n

Claim boundary: research documentation.  NOT evidence for or against the
Riemann Hypothesis in either direction.  NO RH CLAIM.

Usage (from the repo root):
  python rh/verification/run_rh.py            # full run
  python rh/verification/run_rh.py --fast     # integrity + probes + lean
  python rh/verification/run_rh.py --skip-lean
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
INVENTORY = os.path.join(REPO, "rh", "INVENTORY.json")
LEAN_DIR = os.path.join(REPO, "rh", "lean")
VENV_PY = os.path.join(REPO, "experiments", "tfpt-discovery", ".venv",
                       "bin", "python")

# the sealed campaign probe list r250-r464+r467+r468 (order = round order); every
# entry supports --smoke.  This list is frozen with the wave; extend it in
# the same change that extends INVENTORY.json.
# Lean-only rounds (no sealed probe): r305, r310, r310b, r320, C1,
# r362, r373, r376, r380, r384, r397, r406, r412, r426, r430, r434.
PROBES = [
    ("r250", "centered_basefiber_probe.py"),
    ("r251a", "corner_provenance_probe.py"),
    ("r251", "targetreadout_error_probe.py"),
    ("r252", "base_gauge_constant_probe.py"),
    ("r253", "schlesinger_pairing_probe.py"),
    ("r254", "offdiag_gram_probe.py"),
    ("r255", "nodebands_base_probe.py"),
    ("r256", "baseborder_factorial_probe.py"),
    ("r257", "coupledtau_probe.py"),
    ("r258", "budget_anatomy_probe.py"),
    ("r259", "parametrix_pass_probe.py"),
    ("r260", "terminal_crossratio_probe.py"),
    ("r261", "prefix_resummation_probe.py"),
    ("r262", "terminal_triangle_probe.py"),
    ("r263", "cancellation_adjudication_probe.py"),
    ("r264", "quenched_opening_probe.py"),
    ("r265", "s_monotonicity_probe.py"),
    ("r266", "border_resolvent_identity_probe.py"),
    ("r267", "ranktrace_adjudication_probe.py"),
    ("r268", "drive_local_asymptotics_probe.py"),
    ("r269", "phase_bulk_bound_probe.py"),
    ("r270", "kz15_boss_probe.py"),
    ("r271", "universal_pair_theorem_probe.py"),
    ("r272", "l2_scaling_anatomy_probe.py"),
    ("r273", "euler_mechanism_probe.py"),
    ("r274", "wronskian_dictionary_probe.py"),
    ("r275", "kyp_memory_probe.py"),
    ("r276", "minimal_firewall_probe.py"),
    ("r277", "maslov_census_probe.py"),
    ("r278", "metric_stability_probe.py"),
    ("r279", "oriented_theorem_probe.py"),
    ("r280", "budget_localization_probe.py"),
    ("r281", "halffilling_pinning_probe.py"),
    ("r282", "representation_contest_probe.py"),
    ("r283", "fullsource_quasidefiniteness_probe.py"),
    ("r284", "lstar_two_measure_probe.py"),
    ("r285", "christoffel_decomposition_probe.py"),
    ("r286", "lstar_margin_scaling_probe.py"),
    ("r287", "l2_deterministic_cancellation_probe.py"),
    ("r288", "destructive_coherence_probe.py"),
    ("r289", "arch_kernel_diophantine_probe.py"),
    ("r290", "profile_functional_probe.py"),
    ("r291", "ridge_anatomy_probe.py"),
    ("r292", "curvature_form_probe.py"),
    ("r293", "metric_reconciliation_probe.py"),
    ("r294", "f10_stability_probe.py"),
    ("r295", "f10_sp_hardening_probe.py"),
    ("r296", "dens_identity_probe.py"),
    ("r297", "vdc_chain_provenance_probe.py"),
    ("r298", "window_border_transfer_probe.py"),
    ("r299", "fejer_decay_probe.py"),
    ("r300", "diag_target_probe.py"),
    ("r301", "neff_target_probe.py"),
    ("r302", "unif_target_probe.py"),
    ("r303", "atom_target_probe.py"),
    ("r304", "shortrange_law_probe.py"),
    ("r306", "renyi3_probe.py"),
    ("r307", "fixed_head_probe.py"),
    ("r308", "block_green_probe.py"),
    ("r309", "paired_cone_probe.py"),
    ("r311", "blockgreen_nontriviality_probe.py"),
    ("r312", "blockgreen_membership_probe.py"),
    ("r313", "renyi3_proof_fork_probe.py"),
    ("r314", "signed_cubic_flux_probe.py"),
    ("r315", "phi3_functional_probe.py"),
    ("r316", "two_regime_bound_probe.py"),
    ("r317", "exception_families_probe.py"),
    ("r318", "indefinite_fork_probe.py"),
    ("r321", "continuous_coordinate_probe.py"),
    ("r322", "antiphase_sign_law_probe.py"),
    ("r324pre", "fa_provenance_probe.py"),
    ("r324", "qmax_m2_origin_probe.py"),
    ("r325", "extraction_order_probe.py"),
    ("r327", "group_mass_cap_probe.py"),
    ("r329", "ext3_fresh_anchors_probe.py"),
    ("r330", "dirichlet_secondworld_probe.py"),
    ("r331", "twin_resolution_probe.py"),
    ("r333", "companion_orbit_packing_probe.py"),
    ("r334", "fold_capacity_probe.py"),
    ("r335", "edge_packing_dichotomy_probe.py"),
    ("r336", "lstar_parity_section_probe.py"),
    ("r337", "fold_martingale_probe.py"),
    ("r339", "fold_density_dictionary_probe.py"),
    ("r340", "cauchybinet_hall_probe.py"),
    ("r341", "fold_bellman_reverse_holder_probe.py"),
    ("r342", "pair_extremal_probe.py"),
    ("r343", "pair_coupling_probe.py"),
    ("r344", "fold_two_scale_balance_probe.py"),
    ("r345", "gap_ratio_primary_probe.py"),
    ("r346", "fold_cover_canonization_probe.py"),
    ("r347", "delta_alpha_closure_probe.py"),
    ("r348", "delta_source_anatomy_probe.py"),
    ("r349", "thirdarm_spike_law_probe.py"),
    ("r350", "alpha_source_anatomy_probe.py"),
    ("r351", "qmax_growth_law_probe.py"),
    ("r352", "rhor_source_anatomy_probe.py"),
    ("r353", "second_family_erosion_probe.py"),
    ("r354", "phi_wander_anatomy_probe.py"),
    ("r355", "k2_source_formula_probe.py"),
    ("r356", "borodin_dual_hole_probe.py"),
    ("r357", "dirichlet_matched_frame_probe.py"),
    ("r358", "local_gap_carleson_probe.py"),
    ("r359", "schur_wronskian_dual_probe.py"),
    ("r360", "critical_saturation_probe.py"),
    ("r361", "mean_sieve_floor_probe.py"),
    ("r362", "augmented_borodin_duality_probe.py"),
    ("r363", "canonical_sturm_induction_probe.py"),
    ("r366", "edge_gap_ms_probe.py"),
    ("r367", "final_two_rank_inertia_probe.py"),
    ("r368", "weighted_l2_t1_probe.py"),
    ("r369", "mixed_haynsworth_probe.py"),
    ("r370", "matrix_weyl_index_probe.py"),
    ("r371", "compound_cd_wedge_probe.py"),
    ("r372", "source_prufer_one_defect_probe.py"),
    ("r377", "postcap_pivots_probe.py"),
    ("r381", "g_eps_lemma_probe.py"),
    ("r382", "pivot_entry_lemma_probe.py"),
    ("r383", "compose_premises_probe.py"),
    ("r385", "christoffel_quiet_probe.py"),
    ("r386", "compose_premises2_probe.py"),
    ("r387", "coherence_assist_probe.py"),
    ("r388", "delta_deformation_probe.py"),
    ("r389", "weyl_energy_probe.py"),
    ("r390", "g_eps_mu_probe.py"),
    ("r391", "construction_pure_rl_probe.py"),
    ("r392", "deletion_transform_probe.py"),
    ("r393", "tau_field_probe.py"),
    ("r394", "sign_schur_probe.py"),
    ("r395", "three_gap_mask_probe.py"),
    ("r396", "isolation_lemma_probe.py"),
    ("r398", "high_moment_inertia_probe.py"),
    ("r399", "source_weyl_energy_probe.py"),
    ("r400", "bulk_one_defect_probe.py"),
    ("r401", "edge_signature_probe.py"),
    ("r403", "p1_construction_probe.py"),
    ("r404", "one_defect_gram_probe.py"),
    ("r405", "edge_contractive_lift_probe.py"),
    ("r407", "dual_intertwiner_probe.py"),
    ("r408", "c_threshold_probe.py"),
    ("r409", "borodin_birkhoff_intertwiner_probe.py"),
    ("r410", "hole_nyquist_probe.py"),
    ("r411", "threshold_identity_probe.py"),
    ("r413", "hole_top_mode_probe.py"),
    ("r415", "top_mode_edge_probe.py"),
    ("r416", "debranges_index_probe.py"),
    ("r417", "source_sch_sign_probe.py"),
    ("r418", "phi_bb_sign_probe.py"),
    ("r419", "vacuous_overflow_probe.py"),
    ("r420", "cj_sigma_probe.py"),
    ("r421", "reserve_limit_probe.py"),
    ("r422", "sigma_limit_probe.py"),
    ("r423", "den_limit_probe.py"),
    ("r424", "gamma_chain_probe.py"),
    ("r425", "cross_chain_gamma_probe.py"),
    ("r427", "campaign_audit_probe.py"),
    ("r428", "qn_reopened_probe.py"),
    ("r429", "zloc_head_probe.py"),
    ("r431", "source_potapov_probe.py"),
    ("r431-audit", "r431_audit_probe.py"),
    ("r433", "edge_redheffer_probe.py"),
    ("r435", "p1_overload_probe.py"),
    ("r436", "p2_determinant_probe.py"),
    ("r438", "evolutionary_certificate_probe.py"),
    ("r439", "residual_loewner_probe.py"),
    ("r440", "mean_tau_index_probe.py"),
    ("r441", "diag_lifts_loewner_probe.py"),
    ("r442", "block_mean_probe.py"),
    ("r443", "delta_floor_probe.py"),
    ("r444", "signed_border_mean_probe.py"),
    ("r445", "deep_builder_probe.py"),
    ("r446", "deep_abd_probe.py"),
    ("r447", "exact_atom_probe.py"),
    ("r448", "exact_band_probe.py"),
    ("r449", "flip_vs_stab_probe.py"),
    ("r450", "prefix_mincut_probe.py"),
    ("r451", "nstab_transition_probe.py"),
    ("r452", "plateau_theorem_probe.py"),
    ("r453", "border_mass_probe.py"),
    ("r454", "limit_object_probe.py"),
    ("r455", "arch_chain_probe.py"),
    ("r456", "vacuity_redteam_probe.py"),
    ("r457", "jp_increment_probe.py"),
    ("r458", "cofinal_family_probe.py"),
    ("r459", "fullcomb_cleanup_probe.py"),
    ("r460", "race_proof_probe.py"),
    ("r461", "narrowband_weil_probe.py"),
    # r462 gapmap_probe is a historical diagnostic whose success
    # condition is intentionally invalidated by the r463 repair.
    ("r463", "lean_fidelity_probe.py"),
    ("r464", "inner_bridges_probe.py"),
    ("r467", "overlap_mechanism_probe.py"),
    ("r468", "octave_renorm_probe.py"),
    ("r470", "quadrep_probe.py"),
    ("r471", "classical_cert_probe.py"),
    ("r473", "extraction_joint_probe.py"),
    ("r474", "modulus_upgrade_probe.py"),
    ("r475", "arch_rate_probe.py"),
    ("r476", "crossterm_probe.py"),
    ("r477", "highmode_probe.py"),
    ("r478", "endtoend_fixedl_probe.py"),
]

MODULES = [
    "v955_tau_iiks_toda_dictionary.py",
    "v956_signedmoment_halffilling_duality.py",
    "v958_bordered_tau_readout_dictionary.py",
    "v959_coupledtau_terminal_dictionary.py",
    "v960_terminal_surface_closure.py",
    "v961_midpoint_orientation_dictionary.py",
    "v962_halffilling_pinning_theory.py",
    "v963_lstar_reduction_dictionary.py",
    "v964_lstar_coherence_census.py",
    "v965_lstar_curvature_arc.py",
    "v966_l2_reduction_chain.py",
    "v967_l2_cascade_closure.py",
    "v968_architecture_adjudication.py",
    "v969_forks_and_redteam.py",
    "v970_extraction_and_composition.py",
    "v978_terminal_density_martingale.py",
    "v979_cover_growth_k2.py",
    "v980_lstar_margin_chain.py",
    "v981_lstar_borodin_duality.py",
    "v982_dirichlet_matched_frame.py",
]

PROBE_TIMEOUT = 900        # s per probe (smoke runs are seconds)
MODULE_TIMEOUT = 1800      # s per module (v959 ~210 s)
LEAN_TIMEOUT = 1800


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Suite:
    def __init__(self):
        self.results = []          # (name, ok, detail)

    def record(self, name, ok, detail=""):
        self.results.append((name, bool(ok), detail))
        print("  [%s] %-52s %s" % ("PASS" if ok else "FAIL", name, detail),
              flush=True)

    def info(self, msg):
        print("  [INFO] %s" % msg, flush=True)

    @property
    def failures(self):
        return [r for r in self.results if not r[1]]


def python_bin():
    if os.path.isfile(VENV_PY):
        return VENV_PY
    return sys.executable


def run_cmd(cmd, cwd, timeout):
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                              timeout=timeout)
        return proc.returncode, proc.stdout + proc.stderr, time.time() - t0
    except subprocess.TimeoutExpired:
        return -1, "TIMEOUT after %ds" % timeout, time.time() - t0


def section(title):
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74, flush=True)


def check_integrity(suite):
    section("(1) INTEGRITY -- SHA-256 vs rh/INVENTORY.json (drift detector)")
    if not os.path.isfile(INVENTORY):
        suite.record("inventory-present", False, "rh/INVENTORY.json missing")
        return
    with open(INVENTORY, encoding="utf-8") as f:
        inv = json.load(f)
    n_pin_ok = 0
    for e in inv["entries"]:
        p = os.path.join(REPO, e["path"])
        name = "sha256:" + os.path.basename(e["path"])
        if not os.path.isfile(p):
            if e["pin"]:
                suite.record(name, False, "file missing")
            else:
                suite.info("%s: unpinned file missing (living document)"
                           % e["path"])
            continue
        cur = sha256(p)
        if e["pin"]:
            ok = (cur == e["sha256"])
            if ok:
                n_pin_ok += 1
            else:
                suite.record(name, False,
                             "DRIFT pinned file changed (%s.. != %s..)"
                             % (cur[:12], e["sha256"][:12]))
        else:
            if cur != e["sha256"]:
                suite.info("%s: drift in unpinned living document (expected)"
                           % e["path"])
    n_pin = sum(1 for e in inv["entries"] if e["pin"])
    suite.record("inventory-pinned-files",
                 n_pin_ok == n_pin,
                 "%d/%d pinned entries byte-identical" % (n_pin_ok, n_pin))


def check_probes(suite):
    section("(2) PROBES -- sealed campaign r250-r459 + r431-audit, --smoke mode")
    py = python_bin()
    cwd = os.path.join(REPO, "experiments", "tfpt-discovery")
    for rnd, probe in PROBES:
        path = os.path.join(cwd, probe)
        if not os.path.isfile(path):
            suite.record("%s %s" % (rnd, probe), False, "missing")
            continue
        rc, out, dt = run_cmd([py, probe, "--smoke"], cwd, PROBE_TIMEOUT)
        ok = (rc == 0)
        tail = ""
        if not ok:
            lines = [ln for ln in out.strip().splitlines() if ln.strip()]
            tail = (lines[-1][:80] if lines else "no output")
        suite.record("%s %s" % (rnd, probe), ok,
                     "%.1f s%s" % (dt, ("  " + tail) if tail else ""))


def check_modules(suite):
    section("(3) MODULES -- v955-v970 by path from verification/")
    py = python_bin()
    cwd = os.path.join(REPO, "verification")
    for mod in MODULES:
        path = os.path.join(cwd, mod)
        if not os.path.isfile(path):
            suite.record(mod, False, "missing")
            continue
        rc, out, dt = run_cmd([py, mod], cwd, MODULE_TIMEOUT)
        ok = (rc == 0)
        tail = ""
        if not ok:
            lines = [ln for ln in out.strip().splitlines() if ln.strip()]
            tail = (lines[-1][:80] if lines else "no output")
        suite.record(mod, ok,
                     "%.1f s%s" % (dt, ("  " + tail) if tail else ""))


def check_lean(suite):
    section("(4) LEAN -- lake build in rh/lean/")
    lake = shutil.which("lake")
    if lake is None:
        suite.info("no `lake` on PATH -- Lean gate skipped "
                   "(install elan, see rh/lean/README.md)")
        return
    if not os.path.isfile(os.path.join(LEAN_DIR, "lakefile.toml")):
        suite.record("lean-project-present", False, "rh/lean missing")
        return
    rc, out, dt = run_cmd([lake, "build"], LEAN_DIR, LEAN_TIMEOUT)
    ok = (rc == 0) and ("Build completed successfully" in out)
    n_sorry = out.count("declaration uses `sorry`")
    suite.record("lake build (RH)", ok,
                 "%.1f s, %d intentional `sorry` warnings" % (dt, n_sorry))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fast", action="store_true",
                    help="skip the v9xx modules (minutes each)")
    ap.add_argument("--skip-lean", action="store_true",
                    help="skip the lake build gate")
    args = ap.parse_args()

    t0 = time.time()
    print("RH SUITE -- rh/ workspace  (research documentation, NO RH CLAIM)")
    print("repo: %s" % REPO)
    print("python: %s" % python_bin())

    suite = Suite()
    check_integrity(suite)
    check_probes(suite)
    if args.fast:
        section("(3) MODULES -- SKIPPED (--fast)")
    else:
        check_modules(suite)
    if args.skip_lean:
        section("(4) LEAN -- SKIPPED (--skip-lean)")
    else:
        check_lean(suite)

    dt = time.time() - t0
    npass = sum(1 for r in suite.results if r[1])
    print("\n" + "-" * 74)
    print("items: %d pass / %d total, wall %.1f s"
          % (npass, len(suite.results), dt))
    if suite.failures:
        print("RH SUITE: FAILURES %d" % len(suite.failures))
        for name, _, detail in suite.failures:
            print("  FAILED: %s  %s" % (name, detail))
        return 1
    print("RH SUITE: ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
