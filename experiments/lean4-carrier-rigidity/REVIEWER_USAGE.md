# Reviewer Usage Instructions

This package is the Lean 4 companion project for the TFPT carrier-rigidity audit surface. It is intended for local reproduction of the formal checks, not for Overleaf compilation.

## 1. Scope of the Package

The package verifies the algebraic carrier core relative to explicit Lean input interfaces. In particular, the checked surface includes:

- the carrier polynomial statement `sixY^2 - sixY - 6 = 0`;
- the primitive integer rigidity result `(q_-, q_+) = (-2, 3)`;
- the concrete hypercharge readouts `Tr_E Y = 0` and `6Y^2 - Y - 1 = 0`;
- upstream typed interfaces such as `BoundaryPolarization`, `BoundaryYukawaKernelInterface`, and `SeamWindingInterface`;
- audit modules that check theorem names, theorem signatures, absence of deferred proofs, and axiom dependencies.

This package does not claim unconditional closure of the full TFPT-origin gate. It should be read as the Lean-verified algebraic carrier core `[L]`, plus typed Lean interfaces `[I]` for upstream inputs, with manuscript-level and standard-mathematics components kept separate.

## 2. Contents

The important files are:

- `lean-toolchain`: pins Lean to `leanprover/lean4:v4.29.1`;
- `lakefile.lean`: declares the Lake package and pins Mathlib to `v4.29.1`;
- `lake-manifest.json`: records exact dependency revisions;
- `TfptCarrier.lean`: top-level import file;
- `TfptCarrier/*.lean`: Lean source modules;
- `scripts/audit.sh`: the hard audit script;
- `scripts/build_wall_ladder.sh`: the bounded-memory builder for the generated wall-ladder rungs (see §6.1);
- `README.md`: project overview and theorem map.

The archive intentionally excludes `.lake/`, compiled `.olean` files, local build caches, and generated LaTeX artifacts. These are recreated locally.

## 3. Prerequisites

Install `elan`, Lean's toolchain manager:

```bash
curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh -s -- -y --default-toolchain none
export PATH="$HOME/.elan/bin:$PATH"
```

You also need:

- `git`;
- a POSIX shell such as `bash` or `zsh`;
- internet access for the first dependency download;
- several GB of free disk space for Mathlib cache files.

## 4. Unpack and Enter the Project

From the directory where you received the zip file:

```bash
unzip lean4-carrier-rigidity-reviewer.zip
cd lean4-carrier-rigidity
```

Check the pinned Lean version:

```bash
cat lean-toolchain
```

Expected output:

```text
leanprover/lean4:v4.29.1
```

## 5. Fetch Dependencies

Run:

```bash
lake update
lake exe cache get
```

`lake update` fetches the pinned dependencies recorded by the Lake project. `lake exe cache get` downloads prebuilt Mathlib artifacts, which avoids compiling all of Mathlib from source.

The first run can take several minutes and may download a large cache.

## 6. Run the Full Audit

From the project root:

```bash
./scripts/audit.sh
```

If the script is not executable after unzipping, run:

```bash
chmod +x scripts/audit.sh
./scripts/audit.sh
```

Expected final output:

```text
=== Final verdict ===
  AUDIT: PASS
```

The audit script checks:

1. `lake build` succeeds.
2. No `sorry` or `admit` occurs in `TfptCarrier/` or `TfptCarrier.lean`.
3. No domain-specific `axiom` or `constant` declarations occur.
4. No `unsafe`, `opaque`, or `partial` declarations occur.
5. `#print axioms` reports only Lean's standard axioms `propext`, `Classical.choice`, and `Quot.sound`.
6. `AuditCheck.lean` elaborates all headline theorem names.
7. `AuditContract.lean` elaborates exact theorem signature locks.

### 6.1 Prerequisite for check (1): the generated wall-ladder rungs

`TfptCarrier/WallLadder/RungKz*.lean` are generated certificate modules, each a single kernel `decide` over a packed Cholesky witness. They are the expensive part of check (1), and `lakefile.lean` pins a low default per-process memory ceiling (`leanMemoryMb = 12288`, passed to Lean as `-M`) so that an uncapped Lake fan-out cannot exhaust physical RAM. Consequence: on a fresh checkout the rungs **fail at the default ceiling by design**, and `./scripts/audit.sh` cannot report `AUDIT: PASS` until their `.olean` files exist.

Build them first, with a raised ceiling and bounded concurrency:

```bash
BATCH=2 MEM_MB=163840 ./scripts/build_wall_ladder.sh
```

`BATCH` is how many rungs one `lake build` invocation may run at once, `MEM_MB` is the per-process `-M` ceiling in MB, and already-built rungs are skipped, so an interrupted run can be restarted. Each batch prints its wall time and the peak summed RSS of all live `lean` processes.

Budget honestly. Cost rises steeply with the size of the generated rung source — time roughly quadratically, memory faster than linearly. Measured on a 512 GB Apple-silicon machine: `RungKz12` (0.12 MB source) 189 s, `RungKz20` (0.15 MB) 253 s, `RungKz14` (0.17 MB) 363 s, `RungKz32` (0.32 MB) 1011 s at ~93 GB, the pair `RungKz39`+`RungKz46` (0.37 + 0.39 MB) 1598 s at 360.3 GB summed, and `RungKz27` (0.41 MB) 1700 s at **200.4 GB in a single process**. `RungKz14` fails outright at a 48 GB ceiling (`(kernel) excessive memory consumption detected` after 987 s), so a too-small ceiling looks like a build error rather than a resource warning.

Size the run to your machine: `BATCH=2` with `MEM_MB=163840` is comfortable up to ~0.35 MB rungs on 512 GB, while from ~0.4 MB upward a single rung wants 200 GB or more and must be built with `BATCH=1` and a several-hundred-GB ceiling. On a workstation with 64 GB or less only the smaller rungs are reachable at all — in that case report the rungs you could not build rather than treating the audit as passed, since check (1) is then not satisfied.

Because `-M` is passed through `weakLeanArgs`, it is not part of Lake's build trace: an `.olean` produced under a raised ceiling is accepted unchanged by the subsequent default-ceiling `lake build` inside `./scripts/audit.sh`. Building the rungs this way is a memory-scheduling decision, not a change to what is checked.

## 7. Useful Focused Commands

Build only the Lean library:

```bash
lake build TfptCarrier
```

Check theorem-name elaboration:

```bash
lake build TfptCarrier.AuditCheck
```

Check exact theorem signatures:

```bash
lake build TfptCarrier.AuditContract
```

Run static checks without rebuilding:

```bash
./scripts/audit.sh --quick
```

`--quick` is only a convenience mode. The full audit is the authoritative reproduction command.

## 8. How to Read the Status

For the audit map, use the following separation:

- `[L]`: Lean-verified statements in this package;
- `[I]`: typed Lean interface targets imported as upstream inputs;
- `[M]`: manuscript-level proof or discussion outside this Lean package;
- `[S]`: standard external mathematics cited outside this Lean package;
- `[C]`: numerical or computational results requiring a reproducible calculation ledger;
- `[P]`: programmatic or phenomenological pipeline status.

The intended Paper 2/D-03 reading is:

```text
carrier core [L] / manuscript discharge [M] / TFPT-origin gate [I]-conditional
```

This means the carrier arithmetic is machine-checked here, while the full TFPT-origin closure remains visibly conditional on upstream formal interfaces and manuscript-level inputs.

## 9. Troubleshooting

If `lake` is not found, ensure elan is installed and on your `PATH`:

```bash
export PATH="$HOME/.elan/bin:$PATH"
lake --version
```

If dependency fetching fails, retry after checking network access:

```bash
lake update
lake exe cache get
```

If the build fails after a Lean or Mathlib version change, restore the pinned files from the archive:

```text
lean-toolchain
lakefile.lean
lake-manifest.json
```

Do not silently upgrade Lean or Mathlib when reproducing the audit. Any version bump should be treated as a new audit run.

