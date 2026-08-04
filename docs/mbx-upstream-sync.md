# RASPA3-MBX upstream synchronization report

This document records the 2026-08 integration of the Chung Research Group
RASPA3-MBX fork with current RASPA3. It is both a change checklist for this
update and the maintenance policy for later updates.

## Repository baseline

| Item | Commit or location |
|---|---|
| Original MBX working tree | Historical local `RASPA3.0.21-MBX` checkout |
| Original MBX branch | `E_print` at `186f3b8cc6c4bbc28cbdbc5954ddd05b63357c4f` |
| Upstream release on which the fork was named | v3.0.21 at `729be4b40ca6b011488fd83d2a3d94eb15e41a66` |
| True Git merge base | `6b8b25f2c944bbea751fcb63eb79668040c87042` |
| Imported upstream | `a8a70ca7` (`feat: structurekit`, current `upstream/main` on 2026-08-03) |
| Latest upstream merge commit | `afab6925` |
| Integration working tree | This `RASPA3-MBX` repository |
| Integration branch | `main` (tested publication branch) |

The merge base is 23 upstream commits newer than the v3.0.21 release tag.
From that true base, the fork contains 57 commits and upstream contains 248
commits. Current upstream is 271 commits beyond the v3.0.21 release. Using the
true merge base is important: comparing two directory snapshots alone makes
upstream refactors look like local MBX changes and leads to the wrong conflict
choices.

The historical checkout was not modified during integration. This repository
was created as a compact Git clone without the approximately 3.7 GB ignored
build directory. The original untracked MOF-74 Morse inputs and structure were
copied into this working tree so that scientific inputs were not lost.

The integration copy has two remotes:

- `origin`: the Chung Research Group RASPA3-MBX repository.
- `upstream`: the official iRASPA/RASPA3 repository.

## Required forward ports and conflict decisions

| Area | Required modification and decision |
|---|---|
| Build system | Keep the current upstream CMake target layout and C++ module architecture. Make MBX optional through `BUILD_MBX`; discover its headers, generated configuration header, library, and FFTW through cache variables rather than hard-coded Carbon/Oxygen paths. Add optional installed HDF5 and GoogleTest paths for offline/reproducible builds. |
| System state | Add `useMBX`, the canonical MBX settings path, `writeEnergyLog`, the zero-loading Widom heat control, the precomputed bare-framework permanent energy, and logger state to the current upstream `System`. Preserve all new upstream thermobarostat, group, flexible-framework, interpolation, and property state. |
| Input reader | Parse and validate `UseMBX`, `MBXSettingsFile`, canonical `PrintEnergyTerms` (with `WriteEnergyLog` as a legacy alias), `ComputeZeroLoadingHeatOfAdsorption`, `WriteRestartEvery`, and `EnergyEvaluation`. Resolve relative settings and coordinate snapshots against the simulation input directory, validate restart coordinates, and fail clearly if the executable was built without MBX. |
| MBX adapter | Port `interactions_mbx` to the new module interfaces. Use stack ownership/RAII and spans, validate component/atom layouts, keep atom tags within MBX's integer range, and avoid the old raw-allocation and `noexcept` failure modes. |
| Total-energy contract | Prevent classical RASPA guest VDW, guest charge, Ewald, polarization, and intramolecular terms from being added on top of MBX. Retain external-field energy and classical framework–guest VDW plus its framework-tail correction. |
| Running energy | Start from current upstream `RunningEnergy`, including thermobarostat energy and per-group `dU/dlambda` arrays, and add one `mbxEnergy` field to totals, arithmetic, clearing, text, JSON, and archives. Do not restore the fork's obsolete split-tail fields. |
| Energy statistics | Add MBX-aware `EnergyStatus` arithmetic, block statistics, text output, and archives. An MBX total is external field + framework–guest VDW + framework tail + MBX. Classical component decompositions remain available for the pressure estimator but do not enter the MBX total. |
| MC moves | Forward-port MBX acceptance and bookkeeping to translation, rotation, conventional insertion/deletion, CBMC insertion/deletion, CBMC reinsertion, Widom, and isotropic framework-free volume changes. Preserve current upstream TMMC updates, dual cutoffs, Ewald trial caches, polarization caches, CBMC Rosenbluth handling, and timing. |
| Drift correction | Recompute the old and proposed full MBX energy for an MBX move. Use their physical difference in acceptance and the exact accepted MBX total in diagnostic output. This incorporates the important stale-bookkeeping fix from fork commit `e48115b0`. |
| Energy logs | Replace scattered `std::cerr` fragments with persistent per-system CSV streams under `output/`. Accepted state-changing moves use `energy_terms.s<id>.csv`; completed conventional Widom ghost insertions use the distinct `widom_energy_terms.s<id>.csv` schema and never emit an accepted row. `PrintEnergyTerms` defaults to `true` and can be disabled explicitly. |
| Zero-loading Widom heat | Add the opt-in rigid-host/rigid-adsorbate, fixed-volume estimator `Qst^0 = kBT - <W DeltaU>/<W>`, block uncertainty, K and kJ/mol output, and hard guards that guest loading remains zero. |
| Restart cadence | Replace the hard-coded 5000-cycle coordinate-snapshot interval with top-level `WriteRestartEvery`; `0` disables periodic writes. Count completed cycles, preserve the value in version-2 Monte Carlo binary archives, and retain version-1 read compatibility. |
| One-shot energy | Add `SimulationType: EnergyEvaluation` to evaluate an input configuration or coordinate restart once, without MC moves. Emit a human-readable decomposition and `output/energy_evaluation.json`, including all seven raw MBX terms when enabled. |
| CPU timing | Add a dedicated `Move::Timing::MBX` entry without renumbering existing timing values. Bump the CPU-time archive version because the serialized timing row size changes. |
| Output | Include MBX in current MC/MD/repr/JSON energy output, print the MBX settings/status in the normal system report, and expose the MBX settings and current MBX energy in system JSON. |
| Swap bookkeeping | Fix current-upstream conventional insertion bookkeeping in pre-initialization and equilibration: an accepted insertion must be added to `runningEnergies`, as production already does. |
| Tests | Restore the fork's MBX unit tests, adapt them to the current targets, and gate them with `BUILD_TESTING_MBX`. Add a multi-component global-ID regression. Preserve upstream tests and use current native Morse-potential tests as the baseline for Morse behavior. Restore the `tests/test_support.hpp` RAII helper referenced by current upstream tests but absent from the downloaded upstream tree. |
| Relocated drivers | Use the new upstream files in `src/raspakit/drivers/`. Do not revive the deleted top-level Monte Carlo, TMMC, MD, or old VDW module files. |
| Multidimensional TMMC | Preserve `transition_matrix_nd.cpp/.ixx` and its calculation examples, but do not silently graft the old multidimensional implementation into the substantially changed current TMMC driver. It remains an unwired feature that should be rebased and tested on its own branch. |
| Generated data | Do not carry the old 38,000-line debug output or build products into the update. Scientific input files are retained; generated output belongs in ignored run directories. |

## MBX energy definition

The retained total potential energy in an MBX Monte Carlo run is

```text
U = U_external-field + U_framework-guest-VDW + U_framework-tail + U_MBX
```

`U_MBX` contains the following MBX terms:

```text
2B + 3B + 4B + dispersion
   + (permanent electrostatics - bare-framework permanent electrostatics)
   + induced electrostatics
```

The MBX 1B value is recorded as a diagnostic raw term but is deliberately not
included in the retained total, matching the established RASPA3-MBX model.
MBX Buckingham and Lennard-Jones energies are also excluded. Classical RASPA
supplies only the framework–guest VDW contribution. This division is explicit
so that guest electrostatics, Ewald, and polarization are not counted twice.

The accepted-move MBX term order is

```text
[1B, 2B, 3B, 4B, dispersion, permanent electrostatics, induced electrostatics]
```

Raw MBX values arrive in kcal/mol and are converted to the RASPA internal
energy unit before they are combined with `RunningEnergy`. The historical CSV
column order is preserved.

## Supported MBX execution surface

MBX trajectory generation currently supports the regular single-system
`MonteCarlo` driver only. The non-mutating `EnergyEvaluation` analysis mode is
also supported.
The input reader rejects unported combinations instead of allowing a classical
move to corrupt an MBX trajectory.

| Move or mode | Status |
|---|---|
| Translation | Supported |
| Rotation | Supported |
| Reinsertion CBMC | Supported |
| Conventional insertion/deletion (`SwapConventionalProbability`) | Supported for rigid components; flexible components must use CBMC |
| CBMC insertion/deletion (`SwapProbability`) | Supported |
| Widom | Supported |
| Isotropic volume, framework-free box without an external field | Supported |
| Isotropic volume with a framework | Rejected |
| Isotropic volume with an external field | Rejected; an accepted box change would invalidate the field grid |
| Random translation / random rotation | Rejected; these are separate new upstream implementations, not wrappers |
| Smart/gradient moves and hybrid MC | Rejected |
| Partial reinsertion and identity change | Rejected |
| Anisotropic volume | Rejected |
| CFCMC, pair/group swap, Gibbs, reactions | Rejected |
| Transition-matrix and parallel-replica drivers | Rejected |
| Molecular dynamics, minimization, thermodynamic integration | Rejected; MBX gradients are not forward-ported |
| One-shot `EnergyEvaluation` | Supported for JSON coordinate snapshots; no moves or MBX gradients are used |
| Classical RASPA polarization together with MBX | Rejected; MBX supplies the retained induced term |
| Force-based RDF sampling | Rejected; it invokes the classical full-gradient path and would overwrite MBX running energy |

NPT acceptance for the supported framework-free isotropic-volume move uses the
recomputed MBX energy. The pressure tensor reported by current RASPA3-MBX is
still the classical molecular-pressure estimate with MBX energy attached for
reporting; it is not an MBX virial. Do not interpret it as a many-body MBX
pressure until MBX strain derivatives/forces are integrated explicitly.

This is intentionally conservative. A new move should be enabled only after it
has both an MBX energy-difference implementation and a test showing that
incremental `runningEnergies` agrees with a full recomputation.

## Input and local build configuration

An MBX system uses the following current input fields:

```json
{
  "UseMBX": true,
  "MBXSettingsFile": "mbx.json",
  "PrintEnergyTerms": true,
  "ComputeZeroLoadingHeatOfAdsorption": false
}
```

`MBXSettingsFile` may be relative to `simulation.json`. `PrintEnergyTerms` is
optional and defaults to `true` for compatibility with the fork. Set it to
`false` for production without move/trial energy rows. The legacy
`WriteEnergyLog` alias remains accepted; conflicting values are rejected.
Regular Monte Carlo writes accepted state changes to
`output/energy_terms.s0.csv`; when conventional Widom moves are enabled it
also writes completed ghost trials to `output/widom_energy_terms.s0.csv`.
Both truncate on a fresh run and append on binary resume without a duplicate
header. The CSV and binary checkpoint are separate streams: after an abrupt
crash, rows written after the last binary checkpoint can be replayed and
duplicated in the appended tail. Treat them as diagnostic output and split or
deduplicate that final segment when resuming a checkpoint.

Coordinate restart cadence is controlled separately:

```json
{
  "WriteRestartEvery": 5000
}
```

For regular Monte Carlo and thermodynamic integration, `0` disables periodic
coordinate JSON writes, but the drivers' existing initial and final snapshot
writes remain. This is distinct from `WriteBinaryRestartEvery`, which controls
full crash-recovery checkpoints.

For a fixed snapshot comparison use:

```json
{
  "SimulationType": "EnergyEvaluation",
  "Systems": [{"RestartFileName": "output/restart_298_100000.s0.json"}],
  "Components": [{"Name": "co2", "CreateNumberOfMolecules": 0}]
}
```

The abbreviated fragment above must be combined with the normal force-field,
system, and component definitions. The result is written to
`output/energy_evaluation.json`; no accepted-move CSV or new restart snapshot is
created.

Machine-specific paths live in ignored `CMakeUserPresets.json`, not in the
shared `CMakePresets.json`. A portable template is provided as
`CMakeUserPresets.json.example`. The prepared local workflow is:

```bash
cmake --preset mbx-local --fresh
cmake --build --preset mbx-local --parallel 4
ctest --preset mbx-local
```

The audited Carbon build used CMake 4.0.3, Clang 20, public MBX v1.4.0 at
`0e01b75b47611d7d51f27a34b112bdc5e2090a50`, and Conda
FFTW/HDF5/GoogleTest packages. The root CMake file also retains the token
required by CMake 4.3 and newer. MBX and RASPA were compiled against LLVM
OpenMP (`libomp`); linking with GNU `libgomp` as a second runtime is unsafe.
The public dependency comparison and clean-build recipe are recorded in the
[public-MBX migration audit](mbx-public-migration.md).

## Verification performed

The integration was validated in the prepared Carbon environment as follows:

| Check | Result |
|---|---|
| Public MBX dependency | Fresh v1.4.0 build and install passed 34/34 upstream tests and resolved one LLVM OpenMP runtime |
| Current MBX-enabled build | Passed against the public MBX prefix after the StructureKit merge, including `libraspakit_base.a`, `app/raspa3`, and `unit_tests_mbx` |
| Complete MBX-disabled build | Passed from the current sources in a separate configuration, with no MBX headers or library in the target |
| User-facing executables | `app/raspa3` and `cli/raspa3-cli` linked successfully and both passed help/startup checks |
| Focused MBX/diagnostic tests | 15/15 passed: four MBX energy/bookkeeping tests, five accepted/Widom energy-log tests, two zero-loading Widom-heat tests, one restart-cadence/archive test, and three fixed-snapshot evaluation tests |
| Broader current-source test suite | 532/535 enabled tests passed, with 6 additional tests disabled (541 registered); the three failures are in unchanged hybrid-MC, second-order Taylor-shifted Lennard-Jones, and exact-sphere-pruning paths described below |
| Native Morse tests | 7/7 relevant tests passed: reference energy, shifted cutoff, fractional-lambda finiteness, lambda derivative, spatial derivatives, bonded Morse, and Urey-Bradley Morse |
| Two-CO2 MBX application smoke test | Passed through the real input reader and `raspa3` executable with initialization and 20 production MC steps |
| Retained MOF-74 Morse input | Passed with one CO2, MBX enabled, the two Mg-CO2 Morse pairs active, and 20 translation/rotation MC steps; the reported incremental-versus-recomputed MBX drift was zero at output precision |
| Accepted-energy file smoke test | A fresh two-cycle Monte Carlo run configured with `WriteRestartEvery: 1` wrote one-header `output/energy_terms.s0.csv` and produced no energy rows on standard error |
| Classical snapshot evaluation | The MBX-disabled executable loaded a relative two-molecule coordinate restart, evaluated it once without moves, and created only `output/energy_evaluation.json` |
| MBX + Morse snapshot evaluation | A nonzero one-CO2 MOF-74 snapshot produced the classical Morse framework term, retained MBX total, and all seven raw MBX terms without changing the snapshot |
| Private/public fixed-snapshot parity | All seven raw terms agreed; the largest difference was `8.731e-11` kcal/mol and the retained/total difference was `4.394e-8` K |
| Merge hygiene | No unresolved paths or conflict markers; the original source repository remains unchanged |

The native RASPAKit test executable itself also builds after restoring the
missing upstream test-support header. Three broader tests in files unchanged by
this diagnostic feature set fail reproducibly:

- `hybrid_mc.flexible_framework_only_does_not_throw_and_restores_on_reject`
  expects its deterministic proposal to be rejected, but it is accepted;
- `vdw_potentials.second_order_taylor_shifted_spatial_derivatives_match_finite_difference`
  misses its finite-difference tolerance at `r = 5` by about `1.75e-4` beyond
  the allowed error; and
- `exact_sphere_sweep.pruning_keeps_one_of_a_pair_of_equals_and_the_order_of_the_rest`
  retains three circles where the test expects two.

The relevant Morse-bearing tests and all focused MBX/new-control tests pass.
These three failures should be investigated separately instead of being hidden
by changing tolerances or expectations during the MBX synchronization.

The accepted-probability formulas for CBMC insertion/deletion, reinsertion,
and Widom were also checked algebraically. Their MBX correction is evaluated
as one log-space expression, avoiding the former `0 * infinity -> NaN` case
when large native and MBX exponents cancel. The classical paths are unchanged.

## Morse potential audit

Do not submit a pull request that merely adds Morse potential support. Current
upstream already contains a much more complete native implementation than the
fork's final `add morse related terms` commit.

Current upstream supports:

- the parameter order `[D, a, r0]` for self and binary interactions;
- the shifted potential
  `D * ((1 - exp(-a * (r - r0)))^2 - 1)`;
- lambda scaling, gradients, strain derivatives, and Hessians;
- energy and pressure tail corrections;
- serialization and direct CPU potential tests.

The old fork patch was incomplete: its self-interaction path mishandled `r0`,
and it did not supply the current derivative, Hessian, tail, serialization, or
test coverage. The integration therefore keeps the upstream implementation.

### Scientific compatibility warning

The retained MOF-74 Morse force field sets `"TailCorrections": true`.
Current upstream evaluates a Morse tail correction, whereas the old fork did
not. Results may therefore differ even with identical `[D, a, r0]` values.
For a controlled comparison, first run old and new executables with tail
corrections disabled, then enable the new correction and report that change as
a model change rather than a merge effect.

Morse is supported in direct CPU energy/force routes. Grid interpolation and
some OpenCL kernels remain Lennard-Jones-specific. Do not enable grids/OpenCL
for a Morse pair until the route either implements Morse or validates and
rejects it explicitly.

### Useful upstream Morse pull request

A focused upstream contribution would complete Morse I/O and reporting rather
than reimplementing the potential:

1. Add Morse entries to `ForceField::jsonForceFieldStatus()`; it currently
   emits detailed interaction JSON only for Lennard-Jones.
2. Include shift/tail status in the Morse text report, consistently with other
   VDW potentials.
3. Document the formula, `[D, a, r0]` order, units, binary-interaction example,
   and the fact that arbitrary mixing is not defined for Morse overrides.
4. Add parser tests for an exact three-parameter count, missing `r0`, self and
   binary interactions, and energy-unit conversion.
5. Add JSON/text/archive and tail-integral regression tests.
6. Decide explicitly whether extra parameters should be rejected; the current
   parser can accept more than the declared metadata count silently.

That work should be a small standalone branch based directly on upstream
`main`, with no MBX dependency. A later, separate pull request can add
validation or implementations for Morse in interpolation/OpenCL routes.

## Remaining constraints and hardening work

- Normal input-driven runs, restart-driven `MonteCarlo::run()`, stages, and
  `performCycle()` are centrally validated. Calling an individual low-level
  move function directly can still bypass that contract; such functions are
  not a supported MBX application interface.
- MBX trial safety currently relies on the normal classical proposal
  interactions to reject severe overlaps before MBX is called. Do not use
  `omitInterInteractions` or zero/non-repulsive proposal pair potentials in an
  MBX run. A future hardening change should add a geometry/model-independent
  close-contact screen before every MBX evaluation.
- Final CBMC/Widom reweighting is stable in log space, but the upstream CBMC
  generator still accumulates Rosenbluth weights in ordinary floating point.
  Extremely large trial weights can therefore overflow before the final ratio;
  a complete solution requires log-sum-exp accumulation inside CBMC itself.
- The reported pressure tensor is still a classical virial with the MBX energy
  attached. MBX force/strain derivatives are required for a true MBX virial.
- The current Morse text status omits `r0`, and detailed force-field JSON is
  Lennard-Jones-centric. These are reporting gaps, not energy-path gaps, and
  are the best scope for a small upstream Morse pull request.

## Restart compatibility

`System`, `RunningEnergy`, `EnergyStatus`, MC timing, and Monte Carlo driver
archive versions are bumped where their layouts changed. Monte Carlo archive
version 2 stores `WriteRestartEvery`; version-1 checkpoints retain the input or
default 5000-cycle value when read. Current upstream version-1 `System` objects
can be read with MBX state defaulted off. Old RASPA3.0.21-MBX restart files are
not a safe interchange format: that fork reused archive version 1 for a
different layout and its `RunningEnergy` writer and reader were themselves
inconsistent. Start a fresh run for the updated executable rather than
resuming an old MBX binary checkpoint.

## Long-term branch and update plan

### One canonical repository

Use this `RASPA3-MBX` repository as the sole active checkout. It contains the
original MBX/Morse history plus the current upstream history, and the old local
`RASPA3.0.21-MBX` checkout has no tracked working-tree changes that need to be
merged in place. Updating the old folder would retain its obsolete
multi-gigabyte build tree and would make the repository ancestry harder to
audit.

The regression matrix is complete and the tested integration branch is
`main`. The safe one-time publication or recovery sequence is:

1. Create a Git bundle of the old repository and preserve any deliberately
   untracked research files before removing anything.
2. Review the remote diff, then push `main` to the Chung Research Group fork
   and make it the default branch. The old `E_print` and `development` tips may
   remain as archival remote branches until the first tagged release.
3. Only after cloning the promoted remote into a temporary directory and
   rerunning a smoke test should the old checkout be archived or deleted.

No source merge back into the old folder is needed. Future updates happen in
this one checkout by fetching and merging `upstream/main`; build directories
are disposable artifacts, not repository versions. Remote pushes, default
branch changes, and removal of the old checkout remain explicit user actions.

### Upstream synchronization workflow

Keep the fork as a thin set of reviewable layers instead of one long-lived
mixed feature branch:

```text
upstream/main                    exact mirror of iRASPA/RASPA3
main                             tested RASPA3-MBX release integration
feature/mbx-core                 adapter, System fields, input validation
feature/accepted-energy-log      optional diagnostic reporting
feature/tmmc-nd                  multidimensional TMMC, independently tested
sync/upstream-YYYYMMDD-SHA       temporary upstream-integration branch
```

The feature branches describe ownership; they need not be rebased and merged
manually for every release if the same separation is maintained as a short,
logical commit series on `main`. The crucial point is that MBX, logging, and
multidimensional TMMC are not interleaved in the same commits.

Recommended cadence and process:

1. Fetch upstream at least monthly and immediately for every upstream release
   or force-field/MC-move change.
2. Record the exact upstream SHA, create `sync/upstream-YYYYMMDD-SHA` from the
   last tested MBX `main`, and merge `upstream/main` with `--no-ff`.
3. Resolve architecture/build conflicts in favor of upstream first. Reapply the
   small MBX interfaces according to the energy contract in this document.
4. Keep `git rerere` enabled to reuse only identical mechanical conflict
   resolutions; always review energy and archive conflicts manually.
5. Run both an upstream-only configuration (`BUILD_MBX=OFF`) and the MBX CI
   configuration (`BUILD_MBX=ON`).
6. Run upstream tests plus the MBX absolute-energy, accepted-move drift,
   accepted/Widom logger, zero-loading Widom-heat, and Morse direct-CPU
   regressions.
7. Compare one short fixed-seed MOF-74 trajectory against the last released
   fork and explain expected model changes such as new Morse tails.
8. Merge the sync branch into `main`, tag it as
   `raspa3-mbx-<upstream-version>.<fork-revision>`, and record both SHAs in the
   release notes.

Recommended CI matrix:

| Job | Purpose |
|---|---|
| Upstream-only, MBX off | Proves the optional integration does not make normal RASPA depend on MBX headers/libraries |
| MBX unit tests | Absolute MBX decomposition and unit conversion |
| MBX supported-move drift tests | Each accepted move agrees with `computeTotalEnergies()` within tolerance |
| Input rejection tests | Every unsupported driver/move fails before a simulation starts |
| Native Morse CPU tests | Energy, gradient, Hessian, and tail regressions |
| Short MOF-74 smoke run | Exercises real settings, framework charges, logging, and insertion/deletion |

For reproducibility, provide MBX through a pinned package/container or a CI
artifact built from a recorded MBX commit. Never commit local build trees,
absolute dependency presets, production output, or downloaded dependencies.

## Pull-request policy

Do not open or push a pull request automatically from the integration branch.
The RASPA3-MBX integration update belongs in the fork repository after local
review. The Morse I/O/reporting work described above should be prepared
separately from a clean official-upstream branch. Publishing either branch
requires an explicit review of the final diff, commit message, target
repository, and target branch.
