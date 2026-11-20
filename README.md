# RASPA3-MBX

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
[![Upstream: iRASPA/RASPA3](https://img.shields.io/badge/upstream-iRASPA%2FRASPA3-2f80ed.svg)](https://github.com/iRASPA/RASPA3)
[![MBX: optional](https://img.shields.io/badge/MBX-optional-6f42c1.svg)](#building)

RASPA3-MBX is the Chung Research Group integration of
[RASPA3](https://github.com/iRASPA/RASPA3) with the MBX many-body
energy model. It retains the normal classical RASPA3 build while adding an
optional MBX path for supported Monte Carlo calculations, Morse-potential
inputs, accepted-move and Widom-trial energy diagnostics, configurable coordinate restarts,
and one-shot evaluation of coordinate snapshots.

Use `main` or `mbx` for the combined rigid-MBX and two-component TMMC implementation.
The [standalone 2d-tmmc branch](https://github.com/Chung-Research-Group/RASPA3-MBX/tree/2d-tmmc)
provides the force-field implementation on the v3.0.30 release base. For the Linux
executable, start with the [branch installation guide](INSTALLATION.md). The version printed in simulation output identifies the upstream release. Record the fork commit (`git rev-parse HEAD`) alongside it to identify the MBX and TMMC code used.

## What this fork adds

| Capability | RASPA3-MBX behavior |
|---|---|
| MBX guest energy | Enabled per system with `"UseMBX": true` and an `MBXSettingsFile` |
| Two-component TMMC | Set per-component population windows for force-field or rigid MBX exchange. See the [input reference](docs/manual/commands.md#joint-particle-number-windows-for-transition-matrix-monte-carlo). |
| Morse potential | Supported by the current native RASPA3 force-field implementation and retained in the MOF-74 inputs |
| Energy-term logs | Opt in with `PrintEnergyTerms: true` to write accepted moves to `output/energy_terms.s<id>.csv` and conventional Widom trials to `output/widom_energy_terms.s<id>.csv` |
| Zero-loading heat | `ComputeZeroLoadingHeatOfAdsorption` reports energy-weighted conventional-Widom (Q_{st}^0) statistics |
| Coordinate restart cadence | Top-level `WriteRestartEvery` replaces the previous hard-coded interval of 5000 cycles |
| Snapshot comparison | `SimulationType: EnergyEvaluation` reads a coordinate-restart JSON, evaluates it once, and exits without a Monte Carlo move |
| Upstream maintenance | The official repository is kept as the `upstream` remote and integrated through reviewable synchronization branches |

The retained MBX total is

```text
U = U_external-field + U_framework-guest-VDW + U_framework-tail + U_MBX
```

Raw MBX values are received in kcal/mol and converted to RASPA's internal
kelvin energy convention. The raw MBX one-body term is reported for diagnosis
but is deliberately excluded from the retained RASPA3-MBX total, matching the
established model in this fork. See the
[energy-definition section](docs/mbx-upstream-sync.md#mbx-energy-definition)
before comparing results with another program.

## Supported MBX calculations

Rigid MBX trajectory generation supports the single-system `MonteCarlo` and
`MonteCarloTransitionMatrix` drivers. Supported state-changing moves include translation,
rotation, CBMC reinsertion, conventional or CBMC insertion/deletion, and a
framework-free isotropic volume move without an external field. All guest components must be rigid, including for CBMC and `EnergyEvaluation`.
For TMMC use a fixed volume and a fixed host. The documented joint particle-number
windows support conventional and CBMC exchange of rigid guests.
See [component and system macrostate bounds](docs/manual/commands.md#joint-particle-number-windows-for-transition-matrix-monte-carlo).
Flexible-guest targets are maintained on development branches. Widom evaluation and the non-mutating `EnergyEvaluation` mode are supported.

The input reader rejects unported MBX combinations before starting a run. In
particular, do not enable MBX for molecular dynamics, minimization,
thermodynamic integration, CFCMC, Gibbs/reaction moves or
parallel-replica drivers, smart/gradient moves, anisotropic volume changes, or
flexible-framework volume changes. The full support matrix and pressure/virial
limitations are in
[Supported MBX execution surface](docs/mbx-upstream-sync.md#supported-mbx-execution-surface).

Classical RASPA3 calculations remain available by building with
`BUILD_MBX=OFF` or by leaving `UseMBX` false.

Known MBX monomers require database site order: `co2` is C,O,O and `h2o` is O,H,H.
The component definition and coordinate restart must use the same order; the
adapter does not sort atoms. Known definitions are checked using pseudo-atom
atomic numbers. Other monomers remain subject to the installed MBX model's
site count, order and supported interaction combinations.

The current maintenance and branch policy is in [Maintaining the MBX interface](docs/mbx-maintenance.md).
The archived synchronization report records earlier verification and must not
be read as the test receipt for a new build.

## Building

The [installation guide](INSTALLATION.md) walks through cloning the right branch,
creating the compiler environment, building public MBX when needed, installing
the executable and checking it with a calculation.

| Linux x86-64 application | Configure preset | Build directory |
| --- | --- | --- |
| Rigid MBX and force-field calculations, including 2D-TMMC | `raspa-mbx` | `build-mbx-app` |
| Force-field calculations, including 2D-TMMC | `raspa-classical` | `build-classical` |

Both presets are available in a fresh checkout. For developer tests and the
optional AVX2 and server profiles, see [additional build profiles](INSTALLATION.md#additional-build-profiles).

### Use the executable from the build tree

RASPA3 reads a file named `simulation.json` from the current working directory;
it does not take the JSON path as a positional argument. Keep the Conda runtime
libraries active, enter the calculation directory, and run the executable:

```bash
conda activate raspa3-mbx
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /path/to/calculation
/path/to/RASPA3-MBX/build-mbx-app/app/raspa3
```

Useful checks are:

```bash
raspa3 --help
raspa3 --opencl
```

The simulation writes its normal `output/`, restart, and other result files
relative to the calculation directory. An MBX settings path supplied in
`simulation.json` is also resolved relative to that input file.

### Install and expose RASPA3-MBX as a server module

For a shared server installation, install a tested build into a versioned
application directory. Do not expose only the build directory without its
matching Conda runtime libraries:

```bash
cmake --preset raspa-mbx
cmake --build build-mbx-app --target raspa3 --parallel 4
cmake --install build-mbx-app --component app --prefix /shared/apps/raspa3-mbx/REPLACE_WITH_VERSION
cmake --install build-mbx-app --component examples --prefix /shared/apps/raspa3-mbx/REPLACE_WITH_VERSION
```

An Lmod module template is provided at
[`packaging/modulefiles/raspa3-mbx.lua.example`](packaging/modulefiles/raspa3-mbx.lua.example).
Copy it into the server's module tree and edit its `root` and `deps` paths. The
module performs three required operations:

- adds the installed `bin` directory to `PATH`;
- adds the matching Conda `lib` directory to `LD_LIBRARY_PATH`;
- sets `RASPA_DIR` to the installation prefix for force-field lookup. Installed molecule definitions are available to copy from `$RASPA_DIR/share/raspa3/molecules`.

After the administrator installs the module file, users can run:

```bash
module load raspa3-mbx/REPLACE_WITH_VERSION
command -v raspa3
raspa3 --help

cd /path/to/calculation-containing-simulation.json
raspa3
```

For a single-process Slurm calculation, a minimal job body is:

```bash
module purge
module load raspa3-mbx/REPLACE_WITH_VERSION
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
srun --cpu-bind=cores raspa3
```

The module must refer to the same Conda dependency environment used when the
executable was linked. Without that library path, the system can select an
incompatible `libstdc++`, HDF5, OpenBLAS, or OpenMP library and fail before the
simulation starts.

### Classical build without MBX

On this branch, select the `raspa-classical` preset and follow the
[force-field build](INSTALLATION.md#build-and-install-raspa3).

## Running an MBX calculation

`raspa3` reads `simulation.json` from the current working directory.
[Self-contained MBX examples](examples/raspa3_mbx) include water NPT MC,
CO2/Mg-MOF-74 GCMC, conventional Widom zero-loading Qst, energy logging,
coordinate restarts and one-shot snapshot evaluation. Every numbered directory
contains its own force field, molecular definition, MBX settings and any
framework or snapshot file. Each may be copied and run independently.

For example:

```bash
RASPA3_MBX_ROOT=/path/to/RASPA3-MBX
cd "$RASPA3_MBX_ROOT/examples/raspa3_mbx/5_water_npt"
RASPA3_EXECUTABLE="$RASPA3_MBX_ROOT/build-mbx-app/app/raspa3" ./run.sh
```

For adsorption, use `6_gcmc_co2_mg_mof74` instead. The default inputs are short
functional checks. Both directories also include a longer
`simulation.production.json` and instructions for selecting it in a fresh copy.

Add the following fields to a normal Monte Carlo input:

```json
{
  "SimulationType": "MonteCarlo",
  "WriteRestartEvery": 1000,
  "Systems": [
    {
      "UseMBX": true,
      "MBXSettingsFile": "mbx.json",
      "PrintEnergyTerms": true
    }
  ]
}
```

This fragment illustrates the RASPA3-MBX controls; a real input still needs the
normal system, force-field, component, temperature, pressure, and move fields.
Relative `MBXSettingsFile` paths are resolved from the directory containing
`simulation.json`.

### RASPA3-MBX input controls

| Input field | Location | Default | Meaning |
|---|---|---:|---|
| `UseMBX` | Per system | `false` | Use MBX for the supported guest-energy calculation |
| `MBXSettingsFile` | Per system | — | MBX JSON settings file; required when `UseMBX` is true |
| `PrintEnergyTerms` | Per system | `false` | Set true to write accepted-move and conventional-Widom diagnostic CSVs under `output/` |
| `ComputeZeroLoadingHeatOfAdsorption` | Per system | `false` | Accumulate and print (Q_{st}^0) from conventional Widom trials; requires zero guest loading, a positive `WidomProbability`, a rigid host/test component, and fixed volume |
| `WriteRestartEvery` | Top level | `5000` | Periodic coordinate-restart interval for regular MC and thermodynamic integration; `0` disables periodic writes |
| `SimulationType: EnergyEvaluation` | Top level | — | Evaluate a configuration once without changing it |

`WriteEnergyLog` remains a legacy alias for `PrintEnergyTerms`. Supplying both
with different values is an input error.

### Energy diagnostic CSVs

Energy diagnostic CSVs are disabled by default; set `PrintEnergyTerms` to
`true` when move/trial records are needed. Normal scientific summaries,
warnings, errors, and explicitly requested `EnergyEvaluation` results remain
available. Binary restarts that store the logging flag preserve its saved
value; older version-1 system archives without the flag resume with logging
disabled.

With logging enabled, a fresh run truncates `output/energy_terms.s<id>.csv`. A binary-resumed run
appends without duplicating the header. Each row represents an accepted move
and contains the move type, component, loading, retained total, framework-guest
VDW and tail terms, the energy difference, and the acceptance factor. Classical
rows additionally contain guest-guest VDW, framework-guest and guest-guest
charge, and combined Ewald terms. MBX rows instead contain the retained MBX
total plus its 2B, 3B, 4B, dispersion, permanent-electrostatic, and
induced-electrostatic terms. The raw MBX 1B diagnostic is not part of this
historical CSV schema.

Conventional Widom insertions never enter the live system and therefore do not
produce accepted-move rows. When `PrintEnergyTerms` is enabled, each completed
ghost insertion is instead written to
`output/widom_energy_terms.s<id>.csv`. Its row retains the actual loading `N`,
labels the hypothetical loading as `trial_N`, and contains the current and
hypothetical total energies, the hypothetical energy decomposition, insertion
energy `E_insert`, Widom/Rosenbluth `weight`, and `weighted_E_insert =
weight*E_insert`. Failed growths and blocked-pocket trials have zero statistical
weight but no well-defined complete energy decomposition, so they are not CSV
rows.

CSV rows are flushed independently of binary checkpoints. After an abrupt
crash, replaying moves from an older binary checkpoint can therefore duplicate
the CSV rows written after that checkpoint.

### Zero-loading heat of adsorption from Widom trials

Set `ComputeZeroLoadingHeatOfAdsorption` to `true` on a system and give the
rigid component a positive `WidomProbability`. Use a rigid host at fixed volume,
start with zero guest molecules, and do not enable insertion moves that can
change that loading. For example, the relevant parts of a single-component
input are:

```json
{
  "Systems": [
    {
      "ComputeZeroLoadingHeatOfAdsorption": true
    }
  ],
  "Components": [
    {
      "Name": "CO2",
      "CreateNumberOfMolecules": 0,
      "WidomProbability": 1.0
    }
  ]
}
```

RASPA reports block values, an average, and an uncertainty in K and kJ/mol using

```text
Qst^0 = kB*T - <weight * E_insert> / <weight>.
```

The run stops instead of silently relabeling a finite-loading result if any
guest molecule is present initially or appears during sampling. The current
estimator is intentionally restricted to a rigid host and rigid Widom
components at fixed volume; flexible adsorbates require an additional ideal-gas
intramolecular contribution, while flexible/NPT hosts require host-energy or
volume fluctuation terms.

## Evaluating a saved snapshot

Use `EnergyEvaluation` instead of forcing a zero-displacement translation:

```json
{
  "SimulationType": "EnergyEvaluation",
  "ForceField": ".",
  "Systems": [
    {
      "Type": "Framework",
      "Name": "Mg_MOF74_pacman",
      "NumberOfUnitCells": [2, 2, 5],
      "ExternalTemperature": 298.0,
      "ChargeMethod": "Ewald",
      "RestartFileName": "co2_snapshot.json",
      "UseMBX": true,
      "MBXSettingsFile": "mbx.json"
    }
  ],
  "Components": [
    {
      "Name": "co2",
      "MoleculeDefinition": ".",
      "CreateNumberOfMolecules": 0
    }
  ]
}
```

The program prints a human-readable decomposition to standard output and
writes `output/energy_evaluation.json`. The JSON contains the simulation box,
molecule counts, the RASPA energy decomposition in kelvin, and all seven raw
MBX terms in kcal/mol when MBX is enabled.

This mode reads lightweight coordinate-restart JSON, not a binary driver
checkpoint. It requires whole-molecule coordinates and rejects generated extra
molecules, CFCMC/fixed-lambda state, and ambiguous component mappings. For a
framework system, the framework coordinates and cell come from the configured
CIF; a coordinate snapshot containing `SimulationBox` is rejected because the
JSON does not contain matching framework coordinates. See
[restart and snapshot evaluation](docs/manual/restart.md#evaluating-the-energy-of-a-coordinate-snapshot)
for the full format and restrictions.

## Morse-potential inputs

Morse interactions are selected in `force_field.json`, for example:

```json
{
  "names": ["Mg", "Cco2"],
  "type": "morse",
  "parameters": [19.7068, 1.3677, 4.5000],
  "source": "Fit"
}
```

The retained MOF-74 force field contains explicit Mg-CO2 Morse pairs. Use the
direct CPU energy routes documented in the synchronization report; do not
assume that every interpolation or OpenCL path implements a Morse override.

## Testing

After configuring and building the [developer test profile](INSTALLATION.md#additional-build-profiles), run its suite with:

```bash
ctest --preset mbx-local
```

Run only the MBX and new diagnostic controls with:

```bash
ctest --preset mbx-local \
  -R '^(mbx_static_energy|mbx_energy_log|run_control|energy_evaluation)' \
  --output-on-failure
```

The combined branch validation on `e6659c5d` is recorded in the [hosted build and calculation results](https://github.com/Chung-Research-Group/RASPA3-MBX/actions/runs/44571393370). Earlier MBX-specific repair results are recorded in
[Maintained verification](docs/validation-current.md). Earlier synchronization
counts are preserved as historical evidence, not as this edition's results.

## Keeping the fork synchronized

The daily release check stages each new published stable release in
`sync-to-upstream/<tag>`. After validation with MBX enabled and disabled,
the job advances `main` and `mbx` together. The independent `2d-tmmc`
branch retains its release base. See the [release synchronization guide](docs/mbx-upstream-sync.md#published-release-synchronization-19-november-2026).

## Documentation

- [RASPA3-MBX synchronization and compatibility report](docs/mbx-upstream-sync.md)
- [Public MBX replacement and validation report](docs/mbx-public-migration.md)
- [Input-command reference](docs/manual/commands.md)
- [Coordinate and binary restarts](docs/manual/restart.md)
- [Full RASPA3 manual](docs/manual/manual.md)
- [Official upstream RASPA3 documentation](https://iraspa.github.io/RASPA3/)

## Citation and credit

When using this fork, cite RASPA3:

Y. A. Ran, S. Sharma, S. R. G. Balestra, Z. Li, S. Calero, T. J. H. Vlugt,
R. Q. Snurr, and D. Dubbeldam, “RASPA3: A Monte Carlo code for computing
adsorption and diffusion in nanoporous materials and thermodynamics properties
of fluids,” *Journal of Chemical Physics* **161**, 114106 (2024),
[doi:10.1063/5.0226249](https://doi.org/10.1063/5.0226249).

Also cite the MBX/MB-nrg model papers appropriate to the species and parameter
set used by your `mbx.json`, including:

M. Riera, C. Knight, E. F. Bull-Vulpe, X. Zhu, H. Agnew, D. G. A. Smith,
A. C. Simmonett, and F. Paesani, “MBX: A many-body energy and force calculator
for data-driven many-body simulations,” *Journal of Chemical Physics* **159**,
054802 (2023),
[doi:10.1063/5.0156036](https://doi.org/10.1063/5.0156036).

RASPA3 authors and contributors remain credited in the upstream repository and
source history.

## License

RASPA3-MBX source code is distributed under the [MIT License](LICENSE.md).
MBX is a separately obtained dependency with separate licensing terms; this
repository does not relicense MBX.
