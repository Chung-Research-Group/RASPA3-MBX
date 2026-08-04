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

This is a research fork, not an official iRASPA release. The current
integration uses the RASPA3 3.1.0 source layout and contains the audited
upstream history through commit `a8a70ca7`. The exact merge decisions,
scientific energy convention, verification record, and future update process
are documented in [the synchronization report](docs/mbx-upstream-sync.md).

## What this fork adds

| Capability | RASPA3-MBX behavior |
|---|---|
| MBX guest energy | Enabled per system with `"UseMBX": true` and an `MBXSettingsFile` |
| Morse potential | Supported by the current native RASPA3 force-field implementation and retained in the MOF-74 inputs |
| Energy-term logs | `PrintEnergyTerms` writes accepted moves to `output/energy_terms.s<id>.csv` and conventional Widom trials to `output/widom_energy_terms.s<id>.csv` |
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

MBX trajectory generation currently supports the regular, single-system
`MonteCarlo` driver. Supported state-changing moves include translation,
rotation, CBMC reinsertion, conventional or CBMC insertion/deletion, and a
framework-free isotropic volume move without an external field. Conventional
insertion/deletion is limited to rigid components; flexible components must use
CBMC. Widom evaluation and the non-mutating `EnergyEvaluation` mode are also
supported.

The input reader rejects unported MBX combinations before starting a run. In
particular, do not enable MBX for molecular dynamics, minimization,
thermodynamic integration, CFCMC, Gibbs/reaction moves, transition-matrix or
parallel-replica drivers, smart/gradient moves, anisotropic volume changes, or
flexible-framework volume changes. The full support matrix and pressure/virial
limitations are in
[Supported MBX execution surface](docs/mbx-upstream-sync.md#supported-mbx-execution-surface).

Classical RASPA3 calculations remain available by building with
`BUILD_MBX=OFF` or by leaving `UseMBX` false.

## Building

### Requirements

The current Linux build requires:

- CMake 4.0.3 through 4.4.x, Ninja, and a compiler with the required C++26
  module support; the maintained build was tested with Clang 20. `import std`
  is still an experimental CMake feature whose activation token can change
  between CMake releases, so newer CMake versions must be reviewed before use;
- BLAS, LAPACK, OpenMP, OpenCL, zlib, and HDF5;
- for an MBX-enabled build, MBX headers, the generated `mbx_version.h`, the MBX
  library, and FFTW3.

RASPA3-MBX and MBX must use a compatible OpenMP runtime. For example, do not
mix an MBX library built with LLVM `libomp` into a RASPA3 executable configured
for GNU `libgomp`.

MBX is an external dependency and is not vendored here. The maintained public
dependency is [paesanilab/MBX](https://github.com/paesanilab/MBX), pinned to
release `v1.4.0` at commit
`0e01b75b47611d7d51f27a34b112bdc5e2090a50`. It replaces the former private
MBX 1.3.3 fork: the private fork contained no private potential, coefficient,
or RASPA API, and its one compile-only OpenMP patch is superseded in v1.4.0.
Compatibility with arbitrary moving MBX branches is not implied. See the
[public-MBX migration audit](docs/mbx-public-migration.md) for the source,
API, numerical, and licensing comparison.

This is a source-compatible migration, not a binary hot-swap. Do not relink
old RASPA objects against the new archive: public v1.4.0 also updates MBX's
bundled JSON implementation. Reconfigure from a fresh RASPA build tree using
the public headers and library together.

Record the exact MBX source revision alongside each RASPA3-MBX release. MBX
also has licensing terms separate from this repository's MIT license, so
consult the license supplied with your MBX build.

`BUILD_MBX` defaults to `ON` in this fork. Pass `-DBUILD_MBX=OFF` explicitly
when configuring a classical-only build without the external MBX dependency.

### Build the audited public MBX dependency

Use MBX's basic, non-MPI installation. The following is the validated
Clang/Conda workflow; choose a permanent versioned path for `MBX_ROOT`:

```bash
git clone https://github.com/paesanilab/MBX.git
cd MBX
git checkout --detach 0e01b75b47611d7d51f27a34b112bdc5e2090a50

export MBX_ROOT=/path/to/dependencies/mbx-1.4.0-0e01b75b
export FFTW_HOME="$CONDA_PREFIX"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

autoreconf -fi
CC="$CONDA_PREFIX/bin/clang" \
  CXX="$CONDA_PREFIX/bin/clang++" \
  LDFLAGS="-L$CONDA_PREFIX/lib" \
  ./configure --prefix="$MBX_ROOT"

make -j8
make check
make install
```

Passing the Conda library directory in `LDFLAGS` during `configure` is
important. In the maintained Clang environment it lets MBX select the
`libiomp5` compatibility name for LLVM `libomp`, instead of falling back to GNU
`libgomp`. Verify the installed executable before building RASPA3-MBX:

```bash
ldd "$MBX_ROOT/bin/single_point" | grep -E 'fftw|[ig]?omp'
```

The maintained build shows Conda `libfftw3.so` and one `libomp.so`; it does not
load both `libgomp` and `libomp`. All 34 MBX v1.4.0 tests must pass before the
installation is used.

### Build RASPA3-MBX against MBX

The most reproducible local setup is an ignored `CMakeUserPresets.json` that
inherits the tracked `linux_conda` preset and supplies machine-specific MBX,
FFTW, and matching OpenMP paths. Start from the provided template; absolute
dependency paths should not be committed:

```bash
cp CMakeUserPresets.json.example CMakeUserPresets.json
export MBX_ROOT=/path/to/installed/mbx

cmake --preset mbx-local --fresh
cmake --build --preset mbx-local --parallel
ctest --preset mbx-local
```

Activate the Conda dependency environment first so `CONDA_PREFIX` is defined.
The template expects `libmbx.a` under `$MBX_ROOT/lib`; edit the ignored local
preset if your installation layout differs. It uses LLVM `libomp`, matching the
maintained build. Select a different but consistent compiler/OpenMP combination
when rebuilding both RASPA3-MBX and MBX together.

Before configuring, confirm that the active `cmake` is the one from the Conda
environment. On the maintained workstation, `/usr/local/bin/cmake` is an older
3.28 installation and must not be used for this source tree:

```bash
conda activate raspa3-mbx_2
hash -r
command -v cmake
cmake --version
```

The user-preset template provides three roles without creating separate source
repositories:

| Preset | Build directory | Intended use |
|---|---|---|
| `mbx-local` | `build-mbx` | Portable validation build with focused MBX tests |
| `mbx-release-avx2` | `build-mbx-release-avx2` | Production build for the current AVX2 workstation |
| `mbx-server` | `build-mbx-server` | Portable x86-64 server build without workstation-specific `-march` flags |

For the optimized workstation executable, use:

```bash
cmake --preset mbx-release-avx2 --fresh
cmake --build --preset mbx-release-avx2 --parallel
```

For a server, clone this same repository, install the matching Conda environment
and audited MBX dependency on that server, set `MBX_ROOT`, copy the user-preset
template, and run:

```bash
cmake --preset mbx-server --fresh
cmake --build --preset mbx-server --parallel
```

Build on the server itself unless its compiler, standard library, and CPU are
known to match the workstation. C++ module artifacts are toolchain-specific;
copying the source repository and rebuilding is safer than copying a build
directory. The old `MBX_INCLUDES` spelling is intentionally not used: this
version expects `MBX_INCLUDE_DIR`, `MBX_CONFIG_INCLUDE_DIR`, and `MBX_LIBRARY`.

If configuration reports that `CXX_MODULE_STD` lacks toolchain support, first
run the configure command with `--fresh`. Then verify the CMake version as
shown above. RASPA3-MBX contains the experimental activation tokens for CMake
4.0.3--4.4.x; a version outside that range now stops with a direct diagnostic
instead of failing later while generating `raspakit_base`.

The main executable is:

```text
build-mbx/app/raspa3
```

### Create the dependency environment

The tracked [`env.yml`](env.yml) defines the validated Linux x86-64 direct
dependencies: Clang 20, CMake 4.0.3, OpenBLAS, HDF5, FFTW, the OpenCL loader,
and LLVM OpenMP. It also includes Autotools and related utilities needed to
build the audited public MBX dependency:

```bash
conda env create -f env.yml
conda activate raspa3-mbx
```

The compiler and CMake versions are intentionally pinned because C++ standard
library module artifacts are compiler-specific and CMake's `import std` gate is
version-sensitive. The YAML has no `prefix`, so each user can install it in a
normal writable Conda environment. It is intended for Linux x86-64 servers;
other operating systems need a platform-specific environment.

MBX itself is not installed by this YAML because the audited MBX revision is
not distributed as the required Conda package. Build and install public MBX
with this activated environment using the preceding instructions, then set
`MBX_ROOT` to that installation before configuring RASPA3-MBX.

`ocl-icd` supplies the OpenCL loader. A server still needs an appropriate
vendor OpenCL implementation if OpenCL calculations are requested; the NVIDIA,
AMD, or Intel driver is normally managed outside Conda by the cluster.

### Use the executable from the build tree

RASPA3 reads a file named `simulation.json` from the current working directory;
it does not take the JSON path as a positional argument. Keep the Conda runtime
libraries active, enter the calculation directory, and run the executable:

```bash
conda activate raspa3-mbx
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /path/to/calculation
/path/to/RASPA3-MBX/build-mbx/app/raspa3
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
cmake --build --preset mbx-server --parallel
cmake --install build-mbx-server --prefix /shared/apps/raspa3-mbx/REPLACE_WITH_VERSION
```

An Lmod module template is provided at
[`packaging/modulefiles/raspa3-mbx.lua.example`](packaging/modulefiles/raspa3-mbx.lua.example).
Copy it into the server's module tree and edit its `root` and `deps` paths. The
module performs three required operations:

- adds the installed `bin` directory to `PATH`;
- adds the matching Conda `lib` directory to `LD_LIBRARY_PATH`;
- sets `RASPA_DIR` so installed force fields and molecule definitions can be
  found under `share/raspa3`.

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

Use a separate disposable build directory:

```bash
cmake --preset linux_conda -B build-classical \
  -DBUILD_MBX=OFF \
  -DBUILD_TESTING_MBX=OFF

cmake --build build-classical --parallel
```

This build must not require MBX headers or libraries.

## Running an MBX calculation

`raspa3` reads `simulation.json` from the current working directory. A complete
MOF-74/CO2 example with MBX and Morse interactions is stored under
[`calculations/MOF74-CO2/morse`](calculations/MOF74-CO2/morse).
Short, runnable examples for energy logging, restart cadence, conventional
Widom zero-loading Qst, and one-shot snapshot evaluation are collected under
[`examples/raspa3_mbx`](examples/raspa3_mbx).

For example:

```bash
RASPA3_MBX_ROOT=/path/to/RASPA3-MBX
cd "$RASPA3_MBX_ROOT/calculations/MOF74-CO2/morse/calculation/Mg_MOF74_pacman/298K/1e5Pa"
"$RASPA3_MBX_ROOT/build-mbx/app/raspa3"
```

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
| `PrintEnergyTerms` | Per system | `true` | Write accepted-move and conventional-Widom diagnostic CSVs under `output/`; set false to disable |
| `ComputeZeroLoadingHeatOfAdsorption` | Per system | `false` | Accumulate and print (Q_{st}^0) from conventional Widom trials; requires zero guest loading, a positive `WidomProbability`, a rigid host/test component, and fixed volume |
| `WriteRestartEvery` | Top level | `5000` | Periodic coordinate-restart interval for regular MC and thermodynamic integration; `0` disables periodic writes |
| `SimulationType: EnergyEvaluation` | Top level | — | Evaluate a configuration once without changing it |

`WriteEnergyLog` remains a legacy alias for `PrintEnergyTerms`. Supplying both
with different values is an input error.

### Energy diagnostic CSVs

A fresh run truncates `output/energy_terms.s<id>.csv`. A binary-resumed run
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

Run the configured suite with:

```bash
ctest --preset mbx-local
```

Run only the MBX and new diagnostic controls with:

```bash
ctest --preset mbx-local \
  -R '^(mbx_static_energy|mbx_energy_log|run_control|energy_evaluation)' \
  --output-on-failure
```

The latest audited results, including known failures in unchanged upstream
code, are recorded under
[Verification performed](docs/mbx-upstream-sync.md#verification-performed).

## Keeping the fork synchronized

Use one active checkout with two remotes:

```text
origin    https://github.com/Chung-Research-Group/RASPA3-MBX.git
upstream  https://github.com/iRASPA/RASPA3.git
```

Do not make a new directory for each upstream update, and do not merge
`upstream/main` blindly into the release branch. Use a temporary integration
branch:

```bash
git fetch upstream
git switch main
git switch -c sync/upstream-YYYYMMDD-SHA
git merge --no-ff upstream/main
```

Resolve build architecture in favor of upstream, then reapply and review the
small MBX interfaces. Build with MBX both enabled and disabled, run the MBX,
Morse, snapshot, logger, and upstream tests, and only then merge the sync branch
into `main` and create a release tag. The detailed process and CI matrix are in
[Long-term branch and update plan](docs/mbx-upstream-sync.md#long-term-branch-and-update-plan).

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
