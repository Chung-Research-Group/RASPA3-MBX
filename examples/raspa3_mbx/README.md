# Self-contained RASPA3-MBX examples

Each numbered directory contains **all simulation inputs** and its own launcher.
You may copy any one directory elsewhere and run it without its parent directory,
other examples, an installed force-field database, or downloading a structure.
The executable and its runtime libraries are the only external requirements.
Use a RASPA3-MBX build with `BUILD_MBX=ON`, built against the maintained
[MBX v1.4.0 dependency](../../README.md#building).

| Directory | Calculation |
| --- | --- |
| `1_energy_logging_and_restarts` | CO2/MFI: accepted-move energy CSV and coordinate restarts |
| `2_widom_zero_loading_qst` | CO2/MFI: conventional Widom trials and zero-loading heat of adsorption |
| `3_snapshot_energy_evaluation` | CO2/MFI: one-shot energy evaluation of a supplied snapshot |
| `4_gcmc_co2_adsorption` | CO2/MFI: conventional insertion/deletion GCMC |
| `5_water_npt` | Rigid MB-pol water: NPT MC with translation, rotation, CBMC reinsertion and volume moves |
| `6_gcmc_co2_mg_mof74` | CO2/Mg-MOF-74: GCMC with CBMC insertion/deletion and reinsertion |

## Run directly in an example

After activating your installed RASPA3-MBX environment:

```bash
cd 5_water_npt
./run.sh
```

If `raspa3` is not on `PATH`, select an executable explicitly:

```bash
RASPA3_EXECUTABLE=/absolute/path/to/raspa3 ./run.sh
```

You may also invoke `raspa3` directly in the example directory. Every input
path is local. The launcher fixes `RASPA_DIR` to that same directory, resolves
the executable before changing directory, and defaults `OMP_NUM_THREADS` to 1.
It does not download, copy or regenerate missing inputs. Each local README
explains the files, ensemble, moves and outputs.

The parent launcher remains a convenience, for example:

```bash
./run-example.sh 6_gcmc_co2_mg_mof74
```

## Short examples and longer calculations

The default inputs are **functional examples**, not converged calculations.
Water starts with 16 molecules in a 24 Å cube to keep the execution check short.
It is not an equilibrated liquid. The Mg-MOF-74 example uses the manuscript's
2 × 2 × 5 framework at 298 K and 1000 Pa but only a few cycles.

Directories 5 and 6 also include `simulation.production.json`: the longer
512-water NPT and CO2/Mg-MOF-74 GCMC setups from the manuscript input bundle.
In a fresh copy of the desired example directory, use
`cp simulation.production.json simulation.json` before running. All other
required files are already present. These calculations can be expensive and
require independent equilibration and sampling checks.

MBX identifiers are case-sensitive. The included rigid component definitions
use `h2o` in O,H,H order and `co2` in C,O,O order. MBX supplies guest many-body
energy and guest–guest/framework–guest electrostatics and polarization.
The retained classical framework contribution is framework–guest vdW and its
tail correction. Framework charges are fixed, nonpolarizable input charges.

## Verify directory portability

From the repository root:

```bash
python tests/run_mbx_examples.py --source . \
  --executable /absolute/path/to/raspa3 \
  --output /absolute/path/to/new-example-checks
```

The checker copies **each example directory alone**, runs its unmodified short
input from that directory and checks native outputs and move statistics.
The output path must be new. Passing these checks establishes functional
execution and local input closure, not ensemble convergence.
