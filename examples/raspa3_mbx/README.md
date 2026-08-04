# RASPA3-MBX feature examples

These short inputs exercise the controls added by the RASPA3-MBX fork. The
inputs reuse the MFI framework and force field from
`examples/basic/10_mc_adorption_co2_in_mfi`. The lowercase `co2` molecule
definition is kept here because MBX monomer names are case-sensitive. This
keeps the cases small while allowing them to install with the standard example
collection.

Build RASPA3-MBX first, or load its server module. Then run an example from
this directory:

```bash
./run-example.sh 1_energy_logging_and_restarts
./run-example.sh 2_widom_zero_loading_qst
./run-example.sh 3_snapshot_energy_evaluation
```

The launcher selects the executable in this order:

1. `RASPA3_EXECUTABLE=/absolute/path/to/raspa3`
2. `raspa3` from `PATH`, including a loaded environment module
3. `build-mbx-release-avx2/app/raspa3` or `build-mbx/app/raspa3` in a source checkout

For example, on a module-based server:

```bash
module purge
module load raspa3-mbx/VERSION
./run-example.sh 2_widom_zero_loading_qst
```

For a local build:

```bash
RASPA3_EXECUTABLE=/path/to/RASPA3-MBX/build-mbx/app/raspa3 \
  ./run-example.sh 1_energy_logging_and_restarts
```

## 1. Energy logging and coordinate restarts

`1_energy_logging_and_restarts/simulation.json` enables:

```json
"WriteRestartEvery": 5
```

and, on the system:

```json
"PrintEnergyTerms": true
```

Accepted translations and rotations are written to
`output/energy_terms.s0.csv`. Coordinate snapshots are written every five
completed cycles. Set `PrintEnergyTerms` to `false` to disable diagnostic CSV
output, or set `WriteRestartEvery` to `0` to disable periodic coordinate JSON
snapshots.

## 2. Conventional Widom energy log and zero-loading Qst

`2_widom_zero_loading_qst/simulation.json` starts with zero CO2 molecules and
uses only conventional Widom trials. It enables:

```json
"ComputeZeroLoadingHeatOfAdsorption": true
```

The live system remains at zero loading. Completed ghost insertions are
written to `output/widom_energy_terms.s0.csv`, not the accepted-move CSV. The
final report prints block and averaged values for

```text
Qst^0 = kB*T - <weight*E_insert>/<weight>
```

in K and kJ/mol. The short run is intended as a functional example; increase
the initialization and production cycles for publishable statistics.

## 3. One-shot snapshot energy evaluation

`3_snapshot_energy_evaluation/simulation.json` reads the fixed coordinates in
`co2_snapshot.json`, evaluates them once, and exits without attempting a Monte
Carlo move. It writes the machine-readable decomposition to
`output/energy_evaluation.json` and prints a human-readable decomposition to
standard output.

The example uses MBX for guest energy. The framework-guest interaction remains
the configured classical MFI/CO2 interaction, which is the RASPA3-MBX energy
contract. The retained MOF-74 example with explicit Morse framework-guest pairs
is available under `calculations/MOF74-CO2/morse` in a source checkout.

## Shared requirements

All three inputs set `UseMBX` to `true` and read `mbx.json` from this directory.
They require an executable built with `BUILD_MBX=ON`. Output is created inside
the selected example directory. Remove or archive an old `output/` directory
before comparing fresh runs.
