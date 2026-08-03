# Restart-files `RASPA`
\page restart Restart Files

## Writing restart-files

Regular Monte Carlo and thermodynamic-integration coordinate restart files are
automatically written and, by default, updated every 5000 completed cycles in
each stage. Set the top-level `"WriteRestartEvery"` value to choose another
interval; `0` disables periodic writes, while the drivers' existing initial and
final snapshots are still written. Consider the following Gibbs-ensemble simulation:
```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfProductionCycles" : 25000,
  "NumberOfInitializationCycles" : 10000,
  "PrintEvery" : 1000,
  "WriteRestartEvery" : 2000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 240.0,
      "ChargeMethod" : "Ewald",
      "GibbsVolumeMoveProbability" : 0.01
    },
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 240.0,
      "ChargeMethod" : "Ewald",
      "GibbsVolumeMoveProbability" : 0.01
    }
  ],

  "Components" :
  [
    {
      "Name" : "CO2",
      "MoleculeDefinition" : "ExampleDefinitions",
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "GibbsSwapCBMCProbability" : 1.0,
      "CreateNumberOfMolecules" : [256, 256]
    }
  ]
}
```
At the end of the simulation you will find in the `output`-directory the files
```
output/restart_240_0.s0.json
output/restart_240_0.s1.json
```


----------------------------------------------------------------------------------

## Using restart-files

To start from the saved (and hopefully equilibrated) positions, point each
system's `"RestartFileName"` at the corresponding file. Relative paths are
resolved from the directory containing `simulation.json`, so the files do not
need to be copied out of `output/`. Because a new run can overwrite those
paths, copy them to a separately named snapshot if they must be preserved.

The contents of the restart files look like
```
{
  "CO2": [
    [ 16.640068355882505, 9.764919807611305, 4.333198721377549 ],
    [ 17.15733718188185, 8.82872907970449, 4.75293778631024 ],
    [ 17.67460600788119, 7.892538351797674, 5.172676851242931 ],
    [ 3.815549021105568, 8.879311693833092, 17.616704775297364 ],
    [ 4.677375964586537, 8.549547910646638, 18.30132961182377 ],
    [ 5.539202908067507, 8.219784127460183, 18.985954448350178 ]
  ],
  "SimulationBox": {
    "angle-alpha": 90.0,
    "angle-beta": 90.0,
    "angle-gamma": 90.0,
    "length-a": 31.881622110825266,
    "length-b": 31.881622110825266,
    "length-c": 31.881622110825266
  }
}

```
So, the component-name as the key, and the value is an array of positions,
and information on the simulation-box (during NPT or Gibbs the box changes).
CO<sub>2</sub> molecules are usually modeled as `rigid` and implemented using quaternions.
```json
{
  "CriticalTemperature": 304.1282,
  "CriticalPressure": 7377300.0,
  "AcentricFactor": 0.22394,
  "pseudoAtoms": [
    ["O_co2", [0.0, 0.0, 1.149]],
    ["C_co2", [0.0, 0.0, 0.0]],
    ["O_co2", [0.0, 0.0, -1.149]]
  ]
}
```
The quaternions are computed from the positions using singular-value decompositions.
The advantage is that, even when the positions are input by hand and do not satisfy
the rigid constraints, an optimal mapping or best fit is obtained.
For output restart-files, reading these in, should result in the exact same
positions and energies.



```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfProductionCycles" : 25000,
  "NumberOfInitializationCycles" : 1000,
  "PrintEvery" : 1000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 240.0,
      "ChargeMethod" : "Ewald",
      "GibbsVolumeMoveProbability" : 0.01,
      "RestartFileName" : "output/restart_240_0.s0.json"
    },
    {
      "Type" : "Box",
      "BoxLengths" : [30.0, 30.0, 30.0],
      "ExternalTemperature" : 240.0,
      "ChargeMethod" : "Ewald",
      "GibbsVolumeMoveProbability" : 0.01,
      "RestartFileName" : "output/restart_240_0.s1.json"
    }
  ],

  "Components" :
  [
    {
      "Name" : "CO2",
      "MoleculeDefinition" : "ExampleDefinitions",
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 0.5,
      "GibbsSwapCBMCProbability" : 1.0,
      "CreateNumberOfMolecules" : [0, 0]
    }
  ]
}
```
Note because the positions are already equilibrated in this case, you can lower the amount of initialization cycles.
This is the main point of using restart-file. But if you use the files to run at different conditions (for example
a different temperature) you would still need a significant amount of initialization cycles.

The number of `CreateNumberOfMolecules` should be put to zero, since the molecules are read from the restart-files.
For a normal simulation, additional molecules may still be created on top of
the loaded coordinates. This can be useful for constructing high-density
systems. For a one-shot energy comparison, additional random molecules would
change the snapshot and are therefore rejected, as described below.

----------------------------------------------------------------------------------

## Evaluating the energy of a coordinate snapshot

Use `"SimulationType" : "EnergyEvaluation"` to load a coordinate restart,
evaluate it once, and exit without changing atomic positions or attempting a
Monte Carlo move:

```json
{
  "SimulationType" : "EnergyEvaluation",
  "ForceField" : ".",
  "Systems" : [
    {
      "Type" : "Framework",
      "Name" : "Mg_MOF74_pacman",
      "NumberOfUnitCells" : [2, 2, 5],
      "ExternalTemperature" : 298.0,
      "ChargeMethod" : "Ewald",
      "RestartFileName" : "snapshots/co2_state.json",
      "UseMBX" : true,
      "MBXSettingsFile" : "mbx.json"
    }
  ],
  "Components" : [
    {
      "Name" : "co2",
      "MoleculeDefinition" : ".",
      "CreateNumberOfMolecules" : 0
    }
  ]
}
```

The human-readable decomposition is printed to standard output. The stable,
machine-readable result is written to `output/energy_evaluation.json` and
contains the simulation box, molecule counts, the full RASPA energy
decomposition in kelvin, and—when MBX is enabled—the seven raw MBX terms in
kcal/mol. The MBX 1B term is reported but explicitly excluded from the retained
RASPA3-MBX total.

This mode reads the lightweight JSON coordinate format, not a binary driver
checkpoint. The JSON format stores whole guest-molecule atom positions and, for
a `"Box"` system, its simulation box. A `"Framework"` system obtains its cell
from the matching CIF/unit-cell definition; a framework restart that also
declares `"SimulationBox"` is rejected because no matching framework coordinates
are stored. Component keys must match the configured names one-to-one, including
empty arrays for zero-molecule components. The format does not store velocities,
random-number state, statistics, fractional-molecule state, or flexible-framework
coordinates, so CFCMC/fixed-lambda configurations are rejected. Use a binary
restart for exact trajectory continuation, and use `EnergyEvaluation` for a
fixed-coordinate potential-energy comparison.

----------------------------------------------------------------------------------

## Using binary-restart-files

Binary restart-files are usefull in case of system crashes or limits on running-times in computer clusters.
By continuing from the binary restart-file you would get identical results are running the long job, i.e.
it contains the entire state of the simulation. Note that therefore whatever you change in the `simulation.json` 
is ignored since it uses the state from the binary restart-file.

During the simulation a file named `restart_data.bin` is written out.
The option to continue from this file is
```json
{
  "RestartFromBinaryFile" : true
}
```
Note not everything can be recovered, like written pdb-movies.
The setting
```json
{
  "WriteBinaryRestartEvery" : 5000
}
```
controls how often the file is written out. The default is every 5000 cycles.

Binary restart files are version- and layout-sensitive. Use the same RASPA3-MBX
build (ideally the same commit) to resume them; compatibility with another
RASPA release or fork is not guaranteed. Coordinate JSON is more portable for
whole-molecule positions, but it does not preserve the complete simulation
state described above.
