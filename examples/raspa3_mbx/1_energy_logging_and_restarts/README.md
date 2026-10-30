# Energy logging and coordinate restarts

A fixed-population CO2/MFI Monte Carlo example with translation and rotation
moves (one molecule; 2 initialization and 5 production cycles).
`PrintEnergyTerms` writes accepted-move diagnostics to
`output/energy_terms.s0.csv`. `WriteRestartEvery: 5` requests periodic
coordinate restart snapshots. Both controls can be disabled in
`simulation.json` when these diagnostics are unnecessary.

## Run

This directory contains every simulation input. Only a RASPA3-MBX executable
built with `BUILD_MBX=ON` and its installed runtime libraries are required.
Run directly here:

```sh
RASPA3_EXECUTABLE=/absolute/path/to/raspa3 ./run.sh
```

If `raspa3` is on `PATH`, `./run.sh` is sufficient. The launcher changes to
this directory, sets `RASPA_DIR` here and defaults to one OpenMP thread. It
preserves an explicitly set `OMP_NUM_THREADS` and does not rewrite any input.
You may also run the executable directly from this directory.

The local files are `simulation.json`, `mbx.json`, `co2.json`,
`force_field.json` and `MFI_SI.cif`. CO2 sites are ordered C,O,O, as required
by the case-sensitive MBX monomer name `co2`. Framework and force-field
parameters are copied from the repository's MFI/CO2 example; their sources
are retained in `force_field.json`. This is an illustrative MFI calculation,
not the Mg-MOF-74 calculation in the manuscript.

These intentionally short runs demonstrate functionality, not equilibrated
or converged adsorption properties. Archive an existing `output/` directory
before rerunning if you need to retain its results.
