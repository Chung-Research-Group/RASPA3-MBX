# Conventional Widom insertion at zero loading

An empty MFI framework is sampled with conventional Widom test insertions
(5 initialization and 10 production cycles). Ghost molecules are not retained
in the system. `ComputeZeroLoadingHeatOfAdsorption` enables the zero-loading
heat-of-adsorption estimator; `PrintEnergyTerms` writes trial energies to
`output/widom_energy_terms.s0.csv`. The final report contains the resulting
zero-loading estimate. The few test insertions here are insufficient for a
reliable numerical estimate.

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
