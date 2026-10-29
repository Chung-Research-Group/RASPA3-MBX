# CO2 adsorption in Mg-MOF-74: GCMC with MBX

Run directly in this directory with an MBX-enabled RASPA3 installation
(maintained dependency: MBX v1.4.0):

```bash
./run.sh
```

The launcher selects `RASPA3_EXECUTABLE` or `raspa3` from `PATH`:

```bash
RASPA3_EXECUTABLE=/absolute/path/to/raspa3 ./run.sh
```

With the executable and runtime libraries available, `raspa3` also works
directly here. There are no input dependencies on sibling example directories,
an installed structure/force-field database, or a restart file.

## Included inputs

| File | Purpose |
| --- | --- |
| `simulation.json` | Short GCMC run at 298 K and 1000 Pa (0.01 bar) |
| `Mg_MOF74_pacman.cif` | Framework geometry and fixed PACMAN/DDEC6 atomic charges |
| `co2.json` | Rigid MBX `co2` monomer in required C,O,O site order |
| `force_field.json` | Fitted framework–CO2 LJ pairs and classical proposal parameters |
| `mbx.json` | MBX settings, with CO2 three-body terms enabled |
| `simulation.production.json` | Longer manuscript setup at the same state point |

Both inputs use a fixed **2 × 2 × 5 framework supercell**, an initially empty
pore, and fugacity coefficient 1. Translation, rotation and CBMC reinsertion
have relative weights 1 each. CBMC insertion/deletion has weight 4
(`SwapProbability`). Conventional insertion/deletion is disabled. The
population can change, while framework coordinates and cell remain fixed.

The default has two initialization, two equilibration and five production
cycles. It demonstrates input loading and move execution, **not an equilibrated
isotherm or converged heat of adsorption**. The framework charges are read from
the included CIF. MBX supplies guest–guest many-body interactions and
guest–guest/framework–guest electrostatics and polarization; the retained
classical contribution is framework–guest vdW plus its tail correction.
Framework atoms carry fixed charges and are not polarizable MBX monomers.

## Longer calculation

In a **fresh copy of this directory**:

```bash
cp simulation.production.json simulation.json
./run.sh
```

This selects the supplied manuscript setup with 10000 initialization and 20000
production cycles at 298 K and 1000 Pa. Allow for its much longer runtime and
establish equilibration independently. The files are curated from the
manuscript's complete fitted-LJ/MBX input bundle without changing numerical
force-field values or framework charges.

## Output

Native text/JSON reports and coordinate snapshots appear in `output/`.
The short default also opts into `output/energy_terms.s0.csv`, recording
accepted moves with their energy decomposition, and writes a density grid
and molecular movie. Inspect the final native move statistics for attempts
and acceptance. A few cycles do not provide reliable adsorption averages
or uncertainties. Run each calculation in a fresh directory.
