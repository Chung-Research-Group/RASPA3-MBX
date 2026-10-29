# Rigid-water NPT Monte Carlo with MBX

Run from this directory after installing RASPA3-MBX with `BUILD_MBX=ON`
against the supported MBX v1.4.0 dependency:

```bash
./run.sh
```

The launcher uses `RASPA3_EXECUTABLE` when set, otherwise `raspa3` from `PATH`.
For example:

```bash
RASPA3_EXECUTABLE=/absolute/path/to/raspa3 ./run.sh
```

With the executable and its runtime libraries available, `raspa3` also works
directly here. No other example, installed force-field database, downloaded
structure, or coordinate restart is needed.

## Included inputs

| File | Purpose |
| --- | --- |
| `simulation.json` | Short functional NPT run at 298 K and 100000 Pa |
| `h2o.json` | Rigid water geometry in MBX site order O,H,H |
| `force_field.json` | Fixed classical parameters used for trial proposals |
| `mbx.json` | MBX settings, including 9 Å/5 Å two-/three-body cutoffs |
| `simulation.production.json` | Longer 512-water setup from the manuscript input bundle |

The short default creates **16 waters in a 24 Å cube**, with two initialization,
two equilibration and five production cycles. It starts from dilute water to
make a quick, self-contained execution check. It is **not an equilibrated liquid
or a density benchmark**. Translation, rotation, CBMC reinsertion and isotropic
volume changes are enabled. No swap move is enabled, so the population is fixed.

The case uses MBX `h2o` (original MB-pol) with rigid molecular geometry. It is
not MB-pol(2023), flexible MB-pol MD, or a fluctuating-charge water simulation.
The local classical proposal parameters are fixed: oxygen epsilon/kB = 135 K,
sigma = 3.176 Å, qO = -0.672 e and qH = +0.336 e. These are the parameters in the
project's supplied water inputs. The historical SPC-FQ label does not mean
that Rick's charge-equilibration algorithm is implemented here. MBX supplies
the retained intermolecular energy, including electrostatics and polarization.

## Longer calculation

In a **fresh copy of this directory**, select the supplied 512-water input:

```bash
cp simulation.production.json simulation.json
./run.sh
```

This keeps the 24 Å initial box, 298 K, 1 bar, 10000 initialization and 10000
production cycles of the manuscript input bundle. It is substantially slower.
Check equilibration and sampling before interpreting its results. The imposed
pressure controls NPT acceptance. The printed classical virial/pressure is not
an MBX mechanical-pressure estimator.

## Output

Native reports and coordinate snapshots appear under `output/`, including
`output_*.s0.txt`, `output_*.s0.json` and `restart_*.s0.json`. The short input
opts into `energy_terms.s0.csv` for accepted-move diagnostics. RDF and molecular
movie outputs are also enabled. Short-run averages and error estimates are
not converged. Use a fresh directory for each run to keep outputs separate.
