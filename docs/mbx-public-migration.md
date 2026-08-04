# Public MBX dependency audit

Audit date: 2026-08-04

## Decision

RASPA3-MBX can use the public
[`paesanilab/MBX`](https://github.com/paesanilab/MBX) repository instead of the
private `hassan-azizi/MBX-dev` fork. The release validated here is MBX `v1.4.0`
at commit `0e01b75b47611d7d51f27a34b112bdc5e2090a50`.

There is no private scientific model, coefficient set, data file, or RASPA
interface in the private fork. Public MBX is therefore the maintained
dependency for new builds. Keep the old private checkout read-only until a
public-MBX RASPA3-MBX release and its reference outputs have been archived, but
do not require new users to clone it.

The migration is source-compatible, not binary-compatible with already
compiled RASPA object files. A deliberate old-object/new-archive relink crashed;
recompiling the adapter against the public `system.h` and bundled JSON headers
fixed the failure. Never combine headers, objects, and an archive from different
MBX revisions.

Keep MBX as a pinned external dependency rather than copying its source or
generated library into RASPA3-MBX. This preserves one canonical RASPA3-MBX
repository, keeps MBX's separate license visible, and avoids adding hundreds of
megabytes to every RASPA clone. For repeated deployment, publish an internal
versioned MBX build artifact or container made from the recorded public commit;
do not commit that artifact to Git.

## Repository comparison

| Item | Private MBX checkout | Public MBX pin |
|---|---|---|
| Repository | `hassan-azizi/MBX-dev` | `paesanilab/MBX` |
| Version | 1.3.3 | 1.4.0 |
| Commit | `ba6727eb262c753df16a424a442c3893f64c8188` | `0e01b75b47611d7d51f27a34b112bdc5e2090a50` |
| Shared base | `1a0116cbaec40b2a230dafafc8cccdffca4f0575` | same |
| Commits after the shared base | one private commit | 262 public commits |
| License | UC research/nonprofit license | byte-identical license |

The one private commit changes only
`src/potential/electrostatics/helpme.h`. It replaces an unguarded call to
`omp_get_thread_num()` with an already guarded thread identifier. The assigned
local variable is unused, so the patch has no numerical effect; it only avoids
a compile error in a build where OpenMP is disabled. Public v1.4.0 rewrites
that interpolation loop and removes the problematic declaration entirely.
There is consequently no private patch to carry forward or submit upstream.

The following relevant trees were also compared directly:

- `examples/PEFs` is byte-identical between the private head and public
  v1.4.0, so the checked potential examples and data did not disappear;
- `src/bblock/system.cpp` is byte-identical;
- `src/bblock/system.h` retains the required API and changes only internal
  includes from `json/json.h` to `json/json.hpp`, plus the public constants
  header;
- public v1.4.0 contains later fixes and performance work in dispersion,
  electrostatics, JSON, plugins, and additional potentials. Those changes make
  numerical regression necessary even though the API is source-compatible.

## What RASPA3-MBX actually requires

The adapter includes only `bblock/system.h` and uses the stock
`bblock::System` C++ interface:

- construction and destruction;
- `GetNumMon`, `AddMonomer`, `SetUpFromJson`, `SetPBC`, and the four-argument
  `SetExternalChargesAndPositions` overload;
- `OneBodyEnergy`, `TwoBodyEnergy`, `ThreeBodyEnergy`, `FourBodyEnergy`,
  `Dispersion`, and `Electrostatics`;
- `GetPermanentElectrostaticEnergy` and
  `GetInducedElectrostaticEnergy`.

Every declaration and compiled symbol remains present in public v1.4.0. RASPA
does not depend on MBX's MPI/LAMMPS interface, gradient interface, C API,
Buckingham call, or Lennard-Jones call. Use MBX's basic non-MPI installation,
not its special MBX_MPI/LAMMPS build.

The checked RASPA inputs use stock lowercase MBX monomer identifiers `co2` and
`h2o`. Framework atoms are supplied as external point charges. The Morse
framework-guest interaction is implemented in RASPA3-MBX and is independent of
the MBX fork, so changing the MBX dependency does not remove Morse support.

The raw MBX term order remains:

```text
1B, 2B, 3B, 4B, dispersion, permanent electrostatics, induced electrostatics
```

RASPA3-MBX retains `2B + 3B + 4B + dispersion + permanent + induced` after
subtracting the bare-framework permanent electrostatic energy. The raw 1B term
is diagnostic and is not included in the retained total.

## Build and ABI rules

MBX exposes C++ types such as `std::string` and `std::vector` across the static
library boundary, and RASPA constructs `bblock::System` itself. Treat the MBX
headers and archive as one inseparable installation. They must use compatible
compiler ABI, C++ standard library, architecture flags, FFTW ABI, and OpenMP
runtime.

For every MBX update:

1. check out a tag and full commit hash, not a moving branch;
2. configure and build in a fresh source/worktree and install to a fresh,
   versioned prefix;
3. point `MBX_INCLUDE_DIR`, `MBX_CONFIG_INCLUDE_DIR`, and `MBX_LIBRARY` to that
   same prefix;
4. configure RASPA in a fresh build directory so stale CMake cache paths cannot
   mix versions;
5. inspect the final executable with `ldd` and require one OpenMP runtime.

The former private installation is not a reproducible build baseline. Its
installed archive was produced with Conda Clang 20, the source tree was later
reconfigured with system Clang 14, and its installed `single_point` executable
loads both GNU `libgomp` and LLVM `libomp`. Git reports the checkout as clean
because generated archives are ignored. Do not rebuild incrementally in that
directory.

The validated Clang workflow passes `-L$CONDA_PREFIX/lib` to MBX at configure
time. This allows MBX's `-liomp5` probe to resolve to the Conda symlink for
LLVM `libomp`; without the early library path, configure can fall through to
GNU `-lgomp` while Clang's `-fopenmp` still adds LLVM `libomp`. The complete
commands are in the main [README](../README.md#build-the-audited-public-mbx-dependency).

## Verification results

All checks below used public MBX v1.4.0 at the pinned commit, Conda Clang
20.1.6, `libstdc++`, Conda FFTW 3.3.10, and LLVM `libomp`.

| Check | Result |
|---|---|
| Public MBX configure/build/install | Passed from a fresh worktree and fresh prefix |
| Public MBX upstream suite | 34/34 tests passed |
| Installed MBX runtime libraries | Conda `libfftw3.so` and one `libomp.so`; no `libgomp` |
| Required RASPA C++ symbols | All declarations and archive symbols present |
| CO2 gas example, private vs public | Both `-6.5793057239` kcal/mol at printed precision |
| CO2 periodic example, private vs public | Both `-6.5825805658` kcal/mol at printed precision |
| MB-pol water gas example | Public/private difference about `2.0e-9` kcal/mol |
| MB-pol water periodic example | Public/private difference about `1.8e-9` kcal/mol |
| Public MB-pol periodic, OpenMP 1 vs 8 | Identical at printed precision |
| RASPA3-MBX configure and full compile | Passed against only the public install |
| Focused RASPA MBX/control suite, OpenMP 1 | 12/12 tests passed |
| RASPA static MBX suite, OpenMP 8 | 4/4 tests passed |
| Final RASPA runtime libraries | Conda `libfftw3.so`, `libomp.so`, and `libstdc++.so`; no second OpenMP runtime |
| MOF-74 + one CO2 + Morse Monte Carlo smoke test | Passed; Morse pairs and accepted-energy CSV were active |
| MOF-74 + CO2 fixed-snapshot public/private comparison | All seven raw MBX terms agree; largest difference was `8.731e-11` kcal/mol in permanent electrostatics |
| Same fixed-snapshot retained/total energy | Difference `4.394e-8` K, far below the existing `1e-6` regression tolerance |
| Public MOF-74 incremental vs recomputed MBX energy | Drift about `4.394e-8` K |

The fixed snapshot comparison also confirmed that the 1B, 2B, 3B, and 4B
values were identical for that case, while the tiny difference came from the
updated public electrostatics implementation. Reference values must not be
silently regenerated if a future MBX update produces a material difference.

## Corrections to the January 2026 build guide

| Old instruction | Maintained instruction |
|---|---|
| Clone a private `MBX-dev` fork | Clone public `paesanilab/MBX` and check out the pinned commit |
| Treat `paesanilab/MBX-dev` as the public origin | The public repository is `paesanilab/MBX` |
| Clone the old `RASPA3.0.21-MBX` name and `development` branch | Clone `Chung-Research-Group/RASPA3-MBX` and use its tested `main` branch |
| CMake 3.29 | Current RASPA3-MBX requires CMake 4.0.3 or newer |
| Install MBX into its source directory | Use `./configure --prefix=<versioned-prefix>` and keep source and install provenance separate |
| Allow configure to discover any OpenMP library | Use one compiler/runtime deliberately and verify the final `ldd` output |
| Carry an unpinned MBX checkout | Record the MBX tag, full commit, compiler, standard library, FFTW, and OpenMP runtime |
| Use legacy `NumberOfCycles` in retained RASPA inputs | Use current `NumberOfProductionCycles` |

## Upgrade policy

Do not copy a new MBX archive over the validated prefix. For a later MBX
release, create a dependency-update branch, use a new versioned prefix and new
RASPA build directory, repeat the public MBX suite and RASPA MBX/Morse gates,
and compare all raw terms from fixed snapshots. Merge only after recording the
new pin and reference comparison in this document or a release note.

MBX's license permits copying, modification, and distribution for educational,
research, and nonprofit purposes with its notices retained. Commercial use
requires permission from the University of California. It is separate from
RASPA3-MBX's MIT license.
