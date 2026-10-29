# Maintaining the MBX interface

## Scope and branch ownership

`main` contains the current paper's rigid-guest MBX interface and its maintenance
repairs. `development` identifies the latest preserved MBX-TMMC recovery source;
`feature/mbx-flexible-exchange` preserves the separate flexible-target prototype.
The old remote development tip is retained as `archive/development-pre-20261022`.
These branch roles do not deploy a new executable or alter running calculations.

The two experimental implementations are not interchangeable. The flexible
snapshot excludes TMMC rigid reinsertion that the recovery implementation enables.
Porting one onto the other requires a deliberate conflict review and new tests;
copying `system.cpp` wholesale would undo a scientifically relevant guard change.
Source-preservation commits identify their original directory and date. They do
not pretend to reconstruct the chronological commits of earlier development.

The attribution baseline is original RASPA3 `a8a70ca755b2821728bee05c63ba07d8b10085b5`.
The public pre-maintenance tip is `c46c859cec8446cf4dc3d6f68cec839f56f1df33`.
Future synchronization with a newer upstream revision belongs on its own branch.
Upstream acceptance and the form of an eventual contribution remain the original
maintainers' decisions; a clean local branch is not an accepted upstream branch.

## Retained target and state contract

For rigid guests the retained energy is external-field energy, framework–guest
van der Waals energy and its tail, plus MBX 2B/3B/4B, dispersion, permanent
and induced electrostatics. Bare-framework permanent electrostatics is subtracted
once. MBX 1B is reported but omitted. Classical guest energies guide proposals and
screening and must not be added to the retained target. `U_trial - U_current`
uses the complete retained state energy. Insertion/deletion/CBMC reservoir and
Rosenbluth factors remain in the move-specific acceptance code.

Coordinates are in angstroms. The seven diagnostic terms are in kcal/mol;
`RunningEnergy` uses RASPA internal units. Site order is a model contract:
`co2`: C,O,O; `h2o`: O,H,H. Pseudo-atom names do not reorder positions. Other
monomers require their installed database's order and model coverage. Do not
sort coordinates or introduce generic body terms for unsupported combinations.
The main target rejects flexible guests instead of silently omitting their
conformation-dependent 1B energy.

The reported virial/pressure estimator remains classical. Its output is marked
as a diagnostic, not an MBX mechanical pressure. The imposed pressure used in
NPT acceptance is a separate input. A true MBX virial requires a separately
implemented force/strain interface; relabeling does not compute it.

## Adapter entrypoints

All entrypoints create a fresh MBX evaluator. There is currently no persistent
MBX energy cache. The molecular storage is component-major, with global atom
molecule IDs converted to the selected component's local index.

| Operation | Entry and selected-molecule meaning | Returned quantity |
|---|---|---|
| Entire retained configuration | `computeMBXEnergySystem`, complete molecular span and matching counts | Absolute retained MBX contribution |
| Replace an existing molecule | `computeMBXEnergy`, existing global ID, replacement span, `includeSelectedMoleculeAtoms=true` | Absolute contribution with that molecule replaced |
| Insert a new molecule / Widom ghost | `computeMBXEnergy`, a new ID that identifies no existing molecule, trial span, `true` | Absolute contribution with trial molecule included |
| Delete an existing molecule | `computeMBXEnergy`, valid existing ID/span, `false` | Absolute contribution with that molecule excluded |
| Replacement difference only | `computeMBXEnergyDifference`, nonempty old/new spans identifying the same existing molecule | New minus old retained MBX contribution |

The difference entrypoint must not be used for exchange or Widom. Absolute
energy is not a move acceptance probability. The diagnostic order is
`[1B, 2B, 3B, 4B, dispersion, permanent-minus-bare-host, induced]`.
An optional output span must have exactly seven elements. The existing
vector-output overload is retained for compatibility; a typed decomposition or
named operation wrappers should be a separate, tested API change.

## Settings, restart and dependency identity

Keep the MBX settings file immutable throughout a run and across binary restart.
The adapter rereads the file; the checkpoint stores its path and the bare-host
reference, not a content fingerprint. Replacing it in place can change the target
while retaining an inconsistent host reference. A frozen-settings or fingerprinted
archive design is deferred until copy, restart compatibility and intentional
settings changes have dedicated tests. This limitation is not evidence that an
archived run changed settings.

Record the RASPA commit, compiler and flags, MBX source commit, header/library
hashes, input files and settings hashes. Use the pinned public MBX v1.4.0 commit
`0e01b75b47611d7d51f27a34b112bdc5e2090a50` for the maintained build. Do not infer
historical executable identities from present headers. Never relink old adapter
objects against a different MBX C++ ABI. Preserve archived release executables
and physical results when maintaining source.

## Verification and future maintenance

Build both `BUILD_MBX=ON` and `OFF` in fresh directories. The OFF build must not
include or link MBX. Run the MBX suite through CTest and the copied-example runner;
retain XML, commands and dependency identities. The full RASPAKit suite is
reported separately, including any failures and their available baseline controls. Passing focused tests does not
prove ensemble convergence or arbitrary screening support.

Keep each numbered `examples/raspa3_mbx` directory self-contained. Required
component definitions, MBX settings, force fields, framework structures and
snapshot inputs must be local, without parent-directory links or implicit
installed-database lookup. `tests/run_mbx_examples.py` copies and runs each
directory independently and checks completed production, relevant move
statistics and final energy recomputation. Water NPT and Mg-MOF-74 GCMC also
include longer input profiles; short functional defaults must remain explicitly
distinguished from equilibrated scientific calculations. Changes to example
parameters require a physical rationale, not just a successful process exit.

New moves should declare their proposal/target difference, accepted/rejected
state updates, timing/diagnostic semantics, supported ensembles and restart
behavior. Test an accepted state against independent full recomputation, a
rejected state plus the next move, empty and nonfirst-component identity cases,
and errors before mutation. Keep unsupported combinations rejected centrally.

Follow upstream `.clang-format` and `STYLE_GUIDE.md` for changed code, informed
by adjacent maintained modules. Add comments for units, signs, excluded energies
and state ownership; avoid narrating obvious syntax. Keep fixes and optional
extensions in separate reviewable commits. Do not mass-format inherited code,
weaken a regression tolerance, or silently regenerate golden energies to obtain
a passing result. Cloud CI, deployment and upstream review have distinct receipts.

Preserved experimental snapshots do not automatically include later main maintenance
fixes. Compare and port applicable fixes explicitly before integrating a feature.
