# RASPA3–MBX: six general-suite failures

23 October 2026. R153. This report addresses the author's request to explain the six failures, attempt repairs, and assess whether they must be fixed. The isolated candidate is based on published main `df7a90656e8e2a8ad276ca8310b161ec61f506af`; original evidence and sealed manuscript v42 are preserved.

## What the failures mean

The six failures do not mean six defective MBX energies or six invalid manuscript datasets. Five originate in test references, branch setup or file cleanup. One is a real boundary-case defect in a shared exact-geometry utility. All six were reproduced on the original executable in a fresh NFS directory before final validation. No original test tolerance was relaxed, no expected value was changed to match the implementation, and no failing test was disabled.

The original affected files are unchanged from the merged upstream RASPA3 source attribution base `a8a70ca755b2821728bee05c63ba07d8b10085b5`. This establishes inherited source, not an execution of a separate upstream binary. We nevertheless repaired the bounded problems locally because reproducible tests and a correct geometry utility are valuable for maintenance.

## Failure-by-failure findings

### 1. Flexible-framework HMC rejection fixture

Test: `hybrid_mc.flexible_framework_only_does_not_throw_and_restores_on_reject`; source: `tests/raspakit-tests/integrators.cpp:691`.

The fixture assumed that five steps with a large timestep would almost certainly be rejected. Independent harmonic-dimer propagation reproduces the native final coordinates to about 4e-15 Å. Its positive Hamiltonian change is 57.088781 internal energy units, giving acceptance probability 0.795429. The seeded random draw is 0.287749, so acceptance is the correct outcome for that proposal. The test's expectation of rejection was wrong.

Repair: retain the seed, velocities, physical parameters and integration path; condition only the final uniform acceptance draw to exercise both acceptance and rejection. Keep the original position tolerance and additionally verify restoration of velocities, gradients and timestep. Count the random draws to guard the controlled-fixture assumption. Production RNG and acceptance code are unchanged.

Necessity: needed for meaningful rejection-state regression coverage and a trustworthy green suite. It does not require recalculating this paper's rigid MC data.

### 2–3. Charged and polarizable variable-cell Hessian checks

Tests: `minimization_variable_cell.rigid_charged_mixed_blocks_match_finite_difference` and `minimization_variable_cell.polarization_rigid_molecule_real_space_matches_finite_difference`; helper: `tests/raspakit-tests/minimization_variable_cell.cpp:363`.

Both fail at the same null rotation/cell-shear entry: analytic curvature is zero, while the four-energy finite difference gives about 0.009037 or 0.009039 against tolerance 0.005. Rotation about the linear molecule's own axis should leave the configuration unchanged. The original perturbation 2e-5 yields a 1.78e-15 Å quaternion/regeneration difference, and the second difference amplifies the resulting tiny energy residual. Differentiating the analytic gradient gives zero; a perturbation scan confirms the unstable reference rather than a missing Hessian term.

Repair: use a perturbation of machine-epsilon to the one-quarter power, approximately 1.22e-4, balancing second-difference cancellation and truncation. All original rows, assertions and tolerances remain. The complete 16-test variable-cell group passes against the unchanged production library. This is a fixture-scale numerical choice, not a claim that one dimensional step is optimal for every system.

Necessity: needed for reliable numerical validation; neither failure establishes an error in the paper's energy or MC calculations. No scientific data rerun is required.

### 4. Normal-mode movie cleanup on NFS

Test: `normal_modes.writes_one_pdb_movie_per_mode`; source: `tests/raspakit-tests/normal_modes.cpp:162`.

The test keeps an input stream open while removing its parent directory. NFS preserves an unlinked open file as a temporary `.nfs*` entry, so directory removal throws “Directory not empty.” The same original executable passes on local XFS. The failure is after output-content checks, not a failed eigensolver or mode calculation.

Repair: close the reader after loading the contents and before removing the directory. All file-count, content and frame-count checks remain; the writer and numerical calculation are unchanged.

Necessity: needed for portable tests on this server and other NFS environments. It is not a reason to rerun manuscript simulations.

### 5. Taylor-shifted Lennard–Jones derivative reference

Test: `vdw_potentials.second_order_taylor_shifted_spatial_derivatives_match_finite_difference`; source: `tests/raspakit-tests/vdw_potentials.cpp:327`.

At r=5 Å, the analytic curvature is -41.9691567566572; the original three-point difference is -41.9685619590382. The error 5.95e-4 slightly exceeds the original 4.20e-4 tolerance. A 35-point native radius/step scan and an independently differentiated 80-decimal-digit LJ expression support the analytic result. The original step is too small for reliable subtraction of nearby energies.

Repair: use the same distance-scaled, fourth-order five-point difference already used by the general spatial-derivative test. At r=5 and the selected step, error is about 3.2e-7. Retain all five radii and original tolerances. Production energies and derivatives are unchanged.

Necessity: needed for stable regression testing, not for correcting the paper's force field or rerunning its calculations.

### 6. Duplicate spherical-cap pruning

Test: `exact_sphere_sweep.pruning_keeps_one_of_a_pair_of_equals_and_the_order_of_the_rest`; production source: `src/structurekit/diagrams/exact/sphere_sweep.ixx:77`.

For two identical caps, the old containment threshold evaluates to 1.0000000000000002 while the identical-axis dot product is 1. The duplicate survives, producing three caps instead of the expected two. This is a genuine floating-point boundary-case error in exact geometric analysis, not merely an incorrect test.

Repair: handle identical axes by ordered cosine half-angles; retain distinct equal-radius caps even when their axis dot product rounds to one; handle whole-sphere caps explicitly. The existing survivor order is preserved. No blanket epsilon is added. Three new regressions cover normalized axes, one-ULP radius differences, very close but distinct axes, and whole-sphere cases. Extracted-source probes have 38 cases: original 8 failures, repaired 0, under both release fast-math and strict floating-point controls.

Necessity: this should be fixed before relying on the affected exact-geometry functionality or offering it as a maintained upstream contribution. Direct callers are exact surface-patch, boundary-component and solvent-excluded geometry routines. They are not the MC/MBX energy path. The paper's excess-uptake conversion uses cell determinant, framework mass, recorded helium fraction and an external pore-volume sensitivity; it does not call this utility. Its current data therefore need no replacement.

## Scope and remaining limits

The complete patch changes one production geometry header and five test files. No MBX model, energy expression, force-field parameter, MC proposal, acceptance formula, accumulator, installed module, simulation input, scheduler state or sealed manuscript is changed. Two agents independently reviewed patches they did not author, including the physical-data implications. Their reports and the author dispositions are retained alongside this report.

A green local regression suite is not proof of every available algorithm or platform. Six pre-existing disabled entries remain disabled; they are distinct from the six executed failures investigated here. They are optional geometry timing/export/dataset-comparison utilities and a golden-record generator. Some require the upstream developer's external CIF directory; the record generator explicitly says it is not a test. Enabling these is not necessary to resolve the six failures, and regenerating its goldens would be inappropriate. Their exact names are retained in the final receipt. Hosted CI and a separate whole upstream baseline were not run. A pre-existing flexible-framework HMC expression using the absolute Hamiltonian difference remains outside these repairs and must not be certified by the positive-drift fixture. It is already a separate development concern and is not used by this paper's rigid MC calculations.

For this paper, these six failures are not a reason to discard results or initiate long replacement simulations. For software maintenance, fixing the five unreliable tests is worthwhile and the geometry correction is necessary for the affected functionality. Publication, upstream merging and deployment remain separate from this locally tested repair.

## Executed verification

Final exact runtime results and source identity are recorded in `FINAL_VERIFICATION.json` and `MAINTENANCE_RECORD.md`. Original-failure XML, native candidate XML, standalone diagnostics, independent reviews, source patches and configuration/build records are included in the evidence package. See those receipts for actual executed counts rather than CTest's source-scanned inventory.
