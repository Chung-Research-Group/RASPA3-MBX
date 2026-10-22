# MBX integration regression scope

The stored energy constants in `mbx.cpp` predate this maintenance review and
remain unchanged. The static MFI fixture compares those constants with full-system
and selected-molecule **adapter entrypoints**; it is not an independently executed
MBX driver. No original generation script has been recovered for those constants.
A changed dependency or target must be investigated, not accommodated by silently
regenerating the values. Public MBX's own model tests are a separate check.

`component_validation.cpp` checks site-order and rigid-model boundaries and
pressure-output semantics. `reporting.cpp` checks authoritative MBX statistics
without modifying serialized accumulators. `volume_integration.cpp` checks
accepted/rejected and screened volume proposals, cutoff restoration, next-move
behavior and unchanged RNG when diagnostics are toggled. Its nonfinite injection
is deliberately observed through object bytes under the supported fast-math build.

The command `ctest --test-dir build-mbx -R '^(mbx_|energy_evaluation|run_control|widom_heat)' --output-on-failure`
runs the integration cases registered in a build with `BUILD_TESTING_MBX=ON`.
Run `python tests/run_mbx_examples.py --help` for copied functional examples.
Short fixtures do not establish equilibrium convergence, general CBMC
forward/reverse screening support, or a validated MBX virial. Deterministic
exchange/reinsertion reversibility is a future extension of the native matrix;
existing algebra and archived-event checks must retain their separate provenance.
