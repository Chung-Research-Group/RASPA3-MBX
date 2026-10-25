# Maintained verification record — 26 October 2026

The code in this edition is maintenance commit `281772cbd61c805f1a1b713b5dd290187bd0cff9`, descended from `df7a90656e8e2a8ad276ca8310b161ec61f506af`. Subsequent R155 changes are documentation and evidence only. The full test receipts below belong to that exact tested code, not to every future build. [Source binding](validation/r155/SOURCE_BINDING.json) records the six changed files and their hashes.

| Evidence | Result | Provenance and scope |
|---|---|---|
| Complete RASPAKit suite, MBX enabled | 522/522 executed pass; 6 disabled utilities remain | R153 native execution; 528 registered. |
| Focused MBX suite | 31/31 pass | R153 native execution on the same repaired source. |
| Affected general groups, MBX disabled | 59/59 pass | R153 native execution. |
| Independent affected-group repetition | 59/59 pass with MBX; 59/59 pass without MBX | R155 reviewer execution of hash-bound R153 binaries. Original df7 control reproduces all six failures in 56 affected tests. |
| Independent focused MBX repetition | 31/31 pass | R155 reviewer execution on df7; those MBX sources are unchanged in 281772c. |

[Original full receipts](validation/r155/r153_FINAL_VERIFICATION.json), [independent affected receipts](validation/r155/r155_independent_affected_receipts.json), and adjacent XML preserve commands, binary hashes, counts and failures. R153 full-suite execution used the working patch on df7 before committing it as 281772c; the receipt records both the base and final candidate, and the committed patch/source hashes bind them. R155 did not claim a new full rebuild or rerun of all 522 tests. Six pre-existing disabled utilities were neither newly disabled nor counted as passes.

## Repairs and effect on the paper

Five failures concerned test setup, finite-difference references or open-file cleanup. One was a shared exact-geometry defect: duplicate or concentric spherical caps were not pruned correctly. The repair changes only `src/structurekit/diagrams/exact/sphere_sweep.ixx` in production code. All MBX adapter, MC proposal/acceptance, energy, molecular-model and adsorption-estimator files are byte-identical to df7. No tolerance was relaxed, expected result altered to match an implementation, or failing test disabled. The [failure-by-failure report](validation/r155/r153_REPORT.md) explains each fix and its limits.

These bounded repairs do not require new trajectories for the paper. The four repository examples and the archived snapshot-compatibility checks were executed on df7 and remain historical evidence on unchanged application paths; they are not presented as newly executed R155 examples. Test success is not proof of ensemble convergence or exhaustive correctness. No experimental TMMC/flexible-guest branch is merged, no shared module is deployed, and no remote publication is implied by this local review branch.

## Reproduction

Use the build and CTest instructions in the root README. The [R153 reproduction record](validation/r155/r153_REPRODUCE.md) preserves original local commands; paths and installed dependencies must be adapted on another machine. The [earlier synchronization report](mbx-upstream-sync.md#verification-performed) is historical and reports a different test inventory. Its counts must not be substituted for this edition's receipts.
