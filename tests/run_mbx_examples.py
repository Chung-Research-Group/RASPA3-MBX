#!/usr/bin/env python3
"""Run installed MBX examples in a fresh directory and verify physical/output contracts.

These are functional checks, not ensemble-convergence or accuracy benchmarks.
Source examples and simulation inputs are never modified. Requires Python 3.10+.
"""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


EXAMPLES = (
    "3_snapshot_energy_evaluation",
    "1_energy_logging_and_restarts",
    "2_widom_zero_loading_qst",
    "4_gcmc_co2_adsorption",
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def verify_inputs(examples):
    model = json.loads((examples / "raspa3_mbx/co2.json").read_text())
    atoms = model["pseudoAtoms"]
    require([item[0] for item in atoms] == ["C_co2", "O_co2", "O_co2"],
            "MBX co2 requires C,O,O atom order; executable success cannot validate this contract")
    positions = json.loads((examples / "raspa3_mbx/3_snapshot_energy_evaluation/co2_snapshot.json").read_text())["co2"]
    require(len(positions) == 3, "The fixed snapshot must contain exactly one CO2")
    require(math.dist(positions[0], positions[1]) < 1.3 and math.dist(positions[0], positions[2]) < 1.3,
            "Snapshot carbon must be the central atom, in the same order as the component definition")
    require(math.dist(positions[1], positions[2]) > 2.0, "The CO2 oxygen endpoints must be distinct")


def verify_output(name, directory):
    output = directory / "output"
    if name == "3_snapshot_energy_evaluation":
        data = json.loads((output / "energy_evaluation.json").read_text())
        require(data["simulationType"] == "EnergyEvaluation", "Incorrect execution mode")
        system = data["systems"][0]
        require(system["moleculeCounts"]["total"] == 1, "Snapshot changed the molecule count")
        energy = system["energy"]
        require(all(math.isfinite(value) for value in energy.values()), "Nonfinite snapshot energy")
        raw = system["mbxRawEnergyTerms [kcal/mol]"]
        # The supplied 1.149-Angstrom rigid monomer has a small, nonzero 1B energy.
        # A grossly distorted/misordered molecule must fail even if 1B is omitted from the target.
        require(abs(raw["oneBody (1B)"]) < 1.0, "Unexpected CO2 one-body energy: check monomer atom order")
        return {"raw_one_body_kcal_mol": raw["oneBody (1B)"],
                "total_energy_K": energy["Total potential energy [K]"]}
    require(list(output.glob("output_*.s0.txt")), "Missing native text report")
    require(list(output.glob("output_*.s0.json")), "Missing native JSON report")
    rows = load_rows(output / "energy_terms.s0.csv")
    if name == "2_widom_zero_loading_qst":
        require(not rows, "Widom ghost insertions must not be logged as accepted moves")
        widom = load_rows(output / "widom_energy_terms.s0.csv")
        require(widom, "No completed Widom trial was logged")
        require(all(int(row["N"]) == 0 and int(row["trial_N"]) == 1 for row in widom),
                "Widom must preserve the empty physical system")
        report = "\n".join(path.read_text() for path in output.glob("output_*.s0.txt"))
        require("Qst^0" in report, "The report must include the zero-loading heat estimator")
        return {"widom_rows": len(widom), "accepted_rows": 0}
    require(rows, "No accepted move was logged")
    require(list(output.glob("restart_*.s0.json")), "Missing configured coordinate restart")
    return {"accepted_rows": len(rows)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=600.0, help="Maximum seconds per unchanged example")
    args = parser.parse_args()
    source, executable, output = args.source.resolve(), args.executable.resolve(), args.output.resolve()
    require(executable.is_file(), f"Executable not found: {executable}")
    output.mkdir(parents=True, exist_ok=False)
    examples = output / "examples"
    for relative in ("raspa3_mbx", "basic/10_mc_adorption_co2_in_mfi"):
        shutil.copytree(source / "examples" / relative, examples / relative,
                        ignore=shutil.ignore_patterns("output", "bias_factors", "restart_data.bin*", "__pycache__"))
    evidence = {"executable": str(executable), "sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
                "source": str(source), "input_edits": [], "examples": [],
                "input_sha256": {
                    str(path.relative_to(examples)): hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in sorted(examples.rglob("*")) if path.is_file()
                }}
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", RASPA3_EXECUTABLE=str(executable))
    result_path = output / "results.json"
    try:
        verify_inputs(examples)
    except (ValueError, KeyError, OSError) as error:
        evidence.update(status="FAIL", error=str(error))
        result_path.write_text(json.dumps(evidence, indent=2) + "\n")
        print(error, file=sys.stderr)
        return 1
    for name in EXAMPLES:
        started = time.monotonic()
        command = [str(examples / "raspa3_mbx/run-example.sh"), name]
        row = {"name": name, "command": command}
        evidence["examples"].append(row)
        try:
            with (output / f"{name}.log").open("w") as log:
                result = subprocess.run(command, cwd=examples / "raspa3_mbx", env=environment,
                                        stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout, check=False)
            row["returncode"] = result.returncode
            require(result.returncode == 0, f"{name} exited {result.returncode}")
            row["checks"] = verify_output(name, examples / "raspa3_mbx" / name)
            row["status"] = "PASS"
        except (ValueError, KeyError, OSError, subprocess.TimeoutExpired) as error:
            row.update(status="FAIL", error=str(error))
            print(f"{name}: {error}", file=sys.stderr)
        row["seconds"] = time.monotonic() - started
        result_path.write_text(json.dumps(evidence, indent=2) + "\n")
    passed = sum(row["status"] == "PASS" for row in evidence["examples"])
    evidence["status"] = "PASS" if passed == len(EXAMPLES) else "FAIL"
    result_path.write_text(json.dumps(evidence, indent=2) + "\n")
    print(f"{passed}/{len(EXAMPLES)} functional examples passed; evidence: {result_path}")
    return 0 if evidence["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
