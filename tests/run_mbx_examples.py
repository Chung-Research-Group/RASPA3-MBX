#!/usr/bin/env python3
"""Run each complete MBX example alone and verify native physical/output contracts.

These are bounded functional checks, not ensemble-convergence or accuracy benchmarks.
Each example is copied independently, without parent definitions or other example trees.
Source files and the copied simulation inputs are never edited. Requires Python 3.10+.
"""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

EXAMPLES = (
    "3_snapshot_energy_evaluation",
    "1_energy_logging_and_restarts",
    "2_widom_zero_loading_qst",
    "4_gcmc_co2_adsorption",
    "5_water_npt",
    "6_gcmc_co2_mg_mof74",
)
ENERGY_TOLERANCE_K = 1.0e-5


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def local_file(directory, value, suffix=""):
    """Require a real input within this example, including through symlink resolution."""
    relative = Path(value)
    require(not relative.is_absolute(), f"Absolute input path is not portable: {value}")
    if suffix and relative.suffix != suffix:
        relative = Path(str(relative) + suffix)
    candidate = (directory / relative).resolve()
    require(
        candidate.is_relative_to(directory.resolve()),
        f"Input escapes example directory: {value}",
    )
    require(candidate.is_file(), f"Missing local input: {relative}")
    return candidate


def verify_inputs(directory):
    directory = directory.resolve()
    require((directory / "run.sh").is_file(), "Missing local run.sh")
    dependencies = set()
    for input_path in sorted(directory.glob("simulation*.json")):
        data = json.loads(input_path.read_text())
        ff = local_file(
            directory, str(Path(data.get("ForceField", ".")) / "force_field.json")
        )
        dependencies.add(ff)
        force_field = json.loads(ff.read_text())
        elements = {
            atom["name"]: atom["element"] for atom in force_field["PseudoAtoms"]
        }
        for component in data["Components"]:
            name = component["Name"]
            molecule_path = local_file(directory, name, ".json")
            dependencies.add(molecule_path)
            model = json.loads(molecule_path.read_text())
            require(
                model["Type"].lower() == "rigid",
                "Main-branch MBX examples require rigid guests",
            )
            atoms = model["pseudoAtoms"]
            order = [elements[atom[0]] for atom in atoms]
            expected = {"co2": ["C", "O", "O"], "h2o": ["O", "H", "H"]}
            require(
                name in expected and order == expected[name],
                f"Incorrect MBX {name} atom order: {order}",
            )
            distances = [math.dist(atoms[0][1], atom[1]) for atom in atoms[1:]]
            require(
                all(0.8 < distance < 1.3 for distance in distances),
                f"Distorted {name} template: check central-atom order and bond lengths",
            )
            if name == "co2":
                require(
                    math.dist(atoms[1][1], atoms[2][1]) > 2.0,
                    "The CO2 oxygen endpoints must be distinct",
                )
        for system in data["Systems"]:
            require(system.get("UseMBX") is True, "Example must explicitly enable MBX")
            settings_path = local_file(directory, system["MBXSettingsFile"])
            dependencies.add(settings_path)
            require(
                "MBX" in json.loads(settings_path.read_text()), "Malformed MBX settings"
            )
            if system["Type"].lower() == "framework":
                dependencies.add(local_file(directory, system["Name"], ".cif"))
            if "RestartFileName" in system:
                snapshot_path = local_file(directory, system["RestartFileName"])
                dependencies.add(snapshot_path)
                positions = json.loads(snapshot_path.read_text())["co2"]
                require(
                    len(positions) == 3, "Snapshot example must contain exactly one CO2"
                )
                require(
                    all(
                        0.8 < math.dist(positions[0], point) < 1.3
                        for point in positions[1:]
                    ),
                    "Snapshot carbon must be central and first, as in the component definition",
                )
                require(
                    math.dist(positions[1], positions[2]) > 2.0,
                    "Snapshot oxygen endpoints must be distinct",
                )
        if data["SimulationType"] == "MonteCarlo":
            require(
                data.get("NumberOfProductionCycles", 0)
                >= data.get("NumberOfBlocks", 5),
                "Example must provide at least one production cycle per statistical block",
            )
    require(dependencies, "No simulation input dependencies found")
    return sorted(str(path.relative_to(directory)) for path in dependencies)


def retained_energy(energy):
    """These rigid examples have no external fields or intramolecular target terms."""
    require(
        all(
            isinstance(value, (int, float)) and math.isfinite(value)
            for value in energy.values()
        ),
        "Nonfinite native energy",
    )
    total = energy["Total potential energy [K]"]
    expected = sum(
        energy[key]
        for key in (
            "MBX energy [K]",
            "framework-molecule VDW [K]",
            "Van der Waals (Tail) [K]",
        )
    )
    require(
        math.isclose(total, expected, rel_tol=1.0e-10, abs_tol=ENERGY_TOLERANCE_K),
        "Retained total must equal MBX plus framework-guest vdW and tail",
    )
    return total


def verify_logged_energy(rows):
    for row in rows:
        values = [float(row[key]) for key in ("total", "hg_VDW", "hg_tail", "mbx_tot")]
        require(
            all(math.isfinite(value) for value in values),
            "Nonfinite accepted-state energy",
        )
        require(
            math.isclose(
                values[0], sum(values[1:]), rel_tol=1.0e-10, abs_tol=ENERGY_TOLERANCE_K
            ),
            "Accepted-state target energy is inconsistent",
        )
        terms = [
            float(row[key])
            for key in ("E2b", "E3b", "E4b", "Edisp", "Eelec_perm", "Eelec_ind")
        ]
        require(
            all(math.isfinite(value) for value in terms), "Nonfinite MBX decomposition"
        )
        require(
            math.isclose(
                values[-1], sum(terms), rel_tol=1.0e-10, abs_tol=ENERGY_TOLERANCE_K
            ),
            "MBX decomposition does not sum to its retained energy",
        )


def summed(value):
    return sum(value) if isinstance(value, list) else value


def move_counts(report, component, move):
    record = report["output"]["MCMoveStatistics"]["system"][component][move]
    return {
        "attempted": summed(record["total"]),
        "accepted": summed(record["accepted"]),
    }


def verify_output(name, directory):
    output = directory / "output"
    if name == "3_snapshot_energy_evaluation":
        data = json.loads((output / "energy_evaluation.json").read_text())
        require(
            data["simulationType"] == "EnergyEvaluation", "Incorrect execution mode"
        )
        system = data["systems"][0]
        require(system["model"] == "MBX", "Snapshot did not use MBX")
        require(
            system["moleculeCounts"]["total"] == 1,
            "Snapshot changed the molecule count",
        )
        total = retained_energy(system["energy"])
        raw = system["mbxRawEnergyTerms [kcal/mol]"]
        require(
            all(math.isfinite(value) for value in raw.values()),
            "Nonfinite raw MBX decomposition",
        )
        # 1B is omitted from the rigid target; it still diagnoses a distorted/misordered monomer.
        require(
            abs(raw["oneBody (1B)"]) < 1.0,
            "Unexpected CO2 one-body energy: check atom order",
        )
        return {"raw_one_body_kcal_mol": raw["oneBody (1B)"], "total_energy_K": total}

    text_paths, json_paths = list(output.glob("output_*.s0.txt")), list(
        output.glob("output_*.s0.json")
    )
    require(
        len(text_paths) == len(json_paths) == 1,
        "Expected one native text and JSON report",
    )
    text, report = text_paths[0].read_text(), json.loads(json_paths[0].read_text())
    data = json.loads((directory / "simulation.json").read_text())
    cycles = data["NumberOfProductionCycles"]
    require(
        "Simulation finished!" in text
        and re.search(rf"Final state after {cycles} cycles\b", text),
        "Requested production cycles did not complete",
    )
    require(
        report["initialization"]["initialConditions"]["useMBX"] is True,
        "Native run did not enable MBX",
    )
    running = retained_energy(report["output"]["runningEnergies"])
    recomputed = retained_energy(report["output"]["recomputedEnergies"])
    require(
        math.isclose(running, recomputed, rel_tol=1.0e-10, abs_tol=ENERGY_TOLERANCE_K),
        "Retained final energy disagrees with independent native recomputation",
    )
    # Undefined short-run fluctuation/pressure statistics are deliberately not validated here.
    rows = load_rows(output / "energy_terms.s0.csv")
    verify_logged_energy(rows)
    result = {
        "production_cycles": cycles,
        "accepted_rows": len(rows),
        "final_energy_K": running,
        "recomputed_energy_K": recomputed,
        "energy_difference_K": running - recomputed,
    }
    if name == "2_widom_zero_loading_qst":
        require(not rows, "Widom ghost insertions must not be logged as accepted moves")
        widom = load_rows(output / "widom_energy_terms.s0.csv")
        require(widom, "No completed Widom trial was logged")
        require(
            all(int(row["N"]) == 0 and int(row["trial_N"]) == 1 for row in widom),
            "Widom must preserve the empty physical system",
        )
        require(
            "Qst^0" in text, "The report must include the zero-loading heat estimator"
        )
        result["widom_rows"] = len(widom)
        return result

    require(rows, "No accepted move was logged")
    restarts = sorted(output.glob("restart_*.s0.json"))
    require(restarts, "Missing configured coordinate restart")
    if name == "5_water_npt":
        require(
            all(int(row["N"]) == 16 for row in rows),
            "NPT changed its fixed 16-water population",
        )
        for restart in restarts:
            snapshot = json.loads(restart.read_text())
            require(
                len(snapshot["h2o"]) == 48,
                "NPT coordinate restart must contain exactly 16 three-site waters",
            )
            box = snapshot["SimulationBox"]
            lengths = [box[key] for key in ("length-a", "length-b", "length-c")]
            require(
                all(math.isfinite(value) and value > 0 for value in lengths),
                "NPT coordinate restart must preserve valid cell lengths",
            )
        # Volume is a system move; the component's placeholder record stays zero.
        volume = move_counts(report, "system", "Volume change")
        require(volume["attempted"] > 0, "No NPT volume move was attempted")
        result.update(
            population=16,
            volume=volume,
            restart_box_lengths_A=lengths,
            accepted_volume_rows=sum(row["type"] == "volume" for row in rows),
        )
    elif name == "6_gcmc_co2_mg_mof74":
        swap = move_counts(report, "co2", "Swap (CBMC)")
        reinsertion = move_counts(report, "co2", "Reinsertion (CBMC)")
        require(swap["attempted"] > 0, "No CBMC insertion/deletion was attempted")
        require(reinsertion["attempted"] > 0, "No CBMC reinsertion was attempted")
        result.update(
            cbmc_swap=swap,
            cbmc_reinsertion=reinsertion,
            accepted_populations=sorted({int(row["N"]) for row in rows}),
        )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--timeout",
        type=float,
        default=600.0,
        help="Maximum seconds per unchanged example",
    )
    parser.add_argument(
        "--examples",
        nargs="+",
        choices=EXAMPLES,
        default=list(EXAMPLES),
        help="Run only these named examples (default: all six)",
    )
    args = parser.parse_args()
    source, executable, output = (
        args.source.resolve(),
        args.executable.resolve(),
        args.output.resolve(),
    )
    require(
        executable.is_file() and os.access(executable, os.X_OK),
        f"Executable not found: {executable}",
    )
    require(
        math.isfinite(args.timeout) and args.timeout > 0,
        "Timeout must be finite and positive",
    )
    require(
        len(set(args.examples)) == len(args.examples),
        "Repeated example names are not allowed",
    )
    output.mkdir(parents=True, exist_ok=False)
    examples = output / "examples"
    evidence = {
        "executable": str(executable),
        "sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "source": str(source),
        "input_edits": [],
        "examples": [],
        "input_sha256": {},
        "energy_tolerance_K": ENERGY_TOLERANCE_K,
        "isolation": "Each example copied alone; no shared parent inputs; RASPA_DIR is its own directory",
    }
    environment = os.environ.copy()
    environment.update(
        OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", RASPA3_EXECUTABLE=str(executable)
    )
    result_path = output / "results.json"
    for name in args.examples:
        started = time.monotonic()
        directory = examples / name
        row = {"name": name, "command": ["./run.sh"], "cwd": str(directory)}
        evidence["examples"].append(row)
        try:
            shutil.copytree(
                source / "examples/raspa3_mbx" / name,
                directory,
                ignore=shutil.ignore_patterns(
                    "output", "bias_factors", "restart_data.bin*", "__pycache__"
                ),
            )
            row["local_input_dependencies"] = verify_inputs(directory)
            hashes = {
                str(path.relative_to(directory)): hashlib.sha256(
                    path.read_bytes()
                ).hexdigest()
                for path in sorted(directory.rglob("*"))
                if path.is_file()
            }
            evidence["input_sha256"][name] = hashes
            environment["RASPA_DIR"] = str(directory)
            with (output / f"{name}.log").open("w") as log:
                completed = subprocess.run(
                    ["./run.sh"],
                    cwd=directory,
                    env=environment,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=args.timeout,
                    check=False,
                )
            row["returncode"] = completed.returncode
            require(completed.returncode == 0, f"{name} exited {completed.returncode}")
            require(
                all(
                    hashlib.sha256((directory / relative).read_bytes()).hexdigest()
                    == digest
                    for relative, digest in hashes.items()
                ),
                "Execution changed supplied example files",
            )
            row["checks"] = verify_output(name, directory)
            row["status"] = "PASS"
        except (
            ValueError,
            KeyError,
            TypeError,
            OSError,
            subprocess.TimeoutExpired,
        ) as error:
            row.update(status="FAIL", error=str(error))
            print(f"{name}: {error}", file=sys.stderr)
        row["seconds"] = time.monotonic() - started
        result_path.write_text(json.dumps(evidence, indent=2) + "\n")
    passed = sum(row["status"] == "PASS" for row in evidence["examples"])
    evidence["status"] = "PASS" if passed == len(args.examples) else "FAIL"
    result_path.write_text(json.dumps(evidence, indent=2) + "\n")
    print(
        f"{passed}/{len(args.examples)} functional examples passed; evidence: {result_path}"
    )
    return 0 if evidence["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
