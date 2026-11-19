#!/usr/bin/env python3
"""Run native single- and two-component TMMC tests in new case directories."""
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def read_json(path):
    return json.loads(path.read_text())


def ideal_gas(work):
    atom = {"name": "X", "framework": False, "print_to_output": True,
            "element": "Ar", "print_as": "Ar", "mass": 39.948, "charge": 0.0}
    interaction = {"name": "X", "type": "lennard-jones", "parameters": [0.0, 1.0]}
    ff = {"MixingRule": "Lorentz-Berthelot", "TruncationMethod": "truncated",
          "TailCorrections": False, "CutOffVDW": 6.0, "PseudoAtoms": [atom],
          "SelfInteractions": [interaction]}
    write_json(work / "force_field.json", ff)
    molecule = {"CriticalTemperature": 150.687, "CriticalPressure": 4863000.0,
                "AcentricFactor": 0.0, "Type": "rigid", "pseudoAtoms": [["X", [0., 0., 0.]]]}
    for name in ("A", "B"):
        write_json(work / (name + ".json"), molecule)


def matrix(work):
    path, = (work / "tmmc").glob("*statistics*.txt")
    labels, rows = [], []
    for line in path.read_text().splitlines():
        column = re.match(r"^# column (\d+):\s*(.+)", line)
        if column:
            labels.append(column[2].strip())
        if re.match(r"^\s*\d+\s+", line):
            rows.append(list(map(float, line.split())))
    if not rows or any(len(row) != len(labels) for row in rows):
        raise AssertionError(f"Empty or malformed matrix: {path}")
    return labels, rows


def get_deck(work, kind):
    fixed = kind in ("fixed", "cbmc-fixed", "coordinates")
    one_axis = kind in ("legacy", "cation")
    # z_i = V f_i / (k_B T). Set z_A=1/2, z_B=7/2 by choosing total z=4.
    temperature, length = 300.0, 24.0
    pressure = 4 * 1.380650324e-23 * temperature / (length * 1e-10) ** 3
    system = {"Type": "Box", "BoxLengths": [length] * 3, "ExternalTemperature": temperature,
              "ExternalPressure": pressure, "ChargeMethod": "None",
              "MacroStateUseBias": False}
    ideal_gas(work)
    names = ["A", "B"]
    components = [{"Name": name, "MoleculeDefinition": ".", "SwapConventionalProbability": 1.,
                   "FugacityCoefficient": 1., "MolFraction": fraction,
                   "CreateNumberOfMolecules": number}
                  for name, fraction, number in zip(names, (0.125, 0.875), (1, 2))]
    if "cbmc" in kind:
        for component in components:
            component["SwapProbability"] = component.pop("SwapConventionalProbability")
    if one_axis:
        components = components[:1]
        components[0].update(MolFraction=1.)
        if kind == "cation":
            components[0]["Type"] = "Cation"
        system.update(MacroStateMinimumNumberOfMolecules=0, MacroStateMaximumNumberOfMolecules=3)
    elif fixed:
        for component in components:
            n = component["CreateNumberOfMolecules"]
            component.update(MinMacrostate=n, MaxMacrostate=n)
    else:
        # The first component supplies only its upper bound and inherits the lower bound.
        system.update(MacroStateMinimumNumberOfMolecules=1, MacroStateMaximumNumberOfMolecules=3)
        components[0].update(MaxMacrostate=2)
        components[1].update(MinMacrostate=2, MaxMacrostate=3)
    if kind == "coordinates":
        # Load one molecule of each species, then have RASPA create the second B.
        system["RestartFileName"] = "coordinates.json"
        write_json(work / "coordinates.json", {"A": [[2., 3., 4.]], "B": [[14., 15., 16.]]})
        components[0]["CreateNumberOfMolecules"] = 0
        components[1]["CreateNumberOfMolecules"] = 1
    deck = {"SimulationType": "MonteCarloTransitionMatrix", "ForceField": ".",
            "RandomSeed": 74169, "PrintEvery": 5,
            "NumberOfInitializationCycles": 1, "NumberOfEquilibrationCycles": 1,
            "WriteBinaryRestartEvery": 111, "NumberOfCycles": 200,
            "Systems": [system], "Components": components}
    if kind == "adaptive":
        system["MacroStateUseBias"] = True
        deck["SampleTMMCEvery"] = 4
        deck["WriteTMMCEvery"] = 10
    write_json(work / "simulation.json", deck)


def verify(kind, work):
    if kind == "coordinates":
        kind = "fixed"
    elif kind == "cation":
        kind = "legacy"
    labels, rows = matrix(work)
    if kind == "legacy":
        assert [int(r[0]) for r in rows] == [0, 1, 2, 3], rows
        assert all(sum(r[1:4]) > 0 and all(math.isfinite(v) for v in r[:4]) for r in rows)
        return {"populations": [int(r[0]) for r in rows], "row_sums": [sum(r[1:4]) for r in rows]}
    required = ["N_0", "N_1", "CM", "CM_0^-", "CM_0^+", "CM_1^-", "CM_1^+",
                "histogram_overall", "average_potential_energy"]
    index = [labels.index(key) for key in required]
    grid = {}
    for row in rows:
        a, b, stay, am, ap, bm, bp, hist, energy = [row[i] for i in index]
        assert all(math.isfinite(v) for v in (a, b, stay, am, ap, bm, bp, hist, energy))
        assert all(v >= 0 for v in (stay, am, ap, bm, bp, hist))
        assert hist > 0 and sum((stay, am, ap, bm, bp)) > 0, row
        key = (int(a), int(b))
        assert key not in grid
        grid[key] = {"weights": [stay, am, ap, bm, bp], "histogram": hist, "energy": energy}
    if kind in ("fixed", "cbmc-fixed"):
        assert set(grid) == {(1, 2)}, list(grid)
        node = grid[(1, 2)]
        assert min(node["weights"][1:]) > 0, node
        # Fixed integer N=(1,2), equality of forward and backward proposal rates.
        # For axis i, P+ = min(1,z_i/(N_i+1)), P- = min(1,N_i/z_i).
        # Columns allow independent reconstruction of ALL four attempt counts, including
        # physically allowed probes rejected by the hard windows, plus the expected stay row.
        acceptance = [1., 0.25, 2. / 3.5, 1.]
        directions = node["weights"][1:]
        attempts = [v / p for v, p in zip(directions, acceptance)]
        assert max(abs(v - round(v)) for v in attempts) < 1e-4, attempts
        assert math.isclose(sum(attempts), sum(node["weights"]), rel_tol=1e-8), node
        assert abs(node["energy"]) < 1e-8, node
        assert min(attempts) / sum(attempts) > .15, attempts
        node["analytical_probabilities"] = acceptance
        node["reconstructed_attempts"] = list(map(round, attempts))
    else:
        assert set(grid) == set(itertools.product((1, 2), (2, 3))), list(grid)
        assert all(min(node["weights"][1:]) > 0 for node in grid.values())
    if kind in ("mixed", "adaptive"):
        total_hist = sum(row["histogram"] for row in grid.values())
        grand_weights = {(a, b): .5**a * 3.5**b / (math.factorial(a) * math.factorial(b))
                         for a, b in grid}
        normalization = sum(grand_weights.values())
        for (a, b), node in grid.items():
            expected = grand_weights[a, b] / normalization
            observed = node["histogram"] / total_hist
            # The deterministic stream is long enough to resolve the four truncated Poisson masses.
            if kind == "mixed":
                assert abs(expected - observed) < .055, ((a, b), observed, expected)
            rates = [min(1., a / .5), min(1., .5 / (a + 1)),
                     min(1., b / 3.5), min(1., 3.5 / (b + 1))]
            attempts = [weight / rate for weight, rate in zip(node["weights"][1:], rates)]
            assert max(abs(count - round(count)) for count in attempts) < 1e-4, ((a, b), attempts)
            assert math.isclose(sum(attempts), sum(node["weights"]), rel_tol=1e-8), node
            node["ideal_probability"] = expected
            node["sampled_probability"] = observed
    result = {str(key): value for key, value in sorted(grid.items())}
    if kind == "adaptive":
        path, = (work / "tmmc").glob("*statistics*.txt")
        match = re.search(r"^# bias updated: (\d+) times$", path.read_text(), re.M)
        assert match and int(match[1]) > 0
        biases = [row[labels.index("bias")] for row in rows]
        assert all(map(math.isfinite, biases)) and max(biases) > min(biases)
        result["adaptive"] = {"updates": int(match[1]), "bias_range": max(biases) - min(biases)}
    path, = (work / "output").glob("output_*.s0.json")
    move_stats = read_json(path)["output"]["MCMoveStatistics"]["system"]
    attempts_by_axis = []
    for name in ("A", "B"):
        swap = move_stats[name]["Swap (CBMC)" if "cbmc" in kind else "Swap"]
        attempts = swap["total"][:2]
        accepted = swap["accepted"][:2]
        assert min(attempts) > 0, swap
        if kind.endswith("fixed"):
            assert accepted == [0, 0], swap
        else:
            assert min(accepted) > 0, swap
        result[name] = {"swap_attempts": attempts, "accepted": accepted}
        attempts_by_axis.append(attempts)
    if kind in ("fixed", "cbmc-fixed"):
        # Native swap counters order insertion first, deletion second.
        counts = [int(x) for each in attempts_by_axis for x in reversed(each)]
        assert counts == grid[1, 2]["reconstructed_attempts"], (counts, grid[1, 2])
    return result


def native(binary, work, env):
    with (work / "console.log").open("w") as stream:
        proc = subprocess.run([str(binary)], cwd=work, env=env, stdout=stream, stderr=subprocess.STDOUT)
    if proc.returncode:
        raise RuntimeError(f"{work.name}: native exit {proc.returncode}, see {work / 'console.log'}")


def restart_check(kind, binary, first, env):
    # Resume in a clean directory with a different input RNG seed and only the explicit checkpoint.
    # This catches a silent fallback to fresh sampling, and stale matrix files cannot mask a lost row.
    resumed = first.with_name(first.name + "-resumed")
    resumed.mkdir()
    for src in first.glob("*.json"):
        shutil.copyfile(src, resumed / src.name)
    checkpoint = first / "restart_data.bin"
    assert checkpoint.stat().st_size > 0
    shutil.copyfile(checkpoint, resumed / "source.bin")
    settings = read_json(resumed / "simulation.json")
    settings.update(RestartFromBinaryFile=True, BinaryRestartFileName="source.bin", RandomSeed=50235)
    # The archived A count is authoritative; the input-only count lies below its lower bound.
    settings["Components"][0]["CreateNumberOfMolecules"] = 0
    write_json(resumed / "simulation.json", settings)
    native(binary, resumed, env)
    reference = matrix(first)
    continuation = matrix(resumed)
    assert reference == continuation, f"Resumed matrix differs: {resumed}"
    verify(kind, resumed)
    return {"source_checkpoint_sha256": hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
            "identical_rows": len(reference[1]), "input_seed_changed": True}


def reject_extra_coordinate(binary, first, env):
    bad = first.with_name("coordinates-outside")
    bad.mkdir()
    for path in first.glob("*.json"):
        shutil.copyfile(path, bad / path.name)
    deck = read_json(bad / "simulation.json")
    deck["Components"][0]["CreateNumberOfMolecules"] = 1
    write_json(bad / "simulation.json", deck)
    run = subprocess.run([str(binary)], cwd=bad, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (bad / "console.log").write_text(run.stdout)
    assert run.returncode != 0 and "component A initial population 2" in run.stdout, run.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cases", nargs="+", choices=("legacy", "fixed", "cbmc-fixed", "mixed", "adaptive", "coordinates", "cation"),
                        default=("legacy", "fixed", "cbmc-fixed", "mixed", "adaptive", "coordinates", "cation"))
    args = parser.parse_args()
    binary, source = args.binary.resolve(), args.source.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    results = {"binary": str(binary), "sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
               "commit": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
               "cases": {}}
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    for kind in args.cases:
        work = args.output / kind
        work.mkdir()
        get_deck(work, kind)
        native(binary, work, env)
        results["cases"][kind] = verify(kind, work)
        if kind == "coordinates":
            reject_extra_coordinate(binary, work, env)
            results["cases"][kind]["out_of_range_rejected"] = True
        if kind == "mixed":
            results["cases"][kind]["restart"] = restart_check(kind, binary, work, env)
        write_json(args.output / "results.json", results)
        print(kind, json.dumps(results["cases"][kind]), flush=True)


if __name__ == "__main__":
    main()
