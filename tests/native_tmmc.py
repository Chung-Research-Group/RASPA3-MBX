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


def mbx_gas(work, source):
    # Use the same reviewed molecule coordinates and pseudo atoms as the native MBX examples.
    examples = source / "examples/raspa3_mbx"
    sources = [examples / "5_water_npt", examples / "4_gcmc_co2_adsorption"]
    pseudo_atoms, interactions = [], []
    for name, case in zip(("h2o", "co2"), sources):
        molecule = read_json(case / (name + ".json"))
        write_json(work / (name + ".json"), molecule)
        needed = {atom[0] for atom in molecule["pseudoAtoms"]}
        ff = read_json(case / "force_field.json")
        pseudo_atoms.extend(atom for atom in ff["PseudoAtoms"] if atom["name"] in needed)
        interactions.extend(atom for atom in ff["SelfInteractions"] if atom["name"] in needed)
    write_json(work / "force_field.json", {"MixingRule": "Lorentz-Berthelot",
               "TruncationMethod": "truncated", "TailCorrections": False, "CutOffVDW": 9.0,
               "PseudoAtoms": pseudo_atoms, "SelfInteractions": interactions})
    write_json(work / "mbx.json", read_json(sources[0] / "mbx.json"))


def matrix(work, pattern="*statistics*.txt"):
    path, = (work / "tmmc").glob(pattern)
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


def get_deck(work, source, modern, kind):
    use_mbx = kind.startswith("mbx")
    fixed = kind in ("fixed", "cbmc-fixed", "mbx-fixed", "coordinates")
    one_axis = kind in ("legacy", "cation")
    # z_i = V f_i / (k_B T). Set z_A=1/2, z_B=7/2 by choosing total z=4.
    temperature, length = 300.0, 24.0
    pressure = 4 * 1.380650324e-23 * temperature / (length * 1e-10) ** 3
    system = {"Type": "Box", "BoxLengths": [length] * 3, "ExternalTemperature": temperature,
              "ExternalPressure": pressure, "ChargeMethod": "Ewald" if use_mbx else "None",
              "MacroStateUseBias": False}
    if use_mbx:
        mbx_gas(work, source)
        names = ["h2o", "co2"]
        system.update(UseMBX=True, MBXSettingsFile="mbx.json", PrintEnergyTerms=True)
    else:
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
            "WriteBinaryRestartEvery": 7 if use_mbx else 111,
            "NumberOfProductionCycles" if modern else "NumberOfCycles": 12 if use_mbx else 200,
            "Systems": [system], "Components": components}
    if kind == "adaptive":
        system["MacroStateUseBias"] = True
        deck["SampleTMMCEvery"] = 4
        deck["WriteTMMCEvery"] = 10
    write_json(work / "simulation.json", deck)


def closure(work):
    path, = (work / "output").glob("output_*.s0.json")
    output = read_json(path)["output"]
    total, mbx = "Total potential energy [K]", "MBX energy [K]"
    running, recomputed = output["runningEnergies"], output["recomputedEnergies"]
    tol = 1e-4  # K, about 2e-7 kcal/mol, for three-site molecules and the Ewald solver.
    for key in (total, mbx):
        assert abs(recomputed[key]) > 1e-10, (key, recomputed[key])
        delta = running[key] - recomputed[key]
        assert abs(delta) < tol and abs(output["drift"][key] - delta) < tol, (key, delta)
    return {"running_total_K": running[total], "recomputed_total_K": recomputed[total],
            "mbx_K": recomputed[mbx], "drift_K": output["drift"][total]}


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
    if kind in ("fixed", "cbmc-fixed", "mbx-fixed"):
        assert set(grid) == {(1, 2)}, list(grid)
        node = grid[(1, 2)]
        assert min(node["weights"][1:]) > 0, node
        if kind in ("fixed", "cbmc-fixed"):
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
            assert abs(node["energy"]) > 1e-9, node
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
    names = ("h2o", "co2") if kind.startswith("mbx") else ("A", "B")
    attempts_by_axis = []
    for name in names:
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
    if kind.startswith("mbx"):
        result["energy_check"] = closure(work)
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


def parallel_coordinates(binary, source, modern, work, env):
    get_deck(work, source, modern, "legacy")
    deck = read_json(work / "simulation.json")
    deck.update(SimulationType="ParallelTMMC", NumberOfWindows=2, ReweightingNumberOfPressures=5,
                NumberOfProductionCycles=30, WriteBinaryRestartEvery=8)
    system = deck["Systems"][0]
    temperature = system.pop("ExternalTemperature")
    system.update(RestartFileName="coordinates.json", ExternalTemperatures=[temperature],
                  MacroStateMinimumNumberOfMolecules=2, MacroStateMaximumNumberOfMolecules=4)
    deck["Components"][0]["CreateNumberOfMolecules"] = 0
    write_json(work / "coordinates.json", {"A": [[2., 3., 4.]]})
    write_json(work / "simulation.json", deck)
    native(binary, work, env)
    files = sorted((work / "tmmc").glob("*.parallel_tmmc.txt"))
    assert len(files) == 2, files
    windows = []
    for path, expected in zip(files, ([2, 3], [3, 4])):
        _, rows = matrix(work, path.name)
        assert [row[0] for row in rows] == expected and all(row[-1] > 0 for row in rows), rows
        assert math.isclose(sum(sum(row[1:4]) for row in rows), sum(row[-1] for row in rows), rel_tol=1e-12)
        windows.append(expected)
    bad = work.with_name("parallel-coordinates-outside")
    resumed = work.with_name("parallel-coordinates-resumed")
    deck["Components"][0]["CreateNumberOfMolecules"] = 2  # One coordinate plus two created exceeds the first lower edge.
    for directory in (bad, resumed):
        directory.mkdir()
        for path in work.glob("*.json"):
            shutil.copyfile(path, directory / path.name)
        write_json(directory / "simulation.json", deck)
    proc = subprocess.run([str(binary)], cwd=bad, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (bad / "console.log").write_text(proc.stdout)
    assert proc.returncode != 0 and "initial population (3)" in proc.stdout, proc.stdout
    shutil.copyfile(work / "restart_data.bin", resumed / "source.bin")
    deck.update(RestartFromBinaryFile=True, BinaryRestartFileName="source.bin", RandomSeed=50235)
    write_json(resumed / "simulation.json", deck)
    native(binary, resumed, env)
    for path in files:
        assert path.read_bytes() == (resumed / "tmmc" / path.name).read_bytes()
    return {"walker_windows": windows, "excess_initial_population_rejected": True, "identical_restarted_matrices": 2}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cases", nargs="+", choices=("legacy", "fixed", "cbmc-fixed", "mixed", "adaptive", "coordinates", "cation", "parallel-coordinates",
                                                        "mbx-fixed", "mbx-mixed", "mbx-cbmc-mixed"),
                        default=("legacy", "fixed", "cbmc-fixed", "mixed", "adaptive", "coordinates", "cation", "parallel-coordinates"))
    args = parser.parse_args()
    binary, source = args.binary.resolve(), args.source.resolve()
    modern = '"NumberOfProductionCycles"' in (source / "src/raspakit/input_reader.cpp").read_text()
    args.output.mkdir(parents=True, exist_ok=True)
    results = {"binary": str(binary), "sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
               "commit": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
               "cases": {}}
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    for kind in args.cases:
        work = args.output / kind
        work.mkdir()
        if kind == "parallel-coordinates":
            results["cases"][kind] = parallel_coordinates(binary, source, modern, work, env)
        else:
            get_deck(work, source, modern, kind)
            native(binary, work, env)
            results["cases"][kind] = verify(kind, work)
            if kind == "coordinates":
                reject_extra_coordinate(binary, work, env)
                results["cases"][kind]["out_of_range_rejected"] = True
            if kind in ("mixed", "mbx-mixed"):
                results["cases"][kind]["restart"] = restart_check(kind, binary, work, env)
        write_json(args.output / "results.json", results)
        print(kind, json.dumps(results["cases"][kind]), flush=True)


if __name__ == "__main__":
    main()
