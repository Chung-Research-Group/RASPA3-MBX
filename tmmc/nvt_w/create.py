import numpy as np
import json
import sys
import shutil
import os

limits_0 = [4 * i for i in range(121)]
limits_1 = [4 * i for i in range(21)]

for i0, n0 in enumerate(limits_0):
    for i1, n1 in enumerate(limits_1):

        simdir = f"{n0}_{n1}"
        print(simdir)

        os.makedirs(simdir, exist_ok=True)

        for name in os.listdir("base"):
            src = os.path.join("base", name)
            dst = os.path.join(simdir, name)
            if os.path.isfile(src):
                shutil.copy2(src, dst)

        with open(f"{simdir}/simulation.json") as f:
            settings = json.load(f)

        settings["SimulationType"] = "MonteCarloTransitionMatrix"
        settings["Components"][0]["MinMacrostate"] = n0
        settings["Components"][0]["MaxMacrostate"] = n0
        settings["Components"][0]["CreateNumberOfMolecules"] = n0
        settings["Components"][1]["MinMacrostate"] = n1
        settings["Components"][1]["MaxMacrostate"] = n1
        settings["Components"][1]["CreateNumberOfMolecules"] = n1
        settings["NumberOfCycles"] = 50000
        settings["NumberOfInitializationCycles"] = 10000
        settings["NumberOfEquilibrationCycles"] = 10000
        settings["PrintEvery"] = 1000
        settings['Systems'][0]['ExternalPressure'] = 9.7e4

        with open(f"{simdir}/simulation.json", "w") as f:
            json.dump(settings, f, indent=4)