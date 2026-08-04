#!/usr/bin/env python3
"""Run this input-driven example with an installed or explicitly selected RASPA3."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys


example_directory = Path(__file__).resolve().parent
executable = os.environ.get("RASPA3_EXECUTABLE") or shutil.which("raspa3")
if executable is None:
    raise SystemExit(
        "RASPA3 was not found. Load the raspa3-mbx module or set "
        "RASPA3_EXECUTABLE=/path/to/raspa3."
    )

subprocess.run([executable, *sys.argv[1:]], cwd=example_directory, check=True)
