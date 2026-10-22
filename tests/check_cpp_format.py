#!/usr/bin/env python3
"""Check changed C++ lines against the repository style without rewriting upstream files."""

import argparse
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, help="Reviewed base commit")
    parser.add_argument("--binary", default="clang-format", help="Use the pinned clang-format 20 tool")
    args = parser.parse_args()
    result = subprocess.run(
        ["git", "clang-format", "--binary", args.binary, "--extensions", "cpp,ixx", "--diff", args.base],
        cwd=Path(__file__).resolve().parents[1],
        text=True, capture_output=True, check=False,
    )
    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")
    if "diff --git " in result.stdout:
        print(result.stdout, end="")
        return 1
    if result.returncode != 0:
        return result.returncode or 1
    print(result.stdout.strip() or "Changed C++ lines match repository style.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
