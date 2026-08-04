#!/bin/sh
set -eu

usage()
{
  echo "Usage: $0 EXAMPLE_DIRECTORY [raspa3 arguments...]" >&2
  echo "Example: $0 2_widom_zero_loading_qst" >&2
}

if [ "$#" -lt 1 ]; then
  usage
  exit 2
fi

script_directory=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
example_directory=$script_directory/$1
shift

# Component lookup falls back to RASPA_DIR. MBX's monomer identifier is the
# case-sensitive lowercase name "co2", whose definition is stored here.
export RASPA_DIR=$script_directory

if [ ! -f "$example_directory/simulation.json" ]; then
  echo "No simulation.json found under: $example_directory" >&2
  usage
  exit 2
fi

if [ -n "${RASPA3_EXECUTABLE:-}" ]; then
  raspa3_executable=$RASPA3_EXECUTABLE
elif command -v raspa3 >/dev/null 2>&1; then
  raspa3_executable=$(command -v raspa3)
else
  repository_root=$(CDPATH= cd -- "$script_directory/../.." && pwd)
  if [ -x "$repository_root/build-mbx-release-avx2/app/raspa3" ]; then
    raspa3_executable=$repository_root/build-mbx-release-avx2/app/raspa3
  elif [ -x "$repository_root/build-mbx/app/raspa3" ]; then
    raspa3_executable=$repository_root/build-mbx/app/raspa3
  else
    echo "RASPA3 was not found." >&2
    echo "Load the raspa3-mbx module or set RASPA3_EXECUTABLE=/path/to/raspa3." >&2
    exit 127
  fi
fi

case $raspa3_executable in
  /*) ;;
  *) raspa3_executable=$(CDPATH= cd -- "$(dirname -- "$raspa3_executable")" && pwd)/$(basename -- "$raspa3_executable") ;;
esac

echo "Running $example_directory with $raspa3_executable"
echo "Using shared definitions from $RASPA_DIR"
cd "$example_directory"
exec "$raspa3_executable" "$@"
