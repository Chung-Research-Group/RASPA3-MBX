#!/bin/sh
set -eu

# Resolve a user-supplied executable before changing the working directory.
if [ -n "${RASPA3_EXECUTABLE:-}" ]; then
  raspa3_executable=$RASPA3_EXECUTABLE
elif command -v raspa3 >/dev/null 2>&1; then
  raspa3_executable=$(command -v raspa3)
else
  echo "RASPA3-MBX was not found. Build with BUILD_MBX=ON, then set" >&2
  echo "RASPA3_EXECUTABLE=/absolute/path/to/raspa3 or add raspa3 to PATH." >&2
  exit 127
fi

case $raspa3_executable in
  /*) ;;
  */*) raspa3_executable=$(CDPATH= cd -- "$(dirname -- "$raspa3_executable")" && pwd)/$(basename -- "$raspa3_executable") ;;
  *) raspa3_executable=$(command -v "$raspa3_executable") ;;
esac

if [ ! -x "$raspa3_executable" ]; then
  echo "RASPA3 executable is not executable: $raspa3_executable" >&2
  exit 126
fi

script_directory=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$script_directory"
export RASPA_DIR=$script_directory
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
exec "$raspa3_executable" "$@"
