#!/bin/sh
set -eu

example_directory=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
if [ -n "${RASPA3_EXECUTABLE:-}" ]; then
  executable=$(command -v "$RASPA3_EXECUTABLE") || {
    echo "Cannot find RASPA3_EXECUTABLE: $RASPA3_EXECUTABLE" >&2
    exit 127
  }
elif command -v raspa3 >/dev/null 2>&1; then
  executable=$(command -v raspa3)
else
  echo "Load an MBX-enabled RASPA3 installation or set RASPA3_EXECUTABLE=/absolute/path/to/raspa3." >&2
  exit 127
fi
case $executable in
  /*) ;;
  *) executable=$(CDPATH= cd -- "$(dirname -- "$executable")" && pwd)/$(basename -- "$executable") ;;
esac

# Keep all input lookup inside this example, including after copying it elsewhere.
export RASPA_DIR=$example_directory
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
cd "$example_directory"
exec "$executable" "$@"
