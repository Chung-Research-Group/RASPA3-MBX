#!/bin/sh
set -eu

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 EXAMPLE_DIRECTORY [raspa3 arguments...]" >&2
  exit 2
fi
script_directory=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
case $1 in
  */*|.|..) echo "Choose a numbered example directory in $script_directory." >&2; exit 2 ;;
esac
example_directory=$script_directory/$1
shift
if [ ! -x "$example_directory/run.sh" ] || [ ! -f "$example_directory/simulation.json" ]; then
  echo "No runnable example found: $example_directory" >&2
  exit 2
fi
exec "$example_directory/run.sh" "$@"
